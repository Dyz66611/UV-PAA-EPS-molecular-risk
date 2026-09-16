from __future__ import annotations

import argparse
import itertools
import json
import re
import warnings
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    brier_score_loss,
    f1_score,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

RANDOM_STATE = 42
PROBABILITY_THRESHOLD = 0.5

# Exact identifiers or non-predictor metadata. Comparison is case-insensitive.
BLOCKED_EXACT = {
    "global_pair_id",
    "pair_id",
    "record_id",
    "row_id",
    "chem_key",
    "number",
    "record_number",
    "smiles",
    "canonical_smiles",
    "chemical",
    "chemical_name",
    "name",
    "stage",
    "treatment_stage",
    "treatment",
    "time",
    "time_point",
    "time_interval",
    "interval",
}

# Any feature containing one of these strings is an outcome, an ECOSAR-derived
# variable, or treatment metadata and is blocked to prevent leakage.
BLOCKED_SUBSTRINGS = {
    "risk_label",
    "risk_binary",
    "risk_change",
    "toxicity",
    "toxic_",
    "ecosar",
    "lc50",
    "ec50",
    "chv",
    "effect_concentration",
    "risk_conc",
    "log10_mgl",
    "endpoint",
    "change_direction",
    "class_label",
    "response_label",
    "outcome",
}

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="XGBoost and SHAP analysis for binary product-risk prediction."
    )
    parser.add_argument("--input_file", type=Path, required=True)
    parser.add_argument("--out_dir", type=Path, required=True)
    parser.add_argument(
        "--sheet_name",
        default=None,
        help="Input sheet. Default: model_input when present, otherwise first sheet.",
    )
    parser.add_argument(
        "--label_text_col", default="Target_risk_label_binary"
    )
    parser.add_argument(
        "--label_num_col", default="y_target_risk_binary"
    )
    parser.add_argument(
        "--feature_prefixes",
        default="Source_,Target_",
        help="Comma-separated prefixes used for automatic feature selection.",
    )
    parser.add_argument(
        "--feature_file",
        type=Path,
        default=None,
        help="Optional TXT/CSV/XLSX file containing an explicit feature whitelist.",
    )
    parser.add_argument("--test_size", type=float, default=0.20)
    parser.add_argument("--random_state", type=int, default=RANDOM_STATE)
    parser.add_argument("--threshold", type=float, default=PROBABILITY_THRESHOLD)
    parser.add_argument("--top_n", type=int, default=10)
    parser.add_argument("--interaction_sample_size", type=int, default=200)
    parser.add_argument("--n_response_pairs", type=int, default=8)
    parser.add_argument("--response_grid_size", type=int, default=70)
    parser.add_argument(
        "--shap_scope",
        choices=("all", "test"),
        default="all",
        help="Rows used for SHAP interpretation. Default matches the Methods text.",
    )
    return parser.parse_args()

def read_input(path: Path, sheet_name: str | None) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        xls = pd.ExcelFile(path)
        if sheet_name is None:
            sheet_name = "model_input" if "model_input" in xls.sheet_names else xls.sheet_names[0]
        if sheet_name not in xls.sheet_names:
            raise ValueError(
                f"Sheet '{sheet_name}' not found. Available sheets: {xls.sheet_names}"
            )
        df = pd.read_excel(path, sheet_name=sheet_name)
        print(f"Loaded sheet: {sheet_name}")
    elif suffix == ".csv":
        df = pd.read_csv(path)
    elif suffix in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t")
    else:
        raise ValueError("Input must be an XLSX, XLS, CSV, TSV, or TXT file.")
    df.columns = [str(c).strip() for c in df.columns]
    if df.columns.duplicated().any():
        duplicates = df.columns[df.columns.duplicated()].tolist()
        raise ValueError(f"Duplicate column names are not allowed: {duplicates}")
    print(f"Input shape: {df.shape}")
    return df

def read_feature_whitelist(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Feature whitelist not found: {path}")
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        values = pd.read_excel(path).iloc[:, 0]
    elif suffix == ".csv":
        values = pd.read_csv(path).iloc[:, 0]
    else:
        values = pd.read_csv(path, header=None, comment="#").iloc[:, 0]
    return list(dict.fromkeys(str(x).strip() for x in values.dropna() if str(x).strip()))

def remove_molecule_prefix(name: str, prefixes: tuple[str, ...]) -> str:
    lowered = name.casefold()
    for prefix in prefixes:
        if lowered.startswith(prefix.casefold()):
            return name[len(prefix):]
    return name

def leakage_reason(column: str, prefixes: tuple[str, ...]) -> str | None:
    """Return the reason a column is unsafe, or None when it is permissible."""
    lowered = column.casefold().strip()
    base = remove_molecule_prefix(column, prefixes).casefold().strip()
    if lowered in BLOCKED_EXACT or base in BLOCKED_EXACT:
        return "identifier_or_metadata"
    if lowered.endswith("_id") or base.endswith("_id"):
        return "identifier"
    if lowered.startswith("id_") or base.startswith("id_"):
        return "identifier"
    for token in BLOCKED_SUBSTRINGS:
        if token in lowered:
            return f"blocked_token:{token}"
    # Stage/time fields may be embedded after the Source_/Target_ prefix.
    if re.search(r"(^|_)(stage|treatment|interval|timepoint|time_point)($|_)", base):
        return "treatment_or_stage"
    return None

def coerce_numeric_feature(series: pd.Series, name: str) -> pd.Series:
    """Convert a candidate feature to numeric without silently accepting text."""
    converted = pd.to_numeric(series, errors="coerce")
    original_nonmissing = series.notna() & series.astype(str).str.strip().ne("")
    if original_nonmissing.any():
        success_fraction = converted[original_nonmissing].notna().mean()
        if success_fraction < 0.90:
            raise ValueError(
                f"Feature '{name}' is not reliably numeric "
                f"({success_fraction:.1%} of non-missing values converted)."
            )
    return converted.replace([np.inf, -np.inf], np.nan)

def select_features(
    df: pd.DataFrame,
    prefixes: tuple[str, ...],
    feature_file: Path | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the numeric predictor matrix and a complete selection audit."""
    if feature_file is not None:
        candidates = read_feature_whitelist(feature_file)
        missing = [c for c in candidates if c not in df.columns]
        if missing:
            raise ValueError(f"Whitelisted features absent from input: {missing}")
        selection_source = "explicit_whitelist"
    else:
        candidates = [c for c in df.columns if c.startswith(prefixes)]
        selection_source = "automatic_prefix_selection"
        if not candidates:
            raise ValueError(
                "No candidate features matched the prefixes "
                f"{prefixes}. Supply --feature_file or change --feature_prefixes."
            )

    selected: dict[str, pd.Series] = {}
    audit_rows: list[dict[str, str]] = []
    blocked_whitelist: list[str] = []
    for column in candidates:
        reason = leakage_reason(column, prefixes)
        if reason is not None:
            audit_rows.append(
                {"feature": column, "status": "excluded", "reason": reason,
                 "selection_source": selection_source}
            )
            if feature_file is not None:
                blocked_whitelist.append(column)
            continue
        try:
            numeric = coerce_numeric_feature(df[column], column)
        except ValueError as exc:
            audit_rows.append(
                {"feature": column, "status": "excluded", "reason": str(exc),
                 "selection_source": selection_source}
            )
            if feature_file is not None:
                raise
            continue
        if numeric.notna().sum() == 0:
            audit_rows.append(
                {"feature": column, "status": "excluded", "reason": "all_missing",
                 "selection_source": selection_source}
            )
            continue
        selected[column] = numeric
        audit_rows.append(
            {"feature": column, "status": "selected", "reason": "valid_numeric_feature",
             "selection_source": selection_source}
        )

    if blocked_whitelist:
        raise ValueError(
            "The explicit whitelist contains leakage-prone columns: "
            f"{blocked_whitelist}. Remove them before rerunning."
        )
    if not selected:
        raise ValueError("No valid numerical molecular features remained after screening.")
    return pd.DataFrame(selected, index=df.index), pd.DataFrame(audit_rows)

def prepare_binary_labels(
    df: pd.DataFrame, label_text_col: str, label_num_col: str
) -> tuple[pd.Series, pd.Series]:
    """Create and validate 0/1 labels; Unknown records remain missing."""
    text_map = {"Low-risk": 0, "High-risk": 1}
    if label_text_col not in df.columns and label_num_col not in df.columns:
        raise ValueError(
            f"Neither label column was found: '{label_text_col}' or '{label_num_col}'."
        )

    mapped_text = None
    if label_text_col in df.columns:
        normalized = df[label_text_col].astype("string").str.strip()
        mapped_text = normalized.map(text_map)

    numeric = None
    if label_num_col in df.columns:
        numeric = pd.to_numeric(df[label_num_col], errors="coerce")
        numeric = numeric.where(numeric.isin([0, 1]))

    if mapped_text is not None and numeric is not None:
        both = mapped_text.notna() & numeric.notna()
        mismatch = both & mapped_text.ne(numeric)
        if mismatch.any():
            examples = df.loc[mismatch, [label_text_col, label_num_col]].head().to_dict("records")
            raise ValueError(f"Text and numeric labels disagree. Examples: {examples}")
        labels = mapped_text.combine_first(numeric)
    elif mapped_text is not None:
        labels = mapped_text
    else:
        labels = numeric

    keep = labels.isin([0, 1])
    return labels.astype("Float64"), keep

def fit_model(X_train: pd.DataFrame, y_train: pd.Series, random_state: int) -> XGBClassifier:
    model = XGBClassifier(
        n_estimators=400,
        max_depth=5,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        objective="binary:logistic",
        eval_metric="logloss",
        random_state=random_state,
        n_jobs=-1,
    )
    model.fit(X_train, y_train)
    return model

def normalize_shap_output(values, positive_class: int = 1) -> np.ndarray:
    """Normalize common TreeExplainer output formats to rows x features."""
    if isinstance(values, list):
        values = values[positive_class] if len(values) > 1 else values[0]
    if hasattr(values, "values"):
        values = values.values
    array = np.asarray(values)
    if array.ndim == 3:
        # Newer APIs may return rows x features x classes.
        array = array[:, :, positive_class]
    if array.ndim != 2:
        raise ValueError(f"Unexpected SHAP value shape: {array.shape}")
    return array

def normalize_interaction_output(values, positive_class: int = 1) -> np.ndarray:
    if isinstance(values, list):
        values = values[positive_class] if len(values) > 1 else values[0]
    if hasattr(values, "values"):
        values = values.values
    array = np.asarray(values)
    if array.ndim == 4:
        array = array[:, :, :, positive_class]
    if array.ndim != 3:
        raise ValueError(f"Unexpected SHAP interaction shape: {array.shape}")
    return array

def save_performance_plots(y_test: pd.Series, probabilities: np.ndarray, out_dir: Path) -> None:
    fpr, tpr, _ = roc_curve(y_test, probabilities)
    precision, recall, _ = precision_recall_curve(y_test, probabilities)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.plot(fpr, tpr, color="#3575b5", linewidth=2)
    ax.plot([0, 1], [0, 1], "--", color="0.65", linewidth=1)
    ax.set(xlabel="False-positive rate", ylabel="True-positive rate", title="Test-set ROC curve")
    fig.tight_layout()
    fig.savefig(out_dir / "test_roc_curve.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.plot(recall, precision, color="#d95f70", linewidth=2)
    ax.set(xlabel="Recall", ylabel="Precision", title="Test-set precision-recall curve")
    fig.tight_layout()
    fig.savefig(out_dir / "test_precision_recall_curve.png", dpi=300)
    plt.close(fig)

def save_global_shap_plots(
    X_shap: pd.DataFrame,
    shap_values: np.ndarray,
    summary: pd.DataFrame,
    out_dir: Path,
    top_n: int,
) -> None:
    shown = summary.head(top_n).sort_values("mean_abs_shap")
    fig, ax = plt.subplots(figsize=(7.2, max(4.5, 0.38 * len(shown) + 1.5)))
    ax.barh(shown["feature"], shown["mean_abs_shap"], color="#4c78a8")
    ax.set(xlabel="Mean |SHAP value|", ylabel="", title="Global feature importance")
    fig.tight_layout()
    fig.savefig(out_dir / "shap_global_importance.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    plt.figure(figsize=(8.0, max(5.0, 0.40 * min(top_n, X_shap.shape[1]) + 1.5)))
    shap.summary_plot(
        shap_values,
        X_shap,
        max_display=top_n,
        show=False,
        plot_size=None,
    )
    plt.tight_layout()
    plt.savefig(out_dir / "shap_beeswarm.png", dpi=300, bbox_inches="tight")
    plt.close()

def strongest_pairs_from_interactions(
    interaction_matrix: pd.DataFrame, n_pairs: int
) -> list[tuple[str, str, float]]:
    pairs = []
    names = interaction_matrix.columns.tolist()
    for i, j in itertools.combinations(range(len(names)), 2):
        pairs.append((names[i], names[j], float(interaction_matrix.iloc[i, j])))
    return sorted(pairs, key=lambda item: item[2], reverse=True)[:n_pairs]

def fallback_pairs(
    mean_abs_shap: pd.Series, n_pairs: int
) -> list[tuple[str, str, float]]:
    top = mean_abs_shap.head(min(6, len(mean_abs_shap))).index.tolist()
    pairs = [
        (a, b, float(mean_abs_shap[a] + mean_abs_shap[b]))
        for a, b in itertools.combinations(top, 2)
    ]
    return sorted(pairs, key=lambda item: item[2], reverse=True)[:n_pairs]

def main() -> None:
    args = parse_args()
    if not 0 < args.test_size < 1:
        raise ValueError("--test_size must be between 0 and 1.")
    if not 0 < args.threshold < 1:
        raise ValueError("--threshold must be between 0 and 1.")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    df_original = read_input(args.input_file, args.sheet_name)
    labels, keep = prepare_binary_labels(
        df_original, args.label_text_col, args.label_num_col
    )
    n_excluded = int((~keep).sum())
    df = df_original.loc[keep].copy().reset_index(names="input_row_index")
    y = labels.loc[keep].astype(int).reset_index(drop=True)
    print(f"Retained labeled records: {len(df)}; excluded Unknown/invalid: {n_excluded}")

    class_counts = y.value_counts().sort_index()
    if set(class_counts.index) != {0, 1}:
        raise ValueError(f"Both classes are required. Observed counts: {class_counts.to_dict()}")
    if class_counts.min() < 2:
        raise ValueError(
            "At least two observations per class are required for stratified splitting. "
            f"Observed counts: {class_counts.to_dict()}"
        )

    prefixes = tuple(x.strip() for x in args.feature_prefixes.split(",") if x.strip())
    if not prefixes and args.feature_file is None:
        raise ValueError("Provide at least one feature prefix or an explicit --feature_file.")
    X, feature_audit = select_features(df, prefixes, args.feature_file)

    row_positions = np.arange(len(df))
    train_pos, test_pos = train_test_split(
        row_positions,
        test_size=args.test_size,
        random_state=args.random_state,
        stratify=y,
    )
    X_train_raw = X.iloc[train_pos].copy()
    X_test_raw = X.iloc[test_pos].copy()
    y_train = y.iloc[train_pos].copy()
    y_test = y.iloc[test_pos].copy()

    # Remove features that cannot be learned from the training subset. This
    # decision uses training data only.
    train_all_missing = X_train_raw.columns[X_train_raw.isna().all()].tolist()
    train_nonmissing_nunique = X_train_raw.nunique(dropna=True)
    train_constant = train_nonmissing_nunique[train_nonmissing_nunique <= 1].index.tolist()
    dropped_after_split = sorted(set(train_all_missing + train_constant))
    if dropped_after_split:
        feature_audit.loc[
            feature_audit["feature"].isin(dropped_after_split), ["status", "reason"]
        ] = ["excluded", "all_missing_or_constant_in_training"]
        X_train_raw = X_train_raw.drop(columns=dropped_after_split)
        X_test_raw = X_test_raw.drop(columns=dropped_after_split)
        X = X.drop(columns=dropped_after_split)
    if X_train_raw.shape[1] == 0:
        raise ValueError("No predictors remained after training-subset quality checks.")

    imputer = SimpleImputer(strategy="median")
    feature_names = X_train_raw.columns.tolist()
    X_train = pd.DataFrame(
        imputer.fit_transform(X_train_raw), columns=feature_names, index=train_pos
    )
    X_test = pd.DataFrame(
        imputer.transform(X_test_raw), columns=feature_names, index=test_pos
    )
    X_all = pd.DataFrame(
        imputer.transform(X[feature_names]), columns=feature_names, index=row_positions
    )

    model = fit_model(X_train, y_train, args.random_state)
    test_probability = model.predict_proba(X_test)[:, 1]
    test_prediction = (test_probability >= args.threshold).astype(int)

    metrics = pd.DataFrame(
        [
            {"metric": "ROC-AUC", "value": roc_auc_score(y_test, test_probability)},
            {"metric": "PR-AUC", "value": average_precision_score(y_test, test_probability)},
            {"metric": "Accuracy", "value": accuracy_score(y_test, test_prediction)},
            {
                "metric": "Precision",
                "value": precision_score(y_test, test_prediction, zero_division=0),
            },
            {"metric": "Recall", "value": recall_score(y_test, test_prediction, zero_division=0)},
            {"metric": "F1", "value": f1_score(y_test, test_prediction, zero_division=0)},
            {"metric": "Brier score", "value": brier_score_loss(y_test, test_probability)},
        ]
    )
    metrics.to_excel(args.out_dir / "model_metrics.xlsx", index=False)

    split = np.full(len(df), "train", dtype=object)
    split[test_pos] = "test"
    probabilities = np.full(len(df), np.nan)
    predictions = np.full(len(df), np.nan)
    train_probability = model.predict_proba(X_train)[:, 1]
    probabilities[train_pos] = train_probability
    probabilities[test_pos] = test_probability
    predictions[train_pos] = (train_probability >= args.threshold).astype(int)
    predictions[test_pos] = test_prediction

    id_columns = [
        c for c in ("global_pair_id", "pair_id", "Source_SMILES", "Target_SMILES")
        if c in df.columns
    ]
    prediction_table = df[["input_row_index", *id_columns]].copy()
    prediction_table["data_subset"] = split
    prediction_table["observed_label"] = y.to_numpy()
    prediction_table["observed_risk"] = np.where(y.to_numpy() == 1, "High-risk", "Low-risk")
    prediction_table["probability_high_risk"] = probabilities
    prediction_table["predicted_label"] = predictions.astype(int)
    prediction_table["predicted_risk"] = np.where(
        predictions == 1, "High-risk", "Low-risk"
    )
    prediction_table.to_excel(args.out_dir / "model_predictions.xlsx", index=False)
    save_performance_plots(y_test, test_probability, args.out_dir)

    # Persist the complete fitted preprocessing/model objects and the exact
    # feature list so predictions can be reproduced.
    joblib.dump(
        {
            "imputer": imputer,
            "model": model,
            "feature_names": feature_names,
            "threshold": args.threshold,
            "label_encoding": {"Low-risk": 0, "High-risk": 1},
        },
        args.out_dir / "xgboost_risk_model.joblib",
    )
    feature_audit.to_excel(args.out_dir / "feature_selection_audit.xlsx", index=False)
    pd.DataFrame(
        {"feature": feature_names, "training_median": imputer.statistics_}
    ).to_excel(args.out_dir / "training_imputation_medians.xlsx", index=False)

    # SHAP interpretation is performed on all imputed labeled rows by default,
    # as stated in the revised Methods. --shap_scope test is available for a
    # strictly held-out interpretation sensitivity analysis.
    if args.shap_scope == "all":
        X_shap = X_all.copy()
        shap_row_positions = row_positions
    else:
        X_shap = X_test.copy()
        shap_row_positions = test_pos

    explainer = shap.TreeExplainer(model)
    shap_values = normalize_shap_output(explainer.shap_values(X_shap))
    if shap_values.shape != X_shap.shape:
        raise ValueError(
            f"SHAP matrix shape {shap_values.shape} does not match X {X_shap.shape}."
        )
    mean_abs_shap = pd.Series(
        np.abs(shap_values).mean(axis=0), index=feature_names, name="mean_abs_shap"
    ).sort_values(ascending=False)
    total_importance = float(mean_abs_shap.sum())
    shap_summary = mean_abs_shap.rename_axis("feature").reset_index()
    shap_summary["pct_of_total_shap"] = (
        shap_summary["mean_abs_shap"] / total_importance * 100
        if total_importance > 0
        else 0.0
    )
    shap_summary["rank"] = np.arange(1, len(shap_summary) + 1)
    shap_summary.to_excel(args.out_dir / "shap_global_importance.xlsx", index=False)

    top_features = mean_abs_shap.head(min(args.top_n, len(mean_abs_shap))).index.tolist()
    save_global_shap_plots(
        X_shap, shap_values, shap_summary, args.out_dir, args.top_n
    )

    # Interaction calculations can be large (rows x features x features). Cap
    # the row count to an estimated 512 MB while retaining reproducibility.
    n_features = len(feature_names)
    memory_limited_n = max(1, int(512_000_000 / max(8 * n_features * n_features, 1)))
    interaction_n = min(args.interaction_sample_size, len(X_shap), memory_limited_n)
    if interaction_n < args.interaction_sample_size:
        warnings.warn(
            f"Interaction sample reduced to {interaction_n} rows to limit memory use."
        )
    interaction_sample = X_shap.sample(n=interaction_n, random_state=args.random_state)
    interaction_values = normalize_interaction_output(
        explainer.shap_interaction_values(interaction_sample)
    )
    top_idx = [feature_names.index(feature) for feature in top_features]
    top_interactions = interaction_values[:, top_idx, :][:, :, top_idx]
    interaction_matrix = pd.DataFrame(
        np.abs(top_interactions).mean(axis=0),
        index=top_features,
        columns=top_features,
    )
    interaction_matrix.to_excel(args.out_dir / "top_shap_interaction_matrix.xlsx")

    fig, ax = plt.subplots(figsize=(max(6.0, 0.7 * len(top_features)), max(5.2, 0.62 * len(top_features))))
    image = ax.imshow(interaction_matrix.to_numpy(), cmap="magma", aspect="auto")
    ax.set_xticks(range(len(top_features)), labels=top_features, rotation=45, ha="right")
    ax.set_yticks(range(len(top_features)), labels=top_features)
    ax.set_title("Mean absolute SHAP interaction value")
    fig.colorbar(image, ax=ax, shrink=0.80)
    fig.tight_layout()
    fig.savefig(args.out_dir / "top_shap_interaction_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    strongest_partner = {}
    for feature in top_features:
        scores = interaction_matrix.loc[feature].drop(index=feature, errors="ignore")
        strongest_partner[feature] = scores.idxmax() if not scores.empty else ""

    dependence_frames = []
    for feature in top_features:
        feature_idx = feature_names.index(feature)
        partner = strongest_partner[feature]
        frame = pd.DataFrame(
            {
                "input_row_index": df.loc[shap_row_positions, "input_row_index"].to_numpy(),
                "feature": feature,
                "feature_value": X_shap[feature].to_numpy(),
                "shap_value": shap_values[:, feature_idx],
                "strongest_interaction_feature": partner,
                "interaction_feature_value": (
                    X_shap[partner].to_numpy() if partner else np.nan
                ),
            }
        )
        dependence_frames.append(frame)

        fig, ax = plt.subplots(figsize=(5.4, 4.5))
        if partner:
            scatter = ax.scatter(
                X_shap[feature],
                shap_values[:, feature_idx],
                c=X_shap[partner],
                cmap="coolwarm",
                s=20,
                alpha=0.78,
                edgecolors="none",
            )
            colorbar = fig.colorbar(scatter, ax=ax)
            colorbar.set_label(partner)
        else:
            ax.scatter(X_shap[feature], shap_values[:, feature_idx], s=20, alpha=0.78)
        ax.axhline(0, color="0.6", linewidth=0.8)
        ax.set(xlabel=feature, ylabel="SHAP value", title=f"SHAP dependence: {feature}")
        fig.tight_layout()
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", feature)[:120]
        fig.savefig(args.out_dir / f"shap_dependence_{safe_name}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    pd.concat(dependence_frames, ignore_index=True).to_csv(
        args.out_dir / "top_shap_dependence_long.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # Select the strongest top-feature pairs by absolute SHAP interaction and
    # construct response surfaces with every other feature fixed at its
    # training-set median.
    selected_pairs = strongest_pairs_from_interactions(
        interaction_matrix, args.n_response_pairs
    )
    if not selected_pairs:
        selected_pairs = fallback_pairs(mean_abs_shap, args.n_response_pairs)

    base_row = pd.Series(imputer.statistics_, index=feature_names).to_dict()
    surface_frames = []
    for rank, (feature_1, feature_2, pair_score) in enumerate(selected_pairs, start=1):
        x_low, x_high = X_train[feature_1].quantile([0.02, 0.98])
        y_low, y_high = X_train[feature_2].quantile([0.02, 0.98])
        if x_low == x_high or y_low == y_high:
            continue
        x_grid = np.linspace(x_low, x_high, args.response_grid_size)
        y_grid = np.linspace(y_low, y_high, args.response_grid_size)
        xx, yy = np.meshgrid(x_grid, y_grid)
        grid = pd.DataFrame(
            np.tile([base_row[f] for f in feature_names], (xx.size, 1)),
            columns=feature_names,
        )
        grid[feature_1] = xx.ravel()
        grid[feature_2] = yy.ravel()
        zz = model.predict_proba(grid)[:, 1].reshape(xx.shape)
        surface_frames.append(
            pd.DataFrame(
                {
                    "pair_rank": rank,
                    "feature_1": feature_1,
                    "feature_2": feature_2,
                    "pair_interaction_score": pair_score,
                    "x_value": xx.ravel(),
                    "y_value": yy.ravel(),
                    "predicted_probability_high_risk": zz.ravel(),
                }
            )
        )

        fig, ax = plt.subplots(figsize=(5.6, 4.8))
        contour = ax.contourf(xx, yy, zz, levels=20, cmap="viridis", vmin=0, vmax=1)
        fig.colorbar(contour, ax=ax, label="Predicted probability of High-risk")
        ax.set(xlabel=feature_1, ylabel=feature_2, title=f"Prediction-response surface {rank}")
        fig.tight_layout()
        fig.savefig(args.out_dir / f"response_surface_{rank:02d}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    if surface_frames:
        pd.concat(surface_frames, ignore_index=True).to_csv(
            args.out_dir / "prediction_response_surfaces_long.csv",
            index=False,
            encoding="utf-8-sig",
        )

    run_metadata = {
        "input_file": str(args.input_file.resolve()),
        "n_input_rows": int(len(df_original)),
        "n_labeled_rows": int(len(df)),
        "n_unknown_or_invalid_excluded": n_excluded,
        "class_counts": {str(k): int(v) for k, v in class_counts.items()},
        "n_features": len(feature_names),
        "train_rows": int(len(train_pos)),
        "test_rows": int(len(test_pos)),
        "test_size": args.test_size,
        "random_state": args.random_state,
        "probability_threshold": args.threshold,
        "shap_scope": args.shap_scope,
        "xgboost_parameters": {
            "n_estimators": 400,
            "max_depth": 5,
            "learning_rate": 0.05,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
        },
    }
    with (args.out_dir / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(run_metadata, handle, indent=2, ensure_ascii=False)

    print("Analysis completed.")
    print(metrics.to_string(index=False))
    print(f"Output directory: {args.out_dir.resolve()}")

if __name__ == "__main__":
    main()
