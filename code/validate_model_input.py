from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT
    / "data"
    / "ECOSAR_Toxicity_Valid_Data"
    / "ECOSAR_Toxicity_Valid_Data.xlsx"
)


def main() -> None:
    df = pd.read_excel(INPUT, sheet_name="model_input")
    required = {
        "global_pair_id",
        "stage",
        "Target_risk_label_binary",
        "y_target_risk_binary",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if df["global_pair_id"].duplicated().any():
        raise ValueError("global_pair_id contains duplicates")
    text_to_num = df["Target_risk_label_binary"].map(
        {"Low-risk": 0, "High-risk": 1}
    )
    numeric = pd.to_numeric(df["y_target_risk_binary"], errors="coerce")
    if text_to_num.isna().any() or numeric.isna().any() or not numeric.isin([0, 1]).all():
        raise ValueError("Labels contain missing or invalid values")
    if not text_to_num.eq(numeric).all():
        raise ValueError("Text and numeric labels disagree")
    prefixes = ("Source_", "Target_")
    features = [c for c in df.columns if c.startswith(prefixes)]
    if not features:
        raise ValueError("No Source_/Target_ features found")
    print(f"Rows: {len(df)}")
    print(f"Features before script screening: {len(features)}")
    print(df["Target_risk_label_binary"].value_counts().to_string())
    print("Model input validation passed.")


if __name__ == "__main__":
    main()
