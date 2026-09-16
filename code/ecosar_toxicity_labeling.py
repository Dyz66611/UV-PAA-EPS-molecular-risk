from pathlib import Path
import argparse

import numpy as np
import pandas as pd

RISK_CUTOFF_MG_L = 10.0

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate molecule-level risk labels from ECOSAR output."
    )
    parser.add_argument(
        "--input_file",
        type=Path,
        required=True,
        help="Input ECOSAR Excel workbook.",
    )
    parser.add_argument(
        "--rows_out",
        type=Path,
        required=True,
        help="Output Excel file containing labeled ECOSAR rows.",
    )
    parser.add_argument(
        "--summary_out",
        type=Path,
        required=True,
        help="Output Excel file containing one row per molecule.",
    )
    return parser.parse_args()

def read_excel_auto(path: Path) -> pd.DataFrame:
    """Read ECOSAR's Batch Output sheet, or the first sheet if absent."""
    xls = pd.ExcelFile(path)
    sheet = "Batch Output" if "Batch Output" in xls.sheet_names else xls.sheet_names[0]
    df = pd.read_excel(path, sheet_name=sheet)
    df.columns = [str(column).strip() for column in df.columns]
    print(f"Loaded sheet: {sheet}")
    print(f"Input shape: {df.shape}")
    return df

def find_col(columns, candidates):
    """Locate a column by exact name first and partial name second."""
    columns_lower = {str(column).lower(): column for column in columns}
    for candidate in candidates:
        if candidate.lower() in columns_lower:
            return columns_lower[candidate.lower()]

    for column in columns:
        column_lower = str(column).lower()
        for candidate in candidates:
            if candidate.lower() in column_lower:
                return column
    return None

def clean_text(value):
    if pd.isna(value):
        return ""
    value = str(value).strip()
    return "" if value.lower() in {"nan", "none", "null"} else value

def normalize_term(value):
    """Normalize endpoint and organism terms without changing original columns."""
    return " ".join(clean_text(value).lower().replace("_", " ").split())

def conc_to_log10(concentration):
    if pd.isna(concentration) or concentration <= 0:
        return np.nan
    return np.log10(concentration)

def risk_label_binary(concentration):
    """Assign the binary risk label used in all downstream analyses."""
    if pd.isna(concentration) or concentration <= 0:
        return "Unknown"
    return (
        "Low-risk"
        if concentration > RISK_CUTOFF_MG_L
        else "High-risk"
    )

def min_with_driver(
    subset,
    concentration_col,
    class_col,
    organism_col,
    endpoint_col,
    alert_col,
):
    """Select the lowest positive prediction and retain its provenance."""
    empty_result = {
        "min_conc": np.nan,
        "driver_class": "",
        "driver_org": "",
        "driver_endpoint": "",
        "driver_alert": "",
    }
    if subset.empty:
        return pd.Series(empty_result)

    valid_subset = subset[
        subset[concentration_col].notna() & (subset[concentration_col] > 0)
    ]
    if valid_subset.empty:
        return pd.Series(empty_result)

    row = valid_subset.loc[valid_subset[concentration_col].idxmin()]
    return pd.Series(
        {
            "min_conc": row[concentration_col],
            "driver_class": row[class_col] if class_col else "",
            "driver_org": row[organism_col] if organism_col else "",
            "driver_endpoint": row[endpoint_col] if endpoint_col else "",
            "driver_alert": row[alert_col] if alert_col else "",
        }
    )

def choose_representative(core_result, all_result, core_basis, all_basis):
    """Prefer core freshwater organisms, with all organisms as a fallback."""
    if pd.notna(core_result["min_conc"]):
        return core_result["min_conc"], core_basis, core_result
    return all_result["min_conc"], all_basis, all_result

def main():
    args = parse_args()
    df = read_excel_auto(args.input_file)

    # Detect the standard ECOSAR columns.
    col_number = find_col(df.columns, ["Number"])
    col_smiles = find_col(df.columns, ["SMILES"])
    col_class = find_col(df.columns, ["ECOSAR Class"])
    col_org = find_col(df.columns, ["Organism"])
    col_endpoint = find_col(df.columns, ["End Point", "Endpoint"])
    col_conc = find_col(
        df.columns,
        ["Concentration (mg/L)", "Concentration (mg L-1)", "Concentration"],
    )
    col_alert = find_col(df.columns, ["Alert"])
    col_chemical = find_col(df.columns, ["Chemical"])

    required = {
        "Number": col_number,
        "SMILES": col_smiles,
        "ECOSAR Class": col_class,
        "Organism": col_org,
        "End Point": col_endpoint,
        "Concentration": col_conc,
    }
    missing = [name for name, column in required.items() if column is None]
    if missing:
        raise ValueError(f"Missing required ECOSAR columns: {missing}")

    # Retain only records with a SMILES and a positive predicted effect concentration.
    for column in [col_smiles, col_class, col_org, col_endpoint, col_alert, col_chemical]:
        if column is not None:
            df[column] = df[column].apply(clean_text)

    df[col_conc] = pd.to_numeric(df[col_conc], errors="coerce")
    df[col_number] = pd.to_numeric(df[col_number], errors="coerce")
    df = df[
        df[col_conc].notna()
        & (df[col_conc] > 0)
        & (df[col_smiles].astype(str).str.strip() != "")
    ].copy()
    if df.empty:
        raise ValueError("No valid ECOSAR rows remained after quality filtering.")
    print(f"Shape after quality filtering: {df.shape}")

    # Molecules are identified by the ECOSAR record number-SMILES combination.
    df["chem_key"] = df.apply(
        lambda row: (
            f"{int(row[col_number])}__{row[col_smiles]}"
            if pd.notna(row[col_number])
            else f"NA__{row[col_smiles]}"
        ),
        axis=1,
    )

    endpoint_normalized = df[col_endpoint].apply(normalize_term)
    organism_normalized = df[col_org].apply(normalize_term)

    acute_endpoints = {"lc50", "ec50"}
    chronic_endpoints = {"chv", "chronic value"}
    freshwater_orgs = {
        "fish",
        "daphnid",
        "daphnids",
        "daphnia",
        "green algae",
        "green alga",
    }

    df["endpoint_type"] = np.select(
        [
            endpoint_normalized.isin(acute_endpoints),
            endpoint_normalized.isin(chronic_endpoints),
        ],
        ["acute", "chronic"],
        default="other",
    )
    df["is_freshwater_core"] = organism_normalized.isin(freshwater_orgs)

    # Row-level labels are useful for auditing but are not the final molecule label.
    df["row_log10_conc"] = df[col_conc].apply(conc_to_log10)
    df["row_risk_label_binary"] = df[col_conc].apply(risk_label_binary)

    summary_rows = []
    for chem_key, subset in df.groupby("chem_key", sort=False):
        first_row = subset.iloc[0]

        acute_all = min_with_driver(
            subset[subset["endpoint_type"] == "acute"],
            col_conc,
            col_class,
            col_org,
            col_endpoint,
            col_alert,
        )
        chronic_all = min_with_driver(
            subset[subset["endpoint_type"] == "chronic"],
            col_conc,
            col_class,
            col_org,
            col_endpoint,
            col_alert,
        )
        acute_core = min_with_driver(
            subset[
                (subset["endpoint_type"] == "acute")
                & subset["is_freshwater_core"]
            ],
            col_conc,
            col_class,
            col_org,
            col_endpoint,
            col_alert,
        )
        chronic_core = min_with_driver(
            subset[
                (subset["endpoint_type"] == "chronic")
                & subset["is_freshwater_core"]
            ],
            col_conc,
            col_class,
            col_org,
            col_endpoint,
            col_alert,
        )

        acute_conc, acute_basis, acute_driver = choose_representative(
            acute_core,
            acute_all,
            "acute_min_freshwater_core",
            "acute_min_all_organisms",
        )
        chronic_conc, chronic_basis, chronic_driver = choose_representative(
            chronic_core,
            chronic_all,
            "chronic_min_freshwater_core",
            "chronic_min_all_organisms",
        )

        candidate_concentrations = [
            value
            for value in [acute_conc, chronic_conc]
            if pd.notna(value)
        ]
        if candidate_concentrations:
            final_conc = min(candidate_concentrations)
            if pd.notna(acute_conc) and (
                pd.isna(chronic_conc) or acute_conc <= chronic_conc
            ):
                final_basis = acute_basis
                final_driver = acute_driver
            else:
                final_basis = chronic_basis
                final_driver = chronic_driver
        else:
            final_conc = np.nan
            final_basis = "unknown"
            final_driver = pd.Series(
                {
                    "driver_class": "",
                    "driver_org": "",
                    "driver_endpoint": "",
                    "driver_alert": "",
                }
            )

        # The minimum available acute/chronic value provides a conservative
        # molecule-level estimate. This directly yields the binary target.
        integrated_label_binary = risk_label_binary(final_conc)
        integrated_rule = (
            "minimum_available_acute_or_chronic"
            if pd.notna(final_conc)
            else "both_endpoint_groups_unknown"
        )

        summary_rows.append(
            {
                "chem_key": chem_key,
                "Number": first_row[col_number],
                "Chemical": first_row[col_chemical] if col_chemical else "",
                "SMILES": first_row[col_smiles],
                "n_rows": len(subset),
                "n_ecosar_classes": subset[col_class].nunique(),
                "acute_risk_conc_mgL": acute_conc,
                "acute_risk_log10_mgL": conc_to_log10(acute_conc),
                "acute_risk_basis": acute_basis,
                "acute_driver_class": acute_driver["driver_class"],
                "acute_driver_org": acute_driver["driver_org"],
                "acute_driver_endpoint": acute_driver["driver_endpoint"],
                "acute_driver_alert": acute_driver["driver_alert"],
                "acute_risk_label_binary": risk_label_binary(acute_conc),
                "chronic_risk_conc_mgL": chronic_conc,
                "chronic_risk_log10_mgL": conc_to_log10(chronic_conc),
                "chronic_risk_basis": chronic_basis,
                "chronic_driver_class": chronic_driver["driver_class"],
                "chronic_driver_org": chronic_driver["driver_org"],
                "chronic_driver_endpoint": chronic_driver["driver_endpoint"],
                "chronic_driver_alert": chronic_driver["driver_alert"],
                "chronic_risk_label_binary": risk_label_binary(chronic_conc),
                "final_risk_conc_mgL": final_conc,
                "final_risk_log10_mgL": conc_to_log10(final_conc),
                "final_risk_basis": final_basis,
                "final_driver_class": final_driver["driver_class"],
                "final_driver_org": final_driver["driver_org"],
                "final_driver_endpoint": final_driver["driver_endpoint"],
                "final_driver_alert": final_driver["driver_alert"],
                "final_risk_rule": integrated_rule,
                "risk_label_binary": integrated_label_binary,
            }
        )

    df_summary = pd.DataFrame(summary_rows)

    # Map all molecule-level results back to the original valid ECOSAR rows.
    map_columns = [
        column
        for column in df_summary.columns
        if column not in {"Number", "Chemical", "SMILES", "n_rows", "n_ecosar_classes"}
    ]
    df_rows_labeled = df.merge(
        df_summary[map_columns],
        on="chem_key",
        how="left",
        validate="many_to_one",
    )

    args.rows_out.parent.mkdir(parents=True, exist_ok=True)
    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(args.rows_out, engine="openpyxl") as writer:
        df_rows_labeled.to_excel(writer, sheet_name="rows_labeled", index=False)
    with pd.ExcelWriter(args.summary_out, engine="openpyxl") as writer:
        df_summary.to_excel(writer, sheet_name="summary_risk", index=False)

    print(f"Unique molecules summarized: {len(df_summary)}")
    print("Final binary risk distribution:")
    print(df_summary["risk_label_binary"].value_counts(dropna=False))
    print("Processing completed.")

if __name__ == "__main__":
    main()
