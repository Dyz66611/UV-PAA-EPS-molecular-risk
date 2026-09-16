```python
"""
Match molecular-formula differences between precursor and product formulas.

This script is provided for reproducibility of the molecular transformation
analysis. The input difference matrix and target Formula differ list are not
distributed in this repository because of file size and data-management
constraints. Users should provide their own input files.

Workflow:
1. Read a molecular-formula difference matrix.
2. Convert the matrix into precursor-product molecular pairs.
3. Match the Formula differ values against a target difference list.
4. Count the occurrence of each matched molecular difference.
5. Export matched molecular pairs and difference statistics.

Example
-------
python code/formula_difference_matching.py ^
    --diff_matrix path/to/difference_matrix.xlsx ^
    --target_file path/to/reaction_type_classification.xlsx ^
    --output results/formula_difference_matching.xlsx
"""

from pathlib import Path
import argparse

import pandas as pd


def load_target_differences(
    target_file: Path,
    column_name: str = "Formula differ",
) -> set[str]:
    """Load target molecular-formula differences from an Excel file."""

    if not target_file.exists():
        raise FileNotFoundError(
            f"Target difference file not found: {target_file}"
        )

    target_df = pd.read_excel(target_file)

    if column_name not in target_df.columns:
        raise ValueError(
            f"Column '{column_name}' was not found in {target_file}. "
            f"Available columns: {list(target_df.columns)}"
        )

    return set(
        target_df[column_name]
        .dropna()
        .astype(str)
        .str.strip()
    )


def flatten_difference_matrix(
    diff_matrix_file: Path,
) -> pd.DataFrame:
    """
    Convert a molecular-formula difference matrix into a long-format table.

    The first column is treated as the precursor molecular formula index.
    Columns represent product molecular formulas.
    """

    if not diff_matrix_file.exists():
        raise FileNotFoundError(
            f"Difference matrix not found: {diff_matrix_file}"
        )

    diff_matrix = pd.read_excel(
        diff_matrix_file,
        index_col=0,
    )

    records = []

    for before_formula, row in diff_matrix.iterrows():
        for after_formula, diff_value in row.items():

            # Preserve the original workflow:
            # only non-empty string difference values are retained.
            if isinstance(diff_value, str) and diff_value.strip():

                records.append(
                    {
                        "sumFormula-before": str(before_formula),
                        "sumFormula-after": str(after_formula),
                        "Formula differ": diff_value.strip(),
                    }
                )

    return pd.DataFrame(
        records,
        columns=[
            "sumFormula-before",
            "sumFormula-after",
            "Formula differ",
        ],
    )


def match_formula_differences(
    difference_df: pd.DataFrame,
    target_differences: set[str],
) -> pd.DataFrame:
    """Select molecular pairs whose Formula differ is in the target list."""

    return difference_df[
        difference_df["Formula differ"].isin(target_differences)
    ].copy()


def summarize_matches(
    matched_df: pd.DataFrame,
) -> pd.DataFrame:
    """Count the occurrence of each matched molecular-formula difference."""

    if matched_df.empty:
        return pd.DataFrame(
            columns=["Formula differ", "Count"]
        )

    return (
        matched_df["Formula differ"]
        .value_counts()
        .rename_axis("Formula differ")
        .reset_index(name="Count")
    )


def process_difference_matrix(
    diff_matrix_file: Path,
    target_differences: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Process one difference matrix and return matched pairs and statistics."""

    difference_df = flatten_difference_matrix(diff_matrix_file)

    matched_df = match_formula_differences(
        difference_df,
        target_differences,
    )

    count_df = summarize_matches(matched_df)

    return count_df, matched_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Match molecular-formula differences between "
            "precursor and product molecules."
        )
    )

    parser.add_argument(
        "--diff_matrix",
        required=True,
        type=Path,
        help="Path to the molecular-formula difference matrix (.xlsx).",
    )

    parser.add_argument(
        "--target_file",
        required=True,
        type=Path,
        help=(
            "Path to the Excel file containing the target "
            "'Formula differ' list."
        ),
    )

    parser.add_argument(
        "--target_column",
        default="Formula differ",
        help="Column containing target molecular-formula differences.",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output Excel file.",
    )

    args = parser.parse_args()

    print("Loading target molecular-formula differences...")
    target_differences = load_target_differences(
        args.target_file,
        args.target_column,
    )

    print(
        f"Loaded {len(target_differences)} unique target "
        "molecular-formula differences."
    )

    print("Processing molecular-formula difference matrix...")
    count_df, matched_df = process_difference_matrix(
        args.diff_matrix,
        target_differences,
    )

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with pd.ExcelWriter(args.output) as writer:

        count_df.to_excel(
            writer,
            index=False,
            sheet_name="Difference_Statistics",
        )

        matched_df.to_excel(
            writer,
            index=False,
            sheet_name="Matched_Molecular_Pairs",
        )

    print("Analysis completed.")
    print(f"Matched molecular pairs: {len(matched_df)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
```
