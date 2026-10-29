# UV/PAA Transformation and Molecular Hazard Analysis of EPS

This repository contains formula-to-SMILES mapping summaries, RDKit descriptors, DFT optimized-coordinate data, and Python analysis scripts. Raw FT-ICR-MS time-point workbooks and the final valid model-input workbook with binary risk labels are not included in the current tree. The analysis scripts remain available, but the missing inputs must be provided by the user.

## Package structure

```text
.
├── code/                         ECOSAR labeling, molecular-formula difference matching, model validation, and XGBoost/SHAP scripts
├── config/model_configuration.json
├── data/
│   ├── Smiles_Mapping_Data/       Three formula_smiles_ranked_summary_*min.xlsx files
│   ├── RDKit_Calculation_Data/    Three formula_smiles_rdkit_descriptors_*min.xlsx files
│   └── Supplementary_DFT_Data/Supplementary_DFT.xlsx
├── docs/                          File guide and data dictionary, including historical schemas for removed data
├── environment.yml
├── requirements.txt
└── FILE_MANIFEST.tsv
```


## Data availability and user-provided inputs

On 2026-10-12, `data/ECOSAR_Toxicity_Valid_Data/ECOSAR_Toxicity_Valid_Data.xlsx` was removed in commit [`f21901ee2f387cce0ccc643eb99815dc12824124`](https://github.com/Dyz66611/UV-PAA-EPS-molecular-risk/commit/f21901ee2f387cce0ccc643eb99815dc12824124), and `data/Raw_FT_ICR_MS_Data/` (0, 10, 30, and 60 min workbooks) was removed in commit [`7e4298d9ae3992e91f54051d51a9f32774b12c1d`](https://github.com/Dyz66611/UV-PAA-EPS-molecular-risk/commit/7e4298d9ae3992e91f54051d51a9f32774b12c1d). These paths are historical references, not bundled example inputs.

| Analysis | User-provided input (not included) | CLI arguments |
|---|---|---|
| XGBoost/SHAP | Model-input table containing eligible numerical molecular descriptors and binary target-risk labels | `--input_file` |
| Molecular-formula difference matching | Difference matrix workbook and target difference-list workbook | `--diff_matrix`, `--target_file` |
| ECOSAR labeling | ECOSAR Batch Output workbook | `--input_file` |

The bundled mapping and RDKit workbooks are not a complete labeled source/target model-input table. The scripts do not automatically reconstruct the removed model-input workbook from these files. Raw FT-ICR-MS workbooks are also user-provided if needed for upstream work; neither analysis command below reads that former directory automatically.

## Installation and XGBoost/SHAP analysis

Python 3.10 is the recommended baseline. Dependencies are declared in `requirements.txt` and `environment.yml`. Run from the package root. For Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

After preparing your own model-input table, replace the placeholder path below with its actual location. This command requires user-provided data; it cannot run from the bundled data alone.

```powershell
.venv\Scripts\python.exe code/xgboost_risk_classification.py --input_file "path/to/your/model_input.xlsx" --out_dir "results/reproduced_run"
```

Input requirements follow `code/xgboost_risk_classification.py`:

- Accepted formats: XLSX, XLS, CSV, TSV, or tab-delimited TXT. XLS files additionally require a compatible Excel reader such as `xlrd`, which is not declared in the supplied dependencies. For Excel, the default sheet is `model_input` when present, otherwise the first sheet; use `--sheet_name` to select another sheet.
- Include at least one of `Target_risk_label_binary` (`Low-risk` / `High-risk`) or `y_target_risk_binary` (0 / 1). If both provide valid labels for a row, they must agree. Label-column names can be changed with `--label_text_col` and `--label_num_col`; rows without a valid binary label are excluded.
- Provide eligible numerical predictors under the default `Source_` and `Target_` prefixes, or specify `--feature_prefixes` / `--feature_file`. Label, toxicity, identifier, and other leakage-prone fields are screened out.
- Both classes must remain, with at least two observations per class and enough rows for the chosen stratified train/test split. Predictors must retain usable variation in the training subset.

Defaults are a stratified 80:20 random row split, seed 42, probability threshold 0.5, and SHAP analysis on all labeled rows. `config/model_configuration.json` documents model settings; the script uses its CLI defaults and source-defined estimator parameters rather than loading this JSON.

`code/validate_model_input.py` is a legacy validator with a hard-coded path to the removed workbook and no `--input_file` option. It requires a user-provided workbook at that historical path (sheet `model_input`), or an edited `INPUT` constant. It is not a prerequisite for the command above; the model script performs its own label and feature checks.

## ECOSAR labeling code

`code/ecosar_toxicity_labeling.py` implements the concentration-to-label processing step. It requires a user-provided ECOSAR Batch Output workbook, which is not included. Its required CLI arguments are `--input_file`, `--rows_out`, and `--summary_out`. Its outputs do not by themselves assemble the source/target descriptor table needed for modeling.
The labeling script retains nonempty SMILES strings and positive numeric concentrations. LC50/EC50 are acute endpoints; ChV/chronic value are chronic endpoints. Within each endpoint group, fish, daphnid/daphnia, and green-algae records are preferred, with other organisms used as fallback. The lower available acute/chronic representative concentration is classified as High-risk at ≤10 mg/L and Low-risk at >10 mg/L; unavailable representative concentrations receive Unknown.

## Molecular-formula difference matching code

`code/formula_difference_matching.py` provides the precursor–product molecular-pair matching workflow based on molecular-formula differences. The script converts a molecular-formula difference matrix into precursor–product pair records, matches the resulting `Formula differ` values against a user-provided target difference list, and generates both matched molecular-pair records and difference-frequency statistics.

Both input workbooks must be provided by the user. The first sheet of the difference matrix must have precursor formulas in its first column and product formulas as the remaining column headers; nonempty string cells hold formula differences. The first sheet of the target workbook must contain a `Formula differ` column, or another column selected with `--target_column`. Matching uses exact strings after trimming surrounding whitespace.

PowerShell example, after replacing both placeholder paths with actual user-provided workbooks:

```powershell
.venv\Scripts\python.exe code/formula_difference_matching.py `
    --diff_matrix "path/to/your/difference_matrix.xlsx" `
    --target_file "path/to/your/target_difference_list.xlsx" `
    --output "results/formula_difference_matching.xlsx"
```

The output workbook contains two sheets:

- `Difference_Statistics`: occurrence counts for each matched molecular-formula difference.
- `Matched_Molecular_Pairs`: precursor–product molecular pairs and their corresponding formula differences.

## Outputs generated by the model script

A successful run creates outputs under the chosen `--out_dir`, including:

- `model_metrics.xlsx`: test-set performance metrics.
- `model_predictions.xlsx`: predictions for both training and test records, identified by `data_subset`.
- `xgboost_risk_model.joblib`: fitted imputer, classifier, feature names, and threshold.
- `feature_selection_audit.xlsx` and `training_imputation_medians.xlsx`: predictor selection and preprocessing records.
- SHAP importance, interaction, dependence and response-surface outputs, plots, and `run_metadata.json`.

The output directory is created when the analysis runs.

## Integrity and citation

`FILE_MANIFEST.tsv` records the current version-controlled file paths, byte sizes, and SHA-256 hashes, excluding the manifest itself. Python bytecode caches are excluded; use the `.py` scripts for analysis.

The README paths and CLI arguments were checked against the current repository tree and script source on 2026-10-29. This is a documentation and file-integrity check; it does not establish reproduction of the original analysis with the missing inputs. User-provided inputs and generated results are not part of the tracked-file manifest. Historical schemas in `docs/DATA_DICTIONARY.md` do not imply that the removed workbooks are bundled.

Please cite the associated manuscript when using these data or scripts.
