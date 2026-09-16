# File guide

Paths are relative to the package root. Field definitions and exact column names are provided in `docs/DATA_DICTIONARY.md`.

## Data files

| Location | Sheet | Contents |
|---|---|---|
| `data/Raw_FT_ICR_MS_Data/` | `1` or `Sheet1` | FT-ICR-MS molecular-formula and peak records for the time point identified by each filename. |
| `data/Smiles_Mapping_Data/` | `ranked_summary` | Selected formula-to-SMILES candidate summaries and no-hit records for each transformation interval. |
| `data/RDKit_Calculation_Data/` | `rdkit_descriptors` | Valid mapped structures and their RDKit descriptors for each transformation interval. |
| `data/ECOSAR_Toxicity_Valid_Data/ECOSAR_Toxicity_Valid_Data.xlsx` | `model_input` | Final model-input table containing source and target descriptors, descriptor changes, stage fields, and stored binary target-risk labels. |
| `data/Supplementary_DFT_Data/Supplementary_DFT.xlsx` | `Optimized_Coordinates` | Optimized Cartesian coordinates in angstroms. This is the only sheet in the DFT workbook; no energy table is included. |

## File relationships

- FT-ICR-MS filenames identify the sampling time point.
- SMILES-mapping and RDKit filenames identify the transformation interval.
- Within the same interval, `original_order` links a mapping record to its corresponding RDKit record.
- The ECOSAR toxicity workbook is the supplied input for the XGBoost/SHAP analysis. The model uses its stored text and numeric risk labels.

## Code and supporting files

| File | Purpose |
|---|---|
| `code/ecosar_toxicity_labeling.py` | Generates molecule-level risk labels from an ECOSAR Batch Output workbook. |
| `code/formula_difference_matching.py` | Matches precursor-product molecular-formula differences to a user-supplied target difference list. The required difference matrix and target list are not included. |
| `code/validate_model_input.py` | Checks the supplied model-input workbook for required identifiers, labels, label agreement, and descriptor fields. |
| `code/xgboost_risk_classification.py` | Runs model fitting, prediction, preprocessing audits, and SHAP interpretation. |
| `config/model_configuration.json` | Records the model settings; the modeling script does not load this file automatically. |
| `requirements.txt`, `environment.yml` | Define the Python dependencies. |
| `docs/DATA_DICTIONARY.md` | Defines the workbook fields and exact column names. |
| `FILE_MANIFEST.tsv` | Lists packaged files with byte sizes and SHA-256 hashes, excluding the manifest itself. |

Run the analysis commands from the package root as described in `README.md`.
