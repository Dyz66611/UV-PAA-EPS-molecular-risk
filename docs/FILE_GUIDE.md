# File guide

Paths are relative to the package root. Field definitions and exact column names are provided in `docs/DATA_DICTIONARY.md`.

## Data files

| Location | Sheet | Contents |
|---|---|---|
| `data/Smiles_Mapping_Data/` | `ranked_summary` | Selected formula-to-SMILES candidate summaries and no-hit records for each transformation interval. |
| `data/RDKit_Calculation_Data/` | `rdkit_descriptors` | Valid mapped structures and their RDKit descriptors for each transformation interval. |
| `data/Supplementary_DFT_Data/Supplementary_DFT.xlsx` | `Optimized_Coordinates` | Optimized Cartesian coordinates in angstroms. This is the only sheet in the DFT workbook; no energy table is included. |

## User-provided inputs (not included)

- Raw FT-ICR-MS time-point workbooks: the former `data/Raw_FT_ICR_MS_Data/` directory was removed on 2026-10-12. Historical field definitions remain in `docs/DATA_DICTIONARY.md`.
- Model-input table for `code/xgboost_risk_classification.py --input_file`: provide source/target numerical descriptors and binary target-risk labels as described in the README. The former `data/ECOSAR_Toxicity_Valid_Data/ECOSAR_Toxicity_Valid_Data.xlsx` workbook was removed on 2026-10-12.
- Molecular-formula difference matrix and target difference list for `code/formula_difference_matching.py --diff_matrix ... --target_file ...`.
- ECOSAR Batch Output workbook for `code/ecosar_toxicity_labeling.py --input_file`.

## File relationships

- Historical FT-ICR-MS filenames identify the sampling time point; these raw workbooks are not bundled.
- SMILES-mapping and RDKit filenames identify the transformation interval.
- Within the same interval, `original_order` links a mapping record to its corresponding RDKit record.
- XGBoost/SHAP analysis uses the user-provided model-input table and its stored risk labels. The bundled mapping and RDKit workbooks are not a complete labeled model-input table.

## Code and supporting files

| File | Purpose |
|---|---|
| `code/ecosar_toxicity_labeling.py` | Generates molecule-level risk labels from an ECOSAR Batch Output workbook. |
| `code/formula_difference_matching.py` | Matches precursor-product molecular-formula differences to a user-supplied target difference list. The required difference matrix and target list are not included. |
| `code/validate_model_input.py` | Legacy validator with a hard-coded path to the removed model-input workbook. To use it, provide your workbook at that path with a `model_input` sheet, or edit its `INPUT` constant. It has no `--input_file` option. |
| `code/xgboost_risk_classification.py` | Runs model fitting, prediction, preprocessing audits, and SHAP interpretation. |
| `config/model_configuration.json` | Records the model settings; the modeling script does not load this file automatically. |
| `requirements.txt`, `environment.yml` | Define the Python dependencies. |
| `docs/DATA_DICTIONARY.md` | Defines the workbook fields and exact column names. |
| `FILE_MANIFEST.tsv` | Lists packaged files with byte sizes and SHA-256 hashes, excluding the manifest itself. |

Run the analysis commands from the package root as described in `README.md`.
