# Data dictionary

This dictionary describes the fields actually present in the packaged workbooks. Exact column lists are provided below; original spelling and capitalization are preserved.

## FT-ICR-MS time-point data

Location: `data/Raw_FT_ICR_MS_Data/`. Each row is a formula/peak record at the time indicated by the filename.

| Fields | Description |
|---|---|
| `sumFormula`; `C`, `H`, `N`, `O`, `S` | Recorded molecular formula and elemental counts. |
| `ObservedIntens`, `relativeintens`, `sumintens` | Recorded intensity, relative intensity, and intensity-total fields. |
| `ObservedM_z`, `calc_M_z` | Observed and calculated mass-to-charge values. |
| `errMDa`, `errPpm`, `abserrPpm` | Mass error in mDa, mass error in ppm, and absolute ppm error. |
| `mSigma` | Stored peak-assignment score field. |
| `H1`, `HCraito`, `OCraito`, `CN2`, `dbe_o`, `NSP` | Original elemental/ratio/derived fields, retaining the source export names and values. |
| `avemass`, `aimod`, `dbe`, `nosc` | Stored average-mass, modified aromaticity, double-bond-equivalent, and nominal carbon oxidation-state fields. |
| `kmd*`, `nkm*` | Stored Kendrick mass-defect and nominal Kendrick-mass fields for the suffix-designated reference units. |

## Formula-to-SMILES mapping

Location: `data/Smiles_Mapping_Data/`, sheet `ranked_summary`. One selected candidate summary or no-hit record per formula record. The interval comes from the filename.

| Fields | Description |
|---|---|
| `original_order` | Record-order key; use together with the file interval to link to the RDKit table. |
| `sumFormula`, `formula` | Source formula text and compact formula representation. |
| `top1_cid`, `cid_order` | Selected PubChem candidate identifier and stored candidate retrieval-order field. |
| `top1_smiles`, `top1_inchikey`, `returned_formula` | Selected candidate structure identifiers and returned formula. |
| `rdkit_valid`, `has_smiles` | Stored structure-validity and SMILES-presence flags. |
| `formal_charge`, `has_metal` | Selected candidate formal charge and metal-presence flag. |
| `exact_formula_match`, `neutral_preferred`, `nonmetal_preferred` | Stored formula-match and candidate-preference flags. |
| `rank_score`, `candidate_rank` | Candidate score and rank retained in the summary. |
| `candidate_count_total`, `valid_candidate_count` | Total and valid candidate counts. |
| `matched`, `confidence_level` | Match status (`Yes`/`No`) and stored confidence category (`High`, `Medium`, `No hit`). |

No-hit records retain their original flag/sentinel values, including `cid_order = 999999`, and blank candidate-identification fields. Confidence categories describe the stored candidate assignments.

## RDKit descriptors

Location: `data/RDKit_Calculation_Data/`, sheet `rdkit_descriptors`. The first 21 columns reproduce the corresponding mapping fields. Twenty additional columns contain:

| Fields | Description |
|---|---|
| `canonical_smiles_rdkit` | RDKit canonical structure representation. |
| `MolWt`, `ExactMolWt` | Molecular-weight and exact-molecular-weight descriptors. |
| `MolLogP` | Calculated octanol/water partition descriptor. |
| `TPSA` | Topological polar surface area (square angstroms). |
| `HeavyAtomCount` | Number of non-hydrogen atoms. |
| `NumHAcceptors`, `NumHDonors` | Hydrogen-bond acceptor and donor counts. |
| `NumRotatableBonds` | Rotatable-bond count. |
| `RingCount`, `NumAromaticRings` | Total ring and aromatic-ring counts. |
| `FractionCSP3` | Fraction of carbon atoms with sp3 hybridization. |
| `AtomCount_C`, `AtomCount_N`, `AtomCount_O`, `AtomCount_S`, `AtomCount_P` | Structure-derived elemental counts. |
| `HeteroAtomCount`, `FormalCharge`, `BertzCT` | Heteroatom count, total formal charge, and Bertz complexity. |

`formal_charge` is retained mapping metadata; `FormalCharge` is the RDKit descriptor column. Their case-sensitive names are preserved.

## Final valid data and risk labels

Location: `data/ECOSAR_Toxicity_Valid_Data/ECOSAR_Toxicity_Valid_Data.xlsx`, sheet `model_input`.

| Field/group | Description |
|---|---|
| `global_pair_id` | Unique stage-specific model-record identifier, such as `0-10__1`. |
| `stage` | Transformation interval: `0-10`, `10-30`, or `30-60`. |
| `stage_order` | Interval order, respectively 1, 2, or 3. |
| `Target_risk_label_binary` | Stored product risk label, `Low-risk` or `High-risk`. |
| `y_target_risk_binary` | Numeric encoding: Low-risk = 0, High-risk = 1. |
| `Source_` descriptors | 25 source-molecule elemental, formula-derived, and RDKit descriptors. |
| `Target_` descriptors | 23 target-molecule descriptors; the separately listed target-label column also has this prefix. |
| `delta_` descriptors | 24 stored target-minus-source descriptor-change fields. |
| `stage_0-10`, `stage_10-30`, `stage_30-60` | One-hot stage indicators. |

Elemental suffixes `nC`, `nH`, `nO`, `nN`, `nS` denote counts. `ExactMass_formula` is the formula-derived exact mass. `H_C`, `O_C`, `N_C`, `S_C` denote elemental ratios. `DBE`, `DBE_C`, `AI_mod`, and `NOSC` denote double-bond equivalents, DBE per carbon, modified aromaticity index, and nominal carbon oxidation state. RDKit suffixes follow the definitions above. Consult the exact column list for the descriptors included in each group.

The model uses the stored labels directly. Its default predictors are eligible numerical source and target descriptors; label fields are screened out. Stage and delta columns remain available as accompanying information. Training-set median imputation is performed by the model script when predictor values are blank.

The ECOSAR labeling script implements High-risk for the final representative effect concentration ≤10 mg/L and Low-risk for >10 mg/L. Concentration processing is described in the README; the final valid workbook supplies the risk labels for modeling.

## DFT records

Location: `data/Supplementary_DFT_Data/Supplementary_DFT.xlsx`.

| Sheet/field | Description |
|---|---|
| `Optimized_Coordinates`: `Model class` | Model-family description. |
| `Species ID` | Species identifier within the coordinate dataset. |
| `Atom index`, `Atom` | Atom index within a species and element symbol. |
| `X (Å)`, `Y (Å)`, `Z (Å)` | Cartesian coordinates in angstroms. |

`Optimized_Coordinates` is the only sheet retained in the supplementary DFT workbook; no energy table is included.

## Exact column lists

The following lists follow the column order in the supplied workbooks. Workbooks in the same group have the same column names.

### FT-ICR-MS: all four time points

`mSigma`, `ObservedIntens`, `ObservedM_z`, `calc_M_z`, `errMDa`, `errPpm`, `sumFormula`, `C`, `H`, `N`, `O`, `S`, `H1`, `HCraito`, `OCraito`, `CN2`, `dbe_o`, `NSP`, `abserrPpm`, `sumintens`, `relativeintens`, `avemass`, `aimod`, `dbe`, `nosc`, `kmdch2`, `nkmch2`, `kmdh2`, `nkmh2`, `kmdh2o`, `nkmh2o`, `kmdo`, `nkmo`, `kmdcoo`, `nkmcoo`.

### SMILES mapping: all three intervals

`original_order`, `sumFormula`, `formula`, `top1_cid`, `cid_order`, `top1_smiles`, `returned_formula`, `top1_inchikey`, `rdkit_valid`, `formal_charge`, `has_metal`, `exact_formula_match`, `neutral_preferred`, `nonmetal_preferred`, `has_smiles`, `rank_score`, `candidate_rank`, `candidate_count_total`, `valid_candidate_count`, `matched`, `confidence_level`.

### RDKit: all three intervals

`original_order`, `sumFormula`, `formula`, `top1_cid`, `cid_order`, `top1_smiles`, `returned_formula`, `top1_inchikey`, `rdkit_valid`, `formal_charge`, `has_metal`, `exact_formula_match`, `neutral_preferred`, `nonmetal_preferred`, `has_smiles`, `rank_score`, `candidate_rank`, `candidate_count_total`, `valid_candidate_count`, `matched`, `confidence_level`, `canonical_smiles_rdkit`, `MolWt`, `ExactMolWt`, `MolLogP`, `TPSA`, `HeavyAtomCount`, `NumHAcceptors`, `NumHDonors`, `NumRotatableBonds`, `RingCount`, `NumAromaticRings`, `FractionCSP3`, `AtomCount_C`, `AtomCount_N`, `AtomCount_O`, `AtomCount_S`, `AtomCount_P`, `HeteroAtomCount`, `FormalCharge`, `BertzCT`.

### Final valid data: model_input

`global_pair_id`, `stage`, `Target_risk_label_binary`, `y_target_risk_binary`, `stage_order`, `Source_nC`, `Source_nH`, `Source_nO`, `Source_nN`, `Source_nS`, `Source_ExactMass_formula`, `Source_H_C`, `Source_O_C`, `Source_N_C`, `Source_S_C`, `Source_DBE`, `Source_DBE_C`, `Source_AI_mod`, `Source_NOSC`, `Source_MolLogP`, `Source_TPSA`, `Source_NumHAcceptors`, `Source_NumHDonors`, `Source_NumRotatableBonds`, `Source_RingCount`, `Source_NumAromaticRings`, `Source_FractionCSP3`, `Source_HeteroAtomCount`, `Source_FormalCharge`, `Source_BertzCT`, `Target_nC`, `Target_nH`, `Target_nO`, `Target_nS`, `Target_ExactMass_formula`, `Target_H_C`, `Target_O_C`, `Target_S_C`, `Target_DBE`, `Target_DBE_C`, `Target_AI_mod`, `Target_NOSC`, `Target_MolLogP`, `Target_TPSA`, `Target_NumHAcceptors`, `Target_NumHDonors`, `Target_NumRotatableBonds`, `Target_RingCount`, `Target_NumAromaticRings`, `Target_FractionCSP3`, `Target_HeteroAtomCount`, `Target_FormalCharge`, `Target_BertzCT`, `delta_nC`, `delta_nH`, `delta_nO`, `delta_nN`, `delta_nS`, `delta_ExactMass_formula`, `delta_H_C`, `delta_O_C`, `delta_N_C`, `delta_S_C`, `delta_DBE`, `delta_DBE_C`, `delta_AI_mod`, `delta_NOSC`, `delta_MolLogP`, `delta_TPSA`, `delta_NumHAcceptors`, `delta_NumHDonors`, `delta_NumRotatableBonds`, `delta_RingCount`, `delta_NumAromaticRings`, `delta_FractionCSP3`, `delta_FormalCharge`, `delta_BertzCT`, `stage_0-10`, `stage_10-30`, `stage_30-60`.

### DFT: Optimized_Coordinates

`Model class`, `Species ID`, `Atom index`, `Atom`, `X (Å)`, `Y (Å)`, `Z (Å)`.
