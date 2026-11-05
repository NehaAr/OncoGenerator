# Changelog

## 1.1.0 - 2026-11-05
**Renamed:** the project is now **OncoGenerator** (formerly SynthProteomics; package `synthproteomics` -> `oncogenerator`, repository `Synthetic-omics-data` -> `OncoGenerator`). The v1.0.0 entry below used the old name.
Behaviour changes (generated values differ from 1.0.0 even with the same seed):
- **Seed is a parameter** (`--seed`, `seed=`). v1.0.0 re-seeded with 42 on every call, so all cohorts were identical. Random streams now come from `numpy.random.default_rng`; distributions of the clinical model are unchanged (checked against the v1.0.0 generator on 20,000 patients: marginals agree within sampling error).
- **HPA prognostic parser fixed.** v1.0.0 stripped every non-digit character from strings like `unprognostic (3.21e-1)`, giving 3.211, so no gene was ever "prognostic" and the strong-effect samplers (+/-2.0) never fired. Now the p-value in parentheses is parsed; 3,954 HPA genes are flagged (p <= 0.05); 36 of the 180 default proteins switch to a strong sampler.
- **Clinical factor matching fixed.** In v1.0.0 factors `POSTMENOPAUSAL`, `TYPE 1/2` and `stage1/2/3` (prefix) never matched (menopause was an integer, tumour type was not compared, stage compared exactly). 21 of the 180 default proteins are affected.
- **Consistent Patient_ID** (`EC0000...`) in clinical and protein tables (v1.0.0 used `EC0000` and `Patient0000`).
- **Dictionary alignment made explicit.** The Gene/Regulation/Factor lists had lengths 230/212/203; positions beyond the shortest list raised IndexError for 14 genes (e.g. CAPS, DES, CDH13). Rows are now in `oncogenerator/data/curated_regulation_v1.csv` with `alignment_status`; the default panel uses the 180 unique genes with complete entries; others are treated as uncurated. Genes with conflicting duplicate entries use the first entry (listed in the ground truth).
- **Stage classifier in the app no longer uses `Survival_Outcome`** (it is generated from stage). The tutorial already described the leak-free feature set.
- Added `Late_Stage`, `Menopause_label`; `Batch` in batch scenarios.
New:
- `oncogenerator` package, CLI, Parquet/CSV writer, tests, data release v1.1.0, ground-truth tables, scenarios (noise, imbalance, batch, MCAR, MNAR), benchmark scripts, dataset card, data dictionary.
- Local use: no `/content` or `/tmp` paths; HPA path is relative.
Docs: README/algorithm corrected (menopause threshold is > 51; pathway correlation exists; genes in HPA but not in the dictionary are not neutral).
`--legacy` reproduces the v1.0.0 factor and HPA rules (not its random numbers).

## 1.0.0
Initial release as SynthProteomics (Gradio app).
