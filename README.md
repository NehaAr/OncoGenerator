# OncoGenerator


[![License: MIT](https://img.shields.io/badge/code-MIT-yellow.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-blue.svg)](LICENSE-DATA.txt)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)

**A probabilistic, fuzzy rule-based generator of synthetic clinical and proteomic data for endometrial cancer research, plus a versioned benchmark data release with ground truth.**

> **All data are synthetic.** They contain no real patient information and are for method development, benchmarking, teaching and pipeline testing. They are not evidence about real endometrial cancer biology, and performance on them does not predict performance on real cohorts.

## What you get

| | |
|---|---|
| **Generator** (`oncogenerator/`) | Seedable Python package + command-line tool + Gradio app (`code/app.py`). |
| **Data release v1.1.0** (full: 80 MB, publish on Zenodo; this repo ships the ~11 MB lite bundle in `release/v1.1.0-lite/`) | 3 core cohorts (1k, 10k, 100k patients) and 6 scenario variants of the 10k cohort: noise low/high, class imbalance, batch effects, MCAR and MNAR missing values. 180 proteins each. ~80 MB as `csv.gz`. |
| **Ground truth** | Per-protein generating parameters, per-contrast true effects (with true nulls), noise-free expected values, batch effects actually applied. See [docs/data_dictionary.md](docs/data_dictionary.md). |
| **Benchmarks** | Runnable scripts in `examples/` and results in [benchmarks/REPRODUCED.md](benchmarks/REPRODUCED.md). |
| **Dataset card** | [docs/datacard.md](docs/datacard.md) - contents, intended use, limitations, licensing. |

## Quick start

```bash
git clone https://github.com/NehaAr/OncoGenerator.git
cd OncoGenerator
pip install numpy pandas scipy scikit-learn          # core dependencies
python -m oncogenerator generate --n 1000 --seed 1 --out out/            # csv.gz
python -m oncogenerator generate --n 10000 --seed 2 --scenario batch5 --out out_batch/
python -m oncogenerator list-scenarios
```

Change `--seed` to get independent replicates; the same seed and options always reproduce the same files. For Parquet (recommended for cloud use) install `pyarrow` and add `--format parquet`. To rebuild the whole release: `python scripts/build_release.py --format parquet`; a ~10 MB lite bundle (scenarios at 2,000 patients, no 100k cohort): `python scripts/build_release.py --lite`.

Python API:

```python
import oncogenerator as sp
ds = sp.generate_dataset(1000, seed=1, scenario="clean")
ds.clinical, ds.proteins          # patient tables
ds.protein_truth, ds.contrast_truth, ds.expected   # ground truth
```

Web app (no code needed): `pip install -r requirements.txt && python code/app.py`.

## Optional: Human Protein Atlas input

`data/HPA.tsv` (Human Protein Atlas export; CC BY 4.0, Uhlen et al. Science 2015) supplies each gene's normal endometrial expression level and TCGA-UCEC prognostic p-value. The generator uses it automatically if present; without it every gene outside the curated dictionary is neutral. See [NOTICE.md](NOTICE.md).

## How the simulation works

See [docs/algorithm.md](docs/algorithm.md) for the full description. In short:

1. **Clinical cohort.** Age ~ Normal(60,15) clipped to 30-85, BMI ~ Normal(25,4) clipped to 21-40; menopausal status is deterministic on age (post-menopausal above 51); a high-risk profile (age >= 60, nulliparous, BMI >= 25) raises the Type 2 probability; tumour type shapes grade, grade shapes stage, stage and grade shape survival.
2. **Protein abundance.** Each protein has a literature-derived regulation direction, an optional clinical factor (e.g. only in Grade 3 or post-menopausal patients), and HPA-derived modifiers. These select one of six Gaussian samplers (strong/weak up/down, low, neutral). A Mamdani-style fuzzy factor kappa in [0.5,1.5] scales effects by grade and stage; proteins in the same pathway share a latent factor (correlation `rho` = 0.4); values are clipped to [-3, 3].

## Known limitations

- Because stage and grade directly scale protein effects, stage is easy to predict from the proteins in this data by construction. High classifier accuracy here says nothing about real data.
- `Myometrial_mm`, `Treatment` and `Ethnicity` are independent noise variables; `Survival_Outcome` is generated from stage and grade and must not be used to predict stage.
- Protein pathway correlation is a single shared latent factor per pathway, not a realistic covariance structure.


## Citation

See [CITATION.cff](CITATION.cff). If you use the data release, please also cite the Human Protein Atlas (Uhlen et al., Science 2015) and, once minted, the dataset DOI.

## License

Code: MIT ([LICENSE](LICENSE)). Data release: CC BY 4.0 ([LICENSE-DATA.txt](LICENSE-DATA.txt)). Third-party notices: [NOTICE.md](NOTICE.md).

**Author:** Neha Arora - PhD candidate, Bioinformatics, University of Newcastle, NSW, Australia.
