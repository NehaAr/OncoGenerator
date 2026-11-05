# OncoGenerator — Benchmarks

> Note: OncoGenerator was formerly named SynthProteomics. The manuscript text quoted below still uses the old name.

This document reports the validation and benchmarking results described in the
accompanying manuscript ("SynthProteomics: A Probabilistic Fuzzy Rule-Based
Framework for Synthetic Clinical and Proteomics Data Simulation in Endometrial
Cancer Research"). It is intended as a reference for anyone evaluating whether
SynthProteomics' synthetic data is fit for their own use case.

> **Status:** the figures below are the numbers reported in the manuscript and have **not** been re-run
> against v1.1.0. Runnable scripts are in `examples/` and their current results are in
> [REPRODUCED.md](REPRODUCED.md) (stage classifier, DE benchmark with ground truth, augmentation script;
> the n=47 TCGA-UCEC cohort used for Use Case 3 is not included). Note that the Use Case 1 importance
> ranking below lists FIGO sub-stage as a feature; a stage classifier must not use stage-derived or
> `Survival_Outcome` features (target leakage). See REPRODUCED.md.

---

## Table of contents

- [Reproducibility protocol](#reproducibility-protocol)
- [Use Case 1 — FIGO stage classification (Random Forest)](#use-case-1--figo-stage-classification-random-forest)
- [Use Case 2 — Differential expression tool benchmarking](#use-case-2--differential-expression-tool-benchmarking)
- [Use Case 3 — Small-cohort data augmentation](#use-case-3--small-cohort-data-augmentation)
- [Biological plausibility validation](#biological-plausibility-validation)
- [Cohort size selection rationale](#cohort-size-selection-rationale)
- [A note on the n=500 augmentation citation](#a-note-on-the-n500-augmentation-citation)

---

## Reproducibility protocol

All experiments below share the following methodology, as stated in the manuscript:

| Aspect | Detail |
|---|---|
| Random seed | Fixed (`random_state=42`), used throughout patient sampling, protein-abundance sampling, and model training/splitting |
| Data leakage | Train/test partitioning performed **before** any feature computation; no held-out patient information used at generation or preprocessing time |
| Missing data | None by construction — all data is synthetically generated |
| Categorical encoding | One-hot encoded before model input |
| Continuous features | Used unscaled (Random Forest is scale-invariant) |
| Feature set | No additional feature selection applied — 14 clinical variables plus the requested proteomic panel, small relative to n=1,000 |

---

## Use Case 1 — FIGO stage classification (Random Forest)

**Task:** predict early stage (Stage 1–2) vs. late stage (Stage 3–4) endometrial
cancer from synthetic clinical + proteomic features.

**Setup:** n=1,000 synthetic patients; Random Forest, 100 estimators
(scikit-learn); 80/20 held-out split.

### Held-out test performance (single split)

| Metric | Value |
|---|---|
| Accuracy (synthetic test set) | **82%** |
| Accuracy (real TCGA-UCEC cohort, n=47) | **74%** |
| Precision | 0.80 |
| Recall | 0.83 |
| F1-score | 0.81 |
| Specificity | 0.82 |
| ROC-AUC | 0.88 |

### Stratified 5-fold cross-validation (n=1,000, stratified on FIGO-stage class)

| Metric | Mean ± SD |
|---|---|
| Accuracy | 0.81 ± 0.02 |
| Precision | 0.79 ± 0.03 |
| Recall | 0.80 ± 0.03 |
| F1-score | 0.80 ± 0.02 |
| Specificity | 0.83 ± 0.03 |
| ROC-AUC | 0.87 ± 0.02 |

The cross-validated results are consistent with the single 80/20-split accuracy
(82%), indicating the reported performance is not an artefact of a favorable split.

### Hyperparameter selection

Grid search, 5-fold CV on the training partition only, optimizing mean F1-score:

| Hyperparameter | Grid searched | Selected |
|---|---|---|
| `n_estimators` | {100, 200, 500} | 100 |
| `max_depth` | {5, 10, 20, None} | None (unrestricted) |
| `min_samples_leaf` | {1, 2, 5} | 1 |

**Why Random Forest** (over logistic regression, gradient boosting, SVM):
- Natively models nonlinear interactions among mixed continuous (age, BMI, log2FC)
  and categorical (histology, stage) variables
- Comparatively robust to label/feature noise intrinsic to probabilistically
  generated data
- Built-in Gini feature-importance estimates, used as a plausibility check (below)

Gradient-boosted trees (XGBoost/LightGBM) were tested in preliminary experiments
and gave comparable but not superior held-out accuracy (±1–2 percentage points);
noted in the repository as an available alternative.

### Feature importance (biological plausibility check)

Gini importance ranked, in order: **FIGO sub-stage > tumor grade > p53/PTEN/HER2
log2FC values**. This ordering matches known clinical/molecular correlates of
endometrial cancer staging, offered as evidence the classifier is using
biologically plausible signal rather than simulation artefacts.

---

## Use Case 2 — Differential expression tool benchmarking

**Task:** recover a known, seeded set of differentially expressed proteins from
synthetic proteomic cohorts, using standard DE tools.

**Setup:** n=200 synthetic patients (100 Stage I vs. 100 Stage III), with
ground-truth log2FC values assigned at generation time.

| Tool | % of seeded DE proteins recovered (FDR < 0.05) | Notes |
|---|---|---|
| **limma** | **87%** | Appropriate for continuous log2FC/intensity data |
| limma-trend | 84–87% | Recommended comparator for continuous data |
| Moderated t-test (normalized intensities) | 84–87% | Recommended comparator for continuous data |
| Wilcoxon rank-sum | — (benchmarked; see manuscript) | Non-parametric comparator |
| **DESeq2** | **52%** | Retained only as a **negative control** — assumes discrete, overdispersed count data under a negative-binomial model, which is a statistical mismatch for SynthProteomics' continuous log2FC output |

**Interpretation:** the gap between limma (87%) and DESeq2 (52%) is presented as
expected, given the underlying statistical mismatch between DESeq2's count-based
model and SynthProteomics' continuous intensity output — not as a finding that
DESeq2 is unreliable in general.

---

## Use Case 3 — Small-cohort data augmentation

**Task:** test whether adding synthetic patients to a small real cohort improves
downstream classifier generalization.

**Setup:** real TCGA-UCEC cohort (n=47) augmented with 500 synthetic patients
(n=547 total), compared against the real-only cohort (n=47).

| Configuration | Accuracy | CV variance |
|---|---|---|
| Real data only (n=47) | 61% | ±0.12 |
| Real + synthetic (n=547) | **74%** | **±0.06** |

**Statistical significance:** two-sided paired t-test (fold-level) on the
accuracy improvement: **p = 0.004**. 95% bootstrap confidence interval for the
improvement: **[6.1, 19.8] percentage points**.

**Interpretation:** the accuracy gain and the reduction in cross-validation
variance are both presented as evidence that synthetic augmentation can help
offset the small-sample-size constraints typical of uncommon cancer subtypes.

---

## Biological plausibility validation

Independent of the three downstream use cases, SynthProteomics' output is
checked against real reference distributions:

| Check | Reference | Result |
|---|---|---|
| Synthetic clinical variable distributions (age, BMI, stage) | AIHW 2023 cancer registry statistics | Kolmogorov–Smirnov test, **p > 0.05** for all three (i.e. no significant difference from the real registry distribution) |
| Synthetic proteomic expression ranges | HPA expression quintiles, matched endometrial tissue | Benchmarked qualitatively against quintile ranges |

---

## Cohort size selection rationale

The manuscript states each use case's cohort size was chosen for a specific,
stated reason rather than tuned to maximize the reported result:

| Use case | n | Stated rationale |
|---|---|---|
| 1 (classifier training) | 1,000 | Smallest cohort size at which held-out RF accuracy plateaued, in a preliminary sweep over n ∈ {200, 500, 1000, 2000} |
| 2 (DE benchmarking) | 200 (100 vs. 100) | Matches common design sizes used in DE benchmarking studies, to keep per-group variance comparable across the three tools tested |
| 3 (augmentation) | 500 synthetic added | Chosen to roughly match the order of magnitude of the real TCGA-UCEC cohort, augmented by one order of magnitude |

---


