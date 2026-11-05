# Reproduced benchmarks (generator v1.1.0)

All numbers below were produced by the scripts in `examples/` on the code in this repository and can be re-run with the commands shown. They are **not** the numbers in the manuscript (`benchmark.md`), which were produced with v1.0.0 code and scripts that were not in the repository. Differences are expected: v1.1.0 changed the generator (see CHANGELOG.md).

**Read these results as properties of the simulator, not of real endometrial cancer.** Stage and grade scale protein effects by construction, so stage is easy to predict from the proteins.

## 1. Early vs late stage classification (`examples/01_stage_classifier.py`)

Random Forest (100 trees), stratified 5-fold CV, mean ± SD over folds.

`python examples/01_stage_classifier.py --n 1000 --seed 1`

```
n=1000  late-stage fraction=0.324  (majority-class accuracy 0.676)
                                                               feature set      accuracy       ROC-AUC
                                                    A. clinical, leak-free 0.732 ± 0.015 0.715 ± 0.021
                                 B. clinical + all 180 proteins, leak-free 0.914 ± 0.015 0.979 ± 0.004
C. v1.0.0-style: 11 clinical (incl. Survival_Outcome) + 8 proteins [LEAKY] 0.771 ± 0.019 0.816 ± 0.016
                     D. same as C but WITHOUT Survival_Outcome [leak-free] 0.762 ± 0.017 0.798 ± 0.024
```

`python examples/01_stage_classifier.py --n 10000 --seed 2`

```
n=10000  late-stage fraction=0.316  (majority-class accuracy 0.684)
                                                               feature set      accuracy       ROC-AUC
                                                    A. clinical, leak-free 0.698 ± 0.010 0.674 ± 0.013
                                 B. clinical + all 180 proteins, leak-free 0.960 ± 0.002 0.995 ± 0.001
C. v1.0.0-style: 11 clinical (incl. Survival_Outcome) + 8 proteins [LEAKY] 0.774 ± 0.009 0.810 ± 0.010
                     D. same as C but WITHOUT Survival_Outcome [leak-free] 0.771 ± 0.007 0.797 ± 0.007
```

- Set C is what the v1.0.0 app used: it includes `Survival_Outcome`, which the generator creates *from* stage and grade (target leakage). Removing it (set D) lowers accuracy by about 1 point here, a small effect only because survival depends weakly on stage in this model. The leakage is still wrong in principle.
- Clinical-only features (A) barely beat the majority-class baseline; adding all 180 proteins (B) gives near-perfect AUC because the protein effects are built from stage and grade.
- The manuscript's Use Case 1 (82% accuracy, AUC 0.88) cannot be compared directly: its exact feature set is not in the repository, and its importance ranking lists FIGO sub-stage, which would be leakage for a stage classifier.

## 2. Differential-expression benchmark with ground truth (`examples/02_de_benchmark.py`)

Late vs early stage; truth from `ground_truth_contrasts` (154 proteins with a generated effect, 26 true nulls). Welch t-test and Mann-Whitney U with Benjamini-Hochberg FDR 5%; mean over 20 random designs. Recall is shown overall and by true effect size (small < 0.3, medium 0.3-0.5, large >= 0.5 log2FC). Observed FDR rests on only 26 true nulls, so it is a rough estimate.

`python examples/02_de_benchmark.py --reps 20 --per-group 100`

```
100 early vs 100 late, 20 random designs, BH FDR 5%
truth: 154 proteins with a generated effect, 26 true nulls (FDR estimates rest on few nulls)

      scenario         test recall small medium large observed_FDR AUROC
         clean      Welch t   0.88  0.84   0.95  0.98         0.01  0.96
         clean Mann-Whitney   0.90  0.86   0.99  0.99         0.01  0.97
     noise_low      Welch t   0.95  0.92   0.97  1.00         0.01  0.98
     noise_low Mann-Whitney   0.96  0.94   1.00  1.00         0.01  0.99
    noise_high      Welch t   0.72  0.66   0.82  0.90         0.01  0.90
    noise_high Mann-Whitney   0.71  0.63   0.88  0.91         0.01  0.91
        batch5      Welch t   0.85  0.80   0.90  0.98         0.01  0.95
        batch5 Mann-Whitney   0.86  0.82   0.97  0.99         0.01  0.96
missing_mcar10      Welch t   0.87  0.83   0.94  0.97         0.01  0.96
missing_mcar10 Mann-Whitney   0.89  0.85   1.00  0.99         0.01  0.97
missing_mnar20      Welch t   0.82  0.79   0.99  0.89         0.01  0.93
missing_mnar20 Mann-Whitney   0.83  0.80   0.99  0.88         0.01  0.93
```

`python examples/02_de_benchmark.py --reps 20 --per-group 25 --scenarios clean,noise_high`

```
25 early vs 25 late, 20 random designs, BH FDR 5%
truth: 154 proteins with a generated effect, 26 true nulls (FDR estimates rest on few nulls)

  scenario         test recall small medium large observed_FDR AUROC
     clean      Welch t   0.67  0.59   0.64  0.87         0.01  0.89
     clean Mann-Whitney   0.65  0.58   0.72  0.85         0.01  0.90
noise_high      Welch t   0.17  0.11   0.12  0.32         0.00  0.75
noise_high Mann-Whitney   0.13  0.08   0.11  0.27         0.00  0.74
```

- Higher noise and MNAR missingness lower sensitivity; MCAR missingness and batch effects (not corrected here) lower it slightly; FDR stays controlled in all scenarios tested.
- limma and DESeq2 (the tools in the manuscript's Use Case 2) are R packages and were **not** run. Export `proteins_log2fc` and `ground_truth_contrasts` to evaluate them. The manuscript's 87% (limma) and 52% (DESeq2) are therefore unverified here.

## 3. Augmentation (`examples/03_augmentation.py`)

The script is provided and was smoke-tested only, on a stand-in file made from synthetic data (which is **not** a result). The manuscript's real TCGA-UCEC n=47 cohort is not in this repository, so Use Case 3 (61% -> 74%, p = 0.004) is **not reproduced**. Run it on your own real data with matching protein column names:

`python examples/03_augmentation.py --real real.csv --label Late_Stage --n-synth 500`
