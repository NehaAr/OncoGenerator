# OncoGenerator - Algorithm (v1.1.0)

OncoGenerator generates two linked tables: (1) clinical records and (2) protein log2 fold-change values. This page matches the code in `oncogenerator/core.py`.

## Inputs & outputs

| Symbol | Description |
|---|---|
| `N`, `seed` | Number of patients; random seed (all randomness comes from `numpy.random.default_rng`) |
| `P` | Protein list (default: 180 curated genes) |
| `D` | Curated regulation table (`oncogenerator/data/curated_regulation_v1.csv`) |
| `HPA` | Optional Human Protein Atlas export (`data/HPA.tsv`) |
| `rho` | Pathway correlation, default 0.4 |
| `C`, `A` | Outputs: clinical table; protein log2FC matrix |

## Part 1 - Clinical data

```
for each patient i = 1..N
  Age_i       ~ Normal(60, 15), truncated to integer, clipped to [30, 85]
  BMI_i       ~ Normal(25, 4), clipped to [21, 40]
  Ethnicity_i ~ Categorical(European .55, East Asian .15, South Asian .10, Middle Eastern .05,
                            African .03, Indigenous Australian .03, Other .09)
  Menopause_i = 1 (pre/peri) if Age_i <= 51 else 2 (post)
  Nulliparity_i = yes with prob 0.5 if Age_i <= 35 else 0.6
  HighRisk_i  = Age_i >= 60 AND Nulliparity_i = yes AND BMI_i >= 25
  Tumor_type_i ~ Type2 with prob 0.6 if HighRisk_i else 0.2 (otherwise Type1)
  Grade_i     ~ Type1: (G1 .5, G2 .3, G3 .2);  Type2: (G1 .2, G2 .3, G3 .5)
  Subtype_i   ~ Type1: Endometrioid .9 / Mucinous .1;  Type2: Serous .55 / Clear Cell .35 / Undifferentiated .10
  Stage_i     ~ ten stage labels; Grade 1-2 weights 0.1333 on Stage1..Stage2B and 0.05 on Stage3C1..Stage4B,
                Grade 3 the reverse (weights are normalised)
  Survival_i  ~ Alive with prob .90 (Stage1* and Grade1), .30 (Stage4*), else .65
  Treatment_i ~ uniform over 4 labels;  Myometrial_mm_i ~ Uniform(2, 40)   (independent noise columns)
  Late_Stage_i = 1 if Stage is III/IV (derived label)
```

The stage weights are described in earlier versions as calibrated to AIHW 2023 registry data; that calibration is reported in the accompanying manuscript and was not re-run for v1.1.0.

## Part 2 - Protein abundance

```
for each protein p: look up regulation r_p and factor f_p in D (genes not in D: r_p = UP/DOWN, f_p = NA);
                    look up normal-tissue level e_p and prognostic flag s_p in HPA (if absent: moderate, no)
  sampler_p = 
    r_p = UP:   up_strong N(+2.0, .4) if s_p = yes else up_weak N(+1.0, .3)
    r_p = DOWN: down_strong N(-2.0, .4) if s_p = yes else down_weak N(-1.0, .3)
    r_p = UP/DOWN: follows e_p: up -> as UP, down -> as DOWN, moderate -> neutral
    r_p = LOW:  N(-0.5, .3) clipped to [-1.5, 0]
    otherwise:  neutral N(0, .5) clipped to [-1.5, 1.5]
for each patient j
  kappa_j   = Mamdani-style fuzzy scale in [0.5, 1.5] from grade and stage (higher = larger effects)
  z_jk      ~ Normal(0,1) per pathway k
  for each protein p
     base   = sampler_p draw if the clinical factor f_p applies to patient j, else neutral draw
     if p belongs to pathway k:  base += rho * z_jk + sqrt(1 - rho^2) * eps,   eps ~ Normal(0,1)
     A[j,p] = clip(kappa_j * base, -3, 3)
```

`noise_scale` (scenarios `noise_low`, `noise_high`) multiplies every random term. The **expected** (noise-free) value is `kappa_j * sampler mean` when the factor applies, else 0.

### Clinical factor matching (v1.1.0)
A protein's factor is a list such as `stage3`, `GRADE2`, `POSTMENOPAUSAL`, `TYPE 2`. Stage terms match by prefix (`stage3` matches Stage3C1 and Stage3C2), grade terms match exactly, `POSTMENOPAUSAL` means Menopause = 2, `TYPE n` matches the tumour type. A protein is affected if **any** listed term matches. (In v1.0.0 menopause, tumour-type and plain `stage3` terms never matched.)

### HPA parsing (v1.1.0)
`prognostic = yes` if the p-value in parentheses in the HPA TCGA-UCEC column is <= 0.05. (In v1.0.0 the string was stripped of every non-digit character, so no gene was ever flagged.) Normal expression: > 100 nTPM up, 10-100 moderate, < 10 down.

## Scenarios

| Name | Change |
|---|---|
| `noise_low` / `noise_high` | all random terms x0.5 / x2.0 (same signal) |
| `imbalanced_late10` | exactly 10% late-stage patients |
| `batch5` | 5 batches; per-batch additive shift per protein ~ N(0, 0.3) and per-batch gain ~ N(1, 0.1), applied after clipping |
| `missing_mcar10` | 10% of values missing completely at random |
| `missing_mnar20` | ~20% missing; probability rises as the value falls |

## Ground truth
Per protein: sampler and parameters. Per contrast (late vs early stage, grade 3 vs 1, post vs pre-menopausal, type 2 vs 1): mean expected log2FC difference, a non-zero flag (true null if zero) and a large-effect flag (|effect| >= 0.5). Effects are before noise and clipping.

## References
- Zadeh, L.A. (1965). Fuzzy sets. *Information and Control*, 8(3), 338-353.
- Uhlen, M. et al. (2015). Tissue-based map of the human proteome. *Science*, 347, 1260419.
- Australian Institute of Health and Welfare (2023). *Cancer in Australia 2023*.
