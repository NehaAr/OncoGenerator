#!/usr/bin/env python
"""Differential-expression benchmark with GROUND TRUTH (late vs early stage).

A protein is a true positive if the generator gave it a non-zero expected log2FC
difference (is_nonzero_effect in ground_truth_contrasts, contrast Late_vs_Early_stage);
the others are true nulls. Recall is also reported by effect size (small <0.3,
medium 0.3-0.5, large >=0.5). For each scenario we draw
`--reps` random designs of `--per-group` early vs late patients from a 10,000
patient cohort and test every protein with Welch t and Mann-Whitney U (BH, FDR 5%).
    python examples/02_de_benchmark.py --reps 20 --per-group 100
limma/DESeq2 are R tools and are not run here; export proteins_log2fc.csv.gz and
ground_truth_contrasts.csv.gz to run them yourself.
"""
import argparse
from _common import *
from scipy import stats
from sklearn.metrics import roc_auc_score

def evaluate(prot, truth, early, late, fdr=0.05):
    cols = [c for c in prot.columns if c != "Patient_ID"]
    A, B = prot.loc[early, cols].to_numpy(float), prot.loc[late, cols].to_numpy(float)
    pt = stats.ttest_ind(B, A, equal_var=False, nan_policy="omit").pvalue
    pw = np.array([stats.mannwhitneyu(B[:, j][~np.isnan(B[:, j])], A[:, j][~np.isnan(A[:, j])]).pvalue
                   for j in range(len(cols))])
    tr = truth.set_index("protein").loc[cols]
    t = tr["is_nonzero_effect"].to_numpy(bool)
    eff = tr["true_effect_log2fc"].abs().to_numpy()
    bins = {"recall_small": t & (eff < 0.3), "recall_medium": t & (eff >= 0.3) & (eff < 0.5),
            "recall_large": eff >= 0.5}
    out = {}
    for name, p in (("Welch t", pt), ("Mann-Whitney", pw)):
        p = np.nan_to_num(p, nan=1.0); hit = bh_fdr(p) < fdr
        tp = (hit & t).sum()
        out[name] = dict(recall=tp / max(t.sum(), 1), fdp=(hit & ~t).sum() / max(hit.sum(), 1),
                         auroc=roc_auc_score(t, -np.log10(np.maximum(p, 1e-300))),
                         **{k: (hit & m).sum() / max(m.sum(), 1) for k, m in bins.items()})
    return out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--per-group", type=int, default=100); ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=2)
    ap.add_argument("--scenarios", default="clean,noise_low,noise_high,batch5,missing_mcar10,missing_mnar20")
    a = ap.parse_args(); rng = np.random.default_rng(a.seed)
    rows = []
    for sc in a.scenarios.split(","):
        clin, prot, ct = load_or_generate(None, a.n, a.seed, sc)
        truth = ct[ct.contrast == "Late_vs_Early_stage"]
        res = {"Welch t": [], "Mann-Whitney": []}
        for _ in range(a.reps):
            e = rng.choice(np.where(clin.Late_Stage == 0)[0], a.per_group, replace=False)
            l = rng.choice(np.where(clin.Late_Stage == 1)[0], a.per_group, replace=False)
            r = evaluate(prot, truth, e, l)
            for k in res: res[k].append(r[k])
        for k, v in res.items():
            d = pd.DataFrame(v)
            rows.append(dict(scenario=sc, test=k, recall=f"{d.recall.mean():.2f}",
                             small=f"{d.recall_small.mean():.2f}", medium=f"{d.recall_medium.mean():.2f}",
                             large=f"{d.recall_large.mean():.2f}", observed_FDR=f"{d.fdp.mean():.2f}",
                             AUROC=f"{d.auroc.mean():.2f}"))
    print(f"{a.per_group} early vs {a.per_group} late, {a.reps} random designs, BH FDR 5%")
    print(f"truth: {int(truth.is_nonzero_effect.sum())} proteins with a generated effect, "
          f"{int((~truth.is_nonzero_effect).sum())} true nulls (FDR estimates rest on few nulls)\n")
    print(pd.DataFrame(rows).to_string(index=False))

if __name__ == "__main__":
    main()
