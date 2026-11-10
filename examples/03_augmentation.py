#!/usr/bin/env python
"""Does adding synthetic patients to a small real cohort help a classifier?

Needs YOUR real data: a CSV with a 0/1 label column and numeric feature columns
whose names match proteins in the synthetic panel. Real-only vs real+synthetic is compared with repeated
stratified CV in which synthetic rows are added to TRAIN folds only.
    python examples/03_augmentation.py --real real.csv --label Late_Stage --n-synth 500
"""
import argparse
from _common import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--real", required=True)
    ap.add_argument("--label", default="Late_Stage"); ap.add_argument("--n-synth", type=int, default=500)
    ap.add_argument("--seed", type=int, default=7); ap.add_argument("--repeats", type=int, default=10)
    a = ap.parse_args()
    real = pd.read_csv(a.real)
    clin, prot, _ = load_or_generate(None, a.n_synth, a.seed)
    feats = [c for c in real.columns if c != a.label and c in prot.columns]
    if not feats: raise SystemExit("No feature columns of the real file match the synthetic protein names.")
    Xr, yr = real[feats].to_numpy(float), real[a.label].to_numpy(int)
    Xs, ys = prot[feats].to_numpy(float), clin["Late_Stage"].to_numpy(int)
    print(f"real n={len(real)} features={len(feats)}  synthetic n={len(Xs)}")
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=a.repeats, random_state=a.seed)
    res = {"real only": [], "real + synthetic": []}
    for tr, te in cv.split(Xr, yr):
        for k, (X, y) in {"real only": (Xr[tr], yr[tr]),
                          "real + synthetic": (np.vstack([Xr[tr], Xs]), np.r_[yr[tr], ys])}.items():
            m = RandomForestClassifier(100, random_state=42, n_jobs=-1).fit(X, y)
            res[k].append((m.score(Xr[te], yr[te]), roc_auc_score(yr[te], m.predict_proba(Xr[te])[:, 1])))
    for k, v in res.items():
        v = np.array(v); print(f"{k:18s} accuracy {v[:,0].mean():.3f} ± {v[:,0].std():.3f}   AUC {v[:,1].mean():.3f} ± {v[:,1].std():.3f}")

if __name__ == "__main__":
    main()
