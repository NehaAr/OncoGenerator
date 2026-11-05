#!/usr/bin/env python
"""Early (I/II) vs late (III/IV) stage classification, 5-fold stratified CV.

Shows how much the reported accuracy depends on which features are allowed.
`Survival_Outcome` is generated FROM stage and grade, so using it to predict
stage is target leakage; the v1.0.0 app included it. Run:
    python examples/01_stage_classifier.py --n 1000 --seed 1
"""
import argparse
from _common import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate

CLEAN = ["Ages", "BMI", "Menopause", "Nulliparity", "Tumor_type", "Grade", "Ethnicity"]
V100_CLIN = ["Ages", "Ethnicity", "Menopause", "Grade", "Tumor_type", "Histological_Subtype",
             "Nulliparity", "BMI", "Myometrial_mm", "Treatment", "Survival_Outcome"]
V100_PROT = ["CTNNB", "EGFR", "ERBB2", "PKM2", "MMP9", "CASP3", "HMGB3", "PRDX1"]  # PTEN/TP53 are not in the 180-gene panel

def run(X, y, seed):
    cv = StratifiedKFold(5, shuffle=True, random_state=seed)
    r = cross_validate(RandomForestClassifier(100, random_state=42, n_jobs=-1), X, y, cv=cv,
                       scoring=["accuracy", "roc_auc"])
    return r["test_accuracy"], r["test_roc_auc"]

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--release")
    a = ap.parse_args()
    clin, prot, _ = load_or_generate(a.release, a.n, a.seed)
    y = clin["Late_Stage"].to_numpy()
    P = prot.drop(columns="Patient_ID")
    sets = {
        "A. clinical, leak-free": clin[CLEAN],
        "B. clinical + all 180 proteins, leak-free": pd.concat([clin[CLEAN], P], axis=1),
        "C. v1.0.0-style: 11 clinical (incl. Survival_Outcome) + 8 proteins [LEAKY]":
            pd.concat([clin[V100_CLIN], P[V100_PROT]], axis=1),
        "D. same as C but WITHOUT Survival_Outcome [leak-free]":
            pd.concat([clin[[c for c in V100_CLIN if c != "Survival_Outcome"]], P[V100_PROT]], axis=1),
    }
    print(f"n={len(clin)}  late-stage fraction={y.mean():.3f}  (majority-class accuracy {max(y.mean(),1-y.mean()):.3f})")
    rows = []
    for name, X in sets.items():
        acc, auc = run(pd.get_dummies(X, drop_first=False), y, a.seed)
        rows.append((name, f"{acc.mean():.3f} ± {acc.std():.3f}", f"{auc.mean():.3f} ± {auc.std():.3f}"))
    print(pd.DataFrame(rows, columns=["feature set", "accuracy", "ROC-AUC"]).to_string(index=False))

if __name__ == "__main__":
    main()
