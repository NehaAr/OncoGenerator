# -*- coding: utf-8 -*-
"""
OncoGenerator web app (Gradio front-end for the `oncogenerator` package).

Run locally:   pip install -r requirements.txt && python code/app.py
Run in Colab:  clone the repo, `pip install -r requirements.txt`, then `python code/app.py`.
No hard-coded /content or /tmp paths; downloads use the system temp directory.

Author: Neha Arora, University of Newcastle.  MIT License.
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import gradio as gr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import oncogenerator as sp

CLIN_COLS = ["Patient_ID", "Ages", "Ethnicity", "Menopause", "Menopause_label", "Grade", "Tumor_type",
             "Histological_Subtype", "Stage", "Nulliparity", "BMI", "Myometrial_mm", "Treatment",
             "Survival_Outcome", "Late_Stage"]
DEFAULT_PANEL = "ANXA2, PKM2, ERBB2, EGFR, MMP9, CASP3, HMGB3, PRDX1"
STATE = {"clinical": pd.DataFrame(), "proteins": pd.DataFrame(), "seed": 42}


def _tmp(name):
    return str(Path(tempfile.gettempdir()) / name)


def gen_clinical(n, seed, subtype, stage, cols):
    STATE["seed"] = int(seed)
    try:
        df = sp.generate_clinical(int(n), seed=int(seed), subtype=subtype, stage_prefix=stage)
    except ValueError as e:
        raise gr.Error(str(e))
    STATE["clinical"], STATE["proteins"] = df, pd.DataFrame()
    cols = [c for c in cols if c in df.columns]
    return df[cols] if cols else df


def gen_proteins(n, text):
    clin = STATE["clinical"]
    if clin.empty:
        raise gr.Error("Generate clinical data first (Clinical Data tab).")
    panel = [p.strip() for p in text.replace("\n", ",").split(",") if p.strip()]
    res = sp.generate_proteins(clin.iloc[: int(n)].reset_index(drop=True), panel, seed=STATE["seed"])
    STATE["proteins"] = res.values
    return res.values.set_index("Patient_ID").round(3).to_dict(orient="index")


def dl_clinical():
    if STATE["clinical"].empty:
        return None
    p = _tmp("OncoGenerator_clinical.csv"); STATE["clinical"].to_csv(p, index=False); return p


def dl_proteins():
    if STATE["proteins"].empty:
        return None
    p = _tmp("OncoGenerator_protein.csv"); STATE["proteins"].to_csv(p, index=False); return p


def run_rf():
    """Early (I/II) vs late (III/IV) classifier on leak-free clinical features only.
    Survival_Outcome and Stage are deliberately NOT used (Survival_Outcome is generated from stage)."""
    clin = STATE["clinical"]
    if clin.empty or len(clin) < 50:
        return "Generate at least 50 patients in the Clinical Data tab first."
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import accuracy_score, classification_report
    from sklearn.model_selection import train_test_split
    feats = ["Ages", "BMI", "Menopause", "Nulliparity", "Tumor_type", "Grade"]
    X = pd.get_dummies(clin[feats]); y = clin["Late_Stage"]
    if y.nunique() < 2:
        return "Only one stage class in this cohort; clear the Stage filter."
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
    clf = RandomForestClassifier(n_estimators=100, random_state=42).fit(Xtr, ytr)
    pred = clf.predict(Xte)
    return (f"Random Forest accuracy: {accuracy_score(yte, pred) * 100:.1f}%\n"
            f"(features: {', '.join(feats)}; majority-class baseline "
            f"{max(y.mean(), 1 - y.mean()) * 100:.1f}%)\n\n"
            + classification_report(yte, pred, target_names=["Early (I/II)", "Late (III/IV)"]))


def plots():
    d = STATE["clinical"]
    fig, ax = plt.subplots(figsize=(6, 3))
    if d.empty:
        ax.text(0.5, 0.5, "Generate clinical data first", ha="center", va="center"); ax.axis("off")
        return fig
    plt.close(fig)
    cols = [c for c in ["Ages", "BMI", "Grade", "Stage", "Tumor_type", "Menopause", "Nulliparity",
                        "Treatment", "Survival_Outcome"] if c in d.columns]
    rows = (len(cols) + 1) // 2
    fig, axes = plt.subplots(rows, 2, figsize=(14, rows * 3), squeeze=False)
    for i, c in enumerate(cols):
        ax = axes[i // 2][i % 2]
        if pd.api.types.is_numeric_dtype(d[c]) and d[c].nunique() > 8:
            ax.hist(d[c], bins=15, color="steelblue", edgecolor="white")
        else:
            vc = d[c].value_counts(); ax.bar(vc.index.astype(str), vc.values, color="steelblue")
            ax.tick_params(axis="x", rotation=35)
        ax.set_title(c, fontweight="bold", fontsize=10)
    for i in range(len(cols), rows * 2):
        axes[i // 2][i % 2].set_visible(False)
    plt.tight_layout()
    return fig


with gr.Blocks(title="OncoGenerator") as demo:
    gr.Markdown("# OncoGenerator\nSynthetic clinical & proteomic data for endometrial cancer "
                "method development. **All data are synthetic.** MIT licence.")
    with gr.Tab("Clinical Data"):
        with gr.Row():
            n_c = gr.Number(label="Number of patients", value=100, minimum=10)
            seed = gr.Number(label="Random seed", value=42, precision=0)
            sub = gr.Dropdown(["All", "Endometrioid", "Serous", "Clear Cell", "Mucinous", "Undifferentiated"],
                              value="All", label="Subtype filter")
            stg = gr.Dropdown(["All", "Stage1", "Stage2", "Stage3", "Stage4"], value="All", label="Stage filter")
        cols = gr.CheckboxGroup(CLIN_COLS, value=["Patient_ID", "Ages", "Grade", "Stage", "Tumor_type",
                                                  "Histological_Subtype", "Survival_Outcome"],
                                label="Columns to display")
        b = gr.Button("Generate Clinical Data", variant="primary"); out_c = gr.Dataframe()
        b.click(gen_clinical, [n_c, seed, sub, stg, cols], out_c)
        f = gr.File(label="Clinical CSV"); gr.Button("Download as CSV").click(dl_clinical, [], f)
    with gr.Tab("Protein Abundance Data"):
        gr.Markdown("_Generate clinical data first. Values are log2 fold-changes clipped to [-3, 3]._")
        n_p = gr.Number(label="Patients (first N of the cohort)", value=10, minimum=1)
        txt = gr.Textbox(lines=3, label="Protein list (comma-separated)", value=DEFAULT_PANEL)
        out_p = gr.JSON(); gr.Button("Generate Abundance Data", variant="primary").click(gen_proteins, [n_p, txt], out_p)
        fp = gr.File(label="Protein CSV"); gr.Button("Download as CSV").click(dl_proteins, [], fp)
    with gr.Tab("Distribution Plots"):
        po = gr.Plot(); gr.Button("Generate Plots", variant="primary").click(plots, [], po)
    with gr.Tab("Use Case: Stage Classifier"):
        gr.Markdown("Random Forest, early vs late stage, leak-free clinical features only.")
        ro = gr.Textbox(lines=10); gr.Button("Run Random Forest", variant="primary").click(run_rf, [], ro)

if __name__ == "__main__":
    demo.launch()
