# OncoGenerator — Step-by-Step Tutorial

This tutorial walks you through a complete, worked example: generating a synthetic
endometrial cancer patient cohort, simulating protein abundance data for that cohort,
visualizing the results, and running the built-in machine learning demo — end to end,
with no prior coding experience required.

By the end of this tutorial you will have:
- Generated a synthetic clinical cohort of 50 patients
- Generated log2 fold-change protein abundance data for 5 marker proteins
- Downloaded both as CSV files
- Visualized the clinical variable distributions
- Trained and evaluated a Random Forest stage classifier on your synthetic cohort

Estimated time: **10 minutes**.

---

## Before you start

You need a Google account and a web browser. No installation is required for this
tutorial — we'll run everything in Google Colab.

---

## Step 1 - Launch the app

Colab or any machine:

```bash
git clone https://github.com/NehaAr/OncoGenerator.git
cd OncoGenerator
pip install -r requirements.txt
python code/app.py
```

Open the local address printed in the output (or the public `gradio.live` link in Colab). The Human Protein Atlas file `data/HPA.tsv` is found automatically from the repository folder. No path editing is needed.

You should see a page titled **OncoGenerator** with four tabs: *Clinical Data*, *Protein Abundance Data*, *Distribution Plots*, and *Use Case: Stage Classifier*. Prefer the command line? `python -m oncogenerator generate --n 50 --seed 1 --out my_cohort/`.

---

## Step 2 — Generate a synthetic patient cohort

We'll start on the **Clinical Data** tab, which should already be open.

1. Set **Number of Patients** to `50` and leave **Random seed** at `42` (change it for an independent cohort).
2. Leave **Cancer Subtype Filter** as `All`.
3. Leave **Stage Filter** as `All`.
4. In **Select Columns to Display**, make sure these boxes are checked (they're
   checked by default):
   `Patient_ID`, `Ages`, `Grade`, `Stage`, `Tumor_type`, `Histological_Subtype`,
   `Survival_Outcome`
5. Click **Generate Clinical Data**.

**What you should see:** a table appears below the button with 50 rows, one per
synthetic patient. Each row has a `Patient_ID` like `EC0000`, `EC0001`, ... `EC0049`,
along with an age, tumor grade, FIGO-style stage, tumor type, histological subtype,
and simulated survival outcome.

Try scrolling through a few rows. Notice that:
- Older patients are more likely to appear with Type 2 tumors when combined with
  other risk factors — the simulator isn't just assigning values randomly, it's
  applying the clinical dependency rules described in the README.
- Patients with early stage + Grade 1 disease are heavily weighted toward
  `Alive` in `Survival_Outcome`; Stage 4 patients are weighted toward `Deceased`.

6. Click **Download as CSV**. A file link appears — click it to download
   `OncoGenerator_clinical.csv` to your computer. This file contains **all**
   generated columns, even ones you didn't select for the on-screen table.

You now have a real, saved dataset of 50 synthetic patients.

---

## Step 3 — Try a filtered cohort (optional but instructive)

Let's see the filters in action before moving on.

1. Change **Stage Filter** to `Stage3`.
2. Click **Generate Clinical Data** again.

**What you should see:** the table now only contains patients whose `Stage` value
starts with `Stage3` (e.g. `Stage3C1`, `Stage3C2`). The cohort will likely be smaller
than 50 rows, since only a fraction of the original 50 patients fall into Stage 3.

3. Set **Stage Filter** back to `All` and **click Generate Clinical Data one more
   time** before continuing to Step 4 — the next steps assume you're working with
   the full 50-patient cohort, not the filtered one.

---

## Step 4 — Generate protein abundance data

Switch to the **Protein Abundance Data** tab.

1. Set **Number of Patients** to `10` (this generates protein data for the first 10
   patients in your current cohort — it doesn't need to match the clinical cohort
   size).
2. In **Protein List**, type:
   ```
   ANXA2, PKM2, ERBB2, EGFR, MMP9
   ```
   These are five real proteins with documented roles in endometrial cancer biology,
   already present in OncoGenerator' internal literature-curated dictionary.
3. Click **Generate Abundance Data**.

**What you should see:** a JSON output containing 10 patient blocks (`EC0000`
through `EC0009`), each with five values — one log2 fold-change number per
protein, all falling between **-3 and 3**. For example:

```json
{
  "EC0000": {
    "ANXA2": 1.87,
    "PKM2": 1.42,
    "ERBB2": 0.31,
    "EGFR": -0.22,
    "MMP9": 1.05
  },
  ...
}
```

Your exact numbers will differ slightly from these — they're sampled from
probability distributions, not fixed lookup values.

4. Click **Download as CSV**. This exports the same data as a patient × protein
   table (`OncoGenerator_protein.csv`), which is the format most downstream tools
   (R/Bioconductor, scikit-learn, etc.) expect.

**Try it yourself:** go back and try a protein that is *not* in the curated list,
such as `RANDOMGENE123`. You'll still get a value back for it — OncoGenerator
doesn't fail or skip unrecognized genes. A gene that is in neither the curated dictionary nor the HPA file gets a near-zero neutral value; a gene that is in HPA but not the dictionary follows the HPA-based UP/DOWN rule (see the README's "Known limitations").

---

## Step 5 — Visualize the clinical cohort

Switch to the **Distribution Plots** tab.

1. Click **Generate Plots**.

**What you should see:** a grid of histograms and bar charts, one per clinical
variable (age, BMI, grade, stage, tumor type, menopause, nulliparity, treatment,
survival outcome), summarizing the 50-patient cohort you generated in Step 2.

Look for:
- A roughly bell-shaped age distribution centered around 60
- A grade distribution that isn't perfectly even — this reflects the
  tumor-type-dependent grade weighting described in the README
- A stage distribution skewed toward earlier stages, since most of your 50
  synthetic patients were assigned Type 1 tumors (the more common, lower-grade
  pathway) unless they matched the "high-risk" profile

This tab is a quick sanity check — if you generate a cohort and something looks
obviously wrong (e.g. every patient the same age), this is where you'd spot it.

---

## Step 6 — Run the Random Forest stage classifier

Switch to the **Use Case: Stage Classifier** tab.

1. Click **Run Random Forest Classifier**.

**What you should see:** a text output reporting classification accuracy and a full
precision/recall/F1 breakdown for predicting **early stage** (Stage 1–2) vs.
**late stage** (Stage 3–4) disease, using only `Ages`, `BMI`, `Menopause`,
`Nulliparity`, `Tumor_type`, and `Grade` as input features — deliberately *not*
using `Stage` itself, since that would be the answer.

```
Accuracy: 84.0%

              precision    recall  f1-score   support

 Early Stage       0.83      0.91      0.87         6
  Late Stage       0.86      0.75      0.80         4

    accuracy                           0.84        10
   macro avg       0.84      0.83      0.83        10
weighted avg       0.84      0.84      0.84        10
```

Your exact numbers will differ — with only 50 patients and a fixed 80/20 split,
results can vary meaningfully. This demo is meant to show the workflow end-to-end,
not to be a rigorously validated benchmark (see the accompanying manuscript for the
full 1,000-patient, cross-validated version of this experiment).

> **If you see an error mentioning scikit-learn:** the script installs it
> automatically at the top, but if you skipped that cell or restarted your Colab
> runtime, re-run the very first cell (Step 1) before retrying this tab.

---

## What you've done

You've now run the full OncoGenerator pipeline end-to-end:

| Step | What it demonstrated |
|---|---|
| 2–3 | Rule-based synthetic clinical cohort generation, with subtype/stage filtering |
| 4 | Literature-informed, clinically-gated synthetic protein abundance simulation |
| 5 | Built-in visual QC of a generated cohort |
| 6 | A downstream ML use case trained entirely on synthetic data |

---

## Next steps

- Try generating a much larger cohort (e.g. 1,000 patients) to see how the
  classifier's accuracy stabilizes with more data.
- Try a longer protein list — the tool accepts any number of comma-separated gene
  symbols.
- If you have real Human Protein Atlas expression data, follow the
  "Optional: adding real HPA expression data" section in the README to make
  protein simulation more tissue-specific.
- For a full reference of every input, output, and underlying rule, see
  [`README.md`](README.md) in this repository.

---

## Getting help

If something in this tutorial doesn't match what you see, please open an issue at
[https://github.com/NehaAr/OncoGenerator/issues](https://github.com/NehaAr/OncoGenerator/issues)
with a screenshot and a description of the step where things diverged.
