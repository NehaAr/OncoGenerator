import numpy as np
import pandas as pd
import pytest

import oncogenerator as sp
from oncogenerator import core


def test_seed_reproducible_and_independent():
    a = sp.generate_dataset(300, seed=1)
    b = sp.generate_dataset(300, seed=1)
    c = sp.generate_dataset(300, seed=2)
    pd.testing.assert_frame_equal(a.clinical, b.clinical)
    pd.testing.assert_frame_equal(a.proteins, b.proteins)
    assert not a.clinical.Ages.equals(c.clinical.Ages)
    assert not a.proteins.iloc[:, 1:].equals(c.proteins.iloc[:, 1:])


def test_ids_join_and_bounds():
    d = sp.generate_dataset(500, seed=3)
    assert d.clinical.Patient_ID.tolist() == d.proteins.Patient_ID.tolist()
    v = d.proteins.drop(columns="Patient_ID").to_numpy()
    assert v.min() >= -3 and v.max() <= 3 and not np.isnan(v).any()
    assert d.clinical.Ages.between(30, 85).all() and d.clinical.BMI.between(21, 40).all()


def test_dictionary_alignment_is_explicit():
    cur = sp.load_curated()
    assert (cur.alignment_status == "complete").sum() == 203
    assert len(sp.default_panel()) == 180
    assert set(sp.default_panel()).isdisjoint({"CAPS", "DES", "CDH13", "MLH1"})


def test_strict_factors_fire_but_legacy_do_not():
    clin = sp.generate_clinical(2000, seed=4)
    post = core._factor_mask("POSTMENOPAUSAL", clin, legacy=False)
    assert post.sum() == (clin.Menopause == 2).sum()
    assert core._factor_mask("POSTMENOPAUSAL", clin, legacy=True).sum() == 0
    assert core._factor_mask("TYPE 2", clin, legacy=False).sum() == (clin.Tumor_type == "Type2").sum()
    assert core._factor_mask("TYPE 2", clin, legacy=True).sum() == 0
    s3 = core._factor_mask("stage3", clin, legacy=False)
    assert s3.sum() == clin.Stage.str.startswith("Stage3").sum()
    assert core._factor_mask("stage3", clin, legacy=True).sum() == 0


def test_hpa_parser_fix():
    assert sum(v[1] == "yes" for v in sp.load_hpa(legacy=True).values()) == 0
    assert sum(v[1] == "yes" for v in sp.load_hpa(legacy=False).values()) > 1000


def test_scenarios():
    n = 1000
    imb = sp.generate_dataset(n, seed=5, scenario="imbalanced_late10")
    assert imb.clinical.Late_Stage.sum() == 100
    bat = sp.generate_dataset(n, seed=5, scenario="batch5")
    assert bat.clinical.Batch.nunique() == 5 and bat.batch_truth is not None
    m1 = sp.generate_dataset(n, seed=5, scenario="missing_mcar10").proteins.drop(columns="Patient_ID")
    assert abs(m1.isna().to_numpy().mean() - 0.10) < 0.01
    m2 = sp.generate_dataset(n, seed=5, scenario="missing_mnar20")
    x = m2.proteins.drop(columns="Patient_ID")
    assert abs(x.isna().to_numpy().mean() - 0.20) < 0.02
    clean = sp.generate_dataset(n, seed=5).proteins.drop(columns="Patient_ID")
    missing_vals = clean.to_numpy()[x.isna().to_numpy()]
    assert missing_vals.mean() < clean.to_numpy()[~x.isna().to_numpy()].mean()  # MNAR: lower values go missing
    lo = sp.generate_dataset(n, seed=5, scenario="noise_low")
    hi = sp.generate_dataset(n, seed=5, scenario="noise_high")
    pd.testing.assert_frame_equal(lo.expected, hi.expected)  # signal independent of noise
    assert hi.proteins.iloc[:, 1:].std().mean() > lo.proteins.iloc[:, 1:].std().mean()


def test_ground_truth_recovers_signal():
    d = sp.generate_dataset(4000, seed=6)
    ct = d.contrast_truth
    t = ct[ct.contrast == "Late_vs_Early_stage"].set_index("protein").true_effect_log2fc
    late = d.clinical.Late_Stage.to_numpy() == 1
    obs = d.proteins.drop(columns="Patient_ID")
    emp = obs[late].mean() - obs[~late].mean()
    big = t.abs() >= 0.5
    assert np.corrcoef(emp[big.index[big]], t[big])[0, 1] > 0.9
