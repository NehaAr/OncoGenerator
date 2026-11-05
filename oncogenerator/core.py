# -*- coding: utf-8 -*-
"""
OncoGenerator core: seedable, vectorised synthetic clinical + proteomic data
generator for endometrial cancer research.

Author: Neha Arora (University of Newcastle). MIT License.

Design notes
------------
* Every random draw comes from a ``numpy.random.Generator`` created from the
  ``seed`` argument, so the same (seed, parameters) always give the same data and
  different seeds give independent replicates.
* ``legacy=True`` reproduces the *rules* of v1.0.0 (see CHANGELOG.md): factor
  matching that never fires for menopause/tumour-type/"stage3" terms, and an HPA
  prognostic parser that never returns "yes". Random streams differ from v1.0.0
  (which used global seeding), so legacy data are distributionally, not
  numerically, identical to v1.0.0 output.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

__version__ = "1.1.0"

PKG_DIR = Path(__file__).resolve().parent
REPO_DIR = PKG_DIR.parent
CURATED_CSV = PKG_DIR / "data" / "curated_regulation_v1.csv"
DEFAULT_HPA = REPO_DIR / "data" / "HPA.tsv"

HPA_EXPR_COL = "Tissue RNA - endometrium 1 [nTPM]"
HPA_PROG_COL = "Cancer prognostics - Uterine Corpus Endometrial Carcinoma (TCGA)"

# Pathway membership used for covariance blending (exact, case-sensitive symbols).
PATHWAY_MAP: Dict[str, List[str]] = {
    "PI3K_AKT": ["CTNNB", "EGFR", "ERBB2", "PKM2", "MMP9", "NFKB", "PTEN"],
    "APOPTOSIS": ["CASP3", "HMGB3", "PRDX1", "PRDX6", "SOD1"],
    "GLYCOLYSIS": ["ENO1", "ALDOA", "LDHA", "PKM", "TPI1", "PGAM2"],
    "WNT_BETA_CATENIN": ["CTNNB1", "APC", "AXIN1", "GSK3B"],
    "MAPK_ERK": ["KRAS", "BRAF", "MAPK1", "EGFR", "ERBB2"],
    "DNA_REPAIR": ["MLH1", "MSH2", "MSH6", "PMS2", "POLE", "TP53"],
    "CHROMATIN_REMODELING": ["ARID1A", "SMARCA4", "ARID1B"],
    "TGF_BETA_EMT": ["TGFB1", "SMAD2", "SMAD3", "SMAD4"],
    "ANGIOGENESIS": ["VEGFA", "KDR", "HIF1A", "ANGPT2"],
}

STAGE_CHOICES = ["Stage1", "Stage1A", "Stage1B", "Stage2", "Stage2A", "Stage2B",
                 "Stage3C1", "Stage3C2", "Stage4A", "Stage4B"]
_EARLY_W = np.array([0.1333] * 6 + [0.05] * 4)   # used for Grade1/Grade2
_LATE_W = np.array([0.05] * 6 + [0.1333] * 4)    # used for Grade3
_STAGE_KAPPA_IDX = {"Stage1": 1, "Stage1A": 1, "Stage1B": 2, "Stage2": 3, "Stage2A": 3,
                    "Stage2B": 4, "Stage3C1": 6, "Stage3C2": 7, "Stage4A": 9, "Stage4B": 10}

ETHNICITIES = ["European", "East Asian", "South Asian", "Middle Eastern", "African",
               "Indigenous Australian", "Other"]
ETHNICITY_W = np.array([0.55, 0.15, 0.10, 0.05, 0.03, 0.03, 0.09])

# name: (mean, sd, clip_lo, clip_hi)
SAMPLERS = {
    "up_strong": (2.0, 0.4, 0.0, 3.0),
    "up_weak": (1.0, 0.3, 0.0, 3.0),
    "down_strong": (-2.0, 0.4, -3.0, 0.0),
    "down_weak": (-1.0, 0.3, -3.0, 0.0),
    "low": (-0.5, 0.3, -1.5, 0.0),
    "neutral": (0.0, 0.5, -1.5, 1.5),
}

CONTRASTS = {
    "Late_vs_Early_stage": ("Late_Stage", 1, 0),
    "Grade3_vs_Grade1": ("Grade", "Grade3", "Grade1"),
    "Postmenopausal_vs_Premenopausal": ("Menopause", 2, 1),
    "Type2_vs_Type1": ("Tumor_type", "Type2", "Type1"),
}
LARGE_EFFECT_THRESHOLD = 0.5  # |expected log2FC difference| labelled "large effect"
TRUE_DE_THRESHOLD = LARGE_EFFECT_THRESHOLD  # backwards-compatible alias
NULL_TOL = 1e-6              # |effect| <= NULL_TOL => true null (no generated difference)
_CHUNK = 20000           # patients per protein-generation chunk (fixed => reproducible)


# ───────────────────────────── curated dictionary ─────────────────────────────
def load_curated() -> pd.DataFrame:
    """Curated regulation table (data/curated_regulation_v1.csv). See its
    ``alignment_status`` column: in v1.0.0 the Gene/Regulation/Factor lists had
    different lengths (230/212/203), so only the first 203 rows have all three."""
    return pd.read_csv(CURATED_CSV, keep_default_na=False)


def default_panel(curated: Optional[pd.DataFrame] = None) -> List[str]:
    """Unique (case-insensitive, first occurrence) genes whose Gene/Regulation/
    Factor entries are all present. 180 genes."""
    c = load_curated() if curated is None else curated
    c = c[(c.alignment_status == "complete") & (c.first_occurrence == 1)]
    return c.gene.tolist()


# ───────────────────────────────── HPA input ──────────────────────────────────
def load_hpa(path: Optional[Path] = None, legacy: bool = False) -> Dict[str, Tuple[str, str]]:
    """Return {gene_lower: (normal_tissue_regulation, prognostic_flag)}.

    normal_tissue_regulation: 'up' (>100 nTPM), 'moderate' (10-100), 'down' (<10).
    prognostic_flag: 'yes' if the HPA TCGA-UCEC prognostic p-value is <= 0.05.
    The legacy parser stripped every non-digit character, turning e.g.
    'unprognostic (3.21e-1)' into 3.211, so no gene was ever flagged 'yes'.
    Returns {} (neutral baseline) if the file is missing.
    """
    p = Path(path) if path else DEFAULT_HPA
    if not p.exists():
        return {}
    h = pd.read_csv(p, sep="\t")
    expr = h[HPA_EXPR_COL]
    norm = np.where(expr > 100, "up", np.where(expr >= 10, "moderate", "down"))
    raw = h[HPA_PROG_COL].astype(str).str.strip()
    if legacy:
        pv = pd.to_numeric(raw.str.replace(r"[^\d.]", "", regex=True), errors="coerce")
    else:
        pv = pd.to_numeric(raw.str.extract(r"\(\s*([0-9.]+(?:[eE][-+]?\d+)?)\s*\)")[0],
                           errors="coerce")
    prog = np.where(pv.notna() & (pv <= 0.05), "yes", "no")
    out: Dict[str, Tuple[str, str]] = {}
    for g, n_, p_ in zip(h["Gene"].astype(str).str.lower(), norm, prog):
        out.setdefault(g, (n_, p_))
    return out


# ─────────────────────────────── protein model ────────────────────────────────
def _choose_sampler(reg: str, norm_reg: str, prog: str) -> str:
    r, pr = str(reg).lower().strip(), str(prog).lower().strip()
    if r == "up":
        return "up_strong" if pr == "yes" else "up_weak"
    if r == "down":
        return "down_strong" if pr == "yes" else "down_weak"
    if "up" in r and "down" in r:
        if norm_reg == "up":
            return "up_strong" if pr == "yes" else "up_weak"
        if norm_reg == "down":
            return "down_strong" if pr == "yes" else "down_weak"
        return "neutral"
    if r == "low":
        return "low"
    return "neutral"


def _reg_class(reg: str) -> str:
    r = str(reg).lower().strip()
    if r in ("up", "down", "low"):
        return r.upper()
    if "up" in r and "down" in r:
        return "UP/DOWN"
    return "UNSPECIFIED"


def protein_table(panel: Sequence[str], hpa: Optional[Dict[str, Tuple[str, str]]] = None
                  ) -> pd.DataFrame:
    """Per-protein generating parameters = the per-protein ground truth."""
    cur = load_curated()
    cur = cur[cur.first_occurrence == 1].set_index(cur[cur.first_occurrence == 1].gene.str.lower())
    hpa = hpa or {}
    rows = []
    for prot in panel:
        pc = prot.strip()
        lo = pc.lower()
        in_dict = lo in cur.index and cur.loc[lo, "alignment_status"] == "complete"
        if in_dict:
            reg, fac = cur.loc[lo, "regulation_raw"], cur.loc[lo, "factor_raw"]
            ndist = int(cur.loc[lo, "n_distinct_regulation_values_for_gene"])
            status = "curated"
        else:
            reg, fac, ndist = "UP/DOWN", "NA", 0
            status = ("curated_but_incomplete_entry_treated_as_uncurated" if lo in cur.index
                      else "not_in_dictionary")
        norm_reg, prog = hpa.get(lo, ("moderate", "no"))
        samp = _choose_sampler(reg, norm_reg, prog)
        mean, sd, clo, chi = SAMPLERS[samp]
        pw = next((k for k, g in PATHWAY_MAP.items() if pc in g), "")
        rows.append(dict(protein=pc, curation_status=status, regulation_raw=reg,
                         regulation_class=_reg_class(reg), factor_raw=fac,
                         normal_tissue_regulation_hpa=norm_reg if lo in hpa else "",
                         prognostic_flag_hpa=prog if lo in hpa else "",
                         sampler=samp, sampler_mean=mean, sampler_sd=sd,
                         clip_lo=clo, clip_hi=chi, pathway=pw,
                         n_distinct_regulation_values_in_dictionary=ndist))
    return pd.DataFrame(rows)


def _parse_terms(factor: str) -> List[str]:
    if not factor or str(factor).upper() == "NA":
        return []
    return [p.strip().lower() for p in re.split(r"[/,]", str(factor)) if p.strip()]


def _factor_mask(factor: str, clin: pd.DataFrame, legacy: bool) -> np.ndarray:
    """Boolean per patient: does the protein's clinical factor apply?"""
    terms = _parse_terms(factor)
    n = len(clin)
    if not terms:
        return np.ones(n, dtype=bool)
    stage = clin["Stage"].str.lower().to_numpy()
    grade = clin["Grade"].str.lower().to_numpy()
    meno = clin["Menopause"].to_numpy()
    ttype = clin["Tumor_type"].str.lower().to_numpy()
    m = np.zeros(n, dtype=bool)
    if legacy:  # v1.0.0: exact match against [stage, grade, str(menopause)]
        for t in terms:
            m |= (stage == t) | (grade == t) | (meno.astype(str) == t)
        return m
    for t in terms:
        t = t.replace(" ", "")
        if t.startswith("stage"):
            m |= np.char.startswith(stage.astype(str), t)
        elif t.startswith("grade"):
            m |= grade == t
        elif t == "postmenopausal":
            m |= meno == 2
        elif t == "premenopausal":
            m |= meno == 1
        elif t.startswith("type"):
            m |= ttype == t
    return m


def mamdani_kappa(grade: pd.Series, stage: pd.Series) -> np.ndarray:
    """Vectorised Mamdani-style scale factor kappa in [0.5, 1.5] from grade & stage."""
    g = grade.astype(str).map(lambda x: 1.0 if "1" in x else (2.0 if "2" in x else 3.0)).to_numpy()
    s = stage.astype(str).map(lambda x: float(_STAGE_KAPPA_IDX.get(x, 5))).to_numpy()
    g_n, s_n = (g - 1.0) / 2.0, (s - 1.0) / 9.0
    r_high = np.minimum(g_n, s_n)
    r_low = np.minimum(1 - g_n, 1 - s_n)
    r_med = 1.0 - np.maximum(r_high, r_low)
    kappa = (r_low * 0.7 + r_med * 1.0 + r_high * 1.3) / (r_low + r_med + r_high + 1e-6)
    return np.clip(kappa, 0.5, 1.5)


@dataclass
class ProteinResult:
    values: pd.DataFrame        # Patient_ID + proteins, log2FC, float32
    expected: pd.DataFrame      # noise-free expected value kappa*mean (pre-clipping)
    table: pd.DataFrame         # per-protein generating parameters (ground truth)
    kappa: np.ndarray


def generate_proteins(clinical: pd.DataFrame, panel: Optional[Sequence[str]] = None, seed: int = 0,
                      rho: float = 0.4, hpa_path: Optional[Path] = None, legacy: bool = False,
                      noise_scale: float = 1.0) -> ProteinResult:
    """Simulate log2FC abundance for every (patient, protein).

    noise_scale multiplies every random component (sampler SD and pathway term);
    1.0 reproduces the documented model. The expected (signal) matrix is
    independent of noise_scale.
    """
    panel = list(default_panel() if panel is None else panel)
    hpa = load_hpa(hpa_path, legacy=legacy)
    tab = protein_table(panel, hpa)
    n, P = len(clinical), len(panel)
    rng = np.random.default_rng([seed, 7919])
    kappa_all = mamdani_kappa(clinical["Grade"], clinical["Stage"])

    pathways = list(PATHWAY_MAP)
    pw_idx = np.array([pathways.index(p) if p else -1 for p in tab.pathway])
    mean = tab.sampler_mean.to_numpy(float)
    sd = tab.sampler_sd.to_numpy(float)
    lo = tab.clip_lo.to_numpy(float)
    hi = tab.clip_hi.to_numpy(float)
    matched_all = np.column_stack([_factor_mask(f, clinical, legacy) for f in tab.factor_raw])

    vals = np.empty((n, P), dtype=np.float32)
    exp = np.empty((n, P), dtype=np.float32)
    for a in range(0, n, _CHUNK):
        b = min(a + _CHUNK, n)
        m = b - a
        matched = matched_all[a:b]
        kappa = kappa_all[a:b, None]
        z = rng.standard_normal((m, P))
        base_m = np.clip(mean + noise_scale * sd * z, lo, hi)
        zn = rng.standard_normal((m, P))
        base_n = np.clip(noise_scale * 0.5 * zn, -1.5, 1.5)
        base = np.where(matched, base_m, base_n)
        zk = rng.standard_normal((m, len(pathways)))
        eps = rng.standard_normal((m, P))
        in_pw = pw_idx >= 0
        shared = zk[:, np.where(in_pw, pw_idx, 0)]
        blend = rho * shared + math.sqrt(1.0 - rho ** 2) * eps
        final = base + np.where(in_pw, noise_scale * blend, 0.0)
        vals[a:b] = np.clip(kappa * final, -3.0, 3.0)
        exp[a:b] = kappa * np.where(matched, mean, 0.0)

    ids = clinical["Patient_ID"].to_numpy()
    v = pd.DataFrame(vals, columns=panel); v.insert(0, "Patient_ID", ids)
    e = pd.DataFrame(exp, columns=panel); e.insert(0, "Patient_ID", ids)
    return ProteinResult(v, e, tab, kappa_all)


# ─────────────────────────────── clinical model ───────────────────────────────
def _cat(rng: np.random.Generator, labels: Sequence[str], p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    p = p / p.sum(axis=1, keepdims=True)
    u = rng.random(len(p))[:, None]
    idx = np.minimum((u > np.cumsum(p, axis=1)).sum(axis=1), len(labels) - 1)
    return np.asarray(labels, dtype=object)[idx]


def _clinical_raw(n: int, rng: np.random.Generator) -> pd.DataFrame:
    ages = np.clip(rng.normal(60, 15, n).astype(int), 30, 85)
    bmi = np.round(np.clip(rng.normal(25, 4, n), 21, 40), 2)
    eth = _cat(rng, ETHNICITIES, np.tile(ETHNICITY_W, (n, 1)))
    treat = np.asarray(["Surgery", "Chemotherapy", "Radiotherapy", "Combined"], dtype=object)[
        rng.integers(0, 4, n)]
    myo = np.round(rng.uniform(2, 40, n), 2)
    meno = np.where(ages <= 51, 1, 2)
    nul = np.where(rng.random(n) < np.where(ages <= 35, 0.5, 0.6), "yes", "no").astype(object)
    high = (ages >= 60) & (nul == "yes") & (bmi >= 25)
    t2 = rng.random(n) < np.where(high, 0.6, 0.2)
    ttype = np.where(t2, "Type2", "Type1").astype(object)

    grades = ["Grade1", "Grade2", "Grade3"]
    gp = np.where(t2[:, None], np.array([0.2, 0.3, 0.5]), np.array([0.5, 0.3, 0.2]))
    grade = _cat(rng, grades, gp)
    sub1 = _cat(rng, ["Endometrioid", "Mucinous"], np.tile([0.9, 0.1], (n, 1)))
    sub2 = _cat(rng, ["Serous", "Clear Cell", "Undifferentiated"], np.tile([0.55, 0.35, 0.10], (n, 1)))
    sub = np.where(t2, sub2, sub1).astype(object)
    sp = np.where((grade == "Grade3")[:, None], _LATE_W, _EARLY_W)
    stage = _cat(rng, STAGE_CHOICES, sp)

    st = pd.Series(stage)
    p_alive = np.where(st.str.startswith("Stage1") & (grade == "Grade1"), 0.90,
                       np.where(st.str.startswith("Stage4"), 0.30, 0.65))
    surv = np.where(rng.random(n) < p_alive, "Alive", "Deceased").astype(object)
    late = (~st.str.startswith(("Stage1", "Stage2"))).astype(int).to_numpy()
    return pd.DataFrame(dict(Ages=ages, Ethnicity=eth, Menopause=meno, Grade=grade,
                             Tumor_type=ttype, Histological_Subtype=sub, Stage=stage,
                             Nulliparity=nul, BMI=bmi, Myometrial_mm=myo, Treatment=treat,
                             Survival_Outcome=surv, Late_Stage=late))


def generate_clinical(n: int, seed: int = 0, late_fraction: Optional[float] = None,
                      subtype: str = "All", stage_prefix: str = "All") -> pd.DataFrame:
    """Synthetic clinical cohort of n patients.

    late_fraction: if set, resample so that exactly round(f*n) patients are late
    stage (III/IV) - used for the class-imbalance scenario. subtype / stage_prefix
    filters are applied by oversampling, so n is still the returned size.
    """
    n = int(n)
    rng = np.random.default_rng([seed, 104729])
    if late_fraction is None and subtype == "All" and stage_prefix == "All":
        df = _clinical_raw(n, rng)
    else:
        pool = _clinical_raw(max(8 * n, 2000), rng)
        if subtype != "All":
            pool = pool[pool.Histological_Subtype == subtype]
        if stage_prefix != "All":
            pool = pool[pool.Stage.str.startswith(stage_prefix)]
        if late_fraction is not None:
            n_late = int(round(late_fraction * n))
            late_i, early_i = np.where(pool.Late_Stage == 1)[0], np.where(pool.Late_Stage == 0)[0]
            if len(late_i) < n_late or len(early_i) < n - n_late:
                raise ValueError("Filters/late_fraction too restrictive to draw %d patients" % n)
            take = np.concatenate([rng.choice(late_i, n_late, replace=False),
                                   rng.choice(early_i, n - n_late, replace=False)])
            rng.shuffle(take)
            df = pool.iloc[take]
        else:
            if len(pool) < n:
                raise ValueError("Filters too restrictive to draw %d patients" % n)
            df = pool.iloc[rng.choice(len(pool), n, replace=False)]
        df = df.reset_index(drop=True)
    w = max(4, len(str(n - 1)))
    df.insert(0, "Patient_ID", [f"EC{i:0{w}d}" for i in range(len(df))])
    df.insert(df.columns.get_loc("Menopause") + 1, "Menopause_label",
              np.where(df.Menopause == 1, "pre_or_peri", "post"))
    return df


# ───────────────────────────────── scenarios ──────────────────────────────────
def apply_batch_effects(values: pd.DataFrame, n_batches: int, seed: int, shift_sd: float = 0.3,
                        gain_sd: float = 0.1) -> Tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    """Add per-batch, per-protein additive shifts and per-batch multiplicative gain
    (applied AFTER the +/-3 clip, so values may exceed +/-3). Returns
    (new values, batch labels, ground-truth table of the applied effects)."""
    rng = np.random.default_rng([seed, 31337])
    prot = [c for c in values.columns if c != "Patient_ID"]
    n = len(values)
    batch = rng.integers(0, n_batches, n)
    shift = rng.normal(0, shift_sd, (n_batches, len(prot)))
    gain = 1.0 + rng.normal(0, gain_sd, n_batches)
    x = values[prot].to_numpy(np.float64)
    x = x * gain[batch, None] + shift[batch]
    out = values.copy()
    out[prot] = x.astype(np.float32)
    rows = [dict(batch=f"B{b + 1}", protein=p, additive_shift=float(shift[b, j]),
                 batch_gain=float(gain[b]))
            for b in range(n_batches) for j, p in enumerate(prot)]
    return out, np.array([f"B{b + 1}" for b in batch], dtype=object), pd.DataFrame(rows)


def apply_missingness(values: pd.DataFrame, seed: int, mcar: float = 0.0, mnar: float = 0.0
                      ) -> pd.DataFrame:
    """Set entries to NaN. mcar: fraction missing completely at random.
    mnar: overall fraction missing, with probability increasing as the value
    decreases (logistic, slope 1.5), mimicking low-abundance dropout."""
    rng = np.random.default_rng([seed, 4242])
    prot = [c for c in values.columns if c != "Patient_ID"]
    x = values[prot].to_numpy(np.float32).copy()
    miss = np.zeros(x.shape, dtype=bool)
    if mcar > 0:
        miss |= rng.random(x.shape) < mcar
    if mnar > 0:
        lo_a, hi_a = -20.0, 20.0
        for _ in range(60):
            a = (lo_a + hi_a) / 2
            if (1 / (1 + np.exp(-(a - 1.5 * x)))).mean() > mnar:
                hi_a = a
            else:
                lo_a = a
        p = 1 / (1 + np.exp(-(a - 1.5 * x)))
        miss |= rng.random(x.shape) < p
    x[miss] = np.nan
    out = values.copy()
    out[prot] = x
    return out


@dataclass
class Scenario:
    name: str
    noise_scale: float = 1.0
    late_fraction: Optional[float] = None
    n_batches: int = 1
    batch_shift_sd: float = 0.3
    missing_mcar: float = 0.0
    missing_mnar: float = 0.0
    description: str = ""


SCENARIOS: Dict[str, Scenario] = {s.name: s for s in [
    Scenario("clean", description="Default generator settings."),
    Scenario("noise_low", noise_scale=0.5, description="All random components x0.5."),
    Scenario("noise_high", noise_scale=2.0, description="All random components x2.0."),
    Scenario("imbalanced_late10", late_fraction=0.10,
             description="Only 10% late-stage (III/IV) patients."),
    Scenario("batch5", n_batches=5,
             description="5 batches with additive protein shifts (SD 0.3) and gains (SD 0.1)."),
    Scenario("missing_mcar10", missing_mcar=0.10,
             description="10% of protein values missing completely at random."),
    Scenario("missing_mnar20", missing_mnar=0.20,
             description="~20% missing, more likely for low values (MNAR)."),
]}


@dataclass
class Dataset:
    clinical: pd.DataFrame
    proteins: pd.DataFrame
    expected: pd.DataFrame
    protein_truth: pd.DataFrame
    contrast_truth: pd.DataFrame
    batch_truth: Optional[pd.DataFrame]
    meta: dict = field(default_factory=dict)


def contrast_truth(clinical: pd.DataFrame, expected: pd.DataFrame) -> pd.DataFrame:
    prot = [c for c in expected.columns if c != "Patient_ID"]
    e = expected[prot].to_numpy(np.float64)
    rows = []
    for name, (col, b, a) in CONTRASTS.items():
        mb, ma = (clinical[col] == b).to_numpy(), (clinical[col] == a).to_numpy()
        if mb.sum() == 0 or ma.sum() == 0:
            continue
        eff = e[mb].mean(0) - e[ma].mean(0)
        for p, v in zip(prot, eff):
            rows.append(dict(contrast=name, protein=p, n_group_B=int(mb.sum()), n_group_A=int(ma.sum()),
                             true_effect_log2fc=round(float(v), 5),
                             is_nonzero_effect=bool(abs(v) > NULL_TOL),
                             is_large_effect=bool(abs(v) >= LARGE_EFFECT_THRESHOLD)))
    return pd.DataFrame(rows)


def generate_dataset(n: int, seed: int = 0, scenario: Scenario | str = "clean",
                     panel: Optional[Sequence[str]] = None, rho: float = 0.4,
                     hpa_path: Optional[Path] = None, legacy: bool = False) -> Dataset:
    sc = SCENARIOS[scenario] if isinstance(scenario, str) else scenario
    clin = generate_clinical(n, seed, late_fraction=sc.late_fraction)
    res = generate_proteins(clin, panel, seed, rho, hpa_path, legacy, sc.noise_scale)
    values, batch_truth = res.values, None
    if sc.n_batches > 1:
        values, batch, batch_truth = apply_batch_effects(values, sc.n_batches, seed, sc.batch_shift_sd)
        clin = clin.copy()
        clin["Batch"] = batch
    if sc.missing_mcar > 0 or sc.missing_mnar > 0:
        values = apply_missingness(values, seed, sc.missing_mcar, sc.missing_mnar)
    ct = contrast_truth(clin, res.expected)
    meta = dict(generator_version=__version__, n=int(len(clin)), seed=int(seed), rho=rho,
                legacy=legacy, scenario=sc.name, scenario_params={k: v for k, v in sc.__dict__.items()},
                n_proteins=int(values.shape[1] - 1), large_effect_threshold=LARGE_EFFECT_THRESHOLD, null_tolerance=NULL_TOL)
    return Dataset(clin, values, res.expected, res.table, ct, batch_truth, meta)
