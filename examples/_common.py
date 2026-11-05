"""Shared helpers for the example/benchmark scripts."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
import pandas as pd
import oncogenerator as sp

def load_or_generate(release_dir=None, n=1000, seed=1, scenario="clean"):
    """Return (clinical, proteins, contrast_truth). Reads a release folder if given,
    otherwise generates in memory with the same generator."""
    if release_dir:
        d = Path(release_dir)
        rd = lambda f: pd.read_csv(next(d.glob(f + ".csv*")))
        return rd("clinical"), rd("proteins_log2fc"), rd("ground_truth_contrasts")
    ds = sp.generate_dataset(n, seed=seed, scenario=scenario)
    return ds.clinical, ds.proteins, ds.contrast_truth

def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); r = p[o] * n / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(r[::-1])[::-1]
    return np.clip(q, 0, 1)
