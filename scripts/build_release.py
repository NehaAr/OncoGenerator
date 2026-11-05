#!/usr/bin/env python
"""Rebuild the v1.1.0 data release.   python scripts/build_release.py [--format csv.gz|parquet] [--out release/v1.1.0]
Deterministic: the same generator version + seeds reproduce the same files."""
import argparse, json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import oncogenerator as sp
from oncogenerator.io import write_dataset

COHORTS = [("core_n1000", 1000, 1, "clean", True),
           ("core_n10000", 10000, 2, "clean", True),
           ("core_n100000", 100000, 3, "clean", False)]
SCEN = ["noise_low", "noise_high", "imbalanced_late10", "batch5", "missing_mcar10", "missing_mnar20"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--format", default="csv.gz", choices=["csv.gz", "parquet", "csv"])
    ap.add_argument("--out", default=None)
    ap.add_argument("--lite", action="store_true",
                    help="small bundle: 1k and 10k core cohorts + scenarios at 2,000 patients (~10 MB)")
    a = ap.parse_args()
    a.out = a.out or ("release/v1.1.0-lite" if a.lite else "release/v1.1.0")
    out = Path(a.out); manifest = {"generator_version": sp.__version__, "format": a.format, "datasets": {}}
    jobs = [(n, N, s, sc, ex) for n, N, s, sc, ex in COHORTS if not (a.lite and N > 10000)]
    sn = 2000 if a.lite else 10000
    jobs += [(f"scenario_{sc}_n{sn}", sn, 2, sc, False) for sc in SCEN]
    for name, N, seed, sc, ex in jobs:
        t = time.time()
        ds = sp.generate_dataset(N, seed=seed, scenario=sc)
        sums = write_dataset(ds, out / name, a.format, with_expected=ex)
        manifest["datasets"][name] = dict(n_patients=int(len(ds.clinical)), n_proteins=ds.meta["n_proteins"],
                                          seed=seed, scenario=sc, files=sums)
        print(f"{name}: {time.time()-t:.1f}s", flush=True)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))

if __name__ == "__main__":
    main()
