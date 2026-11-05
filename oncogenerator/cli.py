"""Command line interface:  python -m oncogenerator generate --n 10000 --seed 2 --out out/"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .core import SCENARIOS, Scenario, __version__, default_panel, generate_dataset
from .io import write_dataset


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="oncogenerator", description=__doc__)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate", help="generate one synthetic dataset")
    g.add_argument("--n", type=int, required=True, help="number of patients")
    g.add_argument("--seed", type=int, default=0, help="random seed (change for independent replicates)")
    g.add_argument("--out", required=True, help="output directory")
    g.add_argument("--format", default="csv.gz", choices=["csv.gz", "csv", "parquet"])
    g.add_argument("--scenario", default="clean", choices=sorted(SCENARIOS),
                   help="preset scenario (overridable by the options below)")
    g.add_argument("--noise-scale", type=float)
    g.add_argument("--late-fraction", type=float, help="target fraction of late-stage patients")
    g.add_argument("--batches", type=int, help="number of batches (>1 adds batch effects)")
    g.add_argument("--missing-mcar", type=float)
    g.add_argument("--missing-mnar", type=float)
    g.add_argument("--rho", type=float, default=0.4, help="pathway correlation (0-1)")
    g.add_argument("--proteins", help="comma-separated gene symbols (default: 180 curated genes)")
    g.add_argument("--hpa", help="path to HPA.tsv (default: data/HPA.tsv in the repo if present)")
    g.add_argument("--legacy", action="store_true", help="reproduce v1.0.0 factor/HPA rules")
    g.add_argument("--no-expected", action="store_true", help="skip the noise-free expected matrix")
    sub.add_parser("list-scenarios", help="show preset scenarios")
    sub.add_parser("panel", help="print the default 180-gene panel")
    a = ap.parse_args(argv)

    if a.cmd == "list-scenarios":
        for s in SCENARIOS.values():
            print(f"{s.name:20s} {s.description}")
        return
    if a.cmd == "panel":
        print(",".join(default_panel()))
        return
    base = SCENARIOS[a.scenario]
    sc = Scenario(**{**base.__dict__})
    if a.noise_scale is not None: sc.noise_scale = a.noise_scale
    if a.late_fraction is not None: sc.late_fraction = a.late_fraction
    if a.batches is not None: sc.n_batches = a.batches
    if a.missing_mcar is not None: sc.missing_mcar = a.missing_mcar
    if a.missing_mnar is not None: sc.missing_mnar = a.missing_mnar
    panel = [p.strip() for p in a.proteins.split(",") if p.strip()] if a.proteins else None
    ds = generate_dataset(a.n, a.seed, sc, panel, a.rho, Path(a.hpa) if a.hpa else None, a.legacy)
    sums = write_dataset(ds, a.out, a.format, with_expected=not a.no_expected)
    print(json.dumps({"out": a.out, "files": sums}, indent=2))
