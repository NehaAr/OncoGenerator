"""Writing datasets to disk (csv.gz by default; Parquet if pyarrow is installed)."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict

import pandas as pd

from .core import Dataset


def _write(df: pd.DataFrame, path_noext: Path, fmt: str, round_floats: int = 4) -> Path:
    if fmt == "parquet":
        try:
            import pyarrow  # noqa: F401
        except ImportError as e:  # pragma: no cover
            raise SystemExit("Parquet output needs pyarrow: pip install pyarrow") from e
        p = path_noext.with_suffix(".parquet")
        df.to_parquet(p, index=False, compression="zstd")
        return p
    d = df.copy()
    for c in d.select_dtypes("float").columns:
        d[c] = d[c].round(round_floats)
    if fmt == "csv.gz":
        p = Path(str(path_noext) + ".csv.gz")
        d.to_csv(p, index=False, compression={"method": "gzip", "mtime": 0})
    elif fmt == "csv":
        p = Path(str(path_noext) + ".csv")
        d.to_csv(p, index=False)
    else:
        raise ValueError(f"unknown format {fmt!r}")
    return p


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_dataset(ds: Dataset, outdir, fmt: str = "csv.gz", with_expected: bool = True) -> Dict[str, str]:
    """Write all tables of a Dataset into outdir; returns {filename: sha256}."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    files = [_write(ds.clinical, out / "clinical", fmt),
             _write(ds.proteins, out / "proteins_log2fc", fmt),
             _write(ds.protein_truth, out / "ground_truth_proteins", fmt),
             _write(ds.contrast_truth, out / "ground_truth_contrasts", fmt)]
    if with_expected:
        files.append(_write(ds.expected, out / "expected_log2fc_noise_free", fmt))
    if ds.batch_truth is not None:
        files.append(_write(ds.batch_truth, out / "ground_truth_batch_effects", fmt))
    (out / "run_metadata.json").write_text(json.dumps(ds.meta, indent=2, default=str))
    files.append(out / "run_metadata.json")
    return {f.name: sha256(f) for f in files}
