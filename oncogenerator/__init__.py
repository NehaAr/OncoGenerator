"""OncoGenerator - synthetic clinical + proteomic data for endometrial cancer research."""
from .core import (__version__, SCENARIOS, Scenario, Dataset, generate_clinical,
                   generate_proteins, generate_dataset, default_panel, load_curated,
                   load_hpa, protein_table, contrast_truth, CONTRASTS, TRUE_DE_THRESHOLD,
                   LARGE_EFFECT_THRESHOLD)
from .io import write_dataset

__all__ = ["__version__", "SCENARIOS", "Scenario", "Dataset", "generate_clinical",
           "generate_proteins", "generate_dataset", "default_panel", "load_curated",
           "load_hpa", "protein_table", "contrast_truth", "CONTRASTS", "TRUE_DE_THRESHOLD", "LARGE_EFFECT_THRESHOLD",
           "write_dataset"]
