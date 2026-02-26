"""ICOR dataset loader.

~40 CDS sequences with heterologous expression levels across different codon
optimization strategies. Best for COS subscale validation.

Data assembled from published codon optimization benchmarks.
"""

from __future__ import annotations

import csv
import io
import logging
from pathlib import Path

import requests

from calibration.config import DATASET_DIR, DOWNLOAD_TIMEOUT
from calibration.datasets.base import BenchmarkRecord, BenchmarkDataset

logger = logging.getLogger(__name__)

ICOR_CACHE = DATASET_DIR / "icor"

# ICOR data from Supplementary Table — codon optimization comparison
ICOR_URL = (
    "https://raw.githubusercontent.com/laurencefreeman/ICOR/"
    "main/data/expression_data.csv"
)


def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def download(force: bool = False) -> Path:
    """Download ICOR data to local cache."""
    ICOR_CACHE.mkdir(parents=True, exist_ok=True)
    out_path = ICOR_CACHE / "expression_data.csv"

    if out_path.exists() and not force:
        logger.info("ICOR data already cached at %s", out_path)
        return out_path

    logger.info("Downloading ICOR data from %s", ICOR_URL)
    try:
        resp = requests.get(ICOR_URL, timeout=DOWNLOAD_TIMEOUT)
        resp.raise_for_status()
        out_path.write_bytes(resp.content)
        logger.info("Saved %d bytes to %s", len(resp.content), out_path)
    except requests.RequestException as e:
        logger.warning("Could not download ICOR data: %s", e)
        logger.info("You can manually place expression_data.csv in %s", ICOR_CACHE)
        raise

    return out_path


def load(cache_path: Path | None = None) -> BenchmarkDataset:
    """Load ICOR dataset from cached CSV.

    Expected columns (flexible):
      - gene / name / construct: identifier
      - cds / coding_sequence: CDS sequence (DNA)
      - utr / five_prime_utr: 5' UTR (optional)
      - expression / protein_level / te: measured expression
    """
    if cache_path is None:
        cache_path = ICOR_CACHE / "expression_data.csv"

    if not cache_path.exists():
        raise FileNotFoundError(
            f"ICOR data not found at {cache_path}. "
            "Run `python -m calibration.cli download --dataset icor` first."
        )

    records = []
    text = cache_path.read_text(encoding="utf-8")
    reader = csv.DictReader(io.StringIO(text))

    def _get(row: dict, *keys: str) -> str:
        for k in keys:
            if k in row and row[k]:
                return row[k].strip()
        return ""

    for row in reader:
        name = _get(row, "gene", "name", "construct", "id")
        cds = _get(row, "cds", "coding_sequence", "cds_sequence")
        utr = _get(row, "utr", "five_prime_utr", "utr_sequence", "5utr")
        expr_str = _get(row, "expression", "protein_level", "te",
                        "expression_level", "fluorescence")

        if not expr_str or not cds:
            continue

        try:
            expr_val = float(expr_str)
        except ValueError:
            continue

        cds_rna = _dna_to_rna(cds)
        utr_rna = _dna_to_rna(utr) if utr else ""

        # Strip stop codon
        if len(cds_rna) >= 3 and cds_rna[-3:] in ("UAA", "UAG", "UGA"):
            cds_rna = cds_rna[:-3]

        records.append(BenchmarkRecord(
            record_id=name or f"icor_{len(records)}",
            utr_sequence=utr_rna,
            cds_sequence=cds_rna,
            measured_te=expr_val,
            expression_system="mammalian",
            has_full_cds=True,
        ))

    logger.info("Loaded %d ICOR records", len(records))

    return BenchmarkDataset(
        name="icor",
        records=records,
        te_metric_name="protein_expression",
        description="ICOR: ~40 CDS with heterologous expression levels",
        citation="ICOR codon optimization benchmark",
    )
