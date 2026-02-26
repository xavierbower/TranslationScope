"""PERSIST-seq dataset loader.

Wayment-Steele et al. (2022) — 233 full-length mRNA constructs with ribosome
load measurements. Contains 5' UTR and CDS variants, best direct analog to
TranslationScope usage.

Data source: Stanford Digital Repository / supplementary tables.
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

PERSIST_CACHE = DATASET_DIR / "persist_seq"

# Supplementary table URL (Wayment-Steele 2022, eLife)
# The dataset is in Supplementary file 3 — mRNA-level ribosome load data
PERSIST_URL = (
    "https://raw.githubusercontent.com/eternagame/SARSCoV2_ribosome_loading/"
    "main/data/ribosome_loading_data.csv"
)


def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def download(force: bool = False) -> Path:
    """Download PERSIST-seq data to local cache."""
    PERSIST_CACHE.mkdir(parents=True, exist_ok=True)
    out_path = PERSIST_CACHE / "ribosome_loading_data.csv"

    if out_path.exists() and not force:
        logger.info("PERSIST-seq data already cached at %s", out_path)
        return out_path

    logger.info("Downloading PERSIST-seq data from %s", PERSIST_URL)
    resp = requests.get(PERSIST_URL, timeout=DOWNLOAD_TIMEOUT)
    resp.raise_for_status()
    out_path.write_bytes(resp.content)
    logger.info("Saved %d bytes to %s", len(resp.content), out_path)
    return out_path


def load(cache_path: Path | None = None) -> BenchmarkDataset:
    """Load PERSIST-seq dataset from cached CSV.

    Expected CSV columns (flexible — we adapt to what's available):
      - construct_name or name: identifier
      - utr_sequence or five_prime_utr: 5' UTR
      - cds_sequence or cds: coding sequence
      - ribosome_load or rl or te: measured translation metric
    """
    if cache_path is None:
        cache_path = PERSIST_CACHE / "ribosome_loading_data.csv"

    if not cache_path.exists():
        raise FileNotFoundError(
            f"PERSIST-seq data not found at {cache_path}. "
            "Run `python -m calibration.cli download --dataset persist_seq` first."
        )

    records = []
    text = cache_path.read_text(encoding="utf-8")
    reader = csv.DictReader(io.StringIO(text))

    # Map flexible column names
    def _get(row: dict, *keys: str) -> str:
        for k in keys:
            if k in row and row[k]:
                return row[k].strip()
        return ""

    for row in reader:
        name = _get(row, "construct_name", "name", "id", "construct")
        utr = _get(row, "utr_sequence", "five_prime_utr", "5utr", "utr")
        cds = _get(row, "cds_sequence", "cds", "coding_sequence")
        te_str = _get(row, "ribosome_load", "rl", "te", "translation_efficiency",
                       "mean_ribosome_load")

        if not te_str:
            continue

        try:
            te_val = float(te_str)
        except ValueError:
            continue

        # Convert DNA to RNA if needed
        utr = _dna_to_rna(utr)
        cds = _dna_to_rna(cds)

        # Strip terminal stop codon if present
        if len(cds) >= 3 and cds[-3:] in ("UAA", "UAG", "UGA"):
            cds = cds[:-3]

        if not utr and not cds:
            continue

        records.append(BenchmarkRecord(
            record_id=name or f"persist_{len(records)}",
            utr_sequence=utr,
            cds_sequence=cds,
            measured_te=te_val,
            expression_system="mammalian",
            has_full_cds=bool(cds),
        ))

    logger.info("Loaded %d PERSIST-seq records", len(records))

    return BenchmarkDataset(
        name="persist_seq",
        records=records,
        te_metric_name="ribosome_load",
        description="PERSIST-seq: 233 full-length mRNAs with ribosome load (Wayment-Steele 2022)",
        citation="Wayment-Steele et al. (2022) eLife 11:e79775",
    )
