"""Sample MPRA dataset loader.

Sample et al. (2019) — ~100K designed 5' UTRs with mean ribosome load (MRL)
using a fixed EGFP reporter CDS. Best for UTR subscale calibration.

Data source: GEO GSE114002 (GSM3130443_designed_library.csv.gz)
"""

from __future__ import annotations

import csv
import gzip
import io
import logging
from pathlib import Path

import requests

from calibration.config import DATASET_DIR, DOWNLOAD_TIMEOUT, MPRA_SUBSAMPLE_SIZE
from calibration.datasets.base import BenchmarkRecord, BenchmarkDataset

logger = logging.getLogger(__name__)

MPRA_CACHE = DATASET_DIR / "sample_mpra"

# GEO download URLs for Sample et al. 2019
DESIGNED_LIBRARY_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?acc=GSM3130443&format=file&file=GSM3130443%5Fdesigned%5Flibrary%2Ecsv%2Egz"
)

# Fixed EGFP CDS used as reporter in the MPRA assay.
# This is a standard EGFP; since CDS is constant across all samples,
# CDS-dependent subscales cannot differentiate — used for UTR calibration only.
EGFP_CDS_RNA = (
    "AUGGUGAGCAAGGGCGAGGAGCUGUUCACCGGGGUGGUGCCCAUCCUGGUCGAGCUGGACG"
    "GCGACGUAAACGGCCACAAGUUCAGCGUGUCCGGCGAGGGCGAGGGCGAUGCCACCUACGG"
    "CAAGCUGACCCUGAAGUUCAUCUGCACCACCGGCAAGCUGCCCGUGCCCUGGCCCACCCUC"
    "GUGACCACCCUGACCUACGGCGUGCAGUGCUUCAGCCGCUACCCCGACCACAUGAAGCAGC"
    "ACGACUUCUUCAAGUCCGCCAUGCCCGAAGGCUACGUCCAGGAGCGCACCAUCUUCUUCAA"
    "GGACGACGGCAACUACAAGACCCGCGCCGAGGUGAAGUUCGAGGGCGACACCCUGGUGAAC"
    "CGCAUCGAGCUGAAGGGCAUCGACUUCAAGGAGGACGGCAACAUCCUGGGGCACAAGCUGG"
    "AGUACAACUACAACAGCCACAACGUCUAUAUCAUGGCCGACAAGCAGAAGAACGGCAUCAAG"
    "GUGAACUUCAAGAUCCGCCACAACAUCGAGGACGGCAGCGUGCAGCUCGCCGACCACUACC"
    "AGCAGAACACCCCCAUCGGCGACGGCCCCGUGCUGCUGCCCGACAACCACUACCUGAGCAC"
    "CCAGUCCGCCCUGAGCAAAGACCCCAACGAGAAGCGCGAUCACAUGGUCCUGCUGGAGUUC"
    "GUGACCGCCGCCGGGAUCACUCUCGGCAUGGACGAGCUGUACAAG"
)


def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def download(force: bool = False) -> Path:
    """Download Sample MPRA designed library from GEO."""
    MPRA_CACHE.mkdir(parents=True, exist_ok=True)
    gz_path = MPRA_CACHE / "designed_library.csv.gz"
    csv_path = MPRA_CACHE / "designed_library.csv"

    if csv_path.exists() and not force:
        logger.info("Sample MPRA data already cached at %s", csv_path)
        return csv_path

    logger.info("Downloading Sample MPRA designed library from GEO GSM3130443...")
    resp = requests.get(DESIGNED_LIBRARY_URL, timeout=DOWNLOAD_TIMEOUT, stream=True)
    resp.raise_for_status()

    gz_path.write_bytes(resp.content)

    # Decompress
    with gzip.open(gz_path, "rt", encoding="utf-8") as gz_f:
        csv_path.write_text(gz_f.read(), encoding="utf-8")

    logger.info("Saved %d bytes to %s", csv_path.stat().st_size, csv_path)
    return csv_path


def load(cache_path: Path | None = None, subsample: int | None = None) -> BenchmarkDataset:
    """Load Sample MPRA dataset.

    Actual CSV columns from GEO GSM3130443:
      - utr: 5' UTR sequence (DNA, 50 nt)
      - rl: mean ribosome load (float)
      - library: type (snv, human_utrs, etc.)
    """
    if cache_path is None:
        csv_path = MPRA_CACHE / "designed_library.csv"
        gz_path = MPRA_CACHE / "designed_library.csv.gz"
        if csv_path.exists():
            cache_path = csv_path
        elif gz_path.exists():
            cache_path = gz_path
        else:
            raise FileNotFoundError(
                f"Sample MPRA data not found. "
                "Run `python -m calibration.cli download --dataset sample_mpra` first."
            )

    if not cache_path.exists():
        raise FileNotFoundError(f"Sample MPRA data not found at {cache_path}")

    records = []

    # Handle both .csv and .csv.gz
    if str(cache_path).endswith(".gz"):
        opener = lambda: gzip.open(cache_path, "rt", encoding="utf-8")
    else:
        opener = lambda: open(cache_path, "r", encoding="utf-8")

    with opener() as f:
        reader = csv.DictReader(f)

        for i, row in enumerate(reader):
            utr_seq = row.get("utr", "").strip()
            mrl_str = row.get("rl", "").strip()

            if not utr_seq or not mrl_str:
                continue

            try:
                mrl = float(mrl_str)
            except ValueError:
                continue

            utr_rna = _dna_to_rna(utr_seq)

            records.append(BenchmarkRecord(
                record_id=f"mpra_{i}",
                utr_sequence=utr_rna,
                cds_sequence=EGFP_CDS_RNA,
                measured_te=mrl,
                expression_system="mammalian",
                has_full_cds=False,
            ))

    logger.info("Loaded %d Sample MPRA records", len(records))

    dataset = BenchmarkDataset(
        name="sample_mpra",
        records=records,
        te_metric_name="mean_ribosome_load",
        description=f"Sample MPRA: {len(records)} designed 5' UTRs with MRL (Sample 2019)",
        citation="Sample et al. (2019) Nature Biotechnology 37:803-809",
    )

    n = subsample if subsample is not None else MPRA_SUBSAMPLE_SIZE
    if n and len(records) > n:
        dataset = dataset.subsample(n)

    return dataset
