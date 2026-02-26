"""Zheng RiboNN dataset loader.

Zheng et al. (2025) — Transcriptome-wide TE atlas from ribosome profiling.
Requires joining with GENCODE annotations to get UTR + CDS sequences.

This is the most complex loader as it requires:
1. TE values from RiboNN output
2. GENCODE GTF for transcript annotations
3. GENCODE FASTA for transcript sequences
"""

from __future__ import annotations

import csv
import gzip
import io
import logging
from pathlib import Path

import requests

from calibration.config import DATASET_DIR, DOWNLOAD_TIMEOUT
from calibration.datasets.base import BenchmarkRecord, BenchmarkDataset

logger = logging.getLogger(__name__)

RIBONN_CACHE = DATASET_DIR / "zheng_ribonn"

# Placeholder URLs — actual data requires download from paper supplement
RIBONN_TE_URL = (
    "https://raw.githubusercontent.com/zhenghp/RiboNN/"
    "main/data/hela_te_values.tsv"
)

# GENCODE human transcript sequences (for joining)
GENCODE_FASTA_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_44/gencode.v44.transcripts.fa.gz"
)


def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def download(force: bool = False) -> Path:
    """Download RiboNN TE data to local cache.

    Note: GENCODE annotations are large (~1.5 GB). This downloads only the
    pre-processed TE values. Full sequence joining requires running the
    join step separately.
    """
    RIBONN_CACHE.mkdir(parents=True, exist_ok=True)
    te_path = RIBONN_CACHE / "hela_te_values.tsv"

    if te_path.exists() and not force:
        logger.info("RiboNN TE data already cached at %s", te_path)
        return te_path

    logger.info("Downloading RiboNN TE data from %s", RIBONN_TE_URL)
    try:
        resp = requests.get(RIBONN_TE_URL, timeout=DOWNLOAD_TIMEOUT)
        resp.raise_for_status()
        te_path.write_bytes(resp.content)
        logger.info("Saved %d bytes to %s", len(resp.content), te_path)
    except requests.RequestException as e:
        logger.warning("Could not download RiboNN data: %s", e)
        logger.info(
            "This dataset requires manual preparation. Place a TSV with "
            "columns [transcript_id, te, utr_sequence, cds_sequence] in %s",
            RIBONN_CACHE,
        )
        raise

    return te_path


def load(cache_path: Path | None = None) -> BenchmarkDataset:
    """Load RiboNN dataset.

    Expects a TSV with columns:
      - transcript_id / gene_name: identifier
      - te / translation_efficiency: measured TE (log2 ribo/RNA)
      - utr_sequence / five_prime_utr: 5' UTR
      - cds_sequence / cds: coding sequence

    If the raw TE file doesn't have sequences, it must first be joined with
    GENCODE annotations (see `join_with_gencode()`).
    """
    if cache_path is None:
        cache_path = RIBONN_CACHE / "hela_te_values.tsv"

    if not cache_path.exists():
        raise FileNotFoundError(
            f"RiboNN data not found at {cache_path}. "
            "Run `python -m calibration.cli download --dataset zheng_ribonn` first."
        )

    records = []
    text = cache_path.read_text(encoding="utf-8")
    # Try TSV first, fallback to CSV
    dialect = "excel-tab" if "\t" in text[:1000] else "excel"
    reader = csv.DictReader(io.StringIO(text), dialect=dialect)

    def _get(row: dict, *keys: str) -> str:
        for k in keys:
            if k in row and row[k]:
                return row[k].strip()
        return ""

    for row in reader:
        tx_id = _get(row, "transcript_id", "gene_name", "gene", "name")
        te_str = _get(row, "te", "translation_efficiency", "TE", "log2_te")
        utr = _get(row, "utr_sequence", "five_prime_utr", "5utr", "utr")
        cds = _get(row, "cds_sequence", "cds", "coding_sequence")

        if not te_str:
            continue

        try:
            te_val = float(te_str)
        except ValueError:
            continue

        utr_rna = _dna_to_rna(utr) if utr else ""
        cds_rna = _dna_to_rna(cds) if cds else ""

        # Strip stop codon
        if len(cds_rna) >= 3 and cds_rna[-3:] in ("UAA", "UAG", "UGA"):
            cds_rna = cds_rna[:-3]

        if not utr_rna and not cds_rna:
            continue

        records.append(BenchmarkRecord(
            record_id=tx_id or f"ribonn_{len(records)}",
            utr_sequence=utr_rna,
            cds_sequence=cds_rna,
            measured_te=te_val,
            expression_system="mammalian",
            has_full_cds=bool(cds_rna),
        ))

    logger.info("Loaded %d RiboNN records", len(records))

    return BenchmarkDataset(
        name="zheng_ribonn",
        records=records,
        te_metric_name="log2_te",
        description="RiboNN: transcriptome-wide TE atlas (Zheng 2025)",
        citation="Zheng et al. (2025) RiboNN",
    )


def join_with_gencode(te_path: Path, gencode_fasta_path: Path,
                       output_path: Path | None = None) -> Path:
    """Join TE values with GENCODE transcript sequences.

    This reads the raw TE TSV (transcript_id + te) and the GENCODE FASTA,
    extracts UTR and CDS regions, and writes a combined TSV.
    """
    if output_path is None:
        output_path = RIBONN_CACHE / "hela_te_with_sequences.tsv"

    # Parse GENCODE FASTA
    from Bio import SeqIO

    logger.info("Parsing GENCODE FASTA (this may take a few minutes)...")
    tx_seqs = {}
    open_fn = gzip.open if str(gencode_fasta_path).endswith(".gz") else open
    with open_fn(gencode_fasta_path, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            # GENCODE header format: >ENST...|ENSG...|...
            tx_id = record.id.split("|")[0]
            tx_seqs[tx_id] = str(record.seq).upper()

    logger.info("Loaded %d GENCODE transcripts", len(tx_seqs))

    # Read TE values and join
    te_text = te_path.read_text(encoding="utf-8")
    dialect = "excel-tab" if "\t" in te_text[:1000] else "excel"
    reader = csv.DictReader(io.StringIO(te_text), dialect=dialect)

    rows = []
    for row in reader:
        tx_id = row.get("transcript_id", "").strip()
        te_str = row.get("te", row.get("TE", "")).strip()
        if not tx_id or not te_str:
            continue

        # Try exact match, then try without version
        seq = tx_seqs.get(tx_id) or tx_seqs.get(tx_id.split(".")[0])
        if not seq:
            continue

        # For now, use the full sequence as CDS approximation
        # A proper implementation would use GTF annotations for UTR/CDS boundaries
        rows.append({
            "transcript_id": tx_id,
            "te": te_str,
            "utr_sequence": "",  # Requires GTF parsing
            "cds_sequence": seq,
        })

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["transcript_id", "te", "utr_sequence", "cds_sequence"],
                                 delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    logger.info("Wrote %d joined records to %s", len(rows), output_path)
    return output_path
