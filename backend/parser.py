"""Parse SnapGene .dna and GenBank .gbk files and extract CDS/UTR sequences."""

from __future__ import annotations

import io
import tempfile
from dataclasses import dataclass, field
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from snapgene_reader import snapgene_file_to_seqrecord

MAX_FILE_SIZE = 10 * 1024 * 1024  # 10 MB
MAX_UTR_INFERRED = 300

RESISTANCE_MARKERS = {
    "ampr", "bla", "kanr", "neor", "hygr", "puror", "zeor",
    "cmr", "tetr", "spcr", "aph", "cat", "dhfr",
    "resistance", "selection",
}
REPORTER_MARKERS = {
    "gfp", "rfp", "mcherry", "egfp", "cfp", "yfp",
    "luciferase", "lacz", "laci",
}
BACKBONE_MARKERS = {"rep", "ori"}
EXCLUDE_MARKERS = RESISTANCE_MARKERS | REPORTER_MARKERS | BACKBONE_MARKERS

STOP_CODONS = {"UAA", "UAG", "UGA"}
VALID_RNA = set("AUGC")


def _sanitize_rna(rna: str) -> tuple[str, bool]:
    """Remove non-AUGC characters from RNA. Returns (cleaned, had_ambiguous)."""
    if all(c in VALID_RNA for c in rna):
        return rna, False
    cleaned = "".join(c for c in rna if c in VALID_RNA)
    return cleaned, True


@dataclass
class CdsCandidate:
    index: int
    label: str
    length_nt: int
    length_aa: int
    start: int
    end: int
    strand: int
    filtered_out: bool = False


@dataclass
class ParseResult:
    candidates: list[CdsCandidate] = field(default_factory=list)
    auto_selected: int | None = None
    filtered_candidates: list[CdsCandidate] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    gene_name: str = ""


@dataclass
class ExtractionResult:
    gene_name: str = ""
    utr_sequence: str = ""
    cds_sequence: str = ""
    full_sequence: str = ""
    aug_position: int = 0
    is_circular: bool = False
    utr_source: str = "inferred"
    strand: int = 1
    warnings: list[str] = field(default_factory=list)


def _label(feature) -> str:
    """Extract a human-readable label from a SeqFeature."""
    for key in ("label", "gene", "product", "note"):
        vals = feature.qualifiers.get(key, [])
        if isinstance(vals, str):
            return vals
        if vals:
            return vals[0]
    return "unnamed"


def _should_exclude(label: str) -> bool:
    low = label.lower()
    return any(marker in low for marker in EXCLUDE_MARKERS)


def parse_file(file_bytes: bytes, filename: str = "") -> tuple[Any, list[Any], bool, str]:
    """Return (seqrecord, features, is_circular, full_sequence)."""
    if len(file_bytes) > MAX_FILE_SIZE:
        raise ValueError("File exceeds 10 MB limit.")

    ext = filename.lower().rsplit(".", 1)[-1] if filename else ""

    if ext in ("gbk", "gb", "genbank"):
        try:
            text = file_bytes.decode("utf-8", errors="replace")
            records = list(SeqIO.parse(io.StringIO(text), "genbank"))
            if not records:
                raise ValueError("No records found in GenBank file.")
            record = records[0]
        except ValueError:
            raise
        except Exception as exc:
            raise ValueError(f"Failed to parse GenBank file: {exc}") from exc
    else:
        # SnapGene .dna — requires a file path, not a file object
        try:
            with tempfile.NamedTemporaryFile(suffix=".dna", delete=False) as tmp:
                tmp.write(file_bytes)
                tmp.flush()
                record = snapgene_file_to_seqrecord(tmp.name)
        except Exception as exc:
            raise ValueError(f"Failed to parse .dna file: {exc}") from exc

    is_circular = record.annotations.get("topology", "").lower() == "circular"
    full_seq = str(record.seq).upper()
    return record, list(record.features), is_circular, full_seq


def _effective_type(feat) -> str:
    """Get the effective feature type, handling pLannotate-style GBK files."""
    ft = feat.type
    if ft in ("CDS", "gene"):
        return ft
    # pLannotate stores real type in the 'other' qualifier
    other = feat.qualifiers.get("other", [])
    if isinstance(other, list):
        for val in other:
            if val in ("CDS", "gene"):
                return val
    elif isinstance(other, str) and other in ("CDS", "gene"):
        return other
    return ft


def select_cds_candidates(features: list[Any]) -> ParseResult:
    """Identify CDS candidates, filtering common markers."""
    result = ParseResult()
    all_cds: list[tuple[int, Any]] = []
    for i, feat in enumerate(features):
        if _effective_type(feat) in ("CDS", "gene"):
            all_cds.append((i, feat))

    if not all_cds:
        raise ValueError("No CDS or gene features found in this file.")

    candidates: list[CdsCandidate] = []
    filtered_out: list[CdsCandidate] = []
    for idx, (orig_i, feat) in enumerate(all_cds):
        lab = _label(feat)
        start = int(feat.location.start)
        end = int(feat.location.end)
        strand = feat.location.strand or 1
        length_nt = end - start
        c = CdsCandidate(
            index=orig_i, label=lab, length_nt=length_nt,
            length_aa=length_nt // 3, start=start, end=end, strand=strand,
        )
        if _should_exclude(lab):
            c.filtered_out = True
            filtered_out.append(c)
        else:
            candidates.append(c)

    if len(candidates) == 1:
        result.candidates = candidates
        result.auto_selected = 0
        result.gene_name = candidates[0].label
    elif len(candidates) > 1:
        result.candidates = candidates
    else:
        # All were filtered – return everything with a warning
        for c in filtered_out:
            c.filtered_out = False
        result.candidates = filtered_out
        filtered_out = []
        result.warnings.append(
            "Auto-filtering removed all CDS features. Showing all candidates — please select manually."
        )

    result.filtered_candidates = filtered_out
    return result


def _extract_feature_dna(feat, full_seq: str) -> str:
    """Extract DNA for a feature, handling compound/join locations."""
    if isinstance(feat.location, CompoundLocation):
        # Multi-exon: concatenate each part in order
        parts = []
        for part in feat.location.parts:
            segment = full_seq[int(part.start):int(part.end)]
            if part.strand == -1:
                segment = str(Seq(segment).reverse_complement())
            parts.append(segment)
        cds_dna = "".join(parts)
        # For compound locations on the minus strand, the parts are already
        # individually reverse-complemented; if the overall feature is minus-strand
        # and the parts weren't individually flagged, reverse the whole thing
        if feat.location.strand == -1 and all(p.strand != -1 for p in feat.location.parts):
            cds_dna = str(Seq(cds_dna).reverse_complement())
        return cds_dna
    else:
        start = int(feat.location.start)
        end = int(feat.location.end)
        cds_dna = full_seq[start:end]
        if (feat.location.strand or 1) == -1:
            cds_dna = str(Seq(cds_dna).reverse_complement())
        return cds_dna


def extract_sequences(
    file_bytes: bytes,
    selected_feature_index: int,
    filename: str = "",
) -> ExtractionResult:
    """Extract UTR and CDS sequences for the selected CDS feature."""
    record, features, is_circular, full_seq = parse_file(file_bytes, filename)
    seq_len = len(full_seq)

    feat = features[selected_feature_index]
    lab = _label(feat)
    start = int(feat.location.start)
    end = int(feat.location.end)
    strand = feat.location.strand or 1

    res = ExtractionResult(gene_name=lab, is_circular=is_circular, strand=strand)

    # Extract CDS DNA (handles compound/join locations)
    cds_dna = _extract_feature_dna(feat, full_seq)

    # Transcribe to RNA and sanitize
    cds_rna = cds_dna.replace("T", "U").upper()
    cds_rna, had_ambiguous = _sanitize_rna(cds_rna)
    if had_ambiguous:
        res.warnings.append(
            "Ambiguous bases (N, R, Y, etc.) detected in CDS and removed. "
            "Codon analysis may be slightly affected."
        )
    if len(cds_rna) % 3 != 0:
        res.warnings.append(
            "CDS length is not divisible by 3. Possible annotation error. Codon analysis may be unreliable."
        )

    # Check/strip stop codon
    if len(cds_rna) >= 3:
        last_codon = cds_rna[-3:]
        if last_codon in STOP_CODONS:
            cds_rna = cds_rna[:-3]
        else:
            res.warnings.append("No stop codon detected at the end of the CDS.")

    res.cds_sequence = cds_rna

    # --- UTR extraction ---
    # Determine upstream region relative to the strand
    if strand == 1:
        cds_start_pos = start
    else:
        cds_start_pos = end  # for minus strand, "upstream" is after end in genome coords

    # Look for annotated 5' UTR
    annotated_utr = None
    for f in features:
        _labs = f.qualifiers.get("label", [])
        _notes = f.qualifiers.get("note", [])
        if isinstance(_labs, str):
            _labs = [_labs]
        if isinstance(_notes, str):
            _notes = [_notes]
        f_label = " ".join(_labs + _notes)
        if "5" in f_label.lower() and "utr" in f_label.lower():
            f_start = int(f.location.start)
            f_end = int(f.location.end)
            if strand == 1 and f_end <= cds_start_pos and cds_start_pos - f_start <= 2000:
                annotated_utr = full_seq[f_start:f_end]
                if f.location.strand == -1:
                    annotated_utr = str(Seq(annotated_utr).reverse_complement())
                break
            elif strand == -1 and f_start >= cds_start_pos and f_end - cds_start_pos <= 2000:
                annotated_utr = full_seq[cds_start_pos:f_end]
                annotated_utr = str(Seq(annotated_utr).reverse_complement())
                break

    if annotated_utr is not None:
        utr_dna = annotated_utr
        res.utr_source = "annotated"
    else:
        # Infer up to 300 nt upstream
        if strand == 1:
            desired = min(MAX_UTR_INFERRED, cds_start_pos)
            if is_circular and cds_start_pos < 500:
                wrap_needed = MAX_UTR_INFERRED - cds_start_pos
                if wrap_needed > 0:
                    utr_dna = full_seq[-(wrap_needed):] + full_seq[:cds_start_pos]
                else:
                    utr_dna = full_seq[cds_start_pos - MAX_UTR_INFERRED:cds_start_pos]
            else:
                utr_dna = full_seq[cds_start_pos - desired:cds_start_pos]
        else:
            desired = min(MAX_UTR_INFERRED, seq_len - cds_start_pos)
            utr_dna = full_seq[cds_start_pos:cds_start_pos + desired]
            utr_dna = str(Seq(utr_dna).reverse_complement())

        res.utr_source = "inferred"
        res.warnings.append(
            "No 5' UTR feature annotated. Using up to 300 nt upstream of CDS start. "
            "This may include promoter or backbone sequence — UTR-based scores should be interpreted cautiously."
        )

    utr_rna = utr_dna.replace("T", "U").upper()
    utr_rna, utr_ambiguous = _sanitize_rna(utr_rna)
    if utr_ambiguous and not had_ambiguous:
        res.warnings.append(
            "Ambiguous bases detected in UTR region and removed."
        )
    res.utr_sequence = utr_rna
    res.aug_position = len(utr_rna)
    res.full_sequence = utr_rna + cds_rna

    return res
