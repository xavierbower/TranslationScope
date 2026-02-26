"""Programmatic scorer — bypasses FastAPI, takes raw RNA strings."""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Callable

from calibration.config import FOLD_CACHE_DIR
from calibration.parameters import ParameterSet, DEFAULT_PARAMS
from calibration.datasets.base import BenchmarkRecord

logger = logging.getLogger(__name__)


@dataclass
class FeatureVector:
    """Raw feature values extracted from an mRNA construct."""
    aug_accessibility: float = 0.5
    cos: float = 0.0
    nterm_cos: float = 0.0
    total_uaugs: int = 0
    high_impact_uaugs: int = 0
    kozak_score: int = 0
    au_content: float = 50.0
    gg_frequency: float = 0.0
    utr_length: int = 0
    utr_mfe: float = 0.0
    ires_detected: bool = False
    aug_skipped: bool = False


@dataclass
class ScoringOutput:
    """Result of scoring a single sequence."""
    record_id: str = ""
    features: FeatureVector = field(default_factory=FeatureVector)
    subscores: dict[str, float] = field(default_factory=dict)
    total_score: float = 0.0
    measured_te: float = 0.0


# ---- Fold cache (disk-backed, keyed by sequence SHA256) ----

def _seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()[:16]


def _cache_get(key: str) -> dict | None:
    path = FOLD_CACHE_DIR / f"{key}.json"
    if path.exists():
        return json.loads(path.read_text())
    return None


def _cache_set(key: str, value: dict) -> None:
    path = FOLD_CACHE_DIR / f"{key}.json"
    path.write_text(json.dumps(value))


def _fold_utr_cached(utr_sequence: str):
    """Fold UTR with disk cache."""
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from structure_analyzer import fold_utr

    key = f"utr_{_seq_hash(utr_sequence)}"
    cached = _cache_get(key)
    if cached is not None:
        from structure_analyzer import UtrFoldResult
        return UtrFoldResult(**cached)

    result = fold_utr(utr_sequence)
    _cache_set(key, {"mfe": result.mfe, "structure_dotbracket": result.structure_dotbracket,
                      "skipped": result.skipped, "skip_reason": result.skip_reason})
    return result


def _aug_accessibility_cached(utr_sequence: str, cds_sequence: str):
    """Compute AUG accessibility with disk cache."""
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from structure_analyzer import compute_aug_accessibility

    key = f"aug_{_seq_hash(utr_sequence + '|' + cds_sequence[:20])}"
    cached = _cache_get(key)
    if cached is not None:
        from structure_analyzer import AugAccessibilityResult
        return AugAccessibilityResult(**cached)

    result = compute_aug_accessibility(utr_sequence, cds_sequence)
    _cache_set(key, {
        "mfe": result.mfe, "aug_accessibility": result.aug_accessibility,
        "window_sequence": result.window_sequence,
        "structure_dotbracket": result.structure_dotbracket,
        "window_truncated": result.window_truncated,
        "skipped": result.skipped, "skip_reason": result.skip_reason,
    })
    return result


def extract_features(utr: str, cds: str, system: str = "mammalian",
                      skip_structure: bool = False) -> FeatureVector:
    """Extract all features from raw RNA sequences.

    This calls the same analysis pipeline as main.py but without FastAPI.
    """
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

    from codon_analyzer import (
        analyze_codons, analyze_uorfs, analyze_kozak,
        analyze_utr_composition,
    )

    fv = FeatureVector()

    # Codon analysis (fast)
    if cds:
        codons = analyze_codons(cds, system)
        fv.cos = codons.cos
        fv.nterm_cos = codons.nterm_cos

    # uORF analysis (fast)
    if utr:
        uorfs = analyze_uorfs(utr)
        fv.total_uaugs = uorfs.total_uaugs
        fv.high_impact_uaugs = uorfs.high_impact_uaugs

    # Kozak (fast)
    if utr and cds:
        kozak = analyze_kozak(utr, cds)
        fv.kozak_score = kozak.kozak_score

    # UTR composition (fast)
    if utr:
        utr_comp = analyze_utr_composition(utr)
        fv.au_content = utr_comp.au_content
        fv.gg_frequency = utr_comp.gg_frequency
        fv.utr_length = utr_comp.utr_length

    # Structure analysis (slow — uses ViennaRNA)
    if not skip_structure and utr:
        utr_fold = _fold_utr_cached(utr)
        fv.utr_mfe = utr_fold.mfe if not utr_fold.skipped else 0.0

        if cds:
            aug_result = _aug_accessibility_cached(utr, cds)
            fv.aug_accessibility = aug_result.aug_accessibility if not aug_result.skipped else 0.5
            fv.aug_skipped = aug_result.skipped
    else:
        fv.aug_skipped = True

    return fv


def score_sequence(utr: str, cds: str, system: str = "mammalian",
                    params: ParameterSet | None = None,
                    skip_structure: bool = False) -> ScoringOutput:
    """Score a single mRNA construct.

    Args:
        utr: 5' UTR sequence (RNA alphabet)
        cds: CDS sequence (RNA, includes AUG, no stop)
        system: Expression system name
        params: Parameter set (None for defaults)
        skip_structure: Skip ViennaRNA folding (much faster)

    Returns:
        ScoringOutput with features, subscores, and total score.
    """
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from scorer import compute_score
    from codon_analyzer import analyze_codons

    p = params if params is not None else DEFAULT_PARAMS
    fv = extract_features(utr, cds, system, skip_structure)

    # Get nterm_codons for fix generation (not needed for scoring, but compute_score expects it)
    nterm_codons = []
    if cds:
        codons = analyze_codons(cds, system)
        nterm_codons = codons.nterm_codons

    result = compute_score(
        aug_accessibility=fv.aug_accessibility,
        cos=fv.cos,
        nterm_cos=fv.nterm_cos,
        nterm_codons=nterm_codons,
        total_uaugs=fv.total_uaugs,
        high_impact_uaugs=fv.high_impact_uaugs,
        kozak_score_val=fv.kozak_score,
        au_content=fv.au_content,
        gg_frequency=fv.gg_frequency,
        utr_length=fv.utr_length,
        utr_mfe=fv.utr_mfe,
        ires_detected=fv.ires_detected,
        aug_skipped=fv.aug_skipped,
        params=p,
    )

    return ScoringOutput(
        features=fv,
        subscores={
            "aug_accessibility": result.scores.aug_accessibility_score,
            "codon_optimality": result.scores.codon_optimality_score,
            "nterm_codons": result.scores.nterm_codon_score,
            "uorf_burden": result.scores.uorf_score,
            "kozak_context": result.scores.kozak_score,
            "utr_composition": result.scores.utr_composition_score,
            "utr_structure": result.scores.utr_structure_score,
        },
        total_score=result.scores.total_score,
    )


def score_from_features(fv: FeatureVector, params: ParameterSet | None = None) -> float:
    """Score from pre-extracted features (fast — no ViennaRNA calls).

    Used during optimization when features are pre-computed.
    """
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from scorer import compute_score

    p = params if params is not None else DEFAULT_PARAMS

    result = compute_score(
        aug_accessibility=fv.aug_accessibility,
        cos=fv.cos,
        nterm_cos=fv.nterm_cos,
        nterm_codons=[],
        total_uaugs=fv.total_uaugs,
        high_impact_uaugs=fv.high_impact_uaugs,
        kozak_score_val=fv.kozak_score,
        au_content=fv.au_content,
        gg_frequency=fv.gg_frequency,
        utr_length=fv.utr_length,
        utr_mfe=fv.utr_mfe,
        ires_detected=fv.ires_detected,
        aug_skipped=fv.aug_skipped,
        params=p,
    )
    return result.scores.total_score


def score_dataset(
    records: list[BenchmarkRecord],
    params: ParameterSet | None = None,
    skip_structure: bool = False,
    progress_callback: Callable[[int, int], None] | None = None,
) -> list[ScoringOutput]:
    """Score a batch of benchmark records.

    Args:
        records: List of BenchmarkRecords
        params: Parameter set
        skip_structure: Skip ViennaRNA folding
        progress_callback: Called with (current, total) after each record

    Returns:
        List of ScoringOutput, one per record.
    """
    outputs = []
    total = len(records)

    for i, rec in enumerate(records):
        out = score_sequence(
            utr=rec.utr_sequence,
            cds=rec.cds_sequence,
            system=rec.expression_system,
            params=params,
            skip_structure=skip_structure,
        )
        out.record_id = rec.record_id
        out.measured_te = rec.measured_te
        outputs.append(out)

        if progress_callback:
            progress_callback(i + 1, total)

    return outputs


def extract_features_batch(
    records: list[BenchmarkRecord],
    skip_structure: bool = False,
    progress_callback: Callable[[int, int], None] | None = None,
) -> list[tuple[BenchmarkRecord, FeatureVector]]:
    """Extract features for a batch of records (pre-computation for optimization).

    Returns list of (record, feature_vector) tuples.
    """
    results = []
    total = len(records)

    for i, rec in enumerate(records):
        fv = extract_features(
            utr=rec.utr_sequence,
            cds=rec.cds_sequence,
            system=rec.expression_system,
            skip_structure=skip_structure,
        )
        results.append((rec, fv))

        if progress_callback:
            progress_callback(i + 1, total)

    return results
