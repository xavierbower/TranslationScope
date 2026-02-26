"""Composite Translation Efficiency Score calculator."""

from __future__ import annotations

from dataclasses import dataclass, field

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from calibration.parameters import ParameterSet, DEFAULT_PARAMS


@dataclass
class ScoreBreakdown:
    aug_accessibility_score: float = 0.0
    codon_optimality_score: float = 0.0
    nterm_codon_score: float = 0.0
    uorf_score: float = 0.0
    kozak_score: float = 0.0
    utr_composition_score: float = 0.0
    utr_structure_score: float = 0.0
    total_score: float = 0.0


@dataclass
class ScorerResult:
    scores: ScoreBreakdown = field(default_factory=ScoreBreakdown)
    rating: str = ""
    rating_color: str = ""
    primary_bottleneck: str = ""
    summary: str = ""
    prioritized_fixes: list[str] = field(default_factory=list)


def _score_aug(accessibility: float, p: ParameterSet) -> float:
    if accessibility >= p.aug_bp1:
        return p.w_aug * p.aug_frac1
    if accessibility >= p.aug_bp2:
        return p.w_aug * p.aug_frac2
    if accessibility >= p.aug_bp3:
        return p.w_aug * p.aug_frac3
    if accessibility >= p.aug_bp4:
        return p.w_aug * p.aug_frac4
    return 0.0


def _score_cos(cos: float, p: ParameterSet) -> float:
    if cos >= p.cos_bp1:
        return p.w_cos * p.cos_frac1
    if cos >= p.cos_bp2:
        return p.w_cos * p.cos_frac2
    if cos >= p.cos_bp3:
        return p.w_cos * p.cos_frac3
    if cos >= p.cos_bp4:
        return p.w_cos * p.cos_frac4
    return 0.0


def _score_nterm(nterm_cos: float, p: ParameterSet) -> float:
    if nterm_cos > p.nterm_bp1:
        return p.w_nterm * p.nterm_frac1
    if nterm_cos >= p.nterm_bp2:
        return p.w_nterm * p.nterm_frac2
    if nterm_cos >= p.nterm_bp3:
        return p.w_nterm * p.nterm_frac3
    return 0.0


def _score_uorf(total: int, high_impact: int, p: ParameterSet) -> float:
    if total == 0:
        return p.w_uorf * p.uorf_frac_none
    if total <= 2 and high_impact == 0:
        return p.w_uorf * p.uorf_frac_few_low
    if total <= 2:
        return p.w_uorf * p.uorf_frac_few_high
    if high_impact == 0:
        return p.w_uorf * p.uorf_frac_many_low
    return 0.0


def _score_kozak(kozak_score: int, p: ParameterSet) -> float:
    if kozak_score >= p.kozak_bp1:
        return p.w_kozak * p.kozak_frac1
    if kozak_score >= p.kozak_bp2:
        return p.w_kozak * p.kozak_frac2
    return 0.0


def _score_utr_composition(au_content: float, gg_frequency: float,
                            utr_length: int, p: ParameterSet) -> float:
    if au_content >= p.utr_comp_au_bp1:
        score = p.w_utr_comp * p.utr_comp_au_frac1
    elif au_content >= p.utr_comp_au_bp2:
        score = p.w_utr_comp * p.utr_comp_au_frac2
    elif au_content >= p.utr_comp_au_bp3:
        score = p.w_utr_comp * p.utr_comp_au_frac3
    else:
        score = 0.0
    if gg_frequency > p.utr_comp_gg_thresh:
        score -= p.utr_comp_gg_penalty
    if utr_length > p.utr_comp_len_thresh:
        score -= p.utr_comp_len_penalty
    return max(0.0, score)


def _score_utr_structure(mfe: float, p: ParameterSet) -> float:
    if mfe > p.utr_struct_bp1:
        return p.w_utr_struct * p.utr_struct_frac1
    if mfe >= p.utr_struct_bp2:
        return p.w_utr_struct * p.utr_struct_frac2
    if mfe >= p.utr_struct_bp3:
        return p.w_utr_struct * p.utr_struct_frac3
    if mfe >= p.utr_struct_bp4:
        return p.w_utr_struct * p.utr_struct_frac4
    return 0.0


def _rating(score: float, p: ParameterSet) -> tuple[str, str]:
    if score >= p.rating_high:
        return "High Efficiency", "green"
    if score >= p.rating_moderate:
        return "Moderate Efficiency", "yellow"
    if score >= p.rating_low:
        return "Low Efficiency", "orange"
    return "Poor Efficiency", "red"


def compute_score(
    aug_accessibility: float,
    cos: float,
    nterm_cos: float,
    nterm_codons: list,
    total_uaugs: int,
    high_impact_uaugs: int,
    kozak_score_val: int,
    au_content: float,
    gg_frequency: float,
    utr_length: int,
    utr_mfe: float,
    ires_detected: bool = False,
    aug_skipped: bool = False,
    params: ParameterSet | None = None,
) -> ScorerResult:
    """Compute composite score with IRES/truncation adjustments.

    If params is None, uses the default (literature-based) parameter set.
    """
    p = params if params is not None else DEFAULT_PARAMS
    result = ScorerResult()
    s = result.scores

    # Compute raw subscores
    s.aug_accessibility_score = _score_aug(aug_accessibility, p)
    s.codon_optimality_score = _score_cos(cos, p)
    s.nterm_codon_score = _score_nterm(nterm_cos, p)
    s.uorf_score = _score_uorf(total_uaugs, high_impact_uaugs, p)
    s.kozak_score = _score_kozak(kozak_score_val, p)
    s.utr_composition_score = _score_utr_composition(au_content, gg_frequency, utr_length, p)
    s.utr_structure_score = _score_utr_structure(utr_mfe, p)

    if ires_detected:
        # Suppress UTR-dependent scores, redistribute to codon-based
        suppressed = s.aug_accessibility_score + s.uorf_score + s.utr_structure_score
        s.aug_accessibility_score = 0.0
        s.uorf_score = 0.0
        s.utr_structure_score = 0.0
        codon_total = s.codon_optimality_score + s.nterm_codon_score
        if codon_total > 0:
            scale = (codon_total + suppressed) / codon_total
            s.codon_optimality_score = min(s.codon_optimality_score * scale, p.w_cos * 2.25)
            s.nterm_codon_score = min(s.nterm_codon_score * scale, p.w_nterm * 2.0)

    if aug_skipped and not ires_detected:
        lost = s.aug_accessibility_score
        s.aug_accessibility_score = 0.0
        remaining_total = (s.codon_optimality_score + s.nterm_codon_score +
                           s.uorf_score + s.kozak_score +
                           s.utr_composition_score + s.utr_structure_score)
        if remaining_total > 0:
            scale = (remaining_total + lost) / remaining_total
            s.codon_optimality_score = min(s.codon_optimality_score * scale, p.w_cos * 1.5)
            s.nterm_codon_score = min(s.nterm_codon_score * scale, p.w_nterm * 1.5)
            s.uorf_score = min(s.uorf_score * scale, p.w_uorf * (4 / 3))

    s.total_score = (s.aug_accessibility_score + s.codon_optimality_score +
                     s.nterm_codon_score + s.uorf_score + s.kozak_score +
                     s.utr_composition_score + s.utr_structure_score)

    # Round for display (keep float internally for optimization)
    s.total_score = round(s.total_score, 2)

    result.rating, result.rating_color = _rating(s.total_score, p)

    # Find primary bottleneck
    subscales = {
        "AUG Accessibility": (s.aug_accessibility_score, p.w_aug),
        "Codon Optimality": (s.codon_optimality_score, p.w_cos),
        "N-terminal Codons": (s.nterm_codon_score, p.w_nterm),
        "uORF Burden": (s.uorf_score, p.w_uorf),
        "Kozak Context": (s.kozak_score, p.w_kozak),
        "5' UTR Composition": (s.utr_composition_score, p.w_utr_comp),
        "5' UTR Structure": (s.utr_structure_score, p.w_utr_struct),
    }

    if ires_detected:
        for k in ("AUG Accessibility", "uORF Burden", "5' UTR Structure"):
            subscales.pop(k, None)

    # Lowest fraction of max
    worst_name = ""
    worst_frac = 2.0
    for name, (sc, mx) in subscales.items():
        frac = sc / mx if mx > 0 else 1.0
        if frac < worst_frac:
            worst_frac = frac
            worst_name = name
    result.primary_bottleneck = worst_name

    # Generate summary
    result.summary = _generate_summary(worst_name, subscales, s.total_score)
    result.prioritized_fixes = _generate_fixes(subscales, nterm_codons, cos, aug_accessibility, kozak_score_val, ires_detected)

    return result


def _generate_summary(bottleneck: str, subscales: dict, total: float) -> str:
    explanations = {
        "AUG Accessibility": (
            "The start codon region appears to be sequestered in secondary structure, "
            "which can impede ribosome scanning and initiation. "
            "Opening the structure around the AUG may improve translation initiation rates."
        ),
        "Codon Optimality": (
            "The overall codon composition is suboptimal for the selected expression system. "
            "Codons associated with slower elongation are overrepresented, which correlates "
            "with reduced mRNA translation efficiency and stability."
        ),
        "N-terminal Codons": (
            "The first 10 codons after the start methionine contain unfavorable codons. "
            "N-terminal codons have ~2x the per-residue impact on TE compared to mid-CDS codons."
        ),
        "uORF Burden": (
            "Upstream open reading frames in the 5' UTR may compete with the main CDS for "
            "ribosome initiation. High-impact uORFs that overlap the main ORF are particularly disruptive."
        ),
        "Kozak Context": (
            "The sequence context around the start codon deviates from the consensus Kozak sequence. "
            "Strong Kozak context (especially A/G at -3 and G at +4) improves start codon recognition."
        ),
        "5' UTR Composition": (
            "The 5' UTR sequence composition is suboptimal. AU-rich UTRs generally correlate with "
            "higher translation efficiency, while GG dinucleotides and excessive length reduce it."
        ),
        "5' UTR Structure": (
            "The 5' UTR is predicted to form stable secondary structure under in vitro conditions. "
            "While in-cell helicases can resolve some structures, very stable folds may still impede scanning."
        ),
    }
    explanation = explanations.get(bottleneck, "Multiple factors may be limiting translation efficiency.")
    return f"Primary bottleneck: {bottleneck}. {explanation}"


def _generate_fixes(
    subscales: dict,
    nterm_codons: list,
    cos: float,
    aug_accessibility: float,
    kozak_val: int,
    ires_detected: bool,
) -> list[str]:
    fixes = []
    scored = [(name, sc, mx) for name, (sc, mx) in subscales.items()]
    scored.sort(key=lambda x: x[1] / x[2] if x[2] > 0 else 1.0)

    for name, sc, mx in scored[:3]:
        frac = sc / mx if mx > 0 else 1.0
        if frac >= 0.8:
            continue
        if name == "N-terminal Codons" and nterm_codons:
            bad = [f"position {i+2} ({c.codon})" for i, c in enumerate(nterm_codons) if c.category == "negative"]
            if bad:
                fixes.append(
                    f"Optimize N-terminal codons at {', '.join(bad[:3])} — "
                    "these are in the unfavorable category. Synonymous substitutions here "
                    "may have disproportionate impact given the ~2x per-residue importance of early CDS codons."
                )
            else:
                fixes.append("Consider optimizing N-terminal codon context (codons 2-11).")
        elif name == "Codon Optimality":
            fixes.append(
                f"Improve overall codon optimality (current COS: {cos:+.3f}). "
                "Replace negative-category codons with synonymous positive or neutral alternatives, "
                "prioritizing the most frequently used unfavorable codons."
            )
        elif name == "AUG Accessibility":
            fixes.append(
                f"Improve start codon accessibility (current: {aug_accessibility:.0%}). "
                "Consider modifying the sequence immediately upstream or downstream of the AUG "
                "to reduce predicted secondary structure in this region."
            )
        elif name == "Kozak Context":
            fixes.append(
                f"Strengthen Kozak context (current score: {kozak_val}/7). "
                "Ensure A or G at position -3 and G at position +4 relative to the AUG."
            )
        elif name == "uORF Burden":
            fixes.append(
                "Remove or mutate upstream AUG codons in the 5' UTR. "
                "If AUG removal is not possible, ensure any remaining uORFs have in-frame stop codons "
                "well upstream of the main CDS start."
            )
        elif name == "5' UTR Composition":
            fixes.append(
                "Consider redesigning the 5' UTR to increase AU content and reduce GG dinucleotide frequency. "
                "Shorter UTRs (50-100 nt) generally correlate with higher translation efficiency."
            )
        elif name == "5' UTR Structure":
            fixes.append(
                "Reduce predicted secondary structure in the 5' UTR by introducing synonymous mutations "
                "or shortening the UTR. Focus on disrupting stems near the 5' cap."
            )

    return fixes[:3]
