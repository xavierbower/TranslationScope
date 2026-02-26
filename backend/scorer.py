"""Composite Translation Efficiency Score calculator."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ScoreBreakdown:
    aug_accessibility_score: int = 0
    codon_optimality_score: int = 0
    nterm_codon_score: int = 0
    uorf_score: int = 0
    kozak_score: int = 0
    utr_composition_score: int = 0
    utr_structure_score: int = 0
    total_score: int = 0


@dataclass
class ScorerResult:
    scores: ScoreBreakdown = field(default_factory=ScoreBreakdown)
    rating: str = ""
    rating_color: str = ""
    primary_bottleneck: str = ""
    summary: str = ""
    prioritized_fixes: list[str] = field(default_factory=list)


def _score_aug(accessibility: float) -> int:
    if accessibility >= 0.85:
        return 25
    if accessibility >= 0.70:
        return 19
    if accessibility >= 0.50:
        return 11
    if accessibility >= 0.30:
        return 5
    return 0


def _score_cos(cos: float) -> int:
    if cos >= 0.05:
        return 20
    if cos >= 0.02:
        return 15
    if cos >= -0.02:
        return 10
    if cos >= -0.05:
        return 5
    return 0


def _score_nterm(nterm_cos: float) -> int:
    if nterm_cos > 0:
        return 10
    if nterm_cos >= -0.1:
        return 7
    if nterm_cos >= -0.2:
        return 3
    return 0


def _score_uorf(total: int, high_impact: int) -> int:
    if total == 0:
        return 15
    if total <= 2 and high_impact == 0:
        return 11
    if total <= 2:
        return 5
    if high_impact == 0:
        return 7
    return 0


def _score_kozak(kozak_score: int) -> int:
    if kozak_score >= 6:
        return 10
    if kozak_score >= 4:
        return 6
    return 0


def _score_utr_composition(au_content: float, gg_frequency: float, utr_length: int) -> int:
    if au_content >= 65:
        score = 10
    elif au_content >= 55:
        score = 7
    elif au_content >= 45:
        score = 4
    else:
        score = 0
    if gg_frequency > 5:
        score -= 3
    if utr_length > 500:
        score -= 2
    return max(0, score)


def _score_utr_structure(mfe: float) -> int:
    if mfe > -5:
        return 10
    if mfe >= -15:
        return 7
    if mfe >= -30:
        return 4
    if mfe >= -50:
        return 2
    return 0


def _rating(score: int) -> tuple[str, str]:
    if score >= 85:
        return "High Efficiency", "green"
    if score >= 65:
        return "Moderate Efficiency", "yellow"
    if score >= 40:
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
) -> ScorerResult:
    """Compute composite score with IRES/truncation adjustments."""
    result = ScorerResult()
    s = result.scores

    # Compute raw subscores
    s.aug_accessibility_score = _score_aug(aug_accessibility)
    s.codon_optimality_score = _score_cos(cos)
    s.nterm_codon_score = _score_nterm(nterm_cos)
    s.uorf_score = _score_uorf(total_uaugs, high_impact_uaugs)
    s.kozak_score = _score_kozak(kozak_score_val)
    s.utr_composition_score = _score_utr_composition(au_content, gg_frequency, utr_length)
    s.utr_structure_score = _score_utr_structure(utr_mfe)

    if ires_detected:
        # Suppress UTR-dependent scores, redistribute to codon-based
        suppressed = s.aug_accessibility_score + s.uorf_score + s.utr_structure_score
        s.aug_accessibility_score = 0
        s.uorf_score = 0
        s.utr_structure_score = 0
        # Scale codon scores proportionally
        codon_total = s.codon_optimality_score + s.nterm_codon_score
        max_codon = 30  # 20 + 10
        if codon_total > 0:
            scale = (codon_total + suppressed) / codon_total
            s.codon_optimality_score = min(int(s.codon_optimality_score * scale), 45)
            s.nterm_codon_score = min(int(s.nterm_codon_score * scale), 20)

    if aug_skipped and not ires_detected:
        # Redistribute AUG weight
        lost = s.aug_accessibility_score
        s.aug_accessibility_score = 0
        remaining_total = (s.codon_optimality_score + s.nterm_codon_score +
                           s.uorf_score + s.kozak_score +
                           s.utr_composition_score + s.utr_structure_score)
        if remaining_total > 0:
            scale = (remaining_total + lost) / remaining_total
            s.codon_optimality_score = min(int(s.codon_optimality_score * scale), 30)
            s.nterm_codon_score = min(int(s.nterm_codon_score * scale), 15)
            s.uorf_score = min(int(s.uorf_score * scale), 20)

    s.total_score = (s.aug_accessibility_score + s.codon_optimality_score +
                     s.nterm_codon_score + s.uorf_score + s.kozak_score +
                     s.utr_composition_score + s.utr_structure_score)

    result.rating, result.rating_color = _rating(s.total_score)

    # Find primary bottleneck
    subscales = {
        "AUG Accessibility": (s.aug_accessibility_score, 25),
        "Codon Optimality": (s.codon_optimality_score, 20),
        "N-terminal Codons": (s.nterm_codon_score, 10),
        "uORF Burden": (s.uorf_score, 15),
        "Kozak Context": (s.kozak_score, 10),
        "5' UTR Composition": (s.utr_composition_score, 10),
        "5' UTR Structure": (s.utr_structure_score, 10),
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


def _generate_summary(bottleneck: str, subscales: dict, total: int) -> str:
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
