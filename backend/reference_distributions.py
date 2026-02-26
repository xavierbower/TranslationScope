"""Approximate percentile context vs. human transcriptome."""

from __future__ import annotations

DISTRIBUTIONS: dict[str, dict[int, float]] = {
    "cos": {10: -0.06, 25: -0.03, 50: 0.00, 75: 0.03, 90: 0.06},
    "utr_length": {10: 30, 25: 60, 50: 130, 75: 280, 90: 600},
    "au_content": {10: 48, 25: 53, 50: 59, 75: 65, 90: 72},
    "aug_accessibility": {10: 0.30, 25: 0.45, 50: 0.60, 75: 0.75, 90: 0.87},
    "cds_length": {10: 300, 25: 600, 50: 1200, 75: 2400, 90: 4800},
}


def get_percentile(metric_name: str, value: float) -> float | None:
    """Return approximate percentile (0-100) using linear interpolation."""
    dist = DISTRIBUTIONS.get(metric_name)
    if dist is None:
        return None

    breakpoints = sorted(dist.items())
    pcts = [p for p, _ in breakpoints]
    vals = [v for _, v in breakpoints]

    if value <= vals[0]:
        return max(0, pcts[0] * value / vals[0]) if vals[0] != 0 else pcts[0]
    if value > vals[-1]:
        return min(100, pcts[-1] + (100 - pcts[-1]) * 0.5)

    for i in range(len(vals) - 1):
        if vals[i] <= value <= vals[i + 1]:
            frac = (value - vals[i]) / (vals[i + 1] - vals[i])
            return round(pcts[i] + frac * (pcts[i + 1] - pcts[i]), 1)

    return 50.0
