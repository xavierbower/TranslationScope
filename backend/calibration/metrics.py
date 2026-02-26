"""Evaluation metrics for calibration benchmarking."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from scipy import stats
import numpy as np


@dataclass
class EvalMetrics:
    """Evaluation metrics comparing predicted scores vs. measured TE."""
    spearman_rho: float = 0.0
    spearman_pvalue: float = 1.0
    pearson_r: float = 0.0
    pearson_pvalue: float = 1.0
    precision_top_quartile: float = 0.0
    precision_bottom_quartile: float = 0.0
    ndcg_at_10: float = 0.0
    mean_abs_rank_error: float = 0.0
    n_samples: int = 0


@dataclass
class SubscaleDiagnostics:
    """Per-feature Spearman correlation with measured TE."""
    feature_correlations: dict[str, float] = field(default_factory=dict)
    feature_pvalues: dict[str, float] = field(default_factory=dict)


def spearman(predicted: list[float], measured: list[float]) -> tuple[float, float]:
    """Spearman rank correlation."""
    if len(predicted) < 3:
        return 0.0, 1.0
    rho, pval = stats.spearmanr(predicted, measured)
    return float(rho), float(pval)


def pearson(predicted: list[float], measured: list[float]) -> tuple[float, float]:
    """Pearson correlation."""
    if len(predicted) < 3:
        return 0.0, 1.0
    r, pval = stats.pearsonr(predicted, measured)
    return float(r), float(pval)


def precision_at_quartile(predicted: list[float], measured: list[float],
                           quartile: str = "top") -> float:
    """Fraction of predicted top/bottom quartile that is in actual top/bottom quartile.

    Args:
        quartile: "top" for top 25%, "bottom" for bottom 25%
    """
    n = len(predicted)
    if n < 4:
        return 0.0

    k = n // 4
    pred_arr = np.array(predicted)
    meas_arr = np.array(measured)

    if quartile == "top":
        pred_set = set(np.argsort(pred_arr)[-k:])
        meas_set = set(np.argsort(meas_arr)[-k:])
    else:
        pred_set = set(np.argsort(pred_arr)[:k])
        meas_set = set(np.argsort(meas_arr)[:k])

    overlap = len(pred_set & meas_set)
    return overlap / k if k > 0 else 0.0


def ndcg_at_k(predicted: list[float], measured: list[float], k: int = 10) -> float:
    """Normalized Discounted Cumulative Gain at k.

    Measures how well the top-k predicted items match the actual ranking.
    """
    n = len(predicted)
    if n < k:
        k = n
    if k == 0:
        return 0.0

    pred_arr = np.array(predicted)
    meas_arr = np.array(measured)

    # Rank by predicted score (descending)
    pred_order = np.argsort(pred_arr)[::-1][:k]

    # Relevance: use measured TE as relevance score (normalized to 0-1)
    meas_min = meas_arr.min()
    meas_range = meas_arr.max() - meas_min
    if meas_range == 0:
        return 1.0
    relevance = (meas_arr - meas_min) / meas_range

    # DCG
    dcg = sum(relevance[pred_order[i]] / math.log2(i + 2) for i in range(k))

    # Ideal DCG (sorted by actual relevance)
    ideal_order = np.argsort(relevance)[::-1][:k]
    idcg = sum(relevance[ideal_order[i]] / math.log2(i + 2) for i in range(k))

    return dcg / idcg if idcg > 0 else 0.0


def mean_absolute_rank_error(predicted: list[float], measured: list[float]) -> float:
    """Mean absolute difference between predicted and actual ranks."""
    n = len(predicted)
    if n < 2:
        return 0.0

    pred_ranks = stats.rankdata(predicted)
    meas_ranks = stats.rankdata(measured)
    return float(np.mean(np.abs(pred_ranks - meas_ranks)))


def evaluate(predicted: list[float], measured: list[float]) -> EvalMetrics:
    """Compute all evaluation metrics."""
    m = EvalMetrics(n_samples=len(predicted))

    m.spearman_rho, m.spearman_pvalue = spearman(predicted, measured)
    m.pearson_r, m.pearson_pvalue = pearson(predicted, measured)
    m.precision_top_quartile = precision_at_quartile(predicted, measured, "top")
    m.precision_bottom_quartile = precision_at_quartile(predicted, measured, "bottom")
    m.ndcg_at_10 = ndcg_at_k(predicted, measured, 10)
    m.mean_abs_rank_error = mean_absolute_rank_error(predicted, measured)

    return m


def subscale_diagnostics(
    feature_vectors: list[dict[str, float]],
    measured: list[float],
) -> SubscaleDiagnostics:
    """Compute Spearman of each raw feature vs. measured TE.

    Args:
        feature_vectors: List of dicts with feature name -> value
        measured: Corresponding measured TE values
    """
    diag = SubscaleDiagnostics()

    if not feature_vectors:
        return diag

    feature_names = list(feature_vectors[0].keys())

    for fname in feature_names:
        values = [fv.get(fname, 0.0) for fv in feature_vectors]
        # Skip non-numeric features
        try:
            numeric_vals = [float(v) for v in values]
        except (TypeError, ValueError):
            continue

        rho, pval = spearman(numeric_vals, measured)
        diag.feature_correlations[fname] = rho
        diag.feature_pvalues[fname] = pval

    return diag


def format_metrics(metrics: EvalMetrics, label: str = "") -> str:
    """Format metrics as a readable string."""
    lines = []
    if label:
        lines.append(f"=== {label} ===")
    lines.append(f"  Spearman rho:        {metrics.spearman_rho:+.4f}  (p={metrics.spearman_pvalue:.2e})")
    lines.append(f"  Pearson r:           {metrics.pearson_r:+.4f}  (p={metrics.pearson_pvalue:.2e})")
    lines.append(f"  Precision@top25%:    {metrics.precision_top_quartile:.3f}")
    lines.append(f"  Precision@bottom25%: {metrics.precision_bottom_quartile:.3f}")
    lines.append(f"  NDCG@10:             {metrics.ndcg_at_10:.3f}")
    lines.append(f"  Mean abs rank error: {metrics.mean_abs_rank_error:.1f}")
    lines.append(f"  N samples:           {metrics.n_samples}")
    return "\n".join(lines)
