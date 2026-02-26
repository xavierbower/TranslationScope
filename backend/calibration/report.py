"""HTML report generation with Matplotlib plots."""

from __future__ import annotations

import base64
import io
import logging
from dataclasses import asdict
from pathlib import Path

import numpy as np

from calibration.parameters import ParameterSet, DEFAULT_PARAMS
from calibration.metrics import EvalMetrics, SubscaleDiagnostics
from calibration.config import REPORT_DIR

logger = logging.getLogger(__name__)


def _fig_to_base64(fig) -> str:
    """Convert matplotlib figure to base64 PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def scatter_plot(predicted: list[float], measured: list[float],
                  title: str = "Predicted vs. Measured TE") -> str:
    """Create scatter plot of predicted score vs. measured TE."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    ax.scatter(predicted, measured, alpha=0.5, s=20, edgecolors="none")
    ax.set_xlabel("Predicted TranslationScope Score")
    ax.set_ylabel("Measured Translation Efficiency")
    ax.set_title(title)

    # Add trend line
    if len(predicted) > 2:
        z = np.polyfit(predicted, measured, 1)
        p = np.poly1d(z)
        x_line = np.linspace(min(predicted), max(predicted), 100)
        ax.plot(x_line, p(x_line), "r--", alpha=0.7, linewidth=1)

    fig.tight_layout()
    b64 = _fig_to_base64(fig)
    plt.close(fig)
    return b64


def weight_comparison_chart(baseline: ParameterSet, optimized: ParameterSet) -> str:
    """Bar chart comparing baseline vs. optimized weights."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["AUG", "COS", "N-term", "uORF", "Kozak", "UTR Comp", "UTR Struct"]
    base_w = baseline.weights
    opt_w = optimized.weights

    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.bar(x - width / 2, base_w, width, label="Baseline", color="#4A90D9", alpha=0.8)
    ax.bar(x + width / 2, opt_w, width, label="Optimized", color="#E85D75", alpha=0.8)
    ax.set_ylabel("Weight (max points)")
    ax.set_title("Subscale Weight Comparison")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.legend()
    ax.set_ylim(0, max(max(base_w), max(opt_w)) * 1.2)

    fig.tight_layout()
    b64 = _fig_to_base64(fig)
    plt.close(fig)
    return b64


def subscale_heatmap(diagnostics: dict[str, SubscaleDiagnostics]) -> str:
    """Heatmap of per-feature Spearman across all datasets."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Collect all feature names
    all_features = set()
    for diag in diagnostics.values():
        all_features.update(diag.feature_correlations.keys())
    features = sorted(all_features)
    datasets = list(diagnostics.keys())

    if not features or not datasets:
        return ""

    matrix = np.zeros((len(features), len(datasets)))
    for j, ds in enumerate(datasets):
        for i, f in enumerate(features):
            matrix[i, j] = diagnostics[ds].feature_correlations.get(f, 0.0)

    fig, ax = plt.subplots(1, 1, figsize=(max(6, len(datasets) * 2), max(4, len(features) * 0.5)))
    im = ax.imshow(matrix, cmap="RdYlBu", aspect="auto", vmin=-0.5, vmax=0.5)
    ax.set_xticks(range(len(datasets)))
    ax.set_xticklabels(datasets, rotation=30, ha="right")
    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(features)
    ax.set_title("Per-Feature Spearman Correlation with Measured TE")
    fig.colorbar(im, ax=ax, label="Spearman rho")

    # Annotate cells
    for i in range(len(features)):
        for j in range(len(datasets)):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center",
                    fontsize=8, color="black" if abs(matrix[i, j]) < 0.3 else "white")

    fig.tight_layout()
    b64 = _fig_to_base64(fig)
    plt.close(fig)
    return b64


def cv_boxplot(fold_rhos: dict[str, list[float]]) -> str:
    """Box plot of CV fold Spearman rho values."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not fold_rhos:
        return ""

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    labels = list(fold_rhos.keys())
    data = [fold_rhos[l] for l in labels]
    ax.boxplot(data, labels=labels)
    ax.set_ylabel("Spearman rho")
    ax.set_title("Cross-Validation Stability")
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)

    fig.tight_layout()
    b64 = _fig_to_base64(fig)
    plt.close(fig)
    return b64


def generate_html_report(
    baseline_metrics: dict[str, EvalMetrics],
    calibrated_metrics: dict[str, EvalMetrics] | None = None,
    baseline_params: ParameterSet = DEFAULT_PARAMS,
    calibrated_params: ParameterSet | None = None,
    scatter_data: dict[str, tuple[list[float], list[float]]] | None = None,
    diagnostics: dict[str, SubscaleDiagnostics] | None = None,
    fold_rhos: dict[str, list[float]] | None = None,
    output_path: Path | None = None,
) -> Path:
    """Generate full HTML calibration report.

    Args:
        baseline_metrics: Dataset name -> EvalMetrics for baseline params
        calibrated_metrics: Same for calibrated params (optional)
        baseline_params: Baseline parameter set
        calibrated_params: Calibrated parameter set (optional)
        scatter_data: Dataset name -> (predicted, measured) for scatter plots
        diagnostics: Dataset name -> SubscaleDiagnostics for heatmap
        fold_rhos: Label -> list of CV fold rhos for box plot
        output_path: Where to save HTML

    Returns:
        Path to generated HTML report.
    """
    if output_path is None:
        output_path = REPORT_DIR / "calibration_report.html"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    html_parts = [
        "<!DOCTYPE html>",
        "<html><head>",
        "<title>TranslationScope Calibration Report</title>",
        "<style>",
        "body { font-family: -apple-system, BlinkMacSystemFont, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }",
        "h1 { color: #2C3E50; }",
        "h2 { color: #34495E; border-bottom: 2px solid #3498DB; padding-bottom: 5px; }",
        "table { border-collapse: collapse; width: 100%; margin: 10px 0; }",
        "th, td { border: 1px solid #ddd; padding: 8px; text-align: right; }",
        "th { background-color: #3498DB; color: white; }",
        "tr:nth-child(even) { background-color: #f2f2f2; }",
        "td:first-child { text-align: left; font-weight: bold; }",
        ".improved { color: #27AE60; font-weight: bold; }",
        ".degraded { color: #E74C3C; }",
        "img { max-width: 100%; border: 1px solid #ddd; margin: 10px 0; }",
        ".metric-card { display: inline-block; background: #ECF0F1; padding: 15px; margin: 5px; border-radius: 8px; min-width: 150px; text-align: center; }",
        ".metric-card .value { font-size: 24px; font-weight: bold; color: #2C3E50; }",
        ".metric-card .label { font-size: 12px; color: #7F8C8D; }",
        "</style>",
        "</head><body>",
        "<h1>TranslationScope Calibration Report</h1>",
    ]

    # 1. Metrics comparison table
    html_parts.append("<h2>1. Evaluation Metrics</h2>")
    html_parts.append("<table><tr><th>Dataset</th><th>Metric</th>")
    html_parts.append("<th>Baseline</th>")
    if calibrated_metrics:
        html_parts.append("<th>Calibrated</th><th>Change</th>")
    html_parts.append("</tr>")

    for ds_name, base_m in baseline_metrics.items():
        cal_m = calibrated_metrics.get(ds_name) if calibrated_metrics else None
        for metric_name, base_val, cal_val in [
            ("Spearman rho", base_m.spearman_rho, cal_m.spearman_rho if cal_m else None),
            ("Pearson r", base_m.pearson_r, cal_m.pearson_r if cal_m else None),
            ("Precision@top25%", base_m.precision_top_quartile, cal_m.precision_top_quartile if cal_m else None),
            ("NDCG@10", base_m.ndcg_at_10, cal_m.ndcg_at_10 if cal_m else None),
        ]:
            html_parts.append(f"<tr><td>{ds_name}</td><td>{metric_name}</td>")
            html_parts.append(f"<td>{base_val:.4f}</td>")
            if cal_val is not None:
                delta = cal_val - base_val
                css_class = "improved" if delta > 0 else "degraded" if delta < 0 else ""
                html_parts.append(f"<td>{cal_val:.4f}</td>")
                html_parts.append(f'<td class="{css_class}">{delta:+.4f}</td>')
            html_parts.append("</tr>")

    html_parts.append("</table>")

    # 2. Parameter comparison
    if calibrated_params:
        html_parts.append("<h2>2. Parameter Comparison</h2>")
        html_parts.append(weight_comparison_chart(baseline_params, calibrated_params))
        html_parts.append("<table><tr><th>Parameter</th><th>Baseline</th><th>Calibrated</th><th>Change</th></tr>")
        base_d = baseline_params.to_dict()
        cal_d = calibrated_params.to_dict()
        for key in base_d:
            b = base_d[key]
            c = cal_d[key]
            if isinstance(b, (int, float)):
                delta = c - b
                css = "improved" if abs(delta) > 0.01 else ""
                html_parts.append(f'<tr><td>{key}</td><td>{b:.4f}</td><td>{c:.4f}</td>')
                html_parts.append(f'<td class="{css}">{delta:+.4f}</td></tr>')
        html_parts.append("</table>")

        # Weight chart as image
        chart_b64 = weight_comparison_chart(baseline_params, calibrated_params)
        html_parts.append(f'<img src="data:image/png;base64,{chart_b64}" alt="Weight comparison">')

    # 3. Scatter plots
    if scatter_data:
        html_parts.append("<h2>3. Predicted vs. Measured Scatter Plots</h2>")
        for ds_name, (predicted, measured) in scatter_data.items():
            b64 = scatter_plot(predicted, measured, title=f"{ds_name}: Predicted vs. Measured")
            html_parts.append(f'<img src="data:image/png;base64,{b64}" alt="{ds_name} scatter">')

    # 4. Subscale diagnostics heatmap
    if diagnostics:
        html_parts.append("<h2>4. Subscale Diagnostic Heatmap</h2>")
        b64 = subscale_heatmap(diagnostics)
        if b64:
            html_parts.append(f'<img src="data:image/png;base64,{b64}" alt="Subscale heatmap">')

    # 5. CV stability
    if fold_rhos:
        html_parts.append("<h2>5. Cross-Validation Stability</h2>")
        b64 = cv_boxplot(fold_rhos)
        if b64:
            html_parts.append(f'<img src="data:image/png;base64,{b64}" alt="CV stability">')

    html_parts.append("</body></html>")

    output_path.write_text("\n".join(html_parts))
    logger.info("Report saved to %s", output_path)
    return output_path
