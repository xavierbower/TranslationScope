"""CLI entry point for the calibration pipeline.

Usage:
    python -m calibration.cli download --all
    python -m calibration.cli download --dataset persist_seq
    python -m calibration.cli evaluate --params default
    python -m calibration.cli evaluate --params path/to/params.json
    python -m calibration.cli optimize --dataset persist_seq --mode weights
    python -m calibration.cli optimize --dataset all --mode full --trials 500
    python -m calibration.cli report --baseline default --optimized optimized_params.json
    python -m calibration.cli export --params optimized_params.json
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import os

# Ensure backend is on path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from calibration.config import (
    DATA_DIR, DEFAULT_K_FOLDS, DEFAULT_OPTUNA_TRIALS,
    DEFAULT_PRIOR_LAMBDA, REPORT_DIR,
)
from calibration.parameters import ParameterSet, DEFAULT_PARAMS
from calibration.datasets.base import BenchmarkDataset

logger = logging.getLogger("calibration")


DATASET_LOADERS = {
    "persist_seq": "calibration.datasets.persist_seq",
    "sample_mpra": "calibration.datasets.sample_mpra",
    "icor": "calibration.datasets.icor",
    "zheng_ribonn": "calibration.datasets.zheng_ribonn",
}


def _get_loader(name: str):
    """Dynamically import a dataset loader module."""
    import importlib
    mod = importlib.import_module(DATASET_LOADERS[name])
    return mod


def _load_params(path_or_default: str) -> ParameterSet:
    if path_or_default in ("default", "baseline", None):
        return DEFAULT_PARAMS
    return ParameterSet.from_json(path_or_default)


def _progress(current: int, total: int):
    pct = current / total * 100
    bar = "#" * int(pct // 2.5) + "-" * (40 - int(pct // 2.5))
    print(f"\r  [{bar}] {current}/{total} ({pct:.0f}%)", end="", flush=True)
    if current == total:
        print()


# ---- Commands ----

def cmd_download(args):
    """Download benchmark datasets."""
    datasets = list(DATASET_LOADERS.keys()) if args.all else [args.dataset]

    for ds_name in datasets:
        if ds_name not in DATASET_LOADERS:
            print(f"Unknown dataset: {ds_name}")
            print(f"Available: {', '.join(DATASET_LOADERS.keys())}")
            sys.exit(1)

        print(f"Downloading {ds_name}...")
        loader = _get_loader(ds_name)
        try:
            path = loader.download(force=args.force)
            print(f"  Saved to {path}")
        except Exception as e:
            print(f"  Failed: {e}")
            if not args.all:
                sys.exit(1)


def cmd_evaluate(args):
    """Evaluate current or specified parameters against benchmark datasets."""
    params = _load_params(args.params)
    datasets_to_eval = list(DATASET_LOADERS.keys()) if args.dataset == "all" else [args.dataset]

    from calibration.scoring import score_dataset, extract_features
    from calibration.metrics import evaluate, format_metrics, subscale_diagnostics

    for ds_name in datasets_to_eval:
        if ds_name not in DATASET_LOADERS:
            print(f"Skipping unknown dataset: {ds_name}")
            continue

        print(f"\n{'='*60}")
        print(f"Evaluating on {ds_name}")
        print(f"{'='*60}")

        try:
            loader = _get_loader(ds_name)
            dataset = loader.load()
        except FileNotFoundError as e:
            print(f"  {e}")
            continue

        if not dataset.records:
            print("  No records loaded.")
            continue

        print(f"  Loaded {len(dataset.records)} records")
        print(f"  Scoring (skip_structure={args.skip_structure})...")

        outputs = score_dataset(
            dataset.records, params,
            skip_structure=args.skip_structure,
            progress_callback=_progress,
        )

        predicted = [o.total_score for o in outputs]
        measured = [o.measured_te for o in outputs]

        metrics = evaluate(predicted, measured)
        print(format_metrics(metrics, label=ds_name))

        # Subscale diagnostics
        if args.diagnostics:
            fv_dicts = []
            for o in outputs:
                fv_dicts.append({
                    "aug_accessibility": o.features.aug_accessibility,
                    "cos": o.features.cos,
                    "nterm_cos": o.features.nterm_cos,
                    "kozak_score": o.features.kozak_score,
                    "au_content": o.features.au_content,
                    "utr_mfe": o.features.utr_mfe,
                    "utr_length": o.features.utr_length,
                })
            diag = subscale_diagnostics(fv_dicts, measured)
            print("\n  Per-feature Spearman with measured TE:")
            for fname, rho in sorted(diag.feature_correlations.items(), key=lambda x: -abs(x[1])):
                pval = diag.feature_pvalues[fname]
                sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
                print(f"    {fname:25s} rho={rho:+.4f}  p={pval:.2e} {sig}")


def cmd_optimize(args):
    """Run parameter optimization."""
    from calibration.scoring import extract_features_batch
    from calibration.optimizer import optimize_weights, optimize_full

    datasets_to_use = list(DATASET_LOADERS.keys()) if args.dataset == "all" else [args.dataset]

    # Load and combine records
    all_records = []
    for ds_name in datasets_to_use:
        try:
            loader = _get_loader(ds_name)
            dataset = loader.load()
            all_records.extend(dataset.records)
            print(f"Loaded {len(dataset.records)} records from {ds_name}")
        except FileNotFoundError as e:
            print(f"Skipping {ds_name}: {e}")

    if not all_records:
        print("No records available for optimization.")
        sys.exit(1)

    print(f"\nTotal records: {len(all_records)}")
    print(f"Extracting features (skip_structure={args.skip_structure})...")

    features_and_records = extract_features_batch(
        all_records,
        skip_structure=args.skip_structure,
        progress_callback=_progress,
    )

    if args.mode == "weights":
        print(f"\nRunning weight-only optimization (Nelder-Mead)...")
        result = optimize_weights(
            features_and_records,
            k_folds=args.folds,
            prior_lambda=args.prior_lambda,
            seed=args.seed,
        )
    elif args.mode == "full":
        print(f"\nRunning full optimization (Optuna, {args.trials} trials)...")
        result = optimize_full(
            features_and_records,
            n_trials=args.trials,
            k_folds=args.folds,
            prior_lambda=args.prior_lambda,
            seed=args.seed,
        )
    else:
        print(f"Unknown mode: {args.mode}")
        sys.exit(1)

    print(f"\nOptimization complete!")
    print(f"  Mode: {result.mode}")
    print(f"  Trials: {result.n_trials}")
    print(f"  Best Spearman rho: {result.best_spearman:.4f}")
    print(f"  CV mean ± std: {result.cv_spearman_mean:.4f} ± {result.cv_spearman_std:.4f}")
    print(f"  Per-fold rhos: {[f'{r:.4f}' for r in result.cv_fold_rhos]}")

    # Save optimized parameters
    out_path = args.output or str(DATA_DIR / f"optimized_params_{result.mode}.json")
    result.best_params.to_json(out_path)
    print(f"\nOptimized parameters saved to: {out_path}")

    # Print weight comparison
    print("\nWeight comparison (baseline -> optimized):")
    labels = ["AUG", "COS", "N-term", "uORF", "Kozak", "UTR Comp", "UTR Struct"]
    for label, base_w, opt_w in zip(labels, DEFAULT_PARAMS.weights, result.best_params.weights):
        delta = opt_w - base_w
        print(f"  {label:12s}: {base_w:5.1f} -> {opt_w:5.1f}  ({delta:+.1f})")


def cmd_report(args):
    """Generate HTML calibration report."""
    from calibration.scoring import score_dataset
    from calibration.metrics import evaluate, subscale_diagnostics, SubscaleDiagnostics
    from calibration.report import generate_html_report

    baseline_params = _load_params(args.baseline)
    calibrated_params = _load_params(args.optimized) if args.optimized else None

    baseline_metrics = {}
    calibrated_metrics = {}
    scatter_data = {}
    all_diagnostics = {}
    fold_rhos = {}

    datasets_to_eval = list(DATASET_LOADERS.keys()) if args.dataset == "all" else [args.dataset]

    for ds_name in datasets_to_eval:
        try:
            loader = _get_loader(ds_name)
            dataset = loader.load()
        except FileNotFoundError:
            continue

        if not dataset.records:
            continue

        print(f"Scoring {ds_name} with baseline params...")
        base_outputs = score_dataset(dataset.records, baseline_params,
                                      skip_structure=args.skip_structure,
                                      progress_callback=_progress)
        base_pred = [o.total_score for o in base_outputs]
        measured = [o.measured_te for o in base_outputs]
        baseline_metrics[ds_name] = evaluate(base_pred, measured)
        scatter_data[f"{ds_name} (baseline)"] = (base_pred, measured)

        if calibrated_params:
            print(f"Scoring {ds_name} with calibrated params...")
            cal_outputs = score_dataset(dataset.records, calibrated_params,
                                         skip_structure=args.skip_structure,
                                         progress_callback=_progress)
            cal_pred = [o.total_score for o in cal_outputs]
            calibrated_metrics[ds_name] = evaluate(cal_pred, measured)
            scatter_data[f"{ds_name} (calibrated)"] = (cal_pred, measured)

        # Subscale diagnostics
        fv_dicts = [{
            "aug_accessibility": o.features.aug_accessibility,
            "cos": o.features.cos,
            "nterm_cos": o.features.nterm_cos,
            "kozak_score": o.features.kozak_score,
            "au_content": o.features.au_content,
            "utr_mfe": o.features.utr_mfe,
        } for o in base_outputs]
        all_diagnostics[ds_name] = subscale_diagnostics(fv_dicts, measured)

    output_path = generate_html_report(
        baseline_metrics=baseline_metrics,
        calibrated_metrics=calibrated_metrics if calibrated_metrics else None,
        baseline_params=baseline_params,
        calibrated_params=calibrated_params,
        scatter_data=scatter_data,
        diagnostics=all_diagnostics,
        fold_rhos=fold_rhos if fold_rhos else None,
        output_path=args.output and Path(args.output),
    )
    print(f"\nReport saved to: {output_path}")


def cmd_export(args):
    """Export optimized parameters back into scorer.py."""
    params = _load_params(args.params)

    print("Optimized parameters:")
    for name in ParameterSet.param_names():
        val = getattr(params, name)
        default_val = getattr(DEFAULT_PARAMS, name)
        if val != default_val:
            print(f"  {name}: {default_val} -> {val}")

    # Save as a reference JSON alongside scorer.py
    from calibration.config import BACKEND_DIR
    ref_path = BACKEND_DIR / "calibrated_params.json"
    params.to_json(ref_path)
    print(f"\nCalibrated parameters saved to: {ref_path}")
    print("To use these parameters, pass them to compute_score(params=ParameterSet.from_json(...))")
    print("Or update DEFAULT_PARAMS in calibration/parameters.py")


# ---- Main parser ----

def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    parser = argparse.ArgumentParser(
        prog="calibration",
        description="TranslationScope empirical calibration pipeline",
    )
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # download
    dl = subparsers.add_parser("download", help="Download benchmark datasets")
    dl.add_argument("--dataset", choices=list(DATASET_LOADERS.keys()), help="Dataset to download")
    dl.add_argument("--all", action="store_true", help="Download all datasets")
    dl.add_argument("--force", action="store_true", help="Re-download even if cached")

    # evaluate
    ev = subparsers.add_parser("evaluate", help="Evaluate parameters against benchmarks")
    ev.add_argument("--params", default="default", help="Parameter set (default or path to JSON)")
    ev.add_argument("--dataset", default="all", help="Dataset to evaluate on")
    ev.add_argument("--skip-structure", action="store_true", help="Skip ViennaRNA folding (faster)")
    ev.add_argument("--diagnostics", action="store_true", help="Show per-feature diagnostics")

    # optimize
    opt = subparsers.add_parser("optimize", help="Optimize scoring parameters")
    opt.add_argument("--dataset", default="persist_seq", help="Dataset(s) to optimize on")
    opt.add_argument("--mode", choices=["weights", "full"], default="weights",
                     help="Optimization mode: weights-only or full parameter search")
    opt.add_argument("--trials", type=int, default=DEFAULT_OPTUNA_TRIALS,
                     help="Number of Optuna trials (full mode only)")
    opt.add_argument("--folds", type=int, default=DEFAULT_K_FOLDS, help="Number of CV folds")
    opt.add_argument("--prior-lambda", type=float, default=DEFAULT_PRIOR_LAMBDA,
                     help="Regularization strength toward defaults")
    opt.add_argument("--seed", type=int, default=42, help="Random seed")
    opt.add_argument("--skip-structure", action="store_true", help="Skip ViennaRNA folding")
    opt.add_argument("--output", help="Output path for optimized params JSON")

    # report
    rpt = subparsers.add_parser("report", help="Generate HTML calibration report")
    rpt.add_argument("--baseline", default="default", help="Baseline params")
    rpt.add_argument("--optimized", help="Optimized params JSON path")
    rpt.add_argument("--dataset", default="all", help="Dataset(s) to include")
    rpt.add_argument("--skip-structure", action="store_true", help="Skip ViennaRNA folding")
    rpt.add_argument("--output", help="Output path for HTML report")

    # export
    exp = subparsers.add_parser("export", help="Export calibrated parameters")
    exp.add_argument("--params", required=True, help="Path to optimized params JSON")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    commands = {
        "download": cmd_download,
        "evaluate": cmd_evaluate,
        "optimize": cmd_optimize,
        "report": cmd_report,
        "export": cmd_export,
    }
    commands[args.command](args)


if __name__ == "__main__":
    main()
