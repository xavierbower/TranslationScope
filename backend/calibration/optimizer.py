"""Parameter optimization: Nelder-Mead (weight-only) and Optuna Bayesian (full)."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from calibration.parameters import ParameterSet, DEFAULT_PARAMS
from calibration.scoring import score_from_features, FeatureVector
from calibration.metrics import spearman
from calibration.cross_validation import k_fold_splits, DataSplit
from calibration.datasets.base import BenchmarkRecord
from calibration.config import (
    DEFAULT_K_FOLDS, DEFAULT_NELDER_MEAD_MAXITER,
    DEFAULT_OPTUNA_TRIALS, DEFAULT_PRIOR_LAMBDA,
)

logger = logging.getLogger(__name__)


@dataclass
class OptimizationResult:
    """Result of parameter optimization."""
    best_params: ParameterSet = field(default_factory=ParameterSet)
    best_spearman: float = 0.0
    cv_spearman_mean: float = 0.0
    cv_spearman_std: float = 0.0
    cv_fold_rhos: list[float] = field(default_factory=list)
    n_trials: int = 0
    mode: str = ""


def _dirichlet_to_weights(raw: np.ndarray, total: float = 100.0) -> np.ndarray:
    """Convert unconstrained values to weights summing to `total` via softmax."""
    exp_vals = np.exp(raw - np.max(raw))  # numerical stability
    return exp_vals / exp_vals.sum() * total


def _weights_to_dirichlet(weights: np.ndarray, total: float = 100.0) -> np.ndarray:
    """Convert weights to unconstrained parameterization (log-proportions)."""
    props = weights / total
    props = np.clip(props, 1e-6, 1.0)
    return np.log(props)


def _evaluate_params_on_split(
    params: ParameterSet,
    features_and_records: list[tuple[BenchmarkRecord, FeatureVector]],
    test_indices: set[int],
) -> float:
    """Compute Spearman rho for test set using given params."""
    predicted = []
    measured = []

    for i, (rec, fv) in enumerate(features_and_records):
        if i in test_indices:
            score = score_from_features(fv, params)
            predicted.append(score)
            measured.append(rec.measured_te)

    if len(predicted) < 3:
        return 0.0

    rho, _ = spearman(predicted, measured)
    return rho


def _cv_objective(
    params: ParameterSet,
    features_and_records: list[tuple[BenchmarkRecord, FeatureVector]],
    folds: list[DataSplit],
    prior_lambda: float = 0.0,
) -> float:
    """Negative mean CV Spearman rho (for minimization).

    Includes optional regularization toward default parameters.
    """
    # Build index sets for each fold
    record_id_to_idx = {rec.record_id: i for i, (rec, _) in enumerate(features_and_records)}

    rhos = []
    for fold in folds:
        test_ids = {rec.record_id for rec in fold.test}
        test_indices = {record_id_to_idx[rid] for rid in test_ids if rid in record_id_to_idx}
        rho = _evaluate_params_on_split(params, features_and_records, test_indices)
        rhos.append(rho)

    mean_rho = np.mean(rhos)

    # Regularization: penalize deviation from defaults
    reg = 0.0
    if prior_lambda > 0:
        default = DEFAULT_PARAMS.to_flat_array()
        current = params.to_flat_array()
        # Normalize by default value to make penalty scale-invariant
        for d, c in zip(default, current):
            if d != 0:
                reg += ((c - d) / d) ** 2
            else:
                reg += (c - d) ** 2
        reg *= prior_lambda

    return -(mean_rho - reg)


def optimize_weights(
    features_and_records: list[tuple[BenchmarkRecord, FeatureVector]],
    k_folds: int = DEFAULT_K_FOLDS,
    max_iter: int = DEFAULT_NELDER_MEAD_MAXITER,
    prior_lambda: float = DEFAULT_PRIOR_LAMBDA,
    seed: int = 42,
) -> OptimizationResult:
    """Optimize weight allocations only (7 weights, fixed breakpoints).

    Uses Dirichlet reparameterization + Nelder-Mead.
    """
    records = [rec for rec, _ in features_and_records]
    folds = k_fold_splits(records, k_folds, seed)

    # Initial weights in Dirichlet space
    initial_weights = np.array(DEFAULT_PARAMS.weights)
    x0 = _weights_to_dirichlet(initial_weights)

    def objective(raw_weights: np.ndarray) -> float:
        weights = _dirichlet_to_weights(raw_weights)
        p = ParameterSet(
            w_aug=weights[0], w_cos=weights[1], w_nterm=weights[2],
            w_uorf=weights[3], w_kozak=weights[4], w_utr_comp=weights[5],
            w_utr_struct=weights[6],
        )
        return _cv_objective(p, features_and_records, folds, prior_lambda)

    logger.info("Starting weight-only optimization (Nelder-Mead, %d iterations)", max_iter)

    result = minimize(
        objective, x0,
        method="Nelder-Mead",
        options={"maxiter": max_iter, "xatol": 1e-4, "fatol": 1e-5},
    )

    best_weights = _dirichlet_to_weights(result.x)
    best_params = ParameterSet(
        w_aug=best_weights[0], w_cos=best_weights[1], w_nterm=best_weights[2],
        w_uorf=best_weights[3], w_kozak=best_weights[4], w_utr_comp=best_weights[5],
        w_utr_struct=best_weights[6],
    )

    # Final CV evaluation
    final_obj = _cv_objective(best_params, features_and_records, folds, 0.0)
    cv_rhos = []
    record_id_to_idx = {rec.record_id: i for i, (rec, _) in enumerate(features_and_records)}
    for fold in folds:
        test_ids = {rec.record_id for rec in fold.test}
        test_indices = {record_id_to_idx[rid] for rid in test_ids if rid in record_id_to_idx}
        rho = _evaluate_params_on_split(best_params, features_and_records, test_indices)
        cv_rhos.append(rho)

    opt_result = OptimizationResult(
        best_params=best_params,
        best_spearman=-result.fun,
        cv_spearman_mean=float(np.mean(cv_rhos)),
        cv_spearman_std=float(np.std(cv_rhos)),
        cv_fold_rhos=cv_rhos,
        n_trials=result.nfev,
        mode="weights",
    )

    logger.info("Weight optimization complete: mean Spearman rho = %.4f ± %.4f",
                opt_result.cv_spearman_mean, opt_result.cv_spearman_std)
    logger.info("Optimized weights: %s", [f"{w:.1f}" for w in best_weights])

    return opt_result


def optimize_full(
    features_and_records: list[tuple[BenchmarkRecord, FeatureVector]],
    n_trials: int = DEFAULT_OPTUNA_TRIALS,
    k_folds: int = DEFAULT_K_FOLDS,
    prior_lambda: float = DEFAULT_PRIOR_LAMBDA,
    seed: int = 42,
) -> OptimizationResult:
    """Full parameter optimization using Optuna Bayesian search.

    Optimizes all ~30 parameters: weights + breakpoints + penalties.
    """
    try:
        import optuna
    except ImportError:
        raise ImportError(
            "Optuna is required for full optimization. "
            "Install with: pip install optuna>=3.4.0"
        )

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    records = [rec for rec, _ in features_and_records]
    folds = k_fold_splits(records, k_folds, seed)
    default = DEFAULT_PARAMS

    best_result = [None]
    best_rho = [-float("inf")]

    def objective(trial: optuna.Trial) -> float:
        # Weights (constrained to sum to 100)
        raw_w = np.array([
            trial.suggest_float("raw_w_aug", -2, 2),
            trial.suggest_float("raw_w_cos", -2, 2),
            trial.suggest_float("raw_w_nterm", -2, 2),
            trial.suggest_float("raw_w_uorf", -2, 2),
            trial.suggest_float("raw_w_kozak", -2, 2),
            trial.suggest_float("raw_w_utr_comp", -2, 2),
            trial.suggest_float("raw_w_utr_struct", -2, 2),
        ])
        weights = _dirichlet_to_weights(raw_w)

        p = ParameterSet(
            w_aug=weights[0], w_cos=weights[1], w_nterm=weights[2],
            w_uorf=weights[3], w_kozak=weights[4], w_utr_comp=weights[5],
            w_utr_struct=weights[6],

            # AUG breakpoints (monotone decreasing)
            aug_bp1=trial.suggest_float("aug_bp1", 0.70, 0.95),
            aug_bp2=trial.suggest_float("aug_bp2", 0.50, 0.85),
            aug_bp3=trial.suggest_float("aug_bp3", 0.30, 0.70),
            aug_bp4=trial.suggest_float("aug_bp4", 0.10, 0.50),

            aug_frac1=1.0,
            aug_frac2=trial.suggest_float("aug_frac2", 0.5, 0.95),
            aug_frac3=trial.suggest_float("aug_frac3", 0.2, 0.7),
            aug_frac4=trial.suggest_float("aug_frac4", 0.05, 0.4),

            # COS breakpoints
            cos_bp1=trial.suggest_float("cos_bp1", 0.02, 0.10),
            cos_bp2=trial.suggest_float("cos_bp2", 0.00, 0.05),
            cos_bp3=trial.suggest_float("cos_bp3", -0.05, 0.01),
            cos_bp4=trial.suggest_float("cos_bp4", -0.10, -0.02),

            cos_frac1=1.0,
            cos_frac2=trial.suggest_float("cos_frac2", 0.5, 0.95),
            cos_frac3=trial.suggest_float("cos_frac3", 0.2, 0.7),
            cos_frac4=trial.suggest_float("cos_frac4", 0.05, 0.4),

            # N-term breakpoints
            nterm_bp1=trial.suggest_float("nterm_bp1", -0.1, 0.1),
            nterm_bp2=trial.suggest_float("nterm_bp2", -0.3, 0.0),
            nterm_bp3=trial.suggest_float("nterm_bp3", -0.5, -0.1),

            nterm_frac1=1.0,
            nterm_frac2=trial.suggest_float("nterm_frac2", 0.4, 0.9),
            nterm_frac3=trial.suggest_float("nterm_frac3", 0.1, 0.5),

            # Kozak breakpoints
            kozak_bp1=trial.suggest_int("kozak_bp1", 5, 7),
            kozak_bp2=trial.suggest_int("kozak_bp2", 3, 5),

            kozak_frac1=1.0,
            kozak_frac2=trial.suggest_float("kozak_frac2", 0.3, 0.8),

            # UTR composition
            utr_comp_au_bp1=trial.suggest_float("utr_comp_au_bp1", 55.0, 75.0),
            utr_comp_au_bp2=trial.suggest_float("utr_comp_au_bp2", 45.0, 65.0),
            utr_comp_au_bp3=trial.suggest_float("utr_comp_au_bp3", 35.0, 55.0),

            utr_comp_gg_thresh=trial.suggest_float("utr_comp_gg_thresh", 2.0, 8.0),
            utr_comp_gg_penalty=trial.suggest_float("utr_comp_gg_penalty", 1.0, 5.0),
            utr_comp_len_thresh=trial.suggest_float("utr_comp_len_thresh", 200.0, 800.0),
            utr_comp_len_penalty=trial.suggest_float("utr_comp_len_penalty", 0.5, 4.0),

            # UTR structure (MFE) breakpoints
            utr_struct_bp1=trial.suggest_float("utr_struct_bp1", -10.0, 0.0),
            utr_struct_bp2=trial.suggest_float("utr_struct_bp2", -25.0, -5.0),
            utr_struct_bp3=trial.suggest_float("utr_struct_bp3", -45.0, -15.0),
            utr_struct_bp4=trial.suggest_float("utr_struct_bp4", -70.0, -30.0),
        )

        # Enforce monotonicity on breakpoints
        if p.aug_bp1 <= p.aug_bp2 or p.aug_bp2 <= p.aug_bp3 or p.aug_bp3 <= p.aug_bp4:
            return float("inf")
        if p.cos_bp1 <= p.cos_bp2 or p.cos_bp2 <= p.cos_bp3 or p.cos_bp3 <= p.cos_bp4:
            return float("inf")
        if p.nterm_bp1 <= p.nterm_bp2 or p.nterm_bp2 <= p.nterm_bp3:
            return float("inf")
        if p.utr_struct_bp1 <= p.utr_struct_bp2 or p.utr_struct_bp2 <= p.utr_struct_bp3 or p.utr_struct_bp3 <= p.utr_struct_bp4:
            return float("inf")
        if p.utr_comp_au_bp1 <= p.utr_comp_au_bp2 or p.utr_comp_au_bp2 <= p.utr_comp_au_bp3:
            return float("inf")

        # Enforce fraction ordering
        if p.aug_frac2 >= 1.0 or p.aug_frac3 >= p.aug_frac2 or p.aug_frac4 >= p.aug_frac3:
            return float("inf")
        if p.cos_frac2 >= 1.0 or p.cos_frac3 >= p.cos_frac2 or p.cos_frac4 >= p.cos_frac3:
            return float("inf")

        obj = _cv_objective(p, features_and_records, folds, prior_lambda)

        if -obj > best_rho[0]:
            best_rho[0] = -obj
            best_result[0] = p

        return obj

    sampler = optuna.samplers.TPESampler(seed=seed)
    study = optuna.create_study(direction="minimize", sampler=sampler)
    study.optimize(objective, n_trials=n_trials)

    best_params = best_result[0] if best_result[0] is not None else DEFAULT_PARAMS

    # Final CV evaluation
    record_id_to_idx = {rec.record_id: i for i, (rec, _) in enumerate(features_and_records)}
    cv_rhos = []
    for fold in folds:
        test_ids = {rec.record_id for rec in fold.test}
        test_indices = {record_id_to_idx[rid] for rid in test_ids if rid in record_id_to_idx}
        rho = _evaluate_params_on_split(best_params, features_and_records, test_indices)
        cv_rhos.append(rho)

    opt_result = OptimizationResult(
        best_params=best_params,
        best_spearman=best_rho[0],
        cv_spearman_mean=float(np.mean(cv_rhos)),
        cv_spearman_std=float(np.std(cv_rhos)),
        cv_fold_rhos=cv_rhos,
        n_trials=n_trials,
        mode="full",
    )

    logger.info("Full optimization complete: mean Spearman rho = %.4f ± %.4f",
                opt_result.cv_spearman_mean, opt_result.cv_spearman_std)

    return opt_result


def stability_analysis(
    params: ParameterSet,
    features_and_records: list[tuple[BenchmarkRecord, FeatureVector]],
    perturbation: float = 0.10,
) -> dict[str, float]:
    """Perturb each parameter ±perturbation% and measure sensitivity.

    Returns dict of param_name -> delta_spearman (larger = more sensitive).
    """
    # Baseline score
    predicted_base = [score_from_features(fv, params) for _, fv in features_and_records]
    measured = [rec.measured_te for rec, _ in features_and_records]
    base_rho, _ = spearman(predicted_base, measured)

    sensitivities = {}
    for name in ParameterSet.param_names():
        val = getattr(params, name)
        if val == 0:
            continue

        # Perturb up
        p_up = ParameterSet.from_dict(params.to_dict())
        setattr(p_up, name, val * (1 + perturbation))
        pred_up = [score_from_features(fv, p_up) for _, fv in features_and_records]
        rho_up, _ = spearman(pred_up, measured)

        # Perturb down
        p_down = ParameterSet.from_dict(params.to_dict())
        setattr(p_down, name, val * (1 - perturbation))
        pred_down = [score_from_features(fv, p_down) for _, fv in features_and_records]
        rho_down, _ = spearman(pred_down, measured)

        delta = max(abs(rho_up - base_rho), abs(rho_down - base_rho))
        sensitivities[name] = delta

    # Sort by sensitivity
    return dict(sorted(sensitivities.items(), key=lambda x: -x[1]))
