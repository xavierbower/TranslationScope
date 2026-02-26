"""K-fold cross-validation and train/test splitting."""

from __future__ import annotations

import logging
import random
from dataclasses import dataclass, field

import numpy as np

from calibration.datasets.base import BenchmarkRecord, BenchmarkDataset
from calibration.config import DEFAULT_K_FOLDS, DEFAULT_TEST_FRACTION

logger = logging.getLogger(__name__)


@dataclass
class DataSplit:
    """A single train/test split."""
    train: list[BenchmarkRecord] = field(default_factory=list)
    test: list[BenchmarkRecord] = field(default_factory=list)
    fold_index: int = 0


def stratified_train_test_split(
    records: list[BenchmarkRecord],
    test_fraction: float = DEFAULT_TEST_FRACTION,
    seed: int = 42,
) -> DataSplit:
    """Split records into train/test, stratified by TE quartile.

    Args:
        records: All records
        test_fraction: Fraction reserved for test set
        seed: Random seed for reproducibility

    Returns:
        DataSplit with train and test sets.
    """
    rng = random.Random(seed)

    # Sort by TE and split into quartiles
    sorted_recs = sorted(records, key=lambda r: r.measured_te)
    n = len(sorted_recs)
    n_quartiles = 4
    quartile_size = n // n_quartiles

    train = []
    test = []

    for q in range(n_quartiles):
        start = q * quartile_size
        end = start + quartile_size if q < n_quartiles - 1 else n
        quartile = sorted_recs[start:end]
        rng.shuffle(quartile)

        n_test = max(1, int(len(quartile) * test_fraction))
        test.extend(quartile[:n_test])
        train.extend(quartile[n_test:])

    rng.shuffle(train)
    rng.shuffle(test)

    logger.info("Split %d records: %d train, %d test (%.0f%%)",
                n, len(train), len(test), test_fraction * 100)

    return DataSplit(train=train, test=test)


def k_fold_splits(
    records: list[BenchmarkRecord],
    k: int = DEFAULT_K_FOLDS,
    seed: int = 42,
) -> list[DataSplit]:
    """Generate K stratified cross-validation folds.

    Records are sorted by TE, divided into k strata, and each stratum
    is split across folds to ensure balanced TE distribution.
    """
    rng = random.Random(seed)
    sorted_recs = sorted(records, key=lambda r: r.measured_te)

    # Assign each record to a fold via round-robin on sorted order
    fold_assignments = [[] for _ in range(k)]
    indices = list(range(len(sorted_recs)))
    rng.shuffle(indices)  # Shuffle within sort-order blocks

    # Stratified: split sorted list into k blocks, shuffle within blocks
    block_size = len(sorted_recs) // k
    for block_start in range(0, len(sorted_recs), block_size):
        block = sorted_recs[block_start:block_start + block_size]
        rng.shuffle(block)
        for i, rec in enumerate(block):
            fold_assignments[(block_start // block_size + i) % k].append(rec)

    splits = []
    for fold_idx in range(k):
        test = fold_assignments[fold_idx]
        train = []
        for j in range(k):
            if j != fold_idx:
                train.extend(fold_assignments[j])
        rng.shuffle(train)
        splits.append(DataSplit(train=train, test=test, fold_index=fold_idx))

    logger.info("Created %d-fold CV: ~%d train, ~%d test per fold",
                k, len(splits[0].train), len(splits[0].test))

    return splits


def cross_validate_score(
    records: list[BenchmarkRecord],
    score_fn,
    k: int = DEFAULT_K_FOLDS,
    seed: int = 42,
) -> list[float]:
    """Run k-fold CV and return per-fold Spearman rho values.

    Args:
        records: All records (train set — held-out test already removed)
        score_fn: Callable(train_records, test_records) -> float (Spearman rho on test)
        k: Number of folds
        seed: Random seed

    Returns:
        List of k Spearman rho values.
    """
    splits = k_fold_splits(records, k, seed)
    rhos = []

    for split in splits:
        rho = score_fn(split.train, split.test)
        rhos.append(rho)
        logger.info("Fold %d: Spearman rho = %.4f", split.fold_index, rho)

    mean_rho = np.mean(rhos)
    std_rho = np.std(rhos)
    logger.info("CV mean Spearman rho: %.4f ± %.4f", mean_rho, std_rho)

    return rhos
