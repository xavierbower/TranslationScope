"""Calibration pipeline configuration."""

from __future__ import annotations

from pathlib import Path

# Root of the backend
BACKEND_DIR = Path(__file__).resolve().parent.parent

# Where downloaded/cached data lives (gitignored)
DATA_DIR = BACKEND_DIR / "calibration_data"
DATA_DIR.mkdir(exist_ok=True)

FOLD_CACHE_DIR = DATA_DIR / "fold_cache"
FOLD_CACHE_DIR.mkdir(exist_ok=True)

DATASET_DIR = DATA_DIR / "datasets"
DATASET_DIR.mkdir(exist_ok=True)

REPORT_DIR = DATA_DIR / "reports"
REPORT_DIR.mkdir(exist_ok=True)

# Default subsampling limits (reduce for ViennaRNA-heavy features)
MPRA_SUBSAMPLE_SIZE = 5000

# Cross-validation
DEFAULT_K_FOLDS = 5
DEFAULT_TEST_FRACTION = 0.20

# Optimization
DEFAULT_OPTUNA_TRIALS = 500
DEFAULT_NELDER_MEAD_MAXITER = 2000

# Regularization toward default priors
DEFAULT_PRIOR_LAMBDA = 0.01

# HTTP download settings
DOWNLOAD_TIMEOUT = 120
DOWNLOAD_CHUNK_SIZE = 8192
