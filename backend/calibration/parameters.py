"""ParameterSet dataclass: all tunable scoring parameters in one place."""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path


@dataclass
class ParameterSet:
    """All tunable parameters for the TranslationScope scorer.

    Weights: max points for each subscale (must sum to 100).
    Breakpoints: thresholds that map raw feature values to scores.
    Penalties: deductions applied in specific conditions.
    """

    # ---- Subscale weights (max points, sum to 100) ----
    # Calibrated against Sample et al. 2019 MPRA (100K 5' UTRs, GEO GSE114002).
    # UTR composition weight increased from 10→14 based on its strong empirical
    # correlation with measured ribosome load (Spearman ~0.45).  Kozak and uORF
    # reduced proportionally; CDS-dependent weights kept at literature defaults
    # (no CDS variation in the calibration dataset).
    w_aug: float = 25.0
    w_cos: float = 20.0
    w_nterm: float = 10.0
    w_uorf: float = 13.0
    w_kozak: float = 8.0
    w_utr_comp: float = 14.0
    w_utr_struct: float = 10.0

    # ---- AUG accessibility breakpoints ----
    # Thresholds (descending): accessibility >= thresh -> fraction of w_aug
    aug_bp1: float = 0.85  # -> 100% of w_aug
    aug_bp2: float = 0.70  # -> 76% of w_aug
    aug_bp3: float = 0.50  # -> 44% of w_aug
    aug_bp4: float = 0.30  # -> 20% of w_aug
    # Below bp4 -> 0

    aug_frac1: float = 1.0
    aug_frac2: float = 0.76
    aug_frac3: float = 0.44
    aug_frac4: float = 0.20

    # ---- Codon optimality (COS) breakpoints ----
    cos_bp1: float = 0.05   # -> 100% of w_cos
    cos_bp2: float = 0.02   # -> 75%
    cos_bp3: float = -0.02  # -> 50%
    cos_bp4: float = -0.05  # -> 25%

    cos_frac1: float = 1.0
    cos_frac2: float = 0.75
    cos_frac3: float = 0.50
    cos_frac4: float = 0.25

    # ---- N-terminal codon breakpoints ----
    nterm_bp1: float = 0.0   # nterm_cos > bp1 -> 100% of w_nterm
    nterm_bp2: float = -0.1  # -> 70%
    nterm_bp3: float = -0.2  # -> 30%

    nterm_frac1: float = 1.0
    nterm_frac2: float = 0.70
    nterm_frac3: float = 0.30

    # ---- uORF burden breakpoints ----
    # These are categorical, not continuous, so we store the score fractions
    uorf_frac_none: float = 1.0       # total_uaugs == 0
    uorf_frac_few_low: float = 0.733  # total <= 2, high_impact == 0
    uorf_frac_few_high: float = 0.333 # total <= 2, high_impact > 0
    uorf_frac_many_low: float = 0.467 # total > 2, high_impact == 0
    # Otherwise -> 0

    # ---- Kozak breakpoints ----
    kozak_bp1: int = 6  # -> 100% of w_kozak
    kozak_bp2: int = 4  # -> 60%

    kozak_frac1: float = 1.0
    kozak_frac2: float = 0.60

    # ---- UTR composition breakpoints ----
    utr_comp_au_bp1: float = 65.0  # -> base 100%
    utr_comp_au_bp2: float = 55.0  # -> base 70%
    utr_comp_au_bp3: float = 45.0  # -> base 40%

    utr_comp_au_frac1: float = 1.0
    utr_comp_au_frac2: float = 0.70
    utr_comp_au_frac3: float = 0.40

    # Penalties (absolute points deducted, then clamped to 0)
    utr_comp_gg_thresh: float = 5.0   # gg_frequency > this -> penalty
    utr_comp_gg_penalty: float = 3.0  # points deducted (on 10-point scale -> 0.3 of w_utr_comp)
    utr_comp_len_thresh: float = 500.0  # utr_length > this -> penalty
    utr_comp_len_penalty: float = 2.0   # points deducted

    # ---- UTR structure (MFE) breakpoints ----
    utr_struct_bp1: float = -5.0   # MFE > bp1 -> 100%
    utr_struct_bp2: float = -15.0  # -> 70%
    utr_struct_bp3: float = -30.0  # -> 40%
    utr_struct_bp4: float = -50.0  # -> 20%

    utr_struct_frac1: float = 1.0
    utr_struct_frac2: float = 0.70
    utr_struct_frac3: float = 0.40
    utr_struct_frac4: float = 0.20

    # ---- Rating thresholds ----
    rating_high: float = 85.0
    rating_moderate: float = 65.0
    rating_low: float = 40.0

    @property
    def weights(self) -> list[float]:
        return [self.w_aug, self.w_cos, self.w_nterm, self.w_uorf,
                self.w_kozak, self.w_utr_comp, self.w_utr_struct]

    @property
    def weight_sum(self) -> float:
        return sum(self.weights)

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, path: str | Path) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_dict(cls, d: dict) -> ParameterSet:
        # Only pass keys that are valid fields
        valid = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in d.items() if k in valid})

    @classmethod
    def from_json(cls, path: str | Path) -> ParameterSet:
        with open(path) as f:
            return cls.from_dict(json.load(f))

    def to_flat_array(self) -> list[float]:
        """Serialize all parameters to a flat array for optimization."""
        return [getattr(self, f.name) for f in self.__dataclass_fields__.values()]

    @classmethod
    def from_flat_array(cls, arr: list[float]) -> ParameterSet:
        """Reconstruct from flat array."""
        names = [f.name for f in cls.__dataclass_fields__.values()]
        return cls(**{n: v for n, v in zip(names, arr)})

    @classmethod
    def param_names(cls) -> list[str]:
        return [f.name for f in cls.__dataclass_fields__.values()]

    @classmethod
    def default(cls) -> ParameterSet:
        return cls()


# Singleton default for backwards compatibility
DEFAULT_PARAMS = ParameterSet()
