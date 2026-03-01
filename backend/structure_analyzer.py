"""RNA secondary structure analysis using ViennaRNA."""

from __future__ import annotations

import subprocess
import re
from dataclasses import dataclass, field

try:
    import RNA as _RNA
    HAS_VIENNARNA = True
except ImportError:
    _RNA = None
    HAS_VIENNARNA = False


@dataclass
class UtrFoldResult:
    mfe: float = 0.0
    structure_dotbracket: str = ""
    skipped: bool = False
    skip_reason: str = ""


@dataclass
class AugAccessibilityResult:
    mfe: float = 0.0
    aug_accessibility: float = 0.0
    window_sequence: str = ""
    structure_dotbracket: str = ""
    window_truncated: bool = False
    skipped: bool = False
    skip_reason: str = ""


@dataclass
class CdsWindow:
    position: int = 0
    mfe: float = 0.0
    is_flagged: bool = False


# --------------- ViennaRNA wrappers ---------------

def _fold(seq: str) -> tuple[str, float]:
    if HAS_VIENNARNA:
        structure, mfe = _RNA.fold(seq)
        return structure, mfe
    return _fold_subprocess(seq)


def _fold_subprocess(seq: str) -> tuple[str, float]:
    proc = subprocess.run(
        ["RNAfold", "--noPS"],
        input=seq, capture_output=True, text=True, timeout=30,
    )
    lines = proc.stdout.strip().split("\n")
    if len(lines) < 2:
        return "", 0.0
    match = re.search(r"\(([^)]+)\)\s*$", lines[-1])
    mfe = float(match.group(1)) if match else 0.0
    struct = lines[-1].split()[0] if lines[-1] else ""
    return struct, mfe


def _compute_bppm_unpaired(seq: str, positions: list[int]) -> list[float]:
    """Compute unpaired probability for given positions via partition function."""
    n = len(seq)
    if HAS_VIENNARNA:
        fc = _RNA.fold_compound(seq)
        fc.pf()
        bppm = fc.bpp()  # returns tuple of tuples, 1-indexed
        unpaired = []
        for pos in positions:
            p1 = pos + 1  # convert 0-indexed to 1-indexed
            paired_sum = 0.0
            for j in range(1, n + 1):
                if j == p1:
                    continue
                i_idx = min(p1, j)
                j_idx = max(p1, j)
                if i_idx < len(bppm) and j_idx < len(bppm[i_idx]):
                    paired_sum += bppm[i_idx][j_idx]
            unpaired.append(max(0.0, 1.0 - paired_sum))
        return unpaired

    # Fallback: use RNAplfold subprocess
    return _compute_bppm_subprocess(seq, positions)


def _compute_bppm_subprocess(seq: str, positions: list[int]) -> list[float]:
    """Fallback: estimate accessibility from RNAplfold."""
    try:
        proc = subprocess.run(
            ["RNAplfold", "-W", str(len(seq)), "-u", "1"],
            input=seq, capture_output=True, text=True, timeout=30,
        )
        # Parse lunp file for unpaired probabilities
        unpaired_map: dict[int, float] = {}
        for line in proc.stdout.strip().split("\n"):
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    pos = int(parts[0])
                    prob = float(parts[1])
                    unpaired_map[pos] = prob
                except ValueError:
                    continue
        return [unpaired_map.get(p + 1, 0.5) for p in positions]
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return [0.5] * len(positions)


# --------------- Public API ---------------

def fold_utr(utr_sequence: str) -> UtrFoldResult:
    """Fold entire 5' UTR and return MFE."""
    if len(utr_sequence) < 10:
        return UtrFoldResult(skipped=True, skip_reason="UTR too short to assess")
    structure, mfe = _fold(utr_sequence)
    return UtrFoldResult(mfe=mfe, structure_dotbracket=structure)


def compute_aug_accessibility(utr_sequence: str, cds_sequence: str) -> AugAccessibilityResult:
    """Compute AUG-proximate accessibility using partition function."""
    utr_len = len(utr_sequence)

    if utr_len < 10:
        return AugAccessibilityResult(
            skipped=True,
            skip_reason="UTR too short for reliable AUG accessibility calculation — score omitted",
        )

    # Build window: up to 30 nt upstream + AUG + 14 nt downstream
    upstream = min(30, utr_len)
    downstream_avail = min(14, max(0, len(cds_sequence) - 3))
    window = utr_sequence[-upstream:] + cds_sequence[:3 + downstream_avail]
    aug_start = upstream  # position of A in AUG within window

    truncated = utr_len < 30

    structure, mfe = _fold(window)
    aug_positions = [aug_start, aug_start + 1, aug_start + 2]
    unpaired_probs = _compute_bppm_unpaired(window, aug_positions)
    accessibility = sum(unpaired_probs) / 3.0

    return AugAccessibilityResult(
        mfe=mfe,
        aug_accessibility=round(accessibility, 4),
        window_sequence=window,
        structure_dotbracket=structure,
        window_truncated=truncated,
    )


MAX_CDS_WINDOWS = 200  # Cap total fold calls to prevent OOM/timeout


def sliding_window_cds(cds_sequence: str, window_size: int = 40, step: int = 10) -> list[CdsWindow]:
    """Sliding window MFE across CDS."""
    cds_len = len(cds_sequence)
    if cds_len < window_size:
        return []

    # Increase step size for very large CDS to stay within window cap
    num_windows = (cds_len - window_size) // step + 1
    if num_windows > MAX_CDS_WINDOWS:
        step = max(step, (cds_len - window_size) // MAX_CDS_WINDOWS + 1)

    results = []
    for pos in range(0, cds_len - window_size + 1, step):
        window = cds_sequence[pos:pos + window_size]
        _, mfe = _fold(window)
        results.append(CdsWindow(
            position=pos,
            mfe=round(mfe, 2),
            is_flagged=mfe < -20,
        ))
    return results
