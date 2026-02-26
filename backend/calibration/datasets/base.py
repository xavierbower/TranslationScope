"""Common dataset format for calibration benchmarks."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class BenchmarkRecord:
    """A single mRNA construct with measured translation efficiency."""
    record_id: str = ""
    utr_sequence: str = ""          # 5' UTR (RNA alphabet: AUGC)
    cds_sequence: str = ""          # CDS (RNA alphabet, includes start AUG, no stop)
    measured_te: float = 0.0        # Measured translation efficiency / expression
    expression_system: str = "mammalian"
    has_full_cds: bool = True       # False for MPRA (fixed reporter CDS)
    metadata: dict = field(default_factory=dict)


@dataclass
class BenchmarkDataset:
    """A collection of benchmark records from one study."""
    name: str = ""
    records: list[BenchmarkRecord] = field(default_factory=list)
    te_metric_name: str = ""        # e.g. "ribosome_load", "MRL", "protein_level"
    description: str = ""
    citation: str = ""

    def __len__(self) -> int:
        return len(self.records)

    def utr_only_records(self) -> list[BenchmarkRecord]:
        """Records that only have UTR variation (fixed CDS)."""
        return [r for r in self.records if not r.has_full_cds]

    def full_cds_records(self) -> list[BenchmarkRecord]:
        """Records with both UTR and CDS variation."""
        return [r for r in self.records if r.has_full_cds]

    def te_values(self) -> list[float]:
        return [r.measured_te for r in self.records]

    def subsample(self, n: int, seed: int = 42, stratify: bool = True) -> BenchmarkDataset:
        """Return a subsampled dataset, optionally stratified by TE decile."""
        import random
        rng = random.Random(seed)

        if n >= len(self.records):
            return self

        if not stratify:
            sampled = rng.sample(self.records, n)
        else:
            # Stratify by TE decile
            sorted_records = sorted(self.records, key=lambda r: r.measured_te)
            n_deciles = 10
            decile_size = len(sorted_records) // n_deciles
            per_decile = n // n_deciles
            remainder = n - per_decile * n_deciles

            sampled = []
            for i in range(n_deciles):
                start = i * decile_size
                end = start + decile_size if i < n_deciles - 1 else len(sorted_records)
                decile_records = sorted_records[start:end]
                take = per_decile + (1 if i < remainder else 0)
                take = min(take, len(decile_records))
                sampled.extend(rng.sample(decile_records, take))

        return BenchmarkDataset(
            name=f"{self.name}_sub{n}",
            records=sampled,
            te_metric_name=self.te_metric_name,
            description=f"Subsampled {n} from {self.name}",
            citation=self.citation,
        )
