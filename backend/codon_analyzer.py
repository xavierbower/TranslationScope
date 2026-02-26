"""Codon optimality, uORF, Kozak, and UTR composition analysis."""

from __future__ import annotations

from dataclasses import dataclass, field

from sequence_screener import CODON_TABLE

# ---- Expression system codon tables ----

CODON_TABLES = {
    "mammalian": {
        "positive": {"GCU", "GGU", "GAU", "AAC", "AAG"},
        "negative": {"AGG", "AGA", "UCA", "UUA", "AGC", "CCC"},
    },
    "ecoli": {
        "positive": {"GAA", "AAA", "GUU", "AUU", "GCU", "GGU", "GAU", "CGU"},
        "negative": {"AGG", "AGA", "CUA", "CGA", "CGG", "AUA", "CCC", "GGG"},
    },
    "yeast": {
        "positive": {"GAA", "AAA", "GUU", "AUU", "UUG", "UCU", "GAU", "AAU"},
        "negative": {"CGU", "CGC", "CGA", "CGG", "AGA", "AGG", "CUC", "CUA"},
    },
    "insect": {
        "positive": {"GAA", "AAA", "GUU", "AUU", "UUU", "GAU", "AAU", "UCA"},
        "negative": {"CGC", "CGG", "CGA", "AGG", "GGC", "CCC", "ACC"},
    },
}

STOP_CODONS = {"UAA", "UAG", "UGA"}


@dataclass
class CodonInfo:
    codon: str = ""
    aa: str = ""
    category: str = "neutral"  # positive / negative / neutral


@dataclass
class CodonResult:
    cos: float = 0.0
    positive_fraction: float = 0.0
    negative_fraction: float = 0.0
    codon_table: list[CodonInfo] = field(default_factory=list)
    nterm_cos: float = 0.0
    nterm_classification: str = ""
    nterm_flag: bool = False
    nterm_codons: list[CodonInfo] = field(default_factory=list)


@dataclass
class UorfInfo:
    position: int = 0
    distance_to_cds: int = 0
    type: str = ""  # contained / overlapping
    impact: str = ""  # high / moderate


@dataclass
class UorfResult:
    total_uaugs: int = 0
    high_impact_uaugs: int = 0
    uaug_list: list[UorfInfo] = field(default_factory=list)


@dataclass
class KozakResult:
    kozak_score: int = 0
    kozak_class: str = ""
    kozak_sequence: str = ""


@dataclass
class UtrCompositionResult:
    au_content: float = 0.0
    gc_content: float = 0.0
    gg_frequency: float = 0.0
    utr_length: int = 0


@dataclass
class LengthResult:
    cds_length: int = 0
    transcript_length: int = 0
    cds_length_log10: float = 0.0
    transcript_length_log10: float = 0.0


def _get_codons(rna: str) -> list[str]:
    return [rna[i:i + 3] for i in range(0, len(rna) - 2, 3)]


def analyze_codons(cds_rna: str, system: str = "mammalian") -> CodonResult:
    """Compute overall and N-terminal codon optimality scores."""
    import math
    table = CODON_TABLES.get(system, CODON_TABLES["mammalian"])
    pos_set = table["positive"]
    neg_set = table["negative"]

    codons = _get_codons(cds_rna)
    if not codons:
        return CodonResult()

    total = len(codons)
    pos_count = sum(1 for c in codons if c in pos_set)
    neg_count = sum(1 for c in codons if c in neg_set)
    pos_frac = pos_count / total
    neg_frac = neg_count / total
    cos = pos_frac - neg_frac

    codon_table = []
    for c in codons:
        aa = CODON_TABLE.get(c, "X")
        if c in pos_set:
            cat = "positive"
        elif c in neg_set:
            cat = "negative"
        else:
            cat = "neutral"
        codon_table.append(CodonInfo(codon=c, aa=aa, category=cat))

    # N-terminal (codons 1-10, skip start Met at position 0)
    nterm = codons[1:11] if len(codons) > 1 else []
    nterm_codons = []
    nterm_pos = 0
    nterm_neg = 0
    for c in nterm:
        aa = CODON_TABLE.get(c, "X")
        if c in pos_set:
            cat = "positive"
            nterm_pos += 1
        elif c in neg_set:
            cat = "negative"
            nterm_neg += 1
        else:
            cat = "neutral"
        nterm_codons.append(CodonInfo(codon=c, aa=aa, category=cat))

    nterm_total = len(nterm) if nterm else 1
    nterm_cos = (nterm_pos - nterm_neg) / nterm_total

    if nterm_cos > 0:
        nterm_class = "N-terminal codon context favorable"
    elif nterm_cos < -0.2:
        nterm_class = "N-terminal codon context poor — high priority fix"
    else:
        nterm_class = "N-terminal codon context neutral"

    nterm_flag = nterm_neg >= 3

    return CodonResult(
        cos=round(cos, 4),
        positive_fraction=round(pos_frac, 4),
        negative_fraction=round(neg_frac, 4),
        codon_table=codon_table,
        nterm_cos=round(nterm_cos, 4),
        nterm_classification=nterm_class,
        nterm_flag=nterm_flag,
        nterm_codons=nterm_codons,
    )


def analyze_uorfs(utr_sequence: str) -> UorfResult:
    """Scan UTR for upstream AUGs and classify uORFs."""
    result = UorfResult()
    utr_len = len(utr_sequence)
    if utr_len < 3:
        return result

    # Find all AUG positions in UTR
    positions = []
    for i in range(utr_len - 2):
        if utr_sequence[i:i + 3] == "AUG":
            positions.append(i)

    for pos in positions:
        distance = utr_len - pos
        # Translate from this uAUG to find stop
        stop_pos = None
        for j in range(pos + 3, utr_len + 300, 3):  # extend into CDS region conceptually
            if j + 3 <= len(utr_sequence):
                codon = utr_sequence[j:j + 3]
            else:
                break
            if codon in STOP_CODONS:
                stop_pos = j
                break

        if stop_pos is not None and stop_pos > utr_len:
            uorf_type = "overlapping"
        else:
            uorf_type = "contained"

        # Overlapping ORFs and close uAUGs are high impact
        if uorf_type == "overlapping" or distance < 50:
            impact = "high"
        else:
            impact = "moderate"

        info = UorfInfo(
            position=pos,
            distance_to_cds=distance,
            type=uorf_type,
            impact=impact,
        )
        result.uaug_list.append(info)

    result.total_uaugs = len(result.uaug_list)
    result.high_impact_uaugs = sum(1 for u in result.uaug_list if u.impact == "high")
    return result


def analyze_kozak(utr_sequence: str, cds_sequence: str) -> KozakResult:
    """Score Kozak context around main AUG."""
    # Build window: -6 to +6 around AUG
    # AUG is at end of UTR / start of CDS
    upstream = utr_sequence[-6:] if len(utr_sequence) >= 6 else utr_sequence
    downstream = cds_sequence[3:9] if len(cds_sequence) >= 9 else cds_sequence[3:]

    # Pad if needed
    upstream = "N" * (6 - len(upstream)) + upstream
    downstream = downstream + "N" * (6 - len(downstream))

    kozak_seq = upstream + "AUG" + downstream  # positions 0-5 = -6 to -1, 6-8 = AUG, 9-14 = +1 to +6

    score = 0
    # Position -3 (index 3): A or G -> +2
    if len(kozak_seq) > 3 and kozak_seq[3] in ("A", "G"):
        score += 2
    # Position -2 (index 4): C -> +1
    if len(kozak_seq) > 4 and kozak_seq[4] == "C":
        score += 1
    # Position -1 (index 5): C -> +1
    if len(kozak_seq) > 5 and kozak_seq[5] == "C":
        score += 1
    # Position +4 (index 9): G -> +2
    if len(kozak_seq) > 9 and kozak_seq[9] == "G":
        score += 2
    # Position +5 (index 10): G -> +1
    if len(kozak_seq) > 10 and kozak_seq[10] == "G":
        score += 1

    if score >= 6:
        klass = "Strong"
    elif score >= 4:
        klass = "Moderate"
    else:
        klass = "Weak"

    return KozakResult(kozak_score=score, kozak_class=klass, kozak_sequence=kozak_seq)


def analyze_utr_composition(utr_sequence: str) -> UtrCompositionResult:
    """Compute UTR sequence composition metrics."""
    length = len(utr_sequence)
    if length == 0:
        return UtrCompositionResult()

    a_count = utr_sequence.count("A")
    u_count = utr_sequence.count("U")
    g_count = utr_sequence.count("G")
    c_count = utr_sequence.count("C")

    au = (a_count + u_count) / length * 100
    gc = (g_count + c_count) / length * 100

    # GG dinucleotide frequency
    gg_count = sum(1 for i in range(length - 1) if utr_sequence[i:i + 2] == "GG")
    gg_freq = gg_count / (length - 1) * 100 if length > 1 else 0

    return UtrCompositionResult(
        au_content=round(au, 2),
        gc_content=round(gc, 2),
        gg_frequency=round(gg_freq, 2),
        utr_length=length,
    )


def compute_lengths(utr_sequence: str, cds_sequence: str) -> LengthResult:
    """Compute CDS and transcript lengths."""
    import math
    cds_len = len(cds_sequence)
    transcript_len = len(utr_sequence) + cds_len
    return LengthResult(
        cds_length=cds_len,
        transcript_length=transcript_len,
        cds_length_log10=round(math.log10(cds_len), 4) if cds_len > 0 else 0,
        transcript_length_log10=round(math.log10(transcript_len), 4) if transcript_len > 0 else 0,
    )
