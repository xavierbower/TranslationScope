"""Pre-analysis sequence screens: IRES, signal peptide, N-terminal tags."""

from __future__ import annotations

from dataclasses import dataclass, field

CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

IRES_KEYWORDS = {"ires", "emcv", "hcv", "poliovirus", "encephalomyocarditis"}

HYDROPHOBIC = set("LIVFMA")
POSITIVE_CHARGE = set("KR")
SMALL_NEUTRAL = set("AGST")

TAG_SIGNATURES = {
    "6xHis": {"type": "codon_pattern"},
    "FLAG": {"seq": "DYKDDDDK"},
    "Strep-tag": {"seq": "WSHPQFEK"},
    "HA": {"seq": "YPYDVPDYA"},
    "Myc": {"seq": "EQKLISEEDL"},
    "MBP": {"prefix": "MKIEEGKLVI"},
    "GST": {"prefix": "MSPILGYWKI"},
}


def translate_rna(rna: str) -> str:
    """Translate RNA sequence to amino acid string."""
    aa = []
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i + 3]
        residue = CODON_TABLE.get(codon, "X")
        if residue == "*":
            break
        aa.append(residue)
    return "".join(aa)


@dataclass
class ScreenResult:
    ires_detected: bool = False
    signal_peptide_suspected: bool = False
    tag_detected: bool = False
    tag_name: str = ""
    tag_length_codons: int = 0
    warnings: list[str] = field(default_factory=list)


def screen_ires(utr_sequence: str, features_labels: list[str] | None = None) -> bool:
    """Check for IRES presence via annotations and heuristic."""
    if features_labels:
        for lab in features_labels:
            low = lab.lower()
            if any(kw in low for kw in IRES_KEYWORDS):
                return True

    # Heuristic: long UTR with low GC
    if len(utr_sequence) > 400:
        gc = (utr_sequence.count("G") + utr_sequence.count("C")) / len(utr_sequence) * 100
        if gc < 45:
            return True
    return False


def screen_signal_peptide(cds_rna: str) -> bool:
    """Check first 30 amino acids for signal peptide characteristics."""
    if len(cds_rna) < 90:
        return False
    aa = translate_rna(cds_rna[:90])
    if len(aa) < 15:
        return False

    # Look for n-region (1-5 aa after Met, may contain K/R) then h-region (7+ hydrophobic)
    # Start scanning from position 1 (after start Met)
    for n_end in range(1, min(7, len(aa))):
        n_region = aa[1:n_end + 1]
        if not any(r in POSITIVE_CHARGE for r in n_region):
            continue
        # Look for hydrophobic stretch starting after n-region
        best_h = 0
        current_h = 0
        for i in range(n_end + 1, min(n_end + 21, len(aa))):
            if aa[i] in HYDROPHOBIC:
                current_h += 1
                best_h = max(best_h, current_h)
            else:
                current_h = 0
        if best_h >= 7:
            return True
    return False


def screen_tags(cds_rna: str) -> tuple[bool, str, int]:
    """Check first 30 codons for known tags."""
    first_90nt = cds_rna[:90]
    aa = translate_rna(first_90nt)

    # 6xHis: >= 4 consecutive CAC or CAU codons
    his_run = 0
    max_his = 0
    for i in range(0, min(len(cds_rna), 90) - 2, 3):
        codon = cds_rna[i:i + 3]
        if codon in ("CAC", "CAU"):
            his_run += 1
            max_his = max(max_his, his_run)
        else:
            his_run = 0
    if max_his >= 4:
        return True, "6xHis", max_his

    # Sequence-based tags
    for tag_name, info in TAG_SIGNATURES.items():
        if tag_name == "6xHis":
            continue
        if "seq" in info:
            if info["seq"] in aa:
                return True, tag_name, len(info["seq"])
        if "prefix" in info:
            if aa[:len(info["prefix"])] == info["prefix"]:
                return True, tag_name, len(info["prefix"])

    return False, "", 0


def run_screens(
    utr_sequence: str,
    cds_sequence: str,
    feature_labels: list[str] | None = None,
) -> ScreenResult:
    """Run all pre-analysis screens."""
    result = ScreenResult()

    # IRES
    result.ires_detected = screen_ires(utr_sequence, feature_labels)
    if result.ires_detected:
        result.warnings.append(
            "IRES sequence detected or suspected. Cap-dependent initiation models "
            "(UTR structure, AUG accessibility, uORF scores) do not apply to "
            "IRES-driven translation and have been hidden. Codon optimality scores remain valid."
        )

    # Signal peptide
    result.signal_peptide_suspected = screen_signal_peptide(cds_sequence)
    if result.signal_peptide_suspected:
        result.warnings.append(
            "A signal peptide may be present in the N-terminal sequence. "
            "Translation of signal-peptide-containing mRNAs involves co-translational "
            "SRP recognition, which affects ribosome behavior in ways not captured by this model. "
            "Scores may be systematically underestimated for secreted or membrane proteins."
        )

    # Tags
    tag_detected, tag_name, tag_len = screen_tags(cds_sequence)
    result.tag_detected = tag_detected
    result.tag_name = tag_name
    result.tag_length_codons = tag_len
    if tag_detected:
        result.warnings.append(
            f"An N-terminal {tag_name} tag was detected. The first ~{tag_len} codons "
            "reflect tag sequence rather than your gene of interest. "
            "N-terminal codon scores may not reflect GOI design."
        )

    return result
