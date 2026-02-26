"""Unit and integration tests for TranslationScope backend."""

import os
import sys
import pytest

# Ensure backend modules are importable
sys.path.insert(0, os.path.dirname(__file__))

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), "..")
DNA_FILE = os.path.join(EXAMPLE_DIR, "..", "4MMZD9_10_pRRM022-C.dna")
GBK_FILE = os.path.join(EXAMPLE_DIR, "..", "K6Z62L_3_pRRM058A_Mini_Retry.gbk")


# ---- parser.py unit tests ----

class TestParser:
    def test_parse_dna_file(self):
        from parser import parse_file
        with open(DNA_FILE, "rb") as f:
            data = f.read()
        record, features, is_circular, seq = parse_file(data, "test.dna")
        assert len(seq) > 0
        assert isinstance(features, list)
        assert len(features) > 0
        assert isinstance(is_circular, bool)

    def test_parse_gbk_file(self):
        from parser import parse_file
        with open(GBK_FILE, "rb") as f:
            data = f.read()
        record, features, is_circular, seq = parse_file(data, "test.gbk")
        assert len(seq) > 0
        assert isinstance(features, list)
        assert len(features) > 0

    def test_select_cds_candidates_dna(self):
        from parser import parse_file, select_cds_candidates
        with open(DNA_FILE, "rb") as f:
            data = f.read()
        _, features, _, _ = parse_file(data, "test.dna")
        result = select_cds_candidates(features)
        assert len(result.candidates) > 0
        for c in result.candidates:
            assert c.length_nt > 0
            assert c.label

    def test_select_cds_candidates_gbk(self):
        from parser import parse_file, select_cds_candidates
        with open(GBK_FILE, "rb") as f:
            data = f.read()
        _, features, _, _ = parse_file(data, "test.gbk")
        result = select_cds_candidates(features)
        assert len(result.candidates) > 0

    def test_extract_sequences_dna(self):
        from parser import parse_file, select_cds_candidates, extract_sequences
        with open(DNA_FILE, "rb") as f:
            data = f.read()
        _, features, _, _ = parse_file(data, "test.dna")
        result = select_cds_candidates(features)
        idx = result.candidates[0].index
        ext = extract_sequences(data, idx, "test.dna")
        assert len(ext.cds_sequence) > 0
        assert all(c in "AUGC" for c in ext.cds_sequence), f"Non-RNA chars in CDS: {set(ext.cds_sequence) - set('AUGC')}"
        assert ext.gene_name

    def test_extract_sequences_gbk(self):
        from parser import parse_file, select_cds_candidates, extract_sequences
        with open(GBK_FILE, "rb") as f:
            data = f.read()
        _, features, _, _ = parse_file(data, "test.gbk")
        result = select_cds_candidates(features)
        idx = result.candidates[0].index
        ext = extract_sequences(data, idx, "test.gbk")
        assert len(ext.cds_sequence) > 0
        assert all(c in "AUGC" for c in ext.cds_sequence)

    def test_file_too_large(self):
        from parser import parse_file
        big = b"\x00" * (11 * 1024 * 1024)
        with pytest.raises(ValueError, match="10 MB"):
            parse_file(big, "test.dna")

    def test_bad_file(self):
        from parser import parse_file
        with pytest.raises(ValueError):
            parse_file(b"not a real file", "test.dna")


# ---- sequence_screener.py unit tests ----

class TestScreener:
    def test_no_ires_short_utr(self):
        from sequence_screener import screen_ires
        assert screen_ires("AUGCAUGCAUGC") is False

    def test_ires_by_annotation(self):
        from sequence_screener import screen_ires
        assert screen_ires("AUGC", ["EMCV IRES"]) is True
        assert screen_ires("AUGC", ["my IRES element"]) is True

    def test_ires_heuristic(self):
        from sequence_screener import screen_ires
        # Long AU-rich UTR triggers heuristic
        utr = "A" * 200 + "U" * 250  # 450 nt, very low GC
        assert screen_ires(utr) is True
        # Long GC-rich UTR does not
        utr_gc = "G" * 250 + "C" * 200
        assert screen_ires(utr_gc) is False

    def test_signal_peptide_detection(self):
        from sequence_screener import screen_signal_peptide, translate_rna
        # Construct: M + R (positive, n-region) + LLLLLLLLL (9 hydrophobic, h-region) + rest
        # M=AUG, R=CGU, L=CUU (x9), then A=GCU padding — must be >= 90 nt
        sp_rna = "AUG" + "CGU" + "CUU" * 9 + "GCU" * 20  # 3+3+27+60 = 93 nt
        aa = translate_rna(sp_rna)
        assert aa[0] == "M"
        assert aa[1] == "R"  # positive charge
        assert all(a == "L" for a in aa[2:11])  # 9 hydrophobic
        assert len(sp_rna) >= 90
        assert screen_signal_peptide(sp_rna) is True

    def test_no_signal_peptide(self):
        from sequence_screener import screen_signal_peptide
        # All glycine — no signal peptide pattern
        rna = "AUG" + "GGU" * 30
        assert screen_signal_peptide(rna) is False

    def test_his_tag_detection(self):
        from sequence_screener import screen_tags
        # 6 consecutive CAC codons
        rna = "AUG" + "CAC" * 6 + "GCU" * 20
        detected, name, length = screen_tags(rna)
        assert detected is True
        assert name == "6xHis"

    def test_flag_tag_detection(self):
        from sequence_screener import screen_tags
        # FLAG = DYKDDDDK -> GAU UAC AAA GAU GAU GAU GAU AAA
        flag_codons = "GAUUACAAAGAUGAUGAUGAUAAA"
        rna = "AUG" + flag_codons + "GCU" * 20
        detected, name, _ = screen_tags(rna)
        assert detected is True
        assert name == "FLAG"

    def test_no_tag(self):
        from sequence_screener import screen_tags
        rna = "AUG" + "GCU" * 30
        detected, _, _ = screen_tags(rna)
        assert detected is False

    def test_run_screens_integration(self):
        from sequence_screener import run_screens
        result = run_screens("AUGCAUGCAUGC", "AUG" + "GCU" * 30)
        assert result.ires_detected is False
        assert isinstance(result.warnings, list)


# ---- codon_analyzer.py unit tests ----

class TestCodonAnalyzer:
    def test_analyze_codons_mammalian(self):
        from codon_analyzer import analyze_codons
        # All GCU (positive in mammalian)
        cds = "GCU" * 30
        result = analyze_codons(cds, "mammalian")
        assert result.cos > 0
        assert result.positive_fraction > 0.9

    def test_analyze_codons_negative(self):
        from codon_analyzer import analyze_codons
        # All AGG (negative in mammalian)
        cds = "AGG" * 30
        result = analyze_codons(cds, "mammalian")
        assert result.cos < 0
        assert result.negative_fraction > 0.9

    def test_analyze_codons_ecoli(self):
        from codon_analyzer import analyze_codons
        cds = "GAA" * 30  # positive in E. coli
        result = analyze_codons(cds, "ecoli")
        assert result.cos > 0

    def test_nterm_flag(self):
        from codon_analyzer import analyze_codons
        # Start codon + 3 negative codons + neutral
        cds = "AUG" + "AGG" * 3 + "GCU" * 27
        result = analyze_codons(cds, "mammalian")
        assert result.nterm_flag is True

    def test_uorf_none(self):
        from codon_analyzer import analyze_uorfs
        result = analyze_uorfs("CCCCCCCCCCCC")
        assert result.total_uaugs == 0

    def test_uorf_found(self):
        from codon_analyzer import analyze_uorfs
        # UTR with one AUG
        utr = "CCCCAUGCCCCUAACCCCC"
        result = analyze_uorfs(utr)
        assert result.total_uaugs == 1

    def test_kozak_strong(self):
        from codon_analyzer import analyze_kozak
        # Perfect Kozak: A at -3, C at -2, C at -1, G at +4, G at +5
        utr = "NNNACCAUG"
        cds = "AUGGGGCCC" + "GCU" * 10
        result = analyze_kozak(utr, cds)
        assert result.kozak_class in ("Strong", "Moderate")
        assert result.kozak_score >= 4

    def test_kozak_weak(self):
        from codon_analyzer import analyze_kozak
        utr = "NNNUUUAUG"
        cds = "AUGUUUCCC" + "GCU" * 10
        result = analyze_kozak(utr, cds)
        assert result.kozak_score < 4

    def test_utr_composition(self):
        from codon_analyzer import analyze_utr_composition
        utr = "AAAAAUUUUU"
        result = analyze_utr_composition(utr)
        assert result.au_content == 100.0
        assert result.gc_content == 0.0
        assert result.utr_length == 10

    def test_compute_lengths(self):
        from codon_analyzer import compute_lengths
        result = compute_lengths("A" * 100, "U" * 900)
        assert result.cds_length == 900
        assert result.transcript_length == 1000


# ---- reference_distributions.py unit tests ----

class TestReferenceDistributions:
    def test_percentile_known(self):
        from reference_distributions import get_percentile
        # COS=0 is the median
        p = get_percentile("cos", 0.0)
        assert p == 50.0

    def test_percentile_high(self):
        from reference_distributions import get_percentile
        p = get_percentile("cos", 0.06)
        assert p >= 89.0  # At or near the 90th percentile boundary

    def test_percentile_interpolated(self):
        from reference_distributions import get_percentile
        p = get_percentile("cos", 0.015)
        assert 50 < p < 75

    def test_unknown_metric(self):
        from reference_distributions import get_percentile
        assert get_percentile("nonexistent", 42) is None


# ---- scorer.py unit tests ----

class TestScorer:
    def test_perfect_score(self):
        from scorer import compute_score
        from codon_analyzer import CodonInfo
        result = compute_score(
            aug_accessibility=0.90,
            cos=0.06,
            nterm_cos=0.5,
            nterm_codons=[],
            total_uaugs=0,
            high_impact_uaugs=0,
            kozak_score_val=7,
            au_content=70,
            gg_frequency=1,
            utr_length=80,
            utr_mfe=-2,
        )
        assert result.scores.total_score == 100
        assert result.rating == "High Efficiency"
        assert result.rating_color == "green"

    def test_poor_score(self):
        from scorer import compute_score
        result = compute_score(
            aug_accessibility=0.1,
            cos=-0.10,
            nterm_cos=-0.3,
            nterm_codons=[],
            total_uaugs=5,
            high_impact_uaugs=3,
            kozak_score_val=1,
            au_content=30,
            gg_frequency=8,
            utr_length=700,
            utr_mfe=-60,
        )
        assert result.scores.total_score < 20
        assert result.rating == "Poor Efficiency"

    def test_ires_suppression(self):
        from scorer import compute_score
        result = compute_score(
            aug_accessibility=0.90,
            cos=0.06,
            nterm_cos=0.5,
            nterm_codons=[],
            total_uaugs=0,
            high_impact_uaugs=0,
            kozak_score_val=7,
            au_content=70,
            gg_frequency=1,
            utr_length=80,
            utr_mfe=-2,
            ires_detected=True,
        )
        assert result.scores.aug_accessibility_score == 0
        assert result.scores.uorf_score == 0
        assert result.scores.utr_structure_score == 0

    def test_prioritized_fixes_generated(self):
        from scorer import compute_score
        result = compute_score(
            aug_accessibility=0.2,
            cos=-0.08,
            nterm_cos=-0.3,
            nterm_codons=[],
            total_uaugs=4,
            high_impact_uaugs=2,
            kozak_score_val=2,
            au_content=40,
            gg_frequency=6,
            utr_length=600,
            utr_mfe=-40,
        )
        assert len(result.prioritized_fixes) > 0
        assert result.primary_bottleneck != ""


# ---- structure_analyzer.py unit tests ----

class TestStructureAnalyzer:
    def test_fold_utr(self):
        from structure_analyzer import fold_utr
        result = fold_utr("GGGAAACCC" * 3)
        assert result.mfe < 0
        assert len(result.structure_dotbracket) > 0

    def test_fold_utr_too_short(self):
        from structure_analyzer import fold_utr
        result = fold_utr("AUGC")
        assert result.skipped is True

    def test_aug_accessibility(self):
        from structure_analyzer import compute_aug_accessibility
        result = compute_aug_accessibility("A" * 30 + "AUGC" * 3, "AUG" + "GCU" * 20)
        assert 0 <= result.aug_accessibility <= 1
        assert result.window_sequence != ""

    def test_aug_accessibility_short_utr(self):
        from structure_analyzer import compute_aug_accessibility
        result = compute_aug_accessibility("AUG", "AUG" + "GCU" * 20)
        assert result.skipped is True

    def test_sliding_window(self):
        from structure_analyzer import sliding_window_cds
        cds = "GGGAAACCC" * 20  # 180 nt
        results = sliding_window_cds(cds)
        assert len(results) > 0
        for w in results:
            assert hasattr(w, 'mfe')
            assert hasattr(w, 'is_flagged')


# ---- Integration tests with real files ----

class TestIntegrationDna:
    """Full pipeline integration test with the .dna example file."""

    @pytest.fixture
    def dna_bytes(self):
        with open(DNA_FILE, "rb") as f:
            return f.read()

    def test_full_pipeline(self, dna_bytes):
        from parser import parse_file, select_cds_candidates, extract_sequences
        from sequence_screener import run_screens
        from structure_analyzer import fold_utr, compute_aug_accessibility, sliding_window_cds
        from codon_analyzer import analyze_codons, analyze_uorfs, analyze_kozak, analyze_utr_composition, compute_lengths
        from reference_distributions import get_percentile
        from scorer import compute_score

        # Parse
        record, features, is_circular, seq = parse_file(dna_bytes, "test.dna")
        assert len(seq) > 100

        # Select CDS
        result = select_cds_candidates(features)
        assert len(result.candidates) > 0
        idx = result.candidates[0].index

        # Extract
        ext = extract_sequences(dna_bytes, idx, "test.dna")
        assert len(ext.cds_sequence) > 10
        assert len(ext.utr_sequence) >= 0

        # Screen
        labels = []
        for f in features:
            for k in ("label", "gene", "product", "note"):
                labels.extend(f.qualifiers.get(k, []))
        screens = run_screens(ext.utr_sequence, ext.cds_sequence, labels)

        # Structure (skip sliding window to keep test fast)
        utr_fold = fold_utr(ext.utr_sequence)
        aug = compute_aug_accessibility(ext.utr_sequence, ext.cds_sequence)

        # Codons
        codons = analyze_codons(ext.cds_sequence, "mammalian")
        assert -1 <= codons.cos <= 1

        uorfs = analyze_uorfs(ext.utr_sequence)
        kozak = analyze_kozak(ext.utr_sequence, ext.cds_sequence)
        utr_comp = analyze_utr_composition(ext.utr_sequence)
        lengths = compute_lengths(ext.utr_sequence, ext.cds_sequence)

        # Percentiles
        cos_pct = get_percentile("cos", codons.cos)
        assert cos_pct is not None

        # Score
        score = compute_score(
            aug_accessibility=aug.aug_accessibility if not aug.skipped else 0.5,
            cos=codons.cos,
            nterm_cos=codons.nterm_cos,
            nterm_codons=codons.nterm_codons,
            total_uaugs=uorfs.total_uaugs,
            high_impact_uaugs=uorfs.high_impact_uaugs,
            kozak_score_val=kozak.kozak_score,
            au_content=utr_comp.au_content,
            gg_frequency=utr_comp.gg_frequency,
            utr_length=utr_comp.utr_length,
            utr_mfe=utr_fold.mfe if not utr_fold.skipped else 0,
            ires_detected=screens.ires_detected,
            aug_skipped=aug.skipped,
        )
        assert 0 <= score.scores.total_score <= 100
        assert score.rating in ("High Efficiency", "Moderate Efficiency", "Low Efficiency", "Poor Efficiency")


class TestIntegrationGbk:
    """Full pipeline integration test with the .gbk example file."""

    @pytest.fixture
    def gbk_bytes(self):
        with open(GBK_FILE, "rb") as f:
            return f.read()

    def test_full_pipeline(self, gbk_bytes):
        from parser import parse_file, select_cds_candidates, extract_sequences
        from sequence_screener import run_screens
        from codon_analyzer import analyze_codons, analyze_kozak
        from scorer import compute_score
        from structure_analyzer import fold_utr, compute_aug_accessibility
        from codon_analyzer import analyze_uorfs, analyze_utr_composition, compute_lengths

        # Parse
        record, features, is_circular, seq = parse_file(gbk_bytes, "test.gbk")
        assert len(seq) > 100

        # Select CDS
        result = select_cds_candidates(features)
        assert len(result.candidates) > 0
        idx = result.candidates[0].index

        # Extract
        ext = extract_sequences(gbk_bytes, idx, "test.gbk")
        assert len(ext.cds_sequence) > 10
        # Verify RNA alphabet
        bad_chars = set(ext.cds_sequence) - set("AUGC")
        assert bad_chars == set(), f"Unexpected chars in CDS RNA: {bad_chars}"

        # Run full analysis
        labels = []
        for f in features:
            for k in ("label", "gene", "product", "note"):
                labels.extend(f.qualifiers.get(k, []))
        screens = run_screens(ext.utr_sequence, ext.cds_sequence, labels)
        utr_fold = fold_utr(ext.utr_sequence)
        aug = compute_aug_accessibility(ext.utr_sequence, ext.cds_sequence)
        codons = analyze_codons(ext.cds_sequence, "mammalian")
        uorfs = analyze_uorfs(ext.utr_sequence)
        kozak = analyze_kozak(ext.utr_sequence, ext.cds_sequence)
        utr_comp = analyze_utr_composition(ext.utr_sequence)
        lengths = compute_lengths(ext.utr_sequence, ext.cds_sequence)

        score = compute_score(
            aug_accessibility=aug.aug_accessibility if not aug.skipped else 0.5,
            cos=codons.cos,
            nterm_cos=codons.nterm_cos,
            nterm_codons=codons.nterm_codons,
            total_uaugs=uorfs.total_uaugs,
            high_impact_uaugs=uorfs.high_impact_uaugs,
            kozak_score_val=kozak.kozak_score,
            au_content=utr_comp.au_content,
            gg_frequency=utr_comp.gg_frequency,
            utr_length=utr_comp.utr_length,
            utr_mfe=utr_fold.mfe if not utr_fold.skipped else 0,
            ires_detected=screens.ires_detected,
            aug_skipped=aug.skipped,
        )
        assert 0 <= score.scores.total_score <= 100


# ---- FastAPI integration tests ----

class TestAPI:
    @pytest.fixture
    def client(self):
        from fastapi.testclient import TestClient
        from main import app
        return TestClient(app)

    def test_health(self, client):
        r = client.get("/health")
        assert r.status_code == 200
        assert r.json()["status"] == "ok"

    def test_parse_dna(self, client):
        with open(DNA_FILE, "rb") as f:
            r = client.post("/parse", files={"file": ("test.dna", f, "application/octet-stream")})
        assert r.status_code == 200
        data = r.json()
        assert "candidates" in data
        assert len(data["candidates"]) > 0

    def test_parse_gbk(self, client):
        with open(GBK_FILE, "rb") as f:
            r = client.post("/parse", files={"file": ("test.gbk", f, "application/octet-stream")})
        assert r.status_code == 200
        data = r.json()
        assert "candidates" in data
        assert len(data["candidates"]) > 0

    def test_parse_bad_extension(self, client):
        r = client.post("/parse", files={"file": ("test.txt", b"hello", "text/plain")})
        assert r.status_code == 400

    def test_analyze_dna(self, client):
        # First parse
        with open(DNA_FILE, "rb") as f:
            r = client.post("/parse", files={"file": ("test.dna", f, "application/octet-stream")})
        data = r.json()
        idx = data["candidates"][0]["index"]

        # Then analyze
        r = client.post("/analyze", data={
            "filename": data["filename"],
            "selected_cds_index": idx,
            "expression_system": "mammalian",
        })
        assert r.status_code == 200
        result = r.json()
        assert "scores" in result
        assert "rating" in result
        assert 0 <= result["scores"]["total_score"] <= 100
        assert result["gene_name"]
        assert result["codons"]["cos"] is not None

    def test_analyze_gbk(self, client):
        with open(GBK_FILE, "rb") as f:
            r = client.post("/parse", files={"file": ("test.gbk", f, "application/octet-stream")})
        data = r.json()
        idx = data["candidates"][0]["index"]

        r = client.post("/analyze", data={
            "filename": data["filename"],
            "selected_cds_index": idx,
            "expression_system": "mammalian",
        })
        assert r.status_code == 200
        result = r.json()
        assert "scores" in result
        assert 0 <= result["scores"]["total_score"] <= 100
