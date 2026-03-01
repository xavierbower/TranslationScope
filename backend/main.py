"""FastAPI backend for TranslationScope."""

from __future__ import annotations

import dataclasses
import logging
import os
import traceback
from pathlib import Path
from typing import Any

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse

from parser import parse_file, select_cds_candidates, extract_sequences
from sequence_screener import run_screens
from structure_analyzer import fold_utr, compute_aug_accessibility, sliding_window_cds
from codon_analyzer import (
    analyze_codons, analyze_uorfs, analyze_kozak,
    analyze_utr_composition, compute_lengths,
)
from reference_distributions import get_percentile
from scorer import compute_score

logger = logging.getLogger("translationscope")

app = FastAPI(title="TranslationScope", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


def _dc_to_dict(obj: Any) -> Any:
    """Recursively convert dataclasses to dicts."""
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return {k: _dc_to_dict(v) for k, v in dataclasses.asdict(obj).items()}
    if isinstance(obj, list):
        return [_dc_to_dict(i) for i in obj]
    if isinstance(obj, dict):
        return {k: _dc_to_dict(v) for k, v in obj.items()}
    return obj


# Store uploaded files temporarily in memory keyed by filename
_file_cache: dict[str, bytes] = {}


@app.get("/health")
async def health():
    return {"status": "ok"}


ALLOWED_EXTENSIONS = {".dna", ".gbk", ".gb", ".genbank"}


@app.post("/parse")
async def parse_dna(file: UploadFile = File(...)):
    """Parse a .dna or .gbk file and return CDS candidates."""
    if not file.filename:
        raise HTTPException(400, "No filename provided.")
    ext = "." + file.filename.lower().rsplit(".", 1)[-1]
    if ext not in ALLOWED_EXTENSIONS:
        raise HTTPException(400, f"Unsupported file type. Accepted: {', '.join(ALLOWED_EXTENSIONS)}")

    contents = await file.read()
    try:
        record, features, is_circular, full_seq = parse_file(contents, file.filename)
        result = select_cds_candidates(features)
    except ValueError as e:
        raise HTTPException(400, str(e))
    except Exception as e:
        logger.error("Unexpected error parsing %s: %s\n%s", file.filename, e, traceback.format_exc())
        raise HTTPException(400, f"Failed to parse file: {e}")

    _file_cache[file.filename] = contents

    # Collect feature labels for IRES screening later
    feature_labels = []
    for f in features:
        for key in ("label", "gene", "product", "note"):
            feature_labels.extend(f.qualifiers.get(key, []))

    return {
        "filename": file.filename,
        "is_circular": is_circular,
        "candidates": [_dc_to_dict(c) for c in result.candidates],
        "filtered_candidates": [_dc_to_dict(c) for c in result.filtered_candidates],
        "auto_selected": result.auto_selected,
        "gene_name": result.gene_name,
        "warnings": result.warnings,
        "feature_labels": feature_labels,
    }


@app.post("/analyze")
async def analyze(
    file: UploadFile = File(None),
    filename: str = Form(None),
    selected_cds_index: int = Form(...),
    expression_system: str = Form("mammalian"),
):
    """Run full analysis pipeline on selected CDS."""
    # Get file bytes
    if file and file.filename:
        contents = await file.read()
        fname = file.filename
        _file_cache[fname] = contents
    elif filename and filename in _file_cache:
        contents = _file_cache[filename]
        fname = filename
    else:
        raise HTTPException(400, "No file provided and no cached file found.")

    # Parse and extract
    try:
        record, features, is_circular, full_seq = parse_file(contents, fname)
        extraction = extract_sequences(contents, selected_cds_index, fname)
    except ValueError as e:
        raise HTTPException(400, str(e))
    except Exception as e:
        logger.error("Unexpected error parsing %s: %s\n%s", fname, e, traceback.format_exc())
        raise HTTPException(400, f"Failed to parse file: {e}")

    try:
        warnings = list(extraction.warnings)

        # Collect feature labels
        feature_labels = []
        for f in features:
            for key in ("label", "gene", "product", "note"):
                feature_labels.extend(f.qualifiers.get(key, []))

        # Run screens
        screens = run_screens(extraction.utr_sequence, extraction.cds_sequence, feature_labels)
        warnings.extend(screens.warnings)

        # Structure analysis
        utr_fold = fold_utr(extraction.utr_sequence)
        aug_result = compute_aug_accessibility(extraction.utr_sequence, extraction.cds_sequence)
        cds_windows = sliding_window_cds(extraction.cds_sequence)

        if aug_result.skipped:
            warnings.append(aug_result.skip_reason)

        # Codon analysis
        codons = analyze_codons(extraction.cds_sequence, expression_system)
        uorfs = analyze_uorfs(extraction.utr_sequence)
        kozak = analyze_kozak(extraction.utr_sequence, extraction.cds_sequence)
        utr_comp = analyze_utr_composition(extraction.utr_sequence)
        lengths = compute_lengths(extraction.utr_sequence, extraction.cds_sequence)

        if codons.nterm_flag:
            warnings.append(
                "Warning: multiple unfavorable codons in the critical N-terminal window (codons 1-10). "
                "This region has disproportionate impact on translation efficiency."
            )

        # Non-mammalian system warning
        if expression_system != "mammalian":
            system_labels = {
                "ecoli": "E. coli",
                "yeast": "Yeast (S. cerevisiae)",
                "insect": "Insect (Sf9)",
            }
            warnings.append(
                f"Codon preferences shown for {system_labels.get(expression_system, expression_system)}. "
                "UTR structure and AUG accessibility scores assume cap-dependent eukaryotic initiation "
                "and may not apply to your expression context."
            )

        # Percentiles
        percentiles = {
            "cos_percentile": get_percentile("cos", codons.cos),
            "utr_length_percentile": get_percentile("utr_length", utr_comp.utr_length),
            "au_content_percentile": get_percentile("au_content", utr_comp.au_content),
            "aug_accessibility_percentile": get_percentile("aug_accessibility", aug_result.aug_accessibility) if not aug_result.skipped else None,
            "cds_length_percentile": get_percentile("cds_length", lengths.cds_length),
        }

        # Scoring
        score_result = compute_score(
            aug_accessibility=aug_result.aug_accessibility,
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
            aug_skipped=aug_result.skipped,
        )
    except Exception as e:
        logger.error("Analysis failed for %s: %s\n%s", fname, e, traceback.format_exc())
        raise HTTPException(400, f"Analysis failed: {e}")

    return {
        "gene_name": extraction.gene_name,
        "expression_system": expression_system,
        "screens": _dc_to_dict(screens),
        "structure": {
            "utr_mfe": utr_fold.mfe,
            "utr_structure_dotbracket": utr_fold.structure_dotbracket,
            "utr_skipped": utr_fold.skipped,
            "utr_skip_reason": utr_fold.skip_reason,
            "aug_mfe": aug_result.mfe,
            "aug_accessibility": aug_result.aug_accessibility,
            "aug_window_sequence": aug_result.window_sequence,
            "aug_structure_dotbracket": aug_result.structure_dotbracket,
            "aug_window_truncated": aug_result.window_truncated,
            "aug_skipped": aug_result.skipped,
            "cds_windows": [_dc_to_dict(w) for w in cds_windows],
        },
        "codons": _dc_to_dict(codons),
        "uorfs": _dc_to_dict(uorfs),
        "kozak": _dc_to_dict(kozak),
        "utr_composition": _dc_to_dict(utr_comp),
        "lengths": _dc_to_dict(lengths),
        "percentiles": percentiles,
        "scores": _dc_to_dict(score_result.scores),
        "rating": score_result.rating,
        "rating_color": score_result.rating_color,
        "primary_bottleneck": score_result.primary_bottleneck,
        "summary": score_result.summary,
        "prioritized_fixes": score_result.prioritized_fixes,
        "warnings": warnings,
    }


# Serve frontend static files in production (when built into ./static)
_static_dir = Path(__file__).parent / "static"
if _static_dir.is_dir():
    app.mount("/assets", StaticFiles(directory=str(_static_dir / "assets")), name="assets")

    @app.get("/{full_path:path}")
    async def serve_spa(full_path: str):
        """Serve the React SPA for any non-API route."""
        file_path = _static_dir / full_path
        if file_path.is_file():
            return FileResponse(file_path)
        return FileResponse(_static_dir / "index.html")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
