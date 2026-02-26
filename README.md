# TranslationScope

Analyze mRNA translation efficiency from SnapGene `.dna` files. TranslationScope uses thermodynamic RNA structure predictions and sequence-composition features to identify potential translation bottlenecks and suggest design improvements.

## Setup

### Backend

```bash
cd backend

# Create virtual environment
python -m venv venv
source venv/bin/activate  # or venv\Scripts\activate on Windows

# Install dependencies
pip install -r requirements.txt

# ViennaRNA: if pip install fails, use conda:
# conda install -c bioconda viennarna

# Start server
uvicorn main:app --reload --port 8000
```

### Frontend

```bash
cd frontend
npm install
npm run dev
```

The frontend runs at `http://localhost:5173` and the backend at `http://localhost:8000`.

## Scores and Evidence

### AUG Accessibility (25 pts) — High confidence

Measures unpaired probability of the AUG start codon region using ViennaRNA partition function. Open structure at the start codon is rate-limiting for initiation — mechanistically direct and well-replicated evidence.

### Codon Optimality (20 pts) — System-dependent

Overall Codon Optimality Score (COS) based on the fraction of positive vs. negative codons for the selected expression system. Mammalian tables derived from RiboNN insertional analysis (Zheng et al., 2025).

### N-terminal Codons (10 pts) — High confidence

First 10 codons after start methionine scored separately. RiboNN found ~2x per-residue impact on TE for N-terminal codons compared to mid-CDS codons.

### uORF Burden (13 pts) — Moderate confidence

Upstream AUGs in the 5' UTR that may compete with the main CDS for ribosome initiation. Overlapping ORFs (oORFs) and close-proximity uAUGs classified as high impact.

### Kozak Context (8 pts) — High confidence

Sequence context around the start codon scored against consensus. Positions -3 (A/G) and +4 (G) are most critical.

### 5' UTR Composition (14 pts) — Empirically calibrated

AU content, GG dinucleotide frequency, and UTR length. AU-rich, shorter UTRs correlate with higher TE.

### 5' UTR Structure (10 pts) — Moderate confidence

Global minimum free energy of the 5' UTR. Note: these are in vitro thermodynamic predictions. In-cell helicases (eIF4A, DDX3) actively unwind structures during scanning.

## Expression System Codon Tables

| System | Source |
|--------|--------|
| Mammalian | RiboNN insertional analysis (Zheng et al., 2025) |
| E. coli | Highly expressed E. coli gene codon usage |
| Yeast | Codon adaptation index for highly expressed yeast genes |
| Insect (Sf9) | Spodoptera frugiperda codon usage tables |

## Pre-analysis Screens

- **IRES Detection**: Annotated IRES features or heuristic (UTR > 400 nt with GC < 45%). Suppresses cap-dependent initiation scores.
- **Signal Peptide**: N-terminal hydrophobic stretch detection. Flags but does not suppress scores.
- **N-terminal Tags**: Detects 6xHis, FLAG, Strep, HA, Myc, MBP, GST tags in the first 30 codons.

## Citation

Zheng et al. (2025) "Predicting the translation efficiency of messenger RNA in mammalian cells." *Nature Biotechnology*. https://doi.org/10.1038/s41587-025-02712-x

## Known Limitations

- Analyzes only one CDS per file (user-selected if ambiguous)
- Structure scores are in vitro thermodynamic predictions
- Codon preferences derived from human/mouse data; non-mammalian tables are approximate
- N-terminal tag detection uses simple heuristics; verify manually
- IRES detection is heuristic; verify by annotation
- Percentile values are approximate, based on published summary statistics rather than a live database query
- Tool is heuristic; not a validated regression model
- Does not model dicodon effects, RNA modifications (m6A, m1Ψ), or RNA-binding protein interactions
- For state-of-the-art mammalian TE modeling, see [RiboNN](https://github.com/Sanofi-Public/RiboNN)
