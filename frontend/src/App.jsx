import React, { useState } from 'react'
import UploadZone from './components/UploadZone'
import SystemSelector from './components/SystemSelector'
import CdsSelector from './components/CdsSelector'
import ScoreCard from './components/ScoreCard'
import CdsWindowChart from './components/CdsWindowChart'
import CodonHeatmap from './components/CodonHeatmap'
import UorfMap from './components/UorfMap'
import CompareView from './components/CompareView'
import RatingBadge from './components/RatingBadge'

const API = import.meta.env.DEV ? 'http://localhost:8000' : ''

async function safeFetch(url, opts) {
  let res
  try {
    res = await fetch(url, opts)
  } catch (e) {
    throw new Error('Network error — the server may be starting up. Please try again in a few seconds.')
  }
  if (!res.ok) {
    let detail
    try {
      const body = await res.json()
      detail = body.detail
    } catch { /* response wasn't JSON */ }
    throw new Error(detail || `Request failed (HTTP ${res.status})`)
  }
  return res.json()
}

const STEPS = [
  'Parsing file...', 'Screening sequences...', 'Folding RNA...',
  'Scoring codons...', 'Analyzing uORFs...', 'Calculating scores...',
]

export default function App() {
  const [state, setState] = useState('upload') // upload | cds_select | loading | results
  const [file1, setFile1] = useState(null)
  const [file2, setFile2] = useState(null)
  const [compareMode, setCompareMode] = useState(false)
  const [system, setSystem] = useState('mammalian')
  const [parseResult, setParseResult] = useState(null)
  const [parseResult2, setParseResult2] = useState(null)
  const [results, setResults] = useState(null)
  const [results2, setResults2] = useState(null)
  const [error, setError] = useState('')
  const [loadingStep, setLoadingStep] = useState(0)

  async function handleAnalyze() {
    if (!file1) return
    setError('')
    setState('loading')

    try {
      // Parse file 1
      setLoadingStep(0)
      const form1 = new FormData()
      form1.append('file', file1)
      const parsed = await safeFetch(`${API}/parse`, { method: 'POST', body: form1 })

      // Parse file 2 if compare mode
      let parsed2 = null
      if (compareMode && file2) {
        const form2 = new FormData()
        form2.append('file', file2)
        parsed2 = await safeFetch(`${API}/parse`, { method: 'POST', body: form2 })
      }

      // Check if CDS selection needed
      if (parsed.auto_selected === null || (parsed2 && parsed2.auto_selected === null)) {
        setParseResult(parsed)
        setParseResult2(parsed2)
        setState('cds_select')
        return
      }

      // Auto-selected — proceed to analyze
      await runAnalysis(parsed, parsed.auto_selected, parsed2, parsed2?.auto_selected)
    } catch (e) {
      setError(e.message)
      setState('upload')
    }
  }

  async function runAnalysis(parsed, cdsIdx, parsed2, cdsIdx2) {
    setLoadingStep(2)
    const form = new FormData()
    form.append('filename', parsed.filename)
    form.append('selected_cds_index', cdsIdx)
    form.append('expression_system', system)

    const data = await safeFetch(`${API}/analyze`, { method: 'POST', body: form })
    setLoadingStep(4)
    setResults(data)

    if (parsed2 && cdsIdx2 !== null && cdsIdx2 !== undefined) {
      const form2 = new FormData()
      form2.append('filename', parsed2.filename)
      form2.append('selected_cds_index', cdsIdx2)
      form2.append('expression_system', system)
      try {
        const data2 = await safeFetch(`${API}/analyze`, { method: 'POST', body: form2 })
        setResults2(data2)
      } catch { /* second file analysis is best-effort */ }
    }

    setLoadingStep(5)
    setState('results')
  }

  async function handleCdsSelect(idx, idx2) {
    setState('loading')
    try {
      await runAnalysis(parseResult, idx, parseResult2, idx2)
    } catch (e) {
      setError(e.message)
      setState('upload')
    }
  }

  function handleReset() {
    setState('upload')
    setFile1(null)
    setFile2(null)
    setResults(null)
    setResults2(null)
    setParseResult(null)
    setParseResult2(null)
    setError('')
    setCompareMode(false)
  }

  function handleExportPdf() {
    window.print()
  }

  // ---- Upload/Configure ----
  if (state === 'upload') {
    return (
      <div className="min-h-screen flex items-center justify-center p-4">
        <div className="bg-white rounded-2xl shadow-lg p-8 w-full max-w-lg">
          <h1 className="text-2xl font-bold text-center mb-1">TranslationScope</h1>
          <p className="text-gray-500 text-center text-sm mb-6">mRNA Translation Efficiency Analyzer</p>

          {error && (
            <div className="bg-red-50 text-red-700 p-3 rounded-lg mb-4 text-sm">{error}</div>
          )}

          <UploadZone file={file1} onFile={setFile1} label="Upload .dna file" />

          <div className="mt-5">
            <SystemSelector value={system} onChange={setSystem} />
          </div>

          <label className="flex items-center gap-2 mt-4 text-sm text-gray-600 cursor-pointer">
            <input
              type="checkbox"
              checked={compareMode}
              onChange={e => { setCompareMode(e.target.checked); if (!e.target.checked) setFile2(null) }}
              className="rounded"
            />
            Compare two constructs side by side
          </label>

          {compareMode && (
            <div className="mt-3">
              <UploadZone file={file2} onFile={setFile2} label="Upload second .dna file" />
            </div>
          )}

          <button
            onClick={handleAnalyze}
            disabled={!file1 || (compareMode && !file2)}
            className="mt-6 w-full bg-blue-600 text-white py-2.5 rounded-lg font-medium
                       hover:bg-blue-700 disabled:bg-gray-300 disabled:cursor-not-allowed transition"
          >
            Analyze
          </button>

          <div className="mt-6 flex items-center justify-center gap-3 text-xs text-gray-400">
            <a href="https://xavbio.com" target="_blank" rel="noopener noreferrer" className="hover:text-gray-600 transition">xavbio.com</a>
            <span>·</span>
            <a href="https://buymeacoffee.com/xavbio" target="_blank" rel="noopener noreferrer" className="hover:text-gray-600 transition">Buy me a coffee ☕</a>
          </div>
        </div>
      </div>
    )
  }

  // ---- CDS Selection ----
  if (state === 'cds_select') {
    return (
      <div className="min-h-screen flex items-center justify-center p-4">
        <div className="bg-white rounded-2xl shadow-lg p-8 w-full max-w-2xl">
          <h2 className="text-xl font-bold mb-2">Select Coding Sequence</h2>
          <CdsSelector
            parseResult={parseResult}
            parseResult2={parseResult2}
            compareMode={compareMode}
            onSelect={handleCdsSelect}
            onBack={() => setState('upload')}
          />
        </div>
      </div>
    )
  }

  // ---- Loading ----
  if (state === 'loading') {
    return (
      <div className="min-h-screen flex items-center justify-center p-4">
        <div className="bg-white rounded-2xl shadow-lg p-8 w-full max-w-md text-center">
          <div className="animate-spin h-10 w-10 border-4 border-blue-500 border-t-transparent rounded-full mx-auto mb-4" />
          <div className="space-y-1">
            {STEPS.map((step, i) => (
              <p key={i} className={`text-sm ${i <= loadingStep ? 'text-blue-600 font-medium' : 'text-gray-400'}`}>
                {i < loadingStep ? '\u2713 ' : i === loadingStep ? '\u25B6 ' : ''}{step}
              </p>
            ))}
          </div>
          <p className="text-xs text-gray-400 mt-4">Estimated time: ~5-15 seconds</p>
        </div>
      </div>
    )
  }

  // ---- Results ----
  if (state === 'results' && results) {
    if (compareMode && results2) {
      return <CompareView results1={results} results2={results2} onReset={handleReset} onExport={handleExportPdf} />
    }

    return <ResultsView data={results} onReset={handleReset} onExport={handleExportPdf} />
  }

  return null
}

function ResultsView({ data, onReset, onExport }) {
  const d = data
  const ires = d.screens.ires_detected

  const SYSTEM_LABELS = {
    mammalian: 'Mammalian', ecoli: 'E. coli', yeast: 'Yeast', insect: 'Insect (Sf9)',
  }

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Sticky top bar */}
      <div className="sticky top-0 z-10 bg-white border-b shadow-sm px-6 py-3 flex items-center justify-between no-print">
        <div className="flex items-center gap-3">
          <span className="font-bold text-lg">{d.gene_name || 'Unknown Gene'}</span>
          <span className="bg-blue-100 text-blue-700 text-xs px-2 py-0.5 rounded-full font-medium">
            {SYSTEM_LABELS[d.expression_system] || d.expression_system}
          </span>
        </div>
        <div className="flex gap-2">
          <button onClick={onExport} className="px-3 py-1.5 text-sm border rounded-lg hover:bg-gray-50">Export PDF</button>
          <button onClick={onReset} className="px-3 py-1.5 text-sm bg-blue-600 text-white rounded-lg hover:bg-blue-700">Analyze Another</button>
        </div>
      </div>

      <div className="max-w-5xl mx-auto p-6 space-y-6">
        {/* Section 1: Score Header */}
        <div className="bg-white rounded-xl shadow p-6 text-center">
          <div className="text-6xl font-bold">{d.scores.total_score}<span className="text-2xl text-gray-400"> / 100</span></div>
          <p className="text-sm text-gray-500 mt-1">Heuristic Translation Efficiency Score</p>
          <div className="flex items-center justify-center gap-3 mt-3">
            <RatingBadge rating={d.rating} color={d.rating_color} />
            <span className={`text-xs px-2 py-0.5 rounded-full font-medium ${
              d.kozak.kozak_class === 'Strong' ? 'bg-green-100 text-green-700' :
              d.kozak.kozak_class === 'Moderate' ? 'bg-yellow-100 text-yellow-700' :
              'bg-red-100 text-red-700'
            }`}>
              {d.kozak.kozak_class} Kozak
            </span>
          </div>
          {d.percentiles.cos_percentile != null && (
            <p className="text-xs text-gray-400 mt-2">
              ~{Math.round(d.percentiles.cos_percentile)}th percentile vs. human transcriptome (approximate)
            </p>
          )}
        </div>

        {/* Warning Banners */}
        {ires && (
          <div className="bg-orange-50 border-l-4 border-orange-400 p-4 rounded-r-lg">
            <p className="text-orange-800 text-sm font-medium">IRES sequence detected or suspected</p>
            <p className="text-orange-700 text-xs mt-1">
              Cap-dependent initiation models (UTR structure, AUG accessibility, uORF scores) do not apply
              to IRES-driven translation and have been hidden. Codon optimality scores remain valid.
            </p>
          </div>
        )}
        {d.screens.signal_peptide_suspected && (
          <div className="bg-yellow-50 border-l-4 border-yellow-400 p-4 rounded-r-lg">
            <p className="text-yellow-800 text-sm">Signal peptide may be present. Scores may be underestimated for secreted/membrane proteins.</p>
          </div>
        )}
        {d.screens.tag_detected && (
          <div className="bg-yellow-50 border-l-4 border-yellow-400 p-4 rounded-r-lg">
            <p className="text-yellow-800 text-sm">
              N-terminal {d.screens.tag_name} tag detected. First ~{d.screens.tag_length_codons} codons reflect tag sequence.
            </p>
          </div>
        )}
        {d.screens.utr_source === 'inferred' && d.warnings.some(w => w.includes('No 5\' UTR')) && (
          <div className="bg-gray-100 border-l-4 border-gray-300 p-4 rounded-r-lg">
            <p className="text-gray-600 text-xs">
              No 5' UTR annotation found. Using up to 300 nt upstream — UTR scores should be interpreted cautiously.
            </p>
          </div>
        )}

        {/* Section 2: Score Breakdown */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          <ScoreCard
            title="AUG Accessibility"
            score={d.scores.aug_accessibility_score} max={25}
            value={ires ? 'N/A' : `${(d.structure.aug_accessibility * 100).toFixed(0)}%`}
            percentile={d.percentiles.aug_accessibility_percentile}
            confidence="High confidence"
            disabled={ires || d.structure.aug_skipped}
          />
          <ScoreCard
            title="Codon Optimality"
            score={d.scores.codon_optimality_score} max={20}
            value={`COS: ${d.codons.cos >= 0 ? '+' : ''}${d.codons.cos.toFixed(3)}`}
            percentile={d.percentiles.cos_percentile}
            confidence={SYSTEM_LABELS[d.expression_system]}
          />
          <ScoreCard
            title="N-terminal Codons"
            score={d.scores.nterm_codon_score} max={10}
            value={`nCOS: ${d.codons.nterm_cos >= 0 ? '+' : ''}${d.codons.nterm_cos.toFixed(3)}`}
            subtitle="First 10 codons"
            confidence="High confidence"
          />
          <ScoreCard
            title="uORF Burden"
            score={d.scores.uorf_score} max={15}
            value={ires ? 'N/A' : `${d.uorfs.total_uaugs} uAUGs (${d.uorfs.high_impact_uaugs} high-impact)`}
            disabled={ires}
            confidence="Moderate confidence"
          />
          <ScoreCard
            title="Kozak Context"
            score={d.scores.kozak_score} max={10}
            value={d.kozak.kozak_sequence}
            confidence="High confidence"
          />
          <ScoreCard
            title="5' UTR Structure"
            score={d.scores.utr_structure_score} max={10}
            value={ires ? 'N/A' : `\u0394G = ${d.structure.utr_mfe.toFixed(1)} kcal/mol`}
            disabled={ires}
            confidence="Moderate confidence"
          />
        </div>

        {/* Section 3: Prioritized Fixes */}
        {d.prioritized_fixes.length > 0 && (
          <div className="bg-white rounded-xl shadow p-6">
            <h3 className="font-bold text-lg mb-3">Suggested Design Improvements</h3>
            <ol className="space-y-3">
              {d.prioritized_fixes.map((fix, i) => (
                <li key={i} className="flex gap-3">
                  <span className="bg-blue-100 text-blue-700 rounded-full w-6 h-6 flex items-center justify-center text-sm font-bold flex-shrink-0">
                    {i + 1}
                  </span>
                  <p className="text-sm text-gray-700">{fix}</p>
                </li>
              ))}
            </ol>
          </div>
        )}

        {/* Section 4: CDS Sliding Window Chart */}
        <div className="bg-white rounded-xl shadow p-6">
          <h3 className="font-bold text-lg mb-3">CDS Structural Complexity</h3>
          <CdsWindowChart
            windows={d.structure.cds_windows}
            tagRegion={d.screens.tag_detected ? d.screens.tag_length_codons * 3 : 0}
          />
          <p className="text-xs text-gray-400 mt-2">
            Flagged regions reflect unusually stable predicted secondary structure under in vitro conditions.
            In-cell helicase activity (eIF4A, DDX3) may resolve many of these. Interpret as exploratory.
          </p>
        </div>

        {/* Section 5: Codon Heatmap */}
        <div className="bg-white rounded-xl shadow p-6">
          <h3 className="font-bold text-lg mb-3">Codon Usage</h3>
          <div className="flex gap-4 text-sm mb-3">
            <span className="text-green-600">Positive: {(d.codons.positive_fraction * 100).toFixed(1)}%</span>
            <span className="text-red-600">Negative: {(d.codons.negative_fraction * 100).toFixed(1)}%</span>
            <span className="text-gray-500">Neutral: {((1 - d.codons.positive_fraction - d.codons.negative_fraction) * 100).toFixed(1)}%</span>
          </div>
          <CodonHeatmap
            codons={d.codons.codon_table}
            ntermCos={d.codons.nterm_cos}
            tagLength={d.screens.tag_detected ? d.screens.tag_length_codons : 0}
          />
          <p className="text-xs text-gray-400 mt-2">
            Codon categories from Zheng et al. (2025, Nature Biotechnology)
          </p>
        </div>

        {/* Section 6: uORF Map */}
        {!ires && (
          <div className="bg-white rounded-xl shadow p-6">
            <h3 className="font-bold text-lg mb-3">Upstream ORF Map</h3>
            <UorfMap uorfs={d.uorfs} utrLength={d.utr_composition.utr_length} />
          </div>
        )}

        {/* Section 7: Summary & Disclaimer */}
        <div className="bg-white rounded-xl shadow p-6">
          <h3 className="font-bold text-lg mb-2">Summary</h3>
          <p className="font-medium text-gray-800">{d.primary_bottleneck}</p>
          <p className="text-sm text-gray-600 mt-2">{d.summary}</p>

          {d.warnings.length > 0 && (
            <div className="bg-yellow-50 rounded-lg p-3 mt-4">
              <p className="text-xs font-medium text-yellow-700 mb-1">Warnings</p>
              <ul className="text-xs text-yellow-700 space-y-1">
                {d.warnings.map((w, i) => <li key={i}>{w}</li>)}
              </ul>
            </div>
          )}

          <div className="mt-4 bg-gray-50 rounded-lg p-4 text-xs text-gray-500 leading-relaxed">
            TranslationScope uses thermodynamic predictions and sequence-composition features to flag
            potential translation inefficiencies. Structure scores reflect in vitro RNA folding and may
            not capture in-cell behavior, where helicases and RNA-binding proteins actively remodel UTR
            conformation. Codon optimality scores are derived from human/mouse ribosome profiling data
            (Zheng et al., 2025) and may not generalize to other organisms or expression systems. This
            tool is intended to guide construct design decisions, not to produce precise TE predictions.
            For state-of-the-art mammalian TE modeling, see RiboNN (github.com/Sanofi-Public/RiboNN).
          </div>
        </div>

        <div className="flex items-center justify-center gap-3 text-xs text-gray-400 mt-6 mb-4 no-print">
          <a href="https://xavbio.com" target="_blank" rel="noopener noreferrer" className="hover:text-gray-600 transition">xavbio.com</a>
          <span>·</span>
          <a href="https://buymeacoffee.com/xavbio" target="_blank" rel="noopener noreferrer" className="hover:text-gray-600 transition">Buy me a coffee ☕</a>
        </div>
      </div>
    </div>
  )
}
