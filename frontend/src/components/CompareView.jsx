import React from 'react'
import ScoreCard from './ScoreCard'
import CdsWindowChart from './CdsWindowChart'
import CodonHeatmap from './CodonHeatmap'
import UorfMap from './UorfMap'
import RatingBadge from './RatingBadge'

const SYSTEM_LABELS = {
  mammalian: 'Mammalian', ecoli: 'E. coli', yeast: 'Yeast', insect: 'Insect (Sf9)',
}

function DeltaBadge({ val1, val2, max }) {
  const delta = val1 - val2
  if (Math.abs(delta) < 5) return null
  return (
    <span className={`text-xs px-1 rounded ${delta > 0 ? 'bg-green-100 text-green-700' : 'bg-red-100 text-red-700'}`}>
      {delta > 0 ? '+' : ''}{delta}
    </span>
  )
}

export default function CompareView({ results1, results2, onReset, onExport }) {
  const r1 = results1
  const r2 = results2

  return (
    <div className="min-h-screen bg-gray-50">
      <div className="sticky top-0 z-10 bg-white border-b shadow-sm px-6 py-3 flex items-center justify-between no-print">
        <div className="flex items-center gap-3">
          <span className="font-bold text-lg">Compare: {r1.gene_name} vs {r2.gene_name}</span>
          <span className="bg-blue-100 text-blue-700 text-xs px-2 py-0.5 rounded-full font-medium">
            {SYSTEM_LABELS[r1.expression_system] || r1.expression_system}
          </span>
        </div>
        <div className="flex gap-2">
          <button onClick={onExport} className="px-3 py-1.5 text-sm border rounded-lg hover:bg-gray-50">Export PDF</button>
          <button onClick={onReset} className="px-3 py-1.5 text-sm bg-blue-600 text-white rounded-lg hover:bg-blue-700">Analyze Another</button>
        </div>
      </div>

      <div className="grid grid-cols-2 gap-4 p-4 max-w-7xl mx-auto">
        {[r1, r2].map((d, col) => (
          <div key={col} className="space-y-4 overflow-y-auto">
            {/* Score header */}
            <div className="bg-white rounded-xl shadow p-6 text-center">
              <p className="text-sm text-gray-500 mb-1">{d.gene_name}</p>
              <div className="text-5xl font-bold">
                {d.scores.total_score}<span className="text-xl text-gray-400"> / 100</span>
              </div>
              <p className="text-xs text-gray-500 mt-1">Heuristic Translation Efficiency Score</p>
              <div className="mt-2">
                <RatingBadge rating={d.rating} color={d.rating_color} />
              </div>
              {col === 1 && (
                <DeltaBadge val1={r1.scores.total_score} val2={r2.scores.total_score} />
              )}
            </div>

            {/* Score cards */}
            <div className="grid grid-cols-2 gap-3">
              <ScoreCard title="AUG Accessibility" score={d.scores.aug_accessibility_score} max={25}
                value={d.screens.ires_detected ? 'N/A' : `${(d.structure.aug_accessibility * 100).toFixed(0)}%`}
                disabled={d.screens.ires_detected} confidence="High confidence" />
              <ScoreCard title="Codon Optimality" score={d.scores.codon_optimality_score} max={20}
                value={`COS: ${d.codons.cos >= 0 ? '+' : ''}${d.codons.cos.toFixed(3)}`}
                confidence={SYSTEM_LABELS[d.expression_system]} />
              <ScoreCard title="N-terminal Codons" score={d.scores.nterm_codon_score} max={10}
                value={`nCOS: ${d.codons.nterm_cos >= 0 ? '+' : ''}${d.codons.nterm_cos.toFixed(3)}`}
                confidence="First 10 codons" />
              <ScoreCard title="uORF Burden" score={d.scores.uorf_score} max={15}
                value={d.screens.ires_detected ? 'N/A' : `${d.uorfs.total_uaugs} uAUGs`}
                disabled={d.screens.ires_detected} confidence="Moderate confidence" />
              <ScoreCard title="Kozak Context" score={d.scores.kozak_score} max={10}
                value={d.kozak.kozak_sequence} confidence="High confidence" />
              <ScoreCard title="5' UTR Structure" score={d.scores.utr_structure_score} max={10}
                value={d.screens.ires_detected ? 'N/A' : `\u0394G = ${d.structure.utr_mfe.toFixed(1)}`}
                disabled={d.screens.ires_detected} confidence="Moderate confidence" />
            </div>

            {/* Fixes */}
            {d.prioritized_fixes.length > 0 && (
              <div className="bg-white rounded-xl shadow p-4">
                <h3 className="font-bold text-sm mb-2">Suggested Improvements</h3>
                <ol className="space-y-2">
                  {d.prioritized_fixes.map((fix, i) => (
                    <li key={i} className="text-xs text-gray-700 flex gap-2">
                      <span className="bg-blue-100 text-blue-700 rounded-full w-5 h-5 flex items-center justify-center text-[10px] font-bold flex-shrink-0">{i+1}</span>
                      {fix}
                    </li>
                  ))}
                </ol>
              </div>
            )}

            {/* CDS chart */}
            <div className="bg-white rounded-xl shadow p-4">
              <h3 className="font-bold text-sm mb-2">CDS Structure</h3>
              <CdsWindowChart windows={d.structure.cds_windows}
                tagRegion={d.screens.tag_detected ? d.screens.tag_length_codons * 3 : 0} />
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}
