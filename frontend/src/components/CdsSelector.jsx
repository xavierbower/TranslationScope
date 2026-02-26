import React, { useState } from 'react'

export default function CdsSelector({ parseResult, parseResult2, compareMode, onSelect, onBack }) {
  const [selected1, setSelected1] = useState(null)
  const [selected2, setSelected2] = useState(null)
  const [showFiltered, setShowFiltered] = useState(false)

  const canProceed = selected1 !== null && (!compareMode || !parseResult2 || selected2 !== null || parseResult2.auto_selected !== null)

  function handleProceed() {
    const idx2 = parseResult2?.auto_selected ?? selected2
    onSelect(selected1, idx2)
  }

  return (
    <div>
      <p className="text-sm text-gray-600 mb-4">
        Multiple coding sequences were found. Select the one you want to analyze.
      </p>

      {parseResult.warnings.map((w, i) => (
        <div key={i} className="bg-yellow-50 text-yellow-700 text-xs p-2 rounded mb-2">{w}</div>
      ))}

      <CdsTable
        label={compareMode ? `File 1: ${parseResult.filename}` : null}
        candidates={parseResult.candidates}
        filteredCandidates={parseResult.filtered_candidates}
        showFiltered={showFiltered}
        onToggleFiltered={() => setShowFiltered(!showFiltered)}
        selected={selected1}
        onSelect={setSelected1}
      />

      {compareMode && parseResult2 && parseResult2.auto_selected === null && (
        <div className="mt-6">
          <CdsTable
            label={`File 2: ${parseResult2.filename}`}
            candidates={parseResult2.candidates}
            filteredCandidates={parseResult2.filtered_candidates}
            showFiltered={showFiltered}
            onToggleFiltered={() => setShowFiltered(!showFiltered)}
            selected={selected2}
            onSelect={setSelected2}
          />
        </div>
      )}

      <div className="flex gap-2 mt-6">
        <button onClick={onBack} className="px-4 py-2 text-sm border rounded-lg hover:bg-gray-50">Back</button>
        <button
          onClick={handleProceed}
          disabled={!canProceed}
          className="px-4 py-2 text-sm bg-blue-600 text-white rounded-lg hover:bg-blue-700 disabled:bg-gray-300"
        >
          Analyze Selected
        </button>
      </div>
    </div>
  )
}

function CdsTable({ label, candidates, filteredCandidates, showFiltered, onToggleFiltered, selected, onSelect }) {
  return (
    <div>
      {label && <h3 className="font-medium text-sm mb-2">{label}</h3>}
      <table className="w-full text-sm">
        <thead>
          <tr className="text-left text-gray-500 border-b">
            <th className="py-2">Feature Label</th>
            <th className="py-2">Length (aa)</th>
            <th className="py-2">Position</th>
            <th className="py-2"></th>
          </tr>
        </thead>
        <tbody>
          {candidates.map((c, i) => (
            <tr key={i} className={`border-b ${selected === c.index ? 'bg-blue-50' : 'hover:bg-gray-50'}`}>
              <td className="py-2">{c.label}</td>
              <td className="py-2">{c.length_aa}</td>
              <td className="py-2 text-gray-500">{c.start}..{c.end} ({c.strand === 1 ? '+' : '-'})</td>
              <td className="py-2">
                <button
                  onClick={() => onSelect(c.index)}
                  className={`px-3 py-1 text-xs rounded ${
                    selected === c.index ? 'bg-blue-600 text-white' : 'bg-gray-100 hover:bg-gray-200'
                  }`}
                >
                  {selected === c.index ? 'Selected' : 'Select'}
                </button>
              </td>
            </tr>
          ))}
        </tbody>
      </table>

      {filteredCandidates && filteredCandidates.length > 0 && (
        <div className="mt-2">
          <button onClick={onToggleFiltered} className="text-xs text-gray-400 hover:text-gray-600">
            Filtered out: {filteredCandidates.map(c => c.label).join(', ')} [{showFiltered ? 'hide' : 'show anyway'}]
          </button>
          {showFiltered && (
            <table className="w-full text-sm mt-2 opacity-60">
              <tbody>
                {filteredCandidates.map((c, i) => (
                  <tr key={i} className="border-b">
                    <td className="py-2">{c.label} <span className="text-xs text-gray-400">(filtered)</span></td>
                    <td className="py-2">{c.length_aa}</td>
                    <td className="py-2 text-gray-500">{c.start}..{c.end}</td>
                    <td className="py-2">
                      <button
                        onClick={() => onSelect(c.index)}
                        className={`px-3 py-1 text-xs rounded ${
                          selected === c.index ? 'bg-blue-600 text-white' : 'bg-gray-100 hover:bg-gray-200'
                        }`}
                      >
                        {selected === c.index ? 'Selected' : 'Select'}
                      </button>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          )}
        </div>
      )}
    </div>
  )
}
