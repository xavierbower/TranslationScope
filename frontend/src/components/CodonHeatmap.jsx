import React from 'react'

const COLORS = {
  positive: 'bg-green-200 text-green-800',
  negative: 'bg-red-200 text-red-800',
  neutral: 'bg-gray-100 text-gray-600',
}

export default function CodonHeatmap({ codons, ntermCos, tagLength }) {
  if (!codons || codons.length === 0) {
    return <p className="text-sm text-gray-400">No codon data.</p>
  }

  // Show first 100 codons in detail, then summary
  const displayCodons = codons.slice(0, 100)
  const remaining = codons.length - 100

  return (
    <div>
      <div className="mb-3 flex gap-4 items-center">
        <span className="text-sm font-medium">
          N-terminal COS: <span className={ntermCos > 0 ? 'text-green-600' : ntermCos < -0.1 ? 'text-red-600' : 'text-gray-600'}>
            {ntermCos >= 0 ? '+' : ''}{ntermCos.toFixed(3)}
          </span>
        </span>
        <div className="flex gap-2 text-xs">
          <span className="flex items-center gap-1"><span className="w-3 h-3 rounded bg-green-200 inline-block" /> Positive</span>
          <span className="flex items-center gap-1"><span className="w-3 h-3 rounded bg-red-200 inline-block" /> Negative</span>
          <span className="flex items-center gap-1"><span className="w-3 h-3 rounded bg-gray-100 inline-block" /> Neutral</span>
        </div>
      </div>

      <div className="flex flex-wrap gap-0.5">
        {displayCodons.map((c, i) => {
          const isNterm = i >= 1 && i <= 10
          const isTag = tagLength > 0 && i < tagLength
          let borderClass = ''
          if (isNterm) borderClass = 'ring-2 ring-blue-400'
          if (isTag) borderClass = 'ring-2 ring-dashed ring-purple-400'

          return (
            <div
              key={i}
              className={`w-10 h-10 flex flex-col items-center justify-center rounded text-[9px] leading-tight
                         ${COLORS[c.category]} ${borderClass}`}
              title={`Position ${i + 1}: ${c.codon} (${c.aa}) - ${c.category}`}
            >
              <span className="font-mono font-medium">{c.codon}</span>
              <span className="opacity-70">{c.aa}</span>
            </div>
          )
        })}
      </div>

      {remaining > 0 && (
        <p className="text-xs text-gray-400 mt-2">+ {remaining} more codons (not shown)</p>
      )}

      <div className="flex gap-3 mt-2 text-xs text-gray-500">
        <span className="flex items-center gap-1"><span className="w-3 h-3 rounded ring-2 ring-blue-400 inline-block" /> N-terminal window (codons 2-11)</span>
        {tagLength > 0 && (
          <span className="flex items-center gap-1"><span className="w-3 h-3 rounded ring-2 ring-dashed ring-purple-400 inline-block" /> Tag region</span>
        )}
      </div>
    </div>
  )
}
