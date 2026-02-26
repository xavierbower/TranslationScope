import React from 'react'

export default function ScoreCard({ title, score, max, value, percentile, confidence, subtitle, disabled }) {
  const pct = max > 0 ? (score / max) * 100 : 0
  const barColor = disabled ? 'bg-gray-200' : pct >= 70 ? 'bg-green-500' : pct >= 40 ? 'bg-yellow-500' : 'bg-red-500'

  if (disabled) {
    return (
      <div className="bg-gray-50 rounded-xl p-4 border border-gray-200">
        <p className="text-sm font-medium text-gray-400">{title}</p>
        <p className="text-2xl font-bold text-gray-300 mt-1">N/A</p>
        <p className="text-xs text-gray-400 mt-1">IRES detected</p>
      </div>
    )
  }

  return (
    <div className="bg-white rounded-xl p-4 border border-gray-200 hover:shadow-md transition">
      <p className="text-sm font-medium text-gray-700">{title}</p>
      <div className="flex items-baseline gap-1 mt-1">
        <span className="text-2xl font-bold">{score}</span>
        <span className="text-sm text-gray-400">/ {max}</span>
      </div>
      <div className="w-full bg-gray-100 rounded-full h-1.5 mt-2">
        <div className={`h-1.5 rounded-full ${barColor}`} style={{ width: `${pct}%` }} />
      </div>
      <p className="text-xs text-gray-600 mt-2 font-mono">{value}</p>
      {subtitle && <p className="text-xs text-gray-400">{subtitle}</p>}
      {percentile != null && (
        <p className="text-xs text-gray-400 mt-1">~{Math.round(percentile)}th percentile</p>
      )}
      {confidence && (
        <p className="text-xs text-gray-400 italic mt-1">{confidence}</p>
      )}
    </div>
  )
}
