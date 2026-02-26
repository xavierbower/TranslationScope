import React from 'react'
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip,
  ReferenceLine, ReferenceArea, ResponsiveContainer, Scatter, ComposedChart,
} from 'recharts'

export default function CdsWindowChart({ windows, tagRegion }) {
  if (!windows || windows.length === 0) {
    return <p className="text-sm text-gray-400">No CDS windows to display.</p>
  }

  const data = windows.map(w => ({
    position: w.position,
    mfe: w.mfe,
    flagged: w.is_flagged ? w.mfe : null,
  }))

  return (
    <ResponsiveContainer width="100%" height={300}>
      <ComposedChart data={data} margin={{ top: 5, right: 20, bottom: 20, left: 20 }}>
        <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
        <XAxis
          dataKey="position"
          label={{ value: 'CDS position (nt)', position: 'bottom', offset: 0, fontSize: 12 }}
          fontSize={11}
        />
        <YAxis
          label={{ value: '\u0394G (kcal/mol)', angle: -90, position: 'insideLeft', fontSize: 12 }}
          fontSize={11}
        />
        <Tooltip
          formatter={(val) => [`${val} kcal/mol`, '\u0394G']}
          labelFormatter={v => `Position: ${v} nt`}
        />
        <ReferenceLine y={-20} stroke="#ef4444" strokeDasharray="5 5" label={{ value: '-20 kcal/mol', fontSize: 10, fill: '#ef4444' }} />
        {tagRegion > 0 && (
          <ReferenceArea x1={0} x2={tagRegion} fill="#bfdbfe" fillOpacity={0.3} label={{ value: 'Tag region', fontSize: 10 }} />
        )}
        <Line type="monotone" dataKey="mfe" stroke="#3b82f6" dot={false} strokeWidth={1.5} />
        <Scatter dataKey="flagged" fill="#ef4444" r={4} />
      </ComposedChart>
    </ResponsiveContainer>
  )
}
