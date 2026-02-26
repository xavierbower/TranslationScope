import React from 'react'

export default function UorfMap({ uorfs, utrLength }) {
  if (uorfs.total_uaugs === 0) {
    return (
      <div className="flex items-center gap-2 text-green-600">
        <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
        </svg>
        <span className="text-sm font-medium">No upstream AUGs detected</span>
      </div>
    )
  }

  const scale = utrLength > 0 ? 100 / utrLength : 1

  return (
    <div>
      <div className="relative h-12 bg-gray-100 rounded-lg overflow-hidden">
        {/* UTR bar */}
        <div className="absolute inset-0 flex items-center">
          <div className="w-full h-4 bg-blue-100 rounded mx-2 relative">
            <span className="absolute left-1 text-[10px] text-blue-500 top-0.5">5' UTR</span>
            <div className="absolute right-0 top-0 h-full w-1 bg-blue-600" title="CDS start (AUG)" />

            {/* uAUG markers */}
            {uorfs.uaug_list.map((u, i) => {
              const left = u.position * scale
              const color = u.impact === 'high' ? '#ef4444' : '#f97316'
              return (
                <div
                  key={i}
                  className="absolute"
                  style={{ left: `${Math.min(left, 95)}%`, top: '-6px' }}
                  title={`uAUG at ${u.position} nt, ${u.distance_to_cds} nt from CDS, ${u.type} (${u.impact} impact)`}
                >
                  <svg width="10" height="12" viewBox="0 0 10 12">
                    <polygon points="5,0 10,12 0,12" fill={color} />
                  </svg>
                </div>
              )
            })}
          </div>
        </div>
      </div>

      {/* Legend */}
      <div className="flex gap-4 mt-2 text-xs text-gray-500">
        <span className="flex items-center gap-1">
          <svg width="8" height="10"><polygon points="4,0 8,10 0,10" fill="#ef4444" /></svg>
          High impact
        </span>
        <span className="flex items-center gap-1">
          <svg width="8" height="10"><polygon points="4,0 8,10 0,10" fill="#f97316" /></svg>
          Moderate impact
        </span>
        <span className="flex items-center gap-1">
          <span className="w-1 h-3 bg-blue-600 inline-block" />
          CDS start
        </span>
      </div>

      {/* Details */}
      <div className="mt-3">
        <table className="text-xs w-full">
          <thead>
            <tr className="text-gray-400">
              <th className="text-left py-1">Position</th>
              <th className="text-left py-1">Distance to CDS</th>
              <th className="text-left py-1">Type</th>
              <th className="text-left py-1">Impact</th>
            </tr>
          </thead>
          <tbody>
            {uorfs.uaug_list.map((u, i) => (
              <tr key={i} className="border-t border-gray-100">
                <td className="py-1">{u.position} nt</td>
                <td className="py-1">{u.distance_to_cds} nt</td>
                <td className="py-1">{u.type === 'overlapping' ? 'Overlapping (oORF)' : 'Contained uORF'}</td>
                <td className="py-1">
                  <span className={`px-1.5 py-0.5 rounded text-[10px] font-medium ${
                    u.impact === 'high' ? 'bg-red-100 text-red-700' : 'bg-orange-100 text-orange-700'
                  }`}>
                    {u.impact}
                  </span>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  )
}
