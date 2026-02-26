import React from 'react'

const SYSTEMS = [
  { value: 'mammalian', label: 'Mammalian' },
  { value: 'ecoli', label: 'E. coli' },
  { value: 'yeast', label: 'Yeast (S. cerevisiae)' },
  { value: 'insect', label: 'Insect (Sf9)' },
]

export default function SystemSelector({ value, onChange }) {
  return (
    <div>
      <label className="block text-sm font-medium text-gray-700 mb-2">Expression System</label>
      <div className="grid grid-cols-2 gap-2">
        {SYSTEMS.map(s => (
          <button
            key={s.value}
            onClick={() => onChange(s.value)}
            className={`px-3 py-2 text-sm rounded-lg border transition ${
              value === s.value
                ? 'bg-blue-600 text-white border-blue-600'
                : 'bg-white text-gray-700 border-gray-300 hover:border-blue-400'
            }`}
          >
            {s.label}
          </button>
        ))}
      </div>
      <p className="text-xs text-gray-400 mt-2">
        Codon optimality scoring depends on the expression host. Structural features assume
        cap-dependent eukaryotic initiation.
      </p>
    </div>
  )
}
