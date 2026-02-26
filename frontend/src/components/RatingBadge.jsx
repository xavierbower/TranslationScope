import React from 'react'

const COLOR_MAP = {
  green: 'bg-green-100 text-green-700 border-green-300',
  yellow: 'bg-yellow-100 text-yellow-700 border-yellow-300',
  orange: 'bg-orange-100 text-orange-700 border-orange-300',
  red: 'bg-red-100 text-red-700 border-red-300',
}

export default function RatingBadge({ rating, color }) {
  const classes = COLOR_MAP[color] || COLOR_MAP.yellow
  return (
    <span className={`px-3 py-1 rounded-full text-sm font-medium border ${classes}`}>
      {rating}
    </span>
  )
}
