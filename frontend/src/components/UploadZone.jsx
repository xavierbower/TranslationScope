import React, { useCallback } from 'react'

export default function UploadZone({ file, onFile, label }) {
  const handleDrop = useCallback((e) => {
    e.preventDefault()
    const f = e.dataTransfer?.files?.[0]
    const exts = ['.dna', '.gbk', '.gb', '.genbank']
    if (f && exts.some(e => f.name.toLowerCase().endsWith(e))) onFile(f)
  }, [onFile])

  const handleChange = useCallback((e) => {
    const f = e.target.files?.[0]
    if (f) onFile(f)
  }, [onFile])

  return (
    <div
      onDrop={handleDrop}
      onDragOver={e => e.preventDefault()}
      className="border-2 border-dashed border-gray-300 rounded-xl p-6 text-center
                 hover:border-blue-400 hover:bg-blue-50/30 transition cursor-pointer"
      onClick={() => document.getElementById(`upload-${label}`)?.click()}
    >
      <input
        id={`upload-${label}`}
        type="file"
        accept=".dna,.gbk,.gb,.genbank"
        className="hidden"
        onChange={handleChange}
      />
      {file ? (
        <div>
          <p className="text-sm text-gray-500">Selected:</p>
          <p className="font-medium text-gray-800">{file.name}</p>
          <p className="text-xs text-gray-400 mt-1">{(file.size / 1024).toFixed(1)} KB</p>
        </div>
      ) : (
        <div>
          <svg className="mx-auto h-8 w-8 text-gray-400 mb-2" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
              d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M15 13l-3-3m0 0l-3 3m3-3v12" />
          </svg>
          <p className="text-sm text-gray-500">{label}</p>
          <p className="text-xs text-gray-400 mt-1">Drag & drop or click to browse (.dna or .gbk)</p>
        </div>
      )}
    </div>
  )
}
