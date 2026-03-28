import { useState } from 'react'
import { api } from '../api/client'

const VERTICALS = ['personal_care', 'home_care', 'industrial', 'pharma', 'food', 'agriculture']

export default function Formulate() {
  const [inputText, setInputText] = useState('')
  const [vertical, setVertical] = useState('personal_care')
  const [batchSize, setBatchSize] = useState(1000)
  const [optMode, setOptMode] = useState('auto')
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)

  async function handleSubmit() {
    setLoading(true)
    setError(null)
    setResult(null)
    try {
      const res = await api.formulate({
        input_text: inputText,
        vertical,
        batch_size: batchSize,
        opt_mode: optMode,
        constraints: {}
      })
      setResult(res.data)
    } catch (err) {
      setError(err.response?.data?.detail || err.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        ⚗️ Formulate
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Describe your formulation goal in plain language
      </p>

      <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
        <textarea
          value={inputText}
          onChange={e => setInputText(e.target.value)}
          placeholder="e.g. Create a sulfate-free shampoo with >80% bio-based content, mild pH, stable foam..."
          rows={4}
          style={{
            background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px',
            color: '#fff', padding: '1rem', fontSize: '0.9rem', resize: 'vertical',
            outline: 'none', fontFamily: 'inherit'
          }}
        />

        <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap' }}>
          <div style={{ flex: 1, minWidth: '160px' }}>
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>
              VERTICAL
            </label>
            <select
              value={vertical}
              onChange={e => setVertical(e.target.value)}
              style={{
                width: '100%', background: '#0D1F3C', border: '1px solid #1e3a5f',
                borderRadius: '6px', color: '#fff', padding: '0.5rem', fontSize: '0.85rem'
              }}
            >
              {VERTICALS.map(v => (
                <option key={v} value={v}>{v.replace('_', ' ')}</option>
              ))}
            </select>
          </div>

          <div style={{ flex: 1, minWidth: '120px' }}>
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>
              BATCH SIZE (kg)
            </label>
            <input
              type="number"
              value={batchSize}
              onChange={e => setBatchSize(Number(e.target.value))}
              style={{
                width: '100%', background: '#0D1F3C', border: '1px solid #1e3a5f',
                borderRadius: '6px', color: '#fff', padding: '0.5rem', fontSize: '0.85rem',
                boxSizing: 'border-box'
              }}
            />
          </div>

          <div style={{ flex: 1, minWidth: '120px' }}>
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>
              OPT MODE
            </label>
            <select
              value={optMode}
              onChange={e => setOptMode(e.target.value)}
              style={{
                width: '100%', background: '#0D1F3C', border: '1px solid #1e3a5f',
                borderRadius: '6px', color: '#fff', padding: '0.5rem', fontSize: '0.85rem'
              }}
            >
              <option value="auto">Auto</option>
              <option value="pareto">Pareto</option>
              <option value="bayesian">Bayesian</option>
              <option value="single">Single</option>
            </select>
          </div>
        </div>

        <button
          onClick={handleSubmit}
          disabled={loading || !inputText.trim()}
          style={{
            background: loading ? '#334155' : '#0D9488',
            color: '#fff', border: 'none', borderRadius: '8px',
            padding: '0.75rem 2rem', fontSize: '0.95rem', fontWeight: 600,
            cursor: loading ? 'not-allowed' : 'pointer', alignSelf: 'flex-start',
            transition: 'background 0.15s'
          }}
        >
          {loading ? 'Optimizing...' : 'Run IntelliForm →'}
        </button>
      </div>

      {error && (
        <div style={{
          marginTop: '1.5rem', padding: '1rem', background: '#450a0a',
          border: '1px solid #7f1d1d', borderRadius: '8px', color: '#fca5a5', fontSize: '0.85rem'
        }}>
          {error}
        </div>
      )}

      {result && (
        <div style={{
          marginTop: '1.5rem', padding: '1.5rem', background: '#0D1F3C',
          border: '1px solid #1e3a5f', borderRadius: '8px'
        }}>
          <div style={{ color: '#0D9488', fontWeight: 700, marginBottom: '1rem' }}>
            ✅ Formulation Result
          </div>
          <pre style={{ color: '#94a3b8', fontSize: '0.8rem', overflow: 'auto', margin: 0 }}>
            {JSON.stringify(result, null, 2)}
          </pre>
        </div>
      )}
    </div>
  )
}
