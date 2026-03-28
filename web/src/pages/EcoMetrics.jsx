import { useState } from 'react'
import { api } from '../api/client'

const CARD = ({ label, value, unit, color }) => (
  <div style={{
    background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px',
    padding: '1.25rem', flex: 1, minWidth: '140px'
  }}>
    <div style={{ color: '#64748b', fontSize: '0.7rem', marginBottom: '4px' }}>{label}</div>
    <div style={{ color: color || '#fff', fontSize: '1.6rem', fontWeight: 700 }}>{value ?? '—'}</div>
    {unit && <div style={{ color: '#475569', fontSize: '0.7rem' }}>{unit}</div>}
  </div>
)

export default function EcoMetrics() {
  const [blend, setBlend] = useState('{"SLS": 12, "Cocamide DEA": 3, "Water": 85}')
  const [vertical, setVertical] = useState('personal_care')
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)

  async function handleRun() {
    setLoading(true)
    setError(null)
    try {
      const parsedBlend = JSON.parse(blend)
      const res = await api.formulate({
        input_text: 'compute eco metrics',
        vertical,
        batch_size: 1000,
        opt_mode: 'single',
        constraints: { blend: parsedBlend }
      })
      setResult(res.data?.eco)
    } catch (err) {
      setError(err.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        🌱 EcoMetrics
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Sustainability scoring against ISO 14040/14044 benchmarks · Spearman ρ = 0.83
      </p>

      <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem', marginBottom: '1.5rem' }}>
        <div>
          <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>
            BLEND (JSON)
          </label>
          <textarea
            value={blend}
            onChange={e => setBlend(e.target.value)}
            rows={4}
            style={{
              width: '100%', background: '#0D1F3C', border: '1px solid #1e3a5f',
              borderRadius: '8px', color: '#fff', padding: '1rem', fontSize: '0.85rem',
              fontFamily: 'monospace', boxSizing: 'border-box', resize: 'vertical'
            }}
          />
        </div>

        <div style={{ display: 'flex', gap: '1rem', alignItems: 'flex-end' }}>
          <div>
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>
              VERTICAL
            </label>
            <select
              value={vertical}
              onChange={e => setVertical(e.target.value)}
              style={{
                background: '#0D1F3C', border: '1px solid #1e3a5f',
                borderRadius: '6px', color: '#fff', padding: '0.5rem', fontSize: '0.85rem'
              }}
            >
              {['personal_care','home_care','industrial','pharma','food','agriculture'].map(v => (
                <option key={v} value={v}>{v.replace('_',' ')}</option>
              ))}
            </select>
          </div>
          <button
            onClick={handleRun}
            disabled={loading}
            style={{
              background: loading ? '#334155' : '#0D9488', color: '#fff',
              border: 'none', borderRadius: '8px', padding: '0.6rem 1.5rem',
              fontSize: '0.9rem', fontWeight: 600, cursor: loading ? 'not-allowed' : 'pointer'
            }}
          >
            {loading ? 'Scoring...' : 'Score Blend →'}
          </button>
        </div>
      </div>

      {error && (
        <div style={{
          padding: '1rem', background: '#450a0a', border: '1px solid #7f1d1d',
          borderRadius: '8px', color: '#fca5a5', fontSize: '0.85rem', marginBottom: '1rem'
        }}>
          {error}
        </div>
      )}

      {result ? (
        <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap' }}>
          <CARD label="ECO SCORE" value={result.score?.toFixed(1)} unit="/ 100" color="#0D9488" />
          <CARD label="BIO-BASED %" value={result.bio_pct?.toFixed(1)} unit="%" color="#D97706" />
          <CARD label="WASTE SCORE" value={result.waste_score?.toFixed(1)} unit="/ 100" />
          <CARD label="CARBON kg CO₂e" value={result.carbon_kg?.toFixed(2)} unit="per batch" />
        </div>
      ) : (
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#334155', textAlign: 'center', fontSize: '0.85rem'
        }}>
          Enter a blend above and click Score Blend
        </div>
      )}
    </div>
  )
}
