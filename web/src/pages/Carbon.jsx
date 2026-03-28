import { useState } from 'react'
import { api } from '../api/client'

export default function Carbon() {
  const [blend, setBlend] = useState('{"SLS": 12, "Cocamide DEA": 3, "Water": 85}')
  const [batchSize, setBatchSize] = useState(1000)
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
        input_text: 'calculate carbon credits',
        vertical,
        batch_size: batchSize,
        opt_mode: 'single',
        constraints: { blend: parsedBlend }
      })
      setResult(res.data?.carbon)
    } catch (err) {
      setError(err.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        ♻️ Carbon Credits
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Estimate carbon savings and verified credit value per batch
      </p>

      <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem', marginBottom: '1.5rem' }}>
        <div>
          <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>BLEND (JSON)</label>
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
        <div style={{ display: 'flex', gap: '1rem', alignItems: 'flex-end', flexWrap: 'wrap' }}>
          <div>
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>VERTICAL</label>
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
          <div>
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>BATCH SIZE (kg)</label>
            <input
              type="number"
              value={batchSize}
              onChange={e => setBatchSize(Number(e.target.value))}
              style={{
                width: '100px', background: '#0D1F3C', border: '1px solid #1e3a5f',
                borderRadius: '6px', color: '#fff', padding: '0.5rem', fontSize: '0.85rem'
              }}
            />
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
            {loading ? 'Calculating...' : 'Calculate →'}
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
          {[
            { label: 'CO₂e SAVED', value: result.co2_saved_kg?.toFixed(2), unit: 'kg per batch', color: '#0D9488' },
            { label: 'CREDIT VALUE', value: result.credit_value ? `$${result.credit_value.toFixed(2)}` : '—', unit: 'USD', color: '#D97706' },
            { label: 'ANNUAL SAVINGS', value: result.annual_co2_saved_kg?.toFixed(0), unit: 'kg CO₂e/yr' },
            { label: 'BASELINE', value: result.baseline_kg?.toFixed(2), unit: 'kg CO₂e baseline' },
          ].map(({ label, value, unit, color }) => (
            <div key={label} style={{
              background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px',
              padding: '1.25rem', flex: 1, minWidth: '140px'
            }}>
              <div style={{ color: '#64748b', fontSize: '0.7rem', marginBottom: '4px' }}>{label}</div>
              <div style={{ color: color || '#fff', fontSize: '1.4rem', fontWeight: 700 }}>{value ?? '—'}</div>
              <div style={{ color: '#475569', fontSize: '0.7rem' }}>{unit}</div>
            </div>
          ))}
        </div>
      ) : (
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#334155', textAlign: 'center', fontSize: '0.85rem'
        }}>
          Enter a blend above and click Calculate
        </div>
      )}
    </div>
  )
}
