import { useState } from 'react'
import { api } from '../api/client'

const RISK_COLOR = { low: '#0D9488', medium: '#D97706', high: '#ef4444' }

export default function Stability() {
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
        input_text: 'predict stability',
        vertical,
        batch_size: 1000,
        opt_mode: 'single',
        constraints: { blend: parsedBlend }
      })
      setResult(res.data?.stability)
    } catch (err) {
      setError(err.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        🧪 Stability
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Predict phase separation, viscosity drift, pH shift and shelf life
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
        <div style={{ display: 'flex', gap: '1rem', alignItems: 'flex-end' }}>
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
          <button
            onClick={handleRun}
            disabled={loading}
            style={{
              background: loading ? '#334155' : '#0D9488', color: '#fff',
              border: 'none', borderRadius: '8px', padding: '0.6rem 1.5rem',
              fontSize: '0.9rem', fontWeight: 600, cursor: loading ? 'not-allowed' : 'pointer'
            }}
          >
            {loading ? 'Predicting...' : 'Predict Stability →'}
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
        <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
          <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap' }}>
            {[
              { label: 'OVERALL RISK', value: result.overall_risk, isRisk: true },
              { label: 'SHELF LIFE', value: result.shelf_life_months ? `${result.shelf_life_months} mo` : '—' },
              { label: 'VISCOSITY DRIFT', value: result.viscosity_risk, isRisk: true },
              { label: 'PHASE SEP RISK', value: result.phase_sep_risk, isRisk: true },
            ].map(({ label, value, isRisk }) => (
              <div key={label} style={{
                background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px',
                padding: '1.25rem', flex: 1, minWidth: '140px'
              }}>
                <div style={{ color: '#64748b', fontSize: '0.7rem', marginBottom: '4px' }}>{label}</div>
                <div style={{
                  color: isRisk ? (RISK_COLOR[value?.toLowerCase()] || '#fff') : '#fff',
                  fontSize: '1.2rem', fontWeight: 700, textTransform: 'capitalize'
                }}>
                  {value ?? '—'}
                </div>
              </div>
            ))}
          </div>
          {result.recommendations?.length > 0 && (
            <div style={{
              padding: '1rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
              borderRadius: '8px'
            }}>
              <div style={{ color: '#D97706', fontWeight: 600, marginBottom: '0.5rem', fontSize: '0.85rem' }}>
                Recommendations
              </div>
              {result.recommendations.map((r, i) => (
                <div key={i} style={{ color: '#94a3b8', fontSize: '0.85rem', marginBottom: '4px' }}>
                  · {r}
                </div>
              ))}
            </div>
          )}
        </div>
      ) : (
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#334155', textAlign: 'center', fontSize: '0.85rem'
        }}>
          Enter a blend above and click Predict Stability
        </div>
      )}
    </div>
  )
}
