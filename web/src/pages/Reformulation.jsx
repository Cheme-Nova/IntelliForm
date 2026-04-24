import { useState } from 'react'
import { api } from '../api/client'

const FAILURE_TYPES = ['viscosity', 'stability', 'pH', 'color', 'odor', 'certification', 'eco_score']

export default function Reformulation() {
  const [blend, setBlend] = useState('{"SLS": 12, "Cocamide DEA": 3, "Water": 85}')
  const [failureType, setFailureType] = useState('viscosity')
  const [vertical, setVertical] = useState('personal_care')
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)

  async function handleRun() {
    setLoading(true)
    setError(null)
    try {
      const parsedBlend = JSON.parse(blend)
      const res = await api.reformulate({
        blend: parsedBlend,
        failure_type: failureType,
        vertical,
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
        🔧 Reformulation Lab
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Diagnose failures and get ranked fix recommendations
      </p>

      <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem', marginBottom: '1.5rem' }}>
        <div>
          <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>CURRENT BLEND (JSON)</label>
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
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>FAILURE TYPE</label>
            <select
              value={failureType}
              onChange={e => setFailureType(e.target.value)}
              style={{
                background: '#0D1F3C', border: '1px solid #1e3a5f',
                borderRadius: '6px', color: '#fff', padding: '0.5rem', fontSize: '0.85rem'
              }}
            >
              {FAILURE_TYPES.map(f => (
                <option key={f} value={f}>{f.replace('_',' ')}</option>
              ))}
            </select>
          </div>
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
              {['personal_care','industrial','agricultural','pharmaceutical','food','fabric_laundry','paint_coatings'].map(v => (
                <option key={v} value={v}>{v.replace('_',' ')}</option>
              ))}
            </select>
          </div>
          <button
            onClick={handleRun}
            disabled={loading}
            style={{
              background: loading ? '#334155' : '#D97706', color: '#fff',
              border: 'none', borderRadius: '8px', padding: '0.6rem 1.5rem',
              fontSize: '0.9rem', fontWeight: 600, cursor: loading ? 'not-allowed' : 'pointer'
            }}
          >
            {loading ? 'Diagnosing...' : 'Diagnose & Fix →'}
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
          {result.root_cause && (
            <div style={{
              padding: '1rem', background: '#1c1107', border: '1px solid #D97706',
              borderRadius: '8px'
            }}>
              <div style={{ color: '#D97706', fontWeight: 600, marginBottom: '4px', fontSize: '0.85rem' }}>Root Cause</div>
              <div style={{ color: '#fbbf24', fontSize: '0.85rem' }}>{result.root_cause.root_cause}</div>
            </div>
          )}
          {(result.suggestions?.length ? result.suggestions : [result.best_suggestion]).filter(Boolean).map((rec, i) => (
            <div key={i} style={{
              background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem',
              display: 'flex', gap: '1rem', alignItems: 'flex-start'
            }}>
              <div style={{
                background: '#0D9488', color: '#fff', borderRadius: '50%',
                width: '24px', height: '24px', display: 'flex', alignItems: 'center',
                justifyContent: 'center', fontSize: '0.75rem', fontWeight: 700, flexShrink: 0
              }}>
                {i + 1}
              </div>
              <div>
                <div style={{ color: '#fff', fontWeight: 600, fontSize: '0.85rem' }}>
                  {rec.action_type} {rec.ingredient} → {rec.suggested_pct}%
                </div>
                {rec.rationale && <div style={{ color: '#64748b', fontSize: '0.8rem', marginTop: '2px' }}>{rec.rationale}</div>}
              </div>
            </div>
          ))}
        </div>
      ) : (
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#334155', textAlign: 'center', fontSize: '0.85rem'
        }}>
          Enter a failing blend and click Diagnose & Fix
        </div>
      )}
    </div>
  )
}
