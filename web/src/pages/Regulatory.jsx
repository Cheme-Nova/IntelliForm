import { useState } from 'react'
import { api } from '../api/client'

export default function Regulatory() {
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
        input_text: 'generate regulatory report',
        vertical,
        batch_size: 1000,
        opt_mode: 'single',
        constraints: { blend: parsedBlend }
      })
      setResult({ reg: res.data?.reg, vreg: res.data?.vreg })
    } catch (err) {
      setError(err.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        📋 Regulatory
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        REACH, SDS flags, restricted substances, and vertical-specific compliance
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
            {loading ? 'Checking...' : 'Run Report →'}
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
          {result.reg?.flags?.length > 0 && (
            <div style={{
              padding: '1rem', background: '#450a0a', border: '1px solid #7f1d1d',
              borderRadius: '8px'
            }}>
              <div style={{ color: '#ef4444', fontWeight: 600, marginBottom: '0.5rem', fontSize: '0.85rem' }}>
                ⚠️ Regulatory Flags
              </div>
              {result.reg.flags.map((f, i) => (
                <div key={i} style={{ color: '#fca5a5', fontSize: '0.85rem', marginBottom: '4px' }}>· {f}</div>
              ))}
            </div>
          )}
          {result.vreg && (
            <div style={{
              padding: '1rem', background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px'
            }}>
              <div style={{ color: '#0D9488', fontWeight: 600, marginBottom: '0.5rem', fontSize: '0.85rem' }}>
                Vertical Compliance — {vertical.replace('_',' ')}
              </div>
              <pre style={{ color: '#94a3b8', fontSize: '0.8rem', overflow: 'auto', margin: 0 }}>
                {JSON.stringify(result.vreg, null, 2)}
              </pre>
            </div>
          )}
        </div>
      ) : (
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#334155', textAlign: 'center', fontSize: '0.85rem'
        }}>
          Enter a blend above and click Run Report
        </div>
      )}
    </div>
  )
}
