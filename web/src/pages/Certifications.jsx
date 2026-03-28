import { useState } from 'react'
import { api } from '../api/client'

const CERTS = ['ECOCERT', 'COSMOS', 'USDA_BIO', 'EU_ECOLABEL', 'CRADLE_TO_CRADLE', 'NSF']

const Badge = ({ name, status }) => {
  const color = status === 'pass' ? '#0D9488' : status === 'fail' ? '#ef4444' : '#D97706'
  const bg = status === 'pass' ? '#0D948822' : status === 'fail' ? '#ef444422' : '#D9770622'
  const label = status === 'pass' ? '✅ Pass' : status === 'fail' ? '❌ Fail' : '⚠️ Review'
  return (
    <div style={{
      background: bg, border: `1px solid ${color}`, borderRadius: '8px',
      padding: '1rem', flex: 1, minWidth: '140px'
    }}>
      <div style={{ color: '#94a3b8', fontSize: '0.7rem', marginBottom: '4px' }}>{name}</div>
      <div style={{ color, fontWeight: 700, fontSize: '1rem' }}>{label}</div>
    </div>
  )
}

export default function Certifications() {
  const [blend, setBlend] = useState('{"SLS": 12, "Cocamide DEA": 3, "Water": 85}')
  const [vertical, setVertical] = useState('personal_care')
  const [bioPct, setBioPct] = useState(80)
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)

  async function handleRun() {
    setLoading(true)
    setError(null)
    try {
      const parsedBlend = JSON.parse(blend)
      const res = await api.formulate({
        input_text: 'check certifications',
        vertical,
        batch_size: 1000,
        opt_mode: 'single',
        constraints: { blend: parsedBlend, bio_pct: bioPct }
      })
      setResult(res.data?.cert)
    } catch (err) {
      setError(err.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        ✅ Certifications
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Check blend eligibility against ECOCERT, COSMOS, USDA, EU Ecolabel and more
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

        <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap', alignItems: 'flex-end' }}>
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
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>BIO-BASED %</label>
            <input
              type="number"
              value={bioPct}
              onChange={e => setBioPct(Number(e.target.value))}
              min={0} max={100}
              style={{
                width: '80px', background: '#0D1F3C', border: '1px solid #1e3a5f',
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
            {loading ? 'Checking...' : 'Check Certs →'}
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
            {CERTS.map(cert => (
              <Badge
                key={cert}
                name={cert.replace(/_/g, ' ')}
                status={result[cert.toLowerCase()]?.status || 'review'}
              />
            ))}
          </div>
          {result.notes && (
            <div style={{
              padding: '1rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
              borderRadius: '8px', color: '#94a3b8', fontSize: '0.85rem'
            }}>
              {result.notes}
            </div>
          )}
        </div>
      ) : (
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#334155', textAlign: 'center', fontSize: '0.85rem'
        }}>
          Enter a blend above and click Check Certs
        </div>
      )}
    </div>
  )
}
