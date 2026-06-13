import { useState, useMemo } from 'react'
import { api } from '../api/client'
import { loadLastRun } from '../lib/session'

const card = { background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }
const input = { width: '100%', background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '6px', color: '#fff', padding: '0.55rem 0.8rem', fontSize: '0.85rem', boxSizing: 'border-box', fontFamily: 'inherit' }
const label = { color: '#64748b', fontSize: '0.72rem', display: 'block', marginBottom: '4px' }
const btn = (disabled) => ({ background: disabled ? '#334155' : '#0D9488', color: '#fff', border: 'none', borderRadius: '7px', padding: '0.6rem 1.4rem', fontSize: '0.88rem', fontWeight: 600, cursor: disabled ? 'not-allowed' : 'pointer' })

const VERTICALS = [
  'personal_care', 'industrial', 'agricultural', 'pharmaceutical', 'food', 'fabric_laundry', 'paint_coatings',
]

export default function SDS() {
  const session = useMemo(() => loadLastRun(), [])
  const lastBlend = session?.request?.blend || null

  const [productName, setProductName] = useState('IntelliForm Formulation')
  const [vertical, setVertical] = useState('personal_care')
  const [version, setVersion] = useState('1.0')
  const [format, setFormat] = useState('json')
  const [blendText, setBlendText] = useState(
    lastBlend ? Object.entries(lastBlend).map(([k, v]) => `${k}: ${v}`).join('\n') : ''
  )
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)

  function parseBlendText(text) {
    const blend = {}
    for (const line of text.split('\n')) {
      const trimmed = line.trim()
      if (!trimmed) continue
      const sep = trimmed.includes(':') ? ':' : trimmed.includes(',') ? ',' : null
      if (!sep) continue
      const [key, val] = trimmed.split(sep, 2)
      const pct = parseFloat(val)
      if (key && !isNaN(pct)) blend[key.trim()] = pct
    }
    return blend
  }

  async function generate() {
    const blend = parseBlendText(blendText)
    if (!Object.keys(blend).length) { setError('Enter at least one ingredient (e.g. "Coco-Glucoside: 35").'); return }
    setLoading(true); setError(null); setResult(null)
    try {
      const res = await api.exportSDS({ blend, product_name: productName, vertical, version, format: 'json' })
      setResult(res.data)
    } catch (err) {
      setError(err.response?.data?.detail || err.message)
    } finally { setLoading(false) }
  }

  const sds = result?.sds

  return (
    <div style={{ maxWidth: '860px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>📄 Safety Data Sheet</h1>
      <p style={{ color: '#64748b', marginBottom: '1.75rem', fontSize: '0.9rem' }}>
        Generate a GHS EU CLP Annex II 16-section SDS for any formulation blend.
      </p>

      {/* Config */}
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem', marginBottom: '1rem' }}>
        <div><span style={label}>Product name</span><input style={input} value={productName} onChange={e => setProductName(e.target.value)} /></div>
        <div>
          <span style={label}>Vertical</span>
          <select style={input} value={vertical} onChange={e => setVertical(e.target.value)}>
            {VERTICALS.map(v => <option key={v} value={v}>{v.replace(/_/g, ' ')}</option>)}
          </select>
        </div>
      </div>
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 3fr', gap: '1rem', marginBottom: '1rem' }}>
        <div><span style={label}>Version</span><input style={input} value={version} onChange={e => setVersion(e.target.value)} placeholder="1.0" /></div>
        <div>
          <span style={label}>Blend (ingredient: %, one per line)</span>
          <textarea
            style={{ ...input, minHeight: '120px', resize: 'vertical', fontFamily: 'monospace', fontSize: '0.8rem' }}
            value={blendText}
            onChange={e => setBlendText(e.target.value)}
            placeholder="Coco-Glucoside: 35&#10;Glycerin: 10&#10;Water: 55"
          />
          {lastBlend && <div style={{ color: '#64748b', fontSize: '0.7rem', marginTop: '4px' }}>Pre-filled from your last formulation run.</div>}
        </div>
      </div>

      <div style={{ display: 'flex', gap: '0.75rem', alignItems: 'center', marginBottom: '1.5rem', flexWrap: 'wrap' }}>
        <button onClick={generate} disabled={loading} style={btn(loading)}>
          {loading ? 'Generating SDS…' : 'Generate SDS →'}
        </button>
        {result?.download_url && (
          <a href={result.download_url} target="_blank" rel="noreferrer"
            style={{ ...btn(false), textDecoration: 'none', background: '#1e3a5f' }}>
            ⬇ Download PDF
          </a>
        )}
      </div>

      {error && <div style={{ ...card, borderColor: '#7f1d1d', color: '#fca5a5', marginBottom: '1rem' }}>{error}</div>}

      {sds && (
        <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem' }}>
          {/* Header strip */}
          <div style={{ ...card, background: '#0a2a40', borderColor: '#0D9488', display: 'flex', gap: '2rem', flexWrap: 'wrap' }}>
            <div>
              <div style={{ color: '#64748b', fontSize: '0.7rem' }}>PRODUCT</div>
              <div style={{ color: '#fff', fontWeight: 700, fontSize: '1rem' }}>{sds.product_name}</div>
            </div>
            <div>
              <div style={{ color: '#64748b', fontSize: '0.7rem' }}>SIGNAL WORD</div>
              <div style={{ color: sds.signal_word === 'Danger' ? '#ef4444' : '#D97706', fontWeight: 700, fontSize: '1rem' }}>{sds.signal_word || '—'}</div>
            </div>
            <div>
              <div style={{ color: '#64748b', fontSize: '0.7rem' }}>REVISION DATE</div>
              <div style={{ color: '#cbd5e1', fontWeight: 600 }}>{sds.revision_date || '—'}</div>
            </div>
            <div>
              <div style={{ color: '#64748b', fontSize: '0.7rem' }}>VERSION</div>
              <div style={{ color: '#cbd5e1', fontWeight: 600 }}>{sds.version || '—'}</div>
            </div>
          </div>

          {/* Hazards */}
          {sds.h_statements?.length > 0 && (
            <div style={card}>
              <div style={{ color: '#D97706', fontWeight: 700, fontSize: '0.78rem', marginBottom: '0.6rem', textTransform: 'uppercase', letterSpacing: '0.07em' }}>Hazard Statements</div>
              {sds.h_statements.map((h, i) => (
                <div key={i} style={{ color: '#cbd5e1', fontSize: '0.85rem', marginBottom: '0.3rem' }}>
                  <span style={{ color: '#D97706', fontFamily: 'monospace', marginRight: '0.5rem' }}>{h.code}</span>{h.text}
                </div>
              ))}
            </div>
          )}

          {/* Pictograms */}
          {sds.pictograms?.length > 0 && (
            <div style={card}>
              <div style={{ color: '#D97706', fontWeight: 700, fontSize: '0.78rem', marginBottom: '0.6rem', textTransform: 'uppercase', letterSpacing: '0.07em' }}>GHS Pictograms</div>
              <div style={{ display: 'flex', gap: '0.6rem', flexWrap: 'wrap' }}>
                {sds.pictograms.map((p, i) => (
                  <span key={i} style={{ padding: '0.3rem 0.7rem', background: '#1e3a5f', borderRadius: '6px', color: '#94a3b8', fontSize: '0.8rem', fontFamily: 'monospace' }}>{p}</span>
                ))}
              </div>
            </div>
          )}

          {/* P-statements */}
          {sds.p_statements?.length > 0 && (
            <div style={card}>
              <div style={{ color: '#D97706', fontWeight: 700, fontSize: '0.78rem', marginBottom: '0.6rem', textTransform: 'uppercase', letterSpacing: '0.07em' }}>Precautionary Statements</div>
              {sds.p_statements.map((p, i) => (
                <div key={i} style={{ color: '#cbd5e1', fontSize: '0.85rem', marginBottom: '0.3rem' }}>
                  <span style={{ color: '#94a3b8', fontFamily: 'monospace', marginRight: '0.5rem' }}>{p.code}</span>{p.text}
                </div>
              ))}
            </div>
          )}

          {/* Ingredients table */}
          {sds.ingredients?.length > 0 && (
            <div style={card}>
              <div style={{ color: '#D97706', fontWeight: 700, fontSize: '0.78rem', marginBottom: '0.75rem', textTransform: 'uppercase', letterSpacing: '0.07em' }}>Section 3 — Composition</div>
              <div style={{ overflowX: 'auto' }}>
                <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.82rem' }}>
                  <thead>
                    <tr>{['Ingredient', 'CAS', 'EC', '% w/w', 'REACH', 'Hazardous'].map(h => (
                      <th key={h} style={{ color: '#64748b', textAlign: 'left', padding: '0.4rem 0.6rem', borderBottom: '1px solid #1e3a5f', whiteSpace: 'nowrap' }}>{h}</th>
                    ))}</tr>
                  </thead>
                  <tbody>
                    {sds.ingredients.map((ing, i) => (
                      <tr key={i} style={{ borderBottom: '1px solid #0f2a4a' }}>
                        <td style={{ color: '#fff', padding: '0.4rem 0.6rem', fontWeight: 500 }}>{ing.name}</td>
                        <td style={{ color: '#94a3b8', padding: '0.4rem 0.6rem', fontFamily: 'monospace', fontSize: '0.78rem' }}>{ing.cas || '—'}</td>
                        <td style={{ color: '#94a3b8', padding: '0.4rem 0.6rem', fontFamily: 'monospace', fontSize: '0.78rem' }}>{ing.ec || '—'}</td>
                        <td style={{ color: '#0D9488', padding: '0.4rem 0.6rem', fontFamily: 'monospace' }}>{ing.pct?.toFixed(1)}%</td>
                        <td style={{ color: '#64748b', padding: '0.4rem 0.6rem', fontSize: '0.78rem' }}>{ing.reach_status || '—'}</td>
                        <td style={{ padding: '0.4rem 0.6rem' }}>
                          <span style={{ color: ing.is_hazardous ? '#fca5a5' : '#0D9488', fontSize: '0.78rem' }}>{ing.is_hazardous ? 'Yes' : 'No'}</span>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          )}

          {/* Full SDS JSON fallback */}
          <details>
            <summary style={{ color: '#64748b', fontSize: '0.8rem', cursor: 'pointer', marginTop: '0.5rem' }}>View full SDS JSON</summary>
            <pre style={{ ...card, color: '#94a3b8', fontSize: '0.72rem', marginTop: '0.5rem', overflowX: 'auto', maxHeight: '400px' }}>
              {JSON.stringify(sds, null, 2)}
            </pre>
          </details>
        </div>
      )}
    </div>
  )
}
