import { useState, useMemo, useEffect } from 'react'
import { api } from '../api/client'
import { loadLastRun } from '../lib/session'
import { BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, Cell } from 'recharts'

const card = { background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }
const label = { color: '#64748b', fontSize: '0.72rem', display: 'block', marginBottom: '4px' }
const btn = (disabled) => ({ background: disabled ? '#334155' : '#0D9488', color: '#fff', border: 'none', borderRadius: '7px', padding: '0.55rem 1.2rem', fontSize: '0.85rem', fontWeight: 600, cursor: disabled ? 'not-allowed' : 'pointer' })

function MetricBox({ label: lbl, value, sub }) {
  return (
    <div style={{ ...card, flex: 1, minWidth: '120px' }}>
      <div style={{ color: '#64748b', fontSize: '0.68rem', marginBottom: '4px', textTransform: 'uppercase', letterSpacing: '0.06em' }}>{lbl}</div>
      <div style={{ color: '#fff', fontSize: '1.25rem', fontWeight: 700 }}>{value ?? '—'}</div>
      {sub && <div style={{ color: '#64748b', fontSize: '0.7rem', marginTop: '2px' }}>{sub}</div>}
    </div>
  )
}

export default function History() {
  const session = useMemo(() => loadLastRun(), [])
  const result = session?.response
  const blend = session?.request?.blend || (result?.blend ?? null)

  const [pricing, setPricing] = useState(null)
  const [risk, setRisk] = useState(null)
  const [loadingMkt, setLoadingMkt] = useState(false)
  const [mktError, setMktError] = useState(null)

  async function loadMarket() {
    if (!blend) return
    setLoadingMkt(true); setMktError(null)
    try {
      const [pRes, rRes] = await Promise.all([
        api.supplyChainPricing(blend),
        api.supplyChainRisk(blend),
      ])
      setPricing(pRes.data.pricing || [])
      setRisk(rRes.data)
    } catch (err) {
      setMktError(err.response?.data?.detail || err.message)
    } finally { setLoadingMkt(false) }
  }

  if (!session) {
    return (
      <div style={{ maxWidth: '800px' }}>
        <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>📊 History & Market Intel</h1>
        <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>Session history and supply chain analysis for your last formulation run.</p>
        <div style={{ ...card, color: '#94a3b8', textAlign: 'center', fontSize: '0.9rem', lineHeight: 1.7 }}>
          Run a formulation first to see session history and market pricing.
        </div>
      </div>
    )
  }

  const savedAt = session.savedAt ? new Date(session.savedAt).toLocaleString() : '—'
  const req = session.request || {}
  const blendEntries = blend ? Object.entries(blend) : []

  const riskColor = { low: '#0D9488', medium: '#D97706', high: '#ef4444' }

  return (
    <div style={{ maxWidth: '860px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>📊 History & Market Intel</h1>
      <p style={{ color: '#64748b', marginBottom: '1.75rem', fontSize: '0.9rem' }}>Last formulation run · {savedAt}</p>

      {/* Summary metrics */}
      {result && (
        <div style={{ display: 'flex', gap: '0.75rem', flexWrap: 'wrap', marginBottom: '1.25rem' }}>
          <MetricBox label="Cost/kg" value={result.cost_per_kg ? `$${result.cost_per_kg.toFixed(2)}` : '—'} />
          <MetricBox label="Bio-based" value={result.bio_pct ? `${result.bio_pct.toFixed(1)}%` : '—'} />
          <MetricBox label="Perf score" value={result.perf_score ? `${result.perf_score}/100` : '—'} />
          <MetricBox label="EcoScore" value={result.eco_score ? `${result.eco_score.toFixed(1)}/100` : '—'} />
          <MetricBox label="Vertical" value={req.vertical?.replace(/_/g, ' ') || '—'} />
        </div>
      )}

      {/* Blend composition */}
      {blendEntries.length > 0 && (
        <div style={{ ...card, marginBottom: '1.25rem' }}>
          <div style={{ color: '#D97706', fontWeight: 700, fontSize: '0.75rem', letterSpacing: '0.07em', textTransform: 'uppercase', marginBottom: '0.75rem' }}>Blend Composition</div>
          {blendEntries.map(([ing, pct]) => (
            <div key={ing} style={{ display: 'flex', alignItems: 'center', gap: '0.75rem', marginBottom: '0.55rem' }}>
              <div style={{ color: '#fff', fontSize: '0.85rem', minWidth: '180px', fontWeight: 500 }}>{ing}</div>
              <div style={{ flex: 1, height: '6px', background: '#1e3a5f', borderRadius: '3px', overflow: 'hidden' }}>
                <div style={{ height: '100%', width: `${Math.min(pct, 100)}%`, background: '#0D9488', borderRadius: '3px' }} />
              </div>
              <div style={{ color: '#0D9488', fontFamily: 'monospace', fontSize: '0.82rem', minWidth: '44px', textAlign: 'right' }}>{typeof pct === 'number' ? pct.toFixed(1) : pct}%</div>
            </div>
          ))}
        </div>
      )}

      {/* Market Intel section */}
      <div style={{ ...card, marginBottom: '1.25rem' }}>
        <div style={{ color: '#D97706', fontWeight: 700, fontSize: '0.75rem', letterSpacing: '0.07em', textTransform: 'uppercase', marginBottom: '0.75rem' }}>Supply Chain & Market Pricing</div>
        {!pricing && !risk && (
          <div style={{ display: 'flex', gap: '1rem', alignItems: 'center' }}>
            <p style={{ color: '#64748b', fontSize: '0.85rem', margin: 0 }}>
              Load live supplier pricing and risk scores for this blend.
            </p>
            <button onClick={loadMarket} disabled={loadingMkt || !blend} style={btn(loadingMkt || !blend)}>
              {loadingMkt ? 'Loading…' : 'Load market data →'}
            </button>
          </div>
        )}
        {mktError && <div style={{ color: '#fca5a5', fontSize: '0.83rem', marginTop: '0.5rem' }}>{mktError}</div>}

        {risk && (
          <div style={{ marginBottom: '1rem' }}>
            <div style={{
              display: 'inline-flex', alignItems: 'center', gap: '0.5rem', padding: '0.4rem 1rem',
              borderRadius: '999px', background: `${riskColor[risk.overall_risk]}22`,
              border: `1px solid ${riskColor[risk.overall_risk]}`,
              color: riskColor[risk.overall_risk], fontWeight: 700, fontSize: '0.82rem',
              marginBottom: '0.75rem',
            }}>
              Overall supply risk: {risk.overall_risk?.toUpperCase()}
            </div>
            <div style={{ display: 'flex', gap: '1.5rem', flexWrap: 'wrap', color: '#94a3b8', fontSize: '0.8rem' }}>
              <span>Uncovered: <strong style={{ color: '#fff' }}>{risk.uncovered_count}/{risk.total_ingredients}</strong> ingredients</span>
              <span>Single-sourced: <strong style={{ color: '#fff' }}>{risk.single_sourced_count}</strong></span>
              {risk.weighted_lead_time && <span>Avg lead: <strong style={{ color: '#fff' }}>{risk.weighted_lead_time}d</strong></span>}
            </div>
          </div>
        )}

        {pricing && pricing.length > 0 && (
          <div style={{ overflowX: 'auto' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.82rem' }}>
              <thead>
                <tr>{['Ingredient', '% Blend', 'Best Price', 'Supplier', 'Lead (days)', 'Suppliers', 'Risk', 'Source'].map(h => (
                  <th key={h} style={{ color: '#64748b', textAlign: 'left', padding: '0.35rem 0.6rem', borderBottom: '1px solid #1e3a5f', whiteSpace: 'nowrap', fontSize: '0.72rem' }}>{h}</th>
                ))}</tr>
              </thead>
              <tbody>
                {pricing.map((pr, i) => {
                  const ir = risk?.ingredient_risks?.find(r => r.ingredient === pr.ingredient)
                  const rc = riskColor[ir?.risk_level] || '#64748b'
                  return (
                    <tr key={i} style={{ borderBottom: '1px solid #0f2a4a' }}>
                      <td style={{ color: '#fff', padding: '0.4rem 0.6rem', fontWeight: 500 }}>{pr.ingredient}</td>
                      <td style={{ color: '#94a3b8', padding: '0.4rem 0.6rem', fontFamily: 'monospace' }}>{pr.pct?.toFixed(1)}%</td>
                      <td style={{ color: '#0D9488', padding: '0.4rem 0.6rem', fontFamily: 'monospace', fontWeight: 600 }}>
                        {pr.best_price_per_kg ? `$${pr.best_price_per_kg.toFixed(2)}` : '—'}
                      </td>
                      <td style={{ color: '#94a3b8', padding: '0.4rem 0.6rem', fontSize: '0.78rem' }}>{pr.best_supplier || 'DB estimate'}</td>
                      <td style={{ color: '#94a3b8', padding: '0.4rem 0.6rem' }}>{pr.lead_time_days ?? '—'}</td>
                      <td style={{ color: '#94a3b8', padding: '0.4rem 0.6rem', textAlign: 'center' }}>{pr.n_suppliers}</td>
                      <td style={{ padding: '0.4rem 0.6rem' }}>
                        {ir && <span style={{ color: rc, fontWeight: 600, fontSize: '0.78rem' }}>{ir.risk_level}</span>}
                      </td>
                      <td style={{ padding: '0.4rem 0.6rem' }}>
                        <span style={{ color: pr.source === 'supplier_listing' ? '#0D9488' : '#64748b', fontSize: '0.72rem' }}>
                          {pr.source === 'supplier_listing' ? 'live' : 'static'}
                        </span>
                      </td>
                    </tr>
                  )
                })}
              </tbody>
            </table>

            {risk?.geo_concentration && Object.keys(risk.geo_concentration).length > 0 && (
              <div style={{ marginTop: '1.25rem' }}>
                <div style={{ color: '#64748b', fontSize: '0.72rem', marginBottom: '0.5rem', textTransform: 'uppercase', letterSpacing: '0.07em' }}>Geographic concentration (%)</div>
                <ResponsiveContainer width="100%" height={140}>
                  <BarChart data={Object.entries(risk.geo_concentration).map(([c, p]) => ({ country: c, pct: p }))} margin={{ top: 0, right: 0, bottom: 20, left: 0 }}>
                    <XAxis dataKey="country" tick={{ fill: '#64748b', fontSize: 11 }} angle={-20} textAnchor="end" />
                    <YAxis tick={{ fill: '#64748b', fontSize: 11 }} unit="%" />
                    <Tooltip contentStyle={{ background: '#0D1F3C', border: '1px solid #1e3a5f', color: '#fff' }} formatter={v => [`${v}%`]} />
                    <Bar dataKey="pct" name="Coverage %" fill="#0D9488" radius={[4, 4, 0, 0]} />
                  </BarChart>
                </ResponsiveContainer>
              </div>
            )}
          </div>
        )}
      </div>

      {/* Prompt used */}
      {req.input_text && (
        <div style={card}>
          <div style={{ color: '#D97706', fontWeight: 700, fontSize: '0.75rem', letterSpacing: '0.07em', textTransform: 'uppercase', marginBottom: '0.5rem' }}>Original prompt</div>
          <div style={{ color: '#94a3b8', fontSize: '0.85rem', lineHeight: 1.7 }}>{req.input_text}</div>
        </div>
      )}
    </div>
  )
}
