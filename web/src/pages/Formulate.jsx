import { useState } from 'react'
import { api } from '../api/client'

const VERTICALS = ['personal_care', 'home_care', 'industrial', 'pharma', 'food', 'agriculture']

const Card = ({ label, value, unit, color }) => (
  <div style={{
    background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px',
    padding: '1.25rem', flex: 1, minWidth: '130px'
  }}>
    <div style={{ color: '#64748b', fontSize: '0.7rem', marginBottom: '4px' }}>{label}</div>
    <div style={{ color: color || '#fff', fontSize: '1.5rem', fontWeight: 700 }}>{value ?? '—'}</div>
    {unit && <div style={{ color: '#475569', fontSize: '0.7rem' }}>{unit}</div>}
  </div>
)

const Section = ({ title, children }) => (
  <div style={{
    background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px',
    padding: '1.25rem', marginBottom: '1rem'
  }}>
    <div style={{ color: '#0D9488', fontWeight: 700, fontSize: '0.85rem', marginBottom: '1rem' }}>
      {title}
    </div>
    {children}
  </div>
)

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
    <div style={{ maxWidth: '860px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        ⚗️ Formulate
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Describe your formulation goal in plain language
      </p>

      <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem', marginBottom: '1.5rem' }}>
        <textarea
          value={inputText}
          onChange={e => setInputText(e.target.value)}
          placeholder="e.g. Create a sulfate-free shampoo with >80% bio-based content, mild pH, stable foam..."
          rows={4}
          style={{
            background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px',
            color: '#fff', padding: '1rem', fontSize: '0.9rem', resize: 'vertical',
            outline: 'none', fontFamily: 'inherit', width: '100%', boxSizing: 'border-box'
          }}
        />
        <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap' }}>
          <div style={{ flex: 1, minWidth: '160px' }}>
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>VERTICAL</label>
            <select value={vertical} onChange={e => setVertical(e.target.value)} style={{
              width: '100%', background: '#0D1F3C', border: '1px solid #1e3a5f',
              borderRadius: '6px', color: '#fff', padding: '0.5rem', fontSize: '0.85rem'
            }}>
              {VERTICALS.map(v => <option key={v} value={v}>{v.replace('_', ' ')}</option>)}
            </select>
          </div>
          <div style={{ flex: 1, minWidth: '120px' }}>
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>BATCH SIZE (kg)</label>
            <input type="number" value={batchSize} onChange={e => setBatchSize(Number(e.target.value))} style={{
              width: '100%', background: '#0D1F3C', border: '1px solid #1e3a5f',
              borderRadius: '6px', color: '#fff', padding: '0.5rem', fontSize: '0.85rem', boxSizing: 'border-box'
            }} />
          </div>
          <div style={{ flex: 1, minWidth: '120px' }}>
            <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>OPT MODE</label>
            <select value={optMode} onChange={e => setOptMode(e.target.value)} style={{
              width: '100%', background: '#0D1F3C', border: '1px solid #1e3a5f',
              borderRadius: '6px', color: '#fff', padding: '0.5rem', fontSize: '0.85rem'
            }}>
              <option value="auto">Auto</option>
              <option value="pareto">Pareto</option>
              <option value="bayesian">Bayesian</option>
              <option value="single">Single</option>
            </select>
          </div>
        </div>
        <button onClick={handleSubmit} disabled={loading || !inputText.trim()} style={{
          background: loading ? '#334155' : '#0D9488', color: '#fff', border: 'none',
          borderRadius: '8px', padding: '0.75rem 2rem', fontSize: '0.95rem', fontWeight: 600,
          cursor: loading ? 'not-allowed' : 'pointer', alignSelf: 'flex-start'
        }}>
          {loading ? 'Optimizing...' : 'Run IntelliForm →'}
        </button>
      </div>

      {error && (
        <div style={{
          padding: '1rem', background: '#450a0a', border: '1px solid #7f1d1d',
          borderRadius: '8px', color: '#fca5a5', fontSize: '0.85rem', marginBottom: '1rem'
        }}>{error}</div>
      )}

      {result && (
        <div>
          {/* Key Metrics */}
          <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap', marginBottom: '1rem' }}>
            <Card label="COST/KG" value={`$${result.result?.cost_per_kg?.toFixed(2)}`} color="#0D9488" />
            <Card label="BIO-BASED" value={`${result.result?.bio_pct?.toFixed(1)}%`} color="#D97706" />
            <Card label="PERF SCORE" value={result.result?.perf_score?.toFixed(1)} unit="/ 100" />
            <Card label="ECO GRADE" value={result.eco?.grade} color="#0D9488" />
            <Card label="ECO SCORE" value={result.eco?.eco_score?.toFixed(1)} unit="/ 100" />
          </div>

          {/* Blend */}
          <Section title="🧪 Optimized Blend">
            <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
              {Object.entries(result.result?.blend || {}).map(([ing, pct]) => (
                <div key={ing} style={{ display: 'flex', alignItems: 'center', gap: '1rem' }}>
                  <div style={{ color: '#94a3b8', fontSize: '0.85rem', minWidth: '220px' }}>{ing}</div>
                  <div style={{ flex: 1, background: '#1e3a5f', borderRadius: '4px', height: '8px' }}>
                    <div style={{ width: `${Math.min(pct, 100)}%`, background: '#0D9488', height: '8px', borderRadius: '4px' }} />
                  </div>
                  <div style={{ color: '#fff', fontSize: '0.85rem', fontWeight: 600, minWidth: '50px', textAlign: 'right' }}>
                    {pct}%
                  </div>
                </div>
              ))}
            </div>
          </Section>

          {/* Carbon */}
          <Section title="♻️ Carbon Impact">
            <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap' }}>
              <Card label="CO₂ DISPLACED" value={result.carbon?.co2_displaced_tonnes?.toFixed(2)} unit="tonnes per batch" color="#0D9488" />
              <Card label="CREDIT VALUE" value={`$${result.carbon?.credit_value_mid?.toFixed(0)}`} unit="mid-market" color="#D97706" />
              <Card label="ANNUAL VALUE" value={`$${result.carbon?.annual_credit_value_mid?.toFixed(0)}`} unit="at 12 batches/yr" />
            </div>
            {result.carbon?.summary && (
              <div style={{ color: '#64748b', fontSize: '0.8rem', marginTop: '1rem', lineHeight: 1.6 }}>
                {result.carbon.summary}
              </div>
            )}
          </Section>

          {/* Certifications */}
          <Section title="✅ Certification Predictions">
            <div style={{ display: 'flex', gap: '0.75rem', flexWrap: 'wrap' }}>
              {Object.entries(result.cert?.predictions || {}).map(([name, pred]) => {
                const pass = pred.pass_probability >= 0.6
                const color = pass ? '#0D9488' : '#ef4444'
                const bg = pass ? '#0D948818' : '#ef444418'
                return (
                  <div key={name} style={{
                    background: bg, border: `1px solid ${color}`, borderRadius: '8px',
                    padding: '0.75rem 1rem', minWidth: '130px'
                  }}>
                    <div style={{ color: '#94a3b8', fontSize: '0.7rem', marginBottom: '2px' }}>
                      {name.replace(/_/g, ' ')}
                    </div>
                    <div style={{ color, fontWeight: 700, fontSize: '0.9rem' }}>{pred.verdict}</div>
                    <div style={{ color: '#475569', fontSize: '0.7rem' }}>
                      {(pred.pass_probability * 100).toFixed(0)}% probability
                    </div>
                  </div>
                )
              })}
            </div>
          </Section>

          {/* Agent Insights */}
          <Section title="🤖 Agent Insights">
            {result.agents?.map((a, i) => (
              <div key={i} style={{
                color: '#94a3b8', fontSize: '0.85rem', marginBottom: '0.5rem',
                paddingBottom: '0.5rem', borderBottom: i < result.agents.length - 1 ? '1px solid #1e3a5f' : 'none'
              }}>
                {a}
              </div>
            ))}
          </Section>

          {/* Stability */}
          <Section title="🧪 Stability">
            <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap' }}>
              <Card label="SHELF LIFE" value={result.stability?.shelf_life_range} />
              <Card label="VISCOSITY" value={result.stability?.viscosity_range} />
              <Card label="pH RANGE" value={`${result.stability?.ph_min}–${result.stability?.ph_max}`} />
              <Card label="RATING" value={result.stability?.overall_rating}
                color={result.stability?.overall_rating === 'Good' ? '#0D9488' : '#D97706'} />
            </div>
          </Section>

          {/* Regulatory */}
          <Section title="📋 Regulatory Status">
            <div style={{ color: result.vreg?.overall_status?.includes('✅') ? '#0D9488' : '#D97706', fontWeight: 600, marginBottom: '0.5rem' }}>
              {result.vreg?.overall_status}
            </div>
            <div style={{ color: '#64748b', fontSize: '0.8rem', marginBottom: '0.75rem' }}>
              {result.vreg?.framework}
            </div>
            {result.reg?.amber_flags?.length > 0 && (
              <div>
                {result.reg.amber_flags.slice(0, 3).map((f, i) => (
                  <div key={i} style={{ color: '#D97706', fontSize: '0.8rem', marginBottom: '2px' }}>⚠️ {f}</div>
                ))}
              </div>
            )}
          </Section>
        </div>
      )}
    </div>
  )
}
