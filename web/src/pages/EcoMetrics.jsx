import { useMemo } from 'react'
import { loadLastRun } from '../lib/session'

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

function fmt(value, digits = 1) {
  return typeof value === 'number' ? value.toFixed(digits) : '—'
}

export default function EcoMetrics() {
  const session = useMemo(() => loadLastRun(), [])
  const result = session?.response?.eco

  if (!result) {
    return (
      <div style={{ maxWidth: '800px' }}>
        <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
          🌱 EcoMetrics
        </h1>
        <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
          Sustainability scoring from the most recent successful formulation run
        </p>
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#94a3b8', textAlign: 'center', fontSize: '0.9rem', lineHeight: 1.7
        }}>
          Run a formulation in the Formulate tab first. EcoMetrics now reads the actual last run instead of recomputing from a placeholder prompt.
        </div>
      </div>
    )
  }

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        🌱 EcoMetrics
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Sustainability scoring for the latest formulation session
      </p>

      <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap', marginBottom: '1rem' }}>
        <CARD label="ECO SCORE" value={fmt(result.eco_score)} unit="/ 100" color="#0D9488" />
        <CARD label="BIODEGRADABILITY" value={fmt(result.biodegradability)} unit="/ 100" color="#D97706" />
        <CARD label="RENEWABILITY" value={fmt(result.renewability)} unit="/ 100" />
        <CARD label="GRADE" value={result.grade} />
      </div>

      <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
        <div style={{ color: '#D97706', fontWeight: 600, marginBottom: '0.75rem', fontSize: '0.85rem' }}>
          Category breakdown
        </div>
        {[
          ['Carbon Footprint', result.carbon_footprint],
          ['Ecotoxicity', result.ecotoxicity],
          ['Regulatory', result.regulatory],
        ].map(([label, value]) => (
          <div key={label} style={{ display: 'flex', justifyContent: 'space-between', color: '#cbd5e1', marginBottom: '0.45rem', fontSize: '0.86rem' }}>
            <span>{label}</span>
            <strong>{fmt(value)}</strong>
          </div>
        ))}
      </div>
    </div>
  )
}
