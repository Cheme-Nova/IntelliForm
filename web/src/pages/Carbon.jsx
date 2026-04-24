import { useMemo } from 'react'
import { loadLastRun } from '../lib/session'

function card(label, value, unit, color = '#fff') {
  return (
    <div key={label} style={{
      background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px',
      padding: '1.25rem', flex: 1, minWidth: '160px'
    }}>
      <div style={{ color: '#64748b', fontSize: '0.7rem', marginBottom: '4px' }}>{label}</div>
      <div style={{ color, fontSize: '1.4rem', fontWeight: 700 }}>{value ?? '—'}</div>
      <div style={{ color: '#475569', fontSize: '0.7rem' }}>{unit}</div>
    </div>
  )
}

export default function Carbon() {
  const session = useMemo(() => loadLastRun(), [])
  const result = session?.response?.carbon

  if (!result) {
    return (
      <div style={{ maxWidth: '800px' }}>
        <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
          ♻️ Carbon
        </h1>
        <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
          Carbon and circularity outputs from the latest formulation session
        </p>
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#94a3b8', textAlign: 'center', fontSize: '0.9rem', lineHeight: 1.7
        }}>
          Run a formulation first to generate carbon and circularity outputs.
        </div>
      </div>
    )
  }

  return (
    <div style={{ maxWidth: '860px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        ♻️ Carbon
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Carbon and circularity outputs from the latest formulation session
      </p>

      <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap', marginBottom: '1rem' }}>
        {[
          card('CO2 DISPLACED', typeof result.co2_displaced_kg === 'number' ? result.co2_displaced_kg.toFixed(1) : '—', 'kg per batch', '#0D9488'),
          card('CREDITS / BATCH', typeof result.credits_per_batch === 'number' ? result.credits_per_batch.toFixed(3) : '—', 'verified credits', '#D97706'),
          card('MID CREDIT VALUE', typeof result.credit_value_mid === 'number' ? `$${result.credit_value_mid.toFixed(2)}` : '—', 'USD'),
          card('CIRCULAR GRADE', result.circular_grade, 'circularity'),
        ]}
      </div>

      <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem', marginBottom: '1rem' }}>
        <div style={{ color: '#D97706', fontWeight: 600, marginBottom: '0.6rem', fontSize: '0.85rem' }}>
          Summary
        </div>
        <div style={{ color: '#cbd5e1', fontSize: '0.88rem', lineHeight: 1.7 }}>{result.summary}</div>
      </div>

      <div style={{ display: 'grid', gap: '1rem', gridTemplateColumns: 'repeat(2, minmax(0, 1fr))' }}>
        <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
          <div style={{ color: '#64748b', fontSize: '0.75rem', marginBottom: '0.35rem' }}>Green formulation CO2e</div>
          <div style={{ color: '#f8fafc', fontWeight: 700 }}>{result.green_co2_per_kg} kg/kg</div>
        </div>
        <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
          <div style={{ color: '#64748b', fontSize: '0.75rem', marginBottom: '0.35rem' }}>Petro baseline CO2e</div>
          <div style={{ color: '#f8fafc', fontWeight: 700 }}>{result.baseline_co2_per_kg} kg/kg</div>
        </div>
      </div>
    </div>
  )
}
