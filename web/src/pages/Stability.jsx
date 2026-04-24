import { useMemo } from 'react'
import { loadLastRun } from '../lib/session'

export default function Stability() {
  const session = useMemo(() => loadLastRun(), [])
  const result = session?.response?.stability

  if (!result) {
    return (
      <div style={{ maxWidth: '800px' }}>
        <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
          🧪 Stability
        </h1>
        <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
          Stability predictions from the latest formulation session
        </p>
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#94a3b8', textAlign: 'center', fontSize: '0.9rem', lineHeight: 1.7
        }}>
          Run a formulation first to populate stability signals.
        </div>
      </div>
    )
  }

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        🧪 Stability
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Stability predictions from the latest formulation session
      </p>

      <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap', marginBottom: '1rem' }}>
        {[
          ['OVERALL', result.overall_rating],
          ['SHELF LIFE', result.shelf_life_range],
          ['VISCOSITY', result.viscosity_range],
          ['PH RANGE', `${result.ph_min ?? '—'} – ${result.ph_max ?? '—'}`],
        ].map(([label, value]) => (
          <div key={label} style={{
            background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px',
            padding: '1.25rem', flex: 1, minWidth: '160px'
          }}>
            <div style={{ color: '#64748b', fontSize: '0.7rem', marginBottom: '4px' }}>{label}</div>
            <div style={{ color: '#fff', fontSize: '1.2rem', fontWeight: 700 }}>{value ?? '—'}</div>
          </div>
        ))}
      </div>

      <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem', marginBottom: '1rem' }}>
        <div style={{ color: '#D97706', fontWeight: 600, marginBottom: '0.6rem', fontSize: '0.85rem' }}>
          Packaging
        </div>
        <div style={{ color: '#cbd5e1', fontSize: '0.86rem' }}>{result.recommended_packaging || '—'}</div>
      </div>

      <div style={{ display: 'grid', gap: '1rem', gridTemplateColumns: 'repeat(2, minmax(0, 1fr))' }}>
        <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
          <div style={{ color: '#D97706', fontWeight: 600, marginBottom: '0.6rem', fontSize: '0.85rem' }}>Stability risks</div>
          {(result.stability_risks || []).map((item) => (
            <div key={item} style={{ color: '#cbd5e1', fontSize: '0.86rem', marginBottom: '0.35rem' }}>· {item}</div>
          ))}
        </div>
        <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
          <div style={{ color: '#D97706', fontWeight: 600, marginBottom: '0.6rem', fontSize: '0.85rem' }}>Boosters</div>
          {(result.stability_boosters || []).map((item) => (
            <div key={item} style={{ color: '#cbd5e1', fontSize: '0.86rem', marginBottom: '0.35rem' }}>· {item}</div>
          ))}
        </div>
      </div>
    </div>
  )
}
