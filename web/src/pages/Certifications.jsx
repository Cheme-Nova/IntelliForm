import { useMemo } from 'react'
import { loadLastRun } from '../lib/session'

const Badge = ({ name, verdict, confidence }) => {
  const isPass = String(verdict || '').includes('Pass')
  const isWarn = String(verdict || '').includes('Partial') || String(verdict || '').includes('Borderline')
  const color = isPass ? '#0D9488' : isWarn ? '#D97706' : '#ef4444'
  const bg = isPass ? '#0D948822' : isWarn ? '#D9770622' : '#ef444422'
  return (
    <div style={{
      background: bg, border: `1px solid ${color}`, borderRadius: '8px',
      padding: '1rem', flex: 1, minWidth: '180px'
    }}>
      <div style={{ color: '#94a3b8', fontSize: '0.7rem', marginBottom: '4px' }}>{name}</div>
      <div style={{ color, fontWeight: 700, fontSize: '1rem' }}>{verdict || 'Review'}</div>
      <div style={{ color: '#64748b', fontSize: '0.75rem', marginTop: '0.25rem' }}>{confidence || '—'} confidence</div>
    </div>
  )
}

export default function Certifications() {
  const session = useMemo(() => loadLastRun(), [])
  const result = session?.response?.cert

  if (!result) {
    return (
      <div style={{ maxWidth: '800px' }}>
        <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
          ✅ Certifications
        </h1>
        <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
          Certification posture from the most recent formulation session
        </p>
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#94a3b8', textAlign: 'center', fontSize: '0.9rem', lineHeight: 1.7
        }}>
          Run a formulation first to populate certification predictions.
        </div>
      </div>
    )
  }

  const predictions = Object.entries(result.predictions || {})

  return (
    <div style={{ maxWidth: '900px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        ✅ Certifications
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Certification predictions for the latest formulation session
      </p>

      <div style={{ marginBottom: '1rem', background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
        <div style={{ color: '#64748b', fontSize: '0.75rem', marginBottom: '0.35rem' }}>TOP CERTIFICATION PATH</div>
        <div style={{ color: '#f8fafc', fontWeight: 700, fontSize: '1.1rem' }}>{result.top_certification || '—'}</div>
        <div style={{ color: '#94a3b8', marginTop: '0.4rem', fontSize: '0.85rem' }}>
          Bio-based: {typeof result.bio_based_pct === 'number' ? `${result.bio_based_pct.toFixed(1)}%` : '—'} · Overall green score: {typeof result.overall_green_score === 'number' ? result.overall_green_score.toFixed(1) : '—'}
        </div>
      </div>

      <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap', marginBottom: '1rem' }}>
        {predictions.map(([name, prediction]) => (
          <Badge key={name} name={name} verdict={prediction.verdict} confidence={prediction.confidence} />
        ))}
      </div>

      {result.recommended_certs?.length ? (
        <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
          <div style={{ color: '#D97706', fontWeight: 600, marginBottom: '0.6rem', fontSize: '0.85rem' }}>
            Recommended next certifications
          </div>
          {result.recommended_certs.map((item) => (
            <div key={item} style={{ color: '#cbd5e1', fontSize: '0.86rem', marginBottom: '0.35rem' }}>
              · {item}
            </div>
          ))}
        </div>
      ) : null}
    </div>
  )
}
