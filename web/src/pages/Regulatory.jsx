import { useMemo } from 'react'
import { loadLastRun } from '../lib/session'

export default function Regulatory() {
  const session = useMemo(() => loadLastRun(), [])
  const reg = session?.response?.reg
  const vreg = session?.response?.vreg

  if (!reg || !vreg) {
    return (
      <div style={{ maxWidth: '860px' }}>
        <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
          📋 Regulatory
        </h1>
        <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
          Regulatory posture from the latest formulation session
        </p>
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#94a3b8', textAlign: 'center', fontSize: '0.9rem', lineHeight: 1.7
        }}>
          Run a formulation first to populate regulatory and vertical-compliance outputs.
        </div>
      </div>
    )
  }

  return (
    <div style={{ maxWidth: '860px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        📋 Regulatory
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Regulatory posture from the latest formulation session
      </p>

      <div style={{ display: 'grid', gap: '1rem', gridTemplateColumns: 'repeat(2, minmax(0, 1fr))', marginBottom: '1rem' }}>
        <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
          <div style={{ color: '#64748b', fontSize: '0.75rem', marginBottom: '0.35rem' }}>General status</div>
          <div style={{ color: '#f8fafc', fontWeight: 700 }}>{reg.overall_status || '—'}</div>
        </div>
        <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
          <div style={{ color: '#64748b', fontSize: '0.75rem', marginBottom: '0.35rem' }}>Vertical framework</div>
          <div style={{ color: '#f8fafc', fontWeight: 700 }}>{vreg.framework || '—'}</div>
        </div>
      </div>

      <div style={{ display: 'grid', gap: '1rem', gridTemplateColumns: 'repeat(2, minmax(0, 1fr))' }}>
        <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
          <div style={{ color: '#ef4444', fontWeight: 600, marginBottom: '0.6rem', fontSize: '0.85rem' }}>Amber / review flags</div>
          {(reg.amber_flags || []).length ? reg.amber_flags.map((item) => (
            <div key={item} style={{ color: '#fca5a5', fontSize: '0.86rem', marginBottom: '0.35rem' }}>· {item}</div>
          )) : (
            <div style={{ color: '#94a3b8', fontSize: '0.86rem' }}>No amber flags</div>
          )}
        </div>
        <div style={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }}>
          <div style={{ color: '#0D9488', fontWeight: 600, marginBottom: '0.6rem', fontSize: '0.85rem' }}>Vertical passes</div>
          {(vreg.passes || []).length ? vreg.passes.map((item) => (
            <div key={item} style={{ color: '#cbd5e1', fontSize: '0.86rem', marginBottom: '0.35rem' }}>· {item}</div>
          )) : (
            <div style={{ color: '#94a3b8', fontSize: '0.86rem' }}>No pass notes available</div>
          )}
        </div>
      </div>
    </div>
  )
}
