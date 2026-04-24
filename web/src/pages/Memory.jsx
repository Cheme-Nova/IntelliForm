import { useState, useEffect } from 'react'
import { api } from '../api/client'
import { useAuth } from '../auth/AuthContext'

export default function Memory() {
  const auth = useAuth()
  const [records, setRecords] = useState([])
  const [loading, setLoading] = useState(false)
  const [n, setN] = useState(20)

  async function fetchMemory() {
    setLoading(true)
    try {
      const res = auth?.user ? await api.projects(n) : await api.memory(n)
      setRecords(res.data || [])
    } catch (err) {
      setRecords([])
    } finally {
      setLoading(false)
    }
  }

  useEffect(() => { fetchMemory() }, [auth?.user])

  function normalizeRecord(rec) {
    if (auth?.user) {
      const payload = {
        application: rec.application,
        blend: rec.blend,
        cost_per_kg: rec.cost_per_kg,
        bio_pct: rec.bio_pct,
        perf_score: rec.perf_score,
        eco_score: rec.eco_score,
        eco_grade: rec.eco_grade,
        optimizer: rec.optimizer,
        parser: rec.parser,
        relaxed: rec.relaxed,
        nl_input: rec.nl_input,
      }
      return {
        label: rec.application || 'project',
        timestamp: rec.created_at || '',
        vertical: rec.application,
        payload,
      }
    }

    return {
      label: rec.type || 'event',
      timestamp: rec.timestamp || '',
      vertical: rec.vertical,
      payload: rec.payload || {},
    }
  }

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        🧠 Memory
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        {auth?.user
          ? 'Your saved formulation runs from the public free account'
          : 'Anonymous local memory feed. Sign in to save account-linked formulation history.'}
      </p>

      <div style={{ display: 'flex', gap: '1rem', alignItems: 'flex-end', marginBottom: '1.5rem' }}>
        <div>
          <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>SHOW LAST</label>
          <input
            type="number"
            value={n}
            onChange={e => setN(Number(e.target.value))}
            min={5} max={100}
            style={{
              width: '80px', background: '#0D1F3C', border: '1px solid #1e3a5f',
              borderRadius: '6px', color: '#fff', padding: '0.5rem', fontSize: '0.85rem'
            }}
          />
        </div>
        <button
          onClick={fetchMemory}
          disabled={loading}
          style={{
            background: loading ? '#334155' : '#0D9488', color: '#fff',
            border: 'none', borderRadius: '8px', padding: '0.6rem 1.5rem',
            fontSize: '0.9rem', fontWeight: 600, cursor: loading ? 'not-allowed' : 'pointer'
          }}
        >
          {loading ? 'Loading...' : 'Refresh →'}
        </button>
      </div>

      {records.length === 0 ? (
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#334155', textAlign: 'center', fontSize: '0.85rem'
        }}>
          No runs recorded yet. Run a formulation to see it here.
        </div>
      ) : (
        <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem' }}>
          {records.slice().reverse().map((rec, i) => {
            const item = normalizeRecord(rec)
            const preview = JSON.stringify(item.payload, null, 2) || ''
            return (
            <div key={i} style={{
              background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem'
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '0.5rem' }}>
                <span style={{
                  background: '#0D948822', color: '#0D9488', padding: '2px 8px',
                  borderRadius: '4px', fontSize: '0.75rem', fontWeight: 600
                }}>
                  {item.label}
                </span>
                <span style={{ color: '#334155', fontSize: '0.75rem' }}>{item.timestamp}</span>
              </div>
              {item.vertical && (
                <div style={{ color: '#64748b', fontSize: '0.75rem', marginBottom: '4px' }}>
                  Vertical: {item.vertical}
                </div>
              )}
              <pre style={{ color: '#475569', fontSize: '0.75rem', overflow: 'auto', margin: 0 }}>
                {preview.slice(0, 400)}
                {preview.length > 400 ? '...' : ''}
              </pre>
            </div>
          )})}
        </div>
      )}
    </div>
  )
}
