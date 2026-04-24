import { useState } from 'react'
import { api } from '../api/client'
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts'

export default function Pareto() {
  const [vertical, setVertical] = useState('agricultural')
  const [nSolutions, setNSolutions] = useState(10)
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)

  const presets = {
    agricultural: { max_cost: 8, min_bio: 70, min_perf: 75 },
    food: { max_cost: 5, min_bio: 80, min_perf: 72 },
    fabric_laundry: { max_cost: 4, min_bio: 75, min_perf: 78 },
    personal_care: { max_cost: 5, min_bio: 70, min_perf: 75 },
    industrial: { max_cost: 6, min_bio: 55, min_perf: 72 },
    pharmaceutical: { max_cost: 12, min_bio: 55, min_perf: 72 },
    paint_coatings: { max_cost: 6, min_bio: 50, min_perf: 72 },
  }

  async function handleRun() {
    setLoading(true)
    setError(null)
    try {
      const res = await api.pareto({ vertical, constraints: presets[vertical] || {}, n_solutions: nSolutions })
      setResult(res.data)
    } catch (err) {
      setError(err.response?.data?.detail || err.message)
    } finally {
      setLoading(false)
    }
  }

  const chartData = result?.frontier?.map((s, i) => ({
    x: s.bio_pct ?? i,
    y: s.cost_per_kg ?? 0,
    name: `Solution ${i + 1}`
  })) || []

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        📊 Pareto Frontier
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Multi-objective optimization — cost vs eco score trade-off surface
      </p>

      <div style={{ display: 'flex', gap: '1rem', alignItems: 'flex-end', marginBottom: '1.5rem', flexWrap: 'wrap' }}>
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
            {['agricultural','food','fabric_laundry','personal_care','industrial','pharmaceutical','paint_coatings'].map(v => (
              <option key={v} value={v}>{v.replace('_',' ')}</option>
            ))}
          </select>
        </div>
        <div>
          <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>SOLUTIONS</label>
          <input
            type="number"
            value={nSolutions}
            onChange={e => setNSolutions(Number(e.target.value))}
            min={5} max={50}
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
          {loading ? 'Optimizing...' : 'Run Pareto →'}
        </button>
      </div>

      {error && (
        <div style={{
          padding: '1rem', background: '#450a0a', border: '1px solid #7f1d1d',
          borderRadius: '8px', color: '#fca5a5', fontSize: '0.85rem', marginBottom: '1rem'
        }}>
          {error}
        </div>
      )}

      {chartData.length > 0 ? (
        <div style={{
          background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1.5rem'
        }}>
          <div style={{ color: '#94a3b8', fontSize: '0.8rem', marginBottom: '1rem' }}>
            Bio-based % vs Cost — {chartData.length} Pareto-optimal solutions
          </div>
          <ResponsiveContainer width="100%" height={300}>
            <ScatterChart>
              <CartesianGrid strokeDasharray="3 3" stroke="#1e3a5f" />
              <XAxis dataKey="x" name="Bio %" stroke="#475569" tick={{ fill: '#64748b', fontSize: 11 }} label={{ value: 'Bio-based %', position: 'insideBottom', offset: -5, fill: '#64748b' }} />
              <YAxis dataKey="y" name="Cost" stroke="#475569" tick={{ fill: '#64748b', fontSize: 11 }} label={{ value: 'Cost', angle: -90, position: 'insideLeft', fill: '#64748b' }} />
              <Tooltip cursor={{ strokeDasharray: '3 3' }} contentStyle={{ background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '6px', color: '#fff' }} />
              <Scatter data={chartData} fill="#0D9488" />
            </ScatterChart>
          </ResponsiveContainer>
          {result?.recommended ? (
            <div style={{ marginTop: '1rem', color: '#94a3b8', fontSize: '0.82rem', lineHeight: 1.7 }}>
              Recommended balance: ${result.recommended.cost_per_kg}/kg · {result.recommended.bio_pct}% bio · {result.recommended.perf_score} performance.
            </div>
          ) : null}
        </div>
      ) : (
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#334155', textAlign: 'center', fontSize: '0.85rem'
        }}>
          Select vertical and click Run Pareto
        </div>
      )}
    </div>
  )
}
