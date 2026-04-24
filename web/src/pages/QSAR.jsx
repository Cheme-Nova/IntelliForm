import { useState } from 'react'
import { api } from '../api/client'

export default function QSAR() {
  const [smiles, setSmiles] = useState('CC(C)CC1=CC=C(C=C1)C(C)C(=O)O\nO=C(O)c1ccccc1')
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)

  async function handleRun() {
    setLoading(true)
    setError(null)
    try {
      const smilesList = smiles.split('\n').map(s => s.trim()).filter(Boolean)
      const res = await api.qsar({ smiles: smilesList, properties: ['biodegradability', 'ecotox', 'performance'] })
      setResult(res.data)
    } catch (err) {
      setError(err.response?.data?.detail || err.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>
        🔬 QSAR Predictor
      </h1>
      <p style={{ color: '#64748b', marginBottom: '2rem', fontSize: '0.9rem' }}>
        Predict biodegradability, ecotoxicity, and performance from SMILES structures
      </p>

      <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem', marginBottom: '1.5rem' }}>
        <div>
          <label style={{ color: '#64748b', fontSize: '0.75rem', display: 'block', marginBottom: '4px' }}>
            SMILES (one per line)
          </label>
          <textarea
            value={smiles}
            onChange={e => setSmiles(e.target.value)}
            rows={5}
            style={{
              width: '100%', background: '#0D1F3C', border: '1px solid #1e3a5f',
              borderRadius: '8px', color: '#fff', padding: '1rem', fontSize: '0.85rem',
              fontFamily: 'monospace', boxSizing: 'border-box', resize: 'vertical'
            }}
          />
        </div>
        <button
          onClick={handleRun}
          disabled={loading}
          style={{
            background: loading ? '#334155' : '#0D9488', color: '#fff',
            border: 'none', borderRadius: '8px', padding: '0.6rem 1.5rem',
            fontSize: '0.9rem', fontWeight: 600, cursor: loading ? 'not-allowed' : 'pointer',
            alignSelf: 'flex-start'
          }}
        >
          {loading ? 'Predicting...' : 'Predict Properties →'}
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

      {result ? (
        <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem' }}>
          {result.predictions?.map((pred, i) => (
            <div key={i} style={{
              background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem'
            }}>
              <div style={{ color: '#D97706', fontSize: '0.75rem', fontFamily: 'monospace', marginBottom: '0.5rem' }}>
                {pred.smiles}
              </div>
              <div style={{ display: 'flex', gap: '1rem', flexWrap: 'wrap' }}>
                {[
                  ['biodegradability', pred.biodegradability],
                  ['ecotoxicity', pred.ecotoxicity],
                  ['performance', pred.performance],
                ].map(([key, val]) => (
                  <div key={key}>
                    <div style={{ color: '#64748b', fontSize: '0.7rem' }}>{key.toUpperCase()}</div>
                    <div style={{ color: '#0D9488', fontWeight: 600 }}>
                      {typeof val === 'number' ? val.toFixed(3) : val}
                    </div>
                  </div>
                ))}
              </div>
              {pred.confidence ? (
                <div style={{ color: '#94a3b8', fontSize: '0.75rem', marginTop: '0.75rem' }}>
                  Confidence: {pred.confidence}
                </div>
              ) : null}
            </div>
          ))}
        </div>
      ) : (
        <div style={{
          padding: '2rem', background: '#0D1F3C', border: '1px solid #1e3a5f',
          borderRadius: '8px', color: '#334155', textAlign: 'center', fontSize: '0.85rem'
        }}>
          Enter SMILES above and click Predict Properties
        </div>
      )}
    </div>
  )
}
