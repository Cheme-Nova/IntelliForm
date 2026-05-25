import { useMemo, useState } from 'react'
import { loadLastRun } from '../lib/session'
import { api } from '../api/client'

const c = {
  ink: '#111827', muted: '#5f6f7d', accent: '#1f8a7a',
  accentSoft: '#eef6f3', steel: '#4c6375', line: '#dce6ee',
  amber: '#b5893a', navy: '#0A1628',
}

const panel = {
  background: '#ffffff', border: `1px solid ${c.line}`,
  borderRadius: '8px', boxShadow: '0 12px 32px rgba(36,49,61,0.08)',
}
const label = { color: c.steel, fontSize: '0.7rem', letterSpacing: '0.08em', textTransform: 'uppercase', fontWeight: 800 }

function MetricCard({ lbl, value, unit, tone = c.ink }) {
  return (
    <div style={{ ...panel, padding: '1rem', flex: '1 1 140px', minWidth: 'min(140px,100%)', boxShadow: 'none' }}>
      <div style={label}>{lbl}</div>
      <div style={{ color: tone, fontSize: '1.4rem', fontWeight: 840, marginTop: '0.35rem' }}>{value ?? '—'}</div>
      {unit && <div style={{ color: c.muted, fontSize: '0.72rem', marginTop: '0.1rem' }}>{unit}</div>}
    </div>
  )
}

// EU ETS reference price — update quarterly
const ETS_RATE = 65.0  // €/tCO2e, March 2026

export default function Carbon() {
  const session = useMemo(() => loadLastRun(), [])
  const result  = session?.response?.carbon
  const blend   = session?.response?.result?.blend || {}
  const batchKg = session?.request?.batchSize || 500

  const [etsPriceEur, setEtsPriceEur] = useState(ETS_RATE)
  const [passportLoading, setPassportLoading] = useState(false)
  const [passportErr, setPassportErr] = useState(null)

  async function downloadPassport() {
    if (!Object.keys(blend).length) return
    setPassportLoading(true)
    setPassportErr(null)
    try {
      const res = await api.carbonPassport({ blend, batch_kg: batchKg, product_name: 'IntelliForm Formulation' })
      const json = res.data?.json || JSON.stringify(res.data, null, 2)
      const blob = new Blob([json], { type: 'application/json' })
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a'); a.href = url; a.download = 'carbon_passport.json'; a.click()
      URL.revokeObjectURL(url)
    } catch (e) {
      setPassportErr(e.response?.data?.detail || e.message)
    } finally {
      setPassportLoading(false)
    }
  }

  if (!result) {
    return (
      <div style={{ maxWidth: '860px' }}>
        <h1 style={{ color: c.accent, fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>Carbon</h1>
        <p style={{ color: c.muted, marginBottom: '2rem', fontSize: '0.9rem' }}>Carbon, CBAM, and circularity outputs from the latest formulation session</p>
        <div style={{ ...panel, padding: '2rem', textAlign: 'center', color: c.muted, fontSize: '0.9rem', lineHeight: 1.7, boxShadow: 'none' }}>
          Run a formulation first to generate carbon and circularity outputs.
        </div>
      </div>
    )
  }

  // CBAM calculation
  const greenCo2PerKg   = result.green_co2_per_kg   ?? 0
  const baselineCo2PerKg = result.baseline_co2_per_kg ?? 0
  const greenCo2Total   = greenCo2PerKg   * batchKg / 1000   // tCO2e
  const baselineCo2Total = baselineCo2PerKg * batchKg / 1000
  const cbamCost        = greenCo2Total * etsPriceEur
  const baselineCbam    = baselineCo2Total * etsPriceEur
  const cbamSaving      = baselineCbam - cbamCost
  const cbamPerKg       = cbamCost / (batchKg || 1)

  return (
    <div style={{ maxWidth: '860px' }}>
      <h1 style={{ color: c.accent, fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>Carbon</h1>
      <p style={{ color: c.muted, marginBottom: '1.5rem', fontSize: '0.9rem' }}>Carbon, CBAM, and circularity outputs from the latest formulation session</p>

      {/* Key metrics */}
      <div style={{ display: 'flex', gap: '0.65rem', flexWrap: 'wrap', marginBottom: '1rem' }}>
        <MetricCard lbl="CO2 Displaced"    value={typeof result.co2_displaced_kg === 'number' ? result.co2_displaced_kg.toFixed(1) : '—'} unit="kg per batch" tone={c.accent} />
        <MetricCard lbl="Credits / Batch"  value={typeof result.credits_per_batch === 'number' ? result.credits_per_batch.toFixed(3) : '—'} unit="verified credits" tone={c.amber} />
        <MetricCard lbl="Credit Value"     value={typeof result.credit_value_mid === 'number' ? `$${result.credit_value_mid.toFixed(2)}` : '—'} unit="mid-price USD" />
        <MetricCard lbl="Circular Grade"   value={result.circular_grade} />
      </div>

      {/* Summary */}
      <div style={{ ...panel, padding: '1rem 1.1rem', marginBottom: '1.1rem', boxShadow: 'none' }}>
        <div style={label}>Summary</div>
        <div style={{ color: '#3d5060', marginTop: '0.5rem', lineHeight: 1.75, fontSize: '0.88rem' }}>{result.summary}</div>
      </div>

      {/* CO2e comparison */}
      <div style={{ display: 'grid', gap: '0.65rem', gridTemplateColumns: 'repeat(2, 1fr)', marginBottom: '1.1rem' }}>
        <div style={{ ...panel, padding: '1rem', boxShadow: 'none' }}>
          <div style={label}>Green formulation CO2e</div>
          <div style={{ color: c.accent, fontWeight: 780, fontSize: '1.1rem', marginTop: '0.3rem' }}>{result.green_co2_per_kg} kg/kg</div>
        </div>
        <div style={{ ...panel, padding: '1rem', boxShadow: 'none' }}>
          <div style={label}>Petro baseline CO2e</div>
          <div style={{ color: c.ink, fontWeight: 780, fontSize: '1.1rem', marginTop: '0.3rem' }}>{result.baseline_co2_per_kg} kg/kg</div>
        </div>
      </div>

      {/* ── CBAM Calculator ─────────────────────────────────────────────────── */}
      <div style={{ ...panel, padding: '1.25rem', marginBottom: '1.1rem' }}>
        <div style={label}>CBAM Calculator</div>
        <h2 style={{ margin: '0.35rem 0 0.5rem', color: c.ink, fontSize: '1.05rem', fontWeight: 800 }}>EU Carbon Border Adjustment Mechanism</h2>
        <div style={{ color: c.muted, fontSize: '0.83rem', marginBottom: '1rem', lineHeight: 1.6 }}>
          EU Regulation 2023/956 — mandatory CO2 cost for chemical imports into the EU from 2026. Adjust the EU ETS price below.
        </div>

        <div style={{ display: 'flex', alignItems: 'center', gap: '1rem', marginBottom: '1rem', flexWrap: 'wrap' }}>
          <div>
            <div style={label}>EU ETS price (€/tCO2e)</div>
            <input
              type="number" min={10} max={200} step={1}
              value={etsPriceEur}
              onChange={e => setEtsPriceEur(Number(e.target.value))}
              style={{ marginTop: '0.4rem', padding: '0.5rem 0.75rem', border: `1px solid ${c.line}`, borderRadius: '6px', fontSize: '0.9rem', width: '120px', color: c.ink }}
            />
          </div>
          <div style={{ color: c.muted, fontSize: '0.78rem', marginTop: '1rem' }}>
            Reference: €{ETS_RATE}/tCO2e (March 2026)
          </div>
        </div>

        <div style={{ display: 'flex', gap: '0.65rem', flexWrap: 'wrap' }}>
          <div style={{ ...panel, padding: '0.85rem 1rem', flex: '1 1 130px', boxShadow: 'none', background: c.accentSoft, borderColor: '#b9d1ca' }}>
            <div style={label}>CBAM cost (green)</div>
            <div style={{ color: c.accent, fontWeight: 780, fontSize: '1.2rem', marginTop: '0.3rem' }}>€{cbamCost.toFixed(2)}</div>
            <div style={{ color: c.muted, fontSize: '0.72rem' }}>per {batchKg} kg batch</div>
          </div>
          <div style={{ ...panel, padding: '0.85rem 1rem', flex: '1 1 130px', boxShadow: 'none', background: '#fffdf0', borderColor: '#e8c96b' }}>
            <div style={label}>CBAM cost (petro baseline)</div>
            <div style={{ color: c.amber, fontWeight: 780, fontSize: '1.2rem', marginTop: '0.3rem' }}>€{baselineCbam.toFixed(2)}</div>
            <div style={{ color: c.muted, fontSize: '0.72rem' }}>per {batchKg} kg batch</div>
          </div>
          <div style={{ ...panel, padding: '0.85rem 1rem', flex: '1 1 130px', boxShadow: 'none' }}>
            <div style={label}>CBAM saving</div>
            <div style={{ color: c.accent, fontWeight: 780, fontSize: '1.2rem', marginTop: '0.3rem' }}>€{cbamSaving.toFixed(2)}</div>
            <div style={{ color: c.muted, fontSize: '0.72rem' }}>vs petrochemical blend</div>
          </div>
          <div style={{ ...panel, padding: '0.85rem 1rem', flex: '1 1 110px', boxShadow: 'none' }}>
            <div style={label}>CBAM per kg</div>
            <div style={{ color: c.ink, fontWeight: 780, fontSize: '1.2rem', marginTop: '0.3rem' }}>€{cbamPerKg.toFixed(4)}</div>
            <div style={{ color: c.muted, fontSize: '0.72rem' }}>€/kg product</div>
          </div>
        </div>
      </div>

      {/* ── Carbon Passport download ─────────────────────────────────────────── */}
      <div style={{ ...panel, padding: '1rem 1.1rem', boxShadow: 'none', background: c.accentSoft, borderColor: '#b9d1ca' }}>
        <div style={label}>ISO 14067 Carbon Passport</div>
        <div style={{ color: '#3d5060', fontSize: '0.84rem', lineHeight: 1.6, marginTop: '0.4rem', marginBottom: '0.85rem' }}>
          Generate a machine-readable, EU CBAM-compliant carbon footprint declaration for this formulation. Download as JSON — ready for submission.
        </div>
        <button
          onClick={downloadPassport}
          disabled={passportLoading || !Object.keys(blend).length}
          style={{ background: passportLoading ? '#9eafa8' : c.accent, color: '#fff', border: 'none', borderRadius: '6px', padding: '0.65rem 1.1rem', fontWeight: 700, cursor: passportLoading ? 'not-allowed' : 'pointer', fontSize: '0.88rem' }}
        >
          {passportLoading ? 'Generating...' : '↓ Download Carbon Passport (JSON)'}
        </button>
        {passportErr && <div style={{ color: '#b45858', fontSize: '0.8rem', marginTop: '0.5rem' }}>{passportErr}</div>}
      </div>
    </div>
  )
}
