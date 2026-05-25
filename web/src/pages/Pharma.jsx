import { useState } from 'react'
import { loadLastRun } from '../lib/session'
import { api } from '../api/client'

const c = {
  ink: '#111827', muted: '#5f6f7d', accent: '#1f8a7a',
  accentSoft: '#eef6f3', steel: '#4c6375', line: '#dce6ee',
  amber: '#b5893a', red: '#b45858', navy: '#0A1628',
}

const panel = {
  background: '#ffffff', border: `1px solid ${c.line}`,
  borderRadius: '8px', boxShadow: '0 12px 32px rgba(36,49,61,0.08)',
}
const lbl = { color: c.steel, fontSize: '0.7rem', letterSpacing: '0.08em', textTransform: 'uppercase', fontWeight: 800 }

const BCS_CLASSES = ['I', 'II', 'III', 'IV']
const DOSAGE_FORMS = [
  { value: 'immediate_release_tablet', label: 'Immediate Release Tablet' },
  { value: 'modified_release_tablet', label: 'Modified Release Tablet' },
  { value: 'hard_gelatin_capsule', label: 'Hard Gelatin Capsule' },
  { value: 'oral_solution', label: 'Oral Solution / Syrup' },
  { value: 'enteric_coated_tablet', label: 'Enteric Coated Tablet' },
  { value: 'topical_cream', label: 'Topical Cream / Ointment' },
]
const MARKETS = ['USA', 'EU', 'Japan', 'India', 'Southeast Asia', 'Sub-Saharan Africa', 'Brazil']

const SEV_COLOR = { Severe: c.red, Moderate: c.amber, Mild: c.steel }

function Badge({ text, color }) {
  return (
    <span style={{
      background: `${color}18`, color, padding: '2px 8px',
      borderRadius: '4px', fontSize: '0.72rem', fontWeight: 700,
    }}>{text}</span>
  )
}

export default function Pharma() {
  const session = loadLastRun()
  const blend = session?.response?.result?.blend || {}

  const [bcsClass, setBcsClass] = useState('I')
  const [dosageForm, setDosageForm] = useState('immediate_release_tablet')
  const [targetMarkets, setTargetMarkets] = useState(['USA', 'EU'])
  const [isGeneric, setIsGeneric] = useState(false)
  const [isPediatric, setIsPediatric] = useState(false)
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)

  const hasBlend = Object.keys(blend).length > 0

  function toggleMarket(m) {
    setTargetMarkets(prev => prev.includes(m) ? prev.filter(x => x !== m) : [...prev, m])
  }

  async function runAnalysis() {
    if (!hasBlend) return
    setLoading(true)
    setError(null)
    setResult(null)
    try {
      const res = await api.pharmaDeepDive({
        blend,
        bcs_class: bcsClass,
        dosage_form: dosageForm,
        target_markets: targetMarkets,
        is_generic: isGeneric,
        is_pediatric: isPediatric,
      })
      setResult(res.data)
    } catch (e) {
      setError(e.response?.data?.detail || e.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={{ maxWidth: '860px' }}>
      <h1 style={{ color: c.accent, fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>Pharma Deep Dive</h1>
      <p style={{ color: c.muted, marginBottom: '1.5rem', fontSize: '0.9rem' }}>
        ICH Q8 / BCS pharmaceutical analysis — excipient compatibility, stability zones, regulatory pathway
      </p>

      {/* Configuration */}
      <div style={{ ...panel, padding: '1.25rem', marginBottom: '1.1rem' }}>
        <div style={lbl}>Analysis Configuration</div>

        {!hasBlend && (
          <div style={{ marginTop: '0.75rem', color: c.muted, fontSize: '0.85rem', background: '#fff8f0', border: `1px solid #f0d090`, borderRadius: '6px', padding: '0.75rem 1rem' }}>
            Run a formulation first to load a blend for analysis.
          </div>
        )}

        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(200px, 1fr))', gap: '1rem', marginTop: '1rem' }}>
          {/* BCS Class */}
          <div>
            <div style={lbl}>BCS Class</div>
            <div style={{ display: 'flex', gap: '0.4rem', marginTop: '0.5rem' }}>
              {BCS_CLASSES.map(cls => (
                <button key={cls} onClick={() => setBcsClass(cls)} style={{
                  flex: 1, padding: '0.45rem', border: `1px solid ${bcsClass === cls ? c.accent : c.line}`,
                  borderRadius: '6px', background: bcsClass === cls ? c.accentSoft : '#fff',
                  color: bcsClass === cls ? c.accent : c.steel, fontWeight: bcsClass === cls ? 700 : 500,
                  cursor: 'pointer', fontSize: '0.85rem',
                }}>
                  {cls}
                </button>
              ))}
            </div>
          </div>

          {/* Dosage Form */}
          <div>
            <div style={lbl}>Dosage Form</div>
            <select
              value={dosageForm}
              onChange={e => setDosageForm(e.target.value)}
              style={{ marginTop: '0.5rem', width: '100%', padding: '0.5rem 0.75rem', border: `1px solid ${c.line}`, borderRadius: '6px', fontSize: '0.85rem', color: c.ink, background: '#fff' }}
            >
              {DOSAGE_FORMS.map(df => <option key={df.value} value={df.value}>{df.label}</option>)}
            </select>
          </div>
        </div>

        {/* Markets */}
        <div style={{ marginTop: '1rem' }}>
          <div style={lbl}>Target Markets</div>
          <div style={{ display: 'flex', gap: '0.4rem', flexWrap: 'wrap', marginTop: '0.5rem' }}>
            {MARKETS.map(m => (
              <button key={m} onClick={() => toggleMarket(m)} style={{
                padding: '0.35rem 0.75rem', border: `1px solid ${targetMarkets.includes(m) ? c.accent : c.line}`,
                borderRadius: '20px', background: targetMarkets.includes(m) ? c.accentSoft : '#fff',
                color: targetMarkets.includes(m) ? c.accent : c.steel,
                fontWeight: targetMarkets.includes(m) ? 700 : 500,
                cursor: 'pointer', fontSize: '0.78rem',
              }}>
                {m}
              </button>
            ))}
          </div>
        </div>

        {/* Toggles */}
        <div style={{ display: 'flex', gap: '1.5rem', marginTop: '1rem' }}>
          {[
            { label: 'Generic (ANDA)', val: isGeneric, set: setIsGeneric },
            { label: 'Pediatric', val: isPediatric, set: setIsPediatric },
          ].map(({ label, val, set }) => (
            <label key={label} style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', cursor: 'pointer', fontSize: '0.85rem', color: c.steel }}>
              <input type="checkbox" checked={val} onChange={e => set(e.target.checked)} style={{ accentColor: c.accent }} />
              {label}
            </label>
          ))}
        </div>

        <button
          onClick={runAnalysis}
          disabled={loading || !hasBlend}
          style={{
            marginTop: '1rem', background: loading || !hasBlend ? '#9eafa8' : c.accent,
            color: '#fff', border: 'none', borderRadius: '6px', padding: '0.65rem 1.5rem',
            fontWeight: 700, cursor: loading || !hasBlend ? 'not-allowed' : 'pointer', fontSize: '0.9rem',
          }}
        >
          {loading ? 'Analysing...' : 'Run Deep Dive →'}
        </button>
        {error && <div style={{ color: c.red, fontSize: '0.8rem', marginTop: '0.5rem' }}>{error}</div>}
      </div>

      {result && <PharmaResult result={result} />}
    </div>
  )
}

function PharmaResult({ result }) {
  const bcs = result.bcs_profile || {}
  const df = result.dosage_form_profile || {}
  const stability = result.stability_profile || {}
  const pathway = result.pathway_profile || {}

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>

      {/* ICM Score */}
      <div style={{ display: 'flex', gap: '0.65rem', flexWrap: 'wrap' }}>
        {[
          { lbl: 'BCS Class', val: result.bcs_class, sub: `${bcs.solubility} Sol · ${bcs.permeability} Perm` },
          { lbl: 'Dosage Form', val: result.recommended_dosage_form, sub: '' },
          { lbl: 'Manufacturing', val: result.manufacturing_route, sub: '' },
          { lbl: 'ICM Score', val: `${result.icm_score?.toFixed(0)}/100`, sub: 'Complexity metric' },
        ].map(({ lbl: l, val, sub }) => (
          <div key={l} style={{ ...panel, padding: '0.85rem 1rem', flex: '1 1 150px', boxShadow: 'none' }}>
            <div style={lbl}>{l}</div>
            <div style={{ color: c.ink, fontWeight: 780, fontSize: '1rem', marginTop: '0.3rem' }}>{val}</div>
            {sub && <div style={{ color: c.muted, fontSize: '0.72rem', marginTop: '0.1rem' }}>{sub}</div>}
          </div>
        ))}
      </div>

      {/* BCS Profile */}
      <div style={{ ...panel, padding: '1.1rem', boxShadow: 'none' }}>
        <div style={lbl}>BCS Class {result.bcs_class} — Formulation Strategy</div>
        <div style={{ color: c.muted, fontSize: '0.8rem', marginTop: '0.35rem', marginBottom: '0.5rem' }}>{bcs.absorption_prediction}</div>
        <div style={{ color: '#3d5060', fontSize: '0.85rem', lineHeight: 1.7 }}>{bcs.formulation_strategy}</div>
        <div style={{ marginTop: '0.75rem' }}>
          <div style={lbl}>Enabling Excipients</div>
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.35rem', marginTop: '0.4rem' }}>
            {(bcs.enabling_excipients || []).map(ex => (
              <span key={ex} style={{ background: c.accentSoft, color: c.accent, border: `1px solid ${c.line}`, borderRadius: '4px', padding: '2px 8px', fontSize: '0.75rem' }}>{ex}</span>
            ))}
          </div>
        </div>
        <div style={{ display: 'flex', gap: '1rem', marginTop: '0.75rem', flexWrap: 'wrap' }}>
          <div style={{ fontSize: '0.8rem', color: c.muted }}>Bioavailability risk: <strong style={{ color: c.ink }}>{bcs.bioavailability_risk}</strong></div>
          <div style={{ fontSize: '0.8rem', color: c.muted }}>IVIVC: <strong style={{ color: c.ink }}>{bcs.ivivc_potential}</strong></div>
        </div>
      </div>

      {/* Compatibility */}
      <div style={{ ...panel, padding: '1.1rem', boxShadow: 'none' }}>
        <div style={lbl}>API-Excipient Compatibility</div>
        <div style={{ color: '#3d5060', fontSize: '0.85rem', marginTop: '0.4rem', marginBottom: result.compatibility_results?.length ? '0.75rem' : 0 }}>
          {result.compatibility_summary}
        </div>
        {(result.compatibility_results || []).map((r, i) => (
          <div key={i} style={{ border: `1px solid ${c.line}`, borderRadius: '6px', padding: '0.75rem 0.9rem', marginBottom: '0.5rem', background: '#fafbfc' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', gap: '0.5rem', flexWrap: 'wrap', marginBottom: '0.35rem' }}>
              <strong style={{ fontSize: '0.85rem', color: c.ink }}>{r.ingredient}</strong>
              <div style={{ display: 'flex', gap: '0.4rem' }}>
                <Badge text={r.severity} color={SEV_COLOR[r.severity] || c.steel} />
                <Badge text={r.interaction_type} color={c.steel} />
              </div>
            </div>
            <div style={{ color: '#3d5060', fontSize: '0.8rem', lineHeight: 1.6, marginBottom: '0.3rem' }}>{r.mechanism}</div>
            <div style={{ color: c.accent, fontSize: '0.78rem' }}>Mitigation: {r.mitigation}</div>
            <div style={{ color: c.muted, fontSize: '0.72rem', marginTop: '0.2rem' }}>Ref: {r.literature_ref}</div>
          </div>
        ))}
      </div>

      {/* ICH Stability */}
      <div style={{ ...panel, padding: '1.1rem', boxShadow: 'none' }}>
        <div style={lbl}>ICH Q1A Stability — Zone {result.stability_zone}</div>
        {stability.zone && <div style={{ color: c.ink, fontWeight: 700, fontSize: '0.9rem', marginTop: '0.35rem' }}>{stability.zone}</div>}
        <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.5rem', marginTop: '0.5rem' }}>
          {[
            { k: 'Long-term', v: stability.long_term },
            { k: 'Intermediate', v: stability.intermediate },
            { k: 'Accelerated', v: stability.accelerated },
          ].filter(x => x.v && x.v !== 'Not required').map(({ k, v }) => (
            <div key={k} style={{ flex: '1 1 180px', background: c.accentSoft, border: `1px solid ${c.line}`, borderRadius: '6px', padding: '0.5rem 0.75rem' }}>
              <div style={lbl}>{k}</div>
              <div style={{ color: c.ink, fontSize: '0.82rem', marginTop: '0.2rem' }}>{v}</div>
            </div>
          ))}
        </div>
        {stability.regions?.length > 0 && (
          <div style={{ marginTop: '0.5rem', fontSize: '0.8rem', color: c.muted }}>Regions: {stability.regions.join(', ')}</div>
        )}
        {result.stability_concerns?.length > 0 && (
          <div style={{ marginTop: '0.65rem' }}>
            <div style={lbl}>Stability Concerns</div>
            {result.stability_concerns.map((s, i) => (
              <div key={i} style={{ color: c.amber, fontSize: '0.8rem', marginTop: '0.3rem' }}>⚠ {s}</div>
            ))}
          </div>
        )}
        <div style={{ marginTop: '0.6rem', fontSize: '0.82rem', color: '#3d5060' }}>
          <strong>Packaging:</strong> {result.packaging_recommendation}
        </div>
      </div>

      {/* Manufacturing */}
      <div style={{ ...panel, padding: '1.1rem', boxShadow: 'none' }}>
        <div style={lbl}>Manufacturing Route</div>
        <div style={{ color: c.accent, fontWeight: 700, fontSize: '1rem', marginTop: '0.35rem' }}>{result.manufacturing_route}</div>
        <div style={{ color: '#3d5060', fontSize: '0.85rem', marginTop: '0.3rem', lineHeight: 1.7 }}>{result.manufacturing_rationale}</div>
      </div>

      {/* Regulatory Pathway */}
      {pathway.pathway && (
        <div style={{ ...panel, padding: '1.1rem', boxShadow: 'none' }}>
          <div style={lbl}>Regulatory Pathway</div>
          <div style={{ color: c.ink, fontWeight: 700, fontSize: '0.95rem', marginTop: '0.35rem' }}>{pathway.pathway}</div>
          {pathway.description && <div style={{ color: '#3d5060', fontSize: '0.85rem', marginTop: '0.3rem', lineHeight: 1.7 }}>{pathway.description}</div>}
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.5rem', marginTop: '0.65rem' }}>
            {pathway.timeline && <div style={{ flex: '1 1 160px', fontSize: '0.8rem', color: c.muted }}>Timeline: <strong style={{ color: c.ink }}>{pathway.timeline}</strong></div>}
            {pathway.clinical_studies && <div style={{ flex: '1 1 200px', fontSize: '0.8rem', color: c.muted }}>Clinical: <strong style={{ color: c.ink }}>{pathway.clinical_studies}</strong></div>}
          </div>
          {pathway.key_requirements?.length > 0 && (
            <div style={{ marginTop: '0.65rem' }}>
              <div style={lbl}>Key Requirements</div>
              <ul style={{ margin: '0.35rem 0 0', paddingLeft: '1.2rem' }}>
                {pathway.key_requirements.map((req, i) => (
                  <li key={i} style={{ color: '#3d5060', fontSize: '0.82rem', marginBottom: '0.2rem' }}>{req}</li>
                ))}
              </ul>
            </div>
          )}
        </div>
      )}

      {/* Alternative Dosage Forms */}
      {result.alternative_forms?.length > 0 && (
        <div style={{ ...panel, padding: '1.1rem', boxShadow: 'none' }}>
          <div style={lbl}>Alternative Dosage Forms</div>
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.4rem', marginTop: '0.5rem' }}>
            {result.alternative_forms.map(f => (
              <span key={f} style={{ background: '#f1f5f8', color: c.steel, border: `1px solid ${c.line}`, borderRadius: '4px', padding: '3px 10px', fontSize: '0.8rem' }}>{f}</span>
            ))}
          </div>
        </div>
      )}

      {/* Risks & Recommendations */}
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.65rem' }}>
        <div style={{ ...panel, padding: '1.1rem', boxShadow: 'none' }}>
          <div style={lbl}>Development Risks</div>
          {result.development_risks?.length === 0 ? (
            <div style={{ color: c.accent, fontSize: '0.85rem', marginTop: '0.4rem' }}>No major risks identified</div>
          ) : (result.development_risks || []).map((r, i) => (
            <div key={i} style={{ color: c.red, fontSize: '0.82rem', marginTop: '0.4rem', lineHeight: 1.5 }}>⚠ {r}</div>
          ))}
        </div>
        <div style={{ ...panel, padding: '1.1rem', boxShadow: 'none' }}>
          <div style={lbl}>Recommendations</div>
          {(result.development_recommendations || []).map((r, i) => (
            <div key={i} style={{ color: '#3d5060', fontSize: '0.82rem', marginTop: '0.4rem', lineHeight: 1.5 }}>→ {r}</div>
          ))}
        </div>
      </div>

    </div>
  )
}
