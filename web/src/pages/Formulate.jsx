import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { PieChart, Pie, Cell, Tooltip as RechartTooltip, ResponsiveContainer } from 'recharts'
import { api, streamFormulate } from '../api/client'
import { PUBLIC_VERTICAL_GUIDES, VERTICAL_OPTIONS } from '../constants'
import { useAuth } from '../auth/AuthContext'
import { saveLastRun } from '../lib/session'

/* ── Design tokens (match luxury-upgrade.css) ──────────────── */
const ink    = '#111827'
const muted  = '#5f6f7d'
const accent = '#1f8a7a'
const accentSoft = '#eef6f3'
const steel  = '#4c6375'
const line   = '#dce6ee'
const amber  = '#b5893a'

const panel = {
  background: '#ffffff',
  border: `1px solid ${line}`,
  borderRadius: '8px',
  boxShadow: '0 12px 32px rgba(36,49,61,0.08)',
}

const label = {
  color: steel,
  fontSize: '0.7rem',
  letterSpacing: '0.08em',
  textTransform: 'uppercase',
  fontWeight: 800,
}

const input = {
  width: '100%',
  background: '#ffffff',
  border: `1px solid ${line}`,
  borderRadius: '6px',
  color: ink,
  padding: '0.9rem 1rem',
  fontSize: '0.92rem',
  boxSizing: 'border-box',
  outline: 'none',
  boxShadow: 'inset 0 1px 0 rgba(17,24,39,0.03)',
}

const DONUT_COLORS = [
  '#1f8a7a','#22BFA0','#6d8f4f','#b5893a','#4c6375',
  '#6DCFBE','#A8E6DC','#34A853','#2980B9','#0D5C4F',
]

/* ── Sub-components ────────────────────────────────────────── */
function MetricCard({ label: lbl, value, unit, tone = ink }) {
  return (
    <div style={{ ...panel, padding: '1rem', flex: '1 1 130px', minWidth: 'min(130px,100%)' }}>
      <div style={label}>{lbl}</div>
      <div style={{ color: tone, fontSize: '1.42rem', fontWeight: 840, marginTop: '0.4rem', overflowWrap: 'anywhere' }}>
        {value ?? '—'}
      </div>
      {unit && <div style={{ color: muted, fontSize: '0.72rem', marginTop: '0.15rem' }}>{unit}</div>}
    </div>
  )
}

function Section({ eyebrow, title, children, aside }) {
  return (
    <section style={{ ...panel, padding: '1.25rem', marginBottom: '1rem' }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', gap: '1rem', alignItems: 'flex-start', marginBottom: '1rem', flexWrap: 'wrap' }}>
        <div style={{ minWidth: 0 }}>
          <div style={label}>{eyebrow}</div>
          <h2 style={{ margin: '0.35rem 0 0', color: ink, fontSize: '1.1rem', fontWeight: 800 }}>{title}</h2>
        </div>
        {aside}
      </div>
      {children}
    </section>
  )
}

/* Live pipeline steps in the hero right column */
function HeroPipeline({ stepFlags, loading }) {
  const STEPS = [
    { key: 'parse',       done: 'parse_done',       label: 'Parse',    desc: 'Brief, vertical, batch size, and claim targets.' },
    { key: 'optimize',    done: 'optimize_done',     label: 'Optimize', desc: 'Pareto or Bayesian search against explicit constraints.' },
    { key: 'proof_stack', done: 'proof_stack_done',  label: 'Screen',   desc: 'Eco, regulatory, stability, carbon, certifications.' },
    { key: 'agents',      done: 'complete',          label: 'Explain',  desc: 'Readable readiness notes and result-level metrics.' },
  ]
  return (
    <div style={{ display: 'grid', gap: '8px' }}>
      {STEPS.map((s) => {
        const isDone = stepFlags[s.done] || stepFlags.complete
        const isActive = loading && stepFlags[s.key] && !isDone
        return (
          <div key={s.key} className={`if-live-step${isDone ? ' done' : isActive ? ' active' : ''}`}>
            <b>{isDone ? '✓ ' : isActive ? '· ' : ''}{s.label}</b>
            <span>{s.desc}</span>
          </div>
        )
      })}
    </div>
  )
}

/* Dark terminal log */
function LiveLog({ events }) {
  const ref = useRef(null)
  useEffect(() => {
    if (ref.current) ref.current.scrollTop = ref.current.scrollHeight
  }, [events])
  if (!events.length) return null
  return (
    <div ref={ref} className="if-terminal">
      {events.map((ev, i) => (
        <div key={i} style={{
          color: ev.step === 'error' ? '#F87171'
            : ev.step === 'complete' ? '#4ADE80'
            : ev.step?.includes('agent') ? '#60A5FA'
            : '#A3E635',
          marginBottom: '0.1rem',
          overflowWrap: 'anywhere',
        }}>
          <span style={{ opacity: 0.38 }}>{String(i + 1).padStart(2, '0')} </span>{ev.message}
        </div>
      ))}
    </div>
  )
}

/* Blend donut chart */
function BlendDonut({ blendEntries }) {
  if (!blendEntries.length) return null
  const data = blendEntries.map(([name, value]) => ({ name, value }))
  return (
    <div style={{ width: '100%', height: 210, marginBottom: '1.1rem' }}>
      <ResponsiveContainer>
        <PieChart>
          <Pie data={data} cx="50%" cy="50%" innerRadius={52} outerRadius={82} paddingAngle={2} dataKey="value">
            {data.map((_, i) => <Cell key={i} fill={DONUT_COLORS[i % DONUT_COLORS.length]} />)}
          </Pie>
          <RechartTooltip
            formatter={(val, name) => [`${val}%`, name]}
            contentStyle={{ fontSize: '0.78rem', borderRadius: '6px', border: `1px solid ${line}`, background: '#fff' }}
          />
        </PieChart>
      </ResponsiveContainer>
    </div>
  )
}

/* Pin & compare strip */
function PinnedCompare({ pinned, onUnpin }) {
  if (!pinned.length) return null
  const metrics = ['cost_per_kg', 'bio_pct', 'perf_score']
  const labls   = { cost_per_kg: 'Cost/kg', bio_pct: 'Bio %', perf_score: 'Perf' }
  const fmt     = { cost_per_kg: (v) => `$${Number(v).toFixed(2)}`, bio_pct: (v) => `${Number(v).toFixed(1)}%`, perf_score: (v) => Number(v).toFixed(1) }
  return (
    <div style={{ ...panel, padding: '1.1rem', marginBottom: '1rem', overflowX: 'auto' }}>
      <div style={label}>Pinned — compare runs</div>
      <div style={{ display: 'flex', gap: '1rem', marginTop: '0.8rem', minWidth: 'min-content' }}>
        {pinned.map((item, idx) => (
          <div key={idx} style={{ flex: '0 0 190px', border: `1px solid ${line}`, borderRadius: '8px', padding: '0.9rem', background: '#f9fafb' }}>
            <div style={{ color: accent, fontSize: '0.68rem', fontWeight: 800, marginBottom: '0.35rem', letterSpacing: '0.06em', textTransform: 'uppercase' }}>
              Run #{idx + 1}
            </div>
            <div style={{ color: ink, fontSize: '0.82rem', fontWeight: 750, marginBottom: '0.5rem', lineHeight: 1.4 }}>
              {item.vertical?.replace(/_/g, ' ')}
            </div>
            {metrics.map((m) => (
              <div key={m} style={{ display: 'flex', justifyContent: 'space-between', fontSize: '0.8rem', marginBottom: '0.25rem' }}>
                <span style={{ color: muted }}>{labls[m]}</span>
                <span style={{ color: ink, fontWeight: 700 }}>{item.result?.[m] != null ? fmt[m](item.result[m]) : '—'}</span>
              </div>
            ))}
            <button onClick={() => onUnpin(idx)} style={{ marginTop: '0.6rem', fontSize: '0.7rem', color: '#b45858', background: 'none', border: 'none', cursor: 'pointer', padding: 0, fontWeight: 650 }}>
              Remove
            </button>
          </div>
        ))}
      </div>
    </div>
  )
}

/* Conversational refinement chat */
function RefineChat({ currentResult, vertical, batchSize, optMode, onRefined }) {
  const [instruction, setInstruction] = useState('')
  const [loading, setLoading] = useState(false)
  const [err, setErr] = useState(null)

  const chips = ['Make it cheaper', 'More bio-based', 'Better performance', 'More sustainable', 'Premium quality']

  async function handleRefine(text) {
    const instr = text || instruction
    if (!instr.trim()) return
    setLoading(true)
    setErr(null)
    try {
      const res = await api.refine({ instruction: instr, current_result: currentResult, vertical, batch_size: batchSize, opt_mode: optMode })
      onRefined(res.data)
      setInstruction('')
    } catch (e) {
      setErr(e.response?.data?.detail || e.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={{ ...panel, padding: '1.1rem', marginBottom: '1rem', borderColor: '#b9d1ca', background: accentSoft }}>
      <div style={label}>Refinement assistant</div>
      <div style={{ color: muted, fontSize: '0.84rem', marginTop: '0.3rem', marginBottom: '0.85rem', lineHeight: 1.55 }}>
        Tell IntelliForm how to adjust this formulation — it re-runs the full pipeline with updated constraints.
      </div>
      <div style={{ display: 'flex', gap: '0.45rem', flexWrap: 'wrap', marginBottom: '0.75rem' }}>
        {chips.map((c) => (
          <button key={c} onClick={() => handleRefine(c)} disabled={loading}
            style={{ fontSize: '0.76rem', padding: '0.3rem 0.75rem', borderRadius: '999px', border: `1px solid ${accent}`, background: '#fff', color: accent, cursor: loading ? 'not-allowed' : 'pointer', fontWeight: 650 }}>
            {c}
          </button>
        ))}
      </div>
      <div style={{ display: 'flex', gap: '0.5rem' }}>
        <input value={instruction} onChange={(e) => setInstruction(e.target.value)} onKeyDown={(e) => e.key === 'Enter' && handleRefine()}
          placeholder="e.g. Reduce cost by 20%, increase bio-based content..."
          disabled={loading} style={{ ...input, flex: 1, padding: '0.7rem 0.9rem', background: '#fff' }} />
        <button onClick={() => handleRefine()} disabled={loading || !instruction.trim()}
          style={{ background: loading ? '#9eafa8' : accent, color: '#fff', border: 'none', borderRadius: '6px', padding: '0.7rem 1rem', fontWeight: 700, cursor: loading || !instruction.trim() ? 'not-allowed' : 'pointer', whiteSpace: 'nowrap', fontSize: '0.88rem' }}>
          {loading ? 'Refining...' : 'Refine'}
        </button>
      </div>
      {err && <div style={{ color: '#b45858', fontSize: '0.8rem', marginTop: '0.6rem' }}>{err}</div>}
    </div>
  )
}

/* ── Interaction Safety panel ──────────────────────────────── */
const SEV_COLOR  = { CRITICAL: '#b45858', HIGH: '#d97706', MEDIUM: '#b5893a', LOW: '#4c6375' }
const SEV_BG     = { CRITICAL: '#fff5f5', HIGH: '#fffbeb', MEDIUM: '#fffdf0', LOW: '#f9fafb' }
const SEV_BORDER = { CRITICAL: '#efb5b5', HIGH: '#fde68a', MEDIUM: '#e8c96b', LOW: '#dce6ee' }

function InteractionPanel({ data }) {
  if (!data) return null
  const { flags = [], critical_count = 0, high_count = 0, medium_count = 0, low_count = 0, safe, summary } = data
  return (
    <Section eyebrow="Interaction Safety" title="Ingredient incompatibility check">
      <div style={{ display: 'flex', gap: '0.65rem', flexWrap: 'wrap', marginBottom: '0.9rem', alignItems: 'stretch' }}>
        <div style={{ ...panel, padding: '0.85rem 1rem', flex: '0 0 auto', boxShadow: 'none', background: safe ? accentSoft : '#fff5f5', borderColor: safe ? '#b9d1ca' : '#efb5b5' }}>
          <div style={label}>Safety</div>
          <div style={{ color: safe ? accent : '#b45858', fontWeight: 780, marginTop: '0.3rem', fontSize: '0.95rem' }}>
            {safe ? '✓ No conflicts' : '⚠ Flags detected'}
          </div>
        </div>
        {critical_count > 0 && <div style={{ ...panel, padding: '0.85rem 1rem', flex: '0 0 auto', boxShadow: 'none' }}><div style={label}>Critical</div><div style={{ color: '#b45858', fontWeight: 780, fontSize: '1.2rem', marginTop: '0.3rem' }}>{critical_count}</div></div>}
        {high_count    > 0 && <div style={{ ...panel, padding: '0.85rem 1rem', flex: '0 0 auto', boxShadow: 'none' }}><div style={label}>High</div><div style={{ color: '#d97706', fontWeight: 780, fontSize: '1.2rem', marginTop: '0.3rem' }}>{high_count}</div></div>}
        {medium_count  > 0 && <div style={{ ...panel, padding: '0.85rem 1rem', flex: '0 0 auto', boxShadow: 'none' }}><div style={label}>Medium</div><div style={{ color: amber, fontWeight: 780, fontSize: '1.2rem', marginTop: '0.3rem' }}>{medium_count}</div></div>}
        {low_count     > 0 && <div style={{ ...panel, padding: '0.85rem 1rem', flex: '0 0 auto', boxShadow: 'none' }}><div style={label}>Low</div><div style={{ color: steel, fontWeight: 780, fontSize: '1.2rem', marginTop: '0.3rem' }}>{low_count}</div></div>}
      </div>
      <div style={{ color: safe ? accent : '#d97706', fontSize: '0.87rem', fontWeight: 650, marginBottom: flags.length ? '0.85rem' : 0 }}>
        {summary}
      </div>
      {flags.map((f, i) => (
        <div key={i} style={{ ...panel, padding: '0.9rem 1rem', marginBottom: '0.55rem', boxShadow: 'none', background: SEV_BG[f.severity], borderColor: SEV_BORDER[f.severity] }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '0.4rem', flexWrap: 'wrap', gap: '0.5rem' }}>
            <div style={{ color: ink, fontWeight: 750, fontSize: '0.88rem', overflowWrap: 'anywhere' }}>
              {f.ingredient_a} <span style={{ color: muted, fontWeight: 400 }}>×</span> {f.ingredient_b}
            </div>
            <span style={{ background: SEV_BG[f.severity], color: SEV_COLOR[f.severity], border: `1px solid ${SEV_BORDER[f.severity]}`, borderRadius: '999px', fontSize: '0.68rem', fontWeight: 800, padding: '0.2rem 0.6rem', letterSpacing: '0.06em', textTransform: 'uppercase', flex: '0 0 auto' }}>
              {f.severity}
            </span>
          </div>
          <div style={{ color: '#3d5060', fontSize: '0.83rem', lineHeight: 1.65, marginBottom: '0.4rem' }}>{f.mechanism}</div>
          <div style={{ color: SEV_COLOR[f.severity], fontSize: '0.8rem', lineHeight: 1.6 }}>
            <strong>Action:</strong> {f.recommendation}
          </div>
        </div>
      ))}
    </Section>
  )
}

/* ── CHEM21 Solvent Greenness panel ────────────────────────── */
const TIER_COLOR  = { 1: '#1f8a7a', 2: '#b5893a', 3: '#d97706', 4: '#b45858' }
const TIER_BG     = { 1: '#eef6f3', 2: '#fffdf0', 3: '#fffbeb', 4: '#fff5f5' }
const TIER_BORDER = { 1: '#b9d1ca', 2: '#e8c96b', 3: '#fde68a', 4: '#efb5b5' }
const TIER_NAMES  = { 1: 'Recommended', 2: 'Problematic', 3: 'Hazardous', 4: 'Highly Hazardous' }

function SolventPanel({ data }) {
  if (!data || !data.solvents_assessed?.length) return null
  const { weighted_score, grade, coverage_pct, solvents_assessed = [], alerts = [], substitution_suggestions = [] } = data
  const scoreColor = weighted_score >= 80 ? accent : weighted_score >= 50 ? amber : '#b45858'
  return (
    <Section eyebrow="Solvent Greenness" title="CHEM21 solvent assessment">
      <div style={{ display: 'flex', gap: '0.65rem', flexWrap: 'wrap', marginBottom: '1rem', alignItems: 'stretch' }}>
        <div style={{ ...panel, padding: '0.85rem 1rem', flex: '1 1 110px', boxShadow: 'none' }}>
          <div style={label}>Greenness Score</div>
          <div style={{ color: scoreColor, fontWeight: 780, fontSize: '1.4rem', marginTop: '0.3rem' }}>
            {weighted_score?.toFixed(0)}<span style={{ fontSize: '0.85rem', color: muted }}>/100</span>
          </div>
        </div>
        <div style={{ ...panel, padding: '0.85rem 1rem', flex: '0 0 auto', boxShadow: 'none' }}>
          <div style={label}>Grade</div>
          <div style={{ color: scoreColor, fontWeight: 780, fontSize: '1.4rem', marginTop: '0.3rem' }}>{grade}</div>
        </div>
        <div style={{ ...panel, padding: '0.85rem 1rem', flex: '1 1 100px', boxShadow: 'none' }}>
          <div style={label}>Coverage</div>
          <div style={{ color: ink, fontWeight: 780, fontSize: '1.1rem', marginTop: '0.3rem' }}>{coverage_pct?.toFixed(0)}%</div>
          <div style={{ color: muted, fontSize: '0.72rem' }}>of blend assessed</div>
        </div>
      </div>
      {alerts.length > 0 && (
        <div style={{ marginBottom: '0.85rem' }}>
          {alerts.map((a, i) => <div key={i} style={{ color: '#d97706', fontSize: '0.83rem', marginBottom: '0.25rem' }}>⚠ {a}</div>)}
        </div>
      )}
      <div style={{ display: 'grid', gap: '0.45rem' }}>
        {solvents_assessed.map((s, i) => (
          <div key={i} style={{ ...panel, padding: '0.75rem 1rem', boxShadow: 'none', background: TIER_BG[s.tier], borderColor: TIER_BORDER[s.tier] }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', flexWrap: 'wrap', gap: '0.5rem' }}>
              <div style={{ color: ink, fontWeight: 680, fontSize: '0.88rem', overflowWrap: 'anywhere' }}>{s.name}</div>
              <span style={{ color: TIER_COLOR[s.tier], background: TIER_BG[s.tier], border: `1px solid ${TIER_BORDER[s.tier]}`, borderRadius: '999px', fontSize: '0.68rem', fontWeight: 800, padding: '0.18rem 0.55rem', letterSpacing: '0.06em', textTransform: 'uppercase', flex: '0 0 auto' }}>
                T{s.tier} · {TIER_NAMES[s.tier]}
              </span>
            </div>
            <div style={{ color: muted, fontSize: '0.79rem', lineHeight: 1.55, marginTop: '0.3rem' }}>{s.rationale}</div>
          </div>
        ))}
      </div>
      {substitution_suggestions.length > 0 && (
        <div style={{ marginTop: '0.85rem', ...panel, padding: '0.9rem 1rem', boxShadow: 'none', background: '#fffdf0', borderColor: '#e8c96b' }}>
          <div style={{ color: amber, fontWeight: 780, fontSize: '0.83rem', marginBottom: '0.4rem' }}>Substitution suggestions</div>
          {substitution_suggestions.map((s, i) => <div key={i} style={{ color: '#3d5060', fontSize: '0.83rem', marginBottom: '0.2rem' }}>→ {s}</div>)}
        </div>
      )}
    </Section>
  )
}

/* ── CompTox / OPERA regulatory screening panel ────────────── */
function CompToxPanel({ data }) {
  if (!data) return null
  const {
    svhc_flags = [], cmr_flags = [], reach_restricted_flags = [],
    ready_biodeg_fraction, avg_log_bcf, avg_log_kow,
    regulatory_citation, coverage = 0,
  } = data

  const isStub = !coverage

  if (isStub) {
    return (
      <Section eyebrow="CompTox Screening" title="EPA OPERA QSAR regulatory screening">
        <div style={{ ...panel, padding: '1rem', background: '#f9fafb', borderColor: '#dce6ee', boxShadow: 'none' }}>
          <div style={{ color: muted, fontSize: '0.83rem', lineHeight: 1.65 }}>
            No CompTox matches returned for this blend. The EPA CCTE API was queried but found no records for these ingredient names.
          </div>
        </div>
      </Section>
    )
  }

  const hasConcerns = svhc_flags.length > 0 || cmr_flags.length > 0 || reach_restricted_flags.length > 0
  return (
    <Section eyebrow="CompTox Screening" title="EPA OPERA QSAR regulatory screening">
      <div style={{ display: 'flex', gap: '0.65rem', flexWrap: 'wrap', marginBottom: '1rem', alignItems: 'stretch' }}>
        <div style={{ ...panel, padding: '0.85rem 1rem', flex: '1 1 110px', boxShadow: 'none' }}>
          <div style={label}>Biodegradable</div>
          <div style={{ color: (ready_biodeg_fraction ?? 0) >= 70 ? accent : amber, fontWeight: 780, fontSize: '1.3rem', marginTop: '0.3rem' }}>{ready_biodeg_fraction?.toFixed(0)}%</div>
          <div style={{ color: muted, fontSize: '0.72rem' }}>OECD 301B (weighted)</div>
        </div>
        {avg_log_bcf != null && (
          <div style={{ ...panel, padding: '0.85rem 1rem', flex: '1 1 100px', boxShadow: 'none' }}>
            <div style={label}>Avg log BCF</div>
            <div style={{ color: avg_log_bcf > 3 ? '#d97706' : ink, fontWeight: 780, fontSize: '1.1rem', marginTop: '0.3rem' }}>{avg_log_bcf?.toFixed(2)}</div>
            <div style={{ color: muted, fontSize: '0.72rem' }}>log L/kg (bioaccum.)</div>
          </div>
        )}
        {avg_log_kow != null && (
          <div style={{ ...panel, padding: '0.85rem 1rem', flex: '1 1 100px', boxShadow: 'none' }}>
            <div style={label}>Avg log Kow</div>
            <div style={{ color: ink, fontWeight: 780, fontSize: '1.1rem', marginTop: '0.3rem' }}>{avg_log_kow?.toFixed(2)}</div>
            <div style={{ color: muted, fontSize: '0.72rem' }}>octanol-water</div>
          </div>
        )}
        <div style={{ ...panel, padding: '0.85rem 1rem', flex: '0 0 auto', boxShadow: 'none' }}>
          <div style={label}>Coverage</div>
          <div style={{ color: ink, fontWeight: 780, fontSize: '1.1rem', marginTop: '0.3rem' }}>{coverage}%</div>
        </div>
      </div>
      {svhc_flags.length > 0 && (
        <div style={{ ...panel, padding: '0.9rem 1rem', background: '#fff5f5', borderColor: '#efb5b5', boxShadow: 'none', marginBottom: '0.55rem' }}>
          <div style={{ color: '#b45858', fontWeight: 780, fontSize: '0.83rem', marginBottom: '0.35rem' }}>⛔ SVHC Candidates ({svhc_flags.length})</div>
          {svhc_flags.map((n, i) => <div key={i} style={{ color: '#3d5060', fontSize: '0.83rem', marginBottom: '0.2rem' }}>· {n}</div>)}
        </div>
      )}
      {cmr_flags.length > 0 && (
        <div style={{ ...panel, padding: '0.9rem 1rem', background: '#fff5f5', borderColor: '#efb5b5', boxShadow: 'none', marginBottom: '0.55rem' }}>
          <div style={{ color: '#b45858', fontWeight: 780, fontSize: '0.83rem', marginBottom: '0.35rem' }}>⚠ CMR Classified ({cmr_flags.length})</div>
          {cmr_flags.map((n, i) => <div key={i} style={{ color: '#3d5060', fontSize: '0.83rem', marginBottom: '0.2rem' }}>· {n}</div>)}
        </div>
      )}
      {reach_restricted_flags.length > 0 && (
        <div style={{ ...panel, padding: '0.9rem 1rem', background: '#fffbeb', borderColor: '#fde68a', boxShadow: 'none', marginBottom: '0.55rem' }}>
          <div style={{ color: '#d97706', fontWeight: 780, fontSize: '0.83rem', marginBottom: '0.35rem' }}>REACH Restricted ({reach_restricted_flags.length})</div>
          {reach_restricted_flags.map((n, i) => <div key={i} style={{ color: '#3d5060', fontSize: '0.83rem', marginBottom: '0.2rem' }}>· {n}</div>)}
        </div>
      )}
      {!hasConcerns && (
        <div style={{ color: accent, fontSize: '0.87rem' }}>✓ No SVHC, CMR, or REACH restrictions detected in this blend.</div>
      )}
      {regulatory_citation && (
        <div style={{ color: muted, fontSize: '0.75rem', marginTop: '0.75rem', lineHeight: 1.55 }}>{regulatory_citation}</div>
      )}
    </Section>
  )
}

const optModes = [
  { value: 'auto',   label: 'Auto' },
  { value: 'pareto', label: 'Pareto' },
  { value: 'bayesian', label: 'Bayesian' },
  { value: 'single', label: 'Single Objective' },
]

function formatVertical(value) {
  const match = VERTICAL_OPTIONS.find((item) => item.value === value)
  return match?.label ?? String(value || 'Unknown').replace(/_/g, ' ')
}

/* ── Main page ─────────────────────────────────────────────── */
export default function Formulate() {
  const auth = useAuth()
  const [inputText, setInputText]     = useState('')
  const [vertical, setVertical]       = useState('personal_care')
  const [batchSize, setBatchSize]     = useState(1000)
  const [optMode, setOptMode]         = useState('auto')
  const [loading, setLoading]         = useState(false)
  const [result, setResult]           = useState(null)
  const [error, setError]             = useState(null)
  const [streamEvents, setStreamEvents] = useState([])
  const [stepFlags, setStepFlags]     = useState({})
  const [agentStream, setAgentStream] = useState([])
  const [pinned, setPinned]           = useState([])
  const [isMobile, setIsMobile]       = useState(() => typeof window !== 'undefined' && window.innerWidth < 960)
  const abortRef = useRef(null)

  useEffect(() => {
    const handler = () => setIsMobile(window.innerWidth < 960)
    window.addEventListener('resize', handler)
    return () => window.removeEventListener('resize', handler)
  }, [])

  const canSubmit  = inputText.trim().length > 12
  const metrics    = result?.result
  const parsed     = result?.parsed
  const meta       = result?.meta
  const blendEntries = useMemo(() => Object.entries(metrics?.blend || {}).sort((a, b) => b[1] - a[1]), [metrics])
  const verticalGuide = PUBLIC_VERTICAL_GUIDES[vertical] || { status: 'beta', label: 'Beta', message: 'No validated public starter prompts configured for this vertical yet.', prompts: [] }

  const handleSubmit = useCallback(() => {
    if (auth?.supabaseEnabled && !auth?.user) { auth.signInWithGoogle(); return }
    if (abortRef.current) abortRef.current()
    setLoading(true); setError(null); setResult(null)
    setStreamEvents([]); setStepFlags({}); setAgentStream([])

    const abort = streamFormulate(
      { input_text: inputText, vertical, batch_size: batchSize, opt_mode: optMode, constraints: {} },
      (event) => {
        const { step, message, result: final } = event
        if (step === 'error')    { setError(message); setLoading(false); return }
        if (step === 'complete' && final) {
          setResult(final)
          saveLastRun({ request: { inputText, vertical, batchSize, optMode }, response: final })
          setStepFlags((f) => ({ ...f, complete: true }))
          setLoading(false); return
        }
        if (step === 'agent_comment') setAgentStream((p) => [...p, event.text ?? message])
        setStepFlags((f) => ({ ...f, [step]: true }))
        if (message) setStreamEvents((p) => [...p, event])
      }
    )
    abortRef.current = abort
  }, [inputText, vertical, batchSize, optMode, auth])

  const displayedAgents = result?.agents?.length ? result.agents : agentStream

  return (
    <div style={{ maxWidth: '1180px', margin: '0 auto', width: '100%', minWidth: 0 }}>

      {/* ── Luxury hero ─────────────────────────────────────── */}
      <div className="if-hero">
        <div>
          <div className="if-kicker">ChemeNova · IntelliForm Public</div>
          <h1 className="if-title">Formulation intelligence,<br />brief to proof stack.</h1>
          <p className="if-copy">
            A public workspace for translating product intent into formulation structure,
            sustainability posture, regulatory readiness, and optimization output that can be reviewed.
          </p>
          <div className="if-capability-strip">
            <div className="if-capability">
              <b>Vertical-aware</b>
              <span>Personal care, industrial, agricultural, pharmaceutical, food, laundry, and coatings.</span>
            </div>
            <div className="if-capability">
              <b>Optimization-ready</b>
              <span>Cost, bio-based share, performance, and certification as measurable constraints.</span>
            </div>
            <div className="if-capability">
              <b>Reviewable</b>
              <span>Results backed by parsed inputs, run metrics, risk notes, and real calculations.</span>
            </div>
          </div>
        </div>

        <aside className="if-pipeline-panel">
          <div className="if-pipeline-panel-title">
            Run architecture
            <span className="if-status-pill">Streaming</span>
          </div>
          <HeroPipeline stepFlags={stepFlags} loading={loading} />
        </aside>
      </div>

      {/* ── Public notice ────────────────────────────────────── */}
      <div style={{ ...panel, padding: '0.8rem 1rem', marginBottom: '1rem', background: '#fafcfb', borderColor: '#c7d4dd', boxShadow: 'none' }}>
        <span style={{ color: muted, fontSize: '0.85rem', lineHeight: 1.6 }}>
          Public edition: free-tier usage limits are enabled to keep IntelliForm broadly accessible. For internal, higher-volume, or enterprise use, keep the Streamlit lab app and the enterprise stack separate.
        </span>
      </div>

      {loading && <LiveLog events={streamEvents} />}

      <PinnedCompare pinned={pinned} onUnpin={(idx) => setPinned((p) => p.filter((_, i) => i !== idx))} />

      {/* ── Two-column work area ──────────────────────────────── */}
      <div style={{ display: 'grid', gridTemplateColumns: isMobile ? '1fr' : 'minmax(310px, 400px) minmax(0, 1fr)', gap: '1rem' }}>

        {/* Left: brief + starters */}
        <div>
          <Section eyebrow="Brief Builder" title="Describe the product target">
            {auth?.supabaseEnabled && !auth?.user ? (
              <div style={{ ...panel, padding: '1rem', background: '#fffdf0', border: `1px solid #e8c96b`, marginBottom: '1rem', boxShadow: 'none' }}>
                <div style={{ color: ink, fontWeight: 780, marginBottom: '0.3rem', fontSize: '0.92rem' }}>Sign in to generate</div>
                <div style={{ color: muted, fontSize: '0.84rem', lineHeight: 1.6, marginBottom: '0.75rem' }}>
                  Live formulation generation is tied to a signed-in free account so usage, history, and future upgrades belong to you.
                </div>
                <button onClick={() => auth.signInWithGoogle()}
                  style={{ background: accent, color: '#fff', border: 'none', borderRadius: '6px', padding: '0.65rem 1rem', fontWeight: 700, cursor: 'pointer' }}>
                  Continue with Google
                </button>
              </div>
            ) : null}

            <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.5rem', ...label }}>
              Natural-language formulation brief
            </label>
            <textarea
              value={inputText}
              onChange={(e) => setInputText(e.target.value)}
              placeholder="e.g. Low-VOC industrial degreaser for heavy equipment with high flash point, strong grease lift, and moderate foam."
              rows={7}
              style={{ ...input, resize: 'vertical', minHeight: '180px', lineHeight: 1.7 }}
            />

            <div style={{ display: 'grid', gap: '0.85rem', gridTemplateColumns: isMobile ? '1fr' : '1fr 1fr', marginTop: '1rem' }}>
              <div>
                <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.4rem', ...label }}>Vertical</label>
                <select value={vertical} onChange={(e) => setVertical(e.target.value)} style={input}>
                  {VERTICAL_OPTIONS.map((o) => <option key={o.value} value={o.value}>{o.label}</option>)}
                </select>
              </div>
              <div>
                <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.4rem', ...label }}>Optimization mode</label>
                <select value={optMode} onChange={(e) => setOptMode(e.target.value)} style={input}>
                  {optModes.map((o) => <option key={o.value} value={o.value}>{o.label}</option>)}
                </select>
              </div>
              <div style={{ gridColumn: '1 / -1' }}>
                <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.4rem', ...label }}>Batch size (kg)</label>
                <input type="number" value={batchSize} onChange={(e) => setBatchSize(Number(e.target.value))} style={input} />
              </div>
            </div>

            <button
              onClick={handleSubmit}
              disabled={loading || !canSubmit}
              style={{
                marginTop: '1rem',
                background: loading ? '#8fa8a0' : accent,
                color: '#fff',
                border: 'none',
                borderRadius: '6px',
                padding: '0.9rem 1.2rem',
                fontSize: '0.95rem',
                fontWeight: 750,
                cursor: loading || !canSubmit ? 'not-allowed' : 'pointer',
                width: '100%',
                boxShadow: loading ? 'none' : '0 8px 22px rgba(31,138,122,0.22)',
                letterSpacing: '0',
              }}
            >
              {loading ? 'Running IntelliForm...' : 'Run IntelliForm'}
            </button>
          </Section>

          <Section eyebrow="Starter Prompts" title={`${formatVertical(vertical)} demo prompts`}>
            <div style={{
              ...panel,
              padding: '0.9rem 1rem',
              marginBottom: '0.9rem',
              boxShadow: 'none',
              background: verticalGuide.status === 'validated' ? accentSoft : '#fffdf0',
              border: `1px solid ${verticalGuide.status === 'validated' ? '#b9d1ca' : '#e8c96b'}`,
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '0.4rem' }}>
                <div style={{ color: ink, fontWeight: 780, fontSize: '0.9rem' }}>{formatVertical(vertical)}</div>
                <span style={{ color: verticalGuide.status === 'validated' ? accent : amber, fontSize: '0.68rem', textTransform: 'uppercase', letterSpacing: '0.07em', fontWeight: 800 }}>
                  {verticalGuide.label}
                </span>
              </div>
              <div style={{ color: muted, fontSize: '0.83rem', lineHeight: 1.6 }}>{verticalGuide.message}</div>
            </div>

            {verticalGuide.prompts.length > 0 ? (
              <div style={{ display: 'grid', gap: '0.75rem' }}>
                {verticalGuide.prompts.map((ex) => (
                  <button key={ex.title} type="button" onClick={() => { setInputText(ex.text); setVertical(ex.vertical) }}
                    style={{ ...panel, padding: '0.9rem 1rem', textAlign: 'left', cursor: 'pointer', width: '100%', overflowWrap: 'anywhere', boxShadow: 'none' }}>
                    <div style={{ color: ink, fontWeight: 780, marginBottom: '0.2rem', fontSize: '0.9rem' }}>{ex.title}</div>
                    <div style={{ color: accent, fontSize: '0.68rem', textTransform: 'uppercase', letterSpacing: '0.07em', marginBottom: '0.35rem', fontWeight: 800 }}>{formatVertical(ex.vertical)}</div>
                    <div style={{ color: '#3d5060', fontSize: '0.83rem', lineHeight: 1.6, marginBottom: ex.note ? '0.4rem' : 0 }}>{ex.text}</div>
                    {ex.note && <div style={{ color: muted, fontSize: '0.74rem', lineHeight: 1.5 }}>{ex.note}</div>}
                  </button>
                ))}
              </div>
            ) : (
              <div style={{ color: muted, fontSize: '0.85rem', lineHeight: 1.6 }}>
                No validated public starter prompts for this vertical yet. For a first successful demo, switch to Agricultural, Food &amp; Beverage, or Fabric &amp; Laundry.
              </div>
            )}
          </Section>
        </div>

        {/* Right: results */}
        <div>
          {error && (
            <div style={{ ...panel, padding: '1rem 1.1rem', color: '#b45858', border: '1px solid #efb5b5', background: '#fff5f5', marginBottom: '1rem', boxShadow: 'none' }}>
              {error}
            </div>
          )}

          {/* Metric strip */}
          <div style={{ display: 'flex', gap: '0.65rem', flexWrap: 'wrap', marginBottom: '1rem', alignItems: 'stretch' }}>
            <MetricCard label="Cost / kg"       value={metrics?.cost_per_kg  ? `$${metrics.cost_per_kg.toFixed(2)}`  : '—'} tone={accent} />
            <MetricCard label="Bio-based"        value={metrics?.bio_pct      ? `${metrics.bio_pct.toFixed(1)}%`      : '—'} tone={amber} />
            <MetricCard label="Perf Score"       value={metrics?.perf_score   ? metrics.perf_score.toFixed(1)         : '—'} unit="/ 100" />
            <MetricCard label="Eco Grade"        value={result?.eco?.grade} tone={accent} />
            <MetricCard label="Ingredient Pool"  value={meta?.ingredient_pool_size ?? '—'} />
          </div>

          {/* Program lens */}
          <Section
            eyebrow="Program Lens"
            title="Parsed brief and optimization posture"
            aside={meta ? (
              <span style={{ padding: '0.3rem 0.75rem', borderRadius: '999px', background: accentSoft, border: '1px solid #b9d1ca', color: accent, fontSize: '0.72rem', fontWeight: 800 }}>
                {parsed?.parser_backend || 'Pending'}
              </span>
            ) : null}
          >
            {parsed ? (
              <div style={{ display: 'grid', gap: '0.85rem', gridTemplateColumns: isMobile ? '1fr' : '1fr 1fr' }}>
                <div style={{ ...panel, padding: '0.95rem 1rem', boxShadow: 'none' }}>
                  <div style={label}>Resolved Vertical</div>
                  <div style={{ color: ink, marginTop: '0.4rem', fontWeight: 780 }}>{formatVertical(meta?.resolved_vertical)}</div>
                  <div style={{ color: muted, marginTop: '0.25rem', fontSize: '0.78rem' }}>
                    Requested: {formatVertical(meta?.requested_vertical)} · Inferred: {formatVertical(meta?.inferred_vertical)}
                  </div>
                </div>
                <div style={{ ...panel, padding: '0.95rem 1rem', boxShadow: 'none' }}>
                  <div style={label}>Constraints Used</div>
                  <div style={{ color: ink, marginTop: '0.4rem', fontWeight: 780 }}>
                    ${meta?.constraints_used?.max_cost}/kg · {meta?.constraints_used?.min_bio}% bio · {meta?.constraints_used?.min_perf} perf
                  </div>
                  <div style={{ color: muted, marginTop: '0.25rem', fontSize: '0.78rem' }}>Mode: {meta?.optimization_mode_requested}</div>
                </div>
                <div style={{ ...panel, padding: '0.95rem 1rem', gridColumn: '1 / -1', boxShadow: 'none' }}>
                  <div style={label}>Parser Reasoning</div>
                  <div style={{ color: '#3d5060', marginTop: '0.5rem', lineHeight: 1.75, fontSize: '0.88rem' }}>{parsed.reasoning}</div>
                </div>
              </div>
            ) : (
              <div style={{ color: muted, lineHeight: 1.7, fontSize: '0.88rem' }}>
                Submit a brief to see how IntelliForm interpreted the request, which vertical it resolved to, and the exact optimization thresholds used.
              </div>
            )}
          </Section>

          {/* Blend architecture */}
          <Section
            eyebrow="Blend Architecture"
            title="Optimized composition"
            aside={blendEntries.length ? (
              <button
                onClick={() => { if (!result || pinned.length >= 3) return; setPinned((p) => [...p, { ...result, vertical }]) }}
                disabled={pinned.length >= 3}
                style={{ fontSize: '0.76rem', padding: '0.3rem 0.85rem', borderRadius: '999px', border: `1px solid ${pinned.length >= 3 ? line : accent}`, background: pinned.length >= 3 ? '#f9fafb' : accentSoft, color: pinned.length >= 3 ? muted : accent, cursor: pinned.length >= 3 ? 'not-allowed' : 'pointer', fontWeight: 680 }}>
                Pin &amp; compare
              </button>
            ) : null}
          >
            {blendEntries.length ? (
              <>
                <BlendDonut blendEntries={blendEntries} />
                <div style={{ display: 'flex', flexDirection: 'column', gap: '0.7rem' }}>
                  {blendEntries.map(([ingredient, pct]) => (
                    <div key={ingredient}>
                      <div style={{ display: 'flex', justifyContent: 'space-between', gap: '0.75rem', marginBottom: '0.3rem', alignItems: 'flex-start' }}>
                        <div style={{ color: '#3d5060', fontSize: '0.86rem', minWidth: 0, overflowWrap: 'anywhere' }}>{ingredient}</div>
                        <div style={{ color: ink, fontWeight: 780, flex: '0 0 auto', fontSize: '0.9rem' }}>{pct}%</div>
                      </div>
                      <div className="if-bar-wrap">
                        <div className="if-bar-fill" style={{ width: `${Math.min(pct, 100)}%` }} />
                      </div>
                    </div>
                  ))}
                </div>
              </>
            ) : (
              <div style={{ color: muted, fontSize: '0.87rem', lineHeight: 1.6 }}>
                No blend yet. The result will appear here after a successful run.
              </div>
            )}
          </Section>

          {/* Refinement chat — appears after a successful run */}
          {result && (
            <RefineChat currentResult={result} vertical={vertical} batchSize={batchSize} optMode={optMode} onRefined={setResult} />
          )}

          {/* Risk & compliance */}
          <Section eyebrow="System Readout" title="Risk, compliance, and readiness">
            <div style={{ display: 'grid', gap: '0.85rem', gridTemplateColumns: isMobile ? '1fr' : '1fr 1fr' }}>
              <div style={{ ...panel, padding: '1rem', boxShadow: 'none' }}>
                <div style={label}>Regulatory Status</div>
                <div style={{ color: result?.vreg?.overall_status?.includes('✅') ? accent : amber, marginTop: '0.4rem', fontWeight: 780 }}>
                  {result?.vreg?.overall_status || 'Pending'}
                </div>
                <div style={{ color: muted, marginTop: '0.3rem', fontSize: '0.78rem', lineHeight: 1.6 }}>
                  {result?.vreg?.framework || 'Framework notes will appear after formulation.'}
                </div>
              </div>
              <div style={{ ...panel, padding: '1rem', boxShadow: 'none' }}>
                <div style={label}>Stability</div>
                <div style={{ color: ink, marginTop: '0.4rem', fontWeight: 780 }}>{result?.stability?.overall_rating || 'Pending'}</div>
                <div style={{ color: muted, marginTop: '0.3rem', fontSize: '0.78rem', lineHeight: 1.6 }}>
                  Shelf life: {result?.stability?.shelf_life_range || '—'} · pH {result?.stability?.ph_min ?? '—'}–{result?.stability?.ph_max ?? '—'}
                </div>
              </div>
            </div>
            {metrics?.warnings?.length ? (
              <div style={{ marginTop: '1rem' }}>
                {metrics.warnings.map((w) => <div key={w} style={{ color: amber, fontSize: '0.83rem', marginBottom: '0.3rem' }}>⚠ {w}</div>)}
              </div>
            ) : null}
            {metrics?.compliance_flags?.length ? (
              <div style={{ marginTop: '0.8rem' }}>
                {metrics.compliance_flags.map((f) => <div key={f} style={{ color: '#b45858', fontSize: '0.83rem', marginBottom: '0.3rem' }}>✕ {f}</div>)}
              </div>
            ) : null}
          </Section>

          {/* Interaction safety */}
          <InteractionPanel data={result?.interactions} />

          {/* CHEM21 solvent greenness */}
          <SolventPanel data={result?.chem21} />

          {/* CompTox OPERA screening */}
          <CompToxPanel data={result?.comptox} />

          {/* Agent commentary */}
          <Section eyebrow="Technical Review" title="Commercial and scientific assessment">
            {displayedAgents.length ? (
              <div style={{ display: 'grid', gap: '0.65rem' }}>
                {displayedAgents.map((agent, i) => (
                  <div key={i} style={{ ...panel, padding: '0.95rem 1rem', color: '#3d5060', lineHeight: 1.75, overflowWrap: 'anywhere', boxShadow: 'none', animation: result?.agents?.length ? 'none' : 'fadeSlide 0.4s ease' }}>
                    {agent}
                  </div>
                ))}
              </div>
            ) : (
              <div style={{ color: muted, fontSize: '0.87rem', lineHeight: 1.7 }}>
                A concise technical review will stream in during the formulation run.
              </div>
            )}
          </Section>
        </div>
      </div>
    </div>
  )
}
