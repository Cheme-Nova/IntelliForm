import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { PieChart, Pie, Cell, Tooltip as RechartTooltip, ResponsiveContainer } from 'recharts'
import { api, streamFormulate } from '../api/client'
import { PUBLIC_VERTICAL_GUIDES, VERTICAL_OPTIONS } from '../constants'
import { useAuth } from '../auth/AuthContext'
import { saveLastRun } from '../lib/session'

const shell = {
  panel: {
    background: '#FFFFFF',
    border: '1px solid #E4E7EC',
    borderRadius: '8px',
    boxShadow: '0 1px 2px rgba(16, 24, 40, 0.04)',
  },
  label: {
    color: '#667085',
    fontSize: '0.72rem',
    letterSpacing: '0.08em',
    textTransform: 'uppercase',
    fontWeight: 750,
  },
  input: {
    width: '100%',
    background: '#FFFFFF',
    border: '1px solid #D0D5DD',
    borderRadius: '8px',
    color: '#101828',
    padding: '0.95rem 1rem',
    fontSize: '0.92rem',
    boxSizing: 'border-box',
    outline: 'none',
  },
}

const ink = '#101828'
const muted = '#667085'
const accent = '#127C67'
const accentSoft = '#EAF7F2'
const line = '#E4E7EC'

const DONUT_COLORS = [
  '#127C67', '#1A9E85', '#22BFA0', '#6DCFBE', '#A8E6DC',
  '#34A853', '#68D391', '#0D5C4F', '#E67E22', '#2980B9',
]

function MetricCard({ label, value, unit, tone = ink }) {
  return (
    <div style={{
      ...shell.panel,
      padding: '0.95rem 1rem',
      minWidth: 'min(138px, 100%)',
      flex: '1 1 138px',
      boxSizing: 'border-box',
    }}>
      <div style={{ color: muted, fontSize: '0.68rem', letterSpacing: '0.06em', textTransform: 'uppercase', fontWeight: 750 }}>
        {label}
      </div>
      <div style={{ color: tone, fontSize: '1.45rem', fontWeight: 820, marginTop: '0.35rem', overflowWrap: 'anywhere' }}>
        {value ?? '—'}
      </div>
      {unit ? <div style={{ color: muted, fontSize: '0.72rem', marginTop: '0.15rem' }}>{unit}</div> : null}
    </div>
  )
}

function Section({ eyebrow, title, children, aside }) {
  return (
    <section style={{ ...shell.panel, padding: '1.25rem', marginBottom: '1rem' }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', gap: '1rem', alignItems: 'flex-start', marginBottom: '1rem', flexWrap: 'wrap' }}>
        <div style={{ minWidth: 0 }}>
          <div style={shell.label}>{eyebrow}</div>
          <h2 style={{ margin: '0.35rem 0 0', color: ink, fontSize: '1.15rem', textAlign: 'left', fontWeight: 780, letterSpacing: 0 }}>{title}</h2>
        </div>
        {aside}
      </div>
      {children}
    </section>
  )
}

// Live step-by-step progress tracker
function AgentRunway({ steps }) {
  const STEP_CONFIG = [
    { key: 'parse', label: 'Parse brief' },
    { key: 'optimize', label: 'Optimize blend' },
    { key: 'proof_stack', label: 'Generate proof stack' },
    { key: 'agents', label: 'Agent commentary' },
  ]
  return (
    <div style={{
      display: 'grid',
      gridTemplateColumns: 'repeat(auto-fit, minmax(128px, 1fr))',
      gap: '0.5rem',
      marginTop: '1rem',
    }}>
      {STEP_CONFIG.map((cfg) => {
        const done = steps[cfg.key + '_done'] || steps.complete
        const active = steps[cfg.key] && !done
        return (
          <div key={cfg.key} style={{
            border: `1px solid ${done ? '#B9E1D3' : active ? '#D0D5DD' : line}`,
            background: done ? accentSoft : active ? '#F9FAFB' : '#F9FAFB',
            borderRadius: '8px',
            padding: '0.7rem',
            minHeight: '66px',
            transition: 'background 0.3s, border-color 0.3s',
          }}>
            <div style={{
              color: done ? accent : active ? '#B54708' : muted,
              fontSize: '0.72rem',
              fontWeight: 800,
              marginBottom: '0.3rem',
              transition: 'color 0.3s',
            }}>
              {done ? 'Complete' : active ? 'Working...' : 'Queued'}
            </div>
            <div style={{ color: ink, fontSize: '0.82rem', fontWeight: 700 }}>{cfg.label}</div>
          </div>
        )
      })}
    </div>
  )
}

// Live streaming agent log
function LiveLog({ events }) {
  const ref = useRef(null)
  useEffect(() => {
    if (ref.current) ref.current.scrollTop = ref.current.scrollHeight
  }, [events])

  if (!events.length) return null
  return (
    <div ref={ref} style={{
      background: '#0D1117',
      borderRadius: '8px',
      padding: '0.9rem 1rem',
      maxHeight: '200px',
      overflowY: 'auto',
      fontFamily: '"JetBrains Mono", "Fira Code", monospace',
      fontSize: '0.78rem',
      lineHeight: 1.65,
      marginTop: '1rem',
      border: '1px solid #21262D',
    }}>
      {events.map((ev, i) => (
        <div key={i} style={{
          color: ev.step === 'error' ? '#F87171'
            : ev.step === 'complete' ? '#4ADE80'
            : ev.step?.includes('agent') ? '#60A5FA'
            : '#A3E635',
          marginBottom: '0.1rem',
          overflowWrap: 'anywhere',
        }}>
          <span style={{ opacity: 0.45 }}>{String(i + 1).padStart(2, '0')} </span>
          {ev.message}
        </div>
      ))}
    </div>
  )
}

// Blend donut chart
function BlendDonut({ blendEntries }) {
  if (!blendEntries.length) return null
  const data = blendEntries.map(([name, value]) => ({ name, value }))
  return (
    <div style={{ width: '100%', height: 220, marginBottom: '1rem' }}>
      <ResponsiveContainer>
        <PieChart>
          <Pie
            data={data}
            cx="50%"
            cy="50%"
            innerRadius={55}
            outerRadius={85}
            paddingAngle={2}
            dataKey="value"
          >
            {data.map((_, i) => (
              <Cell key={i} fill={DONUT_COLORS[i % DONUT_COLORS.length]} />
            ))}
          </Pie>
          <RechartTooltip
            formatter={(val, name) => [`${val}%`, name]}
            contentStyle={{ fontSize: '0.8rem', borderRadius: '6px', border: `1px solid ${line}` }}
          />
        </PieChart>
      </ResponsiveContainer>
    </div>
  )
}

// Pinned formulations comparison strip
function PinnedCompare({ pinned, onUnpin }) {
  if (!pinned.length) return null
  const metrics = ['cost_per_kg', 'bio_pct', 'perf_score']
  const labels = { cost_per_kg: 'Cost/kg', bio_pct: 'Bio %', perf_score: 'Perf' }
  const fmt = { cost_per_kg: (v) => `$${Number(v).toFixed(2)}`, bio_pct: (v) => `${Number(v).toFixed(1)}%`, perf_score: (v) => Number(v).toFixed(1) }
  return (
    <div style={{ ...shell.panel, padding: '1.1rem', marginBottom: '1rem', overflowX: 'auto' }}>
      <div style={shell.label}>Pinned formulations — compare</div>
      <div style={{ display: 'flex', gap: '1rem', marginTop: '0.8rem', minWidth: 'min-content' }}>
        {pinned.map((item, idx) => (
          <div key={idx} style={{ flex: '0 0 200px', border: `1px solid ${line}`, borderRadius: '8px', padding: '0.85rem', background: '#F9FAFB' }}>
            <div style={{ color: accent, fontSize: '0.7rem', fontWeight: 800, marginBottom: '0.4rem', letterSpacing: '0.06em', textTransform: 'uppercase' }}>
              Run #{idx + 1}
            </div>
            <div style={{ color: ink, fontSize: '0.82rem', fontWeight: 700, marginBottom: '0.5rem', lineHeight: 1.4, overflowWrap: 'anywhere' }}>
              {item.vertical?.replace(/_/g, ' ')}
            </div>
            {metrics.map((m) => (
              <div key={m} style={{ display: 'flex', justifyContent: 'space-between', fontSize: '0.8rem', marginBottom: '0.25rem' }}>
                <span style={{ color: muted }}>{labels[m]}</span>
                <span style={{ color: ink, fontWeight: 700 }}>{item.result?.[m] != null ? fmt[m](item.result[m]) : '—'}</span>
              </div>
            ))}
            <button
              onClick={() => onUnpin(idx)}
              style={{ marginTop: '0.6rem', fontSize: '0.72rem', color: '#B42318', background: 'none', border: 'none', cursor: 'pointer', padding: 0 }}
            >
              Remove
            </button>
          </div>
        ))}
      </div>
    </div>
  )
}

// Refinement chat panel
function RefineChat({ currentResult, vertical, batchSize, optMode, onRefined }) {
  const [instruction, setInstruction] = useState('')
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState(null)

  const suggestions = [
    'Make it cheaper',
    'More bio-based',
    'Better performance',
    'More sustainable',
    'Premium quality',
  ]

  async function handleRefine(text) {
    const instr = text || instruction
    if (!instr.trim()) return
    setLoading(true)
    setError(null)
    try {
      const res = await api.refine({
        instruction: instr,
        current_result: currentResult,
        vertical,
        batch_size: batchSize,
        opt_mode: optMode,
      })
      onRefined(res.data)
      setInstruction('')
    } catch (err) {
      setError(err.response?.data?.detail || err.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div style={{ ...shell.panel, padding: '1.1rem', marginBottom: '1rem', borderColor: '#B9E1D3', background: accentSoft }}>
      <div style={shell.label}>Refinement assistant</div>
      <div style={{ color: muted, fontSize: '0.84rem', marginTop: '0.3rem', marginBottom: '0.85rem', lineHeight: 1.55 }}>
        Tell IntelliForm how to adjust this formulation and it will re-run the full pipeline.
      </div>
      <div style={{ display: 'flex', gap: '0.5rem', flexWrap: 'wrap', marginBottom: '0.75rem' }}>
        {suggestions.map((s) => (
          <button
            key={s}
            onClick={() => handleRefine(s)}
            disabled={loading}
            style={{
              fontSize: '0.78rem',
              padding: '0.35rem 0.75rem',
              borderRadius: '999px',
              border: `1px solid ${accent}`,
              background: '#fff',
              color: accent,
              cursor: loading ? 'not-allowed' : 'pointer',
              fontWeight: 650,
            }}
          >
            {s}
          </button>
        ))}
      </div>
      <div style={{ display: 'flex', gap: '0.5rem' }}>
        <input
          value={instruction}
          onChange={(e) => setInstruction(e.target.value)}
          onKeyDown={(e) => e.key === 'Enter' && handleRefine()}
          placeholder="e.g. Reduce cost by 20%, increase bio-based content..."
          disabled={loading}
          style={{ ...shell.input, flex: 1, padding: '0.7rem 0.9rem', background: '#fff' }}
        />
        <button
          onClick={() => handleRefine()}
          disabled={loading || !instruction.trim()}
          style={{
            background: loading ? '#98A2B3' : accent,
            color: '#fff',
            border: 'none',
            borderRadius: '8px',
            padding: '0.7rem 1rem',
            fontWeight: 700,
            cursor: loading || !instruction.trim() ? 'not-allowed' : 'pointer',
            whiteSpace: 'nowrap',
            fontSize: '0.88rem',
          }}
        >
          {loading ? 'Refining...' : 'Refine'}
        </button>
      </div>
      {error ? (
        <div style={{ color: '#B42318', fontSize: '0.82rem', marginTop: '0.6rem' }}>{error}</div>
      ) : null}
    </div>
  )
}

const optModes = [
  { value: 'auto', label: 'Auto' },
  { value: 'pareto', label: 'Pareto' },
  { value: 'bayesian', label: 'Bayesian' },
  { value: 'single', label: 'Single Objective' },
]

function formatVertical(value) {
  const match = VERTICAL_OPTIONS.find((item) => item.value === value)
  return match?.label ?? String(value || 'Unknown').replace(/_/g, ' ')
}

export default function Formulate() {
  const auth = useAuth()
  const [inputText, setInputText] = useState('')
  const [vertical, setVertical] = useState('personal_care')
  const [batchSize, setBatchSize] = useState(1000)
  const [optMode, setOptMode] = useState('auto')
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)
  const [streamEvents, setStreamEvents] = useState([])
  const [stepFlags, setStepFlags] = useState({})
  const [agentStream, setAgentStream] = useState([])
  const [pinned, setPinned] = useState([])
  const [isMobile, setIsMobile] = useState(() => {
    if (typeof window === 'undefined') return false
    return window.innerWidth < 960
  })
  const abortRef = useRef(null)

  useEffect(() => {
    function onResize() { setIsMobile(window.innerWidth < 960) }
    window.addEventListener('resize', onResize)
    return () => window.removeEventListener('resize', onResize)
  }, [])

  const canSubmit = inputText.trim().length > 12
  const metrics = result?.result
  const parsed = result?.parsed
  const meta = result?.meta
  const blendEntries = useMemo(
    () => Object.entries(metrics?.blend || {}).sort((a, b) => b[1] - a[1]),
    [metrics],
  )
  const verticalGuide = PUBLIC_VERTICAL_GUIDES[vertical] || {
    status: 'beta', label: 'Beta',
    message: 'No public starter prompts configured for this vertical yet.',
    prompts: [],
  }

  const handleSubmit = useCallback(() => {
    if (auth?.supabaseEnabled && !auth?.user) {
      auth.signInWithGoogle()
      return
    }
    if (abortRef.current) abortRef.current()
    setLoading(true)
    setError(null)
    setResult(null)
    setStreamEvents([])
    setStepFlags({})
    setAgentStream([])

    const abort = streamFormulate(
      { input_text: inputText, vertical, batch_size: batchSize, opt_mode: optMode, constraints: {} },
      (event) => {
        const { step, message, result: finalResult } = event

        if (step === 'error') {
          setError(message)
          setLoading(false)
          return
        }

        if (step === 'complete' && finalResult) {
          setResult(finalResult)
          saveLastRun({ request: { inputText, vertical, batchSize, optMode }, response: finalResult })
          setStepFlags((f) => ({ ...f, complete: true }))
          setLoading(false)
          return
        }

        if (step === 'agent_comment') {
          setAgentStream((prev) => [...prev, event.text ?? message])
        }

        setStepFlags((f) => ({ ...f, [step]: true }))
        if (message) {
          setStreamEvents((prev) => [...prev, event])
        }
      }
    )
    abortRef.current = abort
  }, [inputText, vertical, batchSize, optMode, auth])

  function loadPrompt(example) {
    setInputText(example.text)
    setVertical(example.vertical)
  }

  function pinResult() {
    if (!result || pinned.length >= 3) return
    setPinned((p) => [...p, { ...result, vertical }])
  }

  function unpinResult(idx) {
    setPinned((p) => p.filter((_, i) => i !== idx))
  }

  const displayedAgents = result?.agents?.length ? result.agents : agentStream

  return (
    <div style={{ maxWidth: '1180px', margin: '0 auto', width: '100%', minWidth: 0 }}>
      <section style={{
        ...shell.panel,
        padding: isMobile ? '1.15rem' : '1.5rem',
        marginBottom: '1rem',
        background: '#FFFFFF',
      }}>
        <div style={{ display: 'flex', gap: '1.25rem', flexWrap: 'wrap', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div style={{ maxWidth: '680px', minWidth: 0, flex: '1 1 320px' }}>
            <div style={{ ...shell.label, color: accent }}>IntelliForm Free</div>
            <h1 style={{ margin: '0.45rem 0 0.6rem', color: ink, fontSize: isMobile ? '2rem' : '2.65rem', textAlign: 'left', lineHeight: 1.04, letterSpacing: 0, fontWeight: 860 }}>
              Formulation intelligence, from brief to proof stack.
            </h1>
            <p style={{ color: muted, textAlign: 'left', lineHeight: 1.65, margin: 0, fontSize: '0.98rem' }}>
              Describe the target product. IntelliForm parses the brief, chooses the vertical, optimizes the blend, and prepares sustainability, regulatory, carbon, stability, and certification screens from the same run.
            </p>
          </div>
          <div style={{ display: 'grid', gap: '0.75rem', minWidth: isMobile ? 0 : '280px', flex: '1 1 280px', width: isMobile ? '100%' : 'auto' }}>
            <div style={{ border: `1px solid ${line}`, borderRadius: '8px', padding: '1rem', background: '#F9FAFB' }}>
              <div style={{ color: muted, fontSize: '0.68rem', letterSpacing: '0.06em', textTransform: 'uppercase', fontWeight: 750 }}>Agentic run plan</div>
              <div style={{ color: ink, marginTop: '0.35rem', fontWeight: 760 }}>Parse → optimize → validate → report</div>
            </div>
            <div style={{ border: `1px solid #B9E1D3`, borderRadius: '8px', padding: '1rem', background: accentSoft }}>
              <div style={{ color: accent, fontSize: '0.68rem', letterSpacing: '0.06em', textTransform: 'uppercase', fontWeight: 750 }}>Streaming edition</div>
              <div style={{ color: ink, marginTop: '0.35rem', fontWeight: 760 }}>Live progress · Refinement chat · Compare runs</div>
            </div>
          </div>
        </div>
        <AgentRunway steps={stepFlags} />
        {loading && <LiveLog events={streamEvents} />}
      </section>

      <PinnedCompare pinned={pinned} onUnpin={unpinResult} />

      <div style={{ display: 'grid', gridTemplateColumns: isMobile ? '1fr' : 'minmax(320px, 420px) minmax(0, 1fr)', gap: '1rem' }}>
        <div>
          <Section eyebrow="Brief Builder" title="Describe the product target">
            {auth?.supabaseEnabled && !auth?.user ? (
              <div style={{
                ...shell.panel,
                padding: '0.95rem 1rem',
                background: '#FFF8EB',
                border: '1px solid #FEDF89',
                marginBottom: '1rem',
              }}>
                <div style={{ color: ink, fontWeight: 760, marginBottom: '0.35rem' }}>Sign in to generate</div>
                <div style={{ color: muted, fontSize: '0.84rem', lineHeight: 1.6, marginBottom: '0.75rem' }}>
                  The public edition stays browseable, but live formulation generation is tied to a signed-in free account so usage, history, and future upgrades belong to you.
                </div>
                <button
                  onClick={() => auth.signInWithGoogle()}
                  style={{ background: accent, color: '#fff', border: 'none', borderRadius: '8px', padding: '0.7rem 1rem', fontWeight: 700, cursor: 'pointer' }}
                >
                  Continue with Google
                </button>
              </div>
            ) : null}
            <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.55rem', color: muted, fontSize: '0.82rem', fontWeight: 650 }}>
              Natural-language formulation brief
            </label>
            <textarea
              value={inputText}
              onChange={(event) => setInputText(event.target.value)}
              placeholder="e.g. Low-VOC industrial degreaser for heavy equipment with high flash point, strong grease lift, and moderate foam."
              rows={7}
              style={{ ...shell.input, resize: 'vertical', minHeight: '180px', lineHeight: 1.7 }}
            />
            <div style={{ display: 'grid', gap: '0.85rem', gridTemplateColumns: isMobile ? '1fr' : 'repeat(2, minmax(0, 1fr))', marginTop: '1rem' }}>
              <div>
                <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.45rem', color: muted, fontSize: '0.76rem', fontWeight: 650 }}>Vertical</label>
                <select value={vertical} onChange={(event) => setVertical(event.target.value)} style={shell.input}>
                  {VERTICAL_OPTIONS.map((option) => (
                    <option key={option.value} value={option.value}>{option.label}</option>
                  ))}
                </select>
              </div>
              <div>
                <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.45rem', color: muted, fontSize: '0.76rem', fontWeight: 650 }}>Optimization mode</label>
                <select value={optMode} onChange={(event) => setOptMode(event.target.value)} style={shell.input}>
                  {optModes.map((option) => (
                    <option key={option.value} value={option.value}>{option.label}</option>
                  ))}
                </select>
              </div>
              <div style={{ gridColumn: '1 / -1' }}>
                <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.45rem', color: muted, fontSize: '0.76rem', fontWeight: 650 }}>Batch size (kg)</label>
                <input
                  type="number"
                  value={batchSize}
                  onChange={(event) => setBatchSize(Number(event.target.value))}
                  style={shell.input}
                />
              </div>
            </div>
            <button
              onClick={handleSubmit}
              disabled={loading || !canSubmit}
              style={{
                marginTop: '1rem',
                background: loading ? '#98A2B3' : accent,
                color: '#fff',
                border: 'none',
                borderRadius: '8px',
                padding: '0.9rem 1.2rem',
                fontSize: '0.95rem',
                fontWeight: 700,
                cursor: loading || !canSubmit ? 'not-allowed' : 'pointer',
                width: '100%',
                boxShadow: loading ? 'none' : '0 8px 18px rgba(18, 124, 103, 0.18)',
              }}
            >
              {loading ? 'Running IntelliForm...' : 'Run IntelliForm'}
            </button>
          </Section>

          <Section eyebrow="Starter Prompts" title={`${formatVertical(vertical)} public demo prompts`}>
            <div style={{
              ...shell.panel,
              padding: '0.95rem 1rem',
              background: verticalGuide.status === 'validated' ? accentSoft : '#FFF8EB',
              border: verticalGuide.status === 'validated' ? '1px solid #B9E1D3' : '1px solid #FEDF89',
              marginBottom: '0.9rem',
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', gap: '1rem', alignItems: 'center', marginBottom: '0.45rem' }}>
                <div style={{ color: ink, fontWeight: 760 }}>{formatVertical(vertical)}</div>
                <div style={{ color: verticalGuide.status === 'validated' ? accent : '#B54708', fontSize: '0.72rem', textTransform: 'uppercase', letterSpacing: '0.06em', fontWeight: 800 }}>
                  {verticalGuide.label}
                </div>
              </div>
              <div style={{ color: muted, fontSize: '0.84rem', lineHeight: 1.6 }}>{verticalGuide.message}</div>
            </div>
            {verticalGuide.prompts.length > 0 ? (
              <div style={{ display: 'grid', gap: '0.8rem' }}>
                {verticalGuide.prompts.map((example) => (
                  <button
                    key={example.title}
                    type="button"
                    onClick={() => loadPrompt(example)}
                    style={{ ...shell.panel, padding: '0.95rem 1rem', textAlign: 'left', cursor: 'pointer', background: '#FFFFFF', width: '100%', overflowWrap: 'anywhere' }}
                  >
                    <div style={{ color: ink, fontWeight: 760, marginBottom: '0.25rem' }}>{example.title}</div>
                    <div style={{ color: accent, fontSize: '0.72rem', textTransform: 'uppercase', letterSpacing: '0.06em', marginBottom: '0.4rem', fontWeight: 800 }}>
                      {formatVertical(example.vertical)}
                    </div>
                    <div style={{ color: '#344054', fontSize: '0.84rem', lineHeight: 1.6, marginBottom: '0.5rem' }}>{example.text}</div>
                    {example.note ? (
                      <div style={{ color: muted, fontSize: '0.76rem', lineHeight: 1.5 }}>{example.note}</div>
                    ) : null}
                  </button>
                ))}
              </div>
            ) : (
              <div style={{ color: muted, fontSize: '0.86rem', lineHeight: 1.6 }}>
                No validated public starter prompts for this vertical yet. For a first successful demo, switch to Agricultural, Food &amp; Beverage, or Fabric &amp; Laundry.
              </div>
            )}
          </Section>
        </div>

        <div>
          {error ? (
            <div style={{ ...shell.panel, padding: '1rem 1.1rem', color: '#B42318', border: '1px solid #FECDCA', background: '#FEF3F2', marginBottom: '1rem' }}>
              {error}
            </div>
          ) : null}

          <div style={{ display: 'flex', gap: '0.75rem', flexWrap: 'wrap', marginBottom: '1rem', alignItems: 'stretch' }}>
            <MetricCard label="Cost / kg" value={metrics?.cost_per_kg ? `$${metrics.cost_per_kg.toFixed(2)}` : '—'} tone={accent} />
            <MetricCard label="Bio-based" value={metrics?.bio_pct ? `${metrics.bio_pct.toFixed(1)}%` : '—'} tone="#B54708" />
            <MetricCard label="Perf Score" value={metrics?.perf_score ? metrics.perf_score.toFixed(1) : '—'} unit="/ 100" />
            <MetricCard label="Eco Grade" value={result?.eco?.grade} tone={accent} />
            <MetricCard label="Ingredient Pool" value={meta?.ingredient_pool_size ?? '—'} />
          </div>

          <Section
            eyebrow="Program Lens"
            title="Parsed brief and optimization posture"
            aside={meta ? (
              <div style={{ borderRadius: '999px', padding: '0.35rem 0.75rem', background: accentSoft, border: '1px solid #B9E1D3', color: accent, fontSize: '0.75rem', fontWeight: 800 }}>
                {parsed?.parser_backend || 'Parser pending'}
              </div>
            ) : null}
          >
            {parsed ? (
              <div style={{ display: 'grid', gap: '0.85rem', gridTemplateColumns: isMobile ? '1fr' : 'repeat(2, minmax(0, 1fr))' }}>
                <div style={{ ...shell.panel, padding: '0.95rem 1rem' }}>
                  <div style={shell.label}>Resolved Vertical</div>
                  <div style={{ color: ink, marginTop: '0.45rem', fontWeight: 760 }}>{formatVertical(meta?.resolved_vertical)}</div>
                  <div style={{ color: muted, marginTop: '0.3rem', fontSize: '0.8rem' }}>
                    Requested: {formatVertical(meta?.requested_vertical)} · Inferred: {formatVertical(meta?.inferred_vertical)}
                  </div>
                </div>
                <div style={{ ...shell.panel, padding: '0.95rem 1rem' }}>
                  <div style={shell.label}>Constraints Used</div>
                  <div style={{ color: ink, marginTop: '0.45rem', fontWeight: 760 }}>
                    ${meta?.constraints_used?.max_cost}/kg · {meta?.constraints_used?.min_bio}% bio · {meta?.constraints_used?.min_perf} perf
                  </div>
                  <div style={{ color: muted, marginTop: '0.3rem', fontSize: '0.8rem' }}>Mode: {meta?.optimization_mode_requested}</div>
                </div>
                <div style={{ ...shell.panel, padding: '0.95rem 1rem', gridColumn: '1 / -1' }}>
                  <div style={shell.label}>Parser Reasoning</div>
                  <div style={{ color: '#344054', marginTop: '0.55rem', lineHeight: 1.7, fontSize: '0.88rem' }}>
                    {parsed.reasoning}
                  </div>
                </div>
              </div>
            ) : (
              <div style={{ color: muted, textAlign: 'left', lineHeight: 1.7 }}>
                Submit a brief to see how IntelliForm interpreted the request, which vertical it resolved to, and the exact optimization thresholds used.
              </div>
            )}
          </Section>

          <Section
            eyebrow="Blend Architecture"
            title="Optimized composition"
            aside={blendEntries.length ? (
              <div style={{ display: 'flex', gap: '0.5rem' }}>
                <button
                  onClick={pinResult}
                  disabled={pinned.length >= 3}
                  title={pinned.length >= 3 ? 'Max 3 pinned' : 'Pin this formulation to compare'}
                  style={{
                    fontSize: '0.78rem',
                    padding: '0.35rem 0.85rem',
                    borderRadius: '999px',
                    border: `1px solid ${pinned.length >= 3 ? line : accent}`,
                    background: pinned.length >= 3 ? '#F9FAFB' : accentSoft,
                    color: pinned.length >= 3 ? muted : accent,
                    cursor: pinned.length >= 3 ? 'not-allowed' : 'pointer',
                    fontWeight: 650,
                  }}
                >
                  Pin & compare
                </button>
              </div>
            ) : null}
          >
            {blendEntries.length ? (
              <>
                <BlendDonut blendEntries={blendEntries} />
                <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem' }}>
                  {blendEntries.map(([ingredient, pct]) => (
                    <div key={ingredient}>
                      <div style={{ display: 'flex', justifyContent: 'space-between', gap: '0.75rem', marginBottom: '0.35rem', alignItems: 'flex-start' }}>
                        <div style={{ color: '#344054', fontSize: '0.88rem', textAlign: 'left', minWidth: 0, overflowWrap: 'anywhere' }}>{ingredient}</div>
                        <div style={{ color: ink, fontWeight: 760, flex: '0 0 auto' }}>{pct}%</div>
                      </div>
                      <div style={{ background: '#EAECF0', borderRadius: '999px', height: '10px', overflow: 'hidden' }}>
                        <div style={{ width: `${Math.min(pct, 100)}%`, height: '10px', background: accent, transition: 'width 0.6s ease' }} />
                      </div>
                    </div>
                  ))}
                </div>
              </>
            ) : (
              <div style={{ color: muted, textAlign: 'left' }}>
                No blend yet. The result will appear here after a successful run.
              </div>
            )}
          </Section>

          {result && (
            <RefineChat
              currentResult={result}
              vertical={vertical}
              batchSize={batchSize}
              optMode={optMode}
              onRefined={(refined) => setResult(refined)}
            />
          )}

          <Section eyebrow="System Readout" title="Risk, compliance, and readiness">
            <div style={{ display: 'grid', gap: '0.9rem', gridTemplateColumns: isMobile ? '1fr' : 'repeat(2, minmax(0, 1fr))' }}>
              <div style={{ ...shell.panel, padding: '1rem' }}>
                <div style={shell.label}>Regulatory Status</div>
                <div style={{ color: result?.vreg?.overall_status?.includes('✅') ? accent : '#B54708', marginTop: '0.45rem', fontWeight: 760 }}>
                  {result?.vreg?.overall_status || 'Pending'}
                </div>
                <div style={{ color: muted, marginTop: '0.35rem', fontSize: '0.8rem', lineHeight: 1.6 }}>
                  {result?.vreg?.framework || 'Framework notes will appear after formulation.'}
                </div>
              </div>
              <div style={{ ...shell.panel, padding: '1rem' }}>
                <div style={shell.label}>Stability</div>
                <div style={{ color: ink, marginTop: '0.45rem', fontWeight: 760 }}>{result?.stability?.overall_rating || 'Pending'}</div>
                <div style={{ color: muted, marginTop: '0.35rem', fontSize: '0.8rem', lineHeight: 1.6 }}>
                  Shelf life: {result?.stability?.shelf_life_range || '—'} · pH {result?.stability?.ph_min ?? '—'}–{result?.stability?.ph_max ?? '—'}
                </div>
              </div>
            </div>
            {metrics?.warnings?.length ? (
              <div style={{ marginTop: '1rem' }}>
                {metrics.warnings.map((w) => (
                  <div key={w} style={{ color: '#B54708', textAlign: 'left', fontSize: '0.84rem', marginBottom: '0.35rem' }}>⚠ {w}</div>
                ))}
              </div>
            ) : null}
            {metrics?.compliance_flags?.length ? (
              <div style={{ marginTop: '0.8rem' }}>
                {metrics.compliance_flags.map((f) => (
                  <div key={f} style={{ color: '#B42318', textAlign: 'left', fontSize: '0.84rem', marginBottom: '0.35rem' }}>✕ {f}</div>
                ))}
              </div>
            ) : null}
          </Section>

          <Section eyebrow="Agent Commentary" title="Commercial and technical review">
            {displayedAgents.length ? (
              <div style={{ display: 'grid', gap: '0.7rem' }}>
                {displayedAgents.map((agent, i) => (
                  <div
                    key={i}
                    style={{
                      ...shell.panel,
                      padding: '0.95rem 1rem',
                      color: '#344054',
                      textAlign: 'left',
                      lineHeight: 1.7,
                      overflowWrap: 'anywhere',
                      animation: result?.agents?.length ? 'none' : 'fadeSlide 0.4s ease',
                    }}
                  >
                    {agent}
                  </div>
                ))}
              </div>
            ) : (
              <div style={{ color: muted, textAlign: 'left' }}>
                Agent commentary will stream in during the formulation run.
              </div>
            )}
          </Section>
        </div>
      </div>

      <style>{`
        @keyframes fadeSlide {
          from { opacity: 0; transform: translateY(6px); }
          to   { opacity: 1; transform: translateY(0); }
        }
      `}</style>
    </div>
  )
}
