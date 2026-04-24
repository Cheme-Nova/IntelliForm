import { useMemo, useState } from 'react'
import { api } from '../api/client'
import { PUBLIC_VERTICAL_GUIDES, VERTICAL_OPTIONS } from '../constants'
import { useAuth } from '../auth/AuthContext'

const shell = {
  panel: {
    background: 'linear-gradient(180deg, rgba(13,31,60,0.94), rgba(7,19,37,0.94))',
    border: '1px solid rgba(71, 104, 146, 0.45)',
    borderRadius: '22px',
    boxShadow: '0 24px 60px rgba(2, 6, 23, 0.35)',
  },
  label: {
    color: '#7dd3c8',
    fontSize: '0.72rem',
    letterSpacing: '0.18em',
    textTransform: 'uppercase',
    fontWeight: 700,
  },
  input: {
    width: '100%',
    background: 'rgba(5, 15, 31, 0.82)',
    border: '1px solid rgba(71, 104, 146, 0.55)',
    borderRadius: '16px',
    color: '#f8fafc',
    padding: '0.95rem 1rem',
    fontSize: '0.92rem',
    boxSizing: 'border-box',
    outline: 'none',
  },
}

function MetricCard({ label, value, unit, tone = '#f8fafc' }) {
  return (
    <div style={{
      ...shell.panel,
      padding: '1rem 1.1rem',
      minWidth: '138px',
      flex: 1,
    }}>
      <div style={{ color: '#7c8aa5', fontSize: '0.68rem', letterSpacing: '0.14em', textTransform: 'uppercase' }}>
        {label}
      </div>
      <div style={{ color: tone, fontSize: '1.5rem', fontWeight: 800, marginTop: '0.35rem' }}>
        {value ?? '—'}
      </div>
      {unit ? <div style={{ color: '#64748b', fontSize: '0.72rem', marginTop: '0.15rem' }}>{unit}</div> : null}
    </div>
  )
}

function Section({ eyebrow, title, children, aside }) {
  return (
    <section style={{ ...shell.panel, padding: '1.4rem', marginBottom: '1rem' }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', gap: '1rem', alignItems: 'flex-start', marginBottom: '1rem' }}>
        <div>
          <div style={shell.label}>{eyebrow}</div>
          <h2 style={{ margin: '0.4rem 0 0', color: '#f8fafc', fontSize: '1.2rem', textAlign: 'left' }}>{title}</h2>
        </div>
        {aside}
      </div>
      {children}
    </section>
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

  const canSubmit = inputText.trim().length > 12
  const metrics = result?.result
  const parsed = result?.parsed
  const meta = result?.meta
  const blendEntries = useMemo(
    () => Object.entries(metrics?.blend || {}).sort((a, b) => b[1] - a[1]),
    [metrics],
  )
  const verticalGuide = PUBLIC_VERTICAL_GUIDES[vertical] || { status: 'beta', label: 'Beta', message: 'No public starter prompts configured for this vertical yet.', prompts: [] }

  async function handleSubmit() {
    if (auth?.supabaseEnabled && !auth?.user) {
      await auth.signInWithGoogle()
      return
    }
    setLoading(true)
    setError(null)
    setResult(null)
    try {
      const res = await api.formulate({
        input_text: inputText,
        vertical,
        batch_size: batchSize,
        opt_mode: optMode,
        constraints: {},
      })
      setResult(res.data)
    } catch (err) {
      setError(err.response?.data?.detail || err.message)
    } finally {
      setLoading(false)
    }
  }

  function loadPrompt(example) {
    setInputText(example.text)
    setVertical(example.vertical)
  }

  return (
    <div style={{ maxWidth: '1180px', margin: '0 auto' }}>
      <section style={{
        ...shell.panel,
        padding: '1.6rem',
        marginBottom: '1rem',
        background: 'radial-gradient(circle at top right, rgba(13,148,136,0.16), transparent 26%), linear-gradient(180deg, rgba(11,28,53,0.98), rgba(6,16,31,0.96))',
      }}>
        <div style={{ display: 'flex', gap: '1.25rem', flexWrap: 'wrap', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div style={{ maxWidth: '680px' }}>
            <div style={shell.label}>IntelliForm Studio</div>
            <h1 style={{ margin: '0.5rem 0 0.6rem', color: '#f8fafc', fontSize: '2.2rem', textAlign: 'left', lineHeight: 1.1 }}>
              Translate a market brief into an optimization-ready formulation program.
            </h1>
            <p style={{ color: '#9fb0c8', textAlign: 'left', lineHeight: 1.7, margin: 0 }}>
              This upgraded workspace applies FormulAI’s prompt-hardening lessons to IntelliForm’s deeper chemistry stack:
              cleaner parsing, canonical vertical mapping, clearer constraint visibility, and more trustworthy optimization output.
            </p>
          </div>
          <div style={{ display: 'grid', gap: '0.75rem', minWidth: '280px', flex: '1 1 280px' }}>
            <div style={{ ...shell.panel, padding: '1rem 1.1rem' }}>
              <div style={{ color: '#7c8aa5', fontSize: '0.68rem', letterSpacing: '0.14em', textTransform: 'uppercase' }}>Parser Reliability</div>
              <div style={{ color: '#f8fafc', marginTop: '0.35rem', fontWeight: 700 }}>Local JSON repair + canonical vertical inference</div>
            </div>
            <div style={{ ...shell.panel, padding: '1rem 1.1rem' }}>
              <div style={{ color: '#7c8aa5', fontSize: '0.68rem', letterSpacing: '0.14em', textTransform: 'uppercase' }}>Showcase Flow</div>
              <div style={{ color: '#f8fafc', marginTop: '0.35rem', fontWeight: 700 }}>Prompt cards, parsed brief panel, and proof-ready result sections</div>
            </div>
          </div>
        </div>
      </section>

      <div style={{ display: 'grid', gridTemplateColumns: 'minmax(320px, 420px) minmax(0, 1fr)', gap: '1rem' }}>
        <div>
          <Section eyebrow="Brief Builder" title="Describe the product target">
            {auth?.supabaseEnabled && !auth?.user ? (
              <div style={{
                ...shell.panel,
                padding: '0.95rem 1rem',
                background: 'rgba(217, 119, 6, 0.10)',
                border: '1px solid rgba(245, 158, 11, 0.28)',
                marginBottom: '1rem',
              }}>
                <div style={{ color: '#f8fafc', fontWeight: 700, marginBottom: '0.35rem' }}>Sign in to generate</div>
                <div style={{ color: '#9fb0c8', fontSize: '0.84rem', lineHeight: 1.6, marginBottom: '0.75rem' }}>
                  The public edition stays browseable, but live formulation generation is tied to a signed-in free account so usage, history, and future upgrades belong to you.
                </div>
                <button
                  onClick={() => auth.signInWithGoogle()}
                  style={{
                    background: '#0D9488',
                    color: '#fff',
                    border: 'none',
                    borderRadius: '10px',
                    padding: '0.7rem 1rem',
                    fontWeight: 700,
                    cursor: 'pointer',
                  }}
                >
                  Continue with Google
                </button>
              </div>
            ) : null}
            <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.55rem', color: '#94a3b8', fontSize: '0.82rem' }}>
              Natural-language formulation brief
            </label>
            <textarea
              value={inputText}
              onChange={(event) => setInputText(event.target.value)}
              placeholder="e.g. Low-VOC industrial degreaser for heavy equipment with high flash point, strong grease lift, and moderate foam."
              rows={7}
              style={{ ...shell.input, resize: 'vertical', minHeight: '180px', lineHeight: 1.7 }}
            />

            <div style={{ display: 'grid', gap: '0.85rem', gridTemplateColumns: 'repeat(2, minmax(0, 1fr))', marginTop: '1rem' }}>
              <div>
                <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.45rem', color: '#94a3b8', fontSize: '0.76rem' }}>Vertical</label>
                <select value={vertical} onChange={(event) => setVertical(event.target.value)} style={shell.input}>
                  {VERTICAL_OPTIONS.map((option) => (
                    <option key={option.value} value={option.value}>{option.label}</option>
                  ))}
                </select>
              </div>
              <div>
                <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.45rem', color: '#94a3b8', fontSize: '0.76rem' }}>Optimization mode</label>
                <select value={optMode} onChange={(event) => setOptMode(event.target.value)} style={shell.input}>
                  {optModes.map((option) => (
                    <option key={option.value} value={option.value}>{option.label}</option>
                  ))}
                </select>
              </div>
              <div style={{ gridColumn: '1 / -1' }}>
                <label style={{ display: 'block', textAlign: 'left', marginBottom: '0.45rem', color: '#94a3b8', fontSize: '0.76rem' }}>Batch size (kg)</label>
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
                background: loading ? '#334155' : 'linear-gradient(135deg, #0D9488, #0f766e)',
                color: '#fff',
                border: 'none',
                borderRadius: '14px',
                padding: '0.9rem 1.2rem',
                fontSize: '0.95rem',
                fontWeight: 700,
                cursor: loading || !canSubmit ? 'not-allowed' : 'pointer',
                width: '100%',
                boxShadow: loading ? 'none' : '0 14px 32px rgba(13, 148, 136, 0.25)',
              }}
            >
              {loading ? 'Running IntelliForm...' : 'Run IntelliForm'}
            </button>
          </Section>

          <Section eyebrow="Starter Prompts" title={`${formatVertical(vertical)} public demo prompts`}>
            <div style={{
              ...shell.panel,
              padding: '0.95rem 1rem',
              background: verticalGuide.status === 'validated' ? 'rgba(13, 148, 136, 0.10)' : 'rgba(217, 119, 6, 0.10)',
              border: verticalGuide.status === 'validated'
                ? '1px solid rgba(125, 211, 200, 0.28)'
                : '1px solid rgba(245, 158, 11, 0.28)',
              marginBottom: '0.9rem',
            }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', gap: '1rem', alignItems: 'center', marginBottom: '0.45rem' }}>
                <div style={{ color: '#f8fafc', fontWeight: 700 }}>{formatVertical(vertical)}</div>
                <div style={{
                  color: verticalGuide.status === 'validated' ? '#7dd3c8' : '#fbbf24',
                  fontSize: '0.72rem',
                  textTransform: 'uppercase',
                  letterSpacing: '0.12em',
                  fontWeight: 700,
                }}>
                  {verticalGuide.label}
                </div>
              </div>
              <div style={{ color: '#9fb0c8', fontSize: '0.84rem', lineHeight: 1.6 }}>
                {verticalGuide.message}
              </div>
            </div>

            {verticalGuide.prompts.length > 0 ? (
              <div style={{ display: 'grid', gap: '0.8rem' }}>
                {verticalGuide.prompts.map((example) => (
                  <button
                    key={example.title}
                    type="button"
                    onClick={() => loadPrompt(example)}
                    style={{
                      ...shell.panel,
                      padding: '0.95rem 1rem',
                      textAlign: 'left',
                      cursor: 'pointer',
                      background: 'rgba(8, 18, 35, 0.86)',
                    }}
                  >
                    <div style={{ color: '#f8fafc', fontWeight: 700, marginBottom: '0.25rem' }}>{example.title}</div>
                    <div style={{ color: '#7dd3c8', fontSize: '0.72rem', textTransform: 'uppercase', letterSpacing: '0.12em', marginBottom: '0.4rem' }}>
                      {formatVertical(example.vertical)}
                    </div>
                    <div style={{ color: '#93a4bd', fontSize: '0.84rem', lineHeight: 1.6, marginBottom: '0.5rem' }}>{example.text}</div>
                    {example.note ? (
                      <div style={{ color: '#64748b', fontSize: '0.76rem', lineHeight: 1.5 }}>
                        {example.note}
                      </div>
                    ) : null}
                  </button>
                ))}
              </div>
            ) : (
              <div style={{ color: '#94a3b8', fontSize: '0.86rem', lineHeight: 1.6 }}>
                No validated public starter prompts for this vertical yet. For a first successful demo, switch to Agricultural, Food & Beverage, or Fabric & Laundry.
              </div>
            )}
          </Section>
        </div>

        <div>
          {error ? (
            <div style={{
              ...shell.panel,
              padding: '1rem 1.1rem',
              color: '#fecaca',
              border: '1px solid rgba(220, 38, 38, 0.45)',
              background: 'rgba(69, 10, 10, 0.35)',
              marginBottom: '1rem',
            }}>
              {error}
            </div>
          ) : null}

          <div style={{ display: 'flex', gap: '0.9rem', flexWrap: 'wrap', marginBottom: '1rem' }}>
            <MetricCard label="Cost / kg" value={metrics?.cost_per_kg ? `$${metrics.cost_per_kg.toFixed(2)}` : '—'} tone="#7dd3c8" />
            <MetricCard label="Bio-based" value={metrics?.bio_pct ? `${metrics.bio_pct.toFixed(1)}%` : '—'} tone="#fbbf24" />
            <MetricCard label="Perf Score" value={metrics?.perf_score ? metrics.perf_score.toFixed(1) : '—'} unit="/ 100" />
            <MetricCard label="Eco Grade" value={result?.eco?.grade} tone="#7dd3c8" />
            <MetricCard label="Ingredient Pool" value={meta?.ingredient_pool_size ?? '—'} />
          </div>

          <Section
            eyebrow="Program Lens"
            title="Parsed brief and optimization posture"
            aside={meta ? (
              <div style={{
                borderRadius: '999px',
                padding: '0.35rem 0.75rem',
                background: 'rgba(13, 148, 136, 0.12)',
                border: '1px solid rgba(13, 148, 136, 0.4)',
                color: '#99f6e4',
                fontSize: '0.75rem',
                fontWeight: 700,
              }}>
                {parsed?.parser_backend || 'Parser pending'}
              </div>
            ) : null}
          >
            {parsed ? (
              <div style={{ display: 'grid', gap: '0.85rem', gridTemplateColumns: 'repeat(2, minmax(0, 1fr))' }}>
                <div style={{ ...shell.panel, padding: '0.95rem 1rem' }}>
                  <div style={shell.label}>Resolved Vertical</div>
                  <div style={{ color: '#f8fafc', marginTop: '0.45rem', fontWeight: 700 }}>{formatVertical(meta?.resolved_vertical)}</div>
                  <div style={{ color: '#64748b', marginTop: '0.3rem', fontSize: '0.8rem' }}>
                    Requested: {formatVertical(meta?.requested_vertical)} · Inferred: {formatVertical(meta?.inferred_vertical)}
                  </div>
                </div>
                <div style={{ ...shell.panel, padding: '0.95rem 1rem' }}>
                  <div style={shell.label}>Constraints Used</div>
                  <div style={{ color: '#f8fafc', marginTop: '0.45rem', fontWeight: 700 }}>
                    ${meta?.constraints_used?.max_cost}/kg · {meta?.constraints_used?.min_bio}% bio · {meta?.constraints_used?.min_perf} perf
                  </div>
                  <div style={{ color: '#64748b', marginTop: '0.3rem', fontSize: '0.8rem' }}>
                    Mode: {meta?.optimization_mode_requested}
                  </div>
                </div>
                <div style={{ ...shell.panel, padding: '0.95rem 1rem', gridColumn: '1 / -1' }}>
                  <div style={shell.label}>Parser Reasoning</div>
                  <div style={{ color: '#cbd5e1', marginTop: '0.55rem', lineHeight: 1.7, fontSize: '0.88rem' }}>
                    {parsed.reasoning}
                  </div>
                </div>
              </div>
            ) : (
              <div style={{ color: '#94a3b8', textAlign: 'left', lineHeight: 1.7 }}>
                Submit a brief to see how IntelliForm interpreted the request, which vertical it resolved to, and the exact optimization thresholds used.
              </div>
            )}
          </Section>

          <Section eyebrow="Blend Architecture" title="Optimized composition">
            {blendEntries.length ? (
              <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem' }}>
                {blendEntries.map(([ingredient, pct]) => (
                  <div key={ingredient}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', gap: '1rem', marginBottom: '0.35rem' }}>
                      <div style={{ color: '#dbe7f5', fontSize: '0.88rem', textAlign: 'left' }}>{ingredient}</div>
                      <div style={{ color: '#f8fafc', fontWeight: 700 }}>{pct}%</div>
                    </div>
                    <div style={{ background: 'rgba(30, 58, 95, 0.85)', borderRadius: '999px', height: '10px', overflow: 'hidden' }}>
                      <div style={{ width: `${Math.min(pct, 100)}%`, height: '10px', background: 'linear-gradient(90deg, #0D9488, #2dd4bf)' }} />
                    </div>
                  </div>
                ))}
              </div>
            ) : (
              <div style={{ color: '#94a3b8', textAlign: 'left' }}>
                No blend yet. The result will appear here after a successful run.
              </div>
            )}
          </Section>

          <Section eyebrow="System Readout" title="Risk, compliance, and readiness">
            <div style={{ display: 'grid', gap: '0.9rem', gridTemplateColumns: 'repeat(2, minmax(0, 1fr))' }}>
              <div style={{ ...shell.panel, padding: '1rem' }}>
                <div style={shell.label}>Regulatory Status</div>
                <div style={{ color: result?.vreg?.overall_status?.includes('✅') ? '#7dd3c8' : '#fbbf24', marginTop: '0.45rem', fontWeight: 700 }}>
                  {result?.vreg?.overall_status || 'Pending'}
                </div>
                <div style={{ color: '#64748b', marginTop: '0.35rem', fontSize: '0.8rem', lineHeight: 1.6 }}>
                  {result?.vreg?.framework || 'Framework notes will appear after formulation.'}
                </div>
              </div>
              <div style={{ ...shell.panel, padding: '1rem' }}>
                <div style={shell.label}>Stability</div>
                <div style={{ color: '#f8fafc', marginTop: '0.45rem', fontWeight: 700 }}>
                  {result?.stability?.overall_rating || 'Pending'}
                </div>
                <div style={{ color: '#64748b', marginTop: '0.35rem', fontSize: '0.8rem', lineHeight: 1.6 }}>
                  Shelf life: {result?.stability?.shelf_life_range || '—'} · pH {result?.stability?.ph_min ?? '—'}–{result?.stability?.ph_max ?? '—'}
                </div>
              </div>
            </div>

            {metrics?.warnings?.length ? (
              <div style={{ marginTop: '1rem' }}>
                {metrics.warnings.map((warning) => (
                  <div key={warning} style={{ color: '#fbbf24', textAlign: 'left', fontSize: '0.84rem', marginBottom: '0.35rem' }}>
                    ⚠ {warning}
                  </div>
                ))}
              </div>
            ) : null}

            {metrics?.compliance_flags?.length ? (
              <div style={{ marginTop: '0.8rem' }}>
                {metrics.compliance_flags.map((flag) => (
                  <div key={flag} style={{ color: '#fecaca', textAlign: 'left', fontSize: '0.84rem', marginBottom: '0.35rem' }}>
                    ✕ {flag}
                  </div>
                ))}
              </div>
            ) : null}
          </Section>

          <Section eyebrow="Agent Commentary" title="Commercial and technical review">
            {result?.agents?.length ? (
              <div style={{ display: 'grid', gap: '0.7rem' }}>
                {result.agents.map((agent) => (
                  <div key={agent} style={{ ...shell.panel, padding: '0.95rem 1rem', color: '#cbd5e1', textAlign: 'left', lineHeight: 1.7 }}>
                    {agent}
                  </div>
                ))}
              </div>
            ) : (
              <div style={{ color: '#94a3b8', textAlign: 'left' }}>
                Agent commentary will appear after IntelliForm completes a formulation run.
              </div>
            )}
          </Section>
        </div>
      </div>
    </div>
  )
}
