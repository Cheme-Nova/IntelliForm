import { useState } from 'react'
import { api } from '../api/client'
import { BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, Cell } from 'recharts'

const card = { background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '8px', padding: '1rem' }
const input = { width: '100%', background: '#0D1F3C', border: '1px solid #1e3a5f', borderRadius: '6px', color: '#fff', padding: '0.55rem 0.8rem', fontSize: '0.85rem', boxSizing: 'border-box', fontFamily: 'inherit' }
const label = { color: '#64748b', fontSize: '0.72rem', display: 'block', marginBottom: '4px' }
const btn = (disabled) => ({ background: disabled ? '#334155' : '#0D9488', color: '#fff', border: 'none', borderRadius: '7px', padding: '0.6rem 1.4rem', fontSize: '0.88rem', fontWeight: 600, cursor: disabled ? 'not-allowed' : 'pointer' })

const CERTS = ['ISO9001', 'ISO14001', 'GMP', 'Ecocert', 'COSMOS', 'REACH', 'Halal', 'Kosher', 'Fair Trade', 'USDA Organic', 'RSB', 'RSPO']
const AVAIL = ['in_stock', 'limited', 'out_of_stock']

function Section({ title, children }) {
  return (
    <div style={{ marginBottom: '1.5rem' }}>
      <div style={{ color: '#D97706', fontWeight: 700, fontSize: '0.75rem', letterSpacing: '0.08em', textTransform: 'uppercase', marginBottom: '0.75rem' }}>{title}</div>
      {children}
    </div>
  )
}

function Tag({ label: lbl, active, onClick }) {
  return (
    <button onClick={onClick} style={{
      padding: '0.28rem 0.65rem', fontSize: '0.72rem', fontWeight: 500,
      borderRadius: '999px', cursor: 'pointer', border: '1px solid',
      background: active ? '#0D9488' : 'transparent',
      color: active ? '#fff' : '#64748b',
      borderColor: active ? '#0D9488' : '#1e3a5f',
    }}>{lbl}</button>
  )
}

// ── Register tab ──────────────────────────────────────────────────────────────
function RegisterTab() {
  const [form, setForm] = useState({ name: '', email: '', country: '', description: '' })
  const [certs, setCerts] = useState([])
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)

  async function handleSubmit(e) {
    e.preventDefault()
    if (!form.name || !form.email || !form.country) { setError('Name, email, and country are required.'); return }
    setLoading(true); setError(null)
    try {
      const res = await api.supplierRegister({ ...form, certifications: certs })
      setResult(res.data)
    } catch (err) {
      setError(err.response?.data?.detail || err.message)
    } finally { setLoading(false) }
  }

  if (result) return (
    <div style={{ ...card, borderColor: '#0D9488' }}>
      <div style={{ color: '#0D9488', fontWeight: 700, fontSize: '1rem', marginBottom: '0.5rem' }}>Application submitted</div>
      <div style={{ color: '#94a3b8', fontSize: '0.88rem', lineHeight: 1.7 }}>
        Your supplier ID is <span style={{ color: '#D97706', fontFamily: 'monospace' }}>{result.id}</span>.<br />
        Status: <strong>{result.status}</strong> — ChemeNova will email your API key once approved (typically 1 business day).
      </div>
    </div>
  )

  return (
    <form onSubmit={handleSubmit} style={{ display: 'flex', flexDirection: 'column', gap: '1rem', maxWidth: '600px' }}>
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem' }}>
        <div><span style={label}>Company name *</span><input style={input} value={form.name} onChange={e => setForm(f => ({ ...f, name: e.target.value }))} required /></div>
        <div><span style={label}>Contact email *</span><input style={input} type="email" value={form.email} onChange={e => setForm(f => ({ ...f, email: e.target.value }))} required /></div>
      </div>
      <div><span style={label}>Country *</span><input style={input} value={form.country} onChange={e => setForm(f => ({ ...f, country: e.target.value }))} placeholder="e.g. United States" required /></div>
      <div>
        <span style={label}>Certifications</span>
        <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.4rem', marginTop: '4px' }}>
          {CERTS.map(c => <Tag key={c} label={c} active={certs.includes(c)} onClick={() => setCerts(cs => cs.includes(c) ? cs.filter(x => x !== c) : [...cs, c])} />)}
        </div>
      </div>
      <div><span style={label}>Company description</span><textarea style={{ ...input, resize: 'vertical', minHeight: '72px' }} value={form.description} onChange={e => setForm(f => ({ ...f, description: e.target.value }))} /></div>
      {error && <div style={{ color: '#fca5a5', fontSize: '0.83rem' }}>{error}</div>}
      <button type="submit" style={{ ...btn(loading), alignSelf: 'flex-start' }} disabled={loading}>
        {loading ? 'Submitting…' : 'Apply to Join →'}
      </button>
    </form>
  )
}

// ── Listings tab ──────────────────────────────────────────────────────────────
function ListingsTab() {
  const [key, setKey] = useState('')
  const [supplierId, setSupplierId] = useState('')
  const [listings, setListings] = useState(null)
  const [authError, setAuthError] = useState(null)
  const [form, setForm] = useState({ ingredient_name: '', price_per_kg: '', moq_kg: '25', lead_time_days: '14', availability: 'in_stock', valid_until: '' })
  const [formCerts, setFormCerts] = useState([])
  const [loading, setLoading] = useState(false)
  const [saveMsg, setSaveMsg] = useState(null)

  async function loadListings() {
    if (!key || !supplierId) return
    setAuthError(null)
    try {
      const res = await api.supplierListings(supplierId, key)
      setListings(res.data)
    } catch (err) {
      setAuthError(err.response?.data?.detail || 'Invalid key or supplier ID.')
      setListings(null)
    }
  }

  async function saveListing(e) {
    e.preventDefault()
    setLoading(true); setSaveMsg(null)
    try {
      await api.supplierSubmitListing(supplierId, key, {
        ...form,
        price_per_kg: parseFloat(form.price_per_kg),
        moq_kg: parseFloat(form.moq_kg),
        lead_time_days: parseInt(form.lead_time_days),
        certifications: formCerts,
        valid_until: form.valid_until || null,
      })
      setSaveMsg('Listing saved.')
      await loadListings()
    } catch (err) {
      setSaveMsg(err.response?.data?.detail || err.message)
    } finally { setLoading(false) }
  }

  async function deleteListing(listingId) {
    try {
      await api.supplierDeleteListing(supplierId, key, listingId)
      await loadListings()
    } catch { /* ignore */ }
  }

  return (
    <div style={{ maxWidth: '720px' }}>
      <Section title="Authenticate">
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.75rem', marginBottom: '0.75rem' }}>
          <div><span style={label}>Supplier ID</span><input style={input} value={supplierId} onChange={e => setSupplierId(e.target.value)} placeholder="e.g. a1b2c3" /></div>
          <div><span style={label}>API key (ifs_…)</span><input style={{ ...input, fontFamily: 'monospace' }} type="password" value={key} onChange={e => setKey(e.target.value)} /></div>
        </div>
        <button style={btn(false)} onClick={loadListings}>Load my listings →</button>
        {authError && <div style={{ color: '#fca5a5', fontSize: '0.83rem', marginTop: '0.5rem' }}>{authError}</div>}
      </Section>

      {listings !== null && (
        <>
          <Section title="Add / Update Listing">
            <form onSubmit={saveListing} style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem' }}>
              <div><span style={label}>Ingredient name *</span><input style={input} value={form.ingredient_name} onChange={e => setForm(f => ({ ...f, ingredient_name: e.target.value }))} placeholder="Must match IntelliForm DB name" required /></div>
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3, 1fr)', gap: '0.75rem' }}>
                <div><span style={label}>Price ($/kg) *</span><input style={input} type="number" min="0.01" step="0.01" value={form.price_per_kg} onChange={e => setForm(f => ({ ...f, price_per_kg: e.target.value }))} required /></div>
                <div><span style={label}>MOQ (kg)</span><input style={input} type="number" min="0.1" step="1" value={form.moq_kg} onChange={e => setForm(f => ({ ...f, moq_kg: e.target.value }))} /></div>
                <div><span style={label}>Lead time (days)</span><input style={input} type="number" min="1" step="1" value={form.lead_time_days} onChange={e => setForm(f => ({ ...f, lead_time_days: e.target.value }))} /></div>
              </div>
              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.75rem' }}>
                <div>
                  <span style={label}>Availability</span>
                  <select style={input} value={form.availability} onChange={e => setForm(f => ({ ...f, availability: e.target.value }))}>
                    {AVAIL.map(a => <option key={a} value={a}>{a.replace('_', ' ')}</option>)}
                  </select>
                </div>
                <div><span style={label}>Valid until (YYYY-MM-DD)</span><input style={input} type="date" value={form.valid_until} onChange={e => setForm(f => ({ ...f, valid_until: e.target.value }))} /></div>
              </div>
              <div>
                <span style={label}>Certifications</span>
                <div style={{ display: 'flex', flexWrap: 'wrap', gap: '0.35rem', marginTop: '4px' }}>
                  {CERTS.map(c => <Tag key={c} label={c} active={formCerts.includes(c)} onClick={() => setFormCerts(cs => cs.includes(c) ? cs.filter(x => x !== c) : [...cs, c])} />)}
                </div>
              </div>
              <div style={{ display: 'flex', gap: '0.75rem', alignItems: 'center' }}>
                <button type="submit" style={btn(loading)} disabled={loading}>{loading ? 'Saving…' : 'Save listing'}</button>
                {saveMsg && <span style={{ color: saveMsg === 'Listing saved.' ? '#0D9488' : '#fca5a5', fontSize: '0.83rem' }}>{saveMsg}</span>}
              </div>
            </form>
          </Section>

          <Section title={`Current listings (${listings.length})`}>
            {listings.length === 0 ? (
              <div style={{ ...card, color: '#94a3b8', textAlign: 'center', fontSize: '0.85rem' }}>No listings yet — add one above.</div>
            ) : (
              <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem' }}>
                {listings.map(l => (
                  <div key={l.id} style={{ ...card, display: 'flex', justifyContent: 'space-between', alignItems: 'center', flexWrap: 'wrap', gap: '0.5rem' }}>
                    <div style={{ minWidth: 0 }}>
                      <div style={{ color: '#fff', fontWeight: 600, fontSize: '0.9rem' }}>{l.ingredient_name}</div>
                      <div style={{ color: '#64748b', fontSize: '0.75rem', marginTop: '2px' }}>
                        ${l.price_per_kg}/kg · MOQ {l.moq_kg}kg · {l.lead_time_days}d lead · {l.availability.replace('_', ' ')}
                        {l.certifications?.length ? ` · ${l.certifications.join(', ')}` : ''}
                      </div>
                    </div>
                    <button onClick={() => deleteListing(l.id)} style={{ ...btn(false), background: '#7f1d1d', padding: '0.35rem 0.75rem', fontSize: '0.75rem' }}>Remove</button>
                  </div>
                ))}
              </div>
            )}
          </Section>
        </>
      )}
    </div>
  )
}

// ── Demand tab ────────────────────────────────────────────────────────────────
function DemandTab() {
  const [key, setKey] = useState('')
  const [supplierId, setSupplierId] = useState('')
  const [signals, setSignals] = useState(null)
  const [error, setError] = useState(null)

  async function load() {
    setError(null)
    try {
      const res = await api.supplierDemand(supplierId, key)
      setSignals(res.data.signals || [])
    } catch (err) {
      setError(err.response?.data?.detail || 'Invalid credentials.')
      setSignals(null)
    }
  }

  const avColor = { in_stock: '#0D9488', limited: '#D97706', out_of_stock: '#ef4444' }

  return (
    <div style={{ maxWidth: '720px' }}>
      <p style={{ color: '#64748b', fontSize: '0.85rem', marginBottom: '1.25rem' }}>
        See how many saved formulations on the platform include each ingredient you supply.
      </p>
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.75rem', marginBottom: '0.75rem' }}>
        <div><span style={label}>Supplier ID</span><input style={input} value={supplierId} onChange={e => setSupplierId(e.target.value)} /></div>
        <div><span style={label}>API key (ifs_…)</span><input style={{ ...input, fontFamily: 'monospace' }} type="password" value={key} onChange={e => setKey(e.target.value)} /></div>
      </div>
      <button style={{ ...btn(false), marginBottom: '1.25rem' }} onClick={load}>Load demand signals →</button>
      {error && <div style={{ color: '#fca5a5', fontSize: '0.83rem', marginBottom: '1rem' }}>{error}</div>}

      {signals !== null && (
        signals.length === 0 ? (
          <div style={{ ...card, color: '#94a3b8', textAlign: 'center' }}>No listings yet — add some in Manage Listings.</div>
        ) : (
          <>
            <ResponsiveContainer width="100%" height={220}>
              <BarChart data={signals} margin={{ top: 0, right: 0, bottom: 30, left: 0 }}>
                <XAxis dataKey="ingredient" tick={{ fill: '#64748b', fontSize: 11 }} angle={-30} textAnchor="end" interval={0} />
                <YAxis tick={{ fill: '#64748b', fontSize: 11 }} allowDecimals={false} />
                <Tooltip contentStyle={{ background: '#0D1F3C', border: '1px solid #1e3a5f', color: '#fff' }} />
                <Bar dataKey="formulation_count" name="Formulations" radius={[4, 4, 0, 0]}>
                  {signals.map((s, i) => <Cell key={i} fill={avColor[s.availability] || '#0D9488'} />)}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
            <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem', marginTop: '1rem' }}>
              {signals.map(s => (
                <div key={s.listing_id} style={{ ...card, display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                  <div>
                    <div style={{ color: '#fff', fontWeight: 500, fontSize: '0.88rem' }}>{s.ingredient}</div>
                    <div style={{ color: '#64748b', fontSize: '0.72rem' }}>${s.price_per_kg}/kg · {s.availability.replace('_', ' ')}</div>
                  </div>
                  <div style={{ textAlign: 'right' }}>
                    <div style={{ color: '#0D9488', fontWeight: 700, fontSize: '1.1rem' }}>{s.formulation_count}</div>
                    <div style={{ color: '#64748b', fontSize: '0.65rem' }}>formulations</div>
                  </div>
                </div>
              ))}
            </div>
          </>
        )
      )}
    </div>
  )
}

// ── Main component ────────────────────────────────────────────────────────────
export default function Suppliers() {
  const [tab, setTab] = useState('register')
  const tabs = [
    { id: 'register', label: '📋 Register' },
    { id: 'listings', label: '📦 Manage Listings' },
    { id: 'demand',   label: '📈 Demand Signals' },
  ]

  return (
    <div style={{ maxWidth: '800px' }}>
      <h1 style={{ color: '#0D9488', fontSize: '1.8rem', fontWeight: 800, marginBottom: '0.25rem' }}>🏭 Supplier Portal</h1>
      <p style={{ color: '#64748b', marginBottom: '1.5rem', fontSize: '0.9rem', lineHeight: 1.6 }}>
        List your ingredients on IntelliForm and reach every formulator on the platform.
        Register below — ChemeNova reviews applications within 1 business day and emails your access key.
      </p>

      <div style={{ display: 'flex', gap: '0.4rem', marginBottom: '1.75rem', borderBottom: '1px solid #1e3a5f', paddingBottom: '0.75rem' }}>
        {tabs.map(t => (
          <button key={t.id} onClick={() => setTab(t.id)} style={{
            background: tab === t.id ? '#0D9488' : 'transparent',
            color: tab === t.id ? '#fff' : '#64748b',
            border: `1px solid ${tab === t.id ? '#0D9488' : '#1e3a5f'}`,
            borderRadius: '7px', padding: '0.45rem 0.9rem', fontSize: '0.82rem',
            fontWeight: tab === t.id ? 600 : 400, cursor: 'pointer',
          }}>{t.label}</button>
        ))}
      </div>

      {tab === 'register' && <RegisterTab />}
      {tab === 'listings' && <ListingsTab />}
      {tab === 'demand'   && <DemandTab />}
    </div>
  )
}
