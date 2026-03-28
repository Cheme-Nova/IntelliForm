const NAV = [
  { id: 'formulate', label: '⚗️ Formulate', },
  { id: 'eco', label: '🌱 EcoMetrics' },
  { id: 'certifications', label: '✅ Certifications' },
  { id: 'stability', label: '🧪 Stability' },
  { id: 'carbon', label: '♻️ Carbon' },
  { id: 'regulatory', label: '📋 Regulatory' },
  { id: 'pareto', label: '📊 Pareto' },
  { id: 'qsar', label: '🔬 QSAR' },
  { id: 'reformulation', label: '🔧 Reformulation' },
  { id: 'memory', label: '🧠 Memory' },
]

export default function Sidebar({ activePage, setActivePage }) {
  return (
    <aside style={{
      width: '220px',
      background: '#0D1F3C',
      borderRight: '1px solid #1e3a5f',
      display: 'flex',
      flexDirection: 'column',
      padding: '1.5rem 0'
    }}>
      <div style={{ padding: '0 1.5rem 2rem' }}>
        <div style={{ color: '#0D9488', fontWeight: 800, fontSize: '1.2rem', letterSpacing: '-0.5px' }}>
          IntelliForm™
        </div>
        <div style={{ color: '#64748b', fontSize: '0.7rem', marginTop: '2px' }}>
          Speed with Scientific Rigor
        </div>
      </div>
      <nav style={{ flex: 1 }}>
        {NAV.map(item => (
          <button
            key={item.id}
            onClick={() => setActivePage(item.id)}
            style={{
              display: 'block',
              width: '100%',
              padding: '0.65rem 1.5rem',
              textAlign: 'left',
              background: activePage === item.id ? '#0D9488' + '22' : 'transparent',
              color: activePage === item.id ? '#0D9488' : '#94a3b8',
              border: 'none',
              borderLeft: activePage === item.id ? '3px solid #0D9488' : '3px solid transparent',
              cursor: 'pointer',
              fontSize: '0.85rem',
              fontWeight: activePage === item.id ? 600 : 400,
              transition: 'all 0.15s'
            }}
          >
            {item.label}
          </button>
        ))}
      </nav>
      <div style={{ padding: '1rem 1.5rem', borderTop: '1px solid #1e3a5f', color: '#334155', fontSize: '0.7rem' }}>
        v2.1.0 · ChemeNova LLC
      </div>
    </aside>
  )
}
