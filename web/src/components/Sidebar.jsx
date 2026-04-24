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

export default function Sidebar({ activePage, setActivePage, auth, isMobile = false }) {
  if (isMobile) {
    return (
      <aside style={{
        background: '#0D1F3C',
        borderBottom: '1px solid #1e3a5f',
        display: 'flex',
        flexDirection: 'column',
        padding: '0.9rem 0 0.8rem'
      }}>
        <div style={{ padding: '0 1rem 0.9rem' }}>
          <div style={{ color: '#0D9488', fontWeight: 800, fontSize: '1.1rem', letterSpacing: '-0.4px', textAlign: 'left' }}>
            IntelliForm™
          </div>
          <div style={{ color: '#64748b', fontSize: '0.68rem', marginTop: '2px', textAlign: 'left' }}>
            Speed with Scientific Rigor
          </div>
        </div>
        <nav style={{
          display: 'flex',
          gap: '0.55rem',
          overflowX: 'auto',
          padding: '0 1rem 0.9rem',
          scrollbarWidth: 'none',
        }}>
          {NAV.map(item => (
            <button
              key={item.id}
              onClick={() => setActivePage(item.id)}
              style={{
                flex: '0 0 auto',
                padding: '0.65rem 0.9rem',
                textAlign: 'left',
                background: activePage === item.id ? '#0D9488' + '22' : 'transparent',
                color: activePage === item.id ? '#0D9488' : '#94a3b8',
                border: `1px solid ${activePage === item.id ? '#0D9488' : '#1e3a5f'}`,
                borderRadius: '999px',
                cursor: 'pointer',
                fontSize: '0.8rem',
                fontWeight: activePage === item.id ? 700 : 500,
                whiteSpace: 'nowrap',
              }}
            >
              {item.label}
            </button>
          ))}
        </nav>
        <div style={{ padding: '0 1rem', color: '#334155', fontSize: '0.7rem' }}>
          {auth?.supabaseEnabled ? (
            <div style={{ marginBottom: '0.65rem' }}>
              {auth.user ? (
                <>
                  <div style={{ color: '#94a3b8', fontSize: '0.72rem', marginBottom: '0.35rem', textAlign: 'left', wordBreak: 'break-word' }}>
                    {auth.user.email}
                  </div>
                  <button
                    onClick={auth.signOut}
                    style={{
                      background: 'transparent',
                      color: '#94a3b8',
                      border: '1px solid #1e3a5f',
                      borderRadius: '8px',
                      padding: '0.55rem 0.75rem',
                      fontSize: '0.72rem',
                      cursor: 'pointer',
                      width: '100%',
                      textAlign: 'left',
                    }}
                  >
                    Sign out
                  </button>
                </>
              ) : (
                <button
                  onClick={auth.signInWithGoogle}
                  style={{
                    background: '#0D9488',
                    color: '#fff',
                    border: 'none',
                    borderRadius: '8px',
                    padding: '0.65rem 0.85rem',
                    fontSize: '0.78rem',
                    cursor: 'pointer',
                    width: '100%',
                    textAlign: 'left',
                  }}
                >
                  Sign in with Google
                </button>
              )}
            </div>
          ) : null}
          <div style={{ textAlign: 'left' }}>v2.1.0 · ChemeNova LLC</div>
        </div>
      </aside>
    )
  }

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
        {auth?.supabaseEnabled ? (
          <div style={{ marginBottom: '0.75rem' }}>
            {auth.user ? (
              <>
                <div style={{ color: '#94a3b8', fontSize: '0.72rem', marginBottom: '0.35rem' }}>{auth.user.email}</div>
                <button
                  onClick={auth.signOut}
                  style={{
                    background: 'transparent',
                    color: '#94a3b8',
                    border: '1px solid #1e3a5f',
                    borderRadius: '8px',
                    padding: '0.45rem 0.7rem',
                    fontSize: '0.72rem',
                    cursor: 'pointer',
                    width: '100%',
                    textAlign: 'left',
                  }}
                >
                  Sign out
                </button>
              </>
            ) : (
              <button
                onClick={auth.signInWithGoogle}
                style={{
                  background: '#0D9488',
                  color: '#fff',
                  border: 'none',
                  borderRadius: '8px',
                  padding: '0.5rem 0.75rem',
                  fontSize: '0.72rem',
                  cursor: 'pointer',
                  width: '100%',
                  textAlign: 'left',
                }}
              >
                Sign in with Google
              </button>
            )}
          </div>
        ) : null}
        v2.1.0 · ChemeNova LLC
      </div>
    </aside>
  )
}
