import {
  BarChart3,
  Beaker,
  Brain,
  CheckCircle2,
  ClipboardCheck,
  FlaskConical,
  Gauge,
  GitBranch,
  Leaf,
  RefreshCcw,
  Sparkles,
} from 'lucide-react'

const NAV = [
  { id: 'formulate', label: 'Formulate', icon: FlaskConical },
  { id: 'eco', label: 'EcoMetrics', icon: Leaf },
  { id: 'certifications', label: 'Certifications', icon: CheckCircle2 },
  { id: 'stability', label: 'Stability', icon: Beaker },
  { id: 'carbon', label: 'Carbon', icon: Gauge },
  { id: 'regulatory', label: 'Regulatory', icon: ClipboardCheck },
  { id: 'pareto', label: 'Pareto', icon: BarChart3 },
  { id: 'qsar', label: 'QSAR', icon: GitBranch },
  { id: 'reformulation', label: 'Reformulation', icon: RefreshCcw },
  { id: 'memory', label: 'Memory', icon: Brain },
]

const colors = {
  ink: '#101828',
  muted: '#667085',
  line: '#E4E7EC',
  panel: '#FFFFFF',
  accent: '#127C67',
  accentSoft: '#EAF7F2',
}

function NavButton({ item, active, onClick, isMobile }) {
  const Icon = item.icon
  return (
    <button
      key={item.id}
      onClick={onClick}
      title={item.label}
      style={{
        flex: isMobile ? '0 0 auto' : 'none',
        width: isMobile ? 'auto' : '100%',
        display: 'flex',
        alignItems: 'center',
        gap: '0.6rem',
        padding: isMobile ? '0.62rem 0.8rem' : '0.62rem 0.85rem',
        textAlign: 'left',
        background: active ? colors.accentSoft : 'transparent',
        color: active ? colors.accent : colors.muted,
        border: `1px solid ${active ? '#B9E1D3' : 'transparent'}`,
        borderRadius: '8px',
        cursor: 'pointer',
        fontSize: '0.84rem',
        fontWeight: active ? 750 : 550,
        whiteSpace: 'nowrap',
      }}
    >
      <Icon size={16} strokeWidth={2} />
      <span>{item.label}</span>
    </button>
  )
}

export default function Sidebar({ activePage, setActivePage, auth, isMobile = false }) {
  if (isMobile) {
    return (
      <aside style={{
        background: 'rgba(255,255,255,0.96)',
        borderBottom: `1px solid ${colors.line}`,
        display: 'flex',
        flexDirection: 'column',
        padding: '0.85rem 0 0.75rem',
        position: 'sticky',
        top: 0,
        zIndex: 20,
        backdropFilter: 'blur(18px)',
      }}>
        <div style={{ padding: '0 1rem 0.75rem', display: 'flex', justifyContent: 'space-between', gap: '1rem', alignItems: 'center' }}>
          <div>
            <div style={{ color: colors.ink, fontWeight: 850, fontSize: '1rem', letterSpacing: '0', textAlign: 'left' }}>
              IntelliForm
            </div>
            <div style={{ color: colors.muted, fontSize: '0.72rem', marginTop: '2px', textAlign: 'left' }}>
              Formulation intelligence
            </div>
          </div>
          <div title="Agentic workflow" style={{
            display: 'inline-flex',
            alignItems: 'center',
            gap: '0.35rem',
            color: colors.accent,
            background: colors.accentSoft,
            border: '1px solid #B9E1D3',
            borderRadius: '999px',
            padding: '0.4rem 0.65rem',
            fontSize: '0.72rem',
            fontWeight: 750,
          }}>
            <Sparkles size={14} /> Live
          </div>
        </div>
        <nav style={{
          display: 'flex',
          gap: '0.45rem',
          overflowX: 'auto',
          padding: '0 1rem',
          scrollbarWidth: 'none',
        }}>
          {NAV.map(item => (
            <NavButton
              key={item.id}
              item={item}
              active={activePage === item.id}
              onClick={() => setActivePage(item.id)}
              isMobile
            />
          ))}
        </nav>
      </aside>
    )
  }

  return (
    <aside style={{
      width: '248px',
      background: colors.panel,
      borderRight: `1px solid ${colors.line}`,
      display: 'flex',
      flexDirection: 'column',
      padding: '1.2rem',
      position: 'sticky',
      top: 0,
      minHeight: '100vh',
      boxSizing: 'border-box',
    }}>
      <div style={{ padding: '0.2rem 0.2rem 1.4rem' }}>
        <div style={{
          width: '36px',
          height: '36px',
          borderRadius: '8px',
          background: colors.accent,
          display: 'grid',
          placeItems: 'center',
          color: '#fff',
          marginBottom: '0.8rem',
        }}>
          <Sparkles size={18} />
        </div>
        <div style={{ color: colors.ink, fontWeight: 850, fontSize: '1.12rem', letterSpacing: '0' }}>
          IntelliForm
        </div>
        <div style={{ color: colors.muted, fontSize: '0.76rem', marginTop: '4px', lineHeight: 1.4 }}>
          Agentic formulation intelligence
        </div>
      </div>
      <nav style={{ flex: 1, display: 'grid', alignContent: 'start', gap: '0.22rem' }}>
        {NAV.map(item => (
          <NavButton
            key={item.id}
            item={item}
            active={activePage === item.id}
            onClick={() => setActivePage(item.id)}
          />
        ))}
      </nav>
      <div style={{
        padding: '0.9rem',
        border: `1px solid ${colors.line}`,
        borderRadius: '8px',
        color: colors.muted,
        fontSize: '0.74rem',
        lineHeight: 1.5,
        background: '#F9FAFB',
      }}>
        {auth?.supabaseEnabled ? (
          <div style={{ marginBottom: '0.75rem' }}>
            {auth.user ? (
              <>
                <div style={{ color: colors.ink, fontSize: '0.72rem', marginBottom: '0.45rem', wordBreak: 'break-word' }}>{auth.user.email}</div>
                <button
                  onClick={auth.signOut}
                  style={{
                    background: '#fff',
                    color: colors.muted,
                    border: `1px solid ${colors.line}`,
                    borderRadius: '8px',
                    padding: '0.48rem 0.65rem',
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
                  background: colors.accent,
                  color: '#fff',
                  border: 'none',
                  borderRadius: '8px',
                  padding: '0.55rem 0.75rem',
                  fontSize: '0.74rem',
                  cursor: 'pointer',
                  width: '100%',
                  textAlign: 'left',
                  fontWeight: 700,
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
