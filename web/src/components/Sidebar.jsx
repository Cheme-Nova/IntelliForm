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
  Microscope,
  RefreshCcw,
  Sparkles,
} from 'lucide-react'

const NAV = [
  { id: 'formulate',      label: 'Formulate',     icon: FlaskConical  },
  { id: 'eco',            label: 'EcoMetrics',     icon: Leaf          },
  { id: 'certifications', label: 'Certifications', icon: CheckCircle2  },
  { id: 'stability',      label: 'Stability',      icon: Beaker        },
  { id: 'carbon',         label: 'Carbon',         icon: Gauge         },
  { id: 'regulatory',     label: 'Regulatory',     icon: ClipboardCheck },
  { id: 'pareto',         label: 'Pareto',         icon: BarChart3     },
  { id: 'qsar',           label: 'QSAR',           icon: GitBranch     },
  { id: 'reformulation',  label: 'Reformulation',  icon: RefreshCcw    },
  { id: 'pharma',         label: 'Pharma',         icon: Microscope    },
  { id: 'memory',         label: 'Memory',         icon: Brain         },
]

const c = {
  ink:        '#111827',
  muted:      '#5f6f7d',
  steel:      '#4c6375',
  line:       '#dce6ee',
  surface:    '#ffffff',
  surface2:   '#f1f5f8',
  teal:       '#1f8a7a',
  tealSoft:   '#eef6f3',
  tealBorder: '#b9d1ca',
}

function NavButton({ item, active, onClick, isMobile }) {
  const Icon = item.icon
  return (
    <button
      onClick={onClick}
      title={item.label}
      style={{
        flex: isMobile ? '0 0 auto' : 'none',
        width: isMobile ? 'auto' : '100%',
        display: 'flex',
        alignItems: 'center',
        gap: '0.55rem',
        padding: isMobile ? '0.55rem 0.75rem' : '0.58rem 0.8rem',
        textAlign: 'left',
        background: active ? c.tealSoft : 'transparent',
        color: active ? c.teal : c.steel,
        border: `1px solid ${active ? c.tealBorder : 'transparent'}`,
        borderRadius: '6px',
        cursor: 'pointer',
        fontSize: '0.83rem',
        fontWeight: active ? 750 : 550,
        whiteSpace: 'nowrap',
        minHeight: '40px',
        minWidth: 0,
        letterSpacing: 0,
      }}
    >
      <Icon size={15} strokeWidth={active ? 2.2 : 1.8} style={{ flex: '0 0 auto', opacity: active ? 1 : 0.75 }} />
      <span>{item.label}</span>
    </button>
  )
}

export default function Sidebar({ activePage, setActivePage, auth, isMobile = false }) {
  if (isMobile) {
    return (
      <aside className="if-sidebar-mobile" style={{
        display: 'flex',
        flexDirection: 'column',
        padding: '0.8rem 0 0.65rem',
        position: 'sticky',
        top: 0,
        zIndex: 20,
      }}>
        <div style={{ padding: '0 0.85rem 0.65rem', display: 'flex', justifyContent: 'space-between', gap: '0.75rem', alignItems: 'center' }}>
          <div style={{ minWidth: 0 }}>
            <div style={{ color: c.ink, fontWeight: 840, fontSize: '1rem', letterSpacing: 0, textAlign: 'left' }}>
              IntelliForm
            </div>
            <div style={{ color: c.muted, fontSize: '0.7rem', marginTop: '1px', textAlign: 'left' }}>
              ChemeNova · Formulation intelligence
            </div>
          </div>
          <div style={{
            display: 'inline-flex', alignItems: 'center', gap: '0.3rem',
            color: c.teal, background: c.tealSoft, border: `1px solid ${c.tealBorder}`,
            borderRadius: '999px', padding: '0.35rem 0.6rem',
            fontSize: '0.7rem', fontWeight: 750, flex: '0 0 auto',
          }}>
            <Sparkles size={12} /> Streaming
          </div>
        </div>
        <nav style={{
          display: 'flex', gap: '0.4rem', overflowX: 'auto',
          padding: '0 0.85rem', scrollbarWidth: 'none',
          WebkitOverflowScrolling: 'touch', touchAction: 'pan-x',
        }}>
          {NAV.map((item) => (
            <NavButton key={item.id} item={item} active={activePage === item.id} onClick={() => setActivePage(item.id)} isMobile />
          ))}
        </nav>
      </aside>
    )
  }

  return (
    <aside className="if-sidebar" style={{
      width: '240px',
      display: 'flex',
      flexDirection: 'column',
      padding: '1.2rem',
      position: 'sticky',
      top: 0,
      minHeight: '100vh',
      boxSizing: 'border-box',
    }}>
      {/* Brand */}
      <div style={{ padding: '0.1rem 0.2rem 1.4rem' }}>
        <div style={{
          width: '34px', height: '34px', borderRadius: '8px',
          background: c.teal, display: 'grid', placeItems: 'center',
          color: '#fff', marginBottom: '0.75rem',
          boxShadow: '0 4px 14px rgba(31,138,122,0.28)',
        }}>
          <Sparkles size={17} />
        </div>
        <div style={{ color: c.ink, fontWeight: 840, fontSize: '1.08rem', letterSpacing: 0 }}>
          IntelliForm
        </div>
        <div style={{ color: c.muted, fontSize: '0.72rem', marginTop: '3px', lineHeight: 1.4 }}>
          ChemeNova · Formulation intelligence
        </div>
      </div>

      {/* Nav */}
      <nav style={{ flex: 1, display: 'grid', alignContent: 'start', gap: '0.18rem' }}>
        {NAV.map((item) => (
          <NavButton key={item.id} item={item} active={activePage === item.id} onClick={() => setActivePage(item.id)} />
        ))}
      </nav>

      {/* Footer */}
      <div style={{
        padding: '0.85rem',
        border: `1px solid ${c.line}`,
        borderRadius: '8px',
        background: c.surface2,
        color: c.muted,
        fontSize: '0.73rem',
        lineHeight: 1.5,
      }}>
        {auth?.supabaseEnabled ? (
          <div style={{ marginBottom: '0.75rem' }}>
            {auth.user ? (
              <>
                <div style={{ color: c.ink, fontSize: '0.71rem', marginBottom: '0.45rem', wordBreak: 'break-word' }}>{auth.user.email}</div>
                <button onClick={auth.signOut} style={{
                  background: c.surface, color: c.steel, border: `1px solid ${c.line}`,
                  borderRadius: '6px', padding: '0.42rem 0.6rem', fontSize: '0.71rem',
                  cursor: 'pointer', width: '100%', textAlign: 'left',
                }}>
                  Sign out
                </button>
              </>
            ) : (
              <button onClick={auth.signInWithGoogle} style={{
                background: c.teal, color: '#fff', border: 'none',
                borderRadius: '6px', padding: '0.5rem 0.7rem', fontSize: '0.73rem',
                cursor: 'pointer', width: '100%', textAlign: 'left', fontWeight: 700,
              }}>
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
