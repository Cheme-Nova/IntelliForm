import {
  BarChart3,
  Beaker,
  Brain,
  CheckCircle2,
  ClipboardCheck,
  FileText,
  FlaskConical,
  Gauge,
  GitBranch,
  History,
  Leaf,
  Microscope,
  RefreshCcw,
  Sparkles,
  Store,
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
  { id: 'sds',            label: 'Safety Sheet',   icon: FileText      },
  { id: 'suppliers',      label: 'Suppliers',      icon: Store         },
  { id: 'history',        label: 'History',        icon: History       },
]

/* Grouped so there's a visual break before the secondary tools */
const NAV_PRIMARY   = NAV.slice(0, 6)
const NAV_SECONDARY = NAV.slice(6)

const c = {
  ink:        '#0E1C2A',
  muted:      '#6B7A8A',
  steel:      '#4E6070',
  line:       '#E2D9CF',
  surface:    '#FFFFFF',
  surface2:   '#F4F1EB',
  teal:       '#1A7C6E',
  tealSoft:   '#EBF4F1',
  tealBorder: '#B5D0C9',
}

const sectionLabel = {
  color: c.muted,
  fontSize: '0.6rem',
  fontWeight: 700,
  textTransform: 'uppercase',
  letterSpacing: '0.14em',
  padding: '0 0.4rem',
  marginBottom: '3px',
  marginTop: '10px',
  display: 'block',
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
        gap: '0.5rem',
        padding: isMobile ? '0.5rem 0.7rem' : '0.52rem 0.7rem',
        textAlign: 'left',
        background: active ? c.tealSoft : 'transparent',
        color: active ? c.teal : c.steel,
        border: `1px solid ${active ? c.tealBorder : 'transparent'}`,
        borderRadius: '7px',
        cursor: 'pointer',
        fontSize: '0.8rem',
        fontWeight: active ? 600 : 450,
        whiteSpace: 'nowrap',
        minHeight: '36px',
        minWidth: 0,
        letterSpacing: 0,
        fontFamily: "'DM Sans', sans-serif",
        transition: 'background 120ms ease, color 120ms ease, border-color 120ms ease',
      }}
      onMouseEnter={e => {
        if (!active) {
          e.currentTarget.style.background = '#F4F1EB'
          e.currentTarget.style.color = c.ink
        }
      }}
      onMouseLeave={e => {
        if (!active) {
          e.currentTarget.style.background = 'transparent'
          e.currentTarget.style.color = c.steel
        }
      }}
    >
      <Icon
        size={14}
        strokeWidth={active ? 2.1 : 1.7}
        style={{ flex: '0 0 auto', opacity: active ? 1 : 0.65 }}
      />
      <span>{item.label}</span>
    </button>
  )
}

/* ── Mobile sidebar ─────────────────────────────────────────────────────── */
function MobileSidebar({ activePage, setActivePage }) {
  return (
    <aside className="if-sidebar-mobile" style={{
      display: 'flex',
      flexDirection: 'column',
      padding: '0.7rem 0 0.55rem',
      position: 'sticky',
      top: 0,
      zIndex: 20,
    }}>
      <div style={{
        padding: '0 0.85rem 0.6rem',
        display: 'flex',
        justifyContent: 'space-between',
        gap: '0.75rem',
        alignItems: 'center',
      }}>
        {/* Brand */}
        <div style={{ display: 'flex', alignItems: 'center', gap: '0.55rem', minWidth: 0 }}>
          <div style={{
            width: '28px', height: '28px', borderRadius: '7px',
            background: c.teal, display: 'grid', placeItems: 'center',
            color: '#fff', flex: '0 0 auto',
            boxShadow: '0 4px 12px rgba(26,124,110,0.3)',
          }}>
            <Sparkles size={13} strokeWidth={2} />
          </div>
          <div style={{ minWidth: 0 }}>
            <div style={{
              fontFamily: "'Fraunces', serif",
              fontOpticalSizing: 'auto',
              color: c.ink,
              fontWeight: 600,
              fontSize: '1rem',
              letterSpacing: '-0.02em',
              lineHeight: 1,
            }}>
              IntelliForm
            </div>
            <div style={{ color: c.muted, fontSize: '0.63rem', marginTop: '1px', letterSpacing: '0.01em' }}>
              ChemeNova
            </div>
          </div>
        </div>

        <div style={{
          display: 'inline-flex', alignItems: 'center', gap: '0.3rem',
          color: c.teal, background: c.tealSoft, border: `1px solid ${c.tealBorder}`,
          borderRadius: '999px', padding: '0.28rem 0.6rem',
          fontSize: '0.62rem', fontWeight: 700, flex: '0 0 auto',
          letterSpacing: '0.08em', textTransform: 'uppercase',
        }}>
          <Sparkles size={10} /> Live
        </div>
      </div>

      <nav style={{
        display: 'flex', gap: '0.3rem', overflowX: 'auto',
        padding: '0 0.85rem', scrollbarWidth: 'none',
        WebkitOverflowScrolling: 'touch', touchAction: 'pan-x',
      }}>
        {NAV.map(item => (
          <NavButton
            key={item.id} item={item}
            active={activePage === item.id}
            onClick={() => setActivePage(item.id)}
            isMobile
          />
        ))}
      </nav>
    </aside>
  )
}

/* ── Desktop sidebar ────────────────────────────────────────────────────── */
export default function Sidebar({ activePage, setActivePage, auth, isMobile = false }) {
  if (isMobile) {
    return <MobileSidebar activePage={activePage} setActivePage={setActivePage} />
  }

  return (
    <aside className="if-sidebar" style={{
      width: '220px',
      display: 'flex',
      flexDirection: 'column',
      padding: '1.25rem 0.85rem',
      position: 'sticky',
      top: 0,
      minHeight: '100vh',
      boxSizing: 'border-box',
      gap: '0',
    }}>

      {/* ── Brand ──────────────────────────────────────────────────────────── */}
      <div style={{ padding: '0.1rem 0.3rem 1.5rem' }}>
        {/* Logo mark */}
        <div style={{
          width: '36px', height: '36px', borderRadius: '9px',
          background: `linear-gradient(135deg, ${c.teal} 0%, #0F5C52 100%)`,
          display: 'grid', placeItems: 'center',
          color: '#fff', marginBottom: '0.9rem',
          boxShadow: '0 6px 18px rgba(26,124,110,0.30), 0 1px 3px rgba(26,124,110,0.20)',
        }}>
          <Sparkles size={16} strokeWidth={2} />
        </div>

        {/* Wordmark */}
        <div style={{
          fontFamily: "'Fraunces', ui-serif, Georgia, serif",
          fontOpticalSizing: 'auto',
          color: c.ink,
          fontWeight: 600,
          fontSize: '1.12rem',
          letterSpacing: '-0.025em',
          lineHeight: 1,
        }}>
          IntelliForm
        </div>
        <div style={{
          color: c.muted,
          fontSize: '0.67rem',
          marginTop: '4px',
          letterSpacing: '0.01em',
          lineHeight: 1.4,
        }}>
          ChemeNova · Formulation AI
        </div>
      </div>

      {/* ── Primary nav ────────────────────────────────────────────────────── */}
      <nav style={{ flex: 1, display: 'flex', flexDirection: 'column', gap: '1px' }}>
        <span style={sectionLabel}>Workspace</span>
        {NAV_PRIMARY.map(item => (
          <NavButton
            key={item.id} item={item}
            active={activePage === item.id}
            onClick={() => setActivePage(item.id)}
          />
        ))}

        <span style={{ ...sectionLabel, marginTop: '14px' }}>Analysis</span>
        {NAV_SECONDARY.map(item => (
          <NavButton
            key={item.id} item={item}
            active={activePage === item.id}
            onClick={() => setActivePage(item.id)}
          />
        ))}
      </nav>

      {/* ── Auth / footer ──────────────────────────────────────────────────── */}
      <div style={{
        marginTop: '1rem',
        padding: '0.85rem',
        border: `1px solid ${c.line}`,
        borderRadius: '9px',
        background: '#F9F7F4',
      }}>
        {auth?.supabaseEnabled ? (
          <div style={{ marginBottom: '0.65rem' }}>
            {auth.user ? (
              <>
                <div style={{
                  color: c.ink, fontSize: '0.7rem', marginBottom: '0.5rem',
                  wordBreak: 'break-word', fontWeight: 500,
                }}>
                  {auth.user.email}
                </div>
                <button
                  onClick={auth.signOut}
                  style={{
                    background: c.surface, color: c.steel,
                    border: `1px solid ${c.line}`, borderRadius: '6px',
                    padding: '0.4rem 0.65rem', fontSize: '0.7rem',
                    cursor: 'pointer', width: '100%', textAlign: 'left',
                    fontWeight: 500,
                  }}
                >
                  Sign out
                </button>
              </>
            ) : (
              <button
                onClick={auth.signInWithGoogle}
                style={{
                  background: c.teal, color: '#fff', border: 'none',
                  borderRadius: '6px', padding: '0.5rem 0.7rem', fontSize: '0.72rem',
                  cursor: 'pointer', width: '100%', textAlign: 'left', fontWeight: 650,
                  letterSpacing: '0.01em',
                }}
              >
                Sign in with Google
              </button>
            )}
          </div>
        ) : null}

        <div style={{ color: c.muted, fontSize: '0.64rem', letterSpacing: '0.02em' }}>
          v2.1.0 · ChemeNova LLC
        </div>
      </div>
    </aside>
  )
}
