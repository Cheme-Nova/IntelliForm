import { useEffect, useState } from 'react'
import { QueryClient, QueryClientProvider } from '@tanstack/react-query'
import Sidebar from './components/Sidebar'
import Formulate from './pages/Formulate'
import EcoMetrics from './pages/EcoMetrics'
import Certifications from './pages/Certifications'
import Stability from './pages/Stability'
import Carbon from './pages/Carbon'
import Regulatory from './pages/Regulatory'
import Pareto from './pages/Pareto'
import QSAR from './pages/QSAR'
import Reformulation from './pages/Reformulation'
import Memory from './pages/Memory'
import { useAuth } from './auth/AuthContext'

const queryClient = new QueryClient()

const PAGES = {
  formulate: Formulate,
  eco: EcoMetrics,
  certifications: Certifications,
  stability: Stability,
  carbon: Carbon,
  regulatory: Regulatory,
  pareto: Pareto,
  qsar: QSAR,
  reformulation: Reformulation,
  memory: Memory,
}

export default function App() {
  const auth = useAuth()
  const [activePage, setActivePage] = useState('formulate')
  const [isMobile, setIsMobile] = useState(() => {
    if (typeof window === 'undefined') return false
    return window.innerWidth < 960
  })
  const PageComponent = PAGES[activePage] || Formulate
  const publicMode = import.meta.env.VITE_PUBLIC_MODE !== '0'

  useEffect(() => {
    function onResize() {
      setIsMobile(window.innerWidth < 960)
    }
    window.addEventListener('resize', onResize)
    return () => window.removeEventListener('resize', onResize)
  }, [])

  return (
    <QueryClientProvider client={queryClient}>
      <div style={{
        display: 'flex',
        flexDirection: isMobile ? 'column' : 'row',
        minHeight: '100vh',
        background: '#F7F8FA',
        color: '#111827',
        fontFamily: 'Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }}>
        <Sidebar activePage={activePage} setActivePage={setActivePage} auth={auth} isMobile={isMobile} />
        <main style={{
          flex: 1,
          minWidth: 0,
          width: '100%',
          padding: isMobile ? '0.85rem' : '2rem',
          overflowY: 'auto',
          overflowX: 'hidden',
          boxSizing: 'border-box',
        }}>
          {publicMode ? (
            <div style={{
              marginBottom: '1rem',
              padding: isMobile ? '0.8rem 0.9rem' : '0.85rem 1rem',
              borderRadius: '8px',
              border: '1px solid #D5E8E2',
              background: '#F1FAF7',
              color: '#275E52',
              fontSize: isMobile ? '0.82rem' : '0.9rem',
              lineHeight: 1.6,
              textAlign: 'left',
              overflowWrap: 'anywhere',
            }}>
              Public edition: free-tier usage limits are enabled to keep IntelliForm broadly accessible. For internal, higher-volume, or enterprise use, keep the Streamlit lab app and the enterprise stack separate.
            </div>
          ) : null}
          <PageComponent />
        </main>
      </div>
    </QueryClientProvider>
  )
}
