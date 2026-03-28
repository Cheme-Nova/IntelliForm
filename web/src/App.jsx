import { useState } from 'react'
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
  const [activePage, setActivePage] = useState('formulate')
  const PageComponent = PAGES[activePage] || Formulate

  return (
    <QueryClientProvider client={queryClient}>
      <div style={{ display: 'flex', minHeight: '100vh', background: '#0A1628', color: '#fff', fontFamily: 'Inter, sans-serif' }}>
        <Sidebar activePage={activePage} setActivePage={setActivePage} />
        <main style={{ flex: 1, padding: '2rem', overflowY: 'auto' }}>
          <PageComponent />
        </main>
      </div>
    </QueryClientProvider>
  )
}
