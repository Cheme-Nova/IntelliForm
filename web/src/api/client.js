import axios from 'axios'

const publicApiKey = import.meta.env.VITE_PUBLIC_API_KEY
const configuredApiUrl = import.meta.env.VITE_API_URL?.trim()
const isLocalHost =
  typeof window !== 'undefined' &&
  ['localhost', '127.0.0.1', '::1'].includes(window.location.hostname)
export const apiBaseUrl = configuredApiUrl || (isLocalHost ? 'http://localhost:8000' : 'https://intelliform-4x48.onrender.com')

const client = axios.create({
  baseURL: apiBaseUrl,
  headers: {
    'Content-Type': 'application/json',
    ...(publicApiKey ? { 'X-API-Key': publicApiKey } : {}),
  }
})

client.interceptors.request.use((config) => {
  const token = localStorage.getItem('intelliform_access_token')
  if (token) {
    config.headers.Authorization = `Bearer ${token}`
  }
  return config
})

/**
 * Stream a formulation run via SSE.
 * Calls onEvent(event) for each server-sent event.
 * Returns a cleanup function that aborts the stream.
 */
export function streamFormulate(data, onEvent) {
  const controller = new AbortController()
  const token = localStorage.getItem('intelliform_access_token')

  const headers = {
    'Content-Type': 'application/json',
    Accept: 'text/event-stream',
    ...(publicApiKey ? { 'X-API-Key': publicApiKey } : {}),
    ...(token ? { Authorization: `Bearer ${token}` } : {}),
  }

  fetch(`${apiBaseUrl}/api/v1/stream/formulate`, {
    method: 'POST',
    headers,
    body: JSON.stringify(data),
    signal: controller.signal,
  })
    .then(async (res) => {
      if (!res.ok) {
        const text = await res.text().catch(() => res.statusText)
        onEvent({ step: 'error', message: text })
        return
      }
      const reader = res.body.getReader()
      const decoder = new TextDecoder()
      let buffer = ''
      while (true) {
        const { done, value } = await reader.read()
        if (done) break
        buffer += decoder.decode(value, { stream: true })
        const parts = buffer.split('\n\n')
        buffer = parts.pop() ?? ''
        for (const part of parts) {
          const line = part.trim()
          if (line.startsWith('data: ')) {
            try {
              onEvent(JSON.parse(line.slice(6)))
            } catch {
              // ignore malformed lines
            }
          }
        }
      }
    })
    .catch((err) => {
      if (err.name !== 'AbortError') {
        onEvent({ step: 'error', message: err.message })
      }
    })

  return () => controller.abort()
}

export const api = {
  health: () => client.get('/health'),
  verticals: () => client.get('/api/v1/verticals'),
  failureTypes: () => client.get('/api/v1/failure-types'),
  memory: (n = 10) => client.get(`/api/v1/memory?n=${n}`),
  me: () => client.get('/api/v1/me'),
  projects: (limit = 25) => client.get(`/api/v1/projects?limit=${limit}`),
  formulate: (data) => client.post('/api/v1/formulate', data),
  pareto: (data) => client.post('/api/v1/optimize/pareto', data),
  bayesian: (data) => client.post('/api/v1/optimize/bayesian', data),
  qsar: (data) => client.post('/api/v1/predict/qsar', data),
  reformulate: (data) => client.post('/api/v1/reformulate', data),
  refine: (data) => client.post('/api/v1/refine', data),
  pubchemEnrich: (names) => client.post('/api/v1/pubchem/enrich', { names }),
  pharmaDeepDive: (data) => client.post('/api/v1/pharma/deep-dive', data),
  carbonPassport: (data) => client.post('/api/v1/export/carbon-passport', data),
}

export default client
