import axios from 'axios'

const publicApiKey = import.meta.env.VITE_PUBLIC_API_KEY

const client = axios.create({
  baseURL: import.meta.env.VITE_API_URL || 'http://localhost:8000',
  headers: {
    'Content-Type': 'application/json',
    ...(publicApiKey ? { 'X-API-Key': publicApiKey } : {}),
  }
})

export const api = {
  health: () => client.get('/health'),
  verticals: () => client.get('/api/v1/verticals'),
  failureTypes: () => client.get('/api/v1/failure-types'),
  memory: (n = 10) => client.get(`/api/v1/memory?n=${n}`),
  formulate: (data) => client.post('/api/v1/formulate', data),
  pareto: (data) => client.post('/api/v1/optimize/pareto', data),
  bayesian: (data) => client.post('/api/v1/optimize/bayesian', data),
  qsar: (data) => client.post('/api/v1/predict/qsar', data),
  reformulate: (data) => client.post('/api/v1/reformulate', data),
}

export default client
