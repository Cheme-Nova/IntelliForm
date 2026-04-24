const LAST_RUN_KEY = 'intelliform_last_run'

export function saveLastRun(run) {
  if (typeof window === 'undefined') return
  try {
    window.localStorage.setItem(LAST_RUN_KEY, JSON.stringify({
      ...run,
      savedAt: new Date().toISOString(),
    }))
  } catch (error) {
    console.warn('Unable to save IntelliForm session', error)
  }
}

export function loadLastRun() {
  if (typeof window === 'undefined') return null
  try {
    const raw = window.localStorage.getItem(LAST_RUN_KEY)
    return raw ? JSON.parse(raw) : null
  } catch (error) {
    console.warn('Unable to load IntelliForm session', error)
    return null
  }
}

export function getLastRunResult() {
  return loadLastRun()?.response ?? null
}

