-- IntelliForm v0.9 — Supabase Database Migration
-- Run this in your Supabase SQL editor: app.supabase.com → SQL Editor → New Query
-- migrations/001_create_tables.sql

CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- ── Projects table ────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS intelliform_projects (
    id          UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    session_id  TEXT NOT NULL,
    user_email  TEXT,
    created_at  TIMESTAMPTZ DEFAULT NOW(),
    application TEXT,
    blend       JSONB,
    cost_per_kg FLOAT,
    bio_pct     FLOAT,
    perf_score  FLOAT,
    eco_score   FLOAT,
    eco_grade   TEXT,
    optimizer   TEXT,
    parser      TEXT,
    relaxed     BOOLEAN DEFAULT FALSE,
    nl_input    TEXT
);

-- ── QSAR active learning feedback ─────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS intelliform_feedback (
    id          UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    session_id  TEXT,
    smiles      TEXT,
    target      TEXT CHECK (target IN ('Biodegradability','Ecotoxicity','Performance')),
    predicted   FLOAT,
    actual      FLOAT,
    created_at  TIMESTAMPTZ DEFAULT NOW()
);

-- ── Pilot bookings ─────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS intelliform_pilot_bookings (
    id          UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    session_id  TEXT,
    user_email  TEXT,
    blend       JSONB,
    cost_per_kg FLOAT,
    batch_kg    INT DEFAULT 500,
    quote_usd   FLOAT,
    application TEXT,
    created_at  TIMESTAMPTZ DEFAULT NOW(),
    status      TEXT DEFAULT 'pending' CHECK (status IN ('pending','confirmed','completed','cancelled'))
);

CREATE TABLE IF NOT EXISTS intelliform_user_usage (
    id          UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id     TEXT NOT NULL,
    user_email  TEXT,
    action      TEXT NOT NULL DEFAULT 'formulate',
    created_at  TIMESTAMPTZ DEFAULT NOW()
);

-- ── Indexes ────────────────────────────────────────────────────────────────────
CREATE INDEX IF NOT EXISTS idx_projects_session   ON intelliform_projects(session_id);
CREATE INDEX IF NOT EXISTS idx_projects_created   ON intelliform_projects(created_at DESC);
CREATE INDEX IF NOT EXISTS idx_feedback_session   ON intelliform_feedback(session_id);
CREATE INDEX IF NOT EXISTS idx_bookings_session   ON intelliform_pilot_bookings(session_id);
CREATE INDEX IF NOT EXISTS idx_bookings_status    ON intelliform_pilot_bookings(status);
CREATE INDEX IF NOT EXISTS idx_usage_user         ON intelliform_user_usage(user_id);
CREATE INDEX IF NOT EXISTS idx_usage_created      ON intelliform_user_usage(created_at DESC);

-- ── Row-level security (enable in Supabase dashboard after running this) ───────
-- ALTER TABLE intelliform_projects      ENABLE ROW LEVEL SECURITY;
-- ALTER TABLE intelliform_feedback      ENABLE ROW LEVEL SECURITY;
-- ALTER TABLE intelliform_pilot_bookings ENABLE ROW LEVEL SECURITY;
