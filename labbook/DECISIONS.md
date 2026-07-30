# Architecture decisions

Append-only. Never edit or delete an entry; supersede with a new one.

---

## 2026-07-30-001 — DuckDB + Parquet + YAML catalog, not Postgres

**Context:** Cohort-scale analytical work, lab-scale team, want zero admin burden.

**Decision:** Storage is Parquet files under `data/`, queried by DuckDB attached over them. Semantic layer is a set of YAML files under `catalog/`.

**Consequences:**
- Zero server / ops overhead.
- SQL access without loading data into memory.
- Parquet + DuckDB is well-supported by modern text-to-SQL tooling if we add an LLM later.
- Single-user semantics; if we ever need concurrent multi-user writes, revisit.

---

## 2026-07-30-002 — Path abstraction via `TX_DATA_ROOT` env var, not text-swap

**Context:** Real data at `/Nemo`, dev at `nemo_mock/`. Need to swap without editing code.

**Decision:** All code resolves paths via `tx_data.paths.data_root()`. It reads `TX_DATA_ROOT` env var; falls back to `config/paths.yml` (default `nemo_mock`). No hardcoded `/Nemo` anywhere.

**Consequences:**
- Same code runs locally and on HPC — `export TX_DATA_ROOT=/Nemo` and go.
- No accidental commit of `/Nemo` paths.
- All ETL scripts must use the helper; enforced by review checklist.

---

## 2026-07-30-003 — Custom YAML catalog now; migrate to LinkML later if we publish

**Context:** LinkML is the closest thing to a bio-community standard for schema modelling. Has a learning curve.

**Decision:** Start with a custom YAML format documented in `agent/prompts/data-model.md`. If the project reaches preprint stage, migrate.

**Consequences:**
- Faster prototyping now.
- Migration cost is real but bounded (our YAML shape is close to LinkML by design).
- We lose LinkML's free validator / doc-generation until migration.

---

## 2026-07-30-004 — Lab book: raw daily notes + curated STATE.md + append-only DECISIONS.md

**Context:** Both the agent and Odysseus need persistent memory across sessions.

**Decision:** Three-file pattern:
- `labbook/YYYY-MM-DD.md` — raw daily notes, written freely.
- `labbook/STATE.md` — curated current state, updated per session.
- `labbook/DECISIONS.md` — append-only architecture decision records.

**Consequences:**
- Mirrors the Homer memory pattern (proven).
- Some content duplication between the daily log and `STATE.md` by design (raw vs curated).
- Requires session-end discipline to update `STATE.md`; enforced by review checklist.
