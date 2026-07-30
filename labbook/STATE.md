# tx_data — current state

Last updated: 2026-07-30 UTC

## Landed

- Repo scaffolded at `/home/homer/.openclaw/workspace/projects/tx_data` (commit `4cef840`).
  - Package layout under `src/tx_data/` with `paths.py` + `catalog.py`.
  - Config-driven data root: env var `TX_DATA_ROOT` overrides `config/paths.yml`.
  - `nemo_mock/` mirrors `/Nemo` structure; checked into git.
  - `data/` and `db/` gitignored (derived, regenerable).
  - Smoke tests pass (path resolution + catalog dir).
- Agent scaffold: `agent/SYSTEM.md` + `agent/prompts/{conventions,data-model,etl-pattern,privacy,review}.md`.
- Lab book scaffold: this file + `DECISIONS.md` + today's dated log.

## In flight

- Nothing yet. Awaiting first table description from Odysseus.

## Open questions

- Which table to start with? Recommended: `clinical` as the anchor for other tables' relationships.
- Governance: confirm TRACERx data-access-committee expectations before any LLM integration touches real data.
- Ownership: who is first author on the eventual TRACERx-Bench preprint, and when do we raise this with McGranahan / the consortium?

## Next actions

1. Odysseus describes first table: `/Nemo` path, format, columns, semantics, primary key, relationships.
2. Build: synthetic mock file + conversion script + catalog entry + test. One commit.
3. Loop for next table.
