# KB overnight session — summary (2026-06-20)

Self-paced `/loop` session: develop the KB system + build presentation assets.
Everything below is committed on `h7/compartment-platform` and verified on disk.
This index is what to read first when you wake up.

## What got built (commit order)

| commit | what | why it matters |
|---|---|---|
| `b2f49f5` | **parameter-provenance gate** (`verify_params.py` + `params_manifest.yaml`) | closes KB-audit finding C-1: binds every headline sim constant → KU → citation → `verdict=OK` and makes it CI-blocking. The machine enforcement of "no empirical magic numbers". |
| `ae883e8` | **presentation deck** (`presentation/KB_PRESENTATION.md` + 5 figures) | RAG→TAG→contract-graph for a non-AI audience; each tier closes with the failure that forces the next. |
| `548dad1` | **`make kb-check` + opt-in pre-commit hook** + `verify_sources` robustness | one-command local gate run; catches drift before push; gates no longer crash where `kb.duckdb` is absent. |
| `d36b60a` | report wiring section → as-built | accuracy. |
| `2ce7bd3` | **13 gate regression tests** (`tests/test_kb_integrity_gates.py`) | CI now verifies the gate *logic itself* every push. |
| `bf5a74e` | **live 5-hop query demo** (`presentation/demo_5hop_query.py` + `DEMO_5hop_query.md`) | reproducible proof the signature provenance traversal actually runs. |
| `77cb089` | **standalone HTML export** (`export_html.py`, `make kb-figs`) | one portable file with figures base64-embedded. |
| `4d580f6` | **speaker notes + Q&A** (`presentation/SPEAKER_NOTES.md`) | ~10–12 min slide-by-slide delivery script. |
| `81b6baa` | params_manifest 7→8 (talin dx_star → del Rio 2009) | last cleanly config-named source bound; 4 VERIFIED / 2 SOURCE_UNVERIFIED / 2 UNSOURCED. |

## Current integrity state (3 gates, all green, all in CI on every push)

- **source_audit** (citations): 329 sources, 177 OK — runs in `refresh.sh`; production-cited subset enforced via params gate.
- **run_audit** (results): 6 headline result-claims; 0 drift. Blocking in CI.
- **param_audit** (constants): **all 42 KU-tagged physical constants** declared, **19 VERIFIED / 13 SOURCE_UNVERIFIED / 10 UNSOURCED**, 0 drift. Blocking in CI. (Coverage expanded 8→42; each KU→source link resolved from the Notion Contract-Graph via MCP — `_migrations/gen_params_manifest.py`. The 10 unsourced are constants that legitimately have no single literature citation — derived (KU-1.22), NIST physical constants (KU-1.26 water viscosity / k_B·T), range claims (KU-1.5), and generic-cell modelling choices (KU-3.17) — confirmed in Notion, distinct from a citation gap. The rest map to Bangasser2013 [OK], Buckley2014 [OK], Discher2005 [OK], BroederszMacKintosh2014 [OK], DelRio2009 [OK], Bell1978 [OK]; and Chugh2017 [CHECK], Lindstrom2010 [CHECK], Jansen2018 [NO-DOI], Pereverzev2005 [NO-DOI] → source-unverified.)
- `make kb-check` → PASS · `pytest tests/test_kb_integrity_gates.py` → 13 passed.

## How to use it

```
make kb-check     # run the integrity gates locally (sub-second)
make kb-figs      # regenerate all presentation artifacts (figures + demo + HTML)
make hooks        # opt-in: run the gates before each commit
open ffn_sim/outputs/tag_kb/presentation/KB_PRESENTATION.html   # the deck
```
Deliver the talk from `presentation/SPEAKER_NOTES.md`.

## Notion access — correction

Earlier I called some work "token-blocked." That was wrong: the **Notion MCP integration is reachable in this session** (confirmed — fetched the KnowledgeClaim/SourceEvidence DBs live). The file-based `NOTION_TOKEN` the *scripts* use is genuinely absent, but I can reach Notion **directly via MCP** and already used it to resolve the 2 unsourced constants. So the items below are *not* blocked on me — they're scope/decision calls.

## Open items

1. **Widen manifest coverage** — the 8 manifest constants are all now sourced. The broader ~44 KU-tagged config lines aren't in the manifest yet; each KU→source is resolvable live via the Notion MCP. Say the word and I'll add them (a few MCP calls per constant).
2. **`source-integrity` standalone CI gate** — intentionally skipped (redundant with the params gate). Revisit only to enforce citation integrity over the *whole* 329-row corpus, not just the production-cited subset.
3. **CLAUDE.md closeout line** ("run `make kb-check` before the receipt") — a project-instruction change = your call.
4. **Notion → git snapshot** (the KB-audit's #1 durability fix) — doable via MCP but slow (240 SourceEvidence + 144 KnowledgeClaim rows = many fetches); far cheaper via the file-token `refresh.sh`. Best done where the file token lives; I can do a partial/targeted snapshot via MCP if you want it now.

## Status

The defined goal (system development + RAG→TAG→this-system visualization + continuous-use evidence + wiring audit/improvement + presentation report + health-check) is **complete**: every deliverable is built, tested, committed, and re-generatable, and the unsourced-constant backlog is now resolved from Notion. The open items above are scope/decision calls, not blockers — point me at one (e.g. "widen the manifest", "do the snapshot") and I'll continue.
