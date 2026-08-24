# TAG modularization + prompt-cache measurement — 2026-07-28

This is an operational measurement, not a KB-answer benchmark. It does not score
the system against answers generated from the same DuckDB.

## Module boundary

The former 445-line `tag_query.py` owned syn + exec + gen + both LLM backends.
It is now a 206-line CLI/compatibility surface over:

- `tag_syn.py` — deterministic schema rendering, schema SHA-256, few-shot cache
  prefix, NL-to-SQL, and SQL repair;
- `tag_exec.py` — SELECT/WITH guard and read-only DuckDB execution;
- `tag_gen.py` — authoritative-chain context, BM25 excerpts, row formatting,
  and grounded answer generation;
- `tag_backend.py` — Anthropic Messages API prompt cache and the explicit
  uncached Claude CLI fallback.

`tag_query.syn`, `exec_sql`, `content_search`, `rows_to_text`, `gen`, `llm`, and
the existing CLI flags remain available. `--cache-ttl {1h,5m}` is additive.

## Cache shape and invalidation

The SDK request places a `cache_control` breakpoint after:

```
SYN_SYS + SCHEMA_SHA256 + schema_text + FEWSHOT
```

The question and optional SQL-repair error are a second, uncached content block.
Categorical values are now sorted; before this change, repeated renders of the
same DuckDB could permute `SELECT DISTINCT` values and defeat exact-prefix cache
matching. Three independent connections after the fix produced the same
6,623-character schema and the same SHA-256
`3a9b899c08617a1305d5022a5432ff35d02f8dadececff52a4c39119866e6037`
before the final refresh. A Notion refresh that changes the rendered schema
changes both the visible hash and the exact prefix bytes.

The default TTL is 1 hour. This subsystem is a one-question, human-triggered CLI;
follow-up questions commonly arrive outside five minutes, while it is not a
high-throughput service refreshing the prefix several times per five-minute
window. Anthropic's prompt-cache guidance recommends 1 hour for that cadence;
`--cache-ttl 5m` remains available for bursty scripted use.

## Actual tokens and latency

Test question:

> How many SourceEvidence rows are in the knowledge base?

Both fallback measurements used Claude Opus 5, the same 5,788-character rendered
schema, 1,703-character few-shot block, 682-character syn system prompt, and
produced the same 34-output-token SQL. Values below are the Claude CLI's provider
usage fields plus `time.perf_counter()` wall time.

| path | input | cache creation | cache read | output | wall time |
|---|---:|---:|---:|---:|---:|
| before: legacy `claude -p` from repo cwd, default agent context | 2 | 34,170 | 15,268 | 34 | 3.761 s |
| after: tool-free, explicit system, empty temporary cwd | 2 | 3,901 | 0 | 34 | 3.983 s |

The fallback's provider-accounted input categories fell, but latency increased
by 0.222 s (5.9%). Therefore this measurement does **not** support a “faster”
claim. The large token change is from removing Claude Code's repo/agent context;
it is not an application CAG hit. The CLI may report its own internal cache
categories, but `tag_query` cannot place or verify a `cache_control` breakpoint
through `claude -p`.

`ANTHROPIC_API_KEY` was unset on this host. Consequently an actual SDK cold
write/warm read pair — the measurement that would show
`cache_creation_input_tokens` followed by `cache_read_input_tokens` for this
specific schema prefix — could not be made. No SDK latency or token saving is
claimed. With a key present, every SDK call logs provider input/create/read/output
tokens, schema-hash prefix, TTL, and wall time so the pair can be measured
externally without changing code.

## Verifier duplication measurement

Method: pairwise `SequenceMatcher` over stripped source lines for the six named
verifiers; count only exact contiguous blocks of at least five lines containing
at least four nonblank lines.

Before:

| verifier | lines in any measured clone block |
|---|---:|
| `verify_sources.py` | 5 / 292 (1.7%) |
| `verify_params.py` | 79 / 363 (21.8%) |
| `verify_runs.py` | 75 / 296 (25.3%) |
| `verify_forces.py` | 0 / 553 |
| `verify_gate_contracts.py` | 0 / 478 |
| `verify_gate_contracts_negcontrol.py` | 0 / 237 |

Only `verify_params.py` ↔ `verify_runs.py` had substantive duplication: nine
blocks / 75 matched lines, longest block 12 lines. Their shared YAML manifest
loading, declared-vs-observed drift classification, verdict counting, optional
DuckDB table replacement, and coverage-ratchet dispatch moved to
`verify_manifest_common.py`. The other four verifiers stayed structurally
independent. Re-running the same raw detector leaves six blocks / 54 lines
between params and runs; those are mostly the shared imports/helper call sites
and their domain-specific CLI/report presentation, not duplicated audit
classification or persistence implementations.

The inspection also exposed a separate regression: the gate-contract negative
control still assumed the pre-`run-record@2` artifact shape. It now unwraps only
the `measurements` payload and its 8/8 controls pass; no gate contract or
threshold was changed.

## Generated DuckDB decision

`kb.duckdb` is 44 MB locally, but `git ls-files` returns no entry and
`.gitignore` explicitly excludes it. That is the correct policy: Notion remains
the source of truth, while `refresh.sh` deterministically rebuilds the local
mirror and the lightweight Notion snapshot is versioned.
