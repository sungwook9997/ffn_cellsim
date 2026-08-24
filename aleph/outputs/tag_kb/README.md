# tag_kb — TAG (Table-Augmented Generation) over the ffn_cellsim KB

> **TAG = Table-Augmented Generation** (Biswal et al. 2024, arXiv:2408.14717),
> *not* "태그/label". Answers NL questions over the Notion Contract-Graph + the
> reference PDF corpus by combining exact SQL (joins/aggregations the vector RAG
> can't do) with LLM semantic reasoning over evidence text.

Design rationale + decisions: [`TAG_DESIGN.md`](TAG_DESIGN.md).

## Pipeline (canonical TAG loop)

```
 NL question ── syn ──► DuckDB SQL ── exec ──► rows ─┐
                (Claude, given live schema +          ├─ gen ─► answer + citations
                 actual edges relation vocab)         │  (Claude)
 references/*.pdf ── BM25/FTS over paper_chunks ──────┘
```

- **syn** — `tag_query.py` gives Claude the live `kb.duckdb` schema **plus the
  real relation vocabulary** (`SELECT DISTINCT src_type, rel, dst_type FROM
  edges`) so generated joins are accurate. `tag_syn.py` owns this stage and
  returns one DuckDB `SELECT`.
- **exec** — read-only, `SELECT`/`WITH`-only guard, single self-repair pass that
  feeds the DuckDB error back to `syn`; implemented by `tag_exec.py`.
- **gen** — Claude composes the answer from the SQL rows (structured ground
  truth) + BM25 excerpts from `paper_chunks` (evidence text), citing
  `citation_key` / `KB-x` / `MC-*` / `VG-*` / DOI / `[citation_key pN]`;
  implemented by `tag_gen.py`.

## Files

| file | role |
|---|---|
| `notion_to_duckdb.py` | materialize Notion 9-DB Contract-Graph → `kb.duckdb` (9 node tables + `edges` relation triples). Reuses `../obsidian_rag_full/notion_to_obsidian.py` plumbing (same token, same DB ids). |
| `references_ingest.py` | `references/` PDFs → `paper_refs` + `paper_chunks` (page-chunked, fitz) + DuckDB BM25 FTS. Links to SourceEvidence by DOI; skips exact byte-dup PDFs. Writes `references/tag_corpus.json`. |
| `tag_query.py` | stable CLI orchestration for syn→exec→gen. |
| `tag_syn.py` / `tag_exec.py` / `tag_gen.py` | stage implementations with one concept per module. |
| `tag_backend.py` | Anthropic SDK + prompt-cache accounting, or explicit uncached `claude -p` fallback. |
| `refresh.sh` | rebuild `kb.duckdb` from Notion + PDFs (Notion = SoT). |
| `kb.duckdb` | generated table+content layer (gitignored — regenerable). |
| `tag_eval.py` | regression guard: gold NL questions → assert `tag_query.py` returns the right substrings. |
| `verify_sources*.py` + `AUDIT_FINDINGS.md` | citation-integrity audit (CrossRef→web). Verdicts in `source_audit` table; 3 confirmed hallucinations flagged. |
| `kb_benchmark.py` + `kb_benchmark_vis.py` + `BENCHMARK_REPORT.md` | KB-architecture ablation benchmark (below). |

## KB-architecture benchmark (`kb_benchmark.py`)

This is a retained diagnostic harness, not verification authority. Its historic
C4 score is non-quotable: the gold answers and C4 context came from the same SQL
over the same DuckDB, and the old run is no longer reproducible (see
`STATE_NONQUOTABLE.md` (c)11). Do not use it to claim that the architecture
improved accuracy.

Methodological crux: `claude -p` is a full agent — left default it would *read*
`kb.duckdb` even under "closed-book", invalidating the ablation. The harness
forces a true single-turn completion (`--tools ""`, empty scratch cwd,
tool-call artifacts stripped), so C1 only has parametric memory.

```bash
python kb_benchmark.py                       # full run (answerer=sonnet, judge=opus)
python kb_benchmark.py --n 3 --no-judge      # quick programmatic-only smoke test
python kb_benchmark_vis.py                   # render figs/ from benchmark_results.json
```

## Data model (`kb.duckdb`)

- 9 node tables: `source_evidence`, `knowledge_claim`, `model_contract`,
  `parameter`, `validation_gate`, `code_mapping`, `run_result`,
  `decision_ledger`, `cell_state` (all columns TEXT; `CAST` in SQL for numerics).
- `edges(src_id, src_type, rel, dst_id, dst_type)` — relation triples; traverse
  by joining a node's `id` to `edges.src_id` / `edges.dst_id`.
- `paper_refs(citation_key, title, doi, path, n_pages, sha1, source, se_uid,
  se_citation_key, n_chunks)` + `paper_chunks(chunk_id, citation_key, se_uid,
  page, text)` + BM25 FTS index.

## Usage

```bash
conda activate ffn_sim
bash refresh.sh                                    # build/rebuild kb.duckdb
python tag_query.py "How many SourceEvidence papers have no DOI?"
python tag_query.py "Which 5 papers support the most KnowledgeClaims?"
python tag_query.py --sql-only "..."               # show SQL only
```

## LLM backend and prompt cache

The default model is `claude-opus-5`. With `ANTHROPIC_API_KEY`, the syn request
is split into a stable cached prefix and a dynamic question:

```
SYN_SYS + SCHEMA_SHA256 + schema_text + FEWSHOT  | cache breakpoint
question + optional SQL-repair error              | uncached suffix
```

The SHA-256 is computed from the complete rendered schema, including categorical
values and relation vocabulary. A refresh that changes that rendering therefore
changes both the visible key and the exact cached bytes. The backend logs
provider `input_tokens`, `cache_creation_input_tokens`,
`cache_read_input_tokens`, output tokens, and elapsed time.

The default TTL is **1 hour**. This is a human-triggered, one-question CLI, not a
service issuing the same prefix more often than every five minutes; Anthropic
recommends the 1-hour TTL for follow-ups likely to arrive after five minutes and
before one hour. Use `--cache-ttl 5m` for bursty scripted use.

Without `ANTHROPIC_API_KEY`, `claude -p` cannot send a `cache_control`
breakpoint. The fallback emits a visible `prompt-cache BYPASS` warning with the
schema hash and resent character count, disables tools, supplies a minimal
system prompt, and runs from an empty temporary directory. It never silently
pretends that the prefix was cached.

## Refresh policy

Notion is the source of truth. `kb.duckdb` is a regenerable mirror — re-run
`refresh.sh` after Notion edits or when new PDFs land in `references/`. The
44-MB local mirror is deliberately gitignored; code and the durable Notion
snapshot are tracked instead.
