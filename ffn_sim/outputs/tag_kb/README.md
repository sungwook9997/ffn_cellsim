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
  edges`) so generated joins are accurate. Returns one DuckDB `SELECT`.
- **exec** — read-only, `SELECT`/`WITH`-only guard, single self-repair pass that
  feeds the DuckDB error back to `syn`.
- **gen** — Claude composes the answer from the SQL rows (structured ground
  truth) + BM25 excerpts from `paper_chunks` (evidence text), citing
  `citation_key` / `KB-x` / `MC-*` / `VG-*` / DOI / `[citation_key pN]`.

## Files

| file | role |
|---|---|
| `notion_to_duckdb.py` | materialize Notion 8-DB Contract-Graph → `kb.duckdb` (8 node tables + `edges` relation triples). Reuses `../obsidian_rag_full/notion_to_obsidian.py` plumbing (same token, same DB ids). |
| `references_ingest.py` | `references/` PDFs → `paper_refs` + `paper_chunks` (page-chunked, fitz) + DuckDB BM25 FTS. Links to SourceEvidence by DOI; skips exact byte-dup PDFs. Writes `references/tag_corpus.json`. |
| `tag_query.py` | the syn→exec→gen engine (CLI). |
| `refresh.sh` | rebuild `kb.duckdb` from Notion + PDFs (Notion = SoT). |
| `kb.duckdb` | generated table+content layer (gitignored — regenerable). |

## Data model (`kb.duckdb`)

- 8 node tables: `source_evidence`, `knowledge_claim`, `model_contract`,
  `parameter`, `validation_gate`, `code_mapping`, `run_result`,
  `decision_ledger` (all columns TEXT; `CAST` in SQL for numerics).
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

## LLM backend (no new secret required)

`tag_query.py` uses the **anthropic SDK** if `ANTHROPIC_API_KEY` is set, else
falls back to **`claude -p`** (Claude Code headless). Override the model with
`--model <id>`.

## Refresh policy

Notion is the source of truth. `kb.duckdb` is a regenerable mirror — re-run
`refresh.sh` after Notion edits or when new PDFs land in `references/`.
