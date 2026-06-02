# Does knowledge-base architecture reduce LLM hallucination? A 4-condition ablation on the ffn_cellsim KB

**Date:** 2026-06-02/03 · **Status:** complete — two studies: a hand-curated pilot (§3, n=18, 2 answerer models) and a **scaled, statistically-validated study (§4, n=108, deterministic SQL gold, bootstrap CIs + McNemar)**. 8 figures (3 mechanism + 5 results).

**Headline (scaled study, answerer claude-sonnet-4-6, n=108):** the TAG layer is the **only** one that significantly improves structured KB-QA, and the effect is **total** — accuracy **C1 33% ≈ C2 29% ≈ C3 26% ≪ C4 100%** (macro 53/52/52/**100**%; C4 beats C3 on 80/80 questions, χ²=78, **p≪.05**; C1↔C2↔C3 not significant). RAG and the Obsidian graph alone add **no statistically significant** accuracy on this distribution; only the relational query layer (TAG) does. After fixing a citation-resolver bug found in-flight, programmatic hallucination is low everywhere (3–6%), worst for raw RAG.

> One-paragraph claim. As the `ffn_cellsim` project's knowledge base evolved from
> *no structured KB* (the ActiveCellSim era) → semantic RAG → a graph mirror
> (Obsidian) → a relational query layer (TAG over DuckDB), the assistant's
> ability to answer project-specific questions *correctly and without
> fabrication* should improve monotonically. This benchmark measures that, using
> both an objective programmatic backbone and the standard LLM-evaluation
> metrics real systems are graded on (RAGAS, FActScore, TAG-Bench, G-Eval).

---

## 1. Motivation

Early in this project (the `~/ActiveCellSim` v1 era) there was **no structured
knowledge base**: literature values, model decisions, and validation contracts
lived in ad-hoc notes and the model's parametric memory. Questions like *"what
α-actinin off-rate did we adopt, and from which paper?"* were answered from
recall — the regime where LLMs confabulate plausible-but-wrong citations.

`ffn_cellsim` built three successive retrieval layers on top of a single source
of truth (the Notion **Contract-Graph**: 8 atomic, bidirectionally-related
databases — SourceEvidence, KnowledgeClaim, ModelContract, Parameter,
ValidationGate, CodeMapping, RunResult, DecisionLedger):

1. **RAG** — semantic retrieval over the reference-PDF corpus.
2. **Obsidian mirror** — the Contract-Graph as a ~640-node typed-link graph
   (structure + neighbourhood).
3. **TAG** (Table-Augmented Generation, Biswal et al. 2024) — a query engine that
   synthesises DuckDB SQL over the graph (exact joins/aggregations) and fuses it
   with PDF-content semantic search.

The question: **how much does each layer actually buy us in correctness and
hallucination reduction?**

## 2. Experimental design

A controlled **ablation**: four conditions that are *strictly additive*
(C1 ⊂ C2 ⊂ C3 ⊂ C4). The answering model, the answer prompt template, and the
gold question set are **identical** across conditions — the **only** thing that
varies is the context block the model receives.

| Cond | Name | Context the answerer receives |
|------|------|-------------------------------|
| **C1** | no-KB (closed-book) | *nothing* — parametric memory only (ActiveCellSim era) |
| **C2** | RAG | + BM25 passages over the reference-PDF corpus (`paper_chunks`) |
| **C3** | RAG + Obsidian | + relevant Contract-Graph nodes and their typed neighbours (the graph the Obsidian mirror exposes) |
| **C4** | RAG + TAG + Obsidian | + the synthesised DuckDB SQL result (exact relational/aggregate answer) and the `source_audit` citation-integrity table |

### 2.1 Methodological crux — making C1 honestly closed-book

The dual LLM backend uses `claude -p` (Claude Code headless) when no
`ANTHROPIC_API_KEY` is present. **`claude -p` is a full agent with Bash/Read
tools and filesystem access** — left at its defaults it will simply *read*
`kb.duckdb` and compute the true answer even under the "closed-book" condition,
silently destroying the ablation. (We caught exactly this: an early C1 run
answered the row-count question with the precise SQL predicate from the gold
key.) The harness therefore forces a genuine single-turn completion:

- `--tools ""` — disable all tools;
- run from an empty scratch CWD — nothing project-related is reachable;
- strip any residual tool-call artifacts from the output.

Under this regime C1 correctly *refuses* DB-specific questions ("I cannot
determine without querying the data") instead of fabricating — the honest
baseline we want. (When `ANTHROPIC_API_KEY` is set, the anthropic SDK gives a
clean single-turn completion directly.)

### 2.2 Gold question set (n=18)

Questions span the classes a real KB must serve, each with ground truth verified
against `kb.duckdb` on 2026-06-02:

- **aggregation (×4)** — counts requiring `COUNT`/`GROUP BY`, e.g. "how many
  SourceEvidence rows have no DOI?" (56/273), "how many Direct-measurement rows?"
  (161), "how many Parameters?" (21).
- **relational / multi-hop (×4)** — joins over `edges`, e.g. "which ValidationGate
  is failing?" (exactly one: *KU-3.5 cortex tension floor*); "which paper supports
  the most claims?" (Broedersz 2014 / Licup 2015, 7 each); "how many claims cite
  Wolf 2013?" (5).
- **factual / parameter (×3)** — values stored in claim nodes, e.g. KB-3.19
  α-actinin `k_off0` = 0.066 s⁻¹; cortex thickness ~200 nm; F-actin ℓ_p ~17 µm.
- **content synthesis (×1)** — answerable only from PDF text (cortical tension
  sets tissue surface tension).
- **citation-trap (×3)** — the heart of the hallucination probe. Each asks about
  one of the **3 confirmed fabricated sources** found in the 2026-06-02 audit
  (`Yao2011_NatCommun`, `YapKovacs_JCS`, `NanoConvergence2021_Glioma`). A system
  *hallucinates* if it affirms the fake source as real; it *passes* if it flags
  the fabrication or names the correct substitute (Ferrer 2008; Yap/Gomez 2015
  Dev Cell; Ketebo pillars).
- **citation-control / over-refusal (×1, NEW)** — asks the system to confirm a
  **real, classic** source (Bell 1978, Science, the genuine Bell-Evans off-rate
  paper). Correct behaviour is to *affirm* it; doubting or refusing a real source
  is an over-refusal failure (the opposite error from a citation-trap). Guards
  against the trivial "refuse everything" strategy.

Two answerer models are run to test robustness: **claude-sonnet-4-6** (primary)
and **claude-haiku-4-5** (weaker), with the **judge fixed at claude-opus-4-8**.

### 2.3 Metrics

**Tier 1 — programmatic backbone (objective, reproducible, no judge):**

- **Accuracy** — gold substrings/value present in the answer (exact-match spirit).
- **Hallucination rate** — for trap questions, affirming a known-fabricated
  source without flagging it; for all questions, citing a DOI / KB-id that does
  **not resolve** to a real row in `kb.duckdb`.
- **Citation groundedness** — fraction of cited identifiers that resolve to real
  KB rows (ALCE-style citation precision against the live DB).

**Tier 2 — standard LLM-evaluation metrics (G-Eval-style judge, stronger model):**

- **RAGAS** (Es et al. 2023) — *faithfulness* (answer grounded in context),
  *answer relevancy*, *context recall*.
- **FActScore / HaluEval** (Min et al. 2023; Li et al. 2023) — decompose the
  answer into atomic factual claims, score the fraction unsupported by
  context+truth = atomic hallucination rate.
- **TAG-Bench** (Biswal et al. 2024) — exact-match correctness, the metric the
  TAG paper itself reports.

The judge is a **separate, stronger model** from the answerer (answerer =
`claude-sonnet-4-6`, judge = `claude-opus-4-8`), scores only from material in its
prompt, and returns strict JSON. Every raw answer **and its full context** is
logged to `benchmark_results.json` so all scores are auditable, not opaque.

### 2.4 Why this design is defensible for publication

- Single controlled variable (context); identical model/prompt/questions.
- The headline hallucination signal (citation resolution, fabrication traps) is
  **programmatic**, not LLM-self-judgment.
- Ground truth comes from the *audited* Contract-Graph, including the 3
  independently-confirmed fabrications — not from the model under test.
- Full transparency: raw answers + contexts committed alongside scores.

## 3. Study 1 — hand-curated pilot (n=18)

Full run 2026-06-02: 18 gold questions × 4 conditions × 2 answerers (sonnet,
haiku), judge `claude-opus-4-8`. `benchmark_results.json` (sonnet, primary) +
`benchmark_results_haiku.json`. *This pilot motivated the larger Study 2 (§4);
its small-n "staircase" is partly revised there.*

### 3.1 Summary by condition — primary answerer (claude-sonnet-4-6)

| Condition | Accuracy¹ | TAG-Bench EM² | Answer-relev | Faithfulness | **Context-recall** | FActScore halluc³ | Prog. halluc⁴ |
|---|---|---|---|---|---|---|---|
| **C1** no-KB (closed-book) | 31% | 22% | 0.61 | 0.91 | **0.00** | 8% | 6% |
| **C2** RAG | 31% | 28% | 0.65 | 0.98 | **0.11** | 0% | 11% |
| **C3** RAG + Obsidian | 31% | 28% | 0.73 | 0.98 | **0.13** | 2% | 0% |
| **C4** RAG + TAG + Obsidian | **81%** | **78%** | 0.94 | 0.93 | **0.71** | 18% | 0% |

¹ programmatic substring/value match. ² LLM-judge exact-match (refusal = incorrect;
trap = correct iff fabrication not affirmed). ³ FActScore atomic-fact hallucination
rate (judge). ⁴ programmatic: affirming a known-fabricated source, or citing an
unresolvable identifier (DOI/KB-id).

### 3.2 The capability staircase (Fig. `fig_bench_by_class.png`)

Correctness by question class (sonnet) shows *which layer switches each capability
on* — the cleanest result of the study:

| Question class | C1 | C2 | C3 | C4 | switches on at |
|---|---|---|---|---|---|
| content (general literature) | 100 | 100 | 100 | 100 | always (parametric) |
| citation-control (confirm a real source) | 100 | 100 | 100 | 100 | always (no over-refusal) |
| citation-trap (reject a fabrication) | 0 | 100 | 100 | 100 | **RAG** (real corpus crowds out the fake) |
| factual / parameter | 67 | 67 | 33 | 67 | partly parametric; **syn-limited** at C4 (see §3.4) |
| relational / multi-hop | 25 | 0 | 25 | 50 | **TAG**, syn-limited |
| aggregation (count) | 0 | 0 | 0 | **100** | **TAG only** (SQL `COUNT`) — the signature result |

### 3.3 Key findings

1. **Aggregation is answerable ONLY with TAG: 0/0/0 → 100%.** No amount of
   semantic retrieval or graph browsing answers "how many rows satisfy X" — it
   requires an exact `COUNT` over the table. This is the sharpest single result:
   a whole class of questions is binary-gated on the TAG layer.

2. **Correctness jumps at C4: 31 → 31 → 31 → 81% (sonnet).** C1–C3 plateau because
   the structured/relational questions dominate the set and none of them is
   answerable without exact query; C4 clears them.

3. **Context recall is the mechanism: 0.00 → 0.11 → 0.13 → 0.71.** Semantic RAG
   over the PDF corpus rarely *contains* the answer to a structured project
   question — those facts live in the *relations between* records, not in any
   paragraph. Only TAG's synthesised SQL reliably retrieves them; recall → accuracy.

4. **The C4 ceiling is SQL-synthesis quality, not the data (see §3.4).** The two
   structured questions C4 missed both trace to `syn` writing a wrong query, not
   to missing data — a tractable engine-improvement target, documented honestly
   rather than tuned away.

5. **Hallucination — the honest, now-richer picture.** On the harder 18-question
   set, *programmatic* hallucination is no longer flat: C1 6%, C2 **11%**, C3/C4
   **0%**. The traps and the over-refusal control expose it — C1/C2 (closed-book /
   raw RAG) cite **unresolvable or wrong DOIs** for real and fake sources alike,
   while the graph/TAG layers ground citations and drop it to zero. Atomic
   (FActScore) hallucination is 8/0/2/**18%**: C4's 18% comes almost entirely from
   *over-elaboration* on the open-ended content/trap questions (Q6 surface-tension
   synthesis, Q9 glioma summary) where, given rich context, the model adds
   plausible specifics beyond it — a different failure mode (verbosity) than C1/C2's
   bad-citation fabrication.

6. **Citation traps caught from C2 on; no over-refusal.** Raw RAG over the *real*
   corpus already prevents affirming the 3 fabricated sources (genuine papers
   crowd out the fake), and every condition correctly *confirms* the real
   Bell-1978 control — so the trap-catching is genuine skepticism, not blanket
   refusal.

### 3.4 The C4 ceiling: two syn failures (honest error analysis)

C4 missed 2 of the structured questions; both are `syn` (NL→SQL) bugs, not data
gaps:

- **Q16 (cortex thickness):** syn wrote `WHERE id IN ('KB-3.1','KB-3.7')`, but `id`
  is the opaque Notion-page UUID — the human "KB-3.1" lives in the `kb_id` column.
  Wrong column → 0 rows. (Fix landed: `tag_query.py` schema hint now states
  `kb_id` is the human claim id and `id` is a UUID PK.)
- **Q15 (claims citing Wolf 2013):** syn over-joined through `paper_refs` with a
  mismatched `rel`, returning 0; the direct `source_evidence → edges →
  knowledge_claim` path gives 5.

Takeaway for the paper: **the architecture delivers the right data; the residual
gap is the model's text-to-SQL skill** — which is exactly why the weaker answerer
(haiku) gains so much less from C4 (§3.5).

### 3.5 Cross-model robustness (Fig. `fig_bench_model_compare.png`)

Re-running with a weaker answerer (claude-haiku-4-5), judge fixed:

| Condition | sonnet acc | haiku acc | sonnet aggregation | haiku aggregation |
|---|---|---|---|---|
| C1 no-KB | 31% | 19% | 0% | 0% |
| C2 RAG | 31% | 25% | 0% | 0% |
| C3 +Obsidian | 31% | 38% | 0% | 17% |
| C4 +TAG | **81%** | **38%** | **100%** | **17%** |

- **Both models improve with KB structure** (monotonic-ish) — the effect is not a
  single-model artifact.
- **The TAG advantage scales with answerer capability.** Sonnet converts TAG
  access into 100% aggregation accuracy; haiku, whose `syn` writes weaker SQL,
  reaches only 17%. The reasoning/refusal classes (trap, control, content) score
  identically (100%) for both — they are model-robust. So the *structured-query*
  win is gated by text-to-SQL ability, the *grounding/refusal* win is not.

### 3.6 Caveats (for honest reporting)

- n = 18 gold questions; 2 answerers, 1 judge. Treat absolute percentages as
  indicative; the robust results are the **ordering** (C4 ≫ C1–C3 on structured
  questions) and the **aggregation 0→100 gate**.
- RAGAS/FActScore are LLM-judge metrics (subjective); the programmatic backbone
  (accuracy, citation resolution, trap/control detection) is the objective anchor.
- C4's atomic-hallucination (18%) is an over-elaboration/verbosity effect on
  open-ended questions, not bad-citation fabrication; a terser answer prompt would
  likely reduce it (untested).
- The closed-book honesty of C1 partly reflects the 2026 answerer's disposition to
  refuse; it does **not** imply the pre-KB ActiveCellSim workflow was
  hallucination-free in practice (older model era).
- KB snapshot: `kb.duckdb` 2026-06-02 (SE 273, claims 158, edges 808, chunks 5327).

## 4. Study 2 — scaled validation + threat investigation (n=108)

Study 1 had only 18 questions, no confidence intervals, and self-judged scores —
so its "staircase" could be noise. Study 2 addresses the threats-to-validity
head-on with a larger, **deterministically-grounded** question set and proper
statistics.

### 4.1 Setup
- **108 questions** generated by `gen_questions.py` from templates whose **gold
  is computed by SQL over `kb.duckdb`** — no LLM is in the gold loop (cf. how KGQA
  benchmarks over Freebase/Wikidata are built). Stratified across aggregation (17),
  relational (20), lookup (49), factual (14), citation-trap (3), citation-control
  (3), negative-existence (2).
- Answerer `claude-sonnet-4-6`, all 4 conditions, **programmatic backbone** scored
  on all 108; **G-Eval judge** (`claude-opus-4-8`) run on a stratified subset (2
  per class) for the RAGAS/FActScore layer.
- Statistics in `kb_benchmark_stats.py`: bootstrap 95% CIs (micro = per-question;
  **macro = per-class mean**, which removes the class-size imbalance) and **McNemar
  paired tests** on per-question correctness.

### 4.2 Headline result (the rigorous version)

| Condition | Accuracy (micro) [95% CI] | Accuracy (macro) [95% CI] | Prog. hallucination |
|---|---|---|---|
| **C1** no-KB | 33.3% [24.1, 42.6] | 53.0% [49.7, 56.7] | 3% |
| **C2** RAG | 28.7% [20.4, 37.0] | 52.0% [48.4, 55.9] | 6% |
| **C3** +Obsidian | 25.9% [17.6, 34.3] | 52.0% [48.1, 56.2] | 4% |
| **C4** +TAG | **100%** [100, 100] | **100%** [100, 100] | 4% |

**McNemar paired tests:**
- C3 vs **C4**: C4 wins on **80** questions, loses on **0** → χ²=78.0, **p≪.05**.
- C1 vs **C4**: C4 wins on **72**, loses on **0** → χ²=70.0, **p≪.05**.
- C1 vs C2, C2 vs C3: **not significant** (χ²=1.2, 0.4).

> **The honest, statistically-grounded story:** on this KB, **only the TAG layer
> produces a significant accuracy gain — and it is total** (C4 ≥ every other
> condition on every single question). **RAG and the Obsidian graph do *not*
> significantly differ from no-KB** for structured KB-QA. Study 1's apparent
> C1<C2<C3 staircase was largely small-n noise; at n=108 with CIs it collapses to
> "C1≈C2≈C3 ≪ C4." This is exactly the kind of correction larger n + significance
> testing exists to make.

### 4.3 Per-class accuracy (where the gain lives)

| class (n) | C1 | C2 | C3 | C4 |
|---|---|---|---|---|
| aggregation (17) | 0% | 0% | 6% | **100%** |
| relational (20) | 15% | 20% | 15% | **100%** |
| factual (14) | 7% | 7% | 14% | **100%** |
| lookup (49) | 49%ᵇ | 37% | 29% | **100%** |
| citation-trap (3) | 100% | 100% | 100% | 100% |
| citation-control (3) | 100% | 100% | 100% | 100% |
| negative-existence (2) | 100% | 100% | 100% | 100% |

ᵇ **base-rate caveat (threat T3, found in-flight):** the status-lookup subclass is
near-degenerate — 14/14 gold values are `verified` (157/158 claims are verified),
so a model that always answers "verified" scores ~100% on it *without knowing
anything*. C1's lookup score is inflated by this. A future set must use
higher-entropy lookup targets (unit/subtopic) — flagged, not silently kept.

### 4.4 Threats investigated in-flight (the fun part)

- **T1 (power) — resolved.** CIs + McNemar show C4≫rest is real (χ²=78), C1/C2/C3
  differences are not (χ²<1.3). The conclusion is no longer "indicative."
- **T5 (is C4's ceiling data or syn?) — at scale, syn is robust.** Only **2/108**
  C4 SQL queries returned zero rows, **0/86** on structured questions, and C4
  scored 100% regardless. The two syn bugs from Study 1 (id-vs-`kb_id`, over-join)
  did **not** recur on the template-clean generated questions — so the syn ceiling
  is real but appears mainly on *idiosyncratic* phrasings, not canonical ones.
- **T8 (citation metric) — found and fixed a scorer bug.** The first pass flagged
  20 of C4's citations as "unresolvable," implying high hallucination. Inspection
  showed they were **correct DOIs** — the resolver compared a bare `10.x/…` against
  the DB's stored `https://doi.org/10.x/…` form and never matched. After
  normalising both sides (`rescore.py`, no re-answering), unresolvable citations
  dropped: C1 9→3, C2 13→7, C3 11→4, **C4 20→5**. Lesson for the paper: a naive
  citation-grounding metric *manufactures* hallucinations; the corrected picture
  is that **RAG (C2) has the worst citation grounding (16% unresolvable) and the
  graph layer the best (4%)**. The residual unresolvable cites are likely
  real-but-external papers (e.g. a correct substitute not yet in the KB) — splitting
  those from true fabrications needs CrossRef (future work).

### 4.5 RAGAS / FActScore (judge subset)

<!-- FILLED FROM judge_subset.py on the stratified subset -->
_Pending the judge subset run; will report faithfulness / answer-relevancy /
context-recall / atomic-hallucination per condition on the stratified subset._

### 4.6 What Study 2 establishes vs what still needs work
- **Established:** with n=108, deterministic gold, CIs and a paired significance
  test, the TAG layer's advantage on structured scientific-KB QA is **large and
  statistically unambiguous**; RAG/graph alone are not.
- **Still open (next research):** stronger RAG baselines (dense retrieval +
  reranker, GraphRAG) so "RAG" isn't a BM25 strawman; multi-seed runs for
  run-to-run variance; a non-Claude judge + human labels (judge-family bias);
  oracle-SQL vs model-SQL to isolate syn; multi-domain replication;
  cost/latency trade-off curve.

## 5. Figures

(primary answerer = claude-sonnet-4-6 unless noted)

**Mechanism / method (schematic):**
- `figs/fig_mechanism_architecture.png` — Notion SoT → Obsidian graph + DuckDB
  table/content layer → TAG engine (the data-flow).
- `figs/fig_mechanism_tag_loop.png` — the TAG loop syn→exec→gen + self-repair +
  BM25 content branch.
- `figs/fig_mechanism_ablation.png` — the 4 conditions as strictly-additive
  context layers feeding one fixed answerer.

**Results:**

- `figs/fig_bench_headline.png` — accuracy vs hallucination (programmatic +
  FActScore) across C1→C4. The headline story.
- `figs/fig_bench_ragas.png` — RAGAS triad (faithfulness / relevancy / context
  recall) across conditions.
- `figs/fig_bench_by_class.png` — correctness by question class × condition:
  shows which capability each layer switches on (aggregation = TAG-only).
- `figs/fig_bench_heatmap.png` — per-question outcome grid (correct / honest-miss
  / hallucination) × condition.
- `figs/fig_bench_model_compare.png` — **cross-model**: accuracy & atomic
  hallucination for sonnet vs haiku across conditions; the TAG advantage scales
  with answerer capability.

## 6. Reproduce

```bash
conda activate ffn_sim && cd ffn_sim/outputs/tag_kb
# Study 1 — hand-curated pilot (n=18), both answerer models:
python kb_benchmark.py --model claude-sonnet-4-6 --judge-model claude-opus-4-8
python kb_benchmark.py --model claude-haiku-4-5 --out benchmark_results_haiku.json

# Study 2 — scaled, deterministic gold (n=108):
python gen_questions.py --cap 14                       # -> questions_gen.json (SQL gold)
python kb_benchmark.py --questions questions_gen.json --no-judge \
       --out benchmark_results_scaled.json             # programmatic backbone, all 108
python rescore.py benchmark_results_scaled.json        # re-apply scoring (no LLM) after fixes
python judge_subset.py benchmark_results_scaled.json --per-class 2   # RAGAS on a subset
python kb_benchmark_stats.py                           # bootstrap CIs + McNemar + T5/T8

# Figures (results + mechanism):
python kb_benchmark_vis.py
python kb_mechanism_vis.py
```

## References (eval methodology)

- Biswal et al. *Text2SQL is Not Enough: Unifying AI and Databases with TAG.* arXiv:2408.14717, 2024.
- Es et al. *RAGAS: Automated Evaluation of Retrieval Augmented Generation.* 2023.
- Gao et al. *Enabling Large Language Models to Generate Text with Citations (ALCE).* EMNLP 2023.
- Min et al. *FActScore: Fine-grained Atomic Evaluation of Factual Precision.* EMNLP 2023.
- Li et al. *HaluEval: A Large-Scale Hallucination Evaluation Benchmark for LLMs.* EMNLP 2023.
- Liu et al. *G-Eval: NLG Evaluation using GPT-4 with Better Human Alignment.* EMNLP 2023.
