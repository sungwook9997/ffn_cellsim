# Speaker notes — "From Search Box to Self-Auditing Knowledge"

Companion to `KB_PRESENTATION.md` (full deck) and `KB_PRESENTATION.html` (one-file
slides). Audience: **biophysics / wet-lab researchers, not AI specialists.** Target
**~10–12 min** + Q&A. Each block = one slide: *what's on screen*, *what to say*
(spoken, ~30–60 s), and the *one number/visual to point at*.

---

## Slide 1 — Title + the problem (~60 s)
**On screen:** title; one line "every constant and every result in our simulator is supposed to trace back to a real paper."
**Say:** "We build a fine-grained cell-mechanics simulator. It has dozens of physical constants — off-rates, stall forces, tensions — and the whole credibility of the model rests on each one coming from real measurements, not from us tuning until it looked right. So the question I want to answer is: *how do you let an AI answer questions about that knowledge without it lying to you* — and prove every answer back to the paper it came from. The answer turned out to have three levels, and that's the talk."
**Point at:** the word "without lying."

## Slide 2 — Tier 1: RAG, the natural first idea (~90 s) · `rag_mechanism.png`
**Say:** "The obvious move with hundreds of PDFs is a smart search box. Modern RAG does exactly this: chop documents into chunks, turn each into a 'meaning coordinate' so similar text sits nearby, grab the nearest chunks to your question, and let the model write an answer. It's a brilliant assistant that's skimmed everything. *But* — it only ever does similarity matching over loose text."
**The failure (land this):** "Ask it 'which simulation parameter feeds a validation test that's currently failing?' and it returns paragraphs that merely *sound* relevant — and it will happily invent a plausible citation to fill a gap. That's fatal for a model whose rule is 'no magic numbers.'"
**Point at:** the red box — *no joins, no provenance, can fabricate a citation.*

## Slide 3 — Tier 2: TAG, make it relational (~75 s) · `tag_mechanism.png`
**Say:** "If the killer questions are about relationships and counts, the knowledge should live in a database, not loose text. TAG — Table-Augmented Generation, a 2024 method — is the recipe: the model turns your plain-English question into a precise database query, the database runs it *exactly*, and the model answers from the rows it got back. Now 'how many parameters trace to Smith 2019?' returns an exact, joined answer."
**The failure:** "But TAG assumes the database is *correct*. It will faithfully report 'this constant is sourced to Yao 2011' even if Yao 2011 never existed. TAG guarantees the *query* is exact; it guarantees nothing about whether the *contents are true*."
**Point at:** syn → exec → gen, then say "exact query, but blind to truth."

## Slide 4 — Tier 3: contract-graph + integrity gates (~90 s) · `contractgraph_mechanism.png`
**Say:** "Two things on top of TAG. First, instead of generic tables, eight purpose-built linked databases — source, claim, contract, parameter, gate, code, run, decision — so one query can walk the whole causal chain. Second, and the real point: three *auditors* that check the knowledge against what's actually on disk. One checks every cited paper exists. One checks every headline result still holds in the output files. One checks every simulation constant matches its config value *and* traces to a verified source."
**Point at:** the integrity-gate band, then the 5-hop arrow — "from a result, back to its paper, in one query."

## Slide 5 — Why each step is needed (~45 s) · `three_tier_comparison.png`
**Say:** "Same table, the whole argument in one view. Each tier keeps the last one's power and adds what it structurally couldn't do. RAG retrieves text. TAG adds exact relational joins. The contract-graph adds *enforced provenance* — and that bottom row is the one that matters: only the third tier can stop a wrong number from ever reaching the simulation."
**Point at:** the bottom row, all-green only in the last column.

## Slide 6 — It's continuous, not a one-time cleanup (~90 s) · `usage_in_simulation.png`
**Say:** "This isn't a one-off audit. Every constant is bound through a closed loop that fires on every code push, every refresh, every query. Don't take my word — the gates went live last week and immediately caught real drift: a forbidden way of measuring cell spreading headlined in five reports — the correct metric actually showed the cell *compacting*, not spreading; 'breakthrough' numbers that existed only in figure titles with no committed data; and an inflated speedup quietly corrected to our own record. These aren't citation typos — they're disk-versus-claim integrity errors that neither RAG nor plain TAG can catch."
**Point at:** the red panel of live catches. **Do NOT** lean on the older 'caught 3 fake papers' story — this is the fresh, ongoing evidence.

## Slide 7 — Live demo (~60 s) · `DEMO_5hop_query.md`
**Say:** "One command, the actual query. 'This run failed its gate — what constant is it testing, which paper does that number come from, and is that citation verified?' One SQL statement walks five hops and returns: failed run → gate → contract → gamma_cortex = 0.5 milli-newton-per-meter → Chugh 2017 → citation verdict CHECK, *not yet confirmed*. A search box cannot produce that row — the answer lives in the relationships, not in any chunk of text."
**Point at:** the single result row. (If presenting live: `python .../demo_5hop_query.py`.)

## Slide 8 — How it keeps being used (~45 s)
**Say:** "It's wired in so it can't rot: the blocking gates run in CI on every push and block the merge on drift; there's a one-command `make kb-check` and an opt-in pre-commit hook so you catch drift before you even push; and 13 tests guard the gate logic itself. The enforcement is automatic, not a discipline someone has to remember."
**Point at:** "blocks the merge."

## Slide 9 — What's genuinely novel + takeaway (~45 s)
**Say:** "To be honest about novelty: the query engine *is* TAG — we didn't reinvent that. What's new is two things: the domain ontology that lets provenance be traversed end-to-end, and the integrity gates that make the system *refuse to call a number validated unless an auditor agrees*. Takeaway: RAG finds what sounds relevant; TAG answers what's exactly related; this answers *what's actually true and still holds today* — and proves it back to the paper. A knowledge base that audits itself, every push."
**Point at:** the takeaway line. Stop there.

---

## Q&A — anticipated, with crisp answers

- **"Isn't this just TAG with extra steps?"** — Query-wise, yes, and I said so. The contribution isn't the engine; it's (1) the contract-graph *ontology* and (2) the *governance gates*. TAG answers questions over a database; it never checks the database is true. That check is the whole point here.
- **"Why not just auto-generate the wiki from PDFs (e.g. an LLM-built wiki)?"** — That's the opposite of what we need. Auto-generated content *introduces* the hallucination risk our gates exist to remove. Our source of truth is human-curated and every entry is citation-audited; an auto-writer would undermine that.
- **"N=1 — it's one project. Is this generalizable research?"** — As a *method paper* it's an integration of known parts, fair. Its real value is (a) a reproducibility/provenance discipline for mechanistic simulation that the field largely lacks, and (b) it makes *this* simulator's headline parameter-credibility claim machine-checkable. That's a methods/reproducibility contribution, not a new ML algorithm — and I'm not claiming otherwise.
- **"Some constants are 'UNSOURCED' / 'CHECK' — isn't that a weakness?"** — It's the opposite: the gate *surfaces* them honestly instead of hiding them. 3 of 7 declared constants are fully verified, the rest are flagged as pending — you can see exactly what still needs a source. A system that reported everything as fine would be the worrying one.
- **"What does it cost to run?"** — The gates are sub-second pure-Python over committed files — no Notion token, no network. They run in CI on every push and locally via one command. Effectively free.
- **"What happens if Notion (your source of truth) goes away?"** — Known gap; the gates already run off committed snapshots (the citation audit report, the config files), and committing a periodic Notion export is the planned fix so the whole chain is reproducible from git alone.

---

*Timing: slides 2/4/6 are the load-bearing ones — spend the time there; 5/8/9 are fast. If short on time, cut slide 5 (the comparison repeats the argument visually) and keep the live demo (slide 7).*
