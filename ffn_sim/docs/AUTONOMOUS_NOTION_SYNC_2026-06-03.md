# Autonomous Notion / storage sync — 2026-06-03 (PI asleep ~10h)

Cleanup/sync session running an hourly `/loop` (cron `652871ef`, every :13, session-only)
while other live Claude sessions keep committing to `layer2/spheroid-cbm`. Goal: keep Notion
+ Obsidian + TAG in sync with what the sessions land, register new verified reference papers
as SourceEvidence, and tidy storage — **without deleting anything from `references/`**
(off-topic papers simply aren't pushed to Notion/Obsidian).

## Standing protocol (each iteration)
1. `git log <last-HEAD>..HEAD` across all sessions → new commits.
2. Post un-Notioned session work as Dev-Log child pages under board `365120daec5d81969e74ffbb757d55c8`.
3. `references_to_se.py` (dry-run → review → `--commit`) for new on-topic refs; SE data_source
   `03821ec4-002d-473d-91a5-eb3dca7afb3c`. Off-topic (Source Type stays BLANK) → flag, don't wire.
4. Refresh TAG (`notion_to_duckdb.py` + `references_ingest.py` in `outputs/tag_kb/`) and Obsidian
   (`notion_to_obsidian.py`, `add_code_nodes.py`, `devlogs_to_obsidian.py` in `outputs/obsidian_rag_full/`).
   ⚠️ `refresh.sh` fails under cron (`python: command not found`) — call the .py steps via `conda run -n ffn_sim`.
5. Tidy storage (dupes/tmp only; never delete `references/`).
6. Append an iteration entry below.

⚠️ **Do NOT commit/edit other sessions' uncommitted shared files** (REPORT.md, proliferation.py,
sweep scripts). Additive/safe ops only. Shared-tree collision is live (see memory layer2 line).

## Key IDs
- Dev Logs board: `365120daec5d81969e74ffbb757d55c8`
- SourceEvidence DS: `03821ec4-002d-473d-91a5-eb3dca7afb3c` · KnowledgeClaim DS: `f0c7baf0-9761-4934-98db-eecbbf7462f4`
- L2.4a/b Notion page: `373120daec5d815e8df0cbdbc6aaa7cd` · KB-v3: `373120daec5d81f19bbfe4086277f5ea`

## State
- **last processed HEAD: `7dcca72`** (update each iteration)
- TAG kb.duckdb rebuilt 2026-06-03 ~01:22 (144 keys / 6292 chunks / 150 PDFs; 94 SE-linked).
- SourceEvidence count before iter-1 SE creation: 273.

## Iteration log

### iter 1 — 2026-06-03 ~01:25 KST (setup + backlog)
- Cron `652871ef` armed (hourly :13). Tracking file created.
- Today's commits audited: L2.4 (3622f0f), L2.4a (08052f5), L2.4b (266af51), my L2.4a-verify
  (7dcca72), Benchmark v3 (a8ef474/e13f4e7/ef344b0), kb_benchmark --questions (0d3f858),
  γ-watch (9a8fe81). Notion coverage: L2.4 ✓ (page 373120…f2dbaa), L2.4a/b ✓ (page …a7cd),
  KB-v3 ✓ (page …f5ea). **Un-Notioned = 9a8fe81 γ-production crash watch** → posting Dev-Log.
- SE registration: `references_to_se.py --commit` = **16 new SE rows created** (proper Crossref
  Author-Year keys — fixes the slug-key problem; SE 273→289), 34 skipped (already in SE incl.
  Rakshit2012, Stephens2017, Murrell2015), 1 no-DOI (Rakshit SI, correctly not created).
  New on-topic (15): Manibog2016 + Lou2007 (**catch-bond, L2.5**); Chugh2017, Murrell2012,
  Ennomani2016, Lenz2012, Stam2017, Linsmeier2016, Banerjee2017, Tam2021, Belmonte2017,
  McFadden2017, Popov2016/MEDYAN, Freedman2017 (**actomyosin-network cluster, H.3**);
  JuelPrtner2025 (**MCF7 nucleus AFM**).
- **Off-topic flagged:** `Balog2007_BiophysicalJournal` (SE-288, yeast PGK) renamed
  `ZZ-OFFTOPIC …` + Anchor Status marks it for manual deletion (MCP can't archive). Kept in
  references/ per PI; only excluded from KB intent. PI: delete this row in Notion UI.
- **γ-production crash recurrence** posted as Dev-Log `373120daec5d818192a4e51080d36e54`.
- TAG kb.duckdb + Obsidian re-refreshed (SE 273→289 now rendered as nodes + PDF-linked by DOI).

#### SE→KnowledgeClaim relation backlog (deliberate human-curation per the pipeline; NOT auto-wired)
`references_to_se.py` leaves Source Type BLANK + no claim relations by design. Suggested wiring
for PI review (high-confidence):
- **KB-4.2 cadherin catch-bond** (`372120daec5d81aca6c2dee3a8e86072`, "E-cadherin trans-bond") ←
  supports ← Manibog2016, Lou2007 (+ existing Rakshit2012). These 3 = the L2.5 anchor set.
- **KB-3.x cortex/myosin** ← Chugh2017 (cortex tension), Murrell2012/Lenz2012/Linsmeier2016/
  Stam2017/Ennomani2016/Banerjee2017/Tam2021 (actomyosin contractility), Popov2016/Belmonte2017/
  Freedman2017/McFadden2017 (network simulators). (Map to specific KB-3.y at curation.)
- **KB-H9 nucleus** ← JuelPrtner2025 (MCF7 nucleus viscoelasticity, AFM).
Set Source Type per row when wiring (most = "model/method" or "measurement").

- **Storage observation (NOT deleted — PI said keep references/):** `references_ingest` skipped
  **22 exact byte-duplicate PDFs** (same file under two names, e.g. `Chen2017_MBE.pdf` ==
  `10.3934_mbe.2018016.pdf`; `Arruda2026_NpjSystBiolAppl.pdf` == `s41540-026-00648-9.pdf`;
  several `… (1).pdf` copies). Harmless to TAG (deduped by content hash) but clutters the folder.
  Candidate for a future de-dup pass (PI approval) — listed in `/tmp/sync_refresh.log`. Left as-is.
- Verified in vault: new SE nodes (Manibog2016, Lou2007, Chugh2017, Murrell2012, JuelPrtner2025…)
  + Dev-Log nodes (γ-crash, L2.4a/b, KB-v3). Vault 489→505 notes; SE-linked PDFs 94→110.

last processed HEAD → still `7dcca72` (no new session commits during iter-1).
**iter-1 COMPLETE.** Next cron fire :13 → re-scan for new session commits + repeat.
