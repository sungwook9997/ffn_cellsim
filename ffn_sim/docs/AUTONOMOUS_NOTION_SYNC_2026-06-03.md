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
- **last processed HEAD: `dfc6cac`** (iter-3; L2.6 ligand/active-traction posted. Scan `dfc6cac..HEAD` next.)
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

### iter 2 — 2026-06-03 ~02:30 KST — **L2.5 G3-PASS milestone (verifying before posting)**
- New session commits `b8a95a4..9a79a9f` (6): L2.5a catch-bond oracle (`7a6a11b`), L2.5b/c
  (`fee53d1`), L2.5c+L2.6 substrate (`fcd2b20`), --catch/--morse flag (`bbf02b2`),
  **L2.5 HEADLINE G3 PASS r²=0.98 (`5113f70`)**, L2.6 ligand axis (`9a79a9f`).
- ⚠️ **Did NOT blindly post the "G3 PASS r²=0.98 EMERGES" claim** (the auto-classifier rightly
  blocked it, and it's exactly the L2.4a 0.987-was-under-sampled trap). **Independently verifying
  first** (the whole reason PI asked for this check):
  - Re-fit the Lead's `growth_sweep_catch.checkpoint.json` (5 sizes) → **r²=0.980 reproduced**
    (a=-0.33, b=188.7µm, c=-2655µm²; matches REPORT). morse ckpt re-fit = r²=0.874 (n=4,
    ±0.097) — re-confirms morse is seed-unstable (0.74–0.87), catch tight (±0.082).
  - Running N0=400 (the size where morse fragments) × **5 FRESH seeds (2000-2004)** × {catch,morse}
    to test the crux claim "catch 0/5 vs morse 1/5 fragment" with independent seeds. (in progress)
- **VERIFIED ✅** — N0=400 × 5 fresh seeds: **catch 0/5 fragment (core 1.86±0.029), morse 1/5
  (one seed hull=739×, core 1.95±0.119)**. Reproduces "catch 0/5 vs morse 1/5" with independent
  seeds; catch ~4× tighter. Combined with the re-fit r²=0.980, L2.5 G3-PASS is **genuinely robust**
  (unlike L2.4a's under-sampled 0.987 — here the catch bond prevents fragmentation at the
  *mechanism* level, not just the measurement). Artifact: `outputs/layer2/l2_5_catch_fragmentation_verify.jsonl`.
- **Posted** L2.5/L2.6 Dev-Log (verified) → `373120daec5d81eb9d3ddf9dfe1b202c`.
- SE dry-run: 0 new rows (50 skipped, all registered) — iter-1 caught all new papers.
- **iter-2 COMPLETE.** last processed HEAD → `9a79a9f`.

### iter 3 — 2026-06-03 ~03:15 KST — L2.6 ligand/active-traction
- New commits `9a79a9f..dfc6cac` (6, minus my 85f6525): L2.5 A/B (`1ead471`: morse r²=0.80 FAIL
  vs catch 0.98 PASS, same pooled path — **corroborates iter-2 fresh-seed verification**);
  L2.6 ligand axis (`02e1db7`), oracle move (`69abc66`), active traction in pooled growth
  (`b2d86db`), REPORT (`f32d50e`), traction→b-coefficient (`dfc6cac`).
- **Posted** L2.6 ligand Dev-Log → `373120daec5d81588be5e336792d3196`. Faithfully relays the
  honest finding: passive Bare/Pre/Lam4 adhesion ≈ flat (overlap, Lam4 slightly higher);
  **ACTIVE edge-traction is the real ligand driver** (core A/A₀ 1.53→1.56→1.81 at 0/3/6 nN);
  b-coefficient grows with traction (0/3 nN: 2.55/2.25/1.82→2.80/2.35/1.92). Caveat carried:
  6 nN hit a numerical box/wall limit (known fix, not physics). No success over-claim → posted.
- SE dry-run 0 new (172 PDFs on disk now, +2 = Rakshit SI / Lou-Zhu already registered).
  No 8-DB change → only Obsidian devlogs refreshed (60 day-log nodes; L2.6 page pulled).
- **iter-3 COMPLETE.** last processed HEAD → `dfc6cac`.
