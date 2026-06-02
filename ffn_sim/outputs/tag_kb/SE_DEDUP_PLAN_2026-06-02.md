# SourceEvidence duplicate-row dedup plan (PI review required)

**Date:** 2026-06-02 · **Status:** PLAN ONLY — no Notion writes performed (PI
decision 2026-06-02: script+plan first, PI reviews before any destructive SoT
edit). **Notion is the source of truth; archiving rows is hard to reverse.**

## Finding

`kb.duckdb` has 273 SourceEvidence rows but only **267 distinct citation_keys** —
**6 papers are double-entered in Notion** under the same citation_key, each as a
*distinct* Notion page with a *distinct* uid (e.g. SE5 **and** SE113). Same DOI in
every pair → genuine duplicates in the SoT (not a materializer bug). Likely the
early SE1–SE7 batch was re-entered later as SE82–SE114.

This does **not** break the TAG engine; it only inflates SE counts (273 vs 267
true papers) and splits a paper's claim-links across two rows.

## Decision rule

Keep the row that is **better integrated in the graph** (more `edges`); break ties
by the **earlier uid** (canonical SE1–SE7 batch). Each *dropped* row carries
**exactly one** KnowledgeClaim link the kept row lacks → **re-point that link to
the kept row first**, then archive the dropped row. (A plain archive without
re-pointing would orphan 6 claim→evidence links.)

## The 6 pairs

| citation_key | KEEP (uid / page id) | ARCHIVE (uid / page id) | re-point claim → keep |
|---|---|---|---|
| Bi2016_NatPhys | SE113 `372120daec5d81ecb9d7cf6749082eb5` | SE5 `371120daec5d815aa58adcb0ec8350af` | **KB-5.1** |
| ChughPaluch2018_JCS | SE83 `372120daec5d817da901eacb383ddb02` | SE2 `371120daec5d8125b1a9f8b15595400c` | **KB-3.5** |
| Salbreux2012_TCB | SE82 `372120daec5d81ceaccef2a21eb7db68` | SE1 `371120daec5d8128a789f08016dd8729` | **KB-3.5** |
| Murrell2015_NRMCB | SE84 `372120daec5d81a2bfd3f2db92746b7e` | SE4 `371120daec5d81979de5c242293a084e` | **KB-3.5** |
| Park2015_NatMater | SE6 `371120daec5d816f9dd9ffeba9acaef4` | SE114 `372120daec5d8135b659ceda2dc7d4df` | **KB-3.15** |
| Bieling2016_Cell | SE7 `371120daec5d8183a4e5f4f2248473ee` | SE97 `372120daec5d814aa515d1879bd9c80e` | **KB-3.7** |

> Note: KB-3.5 (cortex tension/mechanics) is the unique link on three different
> dropped rows — after dedup it will be anchored by the three kept cortex reviews
> (Chugh/Paluch, Salbreux, Murrell), which is the intended state.

## Execution (for the PI / a future session with Notion access)

For each pair, in order, via Notion MCP (`notion-update-page`) or the Notion API
with the integration token used by `notion_to_obsidian.py`:

1. **Re-point the claim.** On the KEEP SourceEvidence page, add the listed
   KnowledgeClaim to its `Claims` relation (and add KEEP to that claim's
   evidence relation). Confirm the claim now lists the KEEP row.
2. **Verify** the KEEP row's `Claims` is now the union of both rows' claims.
3. **Archive** the ARCHIVE page (Notion: set `archived=true` / move to trash —
   reversible within 30 days, then permanent).
4. After all 6: `bash refresh.sh` (or `python notion_to_duckdb.py`) to rebuild
   `kb.duckdb`; assert SE count = **267** and `count(distinct citation_key) =
   count(*)` in `source_evidence`.

## Verification query (post-dedup)

```sql
-- expect zero rows
SELECT citation_key, count(*) c FROM source_evidence GROUP BY 1 HAVING c > 1;
```

## Why not auto-executed

Per 2026-06-02 PI decision and CLAUDE.md (SoT/Notion changes need PI sign-off;
deletes are hard to reverse). All analysis is done and IDs are pinned above — this
is a mechanical 6×(re-point + archive) once approved.
