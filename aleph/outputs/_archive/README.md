# `_archive/` — superseded output lines, moved 2026-08-20

Thirty directories and seven loose files moved here from `aleph/outputs/`, which had grown to
43 top-level entries mixing four eras of the project. Nothing was deleted; every path below is a
rename, so `git log --follow` still reaches the original commits and the bytes are unchanged.

## Why a move and not a delete

The PI's instruction was *"인용된 증거만 남기고 나머지 아카이브 이동"* — keep the cited evidence,
archive the rest. 6.65 GB of what sits under `aleph/outputs/` has never been committed
(`.gitignore:160` swallows every `*.npz`), so a delete would have been unrecoverable rather than
merely inconvenient. An archive move costs nothing and is reversible with one `git mv`.

## What stayed at `aleph/outputs/`

| path | why it stayed |
|---|---|
| `ac/` | the current engine line. Every artifact `state_rows.yaml` cites lives here. |
| `outer/` | the outer layer's own artifacts. |
| `tag_kb/`, `obsidian_rag_full/` | LIVE KB infrastructure — `make kb-check` and `refresh.sh` read these paths, and 80 files reference them. Moving them for tidiness would have broken working machinery. |
| `deck_assets/` | deliverable assets, still referenced by the presentation build. |
| `_historical/` | already an archive marker before this move. |

## What moved

The DCM / FF / H-series eras, superseded by the `ac` engine:

`ac_ecm` `ac_magnitude` `cbm` `compartment_smoke` `cortex_mesh` `cupy_cache` `dcm_compaction`
`engine_reval` `ff` `ff_single_opt` `gpu_signoff` `h1` `h1_baoab_freeze` `h2` `h3` `h4` `h5` `h7`
`h8` `h9` `h_0_2` `h_dcm_active` `h_dcm_gpu` `h_dcm_gpu_lod` `h_dcm_native` `h_dcm_spreading`
`h_dcm_two_stage` `layer2` `mech_hier` `warp_decohesion`

plus `loose/` — seven files that were sitting directly at the `outputs/` root with no directory.

## ⚠ What this breaks, stated rather than discovered later

Historical documents cite these results by their OLD paths (`aleph/outputs/ff/...`,
`aleph/outputs/h3/...`). Those citations now point one directory short. The files are all still
present under `_archive/<same name>/<same subpath>`, so a broken citation is repaired by inserting
`_archive/` — but nothing was rewritten to do that automatically, because editing dozens of
historical records to preserve link cosmetics is a worse trade than a documented offset.

**No `STATE.md` tier-(a) row is affected**: all 20 cite artifacts under `ac/`, which did not move.
