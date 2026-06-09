# Layer-2 seam inputs — CONSUMED H.7 artifacts (read-only)

These JSON files are **consumed copies** of H.7 single-cell fine-grained outputs, vendored here
so the Layer-2 scale-bridge (`ffn_sim/spheroid/h7_traction_seam.py`) is self-contained and
reproducible. They are **values/artifacts only** — the H.7 runtime code is owned by the
`h7/full-cell-integration` session and is never edited from the Layer-2 line.

| file | source (H.7 worktree) | what it carries |
|---|---|---|
| `h7_sf_2c_sarcomeric_gpu.json` | `/Users/sw1/ffn_cellsim/ffn_sim/outputs/h7/production/` | ⭐ DECISIVE single sarcomeric ventral-SF coherent contractile traction = **+131.08 ± 8.2 pN** (engaged heads 49.2 → +2.67 pN/head ≈ native F_stall). loop23, branch `h7/full-cell-integration`. |
| `h7_sf_array_SMOKE.json` | same | SF-array SMOKE (n_sf=3, random mixed-polarity control) → EXPANSILE/slackening (−), confirms only SARCOMERIC organization rectifies. NOT the decisive value. |

Provenance docs (H.7 session, read-only):
- `ffn_sim/docs/v2_audit/H7_MYOSIN_OVERLAP_MECHANISM_DESIGN_2026-06-09.md` §12 (decisive +131 pN).
- `ffn_sim/docs/v2_audit/H7_SF_ARRAY_TRACTION_SCALEUP_2026-06-10.md` (array scale-up + seam design §5).

⚠️ **PENDING from H.7 (provisional until landed):** the SF-array GPU aggregate
`h7_sf_array_n4_gpu.json` and the H.7-supplied seam record `h7_traction_seam.json` are
IN-PROGRESS on gbook. Until they land, the Layer-2 seam consumer computes the per-cell
aggregate PROVISIONALLY as per-SF (+131 pN) × N_SF (FA-pairing anchor). When H.7 writes
`h7_traction_seam.json`, drop it here and the consumer reads it directly (per-SF replaced by the
measured array aggregate).
