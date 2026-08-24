---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# H.7 — myosin generation-parameter audit (2026-06-07)

A 13-agent Workflow (one literature auditor per generation parameter → adversarial citation
verification → synthesis) to narrow the active-γ ceiling gap (0.084 mN/m → the real MCF7 datum
0.27, Hosseini 2020) **WITHOUT gate-chasing**. All 18 cited DOIs/PMIDs verified real; no
hallucinations. Result: **the honest literature-anchored set WIDENS the gap, it does not close it.**

## Verified per-parameter findings

| parameter | current | literature-anchored | citation (verified) | ceiling × | verdict |
|---|---|---|---|---|---|
| F_stall_per_head | 2.0 pN | **2.0 pN (keep)** | Chugh 2017 SI in-ensemble (40 pN/~20 heads); Norstrom 2010 NMIIB 2.2 pN (PMID 20511646); Finer 1994 skeletal 3-4 pN | 1.0 | ACCEPT (fix comment) |
| dipole length ℓ | 700 nm | **301 nm** | Billington 2013 EM contour 301±24 nm (PMID 24072716) | 0.43 | ACCEPT — 700 nm has no EM referent |
| heads/minifilament | 20 (10/side) | **56 (28/side)** | Niederman & Pollard 1975, 28 mol/filament (PMID 240861); Nagy 2012 (PMID 23148220); Billington 2013 | 2.8 | REVISE up — 20 was un-anchored "Hocky convention" |
| areal density | 3.0 /µm² | **~0.6 /µm²** (range 0.1-1) | Nie/Vavylonis 2015 HeLa medial cortex (PMID 25641802) | 0.2 | REVISE down — LOW confidence + **misattribution flag** |
| duty / engagement | ~4.5% engaged | realized-only, NOT ceiling | Kovács 2003/2007; Nagy 2012; Stam 2015 | 1.0 | ACCEPT — geometric throttle, not kinetics |
| active-stress x-check | 420 Pa | sanity only | Joanny-Prost 2009; Fischer-Friedrich 2014; Prost 2015 | 1.0 | ACCEPT |

## Honest combined result (NOT forced to target)
Product of accepted multipliers = 1.0 × 0.43 × 2.8 × 0.2 = **0.241**.
**γ_ceiling = 0.084 × 0.241 = 0.020 mN/m** — **~13× BELOW the MCF7 datum 0.27**. The down-corrections
(density ×0.2, dipole ×0.43) overwhelm the up-correction (heads ×2.8). Most generous defensible
combo (density kept 3.0 pending the foci-vs-minifilament ruling) = 0.084 × 2.8 × 0.43 = **0.101 mN/m**,
still ~3× below. **Even the most generous literature-anchored set keeps the ceiling under the datum.**
⇒ the active cortical tension is GENERATION-LIMITED, and honest anchoring STRENGTHENS that verdict.
This is the opposite of gate-chasing: three of four structural re-anchors push the ceiling DOWN.

## Active-stress cross-check (the physical question)
Corrected ceiling σ_a = 0.020 mN/m ÷ 200 nm ≈ **100 Pa**; MCF7 target 0.27 mN/m ≈ **1350 Pa**. Both
sit inside the measured cortical active-stress window (~50 Pa → ~10³ Pa canonical; Fischer-Friedrich
HeLa interphase ~1000 Pa). So 1350 Pa is physically precedented (canonical ~10³ Pa) — the corrected
model is at the LOW end. The model's σ_a = N·f·ℓ/(A·h) treats the cortex as INDEPENDENT force
dipoles; the real crosslinked/percolated cortex may amplify active stress beyond the dipole sum
(network/contractility amplification) — a FORMULA/network question, not only a parameter question.

## Engagement (realized-γ, not ceiling)
Only ~4.5% of heads engage — but the literature root cause is GEOMETRIC, not kinetic: per-head
zero-load duty ~0.05 (NMIIA) already gives P(≥1 of 56 heads bound) ≈ 0.94. The throttle is the
bipolar-pairing / capture-perp 210 nm / MAX_HEADS_PER_BEAD=3 / 824 nm head-placement geometry
(gate4_bipolar_pass 902/1783). So REALIZED γ would approach the (low) ceiling once the geometric
throttle is fixed — but the ceiling, not engagement, is the binding constraint after the re-anchors.

## Integrity flags (CLAUDE.md citation-integrity)
1. **`areal_density` = "Salbreux 2012, 3.0/µm²" is a CONFIRMED MISATTRIBUTION** — the cited review
   contains no such density. Treat like the project's 3 prior hallucinated sources. The real
   anchor is Nie/Vavylonis 2015 (~0.6/µm², HeLa, LOW confidence — non-MCF7, foci may bundle).
2. `heads_per_side=20` and `backbone_length=700 nm` were un-anchored magic numbers (no Magic-Number
   Block) — the audit supplies the EM-anchored values (Niederman 56 heads; Billington 301 nm).

## Recommended changes — these are MAGIC-NUMBER CONTRACT changes → PI sign-off required
ACCEPTED (verified, apply on PI go): `heads_per_side 10→28` (phase1_h3.yaml:328, phase1_h4.yaml:95;
+ the `/20.0` divisor at active_gel_seam.py:102 → `/(2·n_heads_per_side)`); `backbone_length
700e-9→301e-9` (phase1_h3.yaml:329, phase1_h4.yaml:96); F_stall keep 2 pN, fix the comment.
SURFACE-TO-PI (low confidence): `areal_density 3.0→0.6` (the misattribution + foci-vs-minifilament).

## The decision (PI)
Honest anchoring does NOT close the generation gap — the active cortical-tension ceiling is
~0.02-0.10 mN/m, ~3-13× below the real MCF7 datum 0.27, robustly generation-limited. Forks:
(a) apply the re-anchors + accept the floor + RE-TARGET the gate to traction/spreading;
(b) investigate the FORMULA — does the independent-dipole σ_a undercount the real cortex's
network-amplified active stress (~1000 Pa)? (measure σ_a from the actual percolated MD network, not
the dipole sum); (c) resolve the density foci-vs-minifilament ambiguity (the dominant uncertainty).
