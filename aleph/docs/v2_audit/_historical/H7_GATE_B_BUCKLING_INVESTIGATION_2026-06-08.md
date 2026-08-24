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

# H.7 Gate-B design-investigation #1 — why the backbone does not condense

**Date:** 2026-06-08. Branch `h7/full-cell-integration`. Author: Lead.
Feeds `H7_GATE_B_CONTRACT_2026-06-08.md` §3 (the design-investigation item: isolate
the condensation suppressor before building the relaxed-constraint mode).

## Question

Gate-A (locked REFUTE) shows the real grip_walk contraction (γ_soft 16.2× to 3.06e-3)
does NOT transmit into hoop tension, with `r/r0 = 1.0000` (no buckling/condensation —
the Murrell/Lenz/Miyazaki symmetry-breaking). The backbone has SOFT bending angles
(`connected_mesh.py` cortex-angle harmonic), so buckling is geometrically allowed in
principle. Why doesn't it happen?

## Method — Euler buckling threshold vs the myosin per-head load

Cortex params (resolved from `configs/phase1_h3.yaml`, no tuning):
- L_filament = 3.0 µm, beads_per_filament = 7 ⇒ segment ℓ₀ = L/(7−1) = **500 nm**.
- bending modulus κ_B = ℓ_p·k_BT = 17 µm · 4.28e-21 J = **7.0e-26 N·m²** (config value).
- myosin per-head stall F_head = **2.0 pN** (`mcf7_baseline.yaml`, KU-2.4/2.18, Chugh 2017).

Slender-rod Euler: F_crit = π²κ_B / L², where L is the **free (unsupported) length**
between hard anchors.

| free length L | F_crit | F_head/F_crit | buckles? |
|---|---|---|---|
| whole filament 3 µm | 0.080 pN | 25.1 | YES, 25× over |
| crosslink spacing ~0.45 µm | 3.55 pN | 0.56 | no |
| segment ℓ₀ 500 nm | 2.87 pN | 0.70 | no |
| 250 nm | 11.5 pN | 0.17 | no |

## Finding — buckling suppression is a FREE-LENGTH effect, at the margin

Over its full 3 µm contour a filament would buckle at 0.08 pN — a 2 pN head would buckle
it **25× over**. It does not, because the filament is **pinned at ~segment/crosslink
scale** (M-SHAKE fixed-length segment nodes ℓ₀ = 500 nm + crosslinks ~0.45 µm). At that
free length the Euler threshold (2.9–3.5 pN) **just exceeds** the per-head 2 pN
(ratio 0.56–0.70). The system sits **right at the buckling margin** — not deeply
suppressed, just over the edge.

This confirms the contract's degree-of-freedom correction: the suppressor is **not** axial
compliance (k_axial ≈ 154 pN/nm, correctly inextensible) — it is the **anchor density /
effective free length** (candidates i + ii jointly), with angle_k contributing only through
κ_B in F_crit.

## Two implications for the relaxed-constraint mode

1. **The real lever is the free buckling length.** "Permit buckling" = let the filament
   present a longer unsupported span between hard constraints (relax the M-SHAKE segment
   pinning on the COMPRESSION side), so a ~2 pN load crosses F_crit. NOT relaxing axial
   stretch.
2. **The margin is narrow ⇒ myosin binding throughput is decisive.** Per-segment load is
   ~1 head (2 pN) because myosin under-binds (~12/2000 heads engaged, overnight finding).
   A minifilament places ~56 heads; if several load one free span, the summed compressive
   force clears F_crit immediately. So buckling onset couples to the known binding-throughput
   issue — more engaged heads per span tips the margin.

## Caveats (first-order)

Straight-rod Euler idealization, single-head point load, free length assumed = crosslink
spacing. The real cortex is a pre-stressed network; the actual onset depends on the
constructed mesh's anchor-to-anchor free-length distribution and the local multi-head force
balance. **Next step (probe):** measure, on the built physiological cell, (a) the anchor-to-
anchor free-length distribution along filaments, and (b) the actual per-span myosin load,
to pin the margin quantitatively before sizing the relaxed-mode compression-release.

Figure: `outputs/h7/figs/h7_gate_b_buckling_euler.png` (F_crit vs L with the 2 pN line and
the segment/crosslink scales marked).
