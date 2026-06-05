# KU-3.5 v4 FA-anchored smoke — diagnosis (2026-05-30)

Smoke: `python -m ffn_sim.scripts.h3_ku35_v4_fa --smoke` (n_fil=60, equilibrate=True,
reconcile_dt=True). Raw: `/tmp/ku35_v4_smoke.json` (copy the JSON into this dir to retain).

## Result: builds + equilibrates cleanly, but γ stays on the floor (NOT ready for full production)

γ_total plateau = **3.13e-4 mN/m** (soft 2.1e-6 + rigid 3.11e-4) vs target [0.35, 0.65] → still **~1,600× under**.
This 3e-4 is constraint/construction noise, not contractile tension. Do NOT launch the multi-hour
production run — it would reproduce this floor and waste gbook hours.

## Root cause = TWO mechanism blockers (both independent of B1/B2; B1/B2 work correctly)

B1 verified working: dt reconciled 1.30e-8 → **6.52e-9 s**, binding element = FA `k_int_bare`
(stiffest spring). B2 verified working: equilibration drove max_force 0.186 → 4.6e-12 N.

1. **No substrate load path (FA clutch never loads).** Ligands at z=0, cortex shell at r≈10 µm.
   Equilibration relaxes the cell *in place* — it does NOT translate the cell down onto the
   substrate. So `n_fa_clutch_bonds_at_build = 9 / 52` integrins, and the dynamic Pereverzev catch
   bonds never engage: **`n_integrin_bound_final = 0`**. With no engaged clutch there is no
   substrate→cortex tension transmission → the FA-anchor mechanism the whole α restart depends on
   is geometrically not closed. Needs: a cell-onto-substrate placement step (S5/S6 geometry — bring
   the cortex shell into contact range of the z=0 ligand layer, or seed FAs at the contact zone),
   NOT a dt or code-wiring fix.

2. **Myosin does not step (`step_advances = 0`).** 83 heads engage (`bind_total=83`) but the Hill
   stepping never advances. This is the pre-existing KU-3.5 myosin binding-geometry / binned-ratchet
   defect flagged in MECHANISM_AUDIT_2026-05-30 §A2 — predates v4, independent of FA/B1/B2. Without
   stepping there is no active contractile stress regardless of the clutch.

3. Consequence: **r/r0 = 1.0000** (zero contraction) — consistent with both 1+2.

## Verdict / handoff
B1/B2 + the v4 driver are correct INFRASTRUCTURE; the floor is a MECHANISM-design gap (substrate
contact geometry + myosin stepping), which overlaps the parallel research session's active work.
Hold full production until: (a) cell-onto-substrate contact so the clutch loads, and (b) myosin
stepping fixed. Both are physics-design decisions, not unsupervised code patches.
