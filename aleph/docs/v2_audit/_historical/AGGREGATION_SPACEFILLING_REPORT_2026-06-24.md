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

# Deep report — aggregation & space-filling (PI's 3 questions, 2026-06-24)

Answering the PI's three questions about why the aggregated spheroid's cells stay ~round (asphericity
0.031) instead of filling empty space as deformed polyhedra. Built from a 6-agent investigation
(codebase audit + SimuCell3D research + tissue-mechanics literature) with a code-state verification pass
(all ON/OFF + file:line claims **CONFIRMED** against the actual source). Branch `h7/m1-ipc-contact`.

## TL;DR — the one-sentence answer

The machinery to fill space by **cortical tension + cohesion is already implemented and is the RIGHT
mechanism — but the single force that drives junction *faceting*, the cortical surface tension γ, is
turned OFF by default** (and set ~100× below the MCF7 value). Turn it on at the literature value and
equilibrate at N≥400, and cells deform-pack into polyhedra. Active "space-finding" (Q1) would HURT — it
disperses a confluent cluster. Proliferation (Q3) is a useful *second pass*, and SimuCell3D's key trick
for avoiding the mesh-collision the PI worries about is **remesh every step + a clean cleavage/re-triangulation
pipeline + daughter-volume = parent/2** (not collision detection).

---

## Q1 — Do cells have an ACTIVE ability to FIND/fill empty space? Should we build it?

**Short answer: the active machinery exists, but it is the WRONG tool for filling a confluent spheroid —
it would widen gaps and drive cells to escape, not fill space.**

- We DO have active motility built: `ActiveRimTraction` (coarse-grained lamellipodium + contraction belt
  + FA clutch, basal rim cells pull OUTWARD) and a lamellipodium graft (`cell/dcm_active.py`,
  `scripts/dcm_active_spheroid.py`). These are protrusion/traction mechanisms for **spreading / invasion**,
  not void-filling.
- Mechanistically, an "active space-finder" (protrusion biased toward free surface / low local density) is
  a **dispersal** force: in a confluent tissue it increases inter-cell gaps and drives single-cell escape
  (EMT/invasion), the OPPOSITE of confluent space-filling. The tissue-mechanics literature is clear that
  confluent epithelia fill space by **passive deformation-driven confluence**, not active migration.
- **Recommendation:** do NOT add active space-finding to solve the round-cell problem. Reserve
  lamellipodium/active-rim for EMT/invasion overlays. (Caveat: the dedicated deep-dive agent on the full
  active-mechanism inventory failed its schema retries; the directional answer above comes from the
  synthesis lens. If you want the exhaustive active-mechanism audit, I can re-run just that lens.)

## Q2 — Is space-filling by cortical tension + cohesion already implemented? (and why are cells still round?)

**Yes — MOSTLY implemented and running — EXCEPT the one force that does the faceting (surface tension γ)
is OFF by default and unanchored.** This is the core finding and the direct answer to "why round."

| mechanism | status | detail (verified) |
|---|---|---|
| **Turgor** (osmotic ΔP) | **ON** | `dcm_turgor_warp.py:35-72`; ΔP = dP0 + K_vol·(V0−V)/V0; V/V0≈1.13 (inflates cells ~13%, presses them outward) |
| **Node-node cohesion** (cadherin tent) | **ON** | `dcm_neighbor_warp.py:166-205`; bilinear attraction over [r_contact, c_adh], ω=adh 5e7. BUT it is **radial point-contact**, not a continuous flat interface |
| **Node-face contact** (excluded volume) | **ON** | `dcm_neighbor_warp.py:249-250`; penalty rep 2e8 over c_rep=0.30·mean_edge — keeps cells from overlapping |
| **Cortical surface tension γ** | **OFF ⚠️** | `dcm_neighbor_warp.py:531-576`, `--surface-tension` flag; the area-minimizing force that FLATTENS junctions — **disabled by default**, and default γ=1e-4 N/m is **~100× below** the MCF7 datum (~0.27 mN/m, Hosseini 2020) |
| **Node-face continuous adhesion** (flat interface) | **OFF** | `dcm_neighbor_warp.py:208-266`; the mesh-independent flat-interface adhesion — only fires with `--coupling`, OFF in the aggregation run |

**Why cells stay round (asphericity 0.031), mechanistically:** turgor inflates them symmetrically (~13%),
stiff contact repulsion prevents overlap, node-node adhesion only grabs at **sparse points** (not flat
faces), and **surface tension — the force that would minimise area and facet the junctions — is off.** So
the force balance has nothing pulling the shape toward a polyhedron. (At N=400 the *envelope* asphericity
even drops to ~0.0003 by curvature-averaging of 400 round cells; per-junction flattening Ψ~0.96 does occur
locally, but the cells themselves stay round.)

**The fix is a low-cost toggle, not new physics:** enable γ at the literature value and run a full
equilibration. The audit also surfaced 5 hygiene issues to fix alongside (so the toggle is clean, not a
magic number):
1. **Enable surface tension γ at a physiological value** (~0.5–2.5e-3 N/m, MCF7-scaled) — the primary fix.
2. **Resolve the K_vol inconsistency** (1e3 in the dataclass vs 7.73e5 hardcoded — a 773× split; pin one with a sourced derivation).
3. **Fix double-counted excluded-volume repulsion** (legacy node-node repulsion fires alongside node-face; run cohesion kernel with rep=0 and let node-face own excluded volume).
4. **Make the contact penalty monotone past c_rep** (it currently vanishes for deep penetration → tunnelling risk).
5. **Unify node-FACE adhesion** (continuous flat-interface + catch-slip kinetics) instead of the sparse node-node tent + a distance switch.

> Note: the OLD 2026-06-13 spheroid (`spheroid_nf_n400`) already reached **asphericity 0.013 + V/V0 1.12** —
> i.e. proper deform-packing has been achieved before. Reproducing it (γ ON, lit value, N≥400, long settle,
> remesh active) is the concrete target — and if it lands, the mechanism is confirmed with **no outcome-tuning**.

## Q3 — Proliferation to fill space, + the mesh-collision problem, + how SimuCell3D solves it

### Our proliferation (current state)
- **Fully implemented (C7):** host-managed rim-cell division on a dormant node pool (`cof=-1` sentinel),
  daughters parked 2.4R outward from the mother; works (fixes #1–#9 resolved the implicit-integration /
  rim-detection / binder re-pointing issues).
- **But three gaps that are exactly the PI's worry:**
  1. Division does **not** make cells less round (turgor still outweighs cohesion → daughters are round too).
  2. Mesh collision after a division is handled only **passively** (post-hoc contact penalty + the
     projection/implicit-CG step) — there is **no CCD or active gap-closure during the division event**, so
     a freshly-placed daughter can spike penetration.
  3. **Division and remesh are MUTUALLY EXCLUSIVE** (they share the `cof<0` pool) — so mesh quality is NOT
     maintained while cells proliferate. This is precisely how proliferation → "severe mesh collision."

### How SimuCell3D (Runser et al. 2024) avoids division-collision — the answer the PI remembered
SimuCell3D prevents collision **not** by collision-detection during division, but by keeping the mesh
*well-conditioned* and the contact *face-based*, every step:
1. **Clean cleavage + re-triangulation pipeline:** cleavage plane from centroid + PCA **largest** eigenvector
   → edge-plane walk to find the partition → 5-node face subdivision → **Poisson-disk + 2D Delaunay
   re-triangulation of the septum** on the division plane (production-ready in their `cell_divider.cpp`).
2. **Remesh EVERY step** (not our ~1000-step cadence) — enforces edge length l ∈ [l_min, 3·l_min] and a
   triangle-quality floor S_f>0.2, **eliminating sliver triangles before they accumulate**. *This is the
   single most important change* — slivers + stiff edge springs are what cause the force blow-up.
3. **Daughter target_volume = parent/2 immediately** at the split — avoids a pressure-jump transient that
   would drive the daughter into its mother.
4. **Node-vs-FACE signed-distance penalty contact** (not node-node LJ) — distributes force smoothly over a
   face triplet, removing the point-contact singularities our node-LJ creates.
5. **Surface tension γ + area elasticity** — *essential* to stop daughters from spreading unboundedly (ties
   straight back to Q2: γ is load-bearing for both faceting AND division stability).

### Concrete adoption plan (if we go the proliferation route)
- Port the cleavage/re-triangulation algorithm (delaunator, BSD-3, already in SimuCell3D).
- **Move remesh to the TOP of the step loop and run it every step (or ≤10 steps).**
- Set daughter V = parent/2 before the first force eval.
- De-conflict the pool sentinels (`cof=-1` remesh pool vs `cof=-2` parked daughters) so division + remesh
  co-run.
- Add a post-division soft-start (~100 sub-dt warmup) + a `G4_division` gate (no division-triggered pen >
  0.2·mean_edge, no NaN).
- Enable node-face contact + γ (shared with Q2).

## Synthesis — which route, in what order

For a **confluent tissue spheroid**, ranked by impact / cost / collision-risk:

1. **Route 2 — cortical tension + cohesion + turgor deform-packing — DO THIS FIRST (highest impact, lowest cost, correct mechanism).**
   Real epithelia fill space this way (active-foam / vertex-model / Bi–Manning jamming; shape index q*≈3.81 in
   2D, ~5.4 in 3D, sets the jam↔fluid transition). Our machinery exists; the work is *turning γ on at the lit
   value + full equilibration at N≥400 + remesh active*. Add a shape-index q* = P/√A observable to confirm the
   jamming (polyhedral) regime.
2. **Route 3 — proliferation — SECOND PASS, after confluence,** using the SimuCell3D remesh-every-step +
   clean-cleavage machinery so it does NOT blow up the mesh. Best coupled to a **tension-gated cell cycle**
   (contact-inhibition at q* > critical) — a positive feedback that fills space without escape.
3. **Route 1 — active migration — NOT for space-filling.** It widens gaps / drives escape; reserve for
   EMT/invasion overlays.

**Bottom line for the PI:** the spheroid stays round because the faceting force (surface tension γ) is OFF,
not because the approach is wrong. Turn γ on at the physiological value (Route 2) → cells deform-pack;
then proliferation (Route 3, with SimuCell3D's remesh-every-step) amplifies the fill without the
mesh-collision blowup. Active space-finding is the wrong lever here.
