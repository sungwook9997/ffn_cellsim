# Phase C — M1 contact fix design (inside-mesh / signed-distance contact) — 2026-06-23

**The night's critical finding:** the per-face penalty contact (`contact_grid_kernel`) is the single
blocker gating BOTH clean aggregation AND the spreading test. It fails under the lit cohesion strength:
the proxy-free spreading run hit `pen=3.1` (nodes ~7µm into neighbours). This doc lays out the fix
options + a tested prototype, for PI to ratify (it changes the contact force law → a contract change).

## The bug, precisely
`contact_grid_kernel` repulsion: for EVERY other-cell face with `sign = (p−closestpoint)·n̂ < 0` AND
`min_d < c_rep (≈0.67µm)`, add `force = r_vec·rep·area` (depth-proportional, |r_vec|=min_d).
- The `min_d < c_rep` cap means a node penetrating deeper than c_rep gets ZERO force → **tunnelling**.
- Naively removing the cap **EXPLODES** (verified: NN 1.58→6.1): a node OUTSIDE a cell still sits on the
  inner half-space (`sign<0`) of MANY of that cell's oblique/far triangles within the query radius
  (con_q/c_rep ≈ 5); summing all of them = a huge spurious outward force. The `min_d<c_rep` gate was
  silently MASKING this by keeping only near faces. So the gate is load-bearing and cannot be dropped.
- Root cause: **per-face `sign` is a LOCAL half-space test, not a global inside-closed-mesh test.**

## Options
**A — Winding-number / solid-angle inside test.** Per node, sum the signed solid angle of a candidate
other-cell's closed mesh; ≈1 ⇒ inside ⇒ push to nearest surface ∝ depth. Exact, handles non-convex/
deformed cells. Cost O(node × faces_of_candidate_cell) — heavier; robust. Best correctness, most work.

**B — Nearest-face signed-distance repulsion (RECOMMENDED).** Keep the existing closest-point machinery
but for REPULSION use only the SINGLE NEAREST other-cell face (min `min_d` over candidates), no c_rep
cap. For an OUTSIDE node the nearest face has `sign>0` → no repel (the far oblique faces are never the
nearest, so the explosion can't happen). For a PENETRATING node the nearest face is the entry face,
`sign<0`, and the existing `force = r_vec·rep·area` is already depth-proportional → deep nodes are
pushed out ∝ depth. Adhesion stays multi-face (genuinely multi-contact). This is the standard FEM /
SimuCell3D node-face contact done right; minimal change to the kernel; O(same as now). Risk: if a node
tunnels past the GRID QUERY RADIUS (con_q≈3.4µm) the entry face is lost → pair with a modestly larger
contact query radius or an active-set carry-over. Changes the force law (single vs summed repulsion) →
not byte-identical even on gentle runs → PI ratification + full re-baseline.

**C — Centroid spherical backstop (cheap, approximate).** Per node, nearest other-cell centroid; if
within R_eff push radially out ∝ (R_eff−d). Catches catastrophic tunnelling using global cell shape;
O(node×cell), trivial. But approximate (assumes spherical; a flattened cell isn't) and not mechanistic
— a stopgap, not the real fix. Could layer UNDER B as a guaranteed last-resort anti-merge guard.

## Prototype (Option B) — tested 2026-06-23
`scripts/m1_nearest_face_contact_proto.py` monkeypatches `contact_grid_kernel` with the nearest-face
variant and runs (a) gentle 2-cell (flattening must survive) and (b) strong-adhesion 2-cell (pen must
drop). RESULTS:

**(a) GENTLE 2-cell (adh=5e7) — flattening PRESERVED ✓**
| | NN/R | oblate | Psi | deep/R | pen |
|---|---|---|---|---|---|
| baseline (per-face penalty) | 1.578 | 0.9045 | 0.956 | 0.047 | 0.048 |
| Option B (nearest-face) | 1.677 | 0.9449 | 0.969 | 0.058 | 0.128 |
Option B still flattens (oblate 0.94<1) and stays G2-clean (pen 0.13<0.3), but settles slightly farther
(NN +6%) and flatter-by-less — so it CHANGES the equilibrium (single-face vs summed repulsion is a
different law) → NOT byte-identical even on gentle runs → full re-baseline + PI ratification required.

**(b) STRONG-force validation — INCONCLUSIVE (test mis-designed).** My first strong case used raw
`adh=5e8`, which NUMERICALLY DIVERGES for BOTH baseline and Option B (NN→13–15R, pen→800–2100, vv0→1.2)
— that is solver blow-up at an over-strong adhesion, not the penetration mode, so it tells us nothing.
The correct penetration-inducing test is the actual cadherin ×40 bundle (the regime that gave the N=100
spreading pen=3.1); that run is pending. **So Option B is gentle-validated but the decisive claim — does
it cut the real bundle-regime penetration — is NOT yet established.**

## Recommendation (honest, partial)
Option B (nearest-face signed-distance) is the most tractable proper fix and is **promising** — it
preserves gentle flattening and, by construction (single nearest face, never the far oblique faces),
avoids the explosion that killed the naive gate-removal. BUT (1) it changes the contact law → re-baselines
every prior equilibrium, and (2) its core promise (cutting the real bundle-regime pen 3.1) is still
UNVALIDATED. **Do not adopt without:** the cadherin-bundle A/B (pen must drop), a parity-gated `--nearest-
face-contact` flag (default OFF), and PI ratification. Option A (winding-number) remains the
gold-standard if B's nearest-face proves insufficient for deep multi-cell tunnels (where the nearest face
can be lost past the query radius). Option C (centroid backstop) is a cheap last-resort anti-merge guard
to layer under either. **This is a design + prototype, NOT a shipped fix.**
