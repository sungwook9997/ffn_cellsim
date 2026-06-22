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

**(b) BUNDLE-REGIME validation — Option B REFUTED (it DIVERGES).** Tested the actual failure mode
(N=12, cadherin ×40 bundle, the regime that gave the N=100 spreading pen=3.1):
| | pen_final | vv0 |
|---|---|---|
| baseline (per-face penalty) | 2.92 | 0.9999 (penetrates but STABLE) |
| Option B (nearest-face, no cap) | 3.68 | **175.3 — CATASTROPHIC BLOW-UP** |
Option B does NOT fix it — it makes it WORSE and **diverges** (cells invert, vv0→175). Root cause: the
depth-proportional force `|F| = min_d·rep·area` is UNBOUNDED once the c_rep cap is removed — a node at
depth d feels a force ∝ d, so as the ×40 bundle drives it deeper the force runs away → ejection → cell
inversion. **(My earlier adh=5e8 "strong" test also diverged for both — same uncapped blow-up.)**

## ⭐ The real lesson — the `c_rep` cap is DOUBLY load-bearing
The `min_d < c_rep` gate serves TWO purposes, and tunnelling is the price of BOTH:
1. **Localization** — restricts repulsion to near faces (removing it → far-oblique-face explosion, NN→6.1).
2. **Force-bounding** — caps the depth-proportional penalty so it can't run away (removing it → vv0→175).
So you cannot fix tunnelling by simply extending the penalty past c_rep in EITHER direction (nearest-face
OR all-faces) — both break. **The proper fix needs all three at once:** a true inside-closed-mesh test
(winding number / generalized winding — Option A, so far faces never spuriously fire), a **BOUNDED
barrier** force (e.g. a saturating/log-barrier penalty that stays finite at large depth instead of the
linear `min_d·rep·area`), and **implicit treatment** of that contact stiffness so the large near-contact
forces don't blow the step. This is a genuine contact-mechanics piece (IPC-style barrier + CCD, or
SimuCell3D's actual contact model), NOT a one-kernel tweak.

## Recommendation (honest)
**Option B is REFUTED — do not pursue the nearest-face-no-cap path.** The M1 fix is harder than any of
the three sketched options alone: it requires Option A's inside test + a bounded (non-linear, saturating)
barrier + implicit treatment. Until then the per-face penalty (with its tunnelling) is the least-bad
option, and the standing fact is: **at the lit cohesion strength the contact either tunnels (capped) or
diverges (uncapped) — aggregation/spreading are blocked on a proper barrier contact.** Surface to PI as a
scoped contact-mechanics task (reference: Incremental Potential Contact / IPC barrier methods). The
prototype + this negative result save the PI from the nearest-face dead end. **Design + negative
prototype, NOT a shipped fix.**
