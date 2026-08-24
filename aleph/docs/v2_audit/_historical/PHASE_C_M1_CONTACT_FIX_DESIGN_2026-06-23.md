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

## Third data point — SATURATING barrier (nearest-face, force capped at the c_rep level): STABLE but TOO WEAK
Tested a bounded variant: nearest-face, `|F| = rep·area·min(min_d, c_rep)` — force grows to the c_rep
level then SATURATES (constant for deeper penetration). On the cadherin ×40 bundle (N=12):
| | pen_final | vv0 |
|---|---|---|
| baseline (per-face cap) | 2.921 | 0.9999 |
| saturating nearest-face | **2.876** | **1.0000 (STABLE)** |
It cures the divergence (vv0=1.0, no blow-up) but **does NOT reduce the penetration** (2.92→2.88). Why:
a constant-capped force can't out-push a stronger constant force — the node equilibrates deep where the
×40 bundle (~7nN) balances the capped contact (`rep·area·c_rep`); deeper penetration feels no extra push.

## Fourth data point — can stiffer `rep` salvage the per-face penalty? NO (rep-sweep)
Swept `rep_strength` on the cadherin ×40 bundle (N=12); pen does NOT drop, and stiff values diverge:
| rep | 2e8 | 5e8 | 1e9 | 2e9 |
|---|---|---|---|---|
| pen_final | 2.92 | 2.82 | 3.58 | 3.61 |
| vv0 | 1.00 | **7.2 (diverge)** | 1.00 | **39.3 (diverge)** |
Stiffer rep does NOT reduce penetration (pen 2.9→3.6, slightly UP) and at 5e8/2e9 it DIVERGES. Why: the
deep nodes that cause the pen are PAST c_rep, where the force is ZERO regardless of rep (the tunnelling);
raising rep only stiffens the shallow shell → explicit-step blow-up. **So the per-face penalty cannot be
salvaged by stiffness tuning — the fix must be a force that GROWS at depth, i.e. the IPC barrier.**

## ⭐ The three prototypes BRACKET the fix space → IPC log-barrier is the answer
| contact force vs penetration depth | result |
|---|---|
| linear UNBOUNDED (Option B, no cap) | DIVERGES (vv0→175) — too stiff for the explicit-ish step |
| constant SATURATED (cap at c_rep) | STABLE but TOO WEAK (pen 2.9 — can't beat the bundle) |
| → **must GROW near contact but stay solvable** | **IPC log-barrier** `F ∝ −ln(gap)` → ∞ as gap→0, **handled IMPLICITLY** so the stiffness goes in the Hessian/Newton, not an explicit force spike |
A contact that REDUCES penetration must apply MORE force the deeper a node is — but a force that grows
explicitly blows the step (Option B). The resolution is the **IPC (Incremental Potential Contact)
log-barrier**: an energy barrier whose force →∞ as the gap →0, integrated IMPLICITLY (the existing
`device_cg` operator already carries contact stiffness — the barrier Hessian goes there) + continuous
collision detection so no node crosses a face within a step. This is the only one of the four force-laws
that is simultaneously non-tunnelling, bounded-in-practice (implicit), and strong enough to beat the
bundle. **Scoped to PI as a contact-mechanics task (ref: Li et al. IPC, SIGGRAPH 2020).**

## Fifth data point — naive IPC-barrier KERNEL also fails (pen→69): a barrier force ≠ the IPC METHOD
Prototyped the recommended log-barrier as a drop-in nearest-face kernel (barrier `rep·area·(c_rep/gap−1)`
→∞ as gap→0 for outside nodes + a strong recovery for inside nodes, in the implicit operator;
d_hat=c_rep, kappa=rep, both derived). On the cadherin ×40 bundle (N=12): **pen 2.92 → 69.3, vv0 → 0.86
(BLOWS UP).** Why: the cells START overlapping (pen high) so most contact nodes are already INSIDE
(gap<0), where the barrier is undefined and the "recovery" injects huge ejection forces; and the →∞
barrier's JVP is too stiff for the matrix-free CG (guards return bad steps). **A barrier FORCE in a kernel
is NOT the IPC METHOD** — IPC's whole point is the **CCD-filtered line search that GUARANTEES gap stays
> 0 every step** (so the barrier is never evaluated at gap≤0), plus the barrier's exact gradient/Hessian
in a filtered Newton solve. Without CCD the barrier is worse than the penalty.

## ⭐⭐ Exhaustive conclusion — 5 prototypes, the per-face/kernel framework CANNOT solve M1
| prototype | result |
|---|---|
| linear unbounded (drop c_rep cap) | DIVERGES (far-face explosion NN→6.1) |
| nearest-face unbounded | DIVERGES (vv0→175) |
| nearest-face saturated cap | STABLE but TOO WEAK (pen 2.9→2.9) |
| stiffer rep (sweep) | NO drop (pen 2.9→3.6) + diverges at 5e8/2e9 |
| naive IPC barrier kernel | BLOWS UP (pen→69, no CCD) |
| constraint projection (post-step relocate-to-surface) | BLOWS UP (vv0→244) |
**No kernel-level contact law in the penalty/barrier-force family solves it. Nor does the position-
constraint family (6th prototype): a post-step projection that relocates inside nodes to the surface
makes it WORSE (vv0→244) — the discontinuous relocation is seen by the next implicit step as a huge
displacement, the turgor/edge forces react violently → energy injection → blow-up. A correct constraint
must be solved WITHIN the implicit step (a Lagrange multiplier / KKT system), not bolted on after it —
the same "real method, not a hack" lesson as IPC's CCD-in-the-line-search.** The fix is the FULL IPC
method — **CCD-filtered line search (guarantees no penetration) + barrier energy + filtered Newton** — a
real contact-mechanics implementation (Li et al., *Incremental Potential Contact*, SIGGRAPH 2020),
GPU-portable but a genuine project, NOT a kernel swap. This is the definitive scoping for PI: do not spend
more effort on penalty/barrier kernel variants; the path is a proper IPC (or a constraint-based contact
like SimuCell3D's own relocation-to-average with a hard non-penetration constraint).

## ⚠️ Over-claim correction (night-deliverables verification) — CHEAPER rungs exist before full IPC
The 2026-06-23 verification flagged "the fix needs full IPC" as PREMATURE — the prototypes ruled out only
specific simple forms, and skipped cheaper rungs the analysis points to. Concretely, the verifier found:
- **The contact Hessian is currently DISCARDED.** Contact IS launched inside `stiff_force_into` (it
  nominally goes through the implicit operator), but `device_cg`'s `pAp_diag` floor (implicit.py:183-187)
  THROWS AWAY the finite-difference contact Hessian whenever it goes spuriously negative — exactly what a
  gated/discontinuous contact force does near the c_rep boundary. So **the contact is effectively EXPLICIT
  right now** (its stiffness never resists the bundle implicitly). Wiring an **ANALYTIC contact+cadherin
  diagonal Hessian** into the operator (a stub estimator already exists at implicit.py:395) would give the
  IPC barrier-Hessian benefit **without CCD/line-search** — a genuinely cheaper rung not tested.
- **The ×40 cadherin spring is mislabeled "soft."** It is a 1–2 nN spring in the EXPLICIT RHS; moving it
  into the implicit operator (it is stiff, not soft) may let the implicit solve absorb it.
- **Higher cap on the BOUNDED saturating barrier (untested).** The rep-sweep raised rep on the UNsaturated
  penalty (diverges by stiffening the shallow shell); it never raised the CAP on the bounded saturating
  nearest-face barrier — which the verifier hoped "cannot diverge by construction." **TESTED (2026-06-23):
  rep 2e8 → stable but (converged) too weak; rep 2e9 → vv0=1.73 (DIVERGES).** So the bounded saturating
  cap ALSO diverges at a high-enough cap — the "can't diverge" hope is REFUTED. **But this is the smoking
  gun for the real diagnosis:** a *bounded* force diverging means the instability is NOT the force
  magnitude per se — it is that the contact lives in the EXPLICIT RHS / its Hessian is discarded by the
  pAp_diag floor, so the implicit solve never absorbs it. → the analytic-contact-Hessian-in-the-implicit-
  operator rung (next bullet) is the one the evidence points to.
- **Soft surface-constraint spring** inside `stiff_force_into` (not a hard post-step relocation) avoids the
  energy-injection that killed prototype #6, and is cheaper than a full KKT Lagrange solve.

So the honest standing: **the 6 prototypes rule out the SIMPLE penalty/barrier/projection forms, but the
correct sequence before committing to full IPC is: (1) analytic contact Hessian in the implicit operator,
(2) cadherin spring → implicit, (3) higher bounded-saturating cap — all cheaper, none yet fully tested.**
Full IPC (CCD+barrier+Newton) remains the guaranteed-correct fallback, not the proven-necessary first step.

## Recommendation (honest, corrected)
**The 6 simple prototypes are refuted, but full IPC was an OVER-CLAIM as the immediate fix — test the 3
cheaper rungs first (analytic contact Hessian in the implicit operator is the most promising).** The M1 fix is harder than any of
the three sketched options alone: it requires Option A's inside test + a bounded (non-linear, saturating)
barrier + implicit treatment. Until then the per-face penalty (with its tunnelling) is the least-bad
option, and the standing fact is: **at the lit cohesion strength the contact either tunnels (capped) or
diverges (uncapped) — aggregation/spreading are blocked on a proper barrier contact.** Surface to PI as a
scoped contact-mechanics task (reference: Incremental Potential Contact / IPC barrier methods). The
prototype + this negative result save the PI from the nearest-face dead end. **Design + negative
prototype, NOT a shipped fix.**
