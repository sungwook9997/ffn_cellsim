# Codex (GPT) solver-optimization review — Lead's judgment (2026-07-25)

PI asked for an independent Codex consult on the "optimization system" (the GATE-A /
fine-mesh projected-CG + fiber-arclength multigrid solver). Codex ran read-only (90 s,
no stall) and returned a technically strong 7-point review. This records it + the
Lead's judgment. **Codex = debug/idea input, not verification authority** (project rule).

## The key judgment: Codex's premise is superseded by the Lead's root-cause finding

Codex reviewed the solver IN ISOLATION and framed the fine-mesh plateau as a solver
globalization problem to be beaten (trust-region SQP, better MG). **But the Lead had
already proven by production that the plateau is a GEOMETRY artifact — the
`overlap_free` WCA relaxation transverse-kink — not a solver problem:** building with
smooth arcs (`--no-overlap-free`) collapses the fine-75nm plateau 3.42 → 0.13 with the
SAME solver (commits 21a56220 / 3a1f8866). So the solver-optimization target Codex
aimed at is largely moot — the fix was upstream in the cortex construction. Codex could
not see this (it only read the solver). This is exactly the depth difference: the Lead
localized the real cause by measurement + production; Codex optimized the wrong target.

## Codex's 7 points + Lead's verdict

1. **Trust-region/SQP Newton instead of max-force line search; retract trials onto the
   manifold before evaluating.** Reasonable general robustness idea. NOT the plateau fix
   (kink is). Low priority now.
2. **Trial-state inconsistency: `_accumulate_all(cell,pos,f)` evaluates bending/xl/
   branch/myosin/nucleus/membrane at `pos` but steric+pressure at `cell.state`
   (driver.py:164-167).** ✅ REAL code detail, correctly spotted. Lead's read: this is
   likely a DELIBERATE operator-split (the Biot pressure FIELD is solved on the outer
   physical clock and held frozen during the inner mechanical solve — consistent with
   the architecture; steric contact-set frozen per linearization, which Codex itself
   endorses in point 6). Worth confirming the intended semantics, but **does not affect
   any plateau finding**: the force decomposition was measured at the committed plateau
   state (pos == cell.state → consistent), and the root cause was confirmed by an
   independent build-test (smooth arcs). Follow-up: confirm/annotate the frozen-field
   semantics in the probe's line search.
3. **Add the constraint geometric stiffness (DP)[s]F (SQP reduced Hessian).** Legit
   higher-order Newton term; solver-quality, not plateau cause. Queue.
4. **Cache the line-smoother factorization once per level; exact pentadiagonal LDLᵀ
   instead of row-sum-lumped tridiagonal (fiber_arclength_mg.py:259).** ✅ A real MG
   efficiency win — BUT the MG's purpose (beat the plateau) is now moot, so low priority.
5. **`assembled_coarse` drops the projector + inter-fiber off-diagonals; not a true
   Galerkin (fiber_arclength_mg.py:352,403); the SPD/h-indep test is crosslink-free +
   unprojected.** ✅ Fair — this confirms the Lead's OWN flag on FIX2 ("drops inter-fiber,
   needs a native reduction A/B"). Valid test-coverage critique. Queue if MG is revived.
6. **Cache the WCA contact pair list per linearization (hash-grid query runs every
   apply, 4×/V-cycle); replace FP64 atomics in `_dot_kernel` with block/warp reduction
   (implicit_mechanics.py:347,1329).** ✅ Real GPU-perf wins, generally applicable.
   Queue (perf, not correctness).
7. **`a` (global-max omitted stiffness incl κ/h³) can overdamp + reintroduce mesh
   dependence; the MG omega uses a power-iteration Rayleigh quotient (→λmax from BELOW)
   as if an upper bound (fiber_arclength_mg.py:719) — use Gershgorin/Lanczos + Chebyshev.**
   ⚠️ The omega-upper-bound point is a subtle correctness concern worth a look if the MG
   is used in production; the `a`-overdamping/mesh-dependence Lead traced to the kink, not
   `a`. Queue.

## Disposition

Net: a strong review of a now-lower-value subsystem. **No urgent action** — the fine mesh
converges with the kink fix + the existing MG, so MG perf/robustness (4,5,6,7) is
deferred. Point 2 gets a semantics confirmation note. The Lead's priority stays the
full-cell breadth (interior column, sf_arc) per the PI. Codex added genuine solver-detail
value (points 2,4,6) but missed the actual root cause the Lead had already proven.
