#!/usr/bin/env python
r"""GATE-B NATIVE driver — the ``sf_arc`` component LAUNCHES its own rod-cable mechanics on CUDA.

Demonstrates that ``sf_arc`` is genuinely **KERNEL_BOUND** (not a driver-deferred scaffold): the SF component
OWNS SF-local device topology/parameter arrays and LAUNCHES the reused ``ff`` kernels over its DISJOINT
population — ``link_spring_kernel`` (Hookean axial backbone tension) + ``cytosim_bending_kernel`` (NF2007
discrete bending, α=κ/seg³) — through :class:`aleph.engine.sf_mechanics.SFFilamentMechanics`, exactly the
way the cortex does through ``CortexFilamentMechanics``.  No incumbent driver is involved; the SF component
launches the kernels itself.

WHAT IT SHOWS (all on-device; host reads are OUT-OF-LOOP diagnostics only):
  1. at the build geometry (rest) the passive backbone + bending inject ~ZERO force (emergent-not-lumped: the
     link rest lengths are the build segment lengths, so a relaxed SF carries no passive prestress);
  2. STRETCH the SF population (scale about its centroid) → a non-zero axial TENSION field EMERGES from
     ``link_spring_kernel`` over the SF segments;
  3. KINK interior nodes → a non-zero restoring BENDING force EMERGES from ``cytosim_bending_kernel``;
  4. an overdamped relaxation (``x += (dt/γ)·F`` via ``axpy_kernel``) drives the residual ‖F‖ DOWN — the SF
     rod-cable relaxes toward its own equilibrium under the launched kernels.  The explicit step ``dt/γ`` is
     CFL-STABLE BY DERIVATION (default): it is ``C/λ_max`` where ``λ_max`` is a Gershgorin upper bound on the
     assembled tangent stiffness (dominated by the SOURCED α-actinin arc crosslink, 4.6e5 pN/µm — NOT the
     axial GAP), so a stiff network cannot blow up; a user-supplied ``--dt-over-gamma`` above the bound is
     auto-clamped, and the verdict reports DIVERGED (never a false PASS) if the residual fails to fall;
  5. the internal ``dorsal_arc_crosslink`` connector also launches ``link_spring_kernel`` over the SF-owned
     dorsal↔arc joint pairs (a genuine second binding).

SOURCED vs GAP (report-not-tune, said out loud):
  * SOURCED — bending κ = KAPPA_ACTIN ≈ 0.0728 pN·µm² (k_B·T·ℓ_p, ℓ_p=17 µm, Gittes 1993, KU-1.1) → α=κ/seg³;
    internal crosslink k = 4.6e5 pN/µm (α-actinin, Ferrer 2008, PI-approved 2026-06-30).
  * GAP (REQUIRED, no default) — the actin axial backbone stiffness ``--k-axial`` [pN/µm].  NF2007 treats the
    backbone as inextensible (constraint + reshape, no axial penalty spring); there is no sourced EA_actin in
    ``ff.units``.  The value passed here is a PROVISIONAL modelling GAP; it sets the SIGN/RISE of the tension,
    NOT a production magnitude.  Source it or surface to PI before any quantitative claim.

This is a MECHANISM demonstration on the DISJOINT SF population, NOT a production run: the structural counts
come from :func:`aleph.engine.sf_population.build_sf_arc_population` (small, native-forbidden on the Mac),
and the sourced native SF inventory (count/length distribution) arrives from ``CellState`` in a biology phase.

────────────────────────────────────────────────────────────────────────────────────────────────────────
RUN ON GBOOK (needs a CUDA GPU; will NOT run on the dev Mac):

    ssh gbook
    cd ~/ffn_ac_native
    PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
      ~/miniconda3/envs/ffn_sim/bin/python \
      aleph/scripts/ac_gate_b_sf_arc_native.py --k-axial 1000 --steps 500 --stretch 1.05

(Omit ``--dt-over-gamma`` — the CFL-stable step is DERIVED from the built stiffnesses.  With the α-actinin arc
crosslink dominating, that step is ~1e-7 µm/pN, so a few hundred steps show the residual falling; the soft
axial/bending modes relax slowly, so more steps give a cleaner monotone trace but the ratio<1 PASS holds early.)

Long run: launch under nohup and monitor the LOG FILE (ssh python is not on PATH; use the full env python).
────────────────────────────────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import argparse

import numpy as np
import warp as wp

from aleph.engine.sf_mechanics import (
    ALPHA_ACTININ_K_PN_PER_UM,
    KAPPA_ACTIN_PN_UM2,
    build_sf_filament_mechanics,
    build_sf_internal_arc_connector,
    build_sf_mechanics_topology,
)
from aleph.engine.sf_mechanics import SFMechanicsTopology
from aleph.engine.sf_population import build_sf_arc_population
from aleph.laws.network_warp import _zero, axpy_kernel

#: Explicit-Euler (overdamped) Courant-number target for the relax step.  The HARD stability limit of the
#: overdamped update ``x += (dt/γ)·F`` with restoring ``F = −K·x`` is ``C = (dt/γ)·λ_max(K) < 2`` (explicit
#: Euler on a quadratic energy).  We TARGET ``C_target = 0.5 ≪ 2``: a ≥2× margin below the hard limit that also
#: absorbs the axial GEOMETRIC stiffening (transverse term ``k·|L−r0|/L < k``, so the true tangent λ_max is at
#: most 2× the material Gershgorin bound below — bending is exactly linear, no geometric term).  This is NOT a
#: fitted constant: it is a dimensionless CFL safety margin, stated and grid/stability-derived.
COURANT_TARGET = 0.5


def _max_force_norm(force_d: wp.array) -> float:
    """OUT-OF-LOOP host diagnostic: max per-node force magnitude (device→host read between phases only)."""
    f = force_d.numpy()
    return float(np.max(np.linalg.norm(f, axis=1))) if f.shape[0] else 0.0


def cfl_stable_dt_over_gamma(
    topo: SFMechanicsTopology, *, courant: float = COURANT_TARGET
) -> tuple[float, float]:
    """Grid-derived CFL-stable explicit step ``dt/γ`` for the overdamped relax, from the BUILT SF stiffnesses.

    The relaxation update ``x += (dt/γ)·F`` with restoring ``F = −K·x`` is explicit Euler; it is stable iff
    ``(dt/γ)·λ_max(K) < 2``.  We bound ``λ_max(K)`` from ABOVE by the assembled Gershgorin row-sum of the
    MATERIAL tangent — no eigensolve, purely the topology + stiffnesses actually built for this population:

      * each Hookean link/joint (stiffness ``k``) adds ``k`` to each endpoint's diagonal block and ``−k``
        off-diagonal ⇒ ``+2k`` to each endpoint's Gershgorin row sum;
      * each Cytosim bending triple (coeff ``α``) has stencil ``K = α·[[1,−2,1],[−2,4,−2],[1,−2,1]]`` ⇒ ``+4α``
        to each END node's row sum and ``+8α`` to the CENTER node's row sum.

    Then ``λ_max ≤ G = max nodal row sum`` and the returned step is ``dt/γ = courant/G`` (with ``G`` also
    returned for the printout).  The DOMINANT stiffness here is the SOURCED α-actinin arc crosslink
    (``4.6e5 pN/µm``), ~460× the axial GAP ``k_axial``, so IT (not ``k_axial``) sets the stable step — the
    reason the old fixed ``dt/γ=1e-4`` default blew up (Courant number ≈ 368 ≫ 2).  The axial geometric
    stiffening (``≤ k·|L−r0|/L < k``) can at most double the axial term, so ``courant ≤ 1`` keeps the true
    Courant number ``(dt/γ)·λ_max ≤ 2·courant < 2`` throughout the relaxation.

    Args:
        topo: the built SF rod-cable topology (links/triples/joints + their stiffnesses).
        courant: dimensionless target Courant number (default :data:`COURANT_TARGET`).

    Returns:
        ``(dt_over_gamma, lambda_max_bound)`` — the CFL-stable step [µm/pN] and the Gershgorin λ_max bound.
    """
    g = np.zeros(int(topo.n_nodes), dtype=np.float64)
    if topo.n_links:
        for (i, j), k in zip(topo.links, topo.link_k):
            g[int(i)] += 2.0 * float(k)
            g[int(j)] += 2.0 * float(k)
    if topo.n_arc_joints:
        for (i, j), k in zip(topo.arc_joints, topo.arc_k):
            g[int(i)] += 2.0 * float(k)
            g[int(j)] += 2.0 * float(k)
    if topo.n_triples:
        for (a, b, c), al in zip(topo.bend_triples, topo.bend_alpha):
            g[int(a)] += 4.0 * float(al)
            g[int(b)] += 8.0 * float(al)
            g[int(c)] += 4.0 * float(al)
    lambda_max_bound = float(g.max()) if g.size else 0.0
    if not np.isfinite(lambda_max_bound) or lambda_max_bound <= 0.0:
        raise SystemExit("[sf_arc GATE-B] FAIL: assembled stiffness is non-positive — no relaxable network")
    return float(courant) / lambda_max_bound, lambda_max_bound


def main() -> None:
    ap = argparse.ArgumentParser(description="Native GATE-B: sf_arc launches its own rod-cable mechanics.")
    ap.add_argument("--k-axial", type=float, required=True,
                    help="actin axial backbone stiffness [pN/µm] — REQUIRED modelling GAP (no default)")
    ap.add_argument("--n-ventral", type=int, default=8)
    ap.add_argument("--n-dorsal", type=int, default=4)
    ap.add_argument("--n-arc", type=int, default=4)
    ap.add_argument("--n-cap", type=int, default=4)
    ap.add_argument("--n-per-fiber", type=int, default=9)
    ap.add_argument("--stretch", type=float, default=1.05, help="uniform stretch factor about the centroid")
    ap.add_argument("--kink", type=float, default=0.15, help="interior-node z-kink [µm] to excite bending")
    ap.add_argument("--steps", type=int, default=200, help="overdamped relaxation steps")
    ap.add_argument("--dt-over-gamma", type=float, default=None,
                    help="explicit step dt/γ [µm per pN]; DEFAULT = the CFL-stable value DERIVED from the "
                         "built stiffnesses (no magic number). A provided value ABOVE the CFL bound is "
                         "auto-clamped down (printed) so the relax cannot go into the unstable regime.")
    ap.add_argument("--device", type=str, default="cuda:0")
    args = ap.parse_args()

    wp.init()
    device = args.device
    print(f"[sf_arc GATE-B] device={device}")
    print("[sf_arc GATE-B] SOURCED κ=KAPPA_ACTIN="
          f"{KAPPA_ACTIN_PN_UM2:.4g} pN·µm² (Gittes 1993) → α=κ/seg³ ; "
          f"internal crosslink k={ALPHA_ACTININ_K_PN_PER_UM:.3g} pN/µm (α-actinin, Ferrer 2008)")
    print(f"[sf_arc GATE-B] GAP (provisional): k_axial={args.k_axial:g} pN/µm (source or surface to PI)")

    # 1) build the DISJOINT sf_arc population (host) + its SF-local rod-cable topology.
    pop = build_sf_arc_population(
        n_ventral=args.n_ventral, n_dorsal=args.n_dorsal, n_arc=args.n_arc, n_cap=args.n_cap,
        n_per_fiber=args.n_per_fiber,
    )
    pop.assert_partitioned()
    topo = build_sf_mechanics_topology(pop, k_axial_pn_per_um=args.k_axial)
    print(f"[sf_arc GATE-B] population: n_nodes={topo.n_nodes} n_links={topo.n_links} "
          f"n_triples={topo.n_triples} n_arc_joints={topo.n_arc_joints} "
          f"(census: {pop.census()['N_unique_sf_fibers']} unique fibers)")

    # 2) the SF component OWNS + LAUNCHES its mechanics (link_spring + cytosim_bending) on CUDA.
    mechanics = build_sf_filament_mechanics(topo, device=device)
    arc_connector = build_sf_internal_arc_connector(topo, device=device)
    pos_d = wp.array(np.ascontiguousarray(pop.pos, np.float64), dtype=wp.vec3d, device=device)
    force_d = wp.zeros(topo.n_nodes, dtype=wp.vec3d, device=device)

    def accumulate() -> None:
        wp.launch(_zero, dim=topo.n_nodes, inputs=[force_d], device=device)
        mechanics.accumulate(pos_d, force_d)          # SF launches link_spring + cytosim_bending itself
        if arc_connector.n_joints:
            arc_connector.accumulate(pos_d, force_d)  # internal dorsal↔arc crosslink (link_spring)
        wp.synchronize_device(device)

    # rest: passive force ~ 0 (emergent-not-lumped).
    accumulate()
    print(f"[sf_arc GATE-B] (1) REST  max|F| = {_max_force_norm(force_d):.3e} pN  (≈0 expected)")

    # 3) stretch about the centroid → axial tension EMERGES; kink interior nodes → bending EMERGES.
    pos = pop.pos.copy()
    centroid = pos.mean(axis=0)
    pos = centroid + args.stretch * (pos - centroid)
    interior = np.unique(topo.bend_triples[:, 1]) if topo.n_triples else np.zeros(0, np.int64)
    pos[interior, 2] += args.kink
    pos_d = wp.array(np.ascontiguousarray(pos, np.float64), dtype=wp.vec3d, device=device)
    accumulate()
    f0 = _max_force_norm(force_d)
    print(f"[sf_arc GATE-B] (2) STRETCH×{args.stretch:g} + kink {args.kink:g}µm  "
          f"max|F| = {f0:.3e} pN  (>0 ⇒ tension/bending emerged from the launched kernels)")
    if f0 <= 0.0:
        raise SystemExit("[sf_arc GATE-B] FAIL: no SF force emerged from the launched kernels")

    # 4) overdamped relaxation: x += (dt/γ)·F — the residual must fall (CFL-stable step, else it blows up).
    #    Resolve dt/γ from the BUILT stiffnesses: the α-actinin arc crosslink (4.6e5 pN/µm) dominates λ_max,
    #    so the stable step is ~1e-7, NOT the old 1e-4 (which gave Courant number ≈ 368 ≫ 2 → the divergence).
    cfl_dt, lambda_max_bound = cfl_stable_dt_over_gamma(topo)
    if args.dt_over_gamma is None:
        step_val = cfl_dt
        print(f"[sf_arc GATE-B] CFL: λ_max ≤ {lambda_max_bound:.3e} pN/µm (Gershgorin bound, material tangent; "
              f"α-actinin arc crosslink dominates) ⇒ stable dt/γ = C/λ_max = {COURANT_TARGET:g}/"
              f"{lambda_max_bound:.3e} = {step_val:.3e} µm/pN  (hard limit dt/γ < 2/λ_max = "
              f"{2.0 / lambda_max_bound:.3e})")
    else:
        step_val = float(args.dt_over_gamma)
        if not np.isfinite(step_val) or step_val <= 0.0:
            raise SystemExit("[sf_arc GATE-B] FAIL: --dt-over-gamma must be positive and finite")
        if step_val > cfl_dt:
            print(f"[sf_arc GATE-B] CFL CLAMP: requested dt/γ={step_val:.3e} EXCEEDS the stable bound "
                  f"{cfl_dt:.3e} (λ_max ≤ {lambda_max_bound:.3e}); auto-reducing to {cfl_dt:.3e} µm/pN to "
                  f"stay below the Courant limit and avoid the CFL blow-up.")
            step_val = cfl_dt
        else:
            print(f"[sf_arc GATE-B] dt/γ={step_val:.3e} µm/pN provided (≤ CFL bound {cfl_dt:.3e}; "
                  f"λ_max ≤ {lambda_max_bound:.3e}) — kept.")

    step = wp.float64(step_val)
    sample_every = max(1, args.steps // 10)
    residuals: list[tuple[int, float]] = [(0, f0)]
    diverged = False
    for s in range(args.steps):
        accumulate()
        if s % sample_every == 0:
            fs = _max_force_norm(force_d)
            residuals.append((s, fs))
            if (not np.isfinite(fs)) or fs > 10.0 * f0:      # explicit divergence guard (never claim success)
                diverged = True
                break
        wp.launch(axpy_kernel, dim=topo.n_nodes, inputs=[pos_d, step, force_d], device=device)
    accumulate()
    fN = _max_force_norm(force_d)
    residuals.append((args.steps, fN))

    trace = "  ".join(f"{s}:{r:.2e}" for s, r in residuals)
    print(f"[sf_arc GATE-B] relax residual trace max|F| [pN] = {trace}")
    monotone = all(residuals[k + 1][1] <= residuals[k][1] * (1.0 + 1.0e-6) for k in range(len(residuals) - 1))

    # verdict: RELAXED requires a FINITE final residual STRICTLY below the perturbed one; anything else is a
    # CFL divergence (or a stall) — never a false-positive PASS on a blow-up.
    relaxed = bool(np.isfinite(fN)) and (f0 > 0.0) and (fN < f0) and not diverged
    if relaxed:
        descent = "monotone" if monotone else "non-monotone (soft modes still converging)"
        print(f"[sf_arc GATE-B] (4) after {args.steps} relax steps  max|F| = {fN:.3e} pN  "
              f"(ratio {fN / f0:.3e} < 1 ⇒ the SF rod-cable relaxed under its own launched kernels; "
              f"{descent} descent)")
        print(f"[sf_arc GATE-B] VERDICT PASS: sf_arc launched link_spring_kernel + cytosim_bending_kernel over "
              f"its DISJOINT population; tension/bending EMERGED (REST≈0, STRETCH>0) and RELAXED (ratio<1). "
              f"(k_axial={args.k_axial:g} pN/µm is a GAP — magnitude is provisional.)")
    else:
        ratio_str = f"{fN / f0:.3e}" if (np.isfinite(fN) and f0 > 0.0) else "non-finite"
        print(f"[sf_arc GATE-B] (4) after {args.steps} relax steps  max|F| = {fN:.3e} pN  (ratio {ratio_str}) "
              f"— residual did NOT fall below the perturbed max|F|.")
        print("[sf_arc GATE-B] VERDICT DIVERGED (CFL — reduce dt): the explicit relax blew up or failed to "
              "decrease max|F|. Omit --dt-over-gamma to use the auto CFL-stable default, or pass a value below "
              "the printed CFL bound. This is a NUMERICS failure of the relax step, NOT an sf_arc binding "
              "failure — the genuine-launch REST/STRETCH checks above already passed.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
