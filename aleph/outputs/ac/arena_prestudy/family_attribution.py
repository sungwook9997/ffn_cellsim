#!/usr/bin/env python
r"""WHICH stiffness family leaves the bound SF-motor operator rank-deficient?

WHY THIS EXISTS.  The bound tangent has 40 of 156 directions at the machine floor where the bond graph
predicts none, and the first hypothesis — that they are the 20 bound crossbridges' transverse blindness,
two per head — was REFUTED: those predicted directions carry 100.04 pN/um, which is ``k_head_arm``, so a
bound head is not transversally free at all.  The null space is collective (weight 0.52 head / 0.37
backbone / 0.11 sf, mixed rather than confined).  This asks the operator itself which family is
responsible, instead of guessing again.

Two questions per family, both answered by re-assembling the SAME operator with that family's scale
changed and counting zero modes at the SAME machine floor:

  * ONLY f  -> what does f on its own constrain?  (all other families scaled to 0)
  * WITHOUT f -> what does f uniquely contribute?  (f scaled to 0, others at 1)

A family whose removal does not raise the zero count constrains nothing the others do not already
constrain.  A null space that survives every family at full strength is a property of the MODEL, not of
a disabled term — which is what decides whether the arena needs a constraint row, a projection basis, or
neither.

Sanity Gate: no constant introduced; scales are exactly 0 or 1, never fitted.  The zero threshold is the
observe module's own ``n_dof*eps64*||K||``.  CUDA-only.  Reproduction of the 40/156 baseline is asserted
before any variant runs.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp

import aleph
import aleph.scripts.ac_observe_sf_operator_native as obs
from aleph.components.motor.minifilament_topology import MinifilamentTopology
from aleph.components.motor.segment_motor import SegmentDetachKinetics
from aleph.engine.observe import assemble_dense_operator, stiffness_spectrum
from aleph.engine.sf_mechanics import build_sf_mechanics_topology
from aleph.engine.sf_motor_slice import build_sf_motor_slice, sf_sarcomere_geometry_for
from aleph.engine.sf_population import build_sf_arc_population

FAMILIES = ("sf_axial", "sf_bending", "alpha_actinin_arc", "nmii_backbone",
            "nmii_head_arm", "nmii_angle_backbone", "nmii_angle_arm", "crossbridge")
EXPECT_N_DOF, EXPECT_N_ZERO = 156, 40


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", type=str, default="cuda")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--build-commit", type=str, default="")
    ap.add_argument("--out", type=str, default="")
    a = ap.parse_args()

    t0 = time.perf_counter()
    wp.init()
    resolved = wp.get_device(a.device)
    if not resolved.is_cuda:
        raise SystemExit("I0-A: CUDA only")
    device = str(resolved)
    print(f"[fam] device={device}  aleph={aleph.__file__}  warp={wp.config.version}")

    topology_nmii = MinifilamentTopology(
        n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.2)
    geometry = sf_sarcomere_geometry_for(topology_nmii)
    population = build_sf_arc_population(
        n_ventral=1, n_dorsal=0, n_arc=0, n_cap=0, n_per_fiber=9, **geometry)
    population.assert_partitioned()
    sf_topology = build_sf_mechanics_topology(population, k_axial_pn_per_um=1000.0)

    from aleph.engine.cortex_motor_slice import NMIIMotorParams
    params = NMIIMotorParams(
        v0=0.12, f_stall=0.5, kappa=0.5, k_xb=1000.0, r0_head=0.2, r0_xb=0.0,
        capture_radius=0.21, k_on=50.0, k_off0=0.35, f0=4.28e-3 / 0.6e-3,
        detach_kinetics=SegmentDetachKinetics.SLIP)
    slice_ = build_sf_motor_slice(
        sf_population=population, minifilament_topology=topology_nmii, params=params,
        sf_k_axial_pn_per_um=1000.0, nmii_k_backbone_pn_per_um=1000.0,
        nmii_k_head_spring_pn_per_um=100.0, nmii_backbone_persistence_length_um=1.0,
        base_seed=int(a.seed), device=device)
    sf_owner, actuator, connector = slice_.sf_owner, slice_.actuator_state, slice_.connector
    n_sf = int(sf_owner.n_nodes)
    fa_sites = np.asarray(population.fa_sites, np.int64)

    from aleph.engine.sf_implicit import SFImplicitCG, SFImplicitTopology
    it = SFImplicitTopology(
        device=device, n_sf=n_sf, n_nmii=int(actuator.position_d.shape[0]),
        sf_topology=sf_topology, actuator=actuator, connector_state=connector.state,
        k_xb_pn_per_um=1000.0, r0_xb_um=float(params.r0_xb), pinned_sf_nodes=fa_sites)
    solver_bind = SFImplicitCG(it, max_iterations=200)
    lam = float(np.max(np.abs(sf_topology.link_k)))
    mobility = 1.0 / max(1.0, lam)

    def inner_solve() -> None:
        for _ in range(4):
            slice_.zero_forces(); slice_.accumulate()
            solver_bind.step(sf_position_d=sf_owner.position_d, sf_force_d=sf_owner.force_d,
                             nmii_position_d=actuator.position_d, nmii_force_d=actuator.force_d,
                             mobility_step=mobility)
        slice_.zero_forces(); slice_.accumulate()

    accepted_d = wp.array(np.ones(1, np.int32), dtype=wp.int32, device=device)
    for _ in range(10):
        slice_.step(inner_solve, dt_phys=0.01, accepted_d=accepted_d)
        wp.synchronize_device(device)
    n_bound = int(slice_.bound_head_count())
    positions = obs._concatenated_positions(sf_owner, actuator)
    with wp.ScopedDevice(device):
        positions_d = wp.array(np.ascontiguousarray(positions, np.float64), dtype=wp.vec3d)

    def measure(scales: dict) -> tuple[int, float, float]:
        s = obs._make_variant_solver(
            device=device, base_sf_topology=sf_topology, actuator=actuator, connector=connector,
            sf_owner=sf_owner, k_xb=1000.0, r0_xb=float(params.r0_xb), fa_sites=fa_sites,
            scales=scales)
        m = assemble_dense_operator(obs._make_matvec(s, positions_d, device), int(s.n))
        sp = stiffness_spectrum(m, rigid_body_zero_modes=0)
        return int(sp.n_zero_modes), float(sp.lambda_max), (sp.condition_number or float("nan"))

    base_zero, base_lmax, base_cond = measure({})
    print(f"[fam] BASELINE all families: zero={base_zero}/{EXPECT_N_DOF} "
          f"lambda_max={base_lmax:.6g} cond={base_cond:.4g}  (bound heads {n_bound})")
    if base_zero != EXPECT_N_ZERO:
        raise SystemExit(f"REPRODUCTION FAILED: expected {EXPECT_N_ZERO} zero modes, got {base_zero}")

    rows = []
    for f in FAMILIES:
        only = {g: (1.0 if g == f else 0.0) for g in FAMILIES}
        without = {f: 0.0}
        z_only, lmax_only, _ = measure(only)
        z_wo, _, cond_wo = measure(without)
        rows.append({"family": f, "zero_only": z_only, "lambda_max_only": lmax_only,
                     "zero_without": z_wo, "delta_vs_baseline": z_wo - base_zero,
                     "cond_without": cond_wo})
        print(f"[fam] {f:22s} ONLY-> zero {z_only:3d} (lmax {lmax_only:12.6g})   "
              f"WITHOUT-> zero {z_wo:3d}  delta {z_wo - base_zero:+d}")

    result = {
        "n_dof": EXPECT_N_DOF, "n_bound_heads": n_bound,
        "baseline_zero_modes": base_zero, "baseline_lambda_max_pN_per_um": base_lmax,
        "baseline_condition_number": base_cond, "families": rows,
        "aleph_file": aleph.__file__, "warp_version": wp.config.version, "device": device,
        "build_commit": a.build_commit, "wall_seconds": time.perf_counter() - t0,
    }
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"[fam] wrote {a.out}")


if __name__ == "__main__":
    main()
