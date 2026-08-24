#!/usr/bin/env python
r"""Is the bound SF-motor null space EXACTLY the bound crossbridges' transverse blindness?

WHY THIS EXISTS.  ``ac_observe_sf_operator_native.py`` measured, at 156 DOF, a bound-configuration
tangent with 40 eigenvalues at the machine floor where the bond graph predicts ZERO.  The committed
M-H result says a bound crossbridge has no transverse channel at all, which predicts exactly two zero
directions per bound head.  20 bound heads x 2 = 40 is an arithmetic coincidence until the null
VECTORS are shown to be those directions, and that is what this measures.

It is not a re-derivation: the operator is assembled by the same helpers the observe driver uses, so
the spectrum reproduced here must match the committed run or the setup diverged (checked, and the run
aborts if it does).

WHAT IS DECIDED BY IT.  The arena's array list.  If the null space IS the transverse blindness, it is a
measured physical absence and belongs in a CONSTRAINT row (the NF2007 projector this engine already
uses for backbone inextensibility), not a projection basis and not a new per-head stiffness.

Sanity Gate: no constant is introduced.  The zero threshold is the observe module's own
``n_dof * eps64 * ||K||`` machine floor, not a chosen tolerance.  The predicted basis is built from
geometry alone (the bound segment direction), never fitted.  CUDA-only.  Every host read is after the
operator is assembled; nothing is read inside a solve.

Run:
    PYTHONPATH=$HOME/ffn/ffn_cellsim ~/miniforge3/envs/ffn_sim/bin/python nullspace_probe.py --out <path>
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp

# The observe driver's own helpers, so the operator is assembled by identical code.
import aleph
import aleph.scripts.ac_observe_sf_operator_native as obs
from aleph.components.motor.minifilament_topology import MinifilamentTopology
from aleph.components.motor.segment_motor import SegmentDetachKinetics
from aleph.engine.observe import assemble_dense_operator, stiffness_spectrum
from aleph.engine.sf_mechanics import build_sf_mechanics_topology
from aleph.engine.sf_motor_slice import build_sf_motor_slice, sf_sarcomere_geometry_for
from aleph.engine.sf_population import build_sf_arc_population

# The committed observe-run values this probe must reproduce before its own measurement means anything.
EXPECT_N_ZERO_BOUND = 40
EXPECT_COND_BOUND = 2.284e12
EXPECT_N_DOF = 156


def main() -> None:
    ap = argparse.ArgumentParser(description="Localise the bound-configuration null space onto heads.")
    ap.add_argument("--device", type=str, default="cuda")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--bind-steps", type=int, default=10)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--newton", type=int, default=4)
    ap.add_argument("--cg-iterations", type=int, default=200)
    ap.add_argument("--build-commit", type=str, default="")
    ap.add_argument("--out", type=str, default="")
    a = ap.parse_args()

    t_start = time.perf_counter()
    wp.init()
    resolved = wp.get_device(a.device)
    if not resolved.is_cuda:
        raise SystemExit("I0-A: CUDA only; there is no CPU simulation path in this repo")
    device = str(resolved)
    print(f"[null] device={device}  aleph={aleph.__file__}  warp={wp.config.version}")

    # ── the SAME slice the observe run measured ──────────────────────────────────────────────────
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
    n_nmii = int(actuator.position_d.shape[0])
    n = n_sf + n_nmii
    fa_sites = np.asarray(population.fa_sites, np.int64)
    print(f"[null] {n_sf} SF + {n_nmii} NMII = {n} nodes ({3 * n} DOF)")

    # ── drive the same accepted steps so the heads bind ──────────────────────────────────────────
    from aleph.engine.sf_implicit import SFImplicitCG, SFImplicitTopology
    implicit_topology = SFImplicitTopology(
        device=device, n_sf=n_sf, n_nmii=n_nmii, sf_topology=sf_topology, actuator=actuator,
        connector_state=connector.state, k_xb_pn_per_um=1000.0,
        r0_xb_um=float(params.r0_xb), pinned_sf_nodes=fa_sites)
    implicit_solver = SFImplicitCG(implicit_topology, max_iterations=int(a.cg_iterations))
    lambda_g = float(np.max(np.abs(sf_topology.link_k))) if sf_topology.n_links else 0.0
    mobility = 1.0 / max(1.0, float(np.max(np.abs(sf_topology.arc_k)))
                         if sf_topology.n_arc_joints else lambda_g)

    def inner_solve() -> None:
        for _ in range(int(a.newton)):
            slice_.zero_forces()
            slice_.accumulate()
            implicit_solver.step(
                sf_position_d=sf_owner.position_d, sf_force_d=sf_owner.force_d,
                nmii_position_d=actuator.position_d, nmii_force_d=actuator.force_d,
                mobility_step=mobility)
        slice_.zero_forces()
        slice_.accumulate()

    accepted_d = wp.array(np.ones(1, np.int32), dtype=wp.int32, device=device)
    for _ in range(int(a.bind_steps)):
        slice_.step(inner_solve, dt_phys=float(a.dt), accepted_d=accepted_d)
        wp.synchronize_device(device)
    n_bound = int(slice_.bound_head_count())
    print(f"[null] bound heads: {n_bound}/{int(actuator.n_heads)}")

    # ── assemble the tangent with the observe driver's OWN helpers ───────────────────────────────
    positions = obs._concatenated_positions(sf_owner, actuator)
    solver = obs._make_variant_solver(
        device=device, base_sf_topology=sf_topology, actuator=actuator, connector=connector,
        sf_owner=sf_owner, k_xb=1000.0, r0_xb=float(params.r0_xb), fa_sites=fa_sites, scales={})
    with wp.ScopedDevice(device):
        positions_d = wp.array(np.ascontiguousarray(positions, np.float64), dtype=wp.vec3d)
    matrix = assemble_dense_operator(obs._make_matvec(solver, positions_d, device), int(solver.n))
    spectrum = stiffness_spectrum(matrix, rigid_body_zero_modes=0)
    n_zero = int(spectrum.n_zero_modes)
    cond = spectrum.condition_number
    print(f"[null] REPRODUCTION: n_dof={spectrum.n_dof} zero_modes={n_zero} "
          f"cond={cond:.4g} tol={spectrum.zero_tolerance:.4g}")

    # A probe whose operator does not match the committed run measures a different thing. Halt.
    if spectrum.n_dof != EXPECT_N_DOF or n_zero != EXPECT_N_ZERO_BOUND:
        raise SystemExit(
            f"REPRODUCTION FAILED: expected {EXPECT_N_DOF} DOF / {EXPECT_N_ZERO_BOUND} zero modes, "
            f"got {spectrum.n_dof} / {n_zero}. The setup diverged from the observe run; the "
            f"localisation below would not be about the same operator.")

    # ── the PREDICTION, built from geometry alone ────────────────────────────────────────────────
    # M-H: a bound head's crossbridge enters only through <d, w_hat>, so the two directions transverse
    # to the bound segment carry no stiffness. Predict 2 unit vectors per bound head, embedded at that
    # head's three DOF, orthonormalised against each other. Nothing here is fitted.
    state = connector.state
    bound = np.asarray(state.bound_d.numpy(), np.int64).reshape(-1)
    seg_a = np.asarray(state.seg_a_d.numpy(), np.int64).reshape(-1)
    seg_b = np.asarray(state.seg_b_d.numpy(), np.int64).reshape(-1)
    head_node = np.asarray(actuator.head_node_d.numpy(), np.int64).reshape(-1) + n_sf

    columns: list[np.ndarray] = []
    per_head = []
    for h in range(bound.size):
        if bound[h] == 0:
            continue
        w = positions[seg_b[h]] - positions[seg_a[h]]
        nw = float(np.linalg.norm(w))
        if nw <= 0.0:
            continue
        w = w / nw
        # any two unit vectors spanning the plane transverse to w
        seed_v = np.array([1.0, 0.0, 0.0]) if abs(w[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        t1 = np.cross(w, seed_v); t1 /= np.linalg.norm(t1)
        t2 = np.cross(w, t1); t2 /= np.linalg.norm(t2)
        base = 3 * int(head_node[h])
        for t in (t1, t2):
            col = np.zeros(3 * n, np.float64)
            col[base:base + 3] = t
            columns.append(col)
        per_head.append(int(head_node[h]))
    predicted = np.stack(columns, axis=1) if columns else np.zeros((3 * n, 0))
    print(f"[null] predicted basis: {predicted.shape[1]} columns "
          f"({len(per_head)} bound heads x 2 transverse)")

    # ── the test: is the predicted span INSIDE the measured null space? ──────────────────────────
    # Project each predicted direction onto the measured null space and report the residual. A
    # direction that lies in the null space projects onto itself; the residual norm is then ~0.
    null_basis = spectrum.eigenvectors[:, :n_zero]
    proj = null_basis @ (null_basis.T @ predicted)
    residual = np.linalg.norm(predicted - proj, axis=0)
    # Rayleigh quotient: the stiffness the operator actually assigns each predicted direction.
    rayleigh = np.einsum("ij,jk,ik->i", predicted.T, matrix, predicted.T)

    # Where does the measured null space LIVE? Fraction of each null vector's norm on head DOF.
    head_dof = np.zeros(3 * n, bool)
    for hn in per_head:
        head_dof[3 * hn:3 * hn + 3] = True
    weight_on_heads = np.sum(null_basis[head_dof, :] ** 2, axis=0)

    result = {
        "n_dof": int(spectrum.n_dof),
        "n_bound_heads": int(n_bound),
        "measured_zero_modes": n_zero,
        "zero_tolerance_pN_per_um": float(spectrum.zero_tolerance),
        "condition_number": float(cond) if cond is not None else None,
        "lambda_max_pN_per_um": float(spectrum.lambda_max),
        "predicted_basis_columns": int(predicted.shape[1]),
        "containment_residual_max": float(residual.max()) if residual.size else None,
        "containment_residual_median": float(np.median(residual)) if residual.size else None,
        "predicted_direction_stiffness_max_pN_per_um": float(np.max(np.abs(rayleigh)))
        if rayleigh.size else None,
        "null_weight_on_bound_head_dof_min": float(weight_on_heads.min()),
        "null_weight_on_bound_head_dof_median": float(np.median(weight_on_heads)),
        "null_weight_on_bound_head_dof_total_mean": float(weight_on_heads.mean()),
        "aleph_file": aleph.__file__,
        "warp_version": wp.config.version,
        "device": device,
        "build_commit": a.build_commit,
        "wall_seconds": time.perf_counter() - t_start,
    }
    print("[null] " + json.dumps(
        {k: v for k, v in result.items() if k.startswith(("containment", "predicted_dir", "null_weight"))},
        indent=1))



    # ── SUBTRACT THE BOUNDARY CONDITION BEFORE COUNTING AN ANOMALY ───────────────────────────────
    # The IMEX operator is (aI + M K M + a(I-M)) with M the Dirichlet mask, so MKM zeroes the ROWS
    # AND COLUMNS of every pinned DOF: 3 * n_pinned EXACT zero modes exist BY CONSTRUCTION and are
    # the boundary condition, not a finding. Their eigenvectors are the Cartesian basis at those DOF,
    # which costs nothing to build, so containment is testable with the same projector used above.
    pinned_basis = np.zeros((3 * n, 3 * fa_sites.size), np.float64)
    for j, node in enumerate(fa_sites):
        for axis in range(3):
            pinned_basis[3 * int(node) + axis, 3 * j + axis] = 1.0
    pin_resid = np.linalg.norm(
        pinned_basis - null_basis @ (null_basis.T @ pinned_basis), axis=0)
    n_mask = int(pinned_basis.shape[1])
    print(f"[null] pinned nodes={fa_sites.size} -> {n_mask} mask zero modes by construction; "
          f"containment residual max {pin_resid.max():.3e}")


    # QR in _rigid_body_basis MIXES translations with rotations, so the six cannot be read separately
    # from it. Build them apart: a translation must be an EXACT zero mode of any distance-based
    # tangent, whereas a rotation need not be away from equilibrium — a stressed configuration carries
    # geometric stiffness in rotation, which is correct physics, not a defect. Separating them is what
    # distinguishes "translation invariance is broken" from "this configuration is not relaxed".
    trans = np.zeros((3 * n, 3), np.float64)
    for axis in range(3):
        trans[axis::3, axis] = 1.0
    trans /= np.linalg.norm(trans, axis=0)
    rot_cols = []
    centred = positions - positions.mean(axis=0)
    for axis in range(3):
        unit = np.zeros(3); unit[axis] = 1.0
        rot_cols.append(np.cross(unit, centred).reshape(-1))
    rot = np.stack(rot_cols, axis=1)
    rot -= trans @ (trans.T @ rot)            # orthogonalise against translation
    rot /= np.linalg.norm(rot, axis=0)
    t_stiff = np.einsum("ji,jk,ki->i", trans, matrix, trans)
    r_stiff = np.einsum("ji,jk,ki->i", rot, matrix, rot)
    t_res = np.linalg.norm(trans - null_basis @ (null_basis.T @ trans), axis=0)
    r_res = np.linalg.norm(rot - null_basis @ (null_basis.T @ rot), axis=0)
    print(f"[null] TRANSLATION: |Rayleigh| max {np.max(np.abs(t_stiff)):.3e} pN/um  "
          f"(zero floor {spectrum.zero_tolerance:.3e})  containment residual max {t_res.max():.3e}")
    print(f"[null] ROTATION   : |Rayleigh| max {np.max(np.abs(r_stiff)):.3e} pN/um  "
          f"containment residual max {r_res.max():.3e}")
    result["translation_stiffness_max_pN_per_um"] = float(np.max(np.abs(t_stiff)))
    result["translation_containment_residual_max"] = float(t_res.max())
    result["rotation_stiffness_max_pN_per_um"] = float(np.max(np.abs(r_stiff)))
    result["rotation_containment_residual_max"] = float(r_res.max())
    result["translation_is_exact_zero_mode"] = bool(
        np.max(np.abs(t_stiff)) <= spectrum.zero_tolerance)

    # The mask basis is NOT null (residual ~1 above), and the reason is that _make_matvec applies
    # solver._stiffness directly, NOT solver._operator: the probe deliberately sees raw K, with no
    # Dirichlet mask and no aI. So there is no boundary condition in this operator to subtract --
    # and by the same token the FA anchors do NOT anchor it, so every connected component is FREE and
    # contributes its six rigid-body modes. That is what must be subtracted instead. Six per component.
    rb = obs._rigid_body_basis(positions)
    rb_q, _ = np.linalg.qr(rb)
    rb_resid = np.linalg.norm(rb_q - null_basis @ (null_basis.T @ rb_q), axis=0)
    rb_stiff = np.einsum("ji,jk,ki->i", rb_q, matrix, rb_q)
    print(f"[null] rigid-body basis of the WHOLE system: {rb_q.shape[1]} modes; "
          f"containment residual max {rb_resid.max():.3e}; "
          f"max |Rayleigh| {np.max(np.abs(rb_stiff)):.3e} pN/um")
    result["rigid_body_modes_tested"] = int(rb_q.shape[1])
    result["rigid_body_containment_residual_max"] = float(rb_resid.max())
    result["rigid_body_stiffness_max_pN_per_um"] = float(np.max(np.abs(rb_stiff)))
    result["operator_is_raw_K_unmasked"] = True
    result["anomaly_after_rigid_body"] = int(n_zero - rb_q.shape[1])
    print(f"[null] ANOMALY after subtracting rigid-body: {n_zero} - {rb_q.shape[1]} "
          f"= {n_zero - rb_q.shape[1]} of {3 * n} DOF")

    print(f"[null] ANOMALY after subtracting the boundary condition: "
          f"{n_zero} - {n_mask} = {n_zero - n_mask} of {3 * n} DOF")
    result["n_pinned_nodes"] = int(fa_sites.size)
    result["mask_zero_modes_by_construction"] = n_mask
    result["mask_containment_residual_max"] = float(pin_resid.max())
    result["anomalous_zero_modes_after_mask"] = int(n_zero - n_mask)

    # ── REFUTATION FOLLOW-UP: what ARE the 40 modes, if not per-head transverse? ──────────────────
    # Decompose each null vector's squared norm across the three particle classes the slice owns.
    head_idx = np.asarray(actuator.head_node_d.numpy(), np.int64).reshape(-1) + n_sf
    cls = np.empty(n, dtype="<U8"); cls[:] = "backbone"
    cls[:n_sf] = "sf"
    cls[head_idx] = "head"
    weights = {}
    for name in ("sf", "backbone", "head"):
        mask = np.zeros(3 * n, bool)
        for i in np.nonzero(cls == name)[0]:
            mask[3 * i:3 * i + 3] = True
        weights[name] = np.sum(null_basis[mask, :] ** 2, axis=0)
        print(f"[null] null-space weight on {name:9s}: "
              f"mean {weights[name].mean():.4f}  min {weights[name].min():.4f}  "
              f"max {weights[name].max():.4f}")
    result["null_weight_by_class_mean"] = {k: float(v.mean()) for k, v in weights.items()}
    result["null_weight_by_class_max"] = {k: float(v.max()) for k, v in weights.items()}

    # Do the modes decouple per minifilament, i.e. are they a single-object property that TRANSFERS?
    mf_off = np.asarray(actuator.minifilament_offset_d.numpy(), np.int64).reshape(-1)
    print(f"[null] minifilament offsets: {mf_off.tolist()}")
    nmii_mask = np.zeros(3 * n, bool); nmii_mask[3 * n_sf:] = True
    w_nmii = np.sum(null_basis[nmii_mask, :] ** 2, axis=0)
    result["null_weight_on_nmii_block_mean"] = float(w_nmii.mean())
    result["null_weight_on_nmii_block_min"] = float(w_nmii.min())
    print(f"[null] null-space weight on the WHOLE NMII block: mean {w_nmii.mean():.4f} "
          f"min {w_nmii.min():.4f}  -> {'confined to NMII' if w_nmii.min() > 0.99 else 'MIXED with sf'}")

    # The head-arm spring is what actually resists head transverse motion. Report it explicitly so the
    # refutation carries its cause, not just its verdict.
    result["predicted_transverse_stiffness_is_k_head_arm"] = {
        "measured_pN_per_um": float(np.max(np.abs(rayleigh))) if rayleigh.size else None,
        "k_head_arm_input_pN_per_um": 100.0,
        "note": "a bound head is NOT transversally free: the head-arm spring to the backbone carries "
                "it. M-H's 'no transverse channel' is a property of the CROSSBRIDGE BOND, not of the "
                "head's mobility, so 40 = 20x2 is arithmetic coincidence.",
    }

    verdict = (
        "CONFIRMED" if (residual.size and residual.max() < 1e-8
                        and predicted.shape[1] == n_zero) else
        "PARTIAL" if (residual.size and residual.max() < 1e-8) else "REFUTED")
    result["verdict"] = verdict
    print(f"[null] VERDICT: {verdict} "
          f"(predicted {predicted.shape[1]} of {n_zero} measured zero modes; "
          f"max containment residual {residual.max():.3e})")

    if a.out:
        with open(a.out, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"[null] wrote {a.out}")


if __name__ == "__main__":
    main()
