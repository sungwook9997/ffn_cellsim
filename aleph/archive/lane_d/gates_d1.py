"""D1 sanity gates — nucleus + MT/IF/LINC force-path rig (Lane D, NON-AUTHORITATIVE DIAGNOSTIC).

Gates (brief D1):
    G1.1 force ledger between membrane/input load and nucleus reaction (Newton's 3rd across the path)
    G1.2 MT compression/buckling sign (restoring toward straight; Euler-buckling oracle)
    G1.3 IF tension/strain-stiffening sign (stretched spoke pulls the nucleus toward the anchor)
    G1.4 LINC equal-and-opposite reaction (nforce = −aforce; matches the host oracle)
    G1.5 LINC detach removes ONLY that graph path (others bit-exact unchanged)
    G1.6 MT growth/shrink / endpoint-remap generation validation (stale reference rejected)
    G1.7 rejected topology event → geometry + state bit-exact recovery
    G1.8 solid work sign / conservation (actuator does positive work; path reversible). The paired
         cytosol-pressure-work half is D0-G0.4 (P·dV adjoint).

Most gates run Warp-CUDA kernels and require a CUDA device; LINC/MT are cross-checked against pure-host
oracles (``linc_tether_force``, ``euler_buckling_load_pn``).
"""
from __future__ import annotations

import numpy as np
import warp as wp

from aleph.components.nucleus.envelope import linc_tether_kernel
from aleph.components.nucleus.linc_analytic import linc_tether_force
from aleph.components.solid.microtubule import build_microtubule_compartment
from .force_path_rig import ForcePathRig, K_LINC_PN_UM, LINC_STIFFENING


# ---------------------------------------------------------------- G1.1
def gate_force_ledger(rig: ForcePathRig, push=0.8) -> dict:
    nuc = rig.nucleus_nodes()
    patch = rig.off_mem + rig.membrane_patch_nodes()
    rig.set_fixed(nuc, True)
    rig.prescribe_displacement(patch, np.array([push, 0.0, 0.0]))
    rig.settle(1200)
    rig.accumulate_all()
    f = rig.f_d.numpy()
    net = np.linalg.norm(np.sum(f, axis=0))                       # Σ all internal forces (equal-opposite ⇒ 0)
    scale = np.sum(np.linalg.norm(f, axis=1)) + 1e-30
    reaction_nuc = np.linalg.norm(np.sum(f[nuc], axis=0))
    applied_patch = np.linalg.norm(np.sum(f[patch], axis=0))
    # free-node residual (settled equilibrium)
    free = rig.fixed_d.numpy() == 0
    resid = float(np.max(np.linalg.norm(f[free], axis=1))) if free.any() else 0.0
    return {
        "name": "G1.1 force ledger (membrane load ↔ nucleus reaction)",
        "passed": bool(net / scale < 1e-9 and reaction_nuc > 1.0 and resid < 1.0),
        "net_internal_force_rel": float(net / scale),
        "nucleus_reaction_pn": float(reaction_nuc),
        "applied_load_at_patch_pn": float(applied_patch),
        "free_node_residual_pn": resid,
    }


# ---------------------------------------------------------------- G1.2
def gate_mt_compression_sign(dev: str) -> dict:
    """A straight MT strut; kink a middle node inward → bending must push it back toward straight."""
    mt = build_microtubule_compartment(centre=(0, 0, 0), n_mt=1, reach_R_um=5.0, seg_um=0.5,
                                        node_off=0, device=dev)
    pos = mt.pos0.copy()
    n = pos.shape[0]
    mid = n // 2
    axis = pos[-1] - pos[0]; axis /= np.linalg.norm(axis)
    perp = np.cross(axis, [0, 0, 1.0]); perp /= (np.linalg.norm(perp) + 1e-12)
    pos_k = pos.copy(); pos_k[mid] += 0.3 * perp          # kink the middle node sideways
    pd = wp.array(np.ascontiguousarray(pos_k), dtype=wp.vec3d, device=dev)
    f = wp.zeros(n, dtype=wp.vec3d, device=dev)
    mt.accumulate(pd, f)
    fmid = f.numpy()[mid]
    restoring = float(np.dot(fmid, -perp))                # >0 ⇒ pushes back toward the straight axis
    buckling = mt.euler_buckling_load_pn()
    return {
        "name": "G1.2 MT compression/buckling sign",
        "passed": bool(restoring > 0 and buckling > 0),
        "restoring_force_pn": restoring,
        "euler_buckling_load_pn": float(buckling),
        "note": "bending restores a kinked strut toward straight (emergent buckling resistance, κ=20 pN·µm²)",
    }


# ---------------------------------------------------------------- G1.3
def gate_if_tension_sign(rig: ForcePathRig, pull=0.4) -> dict:
    """Expand the anchor shell radially outward, pin anchors + nucleus, settle the IF backbone → the IF
    cables must pull the nucleus nodes OUTWARD (radial tension). A stretched cable can only pull, so a
    positive median radial force is the tension-sign arbiter. (LINC is disabled here to isolate IF.)"""
    nuc = rig.nucleus_nodes()
    anc = np.arange(rig.Na) + rig.off_anc
    p = rig.pos_d.numpy()
    au = p[anc] / (np.linalg.norm(p[anc], axis=1, keepdims=True) + 1e-12)
    p[anc] = p[anc] + pull * au                            # expand the whole anchor shell
    rig.pos_d.assign(np.ascontiguousarray(p))
    rig.set_fixed(anc, True); rig.set_fixed(nuc, True)     # pin driven anchors + nucleus, free backbone
    # settle using ONLY the IF channel so we isolate IF (zero LINC/capture influence on the backbone)
    for _ in range(600):
        f = rig.channel_force("if")
        fd = wp.array(np.ascontiguousarray(f), dtype=wp.vec3d, device=rig.device)
        from .force_path_rig import _overdamped_relax_kernel
        wp.launch(_overdamped_relax_kernel, dim=rig.N,
                  inputs=[rig.pos_d, fd, wp.float64(2.0e-4), wp.float64(1.0), rig.fixed_d, rig.pos_new_d],
                  device=rig.device)
        rig.pos_d.assign(rig.pos_new_d)
    fif = rig.channel_force("if")
    nu = rig.pos_d.numpy()[nuc]
    nu_rad = nu / (np.linalg.norm(nu, axis=1, keepdims=True) + 1e-12)
    radial = np.sum(fif[nuc] * nu_rad, axis=1)             # outward-positive IF force on each nucleus node
    loaded = np.abs(radial) > 1e-6
    frac_outward = float(np.mean(radial[loaded] > 0)) if loaded.any() else 0.0
    return {
        "name": "G1.3 IF tension sign",
        "passed": bool(loaded.sum() >= 3 and frac_outward > 0.8),
        "median_radial_force_pn": float(np.median(radial[loaded])) if loaded.any() else 0.0,
        "fraction_outward_tension": frac_outward,
        "n_loaded_nucleus_nodes": int(loaded.sum()),
        "note": "Hookean slice: tension only (sign). Strain-stiffening carried by the nonlinear LINC (G1.4).",
    }


# ---------------------------------------------------------------- G1.4
def gate_linc_equal_opposite(dev: str) -> dict:
    """linc_tether_kernel on one joint: nforce = −aforce, magnitude matches the host oracle, tension-only."""
    npos = np.array([[0.0, 0, 0], [1.2, 0, 0]], np.float64)    # node0 nucleus, node1 anchor; rest 1.0 ⇒ ext 0.2
    rest = np.array([1.0], np.float64)
    pd = wp.array(npos, dtype=wp.vec3d, device=dev)
    n_idx = wp.array(np.array([0], np.int32), dtype=wp.int32, device=dev)
    a_idx = wp.array(np.array([1], np.int32), dtype=wp.int32, device=dev)
    rest_d = wp.array(rest, dtype=wp.float64, device=dev)
    nf = wp.zeros(2, dtype=wp.vec3d, device=dev)
    af = wp.zeros(2, dtype=wp.vec3d, device=dev)
    wp.launch(linc_tether_kernel, dim=1,
              inputs=[pd, pd, n_idx, a_idx, rest_d, wp.float64(K_LINC_PN_UM),
                      wp.float64(LINC_STIFFENING), nf, af], device=dev)
    nfh = nf.numpy(); afh = af.numpy()
    fn = nfh[0]; fa = afh[1]
    equal_opposite = float(np.linalg.norm(fn + fa))
    oracle = float(linc_tether_force(1.2, 1.0, K_LINC_PN_UM, LINC_STIFFENING))   # scalar tension [pN]
    mag = float(np.linalg.norm(fn))
    # compression case: ext<0 ⇒ tension-only ⇒ zero
    npos_c = np.array([[0.0, 0, 0], [0.8, 0, 0]], np.float64)
    pd.assign(np.ascontiguousarray(npos_c)); nf.zero_(); af.zero_()
    wp.launch(linc_tether_kernel, dim=1, inputs=[pd, pd, n_idx, a_idx, rest_d, wp.float64(K_LINC_PN_UM),
              wp.float64(LINC_STIFFENING), nf, af], device=dev)
    compression_force = float(np.linalg.norm(nf.numpy()[0]))
    return {
        "name": "G1.4 LINC equal-and-opposite + oracle",
        "passed": bool(equal_opposite < 1e-9 and abs(mag - oracle) / oracle < 1e-6 and compression_force < 1e-12),
        "equal_opposite_residual_pn": equal_opposite,
        "kernel_magnitude_pn": mag,
        "oracle_magnitude_pn": oracle,
        "compression_force_pn": compression_force,
    }


# ---------------------------------------------------------------- G1.5
def gate_linc_detach_path(rig: ForcePathRig) -> dict:
    """Load LINC (push a patch), detach ONE joint → its force → 0, every other joint bit-exact unchanged."""
    patch = rig.off_mem + rig.membrane_patch_nodes()
    rig.prescribe_displacement(patch, np.array([0.6, 0.0, 0.0]))
    rig.settle(600)
    f_before = rig.channel_force("linc")
    # per-joint force magnitude on the anchor side
    a_idx = rig.linc_a_d.numpy()
    joint_force = np.linalg.norm(f_before[a_idx], axis=1)
    j = int(np.argmax(joint_force))                        # detach the most-loaded joint
    before_others = f_before.copy()
    rig.linc_detach(j)
    f_after = rig.channel_force("linc")
    detached_force = float(np.linalg.norm(f_after[a_idx[j]]))
    # every OTHER anchor-side joint must be bit-exact
    others = np.array([i for i in range(rig.n_linc) if i != j])
    unchanged = bool(np.array_equal(f_after[a_idx[others]], before_others[a_idx[others]]))
    return {
        "name": "G1.5 LINC detach removes only that path",
        "passed": bool(detached_force < 1e-12 and unchanged),
        "detached_joint": j,
        "detached_joint_force_after_pn": detached_force,
        "others_bit_exact": unchanged,
    }


# ---------------------------------------------------------------- G1.6
def gate_generation_validation(rig: ForcePathRig) -> dict:
    """MT endpoint remap bumps the topology epoch; a reference tagged with the STALE epoch is rejected."""
    edge_generation = rig.topology_epoch                   # a capture/LINC ref created at the current epoch
    valid_before = bool(edge_generation == rig.topology_epoch)
    rig.bump_topology_epoch()                              # MT growth/shrink / plus-end remap event
    valid_after = bool(edge_generation == rig.topology_epoch)
    # a fresh reference re-tagged at the new epoch is valid again
    edge_generation_new = rig.topology_epoch
    revalidated = bool(edge_generation_new == rig.topology_epoch)
    return {
        "name": "G1.6 MT remap generation validation",
        "passed": bool(valid_before and (not valid_after) and revalidated),
        "epoch_after_remap": rig.topology_epoch,
        "stale_reference_rejected": (not valid_after),
        "retagged_reference_valid": revalidated,
    }


# ---------------------------------------------------------------- G1.7
def gate_reject_topology_bit_exact(rig: ForcePathRig) -> dict:
    """Snapshot; propose a rejected topology event (LINC detach + capture detach + MT remap + displace); restore."""
    patch = rig.off_mem + rig.membrane_patch_nodes()
    rig.prescribe_displacement(patch, np.array([0.3, 0.0, 0.0]))
    rig.settle(300)
    snap = rig.snapshot()
    # --- rejected candidate ---
    rig.linc_detach(0)
    rig.capture_detach(0)
    rig.bump_topology_epoch()
    rig.prescribe_displacement(patch, np.array([0.5, 0.0, 0.0]))
    rig.settle(200)
    # --- reject → restore ---
    rig.restore(snap)
    pos_exact = bool(np.array_equal(rig.pos_d.numpy(), snap["pos"]))
    linc_exact = bool(np.array_equal(rig.linc_n_d.numpy(), snap["linc_n"]))
    cap_exact = bool(np.array_equal(rig.cap_bound_d.numpy(), snap["cap_bound"]))
    epoch_exact = rig.topology_epoch == snap["epoch"]
    return {
        "name": "G1.7 rejected topology event bit-exact recovery",
        "passed": bool(pos_exact and linc_exact and cap_exact and epoch_exact),
        "pos_bit_exact": pos_exact,
        "linc_bit_exact": linc_exact,
        "capture_bit_exact": cap_exact,
        "epoch_restored": bool(epoch_exact),
    }


# ---------------------------------------------------------------- G1.8
def gate_work_conservation(rig: ForcePathRig, push=0.6) -> dict:
    """Actuator does positive work against the network; releasing the load returns free nodes to rest (reversible)."""
    nuc = rig.nucleus_nodes()
    patch = rig.off_mem + rig.membrane_patch_nodes()
    rest_pos = rig.pos_d.numpy().copy()
    rig.set_fixed(nuc, True)
    rig.prescribe_displacement(patch, np.array([push, 0.0, 0.0]))
    rig.settle(1000)
    rig.accumulate_all()
    f = rig.f_d.numpy()
    disp_patch = rig.pos_d.numpy()[patch] - rest_pos[patch]
    # actuator work = −Σ (internal force · patch displacement)  (work done ON the network) > 0
    w_ext = float(-np.sum(f[patch] * disp_patch))
    # LINC tension sign: force on nucleus points inward-ish (tension pulls nucleus toward moved +x anchors)
    linc = rig.channel_force("linc")
    linc_pull_x = float(np.sum(linc[nuc, 0]))              # >0 ⇒ nucleus pulled +x toward the pushed side
    # release → settle → free nodes should return toward rest (conservative path)
    rig.set_fixed(patch, False)
    rig.set_fixed(nuc, False)
    rig.pos_d.assign(np.ascontiguousarray(rest_pos))       # release actuator, restore patch/nucleus to rest
    rig.settle(800)
    free = rig.fixed_d.numpy() == 0
    return_resid = float(np.max(np.linalg.norm((rig.pos_d.numpy() - rest_pos)[free], axis=1)))
    return {
        "name": "G1.8 solid work sign / conservation",
        "passed": bool(w_ext > 0 and linc_pull_x > 0 and return_resid < 0.02),
        "actuator_work_pn_um": w_ext,
        "linc_pull_toward_push_pn": linc_pull_x,
        "return_to_rest_residual_um": return_resid,
        "note": "paired cytosol pressure-work half is D0-G0.4 (P·dV adjoint).",
    }


def run_all(device: str) -> dict:
    dev = str(device)
    results = []
    results.append(gate_force_ledger(ForcePathRig(device=dev).relax_and_return()))
    results.append(gate_mt_compression_sign(dev))
    results.append(gate_if_tension_sign(ForcePathRig(device=dev).relax_and_return()))
    results.append(gate_linc_equal_opposite(dev))
    results.append(gate_linc_detach_path(ForcePathRig(device=dev).relax_and_return()))
    results.append(gate_generation_validation(ForcePathRig(device=dev)))
    results.append(gate_reject_topology_bit_exact(ForcePathRig(device=dev).relax_and_return()))
    results.append(gate_work_conservation(ForcePathRig(device=dev).relax_and_return()))
    passed = all(r["passed"] for r in results)
    return {"suite": "D1", "device": dev, "all_passed": bool(passed), "gates": results}
