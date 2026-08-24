"""B — node-by-node force parity: does the arena's ID-range addressing change the force?

THE CLAIM UNDER TEST, and it is the whole redesign in one sentence.  The incumbent gives every component
its own device arrays and indexes them from zero; the arena gives every population a contiguous ID range
of ONE allocation and indexes them globally.  If that substitution is sound, then accumulating the same
force law over the same geometry must give bit-identical forces — with every index shifted by the claim's
lower bound and nothing landing outside it.  If it is not sound, this is where it shows.

WHAT MAKES THIS A TEST OF THE LAYOUT RATHER THAN OF A RE-IMPLEMENTATION.  Both arms launch the SAME
production kernels — ``laws.network_warp.link_spring_kernel`` and ``laws.forces_warp.cytosim_bending_kernel``
— on the SAME positions and the SAME topology.  Nothing is re-derived here and no new force law is written.
The only difference between the arms is where the numbers live and what indices address them.  A
re-implementation would have tested my arithmetic; this tests the arrangement, which is the thing in doubt.

THE GEOMETRY IS THE INCUMBENT'S, not a builder of mine.  ``build_cell`` produces the full native cortex and
this reads its arrays; comparing forces on two different geometries would prove nothing, however closely
they resembled each other.

WHY THESE TWO CHANNELS AND NOT THE AXIAL BACKBONE.  ``cortex_state.CORTEX_CHANNELS_NOT_BOUND`` records
axial inextensibility as "an NF2007 constraint carried by the solver, not an accumulated force", with the
note that binding a spring there would double-count the backbone.  So the cortex's accumulated force is the
crosslink link-spring and the Cytosim bending triple, and those are exactly the two terms the stiffness row
sums are built from.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — positions [um], stiffness [pN/um], rest length [um], alpha [pN/um], force [pN].
  * boundary — the arena claim is deliberately placed at a NON-ZERO offset, because an offset of zero is
    the one case where a wrong addressing would still agree; a leading pad population is claimed first so
    lo > 0 is guaranteed rather than hoped for.
  * conservation/invariant — both arms accumulate atomically into a zeroed array, so ordering is
    non-deterministic within a launch and exact float equality is NOT assumed; the criterion is stated
    below and is a relative bound derived from that, not a chosen tolerance.
  * CFL/precision — nothing is integrated. float64 throughout, matching the production arrays.
  * sign sense — the residual is signed per component and its max ABSOLUTE value is reported, so a sign
    error cannot cancel against its mirror.
  * measurement protocol — one launch per channel per arm, then one host readback per arm after both
    launches complete. No readback inside any loop.

THE CRITERION, declared before the run.  Atomic accumulation makes summation order non-deterministic, so
the two arms may differ in the last bits even when the addressing is identical. The criterion is therefore
  max |f_arena - f_incumbent| / max |f_incumbent|  <=  n_terms * eps64
with n_terms the largest number of contributions landing on one node, taken from the topology rather than
assumed — the same Higham summation bound the accepted-step balance gate uses. PLUS an exact requirement
that carries no tolerance at all: every node outside the claimed range must be EXACTLY zero, because a
range that leaks is a defect no floating-point argument can excuse.

engine units: force pN, length um.  Runtime: NVIDIA Warp on CUDA only.
"""

from __future__ import annotations

import json
import sys

import numpy as np
import warp as wp


def main() -> int:
    wp.init()
    device = wp.get_device()
    if not device.is_cuda:
        raise RuntimeError(f"CUDA-only: resolved {device!r}")

    import aleph
    from aleph.components.incumbent.assemble import CellConfig, build_cell
    from aleph.laws.forces_warp import cytosim_bending_kernel
    from aleph.laws.network_warp import link_spring_kernel
    from aleph.world.arena import Kind, WorldArena

    cell = build_cell(CellConfig())
    n = int(cell.n_total)
    pos_h = cell.pos_d.numpy()
    xl_h = cell.xl_d.numpy().reshape(-1, 2).astype(np.int32)
    kxl_h = cell.kxl_d.numpy().reshape(-1).astype(np.float64)
    r0_h = cell.r0xl_d.numpy().reshape(-1).astype(np.float64)
    tri_h = cell.tri_d.numpy().reshape(-1, 3).astype(np.int32)
    alp_h = cell.alpha_d.numpy().reshape(-1).astype(np.float64)

    # Largest number of contributions landing on a single node, from the topology itself.
    hits = np.zeros(n, np.int64)
    np.add.at(hits, xl_h.reshape(-1).astype(np.int64), 1)
    np.add.at(hits, tri_h.reshape(-1).astype(np.int64), 1)
    n_terms = int(hits.max())
    tol = n_terms * float(np.finfo(np.float64).eps)

    def launch(pos_d, force_d, xl_d, k_d, r_d, tri_d, a_d) -> None:
        with wp.ScopedDevice(device):
            wp.launch(link_spring_kernel, dim=xl_h.shape[0], inputs=[pos_d, xl_d, k_d, r_d, force_d])
            wp.launch(cytosim_bending_kernel, dim=tri_h.shape[0], inputs=[pos_d, tri_d, a_d, force_d])
        wp.synchronize_device(device)

    # ── ARM A: the incumbent arrangement — a private array indexed from zero ────────────────────
    with wp.ScopedDevice(device):
        pos_a = wp.array(pos_h, dtype=wp.vec3d)
        f_a = wp.zeros(n, dtype=wp.vec3d)
        xl_a = wp.array(xl_h, dtype=wp.int32)
        k_a = wp.array(kxl_h, dtype=wp.float64)
        r_a = wp.array(r0_h, dtype=wp.float64)
        tri_a = wp.array(tri_h, dtype=wp.int32)
        alp_a = wp.array(alp_h, dtype=wp.float64)
    launch(pos_a, f_a, xl_a, k_a, r_a, tri_a, alp_a)
    f_inc = f_a.numpy()

    # ── ARM B: the arena — a claimed range of ONE allocation, indexed globally ──────────────────
    pad = 1_237                                     # a deliberately non-round leading claim
    arena = WorldArena(capacity={Kind.NODE: n + pad + 4_096}, device=str(device))
    arena.claim("pad", Kind.NODE, pad)               # forces lo > 0: offset 0 is the case that always agrees
    cortex = arena.claim("cortex", Kind.NODE, n)
    arena.assert_partitioned()
    arena.upload_nodes(cortex, "position", pos_h)

    with wp.ScopedDevice(device):
        xl_b = wp.array(xl_h.astype(np.int64) + cortex.lo, dtype=wp.int32)
        tri_b = wp.array(tri_h.astype(np.int64) + cortex.lo, dtype=wp.int32)
    launch(arena.node_arrays["position"], arena.node_arrays["force"], xl_b, k_a, r_a, tri_b, alp_a)
    f_world_all = arena.node_arrays["force"].numpy()
    f_world = f_world_all[cortex.lo:cortex.hi]

    scale = float(np.abs(f_inc).max())
    resid = float(np.abs(f_world - f_inc).max())
    rel = resid / scale if scale else 0.0
    outside = np.concatenate([f_world_all[:cortex.lo], f_world_all[cortex.hi:arena.n_live(Kind.NODE)]])
    leak = float(np.abs(outside).max()) if outside.size else 0.0

    # Per-node residual against contribution count: if the residual is summation-order noise it must
    # RISE with the number of terms landing on a node, and if it is an addressing defect it will not.
    per_node = np.abs(f_world - f_inc).max(axis=1)
    by_hits = {}
    for h in range(1, n_terms + 1):
        m = hits == h
        if m.sum():
            by_hits[h] = {"n_nodes": int(m.sum()), "max_resid_pN": float(per_node[m].max()),
                          "mean_resid_pN": float(per_node[m].mean())}
    np.savez_compressed("/tmp/b_parity_pernode.npz", per_node=per_node, hits=hits,
                        f_inc_mag=np.linalg.norm(f_inc, axis=1))

    lo_b, hi_b = arena.byte_span("position", cortex)
    base = int(arena.node_arrays["position"].ptr)

    record = {
        "schema": "diagnostic@1",
        "run_label": "b_force_parity_arena_vs_private",
        "kind": "diagnostic",
        "quantitative_claim_status": "BLOCKED",
        "evidence_basis": ("the SAME production kernels on the SAME geometry, launched over a private "
                           "zero-based array and over a claimed range of one arena; no force law is "
                           "re-implemented here"),
        "stamp": {"aleph_file": aleph.__file__, "warp_version": wp.__version__,
                  "python": sys.version.split()[0], "device": str(device)},
        "census": {"n_nodes": n, "n_crosslinks": int(xl_h.shape[0]), "n_bend_triples": int(tri_h.shape[0]),
                   "max_contributions_per_node": n_terms},
        "channels": ["actin_crosslink_link_spring", "actin_bending_cytosim"],
        "channels_excluded": {
            "axial_inextensibility": "an NF2007 solver constraint, not an accumulated force "
                                     "(cortex_state.CORTEX_CHANNELS_NOT_BOUND)",
        },
        "arena": {"pad_claim": pad, "cortex_lo": cortex.lo, "cortex_hi": cortex.hi,
                  "byte_span": [lo_b, hi_b], "array_base_ptr": base,
                  "claim_offset_bytes": lo_b - base},
        "criterion": {
            "declared_before_run": True,
            "relative_bound": tol,
            "expression": "max|f_arena - f_incumbent| / max|f_incumbent| <= n_terms * eps64",
            "why_not_exact": ("both arms accumulate atomically, so summation order is non-deterministic "
                              "within a launch; the bound is the Higham summation floor, not a chosen "
                              "tolerance"),
            "exact_requirement": "every node outside the claimed range must be EXACTLY zero",
        },
        "measured": {
            "max_abs_force_pN": scale,
            "max_abs_residual_pN": resid,
            "relative_residual": rel,
            "within_bound": bool(rel <= tol),
            "leak_outside_range_pN": leak,
            "no_leak": bool(leak == 0.0),
        },
        "residual_by_contribution_count": by_hits,
        "verdict": "PARITY" if (rel <= tol and leak == 0.0) else "MISMATCH",
        "timing": {"comparable": False,
                   "not_comparable_reason": "two force assemblies; no solve and no accepted step"},
    }
    print(json.dumps(record, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
