r"""PHASE 2 — per-law force parity between the arena and the incumbent, and each law's marginal cost.

WHAT THIS EXTENDS.  ``arena_force_parity.py`` (``78942fc4``) asked one question about the CORTEX's two
channels: does addressing a population as a claimed range of one allocation change the force?  It
answered 1.108e-16 relative, with exactly 0.0 pN outside the claim.  This script asks the same question
once per LAW, over the SURFACE populations, and adds the measurement PHASE 2 exists to take: what each
law costs on its own.

WHY THE SURFACE LAWS AND NOT THE CORTEX'S.  The cortex's two are already answered by the run above.  The
surface laws are where the arena's geometry layer differs from the incumbent's in a way that could
plausibly matter and had not been checked: ``world.surface`` builds its own hinge list, in its own row
order, with its own column convention, whereas the incumbent carries the list
``laws.membrane_surface.build_membrane_hinges`` produced.  Two different hinge lists over the same mesh
must give the same assembled force or the arena's surface builder is not the incumbent's surface.

THREE ARMS, AND WHY THE THIRD EXISTS.
  * **A — the incumbent.**  ``build_cell`` composes the membrane and its ``MembraneCompartment`` already
    launches ``helfrich_bending_kernel`` and ``membrane_area_kernel`` over GLOBAL indices into the
    combined node array.  Nothing is re-implemented for this arm; it is the production path.
  * **B — the arena.**  The same vertices, claimed as a range behind a deliberately non-round pad, with
    ``world.surface`` building the faces and hinges and ``world.laws_bind`` launching the same kernels.
  * **C — the host oracle, for the volume law only.**  ``laws.volume`` folds a gradient that had three
    copies, and the ONE device copy (``components/incumbent/membrane_pressure.py:48``) is fused to a
    Biot pressure trace and cannot be launched without a fluid grid.  So that law's A-arm is the two
    HOST copies it folds — ``uniform_pressure_force_reference`` (the incumbent's own stated analytic
    oracle, written in the consistent-load form) and ``nucleus/geometry.py``'s ``mesh_volume_gradient``
    (written in the own-term form).  Agreeing with BOTH is the claim; agreeing with one would only show
    the fold reproduced whichever form it copied.

WHAT MAY NOT BE CLAIMED FROM THIS RUN.  That any residual is lower, that any solve converges, or that a
step is accepted — none of that is touched here and all of it is PHASE 4.  Parity is a statement about
ADDRESSING, and the cost figures are per-kernel launch costs, not step costs.

THE CRITERION, DECLARED BEFORE THE RUN.  Every arm accumulates atomically, so summation order is not
reproducible even within one launch, and the two arms visit hinges in different orders by construction.
The bound is therefore the Higham summation floor rather than a chosen tolerance:

    max |f_B - f_A| / max |f_A|  <=  n_terms * eps64

with ``n_terms`` the largest number of contributions landing on any single node, counted from the
topology itself.  PLUS an exact requirement that carries no tolerance: every node outside the claimed
range must be EXACTLY zero.  A range that leaks is a defect no floating-point argument excuses.

THE COST PROTOCOL, DECLARED BEFORE THE RUN.  Per law: 3 untimed warm-up launches (Warp compiles on first
use and a first-launch time is a compile time), then 20 timed launches each bracketed by
``wp.synchronize_device``, reported as the MEDIAN with the min and max beside it.  The projection to the
plan's subdivision-7 membrane is linear in the launch dimension and is labelled an extrapolation on the
line that carries it, because these kernels are bandwidth-bound and a projection is not a measurement.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — positions [um]; sigma [pN/um]; kappa_tilde [pN.um]; delta_p [pN/um^2]; forces [pN].
  * boundary — the arena claim sits behind a non-round PAD so its lower bound is strictly positive;
    offset zero is the one case where wrong addressing still agrees.  ``delta_p = 0`` is checked to give
    bit-zero force before any parity number is read.
  * conservation/invariant — the arena surface is built on the incumbent's OWN vertex positions scaled
    by exactly 1.0, so the two arms are the same geometry rather than two similar ones; both arms'
    enclosed volumes are reported so that assumption is visible rather than trusted.
  * CFL/precision — nothing is integrated.  float64 throughout.
  * sign sense — residuals are signed per component and reported as a max ABSOLUTE, so a sign error
    cannot cancel against its mirror; the volume law additionally reports the virial identity.
  * measurement protocol — one launch per channel per arm, then ONE host readback per arm after all
    launches complete.  No readback inside any loop.  Timing launches read nothing back.

engine units: force pN, length um.  Runtime: NVIDIA Warp on CUDA only.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import time

import numpy as np
import warp as wp

#: The plan's subdivision-7 plasma membrane, for the labelled linear projection only.
SUBDIV7 = {"vertices": 163_842, "faces": 327_680, "hinges": 491_520}

WARMUP_LAUNCHES = 3
TIMED_LAUNCHES = 20


def _contributions_per_node(n: int, *index_arrays) -> int:
    """Largest number of atomic contributions landing on one node, counted from the topology."""
    hits = np.zeros(n, np.int64)
    for idx in index_arrays:
        if idx is None or idx.size == 0:
            continue
        np.add.at(hits, np.asarray(idx, np.int64).reshape(-1), 1)
    return int(hits.max()) if hits.size else 0


def _time_launch(fn, device) -> dict[str, float]:
    """Median/min/max wall time of ``fn`` [ms], after warm-up, each launch device-synchronised."""
    for _ in range(WARMUP_LAUNCHES):
        fn()
    wp.synchronize_device(device)
    samples = []
    for _ in range(TIMED_LAUNCHES):
        t0 = time.perf_counter()
        fn()
        wp.synchronize_device(device)
        samples.append((time.perf_counter() - t0) * 1.0e3)
    return {"median_ms": statistics.median(samples), "min_ms": min(samples), "max_ms": max(samples)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="", help="write the record here as well as to stdout")
    args = ap.parse_args()

    wp.init()
    device = wp.get_device()
    if not device.is_cuda:
        raise RuntimeError(f"CUDA-only: resolved {device!r}")

    import aleph
    from aleph.components.incumbent.assemble import CellConfig, build_cell
    from aleph.components.incumbent.membrane_pressure import uniform_pressure_force_reference
    from aleph.components.nucleus.geometry import mesh_volume_gradient
    from aleph.laws.membrane_surface import helfrich_bending_kernel, membrane_area_kernel
    from aleph.laws.volume import mesh_volume, volume_pressure_kernel
    from aleph.world.arena import Kind, WorldArena
    from aleph.world.laws_bind import (accumulate_area_tension, accumulate_helfrich,
                                       accumulate_volume_pressure, upload_surface_topology)
    from aleph.world.surface import build_surface

    cell = build_cell(CellConfig())
    mem = cell.membrane
    if mem is None:
        raise RuntimeError(
            "the assembled cell carries no membrane compartment, so there is no surface to compare. "
            "This script tests the SURFACE laws; the cortex's two are answered by arena_force_parity.py."
        )

    # The incumbent's own live vertices and topology — read once, on the host, before anything launches.
    pos_all = cell.pos_d.numpy()
    lo_inc, hi_inc = int(mem.node_off), int(mem.node_off) + int(mem.n_verts)
    verts = np.ascontiguousarray(pos_all[lo_inc:hi_inc], np.float64)
    faces_local = np.ascontiguousarray(mem.mesh.faces, np.int64)
    hinges_inc = mem.hinges_d.numpy().reshape(-1, 4).astype(np.int64) - lo_inc

    n_terms = _contributions_per_node(int(mem.n_verts), faces_local, hinges_inc)
    tol = n_terms * float(np.finfo(np.float64).eps)

    sigma = float(mem.gamma_mem)          # the compartment's OWN resolved tension; nothing chosen here
    kappa_t = float(mem.kappa_tilde)      # its OWN calibrated dihedral coupling
    delta_p = 40.0                        # a pressure JUMP for the A/B; not a sourced physiological value

    # ── ARM A: the incumbent arrangement — the production kernels over the combined cell array ──────
    def _arm_a(kernel, dim, topo_d, coeff):
        with wp.ScopedDevice(device):
            f = wp.zeros(cell.n_total, dtype=wp.vec3d)
        wp.launch(kernel, dim=dim, inputs=[cell.pos_d, topo_d, wp.float64(coeff), f], device=str(device))
        wp.synchronize_device(device)
        return f.numpy()[lo_inc:hi_inc]

    f_a_area = _arm_a(membrane_area_kernel, mem.n_faces, mem.faces_d, sigma)
    f_a_bend = _arm_a(helfrich_bending_kernel, mem.n_hinges, mem.hinges_d, kappa_t)

    # ── ARM B: the arena — the same vertices, behind a non-round pad, built by world.surface ────────
    pad = 1_237
    arena = WorldArena(
        capacity={Kind.NODE: int(mem.n_verts) + pad + 4_096,
                  Kind.FACE: int(mem.n_faces) + 64,
                  Kind.ANGLE4: int(mem.n_hinges) + 64},
        device=str(device),
    )
    arena.claim("pad", Kind.NODE, pad)
    # radius_um=1.0 and centre=0 make build_surface's transform exactly the identity, so this is the
    # incumbent's geometry rather than a mesh that resembles it.
    surf = build_surface(arena, "membrane", vertices=verts, faces=faces_local, radius_um=1.0)
    arena.assert_partitioned()
    arena.upload_nodes(surf.nodes, "position", surf.position)
    topo = upload_surface_topology(arena, surf)

    def _arm_b(accumulate, coeff):
        arena.node_arrays["force"].zero_()
        accumulate(arena, topo, coeff)
        wp.synchronize_device(device)
        full = arena.node_arrays["force"].numpy()
        outside = np.concatenate([full[:surf.nodes.lo], full[surf.nodes.hi:arena.n_live(Kind.NODE)]])
        return full[surf.nodes.lo:surf.nodes.hi], float(np.abs(outside).max()) if outside.size else 0.0

    f_b_area, leak_area = _arm_b(accumulate_area_tension, sigma)
    f_b_bend, leak_bend = _arm_b(accumulate_helfrich, kappa_t)
    f_b_vol, leak_vol = _arm_b(accumulate_volume_pressure, delta_p)

    # delta_p = 0 must be bit-zero before any pressure number is read.
    f_b_zero, _ = _arm_b(accumulate_volume_pressure, 0.0)

    # ── ARM C: the two host copies the volume law folds ─────────────────────────────────────────────
    f_c_consistent = uniform_pressure_force_reference(verts, faces_local, delta_p)
    f_c_ownterm = delta_p * mesh_volume_gradient(verts, faces_local)

    def _compare(f_b, f_a, leak) -> dict[str, object]:
        scale = float(np.abs(f_a).max())
        resid = float(np.abs(f_b - f_a).max())
        rel = resid / scale if scale else 0.0
        return {"max_abs_force_pN": scale, "max_abs_residual_pN": resid, "relative_residual": rel,
                "within_bound": bool(rel <= tol), "leak_outside_range_pN": leak,
                "no_leak": bool(leak == 0.0)}

    laws = {
        "membrane_area_tension": {
            "kernel": "laws.membrane_surface.membrane_area_kernel", "binding": "BIND",
            "parameter": {"sigma_pN_per_um": sigma, "source": "MembraneCompartment.gamma_mem"},
            "launch_dim": int(mem.n_faces), "dim_kind": "FACE",
            "measured": _compare(f_b_area, f_a_area, leak_area),
        },
        "helfrich_bending": {
            "kernel": "laws.membrane_surface.helfrich_bending_kernel", "binding": "BIND",
            "parameter": {"kappa_tilde_pN_um": kappa_t, "source": "MembraneCompartment.kappa_tilde"},
            "launch_dim": int(mem.n_hinges), "dim_kind": "ANGLE4",
            "note": ("arm A carries build_membrane_hinges' row order and arm B world.surface._hinges' "
                     "own, permuted (1,2,0,3); the i<->j swap leaves n1.n2 and its gradient unchanged"),
            "measured": _compare(f_b_bend, f_a_bend, leak_bend),
        },
        "volume_pressure": {
            "kernel": "laws.volume.volume_pressure_kernel", "binding": "FOLD",
            "folded_from": ["components/nucleus/geometry.py:115 (host, own-term)",
                            "components/incumbent/membrane_pressure.py:48 (device, consistent-load, "
                            "fused to the Biot pressure trace)",
                            "components/incumbent/resting_balance_oracle.py:94 (host, own-term)"],
            "parameter": {"delta_p_pN_per_um2": delta_p, "source": "A/B PROBE VALUE, not physiological"},
            "launch_dim": int(mem.n_faces), "dim_kind": "FACE",
            "vs_consistent_load_form": _compare(f_b_vol, f_c_consistent, leak_vol),
            "vs_own_term_form": _compare(f_b_vol, f_c_ownterm, leak_vol),
            "delta_p_zero_is_bit_zero": bool(np.array_equal(f_b_zero, np.zeros_like(f_b_zero))),
            "net_force_pN": float(np.abs(f_b_vol.sum(axis=0)).max()),
            "virial_identity": {
                "sum_f_dot_x": float(np.einsum("ij,ij->", f_b_vol, verts)),
                "three_delta_p_V": 3.0 * delta_p * mesh_volume(verts, faces_local),
            },
        },
    }

    # ── the measurement PHASE 2 exists to take: each law's marginal cost, alone ─────────────────────
    cost_fns = {
        "membrane_area_tension": lambda: accumulate_area_tension(arena, topo, sigma),
        "helfrich_bending": lambda: accumulate_helfrich(arena, topo, kappa_t),
        "volume_pressure": lambda: accumulate_volume_pressure(arena, topo, delta_p),
    }
    for name, fn in cost_fns.items():
        t = _time_launch(fn, device)
        dim = laws[name]["launch_dim"]
        kind = laws[name]["dim_kind"]
        projected_dim = SUBDIV7["faces"] if kind == "FACE" else SUBDIV7["hinges"]
        laws[name]["cost"] = {
            **t,
            "ns_per_element": t["median_ms"] * 1.0e6 / dim if dim else None,
            "projection_subdiv7_ms": t["median_ms"] * projected_dim / dim if dim else None,
            "projection_is_an_extrapolation": True,
            "projection_basis": (f"linear in the launch dimension ({dim} -> {projected_dim} {kind}); "
                                 "these kernels are bandwidth-bound, so this is a projection and not a "
                                 "measurement of the subdivision-7 membrane"),
        }

    parity_ok = all(
        (law["measured"]["within_bound"] and law["measured"]["no_leak"]) if "measured" in law else
        (law["vs_consistent_load_form"]["within_bound"] and law["vs_own_term_form"]["within_bound"]
         and law["vs_consistent_load_form"]["no_leak"] and law["delta_p_zero_is_bit_zero"])
        for law in laws.values()
    )

    record = {
        "schema": "diagnostic@1",
        "run_label": "arena_law_parity_surface",
        "kind": "diagnostic",
        "quantitative_claim_status": "BLOCKED",
        "phase": "WORLD_PORT_PLAN_2026-08-20 PHASE 2",
        "evidence_basis": (
            "the SAME production kernels on the SAME vertices, launched over the incumbent's combined "
            "array and over a claimed arena range behind a non-round pad; no force law is "
            "re-implemented, and the one law that is new (laws.volume) is compared against BOTH host "
            "copies it folds"),
        "may_not_be_claimed": [
            "that any residual is lower — nothing is solved here and that is PHASE 4",
            "that these are step costs — they are per-kernel launch costs with no solve",
            "the subdivision-7 projections as measurements",
        ],
        "stamp": {"aleph_file": aleph.__file__, "warp_version": wp.__version__,
                  "python": sys.version.split()[0], "device": str(device)},
        "census": {
            "membrane_vertices": int(mem.n_verts), "faces": int(mem.n_faces),
            "hinges": int(mem.n_hinges), "cell_n_total": int(cell.n_total),
            "max_contributions_per_node": n_terms,
            "enclosed_volume_um3_arena": surf.volume0_um3,
            "enclosed_volume_um3_incumbent": mesh_volume(verts, faces_local),
        },
        "arena": {"pad_claim": pad, "membrane_lo": surf.nodes.lo, "membrane_hi": surf.nodes.hi,
                  "byte_span": list(arena.byte_span("position", surf.nodes)),
                  "array_base_ptr": int(arena.node_arrays["position"].ptr)},
        "criterion": {
            "declared_before_run": True,
            "relative_bound": tol,
            "expression": "max|f_arena - f_incumbent| / max|f_incumbent| <= n_terms * eps64",
            "why_not_exact": ("every arm accumulates atomically and the two arms visit hinges in "
                              "different orders by construction; the bound is the Higham summation "
                              "floor, not a chosen tolerance"),
            "exact_requirement": "every node outside the claimed range must be EXACTLY zero",
            "cost_protocol": (f"{WARMUP_LAUNCHES} untimed warm-up launches then {TIMED_LAUNCHES} timed, "
                              "each device-synchronised; median reported with min and max"),
        },
        "laws": laws,
        "verdict": "PARITY" if parity_ok else "MISMATCH",
        "timing": {"comparable": False,
                   "not_comparable_reason": ("per-kernel launch costs on one surface population; no "
                                             "solve, no accepted step, and no other law present")},
    }
    text = json.dumps(record, indent=1)
    print(text)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text + "\n")
    return 0 if parity_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
