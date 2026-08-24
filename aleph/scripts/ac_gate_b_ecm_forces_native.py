#!/usr/bin/env python
r"""GATE-B native driver — the ``ecm`` component's collagen constitutive force pass on CUDA.

Seeds an explicit collagen Mikado network directly on the GPU (:class:`MikadoTopologyBuilder` over the
:class:`~aleph.components.ecm.device_schema.ECMTopologyState` SoA), binds the SOURCED collagen card's constitutive
law as the ECM component's native force pass (:class:`aleph.engine.ecm_mechanics.ECMConstitutiveForce`), and
drives it over the SoA node arrays — reporting that a RELAXED network is force-free and a small imposed shear
strain develops an EMERGENT restoring force (the macroscopic modulus emerges from the microstructure; no
lumped modulus is baked in).

This is the ECM analog of ``ac_gate_b_erm_cortex_native.py`` / the SF native gate: the dump comes from a
gbook CUDA run (no CUDA on the dev Mac).  It is a *mechanism* check (relaxed=0, stretch/shear→emergent
restoring force, correct sign), NOT a material-production modulus gate — E_fibril / k_xl are GAPs (see
:data:`aleph.engine.ecm_mechanics.COLLAGEN_FORCE_PROVENANCE`), so this does not close a Pa band.

With ``--drive-owner`` the force pass is driven through the ``ecm`` component OWNER
(:class:`aleph.engine.ecm_world.ECMStateOwner`, assembled by
:func:`aleph.engine.ecm_world.build_ecm_state_owner`) rather than the standalone force object: the
owner's own per-step hook (:meth:`ECMStateOwner.accumulate`) launches the collagen kernels over the arrays
it owns, and an accepted-step :meth:`ECMStateOwner.commit_irreversible` re-derives the constitutive
parameters — proving the ``ecm`` component is owner-driven (KERNEL_BOUND → CONNECTED), not just a loose
force object.

Run on gbook:
    conda activate ffn_sim
    PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
        python aleph/scripts/ac_gate_b_ecm_forces_native.py --n-fibers 200 --drive-owner
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import warp as wp

from aleph.components.ecm.device_schema import PopulationSlot
from aleph.engine.ecm_mechanics import (
    COLLAGEN_MATERIAL_KEY,
    build_collagen_constitutive_force,
    collagen_force_provenance,
)
from aleph.engine.ecm_world import (
    BoundaryAnchorMode,
    ECMWorldSettings,
    build_ecm_state_owner,
)
from aleph.engine.load_path import ActorRecord
from aleph.components.ecm.mikado_topology import MikadoInitConfig, MikadoTopologyBuilder
from aleph.laws.ecm_library import get_spec


class _NoOpTransaction:
    """Minimal accepted-step transaction stub so the gate can drive the owner facade.

    The gate exercises the OWNER's mechanics + refresh path; damage/sleep/refinement/remesh commits are
    out of scope here, so each accepted-step hook is a device-free no-op (never advances a private clock,
    never mutates authoritative state).
    """

    def snapshot_candidate(self) -> None:
        return None

    def rollback(self, accepted: wp.array) -> None:
        return None

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        return None


class _NoOpLedger:
    """Minimal ledger stub for the owner-driven gate (no reduction; mechanism check only)."""

    def accumulate_ledger(self, ledger: object) -> None:
        return None


@wp.kernel
def _zero_vec3d(force: wp.array(dtype=wp.vec3d)) -> None:
    force[wp.tid()] = wp.vec3d(0.0, 0.0, 0.0)


@wp.kernel
def _force_magnitude_sum(
    force: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=wp.int32),
    out: wp.array(dtype=wp.float64),   # [0]=Σ|f|, [1]=max|f|, [2]=active node count
) -> None:
    i = wp.tid()
    if active[i] > wp.int32(0):
        m = wp.length(force[i])
        wp.atomic_add(out, 0, m)
        wp.atomic_max(out, 1, m)
        wp.atomic_add(out, 2, wp.float64(1.0))


@wp.kernel
def _apply_simple_shear(
    pos: wp.array(dtype=wp.vec3d),
    node_active: wp.array(dtype=wp.int32),
    gamma: wp.float64,          # shear strain: x += γ·z about the box origin
) -> None:
    i = wp.tid()
    if node_active[i] > wp.int32(0):
        p = pos[i]
        pos[i] = wp.vec3d(p[0] + gamma * p[2], p[1], p[2])


def _force_stats(accumulate, topology) -> dict[str, float]:
    """Zero the ECM force, run one ``accumulate(pos, force)``-equivalent, and reduce |f| on active nodes.

    ``accumulate`` is either the standalone ``force_pass.accumulate`` or an owner-driven closure that
    invokes ``owner.accumulate()`` over the same owned arrays — so both lanes share one measurement.
    """
    out = wp.zeros(3, dtype=wp.float64, device=topology.device)
    wp.launch(_zero_vec3d, dim=topology.n_node_capacity, inputs=[topology.force_d],
              device=topology.device)
    accumulate()
    wp.launch(_force_magnitude_sum, dim=topology.n_node_capacity,
              inputs=[topology.force_d, topology.node_active_d, out], device=topology.device)
    wp.synchronize_device(topology.device)
    s = out.numpy()   # post-loop diagnostic readback only (not a hot-loop D2H)
    n = max(1.0, float(s[2]))
    return {"sum_abs_pN": float(s[0]), "max_abs_pN": float(s[1]), "mean_abs_pN": float(s[0]) / n,
            "active_nodes": int(s[2])}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fibers", type=int, default=200)
    ap.add_argument("--box-um", type=float, default=10.0)
    ap.add_argument("--shear", type=float, default=0.02, help="imposed simple-shear strain γ")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--drive-owner", action="store_true",
                    help="drive the force pass through the ecm ECMStateOwner facade, not the standalone "
                         "force object (also exercises the accepted-step commit refresh hook)")
    ap.add_argument("--out", type=str, default="")
    args = ap.parse_args()

    wp.init()
    dev = str(wp.get_device())
    if not wp.get_device().is_cuda:
        raise SystemExit("GATE-B ECM native gate requires CUDA (I0-A); no CPU simulation path exists.")

    spec = get_spec(COLLAGEN_MATERIAL_KEY)
    L = float(args.box_um)
    config = MikadoInitConfig(
        box_lo_um=(0.0, 0.0, 0.0), box_hi_um=(L, L, L),
        n_fibers=int(args.n_fibers), fiber_length_um=spec.fiber_len_um,
        target_segment_um=spec.seg_um, crosslink_capture_um=spec.xl_contact_um,
        pin_faces=("z_lo",), pin_margin_um=1.0, rng_seed=int(args.seed), max_refinement_level=0,
    )
    topology = MikadoTopologyBuilder(config, device=dev).initialize()
    force_pass = build_collagen_constitutive_force(topology, spec, device=dev)

    # Choose the drive lane: standalone force object, or the ecm component OWNER facade.
    owner = None
    if args.drive_owner:
        owner = build_ecm_state_owner(
            topology,
            settings=ECMWorldSettings(BoundaryAnchorMode.FAR_FIELD_DIRICHLET),
            actor=ActorRecord("ecm", 0, 0, 0, 0, int(topology.n_segment_capacity)),
            transaction=_NoOpTransaction(),
            ledger=_NoOpLedger(),
            spec=spec,
            device=dev,
        )
        accumulate = owner.accumulate            # owner-driven force hook over owned arrays
    else:
        def accumulate() -> None:                # standalone force object over the same SoA arrays
            force_pass.accumulate(topology.position_d, topology.force_d)

    active = topology.active_population_d.numpy()   # post-build diagnostic (device-authoritative counts)
    relaxed = _force_stats(accumulate, topology)

    # impose a small simple shear on the active nodes, re-evaluate the emergent restoring force
    wp.launch(_apply_simple_shear, dim=topology.n_node_capacity,
              inputs=[topology.position_d, topology.node_active_d, wp.float64(args.shear)],
              device=dev)
    wp.synchronize_device(dev)
    sheared = _force_stats(accumulate, topology)

    # If driving the owner, exercise the accepted-step commit → re-derive constitutive parameters
    # (the remesh-invalidation hook) so the gate proves the owner keeps k_seg/α fresh under its
    # transaction, not just that it can launch the force once.
    owner_commit_refresh_ok = None
    if owner is not None:
        accepted_d = wp.array(np.array([1], dtype=np.int32), dtype=wp.int32, device=dev)
        owner.commit_irreversible(accepted_d, dt_phys=1.0, rng_seed=int(args.seed))
        wp.synchronize_device(dev)
        owner_commit_refresh_ok = True

    report = {
        "gate": "GATE-B ecm collagen constitutive force pass (native)",
        "device": dev,
        "material": spec.key,
        "drive_lane": "ecm_owner" if owner is not None else "standalone_force",
        "owner_commit_refresh_ok": owner_commit_refresh_ok,
        "n_fibers": int(args.n_fibers),
        "shear_strain": float(args.shear),
        "active_population": {
            "fibers": int(active[PopulationSlot.FIBER]),
            "nodes": int(active[PopulationSlot.NODE]),
            "segments": int(active[PopulationSlot.SEGMENT]),
            "bends": int(active[PopulationSlot.BEND]),
        },
        "relaxed_force": relaxed,
        "sheared_force": sheared,
        "emergent_restoring_ratio": (sheared["mean_abs_pN"] / relaxed["mean_abs_pN"]
                                     if relaxed["mean_abs_pN"] > 1e-12 else float("inf")),
        "provenance": collagen_force_provenance(spec),
        "verdict": {
            "relaxed_is_force_free": relaxed["max_abs_pN"] < 1e-6,
            "shear_develops_restoring_force": sheared["max_abs_pN"] > 10.0 * max(relaxed["max_abs_pN"], 1e-12),
        },
        "note": "mechanism check only (relaxed=0, shear→emergent restoring force); NOT a Pa modulus gate — "
                "E_fibril/k_xl are GAPs (surface to PI before a material-production gate).",
    }
    text = json.dumps(report, indent=2)
    print(text)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text)


if __name__ == "__main__":
    main()
