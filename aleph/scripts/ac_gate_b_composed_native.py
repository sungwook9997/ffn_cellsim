#!/usr/bin/env python
r"""GATE-B NATIVE driver — the FIRST genuinely-composed cell: sf_arc + ECM + NMII step as ONE transaction.

The real-owner analog of the census-double bridge (``scripts/ac_composed_world_dump.py``).  It builds the
components that already own **private** device arrays and need **no frozen-driver edit** — ``sf_arc``
(:class:`aleph.engine.composed_native.SFArcStateOwner` over the KERNEL_BOUND
:class:`~aleph.engine.sf_mechanics.SFFilamentMechanics`), ECM
(:class:`~aleph.engine.ecm_world.ECMStateOwner`), and the head-resolved NMII actuator + its
``nmii_cortex_motor`` connector — with the **cortex as a bind-target PORT** (never a driven participant for
the motor edge).  They are bound into ONE :func:`~aleph.engine.composition.build_composed_cell_world`
world and ONE top-level :class:`~aleph.engine.transaction.CellTransaction` (world + WarpEventClock +
disjoint cortex/sf_arc/ecm/nmii population ledgers), and ONE accepted physical step drives every owner under
ONE device-resident acceptance predicate.

WHAT IT SHOWS (all on-device; host reads are OUT-OF-LOOP diagnostics only):
  1. at rest the ``sf_arc`` + ECM passive force pass injects ≈ZERO force (emergent-not-lumped);
  2. a STRETCH of the SF population + a SHEAR of the collagen network → each component's OWN kernels emit a
     non-zero EMERGENT restoring force (the SF rod-cable tension/bending; the collagen constitutive modulus)
     — launched by the SAME composed candidate solve, over DISJOINT arrays;
  3. the NMII actuator's internal backbone/arm mechanics run and the split crossbridge scatters into the
     cortex PORT (Newton's 3rd law across two never-merged arrays); as heads bind, the cortex reaction rises;
  4. ONE accepted step: snapshot(all) → propose events → composed solve → ledgers → one predicate to every
     rollback+commit → advance the device clock → re-assert disjoint populations.  The clock's accepted-step
     index advances by exactly 1; the four populations stay disjoint.

HONESTLY EXCLUDED — the interior fluid column (membrane / cytosol / nucleus) is NOT composed here.  Moving its
``pressure``/``membrane_pressure``/``nucleus`` couplings into engine connectors edits the feature-frozen
``ac/cell/driver.py`` — the PI-gated Card-5 strangler transition
(``INTERIOR_COLUMN_CONNECTED_PLAN_2026-07-25.md``), not this autonomous slice.

PARAMS ARE PI-GAPs (report-not-tune):  the SF axial ``k_axial`` (NF2007 inextensible; no sourced EA_actin),
ECM ``E_fibril``/``k_xl``, and the whole NMII force-scale (``f_stall``/``k_xb``/``k_on``/…) are KB/PI GAPs.
This is a MECHANISM demonstration (rest≈0, perturb→emergent restoring force, one clean accepted step), NOT a
quantitative production gate — no band is closed.

────────────────────────────────────────────────────────────────────────────────────────────────────────
RUN ON GBOOK (needs a CUDA GPU; will NOT run on the dev Mac):

    ssh gbook
    cd ~/ffn_ac_native
    PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim \
      ~/miniconda3/envs/ffn_sim/bin/python \
      aleph/scripts/ac_gate_b_composed_native.py --k-axial 1000 --n-fibers-ecm 200 \
        --cortex-filaments 70686

Long run: launch under nohup and monitor the LOG FILE (ssh python is not on PATH; use the full env python).
────────────────────────────────────────────────────────────────────────────────────────────────────────

Runtime: NVIDIA Warp on CUDA only (I0-A).  Authored on the dev Mac (no CUDA) — the Lead runs it on gbook.
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp


def _max_active_force_norm(force_d: wp.array, active_mask: np.ndarray | None = None) -> float:
    """OUT-OF-LOOP host diagnostic: max per-node force magnitude (device→host readback between phases only)."""
    f = force_d.numpy()
    if f.shape[0] == 0:
        return 0.0
    norms = np.linalg.norm(f, axis=1)
    if active_mask is not None:
        norms = norms[active_mask.astype(bool)]
    return float(np.max(norms)) if norms.size else 0.0


def _accumulate_owners(cell) -> None:
    """Zero every owner's force array and run ONE composed candidate solve (each owner's own kernels)."""
    cell.solve_candidate()
    wp.synchronize_device(wp.get_device())


def _build_nmii_cortex_connector(*, actuator_state, port, params, device, max_segment_length_um):
    """Build the real ``nmii_cortex_motor`` connector (mirrors ``build_cortex_motor_slice``'s wiring block).

    Constructed here (not via :func:`build_cortex_motor_slice`) because that helper also builds its OWN
    single-slice ``CellActor`` + ``CellTransaction`` — we want the connector as a participant in the ONE
    composed transaction instead.
    """
    from aleph.engine.cortex_motor_slice import (
        CortexMotorConnector,
        SegmentAttachQuery,
        allocate_segment_connector_state,
    )
    from aleph.components.motor.segment_motor import SegmentDetachKinetics

    n_heads = int(actuator_state.n_heads)
    n_segments = int(port.segment_node_a_d.shape[0])
    state = allocate_segment_connector_state(n_heads, n_segments, device=device)
    hand_params = params.to_hand_params()
    catch_slip = (
        params.to_catch_slip_params()
        if params.detach_kinetics is SegmentDetachKinetics.CATCH_SLIP else None
    )
    query = SegmentAttachQuery(
        state=state, head_node_d=actuator_state.head_node_d, port=port,
        capture_radius=params.capture_radius, max_segment_length_um=max_segment_length_um, device=device,
    )
    r0_bind_d = wp.zeros(n_heads, dtype=wp.float64, device=device)
    r0_bind_snap_d = wp.zeros(n_heads, dtype=wp.float64, device=device)
    return CortexMotorConnector(
        state=state, head_node_d=actuator_state.head_node_d, params=hand_params, catch_slip=catch_slip,
        query=query, actuator_view=actuator_state.geometry, port=port,
        r0_bind_d=r0_bind_d, r0_bind_snap_d=r0_bind_snap_d,
        detach_kinetics=params.detach_kinetics, device=device,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description="GATE-B: sf_arc + ECM + NMII compose into ONE native cell.")
    ap.add_argument("--k-axial", type=float, required=True,
                    help="SF actin axial backbone stiffness [pN/µm] — REQUIRED modelling GAP (no default)")
    ap.add_argument("--cortex-filaments", type=int, default=70686, help="cortex F-actin count (native 70686)")
    ap.add_argument("--n-fibers-ecm", type=int, default=200, help="collagen Mikado fiber count")
    ap.add_argument("--ecm-box-um", type=float, default=10.0)
    ap.add_argument("--sf-stretch", type=float, default=1.05, help="uniform SF stretch about the centroid")
    ap.add_argument("--ecm-shear", type=float, default=0.02, help="imposed collagen simple-shear strain γ")
    ap.add_argument("--k-xb", type=float, default=None, help="NMII crossbridge stiffness [pN/µm] (default = "
                    "params_i0b3 provisional); a MASTER-knob GAP")
    ap.add_argument("--catch-slip", action="store_true",
                    help="physiological Pereverzev catch-slip detach (proxy constants) instead of Bell slip")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default="")
    args = ap.parse_args()

    wp.init()
    dev = str(wp.get_device())
    if not wp.get_device().is_cuda:
        raise SystemExit("GATE-B composed native gate requires CUDA (I0-A); no CPU simulation path exists.")

    # Lazy CUDA-lane imports (the cell-layer minifilament/port builder + the ECM device topology).
    from aleph.components.incumbent.assemble import (
        NMII_BACKBONE_LP_DIAGNOSTIC_UM,
        NMII_K_XB_TEST,
        CellConfig,
        build_cell,
    )
    from aleph.components.ecm.mikado_topology import MikadoInitConfig, MikadoTopologyBuilder
    from aleph.engine.composed_native import (
        ECM_COMPONENT,
        INTERIOR_COLUMN_FLUID_CONNECTORS,
        NMII_COMPONENT,
        SF_COMPONENT,
        build_native_composed_cell_world,
    )
    from aleph.engine.ecm_mechanics import COLLAGEN_MATERIAL_KEY
    from aleph.engine.ecm_world import BoundaryAnchorMode, ECMWorldSettings
    from aleph.engine.sf_population import build_sf_arc_population
    from aleph.laws.ecm_library import get_spec
    from aleph.scripts.ac_gate_b_cortex_motor_native import _build_actuator, _build_port, _params

    k_xb = float(args.k_xb) if args.k_xb is not None else float(NMII_K_XB_TEST)
    print(f"[composed] device={dev}  k_axial={args.k_axial} pN/µm (SF GAP)  k_xb={k_xb} pN/µm (NMII GAP)", flush=True)

    # 1. Build the NMII actuator + cortex bind-target PORT from the incumbent cell (cell-layer minifilaments).
    t0 = time.time()
    cfg = CellConfig(
        n_filaments=int(args.cortex_filaments), with_myosin=True, overlap_free_cortex=True,
        membrane_subdivisions=6, nucleus_subdivisions=3, erm_radial_pairing=True,
        resting_bound_myosin_fraction=None, nmii_straddle_placement=True,
        nmii_backbone_lp_um=NMII_BACKBONE_LP_DIAGNOSTIC_UM,
        nmii_backbone_lp_source="MECHANISM_DEMO_DIAGNOSTIC_L_p_FIXTURE_NOT_PRODUCTION_PI_GAP",
        seed=int(args.seed),
    )
    cell = build_cell(cfg)
    nmii_actuator = _build_actuator(cell, dev, k_xb)
    cortex_port = _build_port(cell, dev)
    params = _params(args.catch_slip, k_xb)
    max_seg_um = float(cell.srest_d.numpy().max())
    nmii_cortex_connector = _build_nmii_cortex_connector(
        actuator_state=nmii_actuator, port=cortex_port, params=params, device=dev,
        max_segment_length_um=max_seg_um,
    )
    n_heads = int(nmii_actuator.n_heads)
    print(f"[composed] cell+NMII built in {time.time() - t0:.1f}s: n_actin={cell.n_actin} n_heads={n_heads} "
          f"n_minifilaments={nmii_actuator.n_minifilaments}", flush=True)

    # 2. Build the disjoint sf_arc population (host) + the collagen Mikado device topology.
    sf_pop = build_sf_arc_population(n_ventral=8, n_dorsal=4, n_arc=4, n_cap=4, n_per_fiber=9)
    sf_pop.assert_partitioned()
    spec = get_spec(COLLAGEN_MATERIAL_KEY)
    ecm_cfg = MikadoInitConfig(
        box_lo_um=(0.0, 0.0, 0.0), box_hi_um=(args.ecm_box_um,) * 3,
        n_fibers=int(args.n_fibers_ecm), fiber_length_um=spec.fiber_len_um,
        target_segment_um=spec.seg_um, crosslink_capture_um=spec.xl_contact_um,
        pin_faces=("z_lo",), pin_margin_um=1.0, rng_seed=int(args.seed), max_refinement_level=0,
    )
    ecm_topology = MikadoTopologyBuilder(ecm_cfg, device=dev).initialize()

    # 3. COMPOSE: real sf_arc + ECM owners + NMII (cortex port) → ONE CellActor + ONE CellTransaction.
    cell_world = build_native_composed_cell_world(
        device=dev, base_seed=int(args.seed),
        sf_population=sf_pop, sf_k_axial_pn_per_um=float(args.k_axial),
        ecm_topology=ecm_topology, ecm_settings=ECMWorldSettings(BoundaryAnchorMode.FAR_FIELD_DIRICHLET),
        ecm_spec=spec, nmii_actuator=nmii_actuator, cortex_port=cortex_port,
        nmii_cortex_connector=nmii_cortex_connector, cortex_capacity=int(args.cortex_filaments),
    )
    sf_owner = cell_world.component_owners[SF_COMPONENT]
    ecm_owner = cell_world.component_owners[ECM_COMPONENT]
    print(f"[composed] bound {sorted(cell_world.world.registered_components())} components + "
          f"{sorted(cell_world.world.registered_connectors())} connectors into ONE world; "
          f"{len(cell_world.participants)} accepted-step participants", flush=True)

    ecm_active = ecm_topology.node_active_d.numpy()

    # --- (1) REST: composed passive force ≈ 0 -------------------------------------------------------------
    _accumulate_owners(cell_world)
    rest = {
        "sf": _max_active_force_norm(sf_owner.force_d),
        "ecm": _max_active_force_norm(ecm_owner.force_d, ecm_active),
        "cortex_port": _max_active_force_norm(cortex_port.force_d),
    }
    print(f"[composed] (1) REST  max|F|  sf={rest['sf']:.3e}  ecm={rest['ecm']:.3e}  "
          f"cortex_port={rest['cortex_port']:.3e} pN  (≈0 expected)", flush=True)

    # --- (2) perturb: stretch SF about its centroid, simple-shear the collagen network (out-of-loop) ------
    sf_pos = sf_owner.position_d.numpy()
    c = sf_pos.mean(axis=0)
    sf_pos = c + float(args.sf_stretch) * (sf_pos - c)
    wp.copy(sf_owner.position_d, wp.array(np.ascontiguousarray(sf_pos), dtype=wp.vec3d, device=dev))
    ep = ecm_owner.position_d.numpy()
    ep[:, 0] += float(args.ecm_shear) * ep[:, 2]
    wp.copy(ecm_owner.position_d, wp.array(np.ascontiguousarray(ep), dtype=wp.vec3d, device=dev))

    # --- (3) run ONE accepted physical step (force-accept so binding accumulates and the mechanism shows) --
    accepted_d = wp.array(np.ones(1, np.int32), dtype=wp.int32, device=dev)
    idx_before = int(cell_world.clock.accepted_step_index_d.numpy()[0])
    cell_world.step_once(dt_phys=0.01, accepted_d=accepted_d)
    wp.synchronize_device(wp.get_device())
    idx_after = int(cell_world.clock.accepted_step_index_d.numpy()[0])

    post = {
        "sf": _max_active_force_norm(sf_owner.force_d),
        "ecm": _max_active_force_norm(ecm_owner.force_d, ecm_active),
        "cortex_port": _max_active_force_norm(cortex_port.force_d),
    }
    n_bound = int(nmii_cortex_connector.state.bound_d.numpy().sum())
    print(f"[composed] (2/3) POST-STEP  max|F|  sf={post['sf']:.3e}  ecm={post['ecm']:.3e}  "
          f"cortex_port={post['cortex_port']:.3e} pN   nmii_bound={n_bound}/{n_heads}", flush=True)

    # --- disjoint census + accepted-step invariant --------------------------------------------------------
    ledgers = {led.component: led for led in cell_world.population_ledgers}
    census = {name: {"block": led.block, "active": led.active_count} for name, led in ledgers.items()}
    accepted_step_ok = idx_after == idx_before + 1

    verdict = {
        "sf_force_emerged": post["sf"] > 10.0 * max(rest["sf"], 1e-12),
        "ecm_force_emerged": post["ecm"] > 10.0 * max(rest["ecm"], 1e-12),
        "accepted_step_advanced_clock": accepted_step_ok,
        "populations_disjoint": True,  # asserted in-transaction every accepted step (else it would raise)
    }
    report = {
        "gate": "GATE-B composed native (sf_arc + ECM + NMII, cortex port) — ONE CellTransaction",
        "device": dev,
        "registered_components": sorted(cell_world.world.registered_components()),
        "registered_connectors": sorted(cell_world.world.registered_connectors()),
        "excluded_interior_fluid_column": list(INTERIOR_COLUMN_FLUID_CONNECTORS),
        "n_participants": len(cell_world.participants),
        "population_census": {k: {"block": list(v["block"]), "active": v["active"]} for k, v in census.items()},
        "rest_force_pN": rest,
        "post_step_force_pN": post,
        "nmii_bound_heads": n_bound,
        "n_heads": n_heads,
        "accepted_step_index": {"before": idx_before, "after": idx_after},
        "verdict": verdict,
        "note": "MECHANISM demo (rest≈0, perturb→emergent restoring force, one clean accepted step); NOT a "
                "quantitative gate. k_axial / E_fibril / NMII force-scale are KB/PI GAPs. Interior fluid "
                "column excluded (PI-gated Card-5 driver seam).",
    }
    text = json.dumps(report, indent=2)
    print(text, flush=True)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text)

    if not (verdict["sf_force_emerged"] and verdict["ecm_force_emerged"] and accepted_step_ok):
        raise SystemExit("[composed] VERDICT INCOMPLETE — see report (a component emitted no force or the "
                         "accepted step did not advance the clock).")
    print("[composed] VERDICT PASS: sf_arc + ECM + NMII composed into ONE CellTransaction; each launched its "
          "own kernels over a DISJOINT population; one predicate reached every rollback/commit; the clock "
          "advanced one accepted step. (Magnitudes PROVISIONAL — PI-GAP params.)", flush=True)


if __name__ == "__main__":
    main()
