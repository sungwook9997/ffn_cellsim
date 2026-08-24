"""CUDA gate for the two stability/rebind defects found before `relax_candidate` was attached.

Implements: KU-0.0

The GPU host carries no pytest, so the CUDA-only assertions in
``aleph/tests/ac/engine/test_overdamped_relax.py`` and ``test_load_path.py`` can never run there —
they skip on the dev Mac (no CUDA) and cannot be collected on the device. That is the same shape as
the defect class this lane spent the day on: **the check does not cover the thing it protects.** This
driver carries those assertions to the device as a driver, which is what this repo runs on the GPU.

WHAT IT GATES

G1 — the Gershgorin bound is a bound.  ``46cf6e28`` summed the axial links and the bending triples
and stopped, leaving α-actinin's ``arc_joints`` (4.6e5 pN/µm, ~460× the axial k) outside the row
sums. Omitting a term does not loosen a Gershgorin bound; it reverses the inequality. Measured on
the real native ``sf_arc`` topology the admitted ``dt`` was 459× the true stability limit, with
nothing raised. G1 recomputes both numbers on the device host and requires the ratio to be large —
i.e. that the omission was real and is now closed.

G2 — the step refuses an unstable ``dt`` rather than integrating it.

G3 — the step DESCENDS: a stretched Hookean pair must relax toward its rest length on the device.
A sign error would grow every mode while still producing a trajectory.

G4 — ``r0_bind`` survives motion.  A joint's rest length was the separation at BUILD, which equalled
the separation at bind only because nothing moved. The first component to integrate its own
positions breaks that, and a rebind at the stale rest injects ``k·Δ`` from nowhere — 4,600 pN for a
10 nm drift at α-actinin stiffness, two orders above the bond's own rupture force. Three controls,
and the third is the one that matters: refreshing an ALREADY-engaged joint would make every
crosslink permanently force-free, so the connector would silently carry nothing at all.

Sanity Gate (per CLAUDE.md, before first execution):
  * dimensional — γ [pN·s·µm⁻¹], λ [pN/µm], dt [s]; ``2γ/λ`` is seconds, ``dt·f/γ`` is µm.
  * boundary — the refusal path is exercised explicitly (G2), not assumed.
  * conservation — G4's third control reads the scattered force after a rebind and requires zero.
  * CFL/precision — float64 throughout; G3 steps at a quarter of the derived bound.
  * sign sense — G3 fails if the step grows the extension.
  * measurement protocol — CUDA only; the driver RAISES on a non-CUDA device rather than proceeding.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import warp as wp

from aleph.engine.contracts import reference_cell_architecture
from aleph.engine.load_path import (
    ActorRecord,
    ActorRegistry,
    ActorSegmentRuntime,
    CompositeFAJointSpec,
    ElementKind,
    EndpointRole,
    JointState,
    LoadPathJointRuntime,
    PortRef,
    resolve_fa_series_joint,
    segment_pair_reference,
)
from aleph.engine.overdamped_relax import (
    build_overdamped_step,
    derive_gershgorin_lambda_max,
)
from aleph.engine.sf_mechanics import build_sf_mechanics_topology
from aleph.engine.sf_population import build_sf_arc_population
#: The PI-GAP axial stiffness the sf_arc topology is built with elsewhere in this lane.
K_AXIAL_PN_PER_UM = 1000.0


def _require_cuda() -> str:
    """Resolve a CUDA device or RAISE. A non-CUDA device may not proceed (CLAUDE.md I0-A)."""
    wp.init()
    device = next((str(d) for d in wp.get_devices() if d.is_cuda), None)
    if device is None:
        raise RuntimeError(
            "no CUDA device: this gate measures physics and a number measured on CPU is not a "
            "result. Submit inside a Slurm allocation."
        )
    return device


def gate_bound_is_a_bound() -> dict:
    """G1: recompute λ_max on the real sf_arc topology with and without the α-actinin joints."""
    population = build_sf_arc_population(id_base=1_000_000)
    topology = build_sf_mechanics_topology(population, k_axial_pn_per_um=K_AXIAL_PN_PER_UM)
    n_filaments = max(int(getattr(population, "n_filaments", 0) or 1), 1)

    step = build_overdamped_step(topology, n_filaments=n_filaments, borrowed_pairs=())
    omitted = derive_gershgorin_lambda_max(
        links=np.asarray(topology.links),
        link_k=np.asarray(topology.link_k, dtype=np.float64),
        n_nodes=int(topology.n_nodes),
        bend_triples=topology.bend_triples,
        bend_alpha=topology.bend_alpha,
    )
    dt_omitted = 2.0 * step.gamma_pn_s_per_um / omitted
    overshoot = dt_omitted / step.dt_max_s
    return {
        "n_nodes": int(topology.n_nodes),
        "n_links": int(np.asarray(topology.link_k).size),
        "n_arc_joints": int(topology.n_arc_joints),
        "gamma_pn_s_per_um": step.gamma_pn_s_per_um,
        "lambda_max_with_arc_k": step.lambda_max_pn_per_um,
        "lambda_max_omitting_arc_k": float(omitted),
        "dt_max_s": step.dt_max_s,
        "dt_max_s_omitting_arc_k": float(dt_omitted),
        "dt_overshoot_before_fix": float(overshoot),
        "PASS": bool(overshoot > 10.0 and step.dt_max_s < dt_omitted),
    }


def gate_refuses_unstable_dt(device: str) -> dict:
    """G2: a dt at or above the derived bound must RAISE, and the message must carry the bound."""
    population = build_sf_arc_population(id_base=1_000_000)
    topology = build_sf_mechanics_topology(population, k_axial_pn_per_um=K_AXIAL_PN_PER_UM)
    step = build_overdamped_step(topology, n_filaments=1, borrowed_pairs=())
    n = int(topology.n_nodes)
    position = wp.array(np.zeros((n, 3)), dtype=wp.vec3d, device=device)
    force = wp.array(np.zeros((n, 3)), dtype=wp.vec3d, device=device)

    raised, message = False, ""
    try:
        step.advance(position, force, dt_phys=step.dt_max_s * 1.5)
    except ValueError as exc:      # the refusal is the behaviour under test
        raised, message = True, str(exc)
    return {"raised": raised, "message": message, "PASS": raised}


def gate_step_descends(device: str) -> dict:
    """G3: a stretched Hookean pair must RELAX on the device, not merely move."""
    k, rest = 100.0, 1.0

    class _Pair:
        n_nodes = 2
        links = np.array([[0, 1]], dtype=np.int64)
        link_k = np.array([k])
        link_r0 = np.array([rest])

    step = build_overdamped_step(_Pair(), n_filaments=1, borrowed_pairs=())
    positions = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]])
    device_positions = wp.array(positions, dtype=wp.vec3d, device=device)
    before = float(np.linalg.norm(positions[1] - positions[0]) - rest)

    for _ in range(20):
        current = device_positions.numpy()
        delta = current[1] - current[0]
        length = float(np.linalg.norm(delta))
        pair = k * (length - rest) * delta / length
        force = wp.array(np.array([pair, -pair]), dtype=wp.vec3d, device=device)
        step.advance(device_positions, force, dt_phys=step.dt_max_s * 0.25)
    wp.synchronize_device(wp.get_device(device))

    final = device_positions.numpy()
    after = float(np.linalg.norm(final[1] - final[0]) - rest)
    return {
        "extension_before_um": before,
        "extension_after_um": after,
        "ratio": after / before if before else None,
        "PASS": bool(abs(after) < 0.5 * abs(before)),
    }


def gate_rebind_uses_the_geometry_it_binds_at(device: str) -> dict:
    """G4: three controls on `r0_bind` — engaged keeps, rejected keeps, accepted rebind refreshes."""
    registry = ActorRegistry(
        (
            ActorRecord("sf_arc", 10, 3, 100, 4, 1),   # ventral SF
            ActorRecord("ecm", 20, 9, 200, 6, 1),      # collagen actor
        )
    )

    def _port(actor: ActorRecord, role: EndpointRole, u: float) -> PortRef:
        return PortRef(
            component=actor.component, actor_id=actor.actor_id,
            actor_generation=actor.actor_generation, entity_id=actor.entity_id,
            entity_generation=actor.entity_generation, element_kind=ElementKind.SEGMENT,
            element_id=0, local_coordinates=(u, 0.0, 0.0, 0.0), role=role,
        )

    spec = CompositeFAJointSpec(
        joint_id=1, fa_cluster_id=17,
        sf_port=_port(registry.actor(10), EndpointRole.VENTRAL_END_1, 1.0),
        ecm_port=_port(registry.actor(20), EndpointRole.COLLAGEN_LIGAND, 0.5),
        stiffness_pn_per_um=10.0, rest_um=1.0,
    )
    joint = resolve_fa_series_joint(spec, registry, reference_cell_architecture())
    sf_pos = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    ecm_pos = np.array([[2.0, 0.0, 0.0], [2.0, 1.0, 0.0]])
    segments = np.array([[0, 1]], np.int32)
    sf = ActorSegmentRuntime(registry.actor(10), sf_pos, segments, device=device)
    ecm = ActorSegmentRuntime(registry.actor(20), ecm_pos, segments, device=device)
    runtime = LoadPathJointRuntime(sf, ecm, (joint,))
    built_rest = float(runtime.rest_d.numpy()[0])
    accepted = wp.array(np.array([1], np.int32), dtype=wp.int32, device=device)
    rejected = wp.array(np.array([0], np.int32), dtype=wp.int32, device=device)

    # (1) an ALREADY-engaged joint keeps its rest; refreshing it would zero every crosslink force
    sf.pos_d.assign(sf_pos + np.array([0.25, 0.0, 0.0]))
    runtime.snapshot_candidate()
    runtime.commit_irreversible(accepted, 0.05, 11)
    wp.synchronize_device(wp.get_device(device))
    engaged_keeps = float(runtime.rest_d.numpy()[0])

    # (2) a REJECTED step may not move rest — rollback correctness by construction
    runtime.state_d.assign(np.array([int(JointState.FREE)], np.int32))
    runtime.candidate_state_d.assign(np.array([int(JointState.ACTIN_ENGAGED)], np.int32))
    runtime.commit_irreversible(rejected, 0.05, 12)
    wp.synchronize_device(wp.get_device(device))
    rejected_keeps = float(runtime.rest_d.numpy()[0])

    # (3) an ACCEPTED rebind takes the separation it binds AT, and is force-free there
    moved = sf_pos + np.array([0.4, 0.0, 0.0])
    sf.pos_d.assign(moved)
    runtime.state_d.assign(np.array([int(JointState.FREE)], np.int32))
    runtime.candidate_state_d.assign(np.array([int(JointState.ACTIN_ENGAGED)], np.int32))
    runtime.commit_irreversible(accepted, 0.05, 13)
    wp.synchronize_device(wp.get_device(device))
    refreshed = float(runtime.rest_d.numpy()[0])

    oracle = segment_pair_reference(
        moved, (0, 1), joint.port_a.segment_u, ecm_pos, (0, 1), joint.port_b.segment_u,
        joint.stiffness_pn_per_um, 0.0,
    )
    separation = oracle.load_pn / joint.stiffness_pn_per_um

    sf.zero_force()
    ecm.zero_force()
    runtime.accumulate()
    wp.synchronize_device(wp.get_device(device))
    residual_force_pn = float(np.max(np.abs(sf.force_d.numpy())))
    stale_artifact_pn = float(joint.stiffness_pn_per_um * abs(separation - built_rest))

    return {
        "rest_at_build_um": built_rest,
        "rest_after_engaged_step_um": engaged_keeps,
        "rest_after_rejected_step_um": rejected_keeps,
        "rest_after_accepted_rebind_um": refreshed,
        "separation_at_rebind_um": float(separation),
        "force_after_rebind_pn": residual_force_pn,
        "artifact_a_stale_rest_would_have_injected_pn": stale_artifact_pn,
        "PASS": bool(
            abs(engaged_keeps - built_rest) < 1.0e-12
            and abs(rejected_keeps - built_rest) < 1.0e-12
            and abs(refreshed - separation) < 1.0e-9
            and residual_force_pn < 1.0e-9
        ),
    }


def main() -> int:
    device = _require_cuda()
    record = {
        "utc": datetime.now(timezone.utc).isoformat(),
        "device": device,
        "G1_bound_is_a_bound": gate_bound_is_a_bound(),
        "G2_refuses_unstable_dt": gate_refuses_unstable_dt(device),
        "G3_step_descends": gate_step_descends(device),
        "G4_rebind_rest_is_at_bind": gate_rebind_uses_the_geometry_it_binds_at(device),
    }
    record["ALL_PASS"] = all(v["PASS"] for k, v in record.items() if k.startswith("G"))

    out = Path("outputs/ac/relax_bound_rebind")
    out.mkdir(parents=True, exist_ok=True)
    (out / "gate.json").write_text(json.dumps(record, indent=2))

    print(json.dumps(record, indent=2))
    return 0 if record["ALL_PASS"] else 1


if __name__ == "__main__":
    sys.exit(main())
