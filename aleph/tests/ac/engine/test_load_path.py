"""Gates for the VS-0 SF--FA--ECM Active Load-Path Graph seam."""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorContract,
    ConnectorFamily,
    reference_cell_architecture,
)
from aleph.engine.load_path import (
    ActorRecord,
    ActorRegistry,
    ActorSegmentRuntime,
    CompositeFAJointSpec,
    DorsalArcJointSpec,
    ElementKind,
    EndpointRole,
    JointState,
    LoadPathJointRuntime,
    PortRef,
    resolve_dorsal_arc_joint,
    resolve_fa_series_group,
    resolve_fa_series_joint,
    segment_pair_reference,
)


def _registry() -> ActorRegistry:
    return ActorRegistry(
        (
            ActorRecord("sf_arc", 10, 3, 100, 4, 1),       # ventral SF
            ActorRecord("sf_arc", 11, 5, 101, 2, 1),       # dorsal SF
            ActorRecord("sf_arc", 12, 7, 102, 8, 1),       # transverse arc
            ActorRecord("ecm", 20, 9, 200, 6, 1),          # collagen actor
        )
    )


def _port(actor: ActorRecord, role: EndpointRole, u: float = 0.5) -> PortRef:
    return PortRef(
        component=actor.component,
        actor_id=actor.actor_id,
        actor_generation=actor.actor_generation,
        entity_id=actor.entity_id,
        entity_generation=actor.entity_generation,
        element_kind=ElementKind.SEGMENT,
        element_id=0,
        local_coordinates=(u, 0.0, 0.0, 0.0),
        role=role,
    )


def _fa_spec(registry: ActorRegistry) -> CompositeFAJointSpec:
    return CompositeFAJointSpec(
        joint_id=1,
        fa_cluster_id=17,
        sf_port=_port(registry.actor(10), EndpointRole.VENTRAL_END_1, 1.0),
        ecm_port=_port(registry.actor(20), EndpointRole.COLLAGEN_LIGAND, 0.5),
        stiffness_pn_per_um=10.0,
        rest_um=1.0,
    )


def test_geometryless_fa_edges_resolve_to_one_mechanical_joint() -> None:
    registry = _registry()
    architecture = reference_cell_architecture()
    resolution = resolve_fa_series_group(architecture)
    joint = resolve_fa_series_joint(_fa_spec(registry), registry, architecture)
    assert resolution.mechanical_joint_count == 1
    assert set(joint.semantic_edges) == {"fa_actin_anchor", "integrin_collagen_clutch"}
    assert joint.stiffness_pn_per_um == 10.0   # one stiffness, not two semantic-edge springs


def test_double_stiffness_group_is_rejected() -> None:
    base = reference_cell_architecture()
    duplicate = ConnectorContract(
        "duplicate_fa_actin",
        ConnectorFamily.ACTIN_ANCHOR,
        "sf_arc",
        "focal_adhesion",
        True,
        True,
        chemistry_card="actin_talin_integrin",
        mechanical_group="alpha2beta1_collagen_series",
    )
    architecture = CellArchitecture(base.components, (*base.connectors, duplicate))
    with pytest.raises(ValueError, match="exactly"):
        resolve_fa_series_group(architecture)


@pytest.mark.parametrize(
    "role",
    [
        EndpointRole.VENTRAL_END_0,
        EndpointRole.VENTRAL_END_1,
        EndpointRole.DORSAL_BASAL_END,
    ],
)
def test_only_load_bearing_sf_roles_may_enter_fa(role: EndpointRole) -> None:
    registry = _registry()
    spec = CompositeFAJointSpec(
        joint_id=1,
        fa_cluster_id=0,
        sf_port=_port(registry.actor(10), role),
        ecm_port=_port(registry.actor(20), EndpointRole.COLLAGEN_LIGAND),
        stiffness_pn_per_um=10.0,
        rest_um=0.0,
    )
    assert resolve_fa_series_joint(spec, registry, reference_cell_architecture()).port_a.role is role


@pytest.mark.parametrize(
    "role",
    [EndpointRole.DORSAL_FREE_END, EndpointRole.TRANSVERSE_ARC_MATERIAL],
)
def test_dorsal_free_and_arc_roles_cannot_bind_fa(role: EndpointRole) -> None:
    registry = _registry()
    actor = registry.actor(11 if role is EndpointRole.DORSAL_FREE_END else 12)
    spec = CompositeFAJointSpec(
        joint_id=1,
        fa_cluster_id=0,
        sf_port=_port(actor, role),
        ecm_port=_port(registry.actor(20), EndpointRole.COLLAGEN_LIGAND),
        stiffness_pn_per_um=10.0,
        rest_um=0.0,
    )
    with pytest.raises(ValueError, match="ventral end or dorsal basal"):
        resolve_fa_series_joint(spec, registry, reference_cell_architecture())


def test_stale_actor_generation_is_rejected() -> None:
    registry = _registry()
    spec = _fa_spec(registry)
    stale = PortRef(
        component=spec.sf_port.component,
        actor_id=spec.sf_port.actor_id,
        actor_generation=spec.sf_port.actor_generation + 1,
        entity_id=spec.sf_port.entity_id,
        entity_generation=spec.sf_port.entity_generation,
        element_kind=spec.sf_port.element_kind,
        element_id=spec.sf_port.element_id,
        local_coordinates=spec.sf_port.local_coordinates,
        role=spec.sf_port.role,
    )
    stale_spec = CompositeFAJointSpec(
        joint_id=spec.joint_id,
        fa_cluster_id=spec.fa_cluster_id,
        sf_port=stale,
        ecm_port=spec.ecm_port,
        stiffness_pn_per_um=spec.stiffness_pn_per_um,
        rest_um=spec.rest_um,
    )
    with pytest.raises(ValueError, match="stale actor generation"):
        resolve_fa_series_joint(stale_spec, registry, reference_cell_architecture())


def test_internal_dorsal_arc_connector_has_typed_roles() -> None:
    registry = _registry()
    spec = DorsalArcJointSpec(
        joint_id=9,
        dorsal_port=_port(registry.actor(11), EndpointRole.DORSAL_FREE_END, 1.0),
        arc_port=_port(registry.actor(12), EndpointRole.TRANSVERSE_ARC_MATERIAL, 0.25),
        stiffness_pn_per_um=6.0,
        rest_um=0.2,
    )
    joint = resolve_dorsal_arc_joint(spec, registry, reference_cell_architecture())
    assert joint.semantic_edges == ("dorsal_arc_crosslink",)
    assert joint.port_a.component == joint.port_b.component == "sf_arc"

    bad = DorsalArcJointSpec(
        joint_id=10,
        dorsal_port=spec.dorsal_port,
        arc_port=_port(registry.actor(10), EndpointRole.VENTRAL_END_0),
        stiffness_pn_per_um=6.0,
        rest_um=0.2,
    )
    with pytest.raises(ValueError, match="transverse-arc"):
        resolve_dorsal_arc_joint(bad, registry, reference_cell_architecture())


def test_segment_scatter_closes_force_moment_and_virtual_work() -> None:
    pa = np.array([[0.0, -0.5, 0.0], [0.0, 0.5, 0.0]])
    pb = np.array([[2.0, -1.0, 0.0], [2.0, 1.0, 0.0]])
    da = np.array([[0.01, 0.02, 0.00], [-0.02, 0.01, 0.01]])
    db = np.array([[0.00, -0.01, 0.02], [0.03, 0.02, -0.01]])
    result = segment_pair_reference(
        pa,
        (0, 1),
        0.25,
        pb,
        (0, 1),
        0.75,
        8.0,
        1.0,
        displacement_a=da,
        displacement_b=db,
    )
    total_force = result.force_a.sum(axis=0) + result.force_b.sum(axis=0)
    total_moment = (
        np.cross(pa, result.force_a).sum(axis=0)
        + np.cross(pb, result.force_b).sum(axis=0)
    )
    np.testing.assert_allclose(total_force, 0.0, atol=1.0e-12)
    np.testing.assert_allclose(total_moment, 0.0, atol=1.0e-12)

    eps = 1.0e-7
    plus = segment_pair_reference(
        pa + eps * da, (0, 1), 0.25, pb + eps * db, (0, 1), 0.75, 8.0, 1.0
    )
    minus = segment_pair_reference(
        pa - eps * da, (0, 1), 0.25, pb - eps * db, (0, 1), 0.75, 8.0, 1.0
    )
    d_energy = (plus.energy_pn_um - minus.energy_pn_um) / (2.0 * eps)
    assert result.nodal_virtual_work_pn_um == pytest.approx(-d_energy, rel=5.0e-7, abs=1.0e-9)


def test_two_semantic_edges_do_not_double_force() -> None:
    pa = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    pb = np.array([[2.0, 0.0, 0.0], [2.0, 1.0, 0.0]])
    result = segment_pair_reference(pa, (0, 1), 0.5, pb, (0, 1), 0.5, 10.0, 1.0)
    assert result.load_pn == pytest.approx(10.0)       # k * (2 - 1), not 20 pN
    assert result.energy_pn_um == pytest.approx(5.0)   # 1/2 k extension^2, once


_CUDA_DEVICE = next((device for device in wp.get_devices() if device.is_cuda), None)


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="Warp load-path runtime requires CUDA")
def test_cuda_soa_force_and_accepted_kinetics() -> None:
    registry = _registry()
    joint = resolve_fa_series_joint(
        _fa_spec(registry), registry, reference_cell_architecture()
    )
    device = str(_CUDA_DEVICE)
    sf_pos = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    ecm_pos = np.array([[2.0, 0.0, 0.0], [2.0, 1.0, 0.0]])
    sf_seg = np.array([[0, 1]], np.int32)
    ecm_seg = np.array([[0, 1]], np.int32)
    sf = ActorSegmentRuntime(registry.actor(10), sf_pos, sf_seg, device=device)
    ecm = ActorSegmentRuntime(registry.actor(20), ecm_pos, ecm_seg, device=device)
    runtime = LoadPathJointRuntime(sf, ecm, (joint,))
    sf.zero_force()
    ecm.zero_force()
    runtime.accumulate()
    wp.synchronize_device(device)

    # The device force scatter must reproduce the validated host virtual-work oracle
    # for the *actual* composite-series geometry: the ventral end sits at u=1 (node
    # (0,1,0)) and the collagen ligand at u=0.5 (node (2,0.5,0)), so the spring is not
    # axis-aligned.  Comparing against segment_pair_reference (a proven single-count
    # oracle) checks the interpolated gather/scatter weights and the one-spring energy
    # without hard-coding an axis-aligned magic number.
    oracle = segment_pair_reference(
        sf_pos,
        (int(sf_seg[0, 0]), int(sf_seg[0, 1])),
        joint.port_a.segment_u,
        ecm_pos,
        (int(ecm_seg[0, 0]), int(ecm_seg[0, 1])),
        joint.port_b.segment_u,
        joint.stiffness_pn_per_um,
        joint.rest_um,
    )
    np.testing.assert_allclose(sf.force_d.numpy(), oracle.force_a, atol=1.0e-12)
    np.testing.assert_allclose(ecm.force_d.numpy(), oracle.force_b, atol=1.0e-12)
    # Newton's third law across the two component owners: the ungrounded pair scatters
    # zero net force even though each side's interpolation weights differ.
    net_force = sf.force_d.numpy().sum(axis=0) + ecm.force_d.numpy().sum(axis=0)
    np.testing.assert_allclose(net_force, 0.0, atol=1.0e-12)
    assert runtime.load_d.numpy()[0] == pytest.approx(oracle.load_pn)
    assert runtime.energy_d.numpy()[0] == pytest.approx(oracle.energy_pn_um)

    runtime.candidate_state_d.assign(np.array([int(JointState.FREE)], np.int32))
    rejected = wp.array(np.array([0], np.int32), dtype=wp.int32, device=device)
    runtime.commit_irreversible(rejected, 0.05, 123)
    wp.synchronize_device(device)
    assert runtime.state_d.numpy()[0] == int(JointState.ACTIN_ENGAGED)
    assert runtime.candidate_state_d.numpy()[0] == int(JointState.ACTIN_ENGAGED)
    assert runtime.rng_epoch_d.numpy()[0] == 0
    assert runtime.age_d.numpy()[0] == 0.0

    runtime.candidate_state_d.assign(np.array([int(JointState.FREE)], np.int32))
    accepted = wp.array(np.array([1], np.int32), dtype=wp.int32, device=device)
    runtime.commit_irreversible(accepted, 0.05, 123)
    wp.synchronize_device(device)
    assert runtime.state_d.numpy()[0] == int(JointState.FREE)
    assert runtime.rng_epoch_d.numpy()[0] == 1
    assert runtime.age_d.numpy()[0] == 0.0       # state lifetime resets on transition

    runtime.snapshot_candidate()
    runtime.commit_irreversible(accepted, 0.05, 124)
    wp.synchronize_device(device)
    assert runtime.rng_epoch_d.numpy()[0] == 2
    assert runtime.age_d.numpy()[0] == pytest.approx(0.05)

    sf.zero_force()
    ecm.zero_force()
    runtime.accumulate()
    wp.synchronize_device(device)
    assert np.max(np.abs(sf.force_d.numpy())) == 0.0
    assert np.max(np.abs(ecm.force_d.numpy())) == 0.0


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="Warp load-path runtime requires CUDA")
def test_a_rebind_takes_its_rest_from_the_geometry_it_binds_AT() -> None:
    """`r0_bind` held trivially while nothing moved. It stops holding the day a component integrates.

    Build-time separation was the separation at bind because positions never changed. Once
    `overdamped_relax` advances a component's own nodes the two diverge, and a rebind at the stale
    rest injects `k·Δ` from nowhere — at α-actinin's 4.6e5 pN/µm a 10 nm drift is 4,600 pN, two
    orders above the bond's own rupture force. Nothing raises; the cell is pulled by an artifact.

    Three controls, and the third is the one that matters most: refreshing an ALREADY-engaged joint
    would make every crosslink permanently force-free, i.e. the connector would silently carry
    nothing at all — the vacuous PASS this lane keeps having to design against.
    """
    registry = _registry()
    joint = resolve_fa_series_joint(_fa_spec(registry), registry, reference_cell_architecture())
    device = str(_CUDA_DEVICE)
    sf_pos = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    ecm_pos = np.array([[2.0, 0.0, 0.0], [2.0, 1.0, 0.0]])
    seg = np.array([[0, 1]], np.int32)
    sf = ActorSegmentRuntime(registry.actor(10), sf_pos, seg, device=device)
    ecm = ActorSegmentRuntime(registry.actor(20), ecm_pos, seg, device=device)
    runtime = LoadPathJointRuntime(sf, ecm, (joint,))
    built_rest = float(runtime.rest_d.numpy()[0])
    accepted = wp.array(np.array([1], np.int32), dtype=wp.int32, device=device)
    rejected = wp.array(np.array([0], np.int32), dtype=wp.int32, device=device)

    # (1) an ALREADY-engaged joint keeps its rest — otherwise every crosslink goes force-free
    sf.pos_d.assign(sf_pos + np.array([0.25, 0.0, 0.0]))
    runtime.snapshot_candidate()
    runtime.commit_irreversible(accepted, 0.05, 11)
    wp.synchronize_device(device)
    assert runtime.rest_d.numpy()[0] == pytest.approx(built_rest), "an engaged joint owns its rest"

    # (2) a REJECTED step leaves rest untouched, so rollback stays correct by construction
    runtime.state_d.assign(np.array([int(JointState.FREE)], np.int32))
    runtime.candidate_state_d.assign(np.array([int(JointState.ACTIN_ENGAGED)], np.int32))
    runtime.commit_irreversible(rejected, 0.05, 12)
    wp.synchronize_device(device)
    assert runtime.rest_d.numpy()[0] == pytest.approx(built_rest), "a rejected step may not move rest"

    # (3) an ACCEPTED rebind takes the separation it actually binds at — and is force-free there
    moved = sf_pos + np.array([0.4, 0.0, 0.0])
    sf.pos_d.assign(moved)
    runtime.state_d.assign(np.array([int(JointState.FREE)], np.int32))
    runtime.candidate_state_d.assign(np.array([int(JointState.ACTIN_ENGAGED)], np.int32))
    runtime.commit_irreversible(accepted, 0.05, 13)
    wp.synchronize_device(device)

    expected = segment_pair_reference(
        moved, (0, 1), joint.port_a.segment_u, ecm_pos, (0, 1), joint.port_b.segment_u,
        joint.stiffness_pn_per_um, 0.0,
    )
    separation = expected.load_pn / joint.stiffness_pn_per_um   # k·(L − 0) / k = L
    assert runtime.rest_d.numpy()[0] == pytest.approx(separation)
    assert runtime.rest_d.numpy()[0] != pytest.approx(built_rest), "the geometry did move"

    sf.zero_force()
    ecm.zero_force()
    runtime.accumulate()
    wp.synchronize_device(device)
    assert np.max(np.abs(sf.force_d.numpy())) == pytest.approx(0.0, abs=1.0e-9), (
        "a joint that just bound must start force-free; a stale rest would show k·Δ here"
    )


def test_fa_double_count_guard_can_actually_fire() -> None:
    """A rival connector across a group component pair must trip the double-stiffness guard.

    Regression for a guard that was structurally unable to fail: `mechanical_joint_count` was a
    dataclass field defaulting to 1 that `resolve_fa_series_group` never set, so the downstream
    `!= 1` check could not raise under any architecture. The existing double-stiffness test adds a
    THIRD group member, which trips the earlier `len(connectors) != 2` check instead — so the count
    guard itself had no coverage at all.
    """
    registry = _registry()
    base = reference_cell_architecture()
    # a second, independent sf_arc <-> focal_adhesion mechanical path, declared OUTSIDE the group
    rival = ConnectorContract(
        "rival_fa_actin_path",
        ConnectorFamily.ACTIN_ANCHOR,
        "sf_arc",
        "focal_adhesion",
        True,
        True,
        chemistry_card="actin_talin_integrin",
    )
    architecture = CellArchitecture(base.components, (*base.connectors, rival))

    resolution = resolve_fa_series_group(architecture)
    assert resolution.mechanical_joint_count == 2       # the composite + the rival path
    with pytest.raises(RuntimeError, match="double-count mechanical stiffness"):
        resolve_fa_series_joint(_fa_spec(registry), registry, architecture)


def test_nascent_fa_edges_are_not_counted_as_rival_paths() -> None:
    """Connectors that touch focal_adhesion across a DIFFERENT pair are separate biology, not a double count."""
    architecture = reference_cell_architecture()
    touching = {
        connector.name
        for connector in architecture.connectors
        if "focal_adhesion" in (connector.component_a, connector.component_b)
    }
    assert {"lamellipodium_nascent_fa", "filopodium_nascent_fa"} <= touching
    assert resolve_fa_series_group(architecture).mechanical_joint_count == 1
