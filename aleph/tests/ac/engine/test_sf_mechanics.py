r"""Structural gate for the KERNEL_BOUND ``sf_arc`` rod-cable mechanics (CPU-safe; native forbidden).

Proves the advance from SEAMED → KERNEL_BOUND for ``sf_arc``: the SF component now OWNS SF-local device
topology/parameter arrays and LAUNCHES the reused ``ff`` kernels over them, mirroring the cortex's
:class:`aleph.engine.surface_body.CortexFilamentMechanics`.

Two independent proofs, both CUDA-free:
  * LAUNCH — a recording launcher (metadata-only ``_FakeArray`` doubles) proves
    :meth:`SFFilamentMechanics.accumulate` launches ``link_spring_kernel`` then ``cytosim_bending_kernel`` with
    the byte-identical cortex arg order and dims ``n_links``/``n_triples`` (no kernel executes on the host);
  * PHYSICS — the pure-NumPy reference of the EXACT kernel force law proves a STRETCHED SF develops a non-zero
    axial tension and a BENT SF develops a non-zero restoring bending force over the built SF arrays (i.e. the
    binding is real, not empty scaffolding).  No Warp-CPU production-kernel launch is used.

Plus DISJOINT-population + SOURCED/GAP + connector-honesty assertions.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest
import warp as wp

from aleph.engine.population import PopulationLedger
from aleph.engine.sf_mechanics import (
    ALPHA_ACTININ_K_PN_PER_UM,
    DORSAL_ARC_CROSSLINK,
    KAPPA_ACTIN_PN_UM2,
    SF_COMPONENT,
    SF_CONNECTOR_BINDING_STATUS,
    SFFilamentMechanics,
    SFInternalArcJointConnector,
    bending_force_reference,
    build_sf_mechanics_topology,
    link_spring_force_reference,
)
from aleph.engine.sf_population import build_sf_arc_population


# ── metadata-only CUDA doubles + recording launcher (no host kernel execution) ────────────────────────
@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    ptr: int
    shape: tuple[int, ...] = (4,)
    dtype: object = wp.vec3d
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _LaunchRecorder:
    calls: list[dict[str, object]] = field(default_factory=list)

    def __call__(self, kernel, *, dim, inputs, device) -> None:
        self.calls.append({"kernel": kernel, "dim": dim, "inputs": inputs, "device": device})


_LINK_KERNEL = object()
_BENDING_KERNEL = object()


def _mechanics(recorder: _LaunchRecorder, *, n_links: int = 3, n_triples: int = 2) -> SFFilamentMechanics:
    return SFFilamentMechanics(
        device="cuda:0",
        n_nodes=8,
        links_d=_FakeArray(500, shape=(n_links, 2), dtype=wp.int32),
        link_k_d=_FakeArray(501, shape=(n_links,), dtype=wp.float64),
        link_r0_d=_FakeArray(502, shape=(n_links,), dtype=wp.float64),
        bend_triples_d=_FakeArray(503, shape=(n_triples, 3), dtype=wp.int32),
        bend_alpha_d=_FakeArray(504, shape=(n_triples,), dtype=wp.float64),
        link_kernel=_LINK_KERNEL,
        bending_kernel=_BENDING_KERNEL,
        launch=recorder,
    )


# ── LAUNCH proof: the SF component launches BOTH kernels itself, correct arg order ────────────────────
def test_sf_mechanics_launches_link_and_bending_kernels_on_sf_arrays() -> None:
    recorder = _LaunchRecorder()
    mechanics = _mechanics(recorder, n_links=3, n_triples=2)
    pos = _FakeArray(600)
    force = _FakeArray(601)
    mechanics.accumulate(pos, force)

    assert len(recorder.calls) == 2
    link_call, bend_call = recorder.calls
    assert link_call["kernel"] is _LINK_KERNEL
    assert link_call["dim"] == 3  # one thread per SF axial segment
    assert link_call["inputs"] == [pos, mechanics.links_d, mechanics.link_k_d, mechanics.link_r0_d, force]
    assert link_call["device"] == "cuda:0"
    assert bend_call["kernel"] is _BENDING_KERNEL
    assert bend_call["dim"] == 2  # one thread per SF bending triple
    assert bend_call["inputs"] == [pos, mechanics.bend_triples_d, mechanics.bend_alpha_d, force]
    # every launch targets the SF-owned force array (last input), never a foreign one
    for call in recorder.calls:
        assert call["inputs"][0] is pos
        assert call["inputs"][-1] is force


def test_sf_mechanics_validates_topology_and_rejects_cpu_or_mismatched_arrays() -> None:
    recorder = _LaunchRecorder()
    with pytest.raises(TypeError, match="dtype wp.int32"):
        SFFilamentMechanics(
            device="cuda:0", n_nodes=8,
            links_d=_FakeArray(500, shape=(3, 2), dtype=wp.float64),
            link_k_d=_FakeArray(501, shape=(3,), dtype=wp.float64),
            link_r0_d=_FakeArray(502, shape=(3,), dtype=wp.float64),
            bend_triples_d=_FakeArray(503, shape=(2, 3), dtype=wp.int32),
            bend_alpha_d=_FakeArray(504, shape=(2,), dtype=wp.float64),
            link_kernel=_LINK_KERNEL, bending_kernel=_BENDING_KERNEL, launch=recorder,
        )
    with pytest.raises(ValueError, match="match its topology row count"):
        SFFilamentMechanics(
            device="cuda:0", n_nodes=8,
            links_d=_FakeArray(500, shape=(3, 2), dtype=wp.int32),
            link_k_d=_FakeArray(501, shape=(2,), dtype=wp.float64),  # wrong length vs 3 links
            link_r0_d=_FakeArray(502, shape=(3,), dtype=wp.float64),
            bend_triples_d=_FakeArray(503, shape=(2, 3), dtype=wp.int32),
            bend_alpha_d=_FakeArray(504, shape=(2,), dtype=wp.float64),
            link_kernel=_LINK_KERNEL, bending_kernel=_BENDING_KERNEL, launch=recorder,
        )
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        SFFilamentMechanics(
            device="cuda", n_nodes=8,
            links_d=_FakeArray(500, shape=(3, 2), dtype=wp.int32, device=cpu),
            link_k_d=_FakeArray(501, shape=(3,), dtype=wp.float64, device=cpu),
            link_r0_d=_FakeArray(502, shape=(3,), dtype=wp.float64, device=cpu),
            bend_triples_d=_FakeArray(503, shape=(2, 3), dtype=wp.int32, device=cpu),
            bend_alpha_d=_FakeArray(504, shape=(2,), dtype=wp.float64, device=cpu),
            link_kernel=_LINK_KERNEL, bending_kernel=_BENDING_KERNEL, launch=recorder,
        )


# ── DISJOINT population: unique IDs, cortex-disjoint, indices address only SF nodes ───────────────────
def test_topology_is_built_over_the_disjoint_sf_population() -> None:
    pop = build_sf_arc_population(id_base=1_000_000)
    pop.assert_partitioned()
    cortex = PopulationLedger("cortex", id_base=0, capacity=70_686)
    cortex.seed_active(range(70_686))
    pop.assert_disjoint_from(cortex)  # sf_arc shares no id block with the cortex

    topo = build_sf_mechanics_topology(pop, k_axial_pn_per_um=1000.0)
    n = topo.n_nodes
    assert n == pop.n_nodes
    # every link / triple / joint index addresses an sf_arc node (never a cortex node)
    assert topo.links.min() >= 0 and topo.links.max() < n
    assert topo.bend_triples.min() >= 0 and topo.bend_triples.max() < n
    if topo.n_arc_joints:
        assert topo.arc_joints.min() >= 0 and topo.arc_joints.max() < n
    # links/triples exist for a non-trivial population
    assert topo.n_links > 0 and topo.n_triples > 0


# ── PHYSICS proof (NumPy reference of the exact kernel law): stretch → tension, bend → bending force ──
def test_stretched_sf_develops_axial_tension() -> None:
    pop = build_sf_arc_population()
    topo = build_sf_mechanics_topology(pop, k_axial_pn_per_um=1000.0)

    # at rest (build geometry == rest length), the passive backbone injects ZERO force (emergent-not-lumped)
    rest = link_spring_force_reference(pop.pos, topo.links, topo.link_k, topo.link_r0)
    assert np.abs(rest).max() < 1e-9

    # stretch one segment 50% → a non-zero restoring tension appears along that axis on both endpoints
    pos = pop.pos.copy()
    i, j = int(topo.links[0, 0]), int(topo.links[0, 1])
    axis = pos[j] - pos[i]
    pos[j] = pos[i] + 1.5 * axis
    force = link_spring_force_reference(pos, topo.links, topo.link_k, topo.link_r0)
    assert np.linalg.norm(force[i]) > 1e-6 and np.linalg.norm(force[j]) > 1e-6
    # the tension pulls the stretched endpoints back together (i outward-toward-j, j back-toward-i)
    unit = axis / np.linalg.norm(axis)
    assert np.dot(force[i], unit) > 0.0   # node i pulled toward j
    assert np.dot(force[j], unit) < 0.0   # node j pulled back toward i


def test_bent_sf_develops_restoring_bending_force() -> None:
    pop = build_sf_arc_population()
    topo = build_sf_mechanics_topology(pop, k_axial_pn_per_um=1000.0)

    # a straight fiber has zero discrete curvature → zero bending force at rest
    rest = bending_force_reference(pop.pos, topo.bend_triples, topo.bend_alpha)
    assert np.abs(rest).max() < 1e-9

    # kink the center node of the first triple → a non-zero restoring bending force appears
    pos = pop.pos.copy()
    b = int(topo.bend_triples[0, 1])
    pos[b] = pos[b] + np.array([0.0, 0.0, 0.25])
    force = bending_force_reference(pos, topo.bend_triples, topo.bend_alpha)
    assert np.linalg.norm(force[b]) > 1e-6


# ── SOURCED vs GAP discipline ─────────────────────────────────────────────────────────────────────────
def test_bending_alpha_uses_sourced_actin_kappa() -> None:
    pop = build_sf_arc_population()
    topo = build_sf_mechanics_topology(pop, k_axial_pn_per_um=1000.0, end_correction=False)
    # α = κ/seg³ with the SOURCED actin κ (Gittes 1993); reproduce it from the first fiber's segment length.
    bundle = pop.bundles[0]
    off = np.asarray(bundle.fiber_offsets)
    lo, hi = bundle.node_base + int(off[0]), bundle.node_base + int(off[1])
    seg = np.linalg.norm(np.diff(pop.pos[lo:hi], axis=0), axis=1).mean()
    expected = KAPPA_ACTIN_PN_UM2 / seg**3
    assert topo.bend_alpha[0] == pytest.approx(expected, rel=1e-9)


def test_axial_stiffness_is_a_required_gap_no_default() -> None:
    pop = build_sf_arc_population()
    with pytest.raises(ValueError, match="REQUIRED-PARAM 'k_axial_pn_per_um'"):
        build_sf_mechanics_topology(pop, k_axial_pn_per_um=None)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="REQUIRED-PARAM 'k_axial_pn_per_um'"):
        build_sf_mechanics_topology(pop, k_axial_pn_per_um=0.0)


def test_internal_arc_crosslink_uses_sourced_alpha_actinin() -> None:
    pop = build_sf_arc_population(n_dorsal=2, n_arc=2)
    topo = build_sf_mechanics_topology(pop, k_axial_pn_per_um=1000.0)
    assert topo.n_arc_joints > 0
    assert np.allclose(topo.arc_k, ALPHA_ACTININ_K_PN_PER_UM)


# ── the internal dorsal_arc_crosslink connector genuinely launches link_spring (single SF array) ─────
def test_internal_arc_connector_launches_link_spring_over_sf_joints() -> None:
    recorder = _LaunchRecorder()
    connector = SFInternalArcJointConnector(
        name=DORSAL_ARC_CROSSLINK, component_a=SF_COMPONENT, component_b=SF_COMPONENT, device="cuda:0",
        joints_d=_FakeArray(700, shape=(2, 2), dtype=wp.int32),
        joint_k_d=_FakeArray(701, shape=(2,), dtype=wp.float64),
        joint_r0_d=_FakeArray(702, shape=(2,), dtype=wp.float64),
        link_kernel=_LINK_KERNEL, launch=recorder,
    )
    pos = _FakeArray(800)
    force = _FakeArray(801)
    connector.accumulate(pos, force)
    assert len(recorder.calls) == 1
    call = recorder.calls[0]
    assert call["kernel"] is _LINK_KERNEL and call["dim"] == 2
    assert call["inputs"] == [pos, connector.joints_d, connector.joint_k_d, connector.joint_r0_d, force]


def test_internal_arc_connector_must_be_sf_internal() -> None:
    with pytest.raises(ValueError, match="INTERNAL"):
        SFInternalArcJointConnector(
            name=DORSAL_ARC_CROSSLINK, component_a=SF_COMPONENT, component_b="cortex", device="cuda:0",
            joints_d=_FakeArray(700, shape=(2, 2), dtype=wp.int32),
            joint_k_d=_FakeArray(701, shape=(2,), dtype=wp.float64),
            joint_r0_d=_FakeArray(702, shape=(2,), dtype=wp.float64),
            link_kernel=_LINK_KERNEL, launch=_LaunchRecorder(),
        )


# ── honesty: only the two available-target connectors are (or become) bound; the rest are SEAMED ─────
def test_connector_binding_status_is_honest() -> None:
    assert SF_CONNECTOR_BINDING_STATUS["dorsal_arc_crosslink"].startswith("KERNEL_BOUND")
    assert SF_CONNECTOR_BINDING_STATUS["sf_cortex_transient"].startswith("SEAMED")
    assert SF_CONNECTOR_BINDING_STATUS["if_sf_plectin"].startswith("SEAMED")
    assert SF_CONNECTOR_BINDING_STATUS["mt_sf_spectraplakin"].startswith("SEAMED")
