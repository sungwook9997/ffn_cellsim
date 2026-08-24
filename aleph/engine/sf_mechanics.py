r"""``sf_arc``-owned active rod-cable mechanics — the SF component LAUNCHES its own Warp kernels.

WHAT THIS CLOSES.  :mod:`aleph.engine.sf_population` builds the DISJOINT ``sf_arc`` F-actin population
(unique global filament IDs, cortex-disjoint, host NumPy only) and :mod:`aleph.engine.stress_fiber` is the
composition facade — but until now NOTHING launched SF rod/cable mechanics anywhere: the facade's ``mechanics``
slot was a bare :class:`~aleph.engine.runtime.MechanicsContributor` Protocol filled only by test spies, so
``sf_arc`` sat at SEAMED (a population + connector CONTRACTS, no CUDA_UNIT).  This module is the concrete
delegate that makes ``sf_arc`` genuinely **KERNEL_BOUND**: it OWNS SF-local device topology/parameter arrays and
LAUNCHES the already-audited passive-fiber kernels over them, exactly the way the cortex does through
:class:`aleph.engine.surface_body.CortexFilamentMechanics`.

THE TEMPLATE IT MIRRORS (cortex → SF, byte-identical launch signatures).  :class:`SFFilamentMechanics` is the
SF counterpart of :class:`~aleph.engine.surface_body.CortexFilamentMechanics`.  It launches:

  * ``ff.network_warp.link_spring_kernel`` — Hookean axial backbone tension ``k·(L−r0)/L·d`` on the SF
    segment pairs (arg order ``[pos, links_d, link_k_d, link_r0_d, force]``); and
  * ``ff.forces_warp.cytosim_bending_kernel`` — NF2007 discrete bending ``F = α·(m_{i-1}−2m_i+m_{i+1})`` on the
    SF consecutive-node triples (arg order ``[pos, bend_triples_d, bend_alpha_d, force]``), with
    ``α = κ/seg³`` precomputed ONCE at build (never recomputed in the loop).

against the *SF-owned* ``position_d``/``force_d`` supplied by the state owner.  Both kernels ``wp.atomic_add``
into ``force``; the caller owns force zeroing.  It never allocates or touches the cortex arrays and never
concatenates SF nodes into a whole-cell array — the SF population is a SEPARATE component owning a DISJOINT
filament population (PI 2026-07-22), so this table contains only ``sf_arc``'s own segments/triples; joints to
other components go through the connector graph, never a shared row here.

EMERGENT, NOT LUMPED (CLAUDE.md hard rule).  These are the PASSIVE backbone/bending laws.  The link rest length
``r0`` is the build-time segment length, so a relaxed SF injects ZERO axial force at rest — a **stretched** SF
develops restoring tension and a **bent** SF develops restoring bending force.  The ACTIVE contractile prestress
is NOT a lumped ``k_SF`` baked in here; it is the reaction of the discrete NMII heads, which enters through a
separate MOTOR connector (the SF counterpart of ``nmii_cortex_motor``) exactly as the cortex does.  **That motor
connector LANDED on 2026-07-25** — :mod:`aleph.engine.sf_motor_slice` binds ``nmii_sf_motor`` over a real
straddle-placed NMII population and a live SF bind-target port, native-validated (SF axial tension and inward FA
traction emerge from ``k_on`` binding events, from zero).  See :data:`SF_CONNECTOR_BINDING_STATUS`.

SOURCED vs GAP (report-not-tune).
  * SOURCED — bending modulus ``κ = KAPPA_ACTIN = k_B·T·ℓ_p ≈ 0.0728 pN·µm²`` (actin ℓ_p = 17 µm, Gittes et al.
    1993, KU-1.1); the SF filaments ARE actin, so ``α = κ/seg³`` is sourced, not chosen.  SOURCED — the internal
    dorsal↔arc crosslink stiffness ``= ALPHA_ACTININ.link_k = 4.6e5 pN/µm`` (Ferrer 2008 PNAS AFM, PI-approved
    2026-06-30).
  * GAP (required, no default) — the actin **axial backbone** stretch stiffness ``k_axial_pn_per_um``.  NF2007
    treats the actin backbone as INEXTENSIBLE (a hard length constraint + reshape projection, explicitly
    rejecting an axial penalty spring), so a Hookean axial stiffness is a modelling GAP with no sourced value in
    ``ff.units``.  :func:`build_sf_mechanics_topology` therefore REQUIRES it (no convenient default); the Lead
    supplies a sourced/derived value or surfaces to PI, per the no-magic-number rule.

engine units: length µm, force pN, stiffness pN/µm.  This module imports ``warp`` (like ``surface_body``) so the
structural gate can validate the wiring with a recording launcher on a CUDA-free host; :meth:`bind_native` and
the device builders allocate real Warp arrays only on the gbook A5000.

Sanity Gate (self-tested in tests/ac/engine/test_sf_mechanics.py):
  * launch: :meth:`SFFilamentMechanics.accumulate` launches ``link_spring_kernel`` over the SF link arrays and
    ``cytosim_bending_kernel`` over the SF triple arrays, in that order, with the byte-identical cortex arg
    order and dims ``n_links``/``n_triples`` — verified by a recording launcher.
  * physics (NumPy reference of the exact kernel force law): a STRETCHED SF segment develops a non-zero axial
    tension along its axis; a BENT SF triple develops a non-zero restoring bending force — proving the built SF
    arrays are physically meaningful, not empty scaffolding.
  * ownership/disjoint: the topology is built from the DISJOINT ``sf_arc`` population (unique IDs, cortex-
    disjoint) and every link/triple index addresses an ``sf_arc`` node — no cortex node is ever referenced.
  * SOURCED/GAP: ``α`` uses ``KAPPA_ACTIN`` (sourced); ``k_axial_pn_per_um`` is REQUIRED (build raises on None).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.engine.sf_population import SFArcPopulation
from aleph.laws.hand_kmc import ALPHA_ACTININ
from aleph.laws.units import KAPPA_ACTIN, LP_ACTIN_UM

__all__ = [
    "SF_COMPONENT",
    "SFMechanicsTopology",
    "SFFilamentMechanics",
    "SFInternalArcJointConnector",
    "DORSAL_ARC_CROSSLINK",
    "KAPPA_ACTIN_PN_UM2",
    "ALPHA_ACTININ_K_PN_PER_UM",
    "SF_CONNECTOR_BINDING_STATUS",
    "build_sf_mechanics_topology",
    "build_sf_filament_mechanics",
    "build_sf_internal_arc_connector",
    "link_spring_force_reference",
    "bending_force_reference",
]

SF_COMPONENT = "sf_arc"
DORSAL_ARC_CROSSLINK = "dorsal_arc_crosslink"

#: SOURCED — actin bending modulus κ = k_B·T·ℓ_p (Gittes et al. 1993, KU-1.1); ℓ_p = 17 µm.  α = κ/seg³.
KAPPA_ACTIN_PN_UM2 = float(KAPPA_ACTIN)
#: SOURCED — α-actinin crosslink stiffness (Ferrer 2008 PNAS AFM, 455 pN/nm; PI-approved 2026-06-30).
#: Imported from the single source of truth (``ff.hand_kmc.ALPHA_ACTININ.link_k``) — never re-literal'd here,
#: so a future revision of the sourced value cannot silently drift out of sync (verify:binding hygiene flag).
ALPHA_ACTININ_K_PN_PER_UM = float(ALPHA_ACTININ.link_k)

_EPS = 1.0e-12

#: Honest per-connector binding status for the four ``sf_arc`` graph edges (see module docstring / report).
SF_CONNECTOR_BINDING_STATUS: dict[str, str] = {
    # KERNEL_BOUND here: internal SF↔SF joints live in the single SF pos array, so link_spring_kernel binds
    # directly against SF-owned nodes with the SOURCED α-actinin stiffness.
    "dorsal_arc_crosslink": "KERNEL_BOUND (link_spring_kernel over SF-owned joint pairs, α-actinin k)",
    # KERNEL_BOUND + native-validated 2026-07-25 (ac/engine/sf_motor_slice.py): a real head-resolved NMII
    # population is straddle-placed on the STRAIGHT sarcomeres (curved cap excluded by construction) and the
    # two-array adjoint crossbridge scatters into the SF-owned force array.  Magnitudes remain PI-GAPs (N1-N9).
    "nmii_sf_motor": "KERNEL_BOUND (segment-motor KMC + split crossbridge over a live SF port; magnitudes GAP)",
    # SEAMED: inter-component; needs a live cortex bind-target PORT (two never-merged arrays) + a two-array
    # adjoint Hookean scatter (Newton's 3rd), which the single-array link_spring_kernel does not provide.  The
    # internal joint above proves the launch pattern; wire this once the composed cortex+SF slice supplies the
    # cortex port (mirror nmii_cortex_motor's split-ownership scatter).
    "sf_cortex_transient": "SEAMED (needs live cortex port + two-array adjoint scatter)",
    # SEAMED: needs the intermediate_filament component; the edge is owned by the IF rig, SF is the endpoint.
    "if_sf_plectin": "SEAMED (needs intermediate_filament component; owned by the IF rig)",
    # SEAMED: needs the microtubule component; the edge is owned by the MT rig, SF is the endpoint.
    "mt_sf_spectraplakin": "SEAMED (needs microtubule component; owned by the MT rig)",
}


# ── host-side CUDA-array metadata validators (no device data read; mirror surface_body's contract) ────────
def _device_is_cuda(array: object) -> bool:
    return bool(getattr(getattr(array, "device", None), "is_cuda", False))


def _storage_key(array: object) -> tuple[object, ...]:
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


def _validate_pair_array(array: object, *, label: str, cols: int) -> int:
    """Validate a CUDA ``(rows, cols)`` int32 connectivity array; return the row count."""
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != wp.int32:
        raise TypeError(f"{label} must have dtype wp.int32")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != 2 or int(shape[1]) != int(cols):
        raise ValueError(f"{label} must be a two-dimensional (rows, {cols}) int32 array")
    return int(shape[0])


def _validate_scalar_f64(array: object, *, label: str, rows: int) -> None:
    """Validate a CUDA one-dimensional float64 parameter array of a fixed length."""
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != wp.float64:
        raise TypeError(f"{label} must have dtype wp.float64")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != 1:
        raise ValueError(f"{label} must be a one-dimensional float64 array")
    if int(shape[0]) != int(rows):
        raise ValueError(f"{label} length {int(shape[0])} must match its topology row count {int(rows)}")


# ── NumPy references for the EXACT kernel force laws (parity oracles for the CPU structural gate) ─────────
def link_spring_force_reference(
    pos: npt.NDArray[np.float64],
    links: npt.NDArray[np.int32],
    k: npt.NDArray[np.float64],
    r0: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Pure-NumPy reference for ``ff.network_warp.link_spring_kernel`` (Hookean axial link).

    Matches the kernel term-for-term: ``f = k·(L−r0)/L·d`` added to node ``i`` (toward ``j`` when stretched),
    negated on ``j``.  Used ONLY by the CPU structural gate to prove the built SF arrays carry force — the
    production path launches the Warp kernel on CUDA, never this reference.
    """
    force = np.zeros_like(pos)
    for t in range(links.shape[0]):
        i, j = int(links[t, 0]), int(links[t, 1])
        d = pos[j] - pos[i]
        length = float(np.linalg.norm(d))
        if length > _EPS:
            f = (float(k[t]) * (length - float(r0[t])) / length) * d
            force[i] += f
            force[j] -= f
    return force


def bending_force_reference(
    pos: npt.NDArray[np.float64],
    triples: npt.NDArray[np.int32],
    alpha: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Pure-NumPy reference for ``ff.forces_warp.cytosim_bending_kernel`` (NF2007 discrete bending).

    Matches the kernel: ``d = pos[a] − 2·pos[b] + pos[c]``, ``f = α·d``, applying the triplet ``{−f, +2f, −f}``.
    """
    force = np.zeros_like(pos)
    for t in range(triples.shape[0]):
        a, b, c = int(triples[t, 0]), int(triples[t, 1]), int(triples[t, 2])
        d = pos[a] - 2.0 * pos[b] + pos[c]
        f = float(alpha[t]) * d
        force[a] -= f
        force[b] += 2.0 * f
        force[c] -= f
    return force


# ── host topology builder — SF-disjoint links / triples / α, from the population ─────────────────────────
@dataclass(slots=True)
class SFMechanicsTopology:
    """SF-disjoint passive rod-cable topology + parameters (host NumPy), ready to upload to the CUDA lane.

    All node indices are GLOBAL into :attr:`SFArcPopulation.pos` and address ``sf_arc`` nodes ONLY (the
    population is disjoint from the cortex).  ``links`` are per-fiber consecutive axial-segment pairs;
    ``bend_triples`` are per-fiber consecutive triples; ``bend_alpha = κ/seg³`` (SOURCED κ = KAPPA_ACTIN).
    ``arc_joints`` are the internal dorsal-free↔transverse-arc crosslink pairs (fed to the internal connector,
    NOT the component mechanics, to avoid a double launch).
    """

    n_nodes: int
    links: npt.NDArray[np.int32]          # (L, 2) axial-segment node pairs
    link_k: npt.NDArray[np.float64]       # (L,)  axial stiffness [pN/µm] (GAP — required, no default)
    link_r0: npt.NDArray[np.float64]      # (L,)  rest length = build segment length [µm] (unstretched at rest)
    bend_triples: npt.NDArray[np.int32]   # (T, 3) consecutive node triples
    bend_alpha: npt.NDArray[np.float64]   # (T,)  α = κ/seg³ [pN/µm]  (SOURCED κ = KAPPA_ACTIN)
    arc_joints: npt.NDArray[np.int32]     # (J, 2) internal dorsal-free↔arc crosslink pairs (connector, not mech)
    arc_k: npt.NDArray[np.float64]        # (J,)  crosslink stiffness [pN/µm] (SOURCED α-actinin)
    arc_r0: npt.NDArray[np.float64]       # (J,)  crosslink rest = build joint distance [µm]
    k_axial_pn_per_um: float              # the supplied axial GAP value (for provenance)
    # ── per-SEGMENT identity/polarity metadata, emitted in the SAME order as ``links`` ────────────────────
    # Produced by the one loop that emits the links, so a consumer (the ``nmii_sf_motor`` bind-target port)
    # never has to re-walk the population and risk drifting out of link order.  These describe the AXIAL
    # segments only — the α-actinin ``arc_joints`` are crosslinks, not filament segments, and are excluded.
    seg_filament_id: npt.NDArray[np.int64] = field(default_factory=lambda: np.zeros(0, np.int64))
    seg_polarity: npt.NDArray[np.int32] = field(default_factory=lambda: np.zeros(0, np.int32))
    seg_s0: npt.NDArray[np.float64] = field(default_factory=lambda: np.zeros(0, np.float64))
    seg_s1: npt.NDArray[np.float64] = field(default_factory=lambda: np.zeros(0, np.float64))

    @property
    def n_links(self) -> int:
        return int(self.links.shape[0])

    @property
    def n_triples(self) -> int:
        return int(self.bend_triples.shape[0])

    @property
    def n_arc_joints(self) -> int:
        return int(self.arc_joints.shape[0])


def build_sf_mechanics_topology(
    population: SFArcPopulation,
    *,
    k_axial_pn_per_um: float,
    crosslink_k_pn_per_um: float = ALPHA_ACTININ_K_PN_PER_UM,
    end_correction: bool = True,
) -> SFMechanicsTopology:
    """Build the SF-disjoint passive rod-cable topology from a :class:`SFArcPopulation`.

    Walks each bundle's fibers (``fiber_offsets`` within the bundle, offset by ``node_base`` into the flat
    global ``pos``) and emits consecutive axial-segment links + consecutive bending triples.  ``α = κ/seg³``
    uses the SOURCED ``KAPPA_ACTIN`` with the Cytosim ``p/(p−1)`` end correction (``forces_warp`` convention).
    Link/joint rest lengths are the build-time distances, so a relaxed SF injects zero passive force (prestress
    is the separate motor connector, never a lumped ``k_SF`` here).

    Args:
        population: the built DISJOINT ``sf_arc`` population (host NumPy).
        k_axial_pn_per_um: **REQUIRED, no default** — the actin axial backbone Hookean stiffness [pN/µm].  This
            is a modelling GAP (NF2007 treats the backbone as inextensible; no sourced ``EA_actin`` in
            ``ff.units``).  Supply a sourced/derived value or surface to PI; ``None``/non-positive raises.
        crosslink_k_pn_per_um: internal dorsal↔arc crosslink stiffness [pN/µm]; defaults to the SOURCED
            α-actinin value (Ferrer 2008).
        end_correction: apply the Cytosim p/(p−1) bending end-correction (matches ``ff.forces_warp``).

    Returns:
        The :class:`SFMechanicsTopology` (upload with :func:`build_sf_filament_mechanics`).

    Raises:
        ValueError: if ``k_axial_pn_per_um`` is not a positive finite value (REQUIRED-PARAM discipline).
    """
    if k_axial_pn_per_um is None or not np.isfinite(k_axial_pn_per_um) or k_axial_pn_per_um <= 0.0:
        raise ValueError(
            "REQUIRED-PARAM 'k_axial_pn_per_um' must be a supplied positive-finite value "
            "(actin axial backbone stiffness is a modelling GAP — NF2007 treats the backbone as inextensible, "
            "no sourced EA_actin in ff.units; source it or surface to PI, never default)"
        )
    if crosslink_k_pn_per_um <= 0.0 or not np.isfinite(crosslink_k_pn_per_um):
        raise ValueError("crosslink_k_pn_per_um must be positive and finite")

    pos = np.ascontiguousarray(population.pos, dtype=np.float64)
    n_nodes = int(pos.shape[0])

    links: list[tuple[int, int]] = []
    triples: list[tuple[int, int, int]] = []
    alpha: list[float] = []
    seg_fid: list[int] = []
    seg_pol: list[int] = []
    seg_s0: list[float] = []
    seg_s1: list[float] = []
    for bundle in population.bundles:
        offsets = np.asarray(bundle.fiber_offsets, dtype=np.int64)
        base = int(bundle.node_base)
        n_fibers = offsets.size - 1
        for f in range(n_fibers):
            lo = base + int(offsets[f])
            hi = base + int(offsets[f + 1])
            nodes = list(range(lo, hi))
            if len(nodes) < 2:
                continue
            seg_lengths = [float(np.linalg.norm(pos[nodes[s + 1]] - pos[nodes[s]]))
                           for s in range(len(nodes) - 1)]
            # persistent identity + barbed polarity + material arclength, per axial segment, IN LINK ORDER.
            # ``polarity`` is the fiber's barbed-end flag (+1 = barbed at the LAST node, -1 = at the FIRST);
            # a segment inherits its fiber's flag unchanged because links are emitted node-index-increasing,
            # which is the convention ``refresh_segment_barbed_kernel`` reads (polarity < 0 flips b-a).
            fid = int(bundle.global_fiber_ids[f])
            pol = 1 if int(bundle.polarity[f]) >= 0 else -1
            s_cursor = 0.0
            for s in range(len(nodes) - 1):
                links.append((nodes[s], nodes[s + 1]))
                seg_fid.append(fid)
                seg_pol.append(pol)
                seg_s0.append(s_cursor)
                s_cursor += seg_lengths[s]
                seg_s1.append(s_cursor)
            # α = κ/seg³ (mean segment length of THIS fiber), Cytosim end-correction p/(p−1).
            p = len(nodes) - 1
            seg_mean = float(np.mean(seg_lengths)) if seg_lengths else 1.0
            g = (p / max(p - 1.0, 1.0)) if (end_correction and p >= 2) else 1.0
            a_fiber = (KAPPA_ACTIN_PN_UM2 / max(seg_mean, _EPS) ** 3) * g
            for k in range(1, len(nodes) - 1):
                triples.append((nodes[k - 1], nodes[k], nodes[k + 1]))
                alpha.append(a_fiber)

    links_arr = np.asarray(links, dtype=np.int32).reshape(-1, 2)
    triples_arr = np.asarray(triples, dtype=np.int32).reshape(-1, 3)
    link_r0 = np.array([float(np.linalg.norm(pos[int(j)] - pos[int(i)])) for i, j in links_arr],
                       dtype=np.float64) if links_arr.shape[0] else np.zeros(0, np.float64)
    link_k = np.full(links_arr.shape[0], float(k_axial_pn_per_um), dtype=np.float64)
    alpha_arr = np.asarray(alpha, dtype=np.float64).reshape(-1)

    joints = np.asarray(population.dorsal_arc_joints, dtype=np.int32).reshape(-1, 2)
    arc_r0 = np.array([float(np.linalg.norm(pos[int(j)] - pos[int(i)])) for i, j in joints],
                      dtype=np.float64) if joints.shape[0] else np.zeros(0, np.float64)
    arc_k = np.full(joints.shape[0], float(crosslink_k_pn_per_um), dtype=np.float64)

    return SFMechanicsTopology(
        n_nodes=n_nodes, links=links_arr, link_k=link_k, link_r0=link_r0,
        bend_triples=triples_arr, bend_alpha=alpha_arr,
        arc_joints=joints, arc_k=arc_k, arc_r0=arc_r0,
        k_axial_pn_per_um=float(k_axial_pn_per_um),
        seg_filament_id=np.asarray(seg_fid, np.int64).reshape(-1),
        seg_polarity=np.asarray(seg_pol, np.int32).reshape(-1),
        seg_s0=np.asarray(seg_s0, np.float64).reshape(-1),
        seg_s1=np.asarray(seg_s1, np.float64).reshape(-1),
    )


# ── the KERNEL_BOUND SF mechanics delegate (cortex-analog; launches the reused ff kernels itself) ─────────
@dataclass(frozen=True, slots=True)
class SFFilamentMechanics:
    """``sf_arc`` mechanics delegate that LAUNCHES the ``ff`` link + Cytosim-bending Warp kernels itself.

    The SF-side counterpart of :class:`aleph.engine.surface_body.CortexFilamentMechanics`.  Owns SF-LOCAL
    device topology/parameter arrays (``links_d``/``link_k_d``/``link_r0_d``/``bend_triples_d``/``bend_alpha_d``)
    and, on :meth:`accumulate`, launches ``link_spring_kernel`` then ``cytosim_bending_kernel`` against the
    SF-owned ``pos``/``force`` arrays with the byte-identical cortex arg order.  ``n_nodes`` is exposed so the
    SF state owner can assert the mechanics addresses exactly the owned node count (never a concatenated array).

    ``link_kernel``/``bending_kernel`` and ``launch`` are injected so a CUDA-free host gate substitutes a
    recorder; :meth:`bind_native` wires the real kernels for the gbook CUDA lane.
    """

    device: str
    n_nodes: int
    links_d: wp.array
    link_k_d: wp.array
    link_r0_d: wp.array
    bend_triples_d: wp.array
    bend_alpha_d: wp.array
    link_kernel: object
    bending_kernel: object
    launch: Callable[..., object] = wp.launch
    n_links: int = field(init=False)
    n_triples: int = field(init=False)

    def __post_init__(self) -> None:
        if int(self.n_nodes) <= 0:
            raise ValueError("SF mechanics must own at least one node")
        n_links = _validate_pair_array(self.links_d, label="sf_arc.links_d", cols=2)
        _validate_scalar_f64(self.link_k_d, label="sf_arc.link_k_d", rows=n_links)
        _validate_scalar_f64(self.link_r0_d, label="sf_arc.link_r0_d", rows=n_links)
        n_triples = _validate_pair_array(self.bend_triples_d, label="sf_arc.bend_triples_d", cols=3)
        _validate_scalar_f64(self.bend_alpha_d, label="sf_arc.bend_alpha_d", rows=n_triples)
        if n_links == 0 and n_triples == 0:
            raise ValueError("SF mechanics must bind at least one link or bending triple")
        object.__setattr__(self, "n_links", n_links)
        object.__setattr__(self, "n_triples", n_triples)

    @classmethod
    def bind_native(
        cls,
        *,
        device: str,
        n_nodes: int,
        links_d: wp.array,
        link_k_d: wp.array,
        link_r0_d: wp.array,
        bend_triples_d: wp.array,
        bend_alpha_d: wp.array,
        launch: Callable[..., object] = wp.launch,
    ) -> "SFFilamentMechanics":
        """Wire the production ``ff`` link/bending kernels (imported lazily to keep the seam light)."""
        from aleph.laws.forces_warp import cytosim_bending_kernel
        from aleph.laws.network_warp import link_spring_kernel

        return cls(
            device=device, n_nodes=int(n_nodes),
            links_d=links_d, link_k_d=link_k_d, link_r0_d=link_r0_d,
            bend_triples_d=bend_triples_d, bend_alpha_d=bend_alpha_d,
            link_kernel=link_spring_kernel, bending_kernel=cytosim_bending_kernel, launch=launch,
        )

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the SF axial-backbone and Cytosim-bending kernels onto the SF force array.

        Both kernels ``wp.atomic_add`` into ``force``; the caller owns force zeroing.  Argument order is
        byte-identical to :meth:`CortexFilamentMechanics.accumulate` / the ``ff`` launch sites so the reused
        kernels get their published signatures, never a re-derived variant.
        """
        if self.n_links:
            self.launch(
                self.link_kernel,
                dim=self.n_links,
                inputs=[pos, self.links_d, self.link_k_d, self.link_r0_d, force],
                device=self.device,
            )
        if self.n_triples:
            self.launch(
                self.bending_kernel,
                dim=self.n_triples,
                inputs=[pos, self.bend_triples_d, self.bend_alpha_d, force],
                device=self.device,
            )


# ── the internal dorsal_arc_crosslink connector — genuinely KERNEL_BOUND (single SF array, link_spring) ──
@dataclass(frozen=True, slots=True)
class SFInternalArcJointConnector:
    """``dorsal_arc_crosslink`` internal connector — LAUNCHES ``link_spring_kernel`` over SF-owned joint pairs.

    The dorsal-free↔transverse-arc joints are INTERNAL to ``sf_arc`` (both endpoints are SF-owned nodes in the
    single SF ``pos`` array), so the Hookean joint binds directly with ``link_spring_kernel`` — the same reused
    kernel the component mechanics uses, launched over the ``(J, 2)`` joint table with the SOURCED α-actinin
    stiffness.  This is a REAL kernel binding, not a spy; the inter-component ``sf_cortex_transient`` edge stays
    SEAMED because it needs a two-array adjoint scatter into a live cortex port (see
    :data:`SF_CONNECTOR_BINDING_STATUS`).
    """

    name: str
    component_a: str
    component_b: str
    device: str
    joints_d: wp.array
    joint_k_d: wp.array
    joint_r0_d: wp.array
    link_kernel: object
    launch: Callable[..., object] = wp.launch
    n_joints: int = field(init=False)

    def __post_init__(self) -> None:
        if self.name != DORSAL_ARC_CROSSLINK:
            raise ValueError(f"internal SF connector must be {DORSAL_ARC_CROSSLINK!r}")
        if self.component_a != SF_COMPONENT or self.component_b != SF_COMPONENT:
            raise ValueError("dorsal_arc_crosslink is INTERNAL — both endpoints must be sf_arc")
        n_joints = _validate_pair_array(self.joints_d, label="sf_arc.arc_joints_d", cols=2)
        if n_joints == 0:
            raise ValueError("internal arc connector needs at least one joint")
        _validate_scalar_f64(self.joint_k_d, label="sf_arc.arc_k_d", rows=n_joints)
        _validate_scalar_f64(self.joint_r0_d, label="sf_arc.arc_r0_d", rows=n_joints)
        object.__setattr__(self, "n_joints", n_joints)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the Hookean joint force over the SF-owned dorsal↔arc pairs into the SF force array."""
        self.launch(
            self.link_kernel,
            dim=self.n_joints,
            inputs=[pos, self.joints_d, self.joint_k_d, self.joint_r0_d, force],
            device=self.device,
        )

    # ── transaction API (ADDED 2026-08-09 by the connector runtime census) ───────────────────────────
    #
    # `CellActor.assert_fully_bound` REQUIRES all four of these of EVERY connector, and this class had
    # none — so `dorsal_arc_crosslink` was counted as having a runtime while a `require_complete=True`
    # world would have rejected it with a TypeError.  The census found it; nothing else could have,
    # because the slice drivers that exercise this connector all pass `require_complete=False`.
    #
    # They are structural no-ops for a REASON, not for convenience: this joint is INTERNAL to `sf_arc`,
    # both its endpoints are nodes in the SF owner's single array, and it holds no reversible state of
    # its own.  The SF owner's transaction snapshots and reject-restores those nodes.  Owning a second
    # copy here would be two authorities over one array — the co-location trap in reverse.

    def snapshot_candidate(self) -> None:
        """No reversible connector state: the SF owner snapshots the nodes both endpoints live in."""

    def rollback(self, accepted: wp.array) -> None:
        """No-op: the SF owner's own transaction reject-restores its node positions."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """No-op: ``dorsal_arc_crosslink`` declares no kinetics, so it advances no irreversible state."""

    def accumulate_ledger(self, ledger: object) -> None:
        """No-op: the joint force is already in the SF force array the owner reduces into the ledger.

        Contributing it again here would double-count the internal load — the same reason
        ``FocalAdhesionStateOwner`` deliberately has no ledger hook.
        """


# ── CUDA-lane builders (allocate device memory; called on the gbook A5000, not on the dev Mac) ───────────
def _upload(array: npt.NDArray, dtype: object, device: str | None) -> wp.array:
    return wp.array(np.ascontiguousarray(array), dtype=dtype, device=device)


def build_sf_filament_mechanics(
    topology: SFMechanicsTopology,
    *,
    device: str | None = None,
    launch: Callable[..., object] = wp.launch,
) -> SFFilamentMechanics:
    """Upload an :class:`SFMechanicsTopology` and bind the real ``ff`` kernels (CUDA lane).

    The link/triple/parameter arrays become SF-local Warp device arrays; the returned
    :class:`SFFilamentMechanics` LAUNCHES ``link_spring_kernel`` + ``cytosim_bending_kernel`` over them.
    """
    links_d = _upload(topology.links.reshape(-1, 2), wp.int32, device)
    triples_d = _upload(topology.bend_triples.reshape(-1, 3), wp.int32, device)
    return SFFilamentMechanics.bind_native(
        device=str(device), n_nodes=topology.n_nodes,
        links_d=links_d,
        link_k_d=_upload(topology.link_k, wp.float64, device),
        link_r0_d=_upload(topology.link_r0, wp.float64, device),
        bend_triples_d=triples_d,
        bend_alpha_d=_upload(topology.bend_alpha, wp.float64, device),
        launch=launch,
    )


def build_sf_internal_arc_connector(
    topology: SFMechanicsTopology,
    *,
    device: str | None = None,
    launch: Callable[..., object] = wp.launch,
) -> SFInternalArcJointConnector:
    """Upload the internal dorsal↔arc joints and bind ``link_spring_kernel`` (CUDA lane, genuine binding)."""
    from aleph.laws.network_warp import link_spring_kernel

    return SFInternalArcJointConnector(
        name=DORSAL_ARC_CROSSLINK, component_a=SF_COMPONENT, component_b=SF_COMPONENT, device=str(device),
        joints_d=_upload(topology.arc_joints.reshape(-1, 2), wp.int32, device),
        joint_k_d=_upload(topology.arc_k, wp.float64, device),
        joint_r0_d=_upload(topology.arc_r0, wp.float64, device),
        link_kernel=link_spring_kernel, launch=launch,
    )


# provenance for the report / native gate header.
_SOURCED_NOTE = (
    f"SOURCED: bending κ = KAPPA_ACTIN = {KAPPA_ACTIN_PN_UM2:.4g} pN·µm² "
    f"(k_B·T·ℓ_p, ℓ_p={LP_ACTIN_UM} µm, Gittes 1993, KU-1.1) → α=κ/seg³; "
    f"internal crosslink k = {ALPHA_ACTININ_K_PN_PER_UM:.3g} pN/µm (α-actinin, Ferrer 2008). "
    "GAP (required, no default): axial backbone k_axial_pn_per_um (NF2007 inextensible; no sourced EA_actin)."
)
