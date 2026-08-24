r"""``ecm``-owned collagen constitutive force pass — the ECM component LAUNCHES its own Warp kernels.

WHAT THIS CLOSES.  :mod:`aleph.components.ecm.device_schema` owns the CUDA-resident collagen Mikado
topology SoA (:class:`~aleph.components.ecm.device_schema.ECMTopologyState`) and
:mod:`aleph.components.ecm.mikado_topology` builds it on-GPU, but the schema deliberately carries **no
spring law, modulus, or solver** (``device_schema`` lines 1-6; constitutive-hold gate
``mikado_topology`` Sanity Gate).  Until now the ``ecm`` component's ``mechanics`` slot
(:class:`aleph.engine.ecm_world.ECMStateOwner.mechanics`, a bare
:class:`~aleph.engine.runtime.MechanicsContributor`) had no concrete delegate that launched a
real collagen constitutive law over the SoA — the only injectable was the legacy compartment adapter
(:class:`~aleph.engine.ecm_world.LegacyLocalMikadoMechanicsAdapter`).  This module is the
concrete delegate that makes ``ecm`` genuinely **KERNEL_BOUND** over its *own* device SoA: it OWNS the
derived per-element stiffness/α parameter arrays and LAUNCHES the already-audited passive-fiber kernels
over the SoA topology, exactly the way the cortex and ``sf_arc`` do through
:class:`aleph.engine.surface_body.CortexFilamentMechanics` /
:class:`aleph.engine.sf_mechanics.SFFilamentMechanics`.

THE TEMPLATE IT MIRRORS (sf_arc → ecm, byte-identical launch signatures).
:class:`ECMConstitutiveForce` is the ECM counterpart of :class:`SFFilamentMechanics`.  It launches:

  * ``ff.network_warp.link_spring_kernel`` — Hookean axial backbone tension ``k·(L−r0)/L·d`` on the
    collagen segment pairs ``segments_d``  (arg order ``[pos, segments_d, k_seg_d, seg_rest_d, force]``,
    ``r0 = segment_rest_length_d`` = the build segment length); and
  * ``ff.forces_warp.cytosim_bending_kernel`` — NF2007 discrete bending
    ``F = α·(m_{i-1}−2 m_i + m_{i+1})`` on the collagen consecutive-node triples ``bend_triples_d``
    (arg order ``[pos, bend_triples_d, alpha_d, force]``), with ``α = κ/seg³``.

against the *ECM-owned* ``position_d`` / ``force_d`` supplied by :class:`ECMStateOwner.accumulate`.
Both kernels ``wp.atomic_add`` into ``force``; the caller owns force zeroing.  It binds the SoA
``segments_d`` / ``bend_triples_d`` / ``segment_rest_length_d`` **directly** (index-identical layout,
``device_schema`` lines 79/96/92) — a thin device binding, not a host copy, so there is no per-step
GPU→CPU roundtrip (I0-A GPU-only).

DERIVED PARAMETERS (the schema does not carry them).  ``k_seg = EA/seg_rest`` (per-segment axial
stiffness) and ``α = κ/seg³`` (per-triple bending) are computed **on device** from the SOURCED collagen
card scalars + the SoA rest lengths by :meth:`refresh_parameters` (two derivation kernels, no host
readback).  They must be recomputed on any accepted remesh (watch ``topology_epoch_d``); within a single
candidate step the topology is constant, so :meth:`accumulate` never refreshes (no commit inside
accumulate — accepted-step contract, ``runtime`` lines 51-52).

ACTIVE-FLAG GATING (fixed-capacity SoA, no host count).  The launches cover the full segment/bend
**capacity** and never size a launch from a host active count.  Inactive/dormant slots are gated to a
ZERO contribution *through the parameter arrays*: the derivation kernels write ``k_seg = 0`` /
``α = 0`` wherever ``segment_active_d`` / ``bend_active_d`` is 0, so an inactive element's
``atomic_add`` is exactly zero.  Dormant slots additionally hold ``(node_0, node_0)`` pairs
(zero-length ⇒ the ``L > 1e-12`` guard skips them).  This is the parameter-mask form of the SoA
``if *_active_d > 0`` idiom (``mikado_topology`` lines 525/592) that keeps the reused ``ff`` kernels
unchanged.

EMERGENT, NOT LUMPED (CLAUDE.md hard rule).  These are the PASSIVE stretch/bend laws.  The link rest
length ``r0`` is the build-time segment length (``segment_rest_length_d = s1 − s0``), so a **relaxed**
collagen network injects ZERO force — a **stretched** segment develops restoring tension and a **bent**
triple develops restoring bending force.  The macroscopic modulus is NOT a lumped material constant baked
in here; it EMERGES from this microstructure (fiber density → connectivity, κ, EA).  Collagen crosslinks
are NOT part of this pass: they are intentionally absent from ``ECMTopologyState`` (``device_schema``
lines 56-59) and are graph-owned by the ``ecm_crosslink`` connector, which runs its own
``link_spring_kernel`` over accepted-step bond arrays — never sourced from the topology SoA.

SOURCED vs GAP (report-not-tune; provenance in :data:`COLLAGEN_FORCE_PROVENANCE`).
  * SOURCED — bending modulus ``κ = k_B·T·L_p`` with collagen ``L_p = 17 µm`` (Licup 2015 PNAS,
    ``ecm_library`` line 127), so ``α = κ/seg³`` is sourced, not chosen.
  * GAP (required, no default) — the axial stretch modulus ``EA = E_fibril·π r_f²`` carries the
    single-fibril Young's modulus ``E_fibril_Pa = 1.1e6`` which has **no inline DOI** on the collagen
    card (``ecm_library`` line 114) and the fibril diameter ``0.10 µm`` (no inline DOI).  Both feed
    ``k_seg = EA/seg_rest`` and thus set the axial modulus.  :func:`build_collagen_constitutive_force`
    and :class:`ECMConstitutiveForce` therefore **require** a positive-finite ``EA`` taken from the card
    and never fall back to the ``ecm_library.ECMSpec.k_seg_pN_um`` convenience default (``5.0e4``); the
    GAP is flagged for PI, not silently sourced.

engine units: length µm, force pN, stiffness pN/µm, ``κ`` pN·µm², ``EA`` pN (``1 pN/µm² = 1 Pa``).
This module imports ``warp`` (like ``sf_mechanics`` / ``surface_body``) so the structural gate can
validate the wiring with a recording launcher on a CUDA-free host; :func:`build_collagen_constitutive_force`
allocates real Warp arrays only on the gbook A5000.

Sanity Gate (self-tested in tests/ac/ecm/test_ecm_forces.py):
  * launch: :meth:`ECMConstitutiveForce.accumulate` launches ``link_spring_kernel`` over ``segments_d``
    then ``cytosim_bending_kernel`` over ``bend_triples_d``, dims = segment/bend capacity, with the
    byte-identical cortex/SF arg order — verified by a recording launcher (no host kernel execution).
  * physics (NumPy reference of the EXACT kernel force law + the exact derivation): a RELAXED collagen
    network injects ≈0 force; a STRETCHED segment develops correct-sign restoring tension; a BENT triple
    develops a restoring bending force; the force scatters onto the correct SoA nodes.
  * gating: an inactive segment/bend derives ``k_seg = 0`` / ``α = 0`` ⇒ zero contribution.
  * SOURCED/GAP: ``α`` uses the SOURCED collagen ``κ``; ``EA`` is REQUIRED (build raises on ≤0), never
    the ``ecm_library`` 5e4 fallback.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.components.ecm.device_schema import ECMTopologyState
from aleph.laws.ecm_library import ECMSpec, get_spec

__all__ = [
    "ECM_COMPONENT",
    "COLLAGEN_MATERIAL_KEY",
    "COLLAGEN_FORCE_PROVENANCE",
    "ECMConstitutiveForce",
    "collagen_segment_stiffness_np",
    "collagen_bend_alpha_np",
    "link_spring_force_reference",
    "bending_force_reference",
    "build_collagen_constitutive_force",
    "collagen_force_provenance",
]

ECM_COMPONENT = "ecm"
#: Default collagen card key in :data:`aleph.laws.ecm_library.REGISTRY` (MCF7 first baseline ECM).
COLLAGEN_MATERIAL_KEY = "collagen_I"

_EPS = 1.0e-12

#: Honest per-field provenance of the collagen constitutive inputs this pass consumes from the card.
#: SOURCED = literature-anchored; GAP = present on the card but with no inline DOI (surface to PI before
#: closing a modulus-production gate).  This pass consumes ``κ`` (SOURCED) and ``EA`` (GAP via E_fibril).
COLLAGEN_FORCE_PROVENANCE: dict[str, str] = {
    "kappa_pN_um2": "SOURCED (κ = k_B·T·L_p, collagen L_p = 17 µm — Licup 2015 PNAS, ecm_library:127)",
    "E_fibril_Pa": "GAP (single-fibril Young's modulus 1.1e6, no inline DOI — ecm_library:114; feeds EA→k_seg)",
    "fiber_diameter_um": "GAP (fibril diameter 0.10 µm, no inline DOI — ecm_library:114; sets r_f→EA)",
    "EA_pN": "GAP (EA = E_fibril·π r_f²; carries the two GAPs above — the axial modulus setter)",
    "seg_rest_um": "SOURCED-BY-CONSTRUCTION (per-segment build length s1−s0; r0 ⇒ relaxed network = 0 force)",
    "k_seg_pN_um": "DERIVED (k_seg = EA/seg_rest, per-segment; REQUIRED EA, never the ecm_library 5e4 fallback)",
    "k_xl_pN_um": "NOT BOUND HERE (crosslinks are graph-owned; absent from ECMTopologyState by design)",
}


def collagen_force_provenance(spec: ECMSpec) -> dict[str, object]:
    """Return the SOURCED/GAP provenance of the collagen inputs this pass consumes from ``spec``.

    Reports the concrete card values plus the :data:`COLLAGEN_FORCE_PROVENANCE` verdicts so a caller
    (or the native gate) can surface the GAPs to PI before a modulus-production gate.
    """
    return {
        "material": spec.key,
        "kappa_pN_um2": float(spec.kappa_pN_um2),
        "E_fibril_Pa": float(spec.E_fibril_Pa),
        "fiber_diameter_um": float(spec.fiber_diameter_um),
        "EA_pN": float(spec.EA_pN),
        "verdicts": dict(COLLAGEN_FORCE_PROVENANCE),
    }


# ── host-side CUDA-array metadata validators (no device data read; mirror sf_mechanics's contract) ────────
def _device_is_cuda(array: object) -> bool:
    return bool(getattr(getattr(array, "device", None), "is_cuda", False))


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


def _validate_int_vector(array: object, *, label: str, rows: int) -> None:
    """Validate a CUDA one-dimensional int32 array of a fixed length."""
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != wp.int32:
        raise TypeError(f"{label} must have dtype wp.int32")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != 1 or int(shape[0]) != int(rows):
        raise ValueError(f"{label} must be a one-dimensional int32 array of length {int(rows)}")


def _validate_scalar_f64(array: object, *, label: str, rows: int) -> None:
    """Validate a CUDA one-dimensional float64 parameter array of a fixed length."""
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != wp.float64:
        raise TypeError(f"{label} must have dtype wp.float64")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != 1 or int(shape[0]) != int(rows):
        raise ValueError(f"{label} length must equal its topology row count {int(rows)}")


# ── on-device parameter derivation kernels (fill the SoA-absent k_seg / α from card scalars) ─────────────
@wp.kernel
def _collagen_segment_stiffness_kernel(
    seg_rest: wp.array(dtype=wp.float64),       # (S,) per-segment build rest length [µm]
    seg_active: wp.array(dtype=wp.int32),       # (S,) active-flag gate
    ea_pn: wp.float64,                          # axial stretch modulus EA = E_fibril·π r_f² [pN]
    k_seg: wp.array(dtype=wp.float64),          # (S,) out: k = EA/seg_rest [pN/µm], 0 where inactive
):
    """Per-segment axial spring stiffness ``k = EA/seg_rest``; 0 for inactive/zero-length slots (gating)."""
    s = wp.tid()
    if seg_active[s] > wp.int32(0):
        r0 = seg_rest[s]
        if r0 > wp.float64(_EPS):
            k_seg[s] = ea_pn / r0
        else:
            k_seg[s] = wp.float64(0.0)
    else:
        k_seg[s] = wp.float64(0.0)


@wp.kernel
def _collagen_bend_alpha_kernel(
    seg_rest: wp.array(dtype=wp.float64),       # (S,) per-segment build rest length [µm]
    bend_active: wp.array(dtype=wp.int32),      # (B,) active-flag gate
    bend_left: wp.array(dtype=wp.int32),        # (B,) left segment index of the triple
    bend_right: wp.array(dtype=wp.int32),       # (B,) right segment index of the triple
    kappa: wp.float64,                          # bending modulus κ = k_B·T·L_p [pN·µm²] (SOURCED)
    alpha: wp.array(dtype=wp.float64),          # (B,) out: α = κ/seg³ [pN/µm], 0 where inactive
):
    """Per-triple bending ``α = κ/seg³`` with ``seg`` the local mean of the triple's two segments; 0 inactive."""
    b = wp.tid()
    if bend_active[b] > wp.int32(0):
        ls = bend_left[b]
        rs = bend_right[b]
        seg_local = wp.float64(0.5) * (seg_rest[ls] + seg_rest[rs])
        if seg_local > wp.float64(_EPS):
            alpha[b] = kappa / (seg_local * seg_local * seg_local)
        else:
            alpha[b] = wp.float64(0.0)
    else:
        alpha[b] = wp.float64(0.0)


# ── NumPy references for the EXACT kernel force laws + the exact derivations (CPU structural gate) ────────
def collagen_segment_stiffness_np(
    seg_rest: npt.NDArray[np.float64],
    seg_active: npt.NDArray[np.int32],
    ea_pn: float,
) -> npt.NDArray[np.float64]:
    """Pure-NumPy reference for :func:`_collagen_segment_stiffness_kernel` (``k = EA/seg_rest``, 0 inactive)."""
    k = np.zeros(seg_rest.shape[0], dtype=np.float64)
    active = (np.asarray(seg_active) > 0) & (np.asarray(seg_rest) > _EPS)
    k[active] = float(ea_pn) / np.asarray(seg_rest)[active]
    return k


def collagen_bend_alpha_np(
    seg_rest: npt.NDArray[np.float64],
    bend_active: npt.NDArray[np.int32],
    bend_left: npt.NDArray[np.int32],
    bend_right: npt.NDArray[np.int32],
    kappa: float,
) -> npt.NDArray[np.float64]:
    """Pure-NumPy reference for :func:`_collagen_bend_alpha_kernel` (``α = κ/seg³``, 0 inactive)."""
    alpha = np.zeros(bend_active.shape[0], dtype=np.float64)
    for b in range(bend_active.shape[0]):
        if int(bend_active[b]) <= 0:
            continue
        seg_local = 0.5 * (float(seg_rest[int(bend_left[b])]) + float(seg_rest[int(bend_right[b])]))
        if seg_local > _EPS:
            alpha[b] = float(kappa) / seg_local ** 3
    return alpha


def link_spring_force_reference(
    pos: npt.NDArray[np.float64],
    links: npt.NDArray[np.int32],
    k: npt.NDArray[np.float64],
    r0: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Pure-NumPy reference for ``ff.network_warp.link_spring_kernel`` (Hookean axial link).

    Term-for-term: ``f = k·(L−r0)/L·d`` added to node ``i`` (toward ``j`` when stretched), negated on ``j``.
    Used ONLY by the CPU structural gate; the production path launches the Warp kernel on CUDA.
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

    Term-for-term: ``d = pos[a] − 2·pos[b] + pos[c]``, ``f = α·d``, applying the triplet ``{−f, +2f, −f}``.
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


# ── the KERNEL_BOUND ECM constitutive delegate (sf-analog; launches the reused ff kernels itself) ────────
@dataclass(frozen=True, slots=True)
class ECMConstitutiveForce:
    """``ecm`` mechanics delegate that LAUNCHES the ``ff`` link + Cytosim-bending Warp kernels over the SoA.

    The ECM-side counterpart of :class:`aleph.engine.sf_mechanics.SFFilamentMechanics`.  It binds the
    :class:`~aleph.components.ecm.device_schema.ECMTopologyState` connectivity arrays **directly** (``segments_d``,
    ``bend_triples_d``, ``segment_rest_length_d`` as ``r0``) and OWNS the two derived parameter arrays
    (``k_seg_d`` = EA/seg_rest, ``alpha_d`` = κ/seg³) that the schema deliberately does not carry.  On
    :meth:`accumulate` it launches ``link_spring_kernel`` then ``cytosim_bending_kernel`` against the
    ECM-owned ``pos`` / ``force`` supplied by :class:`ECMStateOwner`, over the full segment/bend capacity,
    gated to zero on inactive slots through the derived parameters.

    ``n_nodes`` is exposed so :class:`ECMStateOwner` can cross-check the mechanics addresses exactly the owned
    node capacity (``ecm_world`` lines 669-674).  The kernels and ``launch`` are injected so a CUDA-free host
    gate substitutes a recorder; :func:`build_collagen_constitutive_force` wires the real kernels for gbook.
    """

    device: str
    n_nodes: int
    # SoA connectivity (bound directly — index-identical layout, not copied)
    segments_d: wp.array          # (S, 2) int32 node pairs  (device_schema:79)
    seg_rest_d: wp.array          # (S,)   float64 build rest length = r0  (device_schema:92)
    segment_active_d: wp.array    # (S,)   int32 active gate  (device_schema:81)
    bend_triples_d: wp.array      # (B, 3) int32 consecutive-node triples  (device_schema:96)
    bend_left_segment_d: wp.array   # (B,) int32 left segment of the triple  (device_schema:99)
    bend_right_segment_d: wp.array  # (B,) int32 right segment of the triple  (device_schema:100)
    bend_active_d: wp.array       # (B,)   int32 active gate  (device_schema:98)
    # owned derived parameters (filled on device by refresh_parameters)
    k_seg_d: wp.array             # (S,)   float64 k = EA/seg_rest [pN/µm]
    alpha_d: wp.array             # (B,)   float64 α = κ/seg³ [pN/µm]
    # SOURCED/GAP card scalars (host config metadata, passed as kernel scalars — provenance)
    ea_pn: float                  # EA = E_fibril·π r_f² [pN]  (GAP — required, no default)
    kappa_pn_um2: float           # κ = k_B·T·L_p [pN·µm²]  (SOURCED)
    # injected kernels / launcher
    link_kernel: object
    bending_kernel: object
    stiffness_kernel: object
    alpha_kernel: object
    launch: Callable[..., object] = wp.launch
    n_segments: int = field(init=False)
    n_bends: int = field(init=False)

    def __post_init__(self) -> None:
        if int(self.n_nodes) <= 0:
            raise ValueError("ECM constitutive force must own at least one collagen node")
        if not np.isfinite(self.ea_pn) or self.ea_pn <= 0.0:
            raise ValueError(
                "REQUIRED-PARAM 'ea_pn' must be a supplied positive-finite axial modulus EA = E_fibril·π r_f² "
                "[pN] from the collagen card (E_fibril is a GAP — no inline DOI, ecm_library:114); source it "
                "or surface to PI, never the ecm_library.ECMSpec.k_seg_pN_um 5e4 fallback"
            )
        if not np.isfinite(self.kappa_pn_um2) or self.kappa_pn_um2 <= 0.0:
            raise ValueError("REQUIRED-PARAM 'kappa_pn_um2' must be a positive-finite SOURCED κ = k_B·T·L_p")
        n_segments = _validate_pair_array(self.segments_d, label="ecm.segments_d", cols=2)
        _validate_scalar_f64(self.seg_rest_d, label="ecm.segment_rest_length_d", rows=n_segments)
        _validate_int_vector(self.segment_active_d, label="ecm.segment_active_d", rows=n_segments)
        _validate_scalar_f64(self.k_seg_d, label="ecm.k_seg_d", rows=n_segments)
        n_bends = _validate_pair_array(self.bend_triples_d, label="ecm.bend_triples_d", cols=3)
        _validate_int_vector(self.bend_left_segment_d, label="ecm.bend_left_segment_d", rows=n_bends)
        _validate_int_vector(self.bend_right_segment_d, label="ecm.bend_right_segment_d", rows=n_bends)
        _validate_int_vector(self.bend_active_d, label="ecm.bend_active_d", rows=n_bends)
        _validate_scalar_f64(self.alpha_d, label="ecm.alpha_d", rows=n_bends)
        if n_segments == 0 and n_bends == 0:
            raise ValueError("ECM constitutive force must bind at least one segment or bending triple")
        object.__setattr__(self, "n_segments", n_segments)
        object.__setattr__(self, "n_bends", n_bends)

    def refresh_parameters(self) -> None:
        """(Re)derive ``k_seg_d`` = EA/seg_rest and ``alpha_d`` = κ/seg³ on device from the SoA rest lengths.

        Fully GPU-resident (no host readback): reads ``segment_rest_length_d`` + the active flags + the
        SOURCED/GAP card scalars.  Call at build and after any accepted remesh that changes segment lengths
        or the active set (``topology_epoch_d`` increment).  MUST NOT be called inside a candidate's
        :meth:`accumulate` — topology is constant within one candidate step (accepted-step contract).
        """
        if self.n_segments:
            self.launch(
                self.stiffness_kernel,
                dim=self.n_segments,
                inputs=[self.seg_rest_d, self.segment_active_d, wp.float64(self.ea_pn), self.k_seg_d],
                device=self.device,
            )
        if self.n_bends:
            self.launch(
                self.alpha_kernel,
                dim=self.n_bends,
                inputs=[
                    self.seg_rest_d,
                    self.bend_active_d,
                    self.bend_left_segment_d,
                    self.bend_right_segment_d,
                    wp.float64(self.kappa_pn_um2),
                    self.alpha_d,
                ],
                device=self.device,
            )

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the collagen axial-backbone and Cytosim-bending kernels onto the ECM force array.

        Both kernels ``wp.atomic_add`` into ``force``; the caller owns force zeroing (``ecm_world`` docstring
        lines 943-945).  Argument order is byte-identical to :meth:`SFFilamentMechanics.accumulate` / the ``ff``
        launch sites so the reused kernels get their published signatures.  Launches over the full segment/bend
        capacity; inactive slots are gated to zero through ``k_seg_d`` / ``alpha_d`` (never a host active count).
        """
        if self.n_segments:
            self.launch(
                self.link_kernel,
                dim=self.n_segments,
                inputs=[pos, self.segments_d, self.k_seg_d, self.seg_rest_d, force],
                device=self.device,
            )
        if self.n_bends:
            self.launch(
                self.bending_kernel,
                dim=self.n_bends,
                inputs=[pos, self.bend_triples_d, self.alpha_d, force],
                device=self.device,
            )


# ── CUDA-lane builder (allocates device memory; called on the gbook A5000, not on the dev Mac) ───────────
def build_collagen_constitutive_force(
    topology: ECMTopologyState,
    spec: ECMSpec | str = COLLAGEN_MATERIAL_KEY,
    *,
    device: str | None = None,
    launch: Callable[..., object] = wp.launch,
) -> ECMConstitutiveForce:
    """Bind the collagen constitutive force pass to an on-GPU :class:`ECMTopologyState` (CUDA lane).

    Allocates the two ECM-owned derived parameter arrays (``k_seg_d`` / ``alpha_d``), binds the SoA
    connectivity directly, wires the real ``ff`` kernels, and derives the parameters once
    (:meth:`ECMConstitutiveForce.refresh_parameters`).  ``EA`` and ``κ`` come from the SOURCED collagen card
    (:data:`aleph.laws.ecm_library.REGISTRY`); ``EA`` is REQUIRED (raises on ≤0 — never the 5e4 fallback).

    Args:
        topology: the on-GPU collagen Mikado SoA (from :class:`MikadoTopologyBuilder`).
        spec: an :class:`ECMSpec` or a registry key/alias (default ``"collagen_I"``).
        device: the CUDA device the topology lives on (defaults to the topology's device).
        launch: injectable launcher (defaults to ``wp.launch``).

    Returns:
        The bound :class:`ECMConstitutiveForce` (parameters already derived).

    Raises:
        ValueError: if ``spec`` is not fibrillar, or its ``EA`` / ``κ`` is not positive-finite.
    """
    card = spec if isinstance(spec, ECMSpec) else get_spec(spec)
    if not card.is_fibrillar:
        raise ValueError(f"{card.key} is a continuum gel; the Mikado constitutive pass needs a fibrillar card")
    dev = topology.device if device is None else str(device)
    with wp.ScopedDevice(dev):
        k_seg_d = wp.zeros(topology.n_segment_capacity, dtype=wp.float64)
        alpha_d = wp.zeros(topology.n_bend_capacity, dtype=wp.float64)
    from aleph.laws.forces_warp import cytosim_bending_kernel
    from aleph.laws.network_warp import link_spring_kernel

    force_pass = ECMConstitutiveForce(
        device=str(dev),
        n_nodes=int(topology.n_node_capacity),
        segments_d=topology.segments_d,
        seg_rest_d=topology.segment_rest_length_d,
        segment_active_d=topology.segment_active_d,
        bend_triples_d=topology.bend_triples_d,
        bend_left_segment_d=topology.bend_left_segment_d,
        bend_right_segment_d=topology.bend_right_segment_d,
        bend_active_d=topology.bend_active_d,
        k_seg_d=k_seg_d,
        alpha_d=alpha_d,
        ea_pn=float(card.EA_pN),
        kappa_pn_um2=float(card.kappa_pN_um2),
        link_kernel=link_spring_kernel,
        bending_kernel=cytosim_bending_kernel,
        stiffness_kernel=_collagen_segment_stiffness_kernel,
        alpha_kernel=_collagen_bend_alpha_kernel,
        launch=launch,
    )
    force_pass.refresh_parameters()
    return force_pass
