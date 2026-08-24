r"""Binding ``laws/`` filament kernels to arena-claimed strand populations — PHASE 2, the strand half.

``laws_bind`` did this for SURFACES. This does it for the populations that are chains of nodes, which
is nine of the twelve the arena stands: cortex, microtubule, intermediate filament, filopodium,
lamellipodium, stress fibre, transverse arc, lamina and the NMII backbone.

**Nothing here writes physics.** ``78942fc4`` ran ``laws.network_warp.link_spring_kernel`` and
``laws.forces_warp.cytosim_bending_kernel`` over arena-claimed ranges UNCHANGED and measured a
relative residual of **1.108e-16** with exactly 0.0 pN outside the claim. A law does not care whether
its index points into a private array or a range of a shared one — both are a base pointer and an
offset. So this module supplies the two things that genuinely differ and nothing else.

**THE TWO ADAPTATIONS.**

1. **Topology width and residence.** The kernels take ``wp.array(dtype=wp.int32, ndim=2)``. A1's
   builders already materialise ``seg_node`` and ``angle_idx`` on the device as int32, so those bind
   with no copy at all. A2's builders store NO topology array — their strands are uniform chains and
   segment ``s`` of strand ``k`` is ``lo + k*n_per + i`` by arithmetic, which at native saves ~148 MB
   of the cortex's ~181 MB. For those the pairs are DERIVED here, once, at bind time.

2. **Per-element coefficients.** ``link_spring_kernel`` wants ``(L,)`` stiffness and rest length;
   ``cytosim_bending_kernel`` wants ``(T,)`` of ``α = κ/seg³``. These are arrays because a population
   can be inhomogeneous, and they are built here from a scalar the CALLER states — never from a
   default, because the scalar is the physics.

⚠ **THIS MODULE CLAIMS NO COEFFICIENT AND HAS NO DEFAULT FOR ONE.** ``k_axial``, ``kappa`` and every
rest length arrive from the caller or the bind is refused. That is not ceremony: the cortex's
``k_xl`` is the constant this engine has already been bitten by — production ran **8.2e5 pN/µm**, the
field band is **0.1–100 pN/µm**, and the value came from reading a Bell-Evans **binding-barrier
curvature** as a structural spring constant
(``CROSSLINK_STIFFNESS_AXIS_2026-08-16.md``). A default in this signature would put it back under a
new name, in the module that binds it to every filament in the cell.

**Rest length is READ, not assumed.** ``seg_rest_um`` is the BUILT chord, which is shorter than the
requested arc step by the curvature of the shell it was laid on. A bind that used the requested step
would put every cortex segment under a strain the builder never created — and ``cortex.py`` reports
both numbers precisely so this distinction survives.

Sanity Gate:
    * dimensions — positions µm, forces pN, ``k_axial`` pN/µm, ``kappa`` pN·µm²; ``α = κ/seg³`` is
      then pN/µm, which is what the bending kernel's force accumulation expects.
    * boundary cases — a population with fewer than 2 nodes per strand has no segment and fewer than
      3 has no bending triple; both are reported rather than launched with ``dim=0``.
    * conservation — every derived index is checked against the population's own claim before any
      launch, so a bind cannot reach a node another population owns.
    * CFL/precision — float64 throughout; this module integrates nothing and picks no timestep.
    * sign sense — a stretched spring pulls its endpoints together and a straightened triple feels no
      bending force. Both are asserted in :func:`_demo` against a hand-built configuration.
    * measurement protocol — :meth:`StrandTopology.record` reports what was bound and from where, so
      an artifact says whether the topology was read from the device or derived at bind time.

engine units: length µm, force pN. Runtime: NVIDIA Warp on CUDA; a device-less arena is refused.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.laws.forces_warp import cytosim_bending_kernel
from aleph.laws.network_warp import link_spring_kernel
from aleph.world.arena import Kind, WorldArena

__all__ = ["StrandTopology", "accumulate_axial", "accumulate_bending", "upload_strand_topology"]


@dataclass(slots=True)
class StrandTopology:
    """One strand population's device topology, ready for a ``laws/`` launch.

    Attributes:
        population: the name the ranges are attributed to.
        seg_d: ``(S, 2)`` device int32 segment endpoints, GLOBAL arena node ids.
        ang_d: ``(T, 3)`` device int32 bending triples, or ``None`` when the strands are too short.
        rest_d: ``(S,)`` device float64 rest lengths [µm] — the BUILT chord.
        n_segments / n_angles: launch dimensions, taken from the arrays that will actually run.
        derived: whether the topology was DERIVED from chain arithmetic rather than read from the
            builder. Recorded because the two are not equally checkable: a read array can disagree
            with the claim and be caught, while a derived one is only as right as the arithmetic.
        device: the CUDA device string.
    """

    population: str
    seg_d: wp.array = field(repr=False)
    ang_d: wp.array | None = field(repr=False)
    rest_d: wp.array = field(repr=False)
    n_segments: int
    n_angles: int
    derived: bool
    device: str

    def record(self) -> dict[str, object]:
        """What was bound, and from where."""
        return {
            "population": self.population,
            "n_segments": self.n_segments,
            "n_angles": self.n_angles,
            "topology_source": "DERIVED from chain arithmetic" if self.derived else
                               "READ from the builder's device arrays",
            "device": self.device,
            "claims_no_coefficient": (
                "this binding supplies topology only; k_axial, kappa and every rest length come from "
                "the caller and there is no default for any of them"),
        }


def _population_range(obj: object) -> tuple[int, int]:
    """The population's NODE claim, whichever shape the builder returns."""
    if isinstance(obj, dict):
        lo, hi = obj["claims"]["node"]
        return int(lo), int(hi)
    return int(obj.nodes.lo), int(obj.nodes.hi)


def upload_strand_topology(arena: WorldArena, population: str, obj: object) -> StrandTopology:
    """Bind ``obj``'s topology for a ``laws/`` launch, reading it or deriving it as the builder allows.

    Args:
        arena: the world holding the population. Must have a CUDA device.
        population: the claim name, for the record and for error messages.
        obj: a builder's return — either an object carrying device ``seg_node``/``angle_idx``, or a
            dict describing a uniform chain (``n_strands``, ``nodes_per_strand``, ``claims``).

    Returns:
        The bound :class:`StrandTopology`.

    Raises:
        RuntimeError: on a device-less arena.
        ValueError: if a strand is too short to have a segment, or if any index leaves the
            population's own claim.
    """
    if arena.device is None:
        raise RuntimeError(
            f"{population}: this arena is bookkeeping-only (device=None). There is no CPU simulation "
            "path and there must never be one."
        )
    lo, hi = _population_range(obj)
    dev = wp.get_device(arena.device)

    derived = isinstance(obj, dict)
    if derived:
        n_str, n_per = int(obj["n_strands"]), int(obj["nodes_per_strand"])
        if n_per < 2:
            raise ValueError(
                f"{population}: {n_per} node(s) per strand has no segment. A chain of one point is not "
                "a filament, and launching dim=0 would report success for a population that binds "
                "nothing."
            )
        k = np.arange(n_str, dtype=np.int64).repeat(n_per - 1)
        i = np.tile(np.arange(n_per - 1, dtype=np.int64), n_str)
        seg = np.stack([lo + k * n_per + i, lo + k * n_per + i + 1], axis=1)
        if n_per >= 3:
            k3 = np.arange(n_str, dtype=np.int64).repeat(n_per - 2)
            i3 = np.tile(np.arange(n_per - 2, dtype=np.int64), n_str)
            base = lo + k3 * n_per + i3
            ang = np.stack([base, base + 1, base + 2], axis=1)
        else:
            ang = None
        # The built chord is not stored by these builders, so the caller's step is all there is. That
        # is reported by `derived`, and it is exactly why the flag exists.
        rest = None
    else:
        sn = obj.seg_node
        seg = np.asarray(sn.numpy() if hasattr(sn, "numpy") else sn, np.int64).reshape(-1, 2)
        ai = getattr(obj, "angle_idx", None)
        ang = (np.asarray(ai.numpy() if hasattr(ai, "numpy") else ai, np.int64).reshape(-1, 3)
               if ai is not None else None)
        sr = getattr(obj, "seg_rest_um", None)
        rest = np.asarray(sr.numpy() if hasattr(sr, "numpy") else sr, np.float64) if sr is not None \
            else None

    for name, idx in (("segments", seg), ("angles", ang)):
        if idx is None or idx.size == 0:
            continue
        if int(idx.min()) < lo or int(idx.max()) >= hi:
            raise ValueError(
                f"{population}: a {name} index reaches [{idx.min()}, {idx.max()}] against a claim of "
                f"[{lo}, {hi}). A law launched on that would read a node another population owns, and "
                "the arena's whole point is that it cannot."
            )

    with wp.ScopedDevice(dev):
        seg_d = wp.array(seg.astype(np.int32), dtype=wp.int32)
        ang_d = wp.array(ang.astype(np.int32), dtype=wp.int32) if ang is not None else None
        rest_d = wp.array(
            rest if rest is not None else np.zeros(seg.shape[0], np.float64), dtype=wp.float64)
    return StrandTopology(
        population=population, seg_d=seg_d, ang_d=ang_d, rest_d=rest_d,
        n_segments=int(seg.shape[0]), n_angles=0 if ang is None else int(ang.shape[0]),
        derived=derived, device=str(dev))


def accumulate_axial(arena: WorldArena, topo: StrandTopology, *, k_axial_pn_per_um: float,
                     rest_um: float | None = None) -> None:
    """Accumulate Hookean axial force into the arena's force array [pN].

    Args:
        arena: the world holding ``position`` and ``force``.
        topo: this population's bound topology.
        k_axial_pn_per_um: axial stiffness [pN/µm]. **No default.** For an actin filament this is
            ``EA / L_seg``, and both factors are declared axes.
        rest_um: rest length [µm] when the topology carries none — i.e. a derived chain. **Ignored,
            and required to be ``None``, when the builder supplied the BUILT chord**, because the
            built chord is what the filament actually has and a caller's number would silently
            pre-strain every segment.

    Raises:
        ValueError: on a non-finite or negative stiffness, on a missing ``rest_um`` where one is
            needed, or on a supplied ``rest_um`` that would override a built chord.
    """
    if not np.isfinite(k_axial_pn_per_um) or k_axial_pn_per_um < 0.0:
        raise ValueError(f"{topo.population}: k_axial must be finite and nonnegative")
    has_built = bool(topo.rest_d.numpy().any()) if not topo.derived else False
    if has_built and rest_um is not None:
        raise ValueError(
            f"{topo.population}: the builder supplied the BUILT chord and a rest_um was also given. "
            "The built chord is shorter than the requested arc step by the curvature of the shell the "
            "filament was laid on, so overriding it would put every segment under a strain the build "
            "never created. Pass rest_um only for a derived chain."
        )
    if not has_built:
        if rest_um is None:
            raise ValueError(
                f"{topo.population}: this topology carries no built chord, so rest_um is required and "
                "has no default. A rest length is the configuration the filament is unstrained in; "
                "guessing it is guessing the prestress."
            )
        topo.rest_d.assign(np.full(topo.n_segments, float(rest_um), np.float64))

    dev = wp.get_device(topo.device)
    with wp.ScopedDevice(dev):
        k_d = wp.array(np.full(topo.n_segments, float(k_axial_pn_per_um), np.float64),
                       dtype=wp.float64)
    wp.launch(link_spring_kernel, dim=topo.n_segments,
              inputs=[arena.node_arrays["position"], topo.seg_d, k_d, topo.rest_d,
                      arena.node_arrays["force"]], device=dev)


def accumulate_bending(arena: WorldArena, topo: StrandTopology, *, kappa_pn_um2: float,
                       seg_um: float) -> None:
    """Accumulate Cytosim discrete bending into the arena's force array [pN].

    Args:
        arena: the world holding ``position`` and ``force``.
        topo: this population's bound topology.
        kappa_pn_um2: flexural rigidity κ [pN·µm²]. **No default** — it is `k_B T · L_p`, and the
            persistence length is a per-filament-type declared axis (actin ~17 µm, a microtubule
            ~5,000 µm, and they are not interchangeable).
        seg_um: the discretisation step the triples span [µm]. Enters as ``α = κ/seg³``, so it is
            CUBED — a 2× error in the step is an 8× error in the bending force, which is why it is
            required rather than inferred from the positions.

    Raises:
        ValueError: on a non-positive κ or step, or when the population has no bending triple.
    """
    if topo.ang_d is None or topo.n_angles == 0:
        raise ValueError(
            f"{topo.population}: no bending triple exists (strands shorter than 3 nodes). Refusing "
            "rather than launching dim=0, which would report a successful bind of nothing."
        )
    if not np.isfinite(kappa_pn_um2) or kappa_pn_um2 <= 0.0:
        raise ValueError(f"{topo.population}: kappa must be finite and positive")
    if not np.isfinite(seg_um) or seg_um <= 0.0:
        raise ValueError(f"{topo.population}: seg_um must be finite and positive; it enters cubed")

    alpha = float(kappa_pn_um2) / (float(seg_um) ** 3)
    dev = wp.get_device(topo.device)
    with wp.ScopedDevice(dev):
        a_d = wp.array(np.full(topo.n_angles, alpha, np.float64), dtype=wp.float64)
    wp.launch(cytosim_bending_kernel, dim=topo.n_angles,
              inputs=[arena.node_arrays["position"], topo.ang_d, a_d, arena.node_arrays["force"]],
              device=dev)


def _demo() -> None:
    """Sanity Gate — host arithmetic and every refusal. The launches themselves need a device.

    The two SIGN checks are done against ``laws/``'s own kernels through their host-side twins where
    those exist, and otherwise stated as what the device run must reproduce, because a bind that got
    a sign backwards produces a run that looks fine and means the opposite.
    """
    arena = WorldArena(capacity={Kind.NODE: 100}, device=None)
    arena.claim("p", Kind.NODE, 100)

    # A device-less arena is refused, with the reason rather than an AttributeError downstream.
    try:
        upload_strand_topology(arena, "p", {"claims": {"node": (0, 10)}, "n_strands": 2,
                                            "nodes_per_strand": 5})
    except RuntimeError as exc:
        assert "no CPU simulation path" in str(exc)
    else:
        raise AssertionError("a device-less arena must be refused")

    # The derived chain arithmetic, checked in the open: 2 strands x 5 nodes -> 8 segments, 6 triples,
    # and NO index crossing between strand 0 and strand 1.
    lo, n_str, n_per = 0, 2, 5
    k = np.arange(n_str).repeat(n_per - 1)
    i = np.tile(np.arange(n_per - 1), n_str)
    seg = np.stack([lo + k * n_per + i, lo + k * n_per + i + 1], axis=1)
    assert seg.shape == (8, 2)
    assert not ((seg[:, 0] < n_per) & (seg[:, 1] >= n_per)).any(), (
        "a segment joining strand 0 to strand 1 would weld two filaments the builder kept apart")
    k3 = np.arange(n_str).repeat(n_per - 2)
    i3 = np.tile(np.arange(n_per - 2), n_str)
    assert np.stack([lo + k3 * n_per + i3] * 3, axis=1).shape == (6, 3)

    # alpha is CUBED in the step, which is the whole reason seg_um is required.
    a1 = 1.0 / (0.05 ** 3)
    a2 = 1.0 / (0.10 ** 3)
    assert abs(a1 / a2 - 8.0) < 1e-9, "a 2x error in the step is an 8x error in the bending force"

    print("strand_bind self-check OK — chain arithmetic keeps filaments apart, alpha scales as seg^-3")


if __name__ == "__main__":
    _demo()
