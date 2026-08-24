r"""Binding ``laws/`` kernels to arena-addressed surfaces — PHASE 2, the surface half.

WHAT A BIND IS HERE, AND WHY IT IS NOT A PORT.  ``78942fc4`` ran two production kernels
(``laws.network_warp.link_spring_kernel``, ``laws.forces_warp.cytosim_bending_kernel``) over
arena-claimed ranges UNCHANGED and measured a relative residual of 1.108e-16 with exactly 0.0 pN
outside the claim.  Nothing about a law needs to change to read a range instead of a private array,
because both are just a base pointer and an index.  So this module writes no physics and re-derives
no formula.  What it does write is the small amount of ADAPTING that is genuinely needed, and each
piece of it exists for a reason recorded below rather than for symmetry.

THE THREE ADAPTATIONS, and that is all of them.

1. **Index width.**  ``laws/`` topology kernels take ``wp.array(dtype=wp.int32, ndim=2)``; the arena's
   builders carry ``int64`` host arrays because a global index into a ten-million-node arena is
   computed in Python.  ``10e6`` is comfortably inside int32, but the cast is where a future
   multi-cell arena would silently wrap, so it is checked rather than performed.

2. **Hinge tuple order.**  ``world.surface`` stores a dihedral as ``(a, b, c, d)`` with ``b-c`` the
   shared edge; ``laws.membrane_surface.helfrich_bending_kernel`` wants ``(i, j, k, l)`` with ``i-j``
   the shared edge.  That is the permutation ``(1, 2, 0, 3)`` and NOTHING ELSE — see the next
   paragraph, because the obvious objection to that sentence is correct-sounding and wrong.

   The objection: ``world.surface._hinges`` records ``(uses[0][1], u, v, uses[1][1])`` and DISCARDS
   which of the two faces carries the directed half-edge ``u -> v``.  The kernel's docstring is
   explicit that ``T1`` must contain ``i -> j``, so it looks as though the orientation has been lost
   and must be recovered from the face list.  **It has been lost, and it does not matter.**  Swapping
   ``i`` and ``j`` replaces ``T1=(i,j,k)`` by ``(j,i,k)`` and ``T2=(j,i,l)`` by ``(i,j,l)``, which
   reverses BOTH discrete normals at once, so ``n1 . n2`` — and therefore the energy
   ``kappa~(1 - n1.n2)`` and its gradient — are unchanged.  Measured on a perturbed subdivision-2
   icosphere: the assembled force under the two orderings differs by 5.7e-16 relative, i.e. round-off.
   :func:`_demo` pins the energy half of that with the module's own host twin, and the force half
   follows because the kernel force is ``-dE/dx`` (its own finite-difference test is the arbiter).

   This matters practically, not only for tidiness.  The alternative — calling
   ``laws.membrane_surface.build_membrane_hinges`` to rebuild the list in the kernel's own order — is
   a Python dict loop over every half-edge, which at the plan's subdivision-7 membrane is ~1e6 host
   dict operations for a permutation of four columns.

3. **Launch extent.**  Every launch is over ONE population's own claim count, never
   ``arena.n_live(Kind.FACE)``.  For a per-face force this is not an optimisation, it is correctness:
   :mod:`aleph.laws.volume` shows the assembled volume gradient is translation-invariant only when a
   surface's vertex links close, so a launch spanning two populations' face ranges — or half of
   one — is a different force, not a partial one.

WHAT THIS MODULE DELIBERATELY DOES NOT HOLD.  No tension, no bending rigidity, no pressure, no rest
area.  Those are parameters with provenance and they arrive as arguments.  A default here would be an
unsourced constant governing the physics under a new name, which is the failure the arena's
``radius_um``/``seg_um`` refusals already exist to prevent; this module refuses in the same way.

It also does not decide WHICH law a population carries.  That is the meaning question PHASE 2 answers
one population at a time, and the answer belongs beside the population, not in a registry here.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — positions [um]; ``sigma`` [pN/um]; ``kappa_tilde`` [pN.um]; ``delta_p`` [pN/um^2];
    accumulated force [pN].  Indices carry no unit.
  * boundary — an OPEN surface is refused for the volume-pressure law and accepted for the others,
    because only that law's invariance depends on the links closing; a surface with no hinges is
    refused for bending rather than launched at dimension zero, since "no bending term" and "a bending
    term that was never launched" are different facts.
  * conservation/invariant — every index uploaded is asserted to lie inside the population's own node
    claim before any launch; a topology that reaches outside its range is the defect the arena exists
    to make visible, and it is caught on the host where it is cheap.
  * CFL/precision — nothing is integrated here.  float64 positions and coefficients throughout; the
    int32 narrowing is range-checked.
  * sign sense — untouched: every kernel launched is the production kernel with its own sign, and this
    module reorders columns and passes pointers.
  * measurement protocol — topology is uploaded ONCE at build time.  No host readback, no allocation
    and no upload happens in any function called per iteration; the accumulate helpers launch and
    return.

engine units: length um, force pN.  Runtime: NVIDIA Warp on CUDA; the host-side adapters are
CPU-importable and :func:`_demo` needs no device.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.laws.membrane_surface import helfrich_bending_kernel, membrane_area_kernel
from aleph.laws.volume import is_closed_manifold, volume_pressure_kernel
from aleph.world.arena import Kind, WorldArena
from aleph.world.surface import Surface

__all__ = ["HINGE_PERMUTATION", "SurfaceTopology", "upload_surface_topology"]

#: ``world.surface`` stores ``(a, b, c, d)`` with ``b-c`` shared; ``laws`` wants ``(i, j, k, l)`` with
#: ``i-j`` shared.  A column permutation, and the module docstring records why no orientation recovery
#: is needed with it.
HINGE_PERMUTATION: tuple[int, int, int, int] = (1, 2, 0, 3)

_INT32_MAX = np.iinfo(np.int32).max


def _closedness(surface: object, face_rows: npt.NDArray[np.int64]) -> bool:
    """Whether the mesh is closed — DERIVED when the surface does not carry the answer.

    ⚠ **There are two `ClosedSurface` classes in `world/` and they are not the same object.**
    ``world/surface.py``'s holds host arrays and a ``closed: bool``; ``world/build/membrane.py``'s
    holds DEVICE arrays and has no such field. This function was written against the first and
    ``build_all`` returns the second, so it raised ``AttributeError`` on the first PHASE 4 run.

    Reading the flag when it exists and deriving it otherwise is the fix, and deriving is the better
    half: a closedness flag is a claim about the mesh, and ``laws.volume.is_closed_manifold`` checks
    the mesh itself — every edge shared by exactly two faces. A surface that says it is closed and is
    not would pass an attribute read and fail here, which is the direction the error should point.
    """
    declared = getattr(surface, "closed", None)
    derived = bool(is_closed_manifold(face_rows)) if face_rows.shape[0] else False
    if declared is not None and bool(declared) != derived:
        raise ValueError(
            f"{getattr(surface, 'population', '?')}: the surface declares closed={bool(declared)} but "
            f"its faces say {derived}. The volume-pressure law's translation invariance rests on this, "
            "so a disagreement is not a labelling question."
        )
    return derived


def _host_rows(idx: object, width: int) -> npt.NDArray[np.int64]:
    """A host ``(N, width)`` int64 view of a topology array, whether it lives on host or device.

    One conversion point on purpose. The alternative — reading the array directly for a shape check
    and converting separately for the upload — is what broke the first PHASE 4 run: the two paths
    disagreed the moment a builder moved its topology to the device.
    """
    host = idx.numpy() if hasattr(idx, "numpy") else idx
    return np.asarray(host, np.int64).reshape(-1, width)


def _as_int32_topology(name: str, idx: npt.NDArray[np.int64], surface: Surface) -> npt.NDArray[np.int32]:
    """Narrow a global-index topology array to int32, refusing rather than wrapping.

    Also asserts every index lies inside the surface's OWN node claim.  A face or hinge that reaches
    into a neighbouring population's range is exactly the failure one shared allocation makes possible
    and private arrays made impossible, so it is checked at the only point where it is still cheap.
    """
    if idx.size and (int(idx.min()) < surface.nodes.lo or int(idx.max()) >= surface.nodes.hi):
        raise ValueError(
            f"{surface.population}.{name} indexes node "
            f"[{int(idx.min())}, {int(idx.max())}] outside its own claim "
            f"[{surface.nodes.lo}, {surface.nodes.hi}). In one shared allocation this reads as a "
            "silent coupling to whichever population owns those ids."
        )
    if idx.size and int(idx.max()) > _INT32_MAX:
        raise ValueError(
            f"{surface.population}.{name} holds index {int(idx.max())}, past int32. The laws' topology "
            "kernels take int32; widen them before building an arena this large rather than wrapping."
        )
    return np.ascontiguousarray(idx, np.int32)


@dataclass(slots=True)
class SurfaceTopology:
    """One surface population's topology, resident on the device, in the layout ``laws/`` kernels read.

    Built once and never rebuilt: the arrays here are the only host->device transfer the surface laws
    need, and they are uploaded at construction so nothing in the inner loop allocates or copies.

    Attributes:
        population: the name the ranges belong to.
        faces_d: ``(F, 3)`` int32 GLOBAL node indices, outward-wound.
        hinges_d: ``(H, 4)`` int32 GLOBAL node indices in the kernel's ``(i, j, k, l)`` order, or
            ``None`` for an open surface with no dihedrals.
        n_faces / n_hinges: the launch dimensions — this population's counts, never the arena's live
            totals.
        closed: whether every edge is shared by exactly two faces.  The volume-pressure law needs it.
    """

    population: str
    faces_d: wp.array
    hinges_d: wp.array | None
    n_faces: int
    n_hinges: int
    closed: bool
    device: str = field(default="")


def upload_surface_topology(arena: WorldArena, surface: Surface) -> SurfaceTopology:
    """Upload ``surface``'s faces and hinges to ``arena``'s device in the kernels' layout.

    Args:
        arena: the world the surface was claimed from; supplies the device.
        surface: a built surface of that arena.

    Returns:
        The resident :class:`SurfaceTopology`.

    Raises:
        RuntimeError: if the arena holds no device allocation.
        ValueError: if a face or hinge indexes outside the surface's own node claim, or if the claim
            counts disagree with the topology arrays.
    """
    if arena.device is None:
        raise RuntimeError(
            "this arena is bookkeeping-only (device=None). Topology upload needs a CUDA device; there "
            "is no CPU simulation path and there must never be one."
        )
    # ⚠ Read the topology through NumPy ONCE, before validating against it. This function used to
    # validate on `surface.hinge_idx.reshape(-1, 4)` directly and then convert a few lines later —
    # which works while the surface holds host arrays and raises `array.reshape() takes 2 positional
    # arguments but 3 were given` the moment it holds device ones, because `wp.array.reshape` takes a
    # tuple. Found 2026-08-21 on the first PHASE 4 run, where the builder had already moved its
    # topology to the device. Converting first means the check and the use see the same object.
    face_rows = _host_rows(surface.face_idx, 3)
    hinge_rows = _host_rows(surface.hinge_idx, 4)

    if surface.faces.count != face_rows.shape[0]:
        raise ValueError(
            f"{surface.population}: FACE claim holds {surface.faces.count} but the topology has "
            f"{face_rows.shape[0]} faces — the launch dimension comes from the claim, so these "
            "disagreeing means one of them is not describing what will run."
        )
    if surface.hinges.count != hinge_rows.shape[0]:
        raise ValueError(
            f"{surface.population}: ANGLE4 claim holds {surface.hinges.count} but the topology has "
            f"{hinge_rows.shape[0]} hinges"
        )

    faces = _as_int32_topology("face_idx", face_rows, surface)
    hinges = (_as_int32_topology("hinge_idx", hinge_rows[:, list(HINGE_PERMUTATION)], surface)
              if hinge_rows.shape[0] else None)

    dev = wp.get_device(arena.device)
    with wp.ScopedDevice(dev):
        faces_d = wp.array(faces, dtype=wp.int32)
        hinges_d = wp.array(hinges, dtype=wp.int32) if hinges is not None else None
    return SurfaceTopology(
        population=surface.population, faces_d=faces_d, hinges_d=hinges_d,
        n_faces=int(faces.shape[0]), n_hinges=0 if hinges is None else int(hinges.shape[0]),
        closed=_closedness(surface, face_rows), device=str(dev),
    )


def accumulate_area_tension(arena: WorldArena, topo: SurfaceTopology, sigma: float) -> None:
    """Accumulate in-plane surface tension ``-sigma * dA/dx_i`` into the arena's force array [pN].

    Args:
        arena: the world holding ``position`` and ``force``.
        topo: this population's resident topology.
        sigma: the tension [pN/um].  **No default** — it is the reservoir-buffered membrane tension or
            an envelope's own, and which one is a per-population statement with provenance.
    """
    wp.launch(membrane_area_kernel, dim=topo.n_faces,
              inputs=[arena.node_arrays["position"], topo.faces_d, wp.float64(float(sigma)),
                      arena.node_arrays["force"]],
              device=topo.device)


def accumulate_helfrich(arena: WorldArena, topo: SurfaceTopology, kappa_tilde: float) -> None:
    """Accumulate discrete Helfrich dihedral bending into the arena's force array [pN].

    Args:
        arena: the world holding ``position`` and ``force``.
        topo: this population's resident topology.
        kappa_tilde: the CALIBRATED dihedral coupling [pN.um] from
            ``laws.membrane_surface.calibrate_kappa_tilde`` — not ``kappa_m``, and not the ideal-lattice
            ``2 kappa/sqrt(3)``, which under-represents an icosphere by ~2.9x.

    Raises:
        ValueError: if the surface has no hinges.  Launching at dimension zero would record "bending
            was applied" for a surface that has no dihedral to apply it over.
    """
    if topo.hinges_d is None or topo.n_hinges == 0:
        raise ValueError(
            f"{topo.population} has no hinges, so there is no dihedral for a bending law to act on. "
            "A zero-dimension launch would be indistinguishable in the record from a bending term "
            "that ran and found nothing."
        )
    wp.launch(helfrich_bending_kernel, dim=topo.n_hinges,
              inputs=[arena.node_arrays["position"], topo.hinges_d, wp.float64(float(kappa_tilde)),
                      arena.node_arrays["force"]],
              device=topo.device)


def accumulate_volume_pressure(arena: WorldArena, topo: SurfaceTopology, delta_p: float) -> None:
    """Accumulate ``delta_p * dV/dx_i`` over the surface's complete face range [pN].

    Args:
        arena: the world holding ``position`` and ``force``.
        topo: this population's resident topology.
        delta_p: inside-minus-outside pressure [pN/um^2 == Pa in engine units].  **No default** — the
            osmotic value, a nucleoplasm bulk term and a Biot pore pressure are three different
            physical quantities multiplying this one geometric gradient.

    Raises:
        ValueError: if the surface is not closed.  The assembled gradient is translation-invariant only
            when every vertex link closes, so on an open sheet this law returns an origin-dependent
            force that looks like a surface traction.
    """
    if not topo.closed:
        raise ValueError(
            f"{topo.population} is an OPEN surface. The volume gradient assembles to a "
            "translation-invariant field only when the vertex links close; on an open sheet the result "
            "depends on where the origin is, while still looking like a traction. See aleph.laws.volume."
        )
    wp.launch(volume_pressure_kernel, dim=topo.n_faces,
              inputs=[arena.node_arrays["position"], topo.faces_d, wp.float64(float(delta_p)),
                      arena.node_arrays["force"]],
              device=topo.device)


def _demo() -> None:
    """Self-check: the permutation, the range guard, and the closure condition.  No device needed."""
    from aleph.laws.membrane_surface import build_membrane_hinges, membrane_bending_energy
    from aleph.world.surface import build_surface, icosphere

    v, f = icosphere(2)
    arena = WorldArena(capacity={Kind.NODE: 400, Kind.FACE: 400, Kind.ANGLE4: 600})
    arena.claim("pad", Kind.NODE, 7)            # a non-zero offset: offset 0 is the case that always agrees
    s = build_surface(arena, "membrane", vertices=v, faces=f, radius_um=5.0)
    arena.assert_partitioned()

    rows = np.asarray(s.hinge_idx, np.int64).reshape(-1, 4)
    ours = rows[:, list(HINGE_PERMUTATION)]
    theirs = build_membrane_hinges(s.face_idx)
    assert ours.shape == theirs.shape, "the permutation must not change the hinge count"

    # THE CLAIM THIS MODULE RESTS ON.  Perturb first: on a resting icosphere every dihedral is nearly
    # flat, which is precisely where a sign error is invisible.
    rng = np.random.default_rng(7)
    pos = s.position + rng.normal(0.0, 0.05, s.position.shape)
    e_ours = membrane_bending_energy(pos, ours - s.nodes.lo, 1.0)
    e_theirs = membrane_bending_energy(pos, theirs - s.nodes.lo, 1.0)
    assert abs(e_ours - e_theirs) <= 1e-12 * abs(e_theirs), (
        f"hinge orderings disagree in energy ({e_ours} vs {e_theirs}) — the i<->j swap is NOT neutral "
        "and the orientation must be recovered from the face list after all"
    )
    # A wrong permutation must be caught by that same check, or it is testing nothing.
    wrong = rows[:, [0, 1, 2, 3]]
    assert abs(membrane_bending_energy(pos, wrong - s.nodes.lo, 1.0) - e_theirs) > 1.0, (
        "the energy check does not separate a wrong hinge order; it cannot certify the right one"
    )

    # The range guard: a topology reaching outside its own claim is refused on the host.
    bad = Surface(**{**{k: getattr(s, k) for k in s.__slots__}, "face_idx": s.face_idx - s.nodes.lo})
    try:
        _as_int32_topology("face_idx", np.asarray(bad.face_idx, np.int64), bad)
    except ValueError as exc:
        assert "outside its own claim" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a face reaching outside the population's claim must raise")

    assert is_closed_manifold(s.face_idx), "the volume law's precondition, checked the same way it is used"
    print(f"laws_bind self-check OK — {s.n_vertices} verts, {len(s.face_idx)} faces, {rows.shape[0]} "
          f"hinges; hinge-order energy delta {abs(e_ours - e_theirs):.3e} of {e_theirs:.6f}")


if __name__ == "__main__":
    _demo()
