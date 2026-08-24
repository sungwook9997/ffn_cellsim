r"""The enclosed volume of a closed surface, and the force its pressure exerts — ONE implementation.

WHY THIS MODULE EXISTS.  ``partial V / partial x`` was written three times in this repository, in two
different algebraic forms, and nothing tied them together:

  * ``components/nucleus/geometry.py:115`` — host NumPy, own-term form, the nucleoplasm constraint.
  * ``components/incumbent/membrane_pressure.py:48`` — the ONLY device version, written as the
    finite-element consistent load and fused into the Biot pressure trace so it cannot be used
    without a fluid grid.
  * ``components/incumbent/resting_balance_oracle.py:94`` — an oracle-local copy of the first.

Three copies is a maintenance observation.  **Two algebraic forms is a correctness question**, and it
is the reason this module leads with the lemma rather than the code.

THE LEMMA, AND THE CONDITION IT NEEDS.  For a closed outward-wound triangulation the discrete enclosed
volume is ``V = (1/6) sum_f x_a . (x_b x x_c)``.  Differentiating gives each face's OWN-TERM
contribution to each of its three vertices::

    dV/dx_a |_f = (1/6)(x_b x x_c),   dV/dx_b |_f = (1/6)(x_c x x_a),   dV/dx_c |_f = (1/6)(x_a x x_b)

The consistent-load form instead gives EVERY vertex of the face the same vector ``A n / 3``, i.e.
``((x_b - x_a) x (x_c - x_a)) / 6``, which expands to exactly the SUM of the three own terms above.
Per face the two disagree; their difference at vertex ``a`` is ``(1/6) x_a x (x_b - x_c)``.

Assembled over every face incident on ``a``, the link of ``a`` in a closed manifold mesh is a CYCLE
``v_1 -> v_2 -> ... -> v_n -> v_1`` and the incident faces are ``(a, v_m, v_{m+1})``, so that difference
sums to ``(1/6) x_a x sum_m (v_m - v_{m+1}) = 0`` — it telescopes.  **The two forms are therefore equal
as assembled nodal fields on a closed manifold mesh, and unequal on anything else.**  An open sheet, a
mesh with a boundary edge, or a launch that covers only some of a surface's faces breaks the telescoping
and leaves a force that depends on where the origin happens to be.  ``_demo`` asserts both halves of
that: equality on the closed icosphere, and inequality on an open patch cut out of it.

WHICH FORM THIS MODULE SHIPS, and why it is the own-term one.  Both are correct when the condition
holds, so the tie is broken on what fails loudly.  The own-term form is ``-dE/dx`` of a stated energy
term by construction, so a caller that hands it an open surface gets an origin-dependent answer that a
translation test catches immediately.  The consistent-load form is origin-dependent in the same case but
*looks* like a surface traction, which is the shape a wrong answer hides in.

WHAT THIS MODULE DELIBERATELY DOES NOT DO.  It carries no pressure.  ``delta_p`` is an argument, never a
constant here: the osmotic value, the nucleoplasm bulk modulus and the Biot pore pressure are three
different physical quantities that all multiply this same geometric gradient, and folding any of them in
is what produced a fused kernel that needs a fluid grid to compute a triangle's normal.  It also reads
no field and interpolates nothing — ``membrane_pressure.py``'s spatial pressure TRACE is a separate
concern that composes with this by supplying ``delta_p``.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — ``pos`` [um]; the gradient [um^2]; ``delta_p`` [pN/um^2] (== Pa in engine units);
    the accumulated force [pN].  ``volume`` [um^3].
  * boundary — ``delta_p = 0`` gives bit-zero force.  Reversing the sign of ``delta_p`` reverses every
    force exactly.  A closed surface under uniform pressure has ZERO net force and ZERO net torque, and
    both are asserted rather than assumed.  An OPEN surface is outside the lemma: the module does not
    silently accept one, it offers :func:`is_closed_manifold` and the demo shows what breaks.
  * conservation/invariant — the assembled field is TRANSLATION-INVARIANT on a closed mesh even though
    each face's term is not, because ``V`` itself is; the demo asserts it under a large offset, which is
    the test that fails first if a launch covers a partial face range.  The virial identity
    ``sum_i f_i . x_i = 3 delta_p V`` arbitrates the sign and the factor together.
  * CFL/precision — nothing is integrated.  float64 throughout, matching the production arrays.  The
    device kernel accumulates atomically, so summation order is not reproducible and exact float
    equality against the host twin is NOT claimed; a Higham bound is.
  * sign sense — positive ``delta_p`` means higher pressure INSIDE, and pushes each vertex OUTWARD along
    the area-weighted outward normal, growing ``V``.  The virial identity is the arbiter, not a picture.
  * measurement protocol — the kernel is one launch per closed surface over that surface's OWN face
    range; the host twin is for acceptance and rendering between accepted steps, never in the loop.

engine units: length um, area um^2, volume um^3, pressure pN/um^2, force pN.
Runtime: NVIDIA Warp on CUDA for the kernel; the host twin is pure NumPy and CPU-importable.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import warp as wp

__all__ = [
    "volume_pressure_kernel",
    "mesh_volume",
    "volume_gradient",
    "pressure_force_reference",
    "consistent_load_reference",
    "is_closed_manifold",
]


@wp.kernel
def volume_pressure_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    delta_p: wp.float64,
    force: wp.array(dtype=wp.vec3d),
) -> None:
    """Accumulate ``f_i += delta_p * dV/dx_i`` over one closed surface's triangles.

    ``faces`` holds GLOBAL node indices, so this composes additively with every other law writing into
    the same ``force`` array and needs no offset argument — which is the whole point of addressing a
    population as a range rather than as a private array.

    One thread per FACE, three atomic adds, each vertex receiving its OWN term.  Launch over the
    surface's complete face range: a partial launch is not a partial force, it is a different one, because
    the translation invariance of the assembled field depends on the vertex links closing.
    """
    t = wp.tid()
    ia = faces[t, 0]
    ib = faces[t, 1]
    ic = faces[t, 2]
    xa = pos[ia]
    xb = pos[ib]
    xc = pos[ic]
    s = delta_p / wp.float64(6.0)
    wp.atomic_add(force, ia, s * wp.cross(xb, xc))
    wp.atomic_add(force, ib, s * wp.cross(xc, xa))
    wp.atomic_add(force, ic, s * wp.cross(xa, xb))


def mesh_volume(pos: npt.ArrayLike, faces: npt.ArrayLike) -> float:
    """Enclosed volume [um^3] of a closed outward-wound triangulation, by the divergence theorem.

    Computed from the polyhedron the mesh actually is and never from the analytic shape it approximates:
    an inscribed mesh encloses LESS than its sphere (~1.4% at subdivision 2), and substituting the
    analytic value produced a phantom nucleus force at t0 in this repository once already.
    """
    x = np.asarray(pos, np.float64).reshape(-1, 3)
    f = np.asarray(faces, np.int64).reshape(-1, 3)
    return float(np.einsum("ij,ij->i", x[f[:, 0]], np.cross(x[f[:, 1]], x[f[:, 2]])).sum() / 6.0)


def volume_gradient(pos: npt.ArrayLike, faces: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """``dV/dx_i`` as an ``(N, 3)`` host array [um^2] — the host twin of :func:`volume_pressure_kernel`.

    The own-term form, assembled by scatter-add exactly as the kernel accumulates it.
    """
    x = np.asarray(pos, np.float64).reshape(-1, 3)
    f = np.asarray(faces, np.int64).reshape(-1, 3)
    xa, xb, xc = x[f[:, 0]], x[f[:, 1]], x[f[:, 2]]
    grad = np.zeros_like(x)
    np.add.at(grad, f[:, 0], np.cross(xb, xc) / 6.0)
    np.add.at(grad, f[:, 1], np.cross(xc, xa) / 6.0)
    np.add.at(grad, f[:, 2], np.cross(xa, xb) / 6.0)
    return grad


def pressure_force_reference(
    pos: npt.ArrayLike, faces: npt.ArrayLike, delta_p: float
) -> npt.NDArray[np.float64]:
    """``f_i = delta_p * dV/dx_i`` [pN] — the acceptance oracle for the kernel."""
    return float(delta_p) * volume_gradient(pos, faces)


def consistent_load_reference(
    pos: npt.ArrayLike, faces: npt.ArrayLike, delta_p: float
) -> npt.NDArray[np.float64]:
    """The OTHER algebraic form: ``delta_p * A n / 3`` given identically to all three vertices [pN].

    This is what ``components/incumbent/membrane_pressure.py:48`` assembles.  It is kept here — not to be
    called by the runtime, but so the lemma in the module docstring is testable rather than asserted, and
    so a parity run against the incumbent has a named thing to compare.  Equal to
    :func:`pressure_force_reference` on a closed manifold mesh; NOT equal on anything else.
    """
    x = np.asarray(pos, np.float64).reshape(-1, 3)
    f = np.asarray(faces, np.int64).reshape(-1, 3)
    area2 = np.cross(x[f[:, 1]] - x[f[:, 0]], x[f[:, 2]] - x[f[:, 0]])
    per_vertex = float(delta_p) * area2 / 6.0
    out = np.zeros_like(x)
    for col in range(3):
        np.add.at(out, f[:, col], per_vertex)
    return out


def is_closed_manifold(faces: npt.ArrayLike) -> bool:
    """Whether every edge is shared by exactly two faces — the condition the lemma needs.

    A caller that cannot answer yes is outside this module's guarantee and will get an origin-dependent
    force from either algebraic form.
    """
    f = np.asarray(faces, np.int64).reshape(-1, 3)
    edges = np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    return bool(np.all(counts == 2))


def _demo() -> None:
    """Self-check: the lemma, its condition, and the invariants the Sanity Gate names.  No device."""
    from aleph.laws.surface_manifold import icosphere

    verts, faces = icosphere(2, 5.0)          # radius in the caller's own length unit; um here
    verts = np.ascontiguousarray(verts, np.float64)
    faces = np.ascontiguousarray(faces, np.int64)
    assert is_closed_manifold(faces), "the icosphere must be a closed manifold or the lemma does not apply"

    dp = 40.0                                  # pN/um^2; a pressure JUMP, not a sourced constant
    own = pressure_force_reference(verts, faces, dp)
    fem = consistent_load_reference(verts, faces, dp)

    # THE LEMMA — the two algebraic forms agree as assembled fields on a closed mesh.
    scale = np.abs(own).max()
    assert np.abs(own - fem).max() <= 64 * np.finfo(np.float64).eps * scale, "closed-mesh forms disagree"

    # Its CONDITION — cut the mesh open and they must NOT agree, or the lemma proves nothing.
    open_faces = faces[: faces.shape[0] // 2]
    assert not is_closed_manifold(open_faces)
    d_open = np.abs(pressure_force_reference(verts, open_faces, dp)
                    - consistent_load_reference(verts, open_faces, dp)).max()
    assert d_open > 1e-6 * scale, "an open patch must separate the two forms; it is the whole condition"

    # Translation invariance of the ASSEMBLED field — V is translation-invariant on a closed mesh, so
    # its gradient must be, even though no single face's term is.  This is the check that fails first
    # when a launch covers only part of a surface's face range.
    shifted = pressure_force_reference(verts + np.array([1e3, -7e2, 3e2]), faces, dp)
    assert np.abs(shifted - own).max() <= 1e-6 * scale, "assembled gradient is not translation-invariant"

    # Uniform pressure on a closed surface: zero net force AND zero net torque.
    assert np.abs(own.sum(axis=0)).max() <= 1e-9 * scale
    assert np.abs(np.cross(verts, own).sum(axis=0)).max() <= 1e-9 * scale * 5.0

    # Sign and factor together: the virial identity sum_i f_i . x_i = 3 delta_p V.
    vol = mesh_volume(verts, faces)
    virial = float(np.einsum("ij,ij->", own, verts))
    assert abs(virial - 3.0 * dp * vol) <= 1e-9 * abs(3.0 * dp * vol), "virial identity fails"
    assert vol < 4.0 / 3.0 * np.pi * 125.0, "an inscribed mesh must enclose LESS than its sphere"

    # delta_p = 0 is bit-zero, and the law is exactly odd in delta_p.
    assert np.array_equal(pressure_force_reference(verts, faces, 0.0), np.zeros_like(verts))
    assert np.array_equal(pressure_force_reference(verts, faces, -dp), -own)

    print(f"laws.volume self-check OK — V={vol:.6f} um^3, max|f|={scale:.6f} pN, "
          f"open-patch separation={d_open / scale:.3e} (relative)")


if __name__ == "__main__":
    _demo()
