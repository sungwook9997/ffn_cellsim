r"""Closed-form exterior-Stokes references for the ``extracellular_medium`` component (T10).

WHY THIS FILE EXISTS.  The 14th component was declared by PI decision D5-A (2026-07-28) with the
declaration explicitly buying none of the physics: until T10 lands an exterior solve, the cell's six
whole-body rigid modes are held by a NUMERICAL regulariser (``aI`` in the implicit paths) rather than by
a medium.  This module is the host-side reference the device runtime in
:mod:`aleph.engine.medium_exterior` is scored against, and it is deliberately pure NumPy so the
whole acceptance argument can be developed on a machine with no CUDA (the Mac lane) while the native
run stays a gbook obligation.

WHY A BOUNDARY INTEGRAL AND NOT A GRID.  Outside the cell the medium is free fluid: ``k -> inf`` and the
Brinkman screening length ``l_B = sqrt(k) -> inf``, so Darcy is not an approximation there, it is simply
wrong (``AC_EXECUTION_PLAN_2026-07-25.md`` §T5).  What remains is Stokes flow, and for an exterior domain
the single-layer integral representation *is* the Stokes solution — not a discretisation of the domain.
The only approximations are the surface quadrature and the blob regularisation, and both are refinable,
which is what makes the convergence gate below meaningful.  A grid solve would additionally carry a
truncated far field and an immersed-boundary blob radius, and at any affordable ``dx`` it could not
resolve the cell--substrate lubrication film where crawl traction actually lives.

THE KERNEL.  Cortez (2001) regularised Stokeslet for the blob
``phi_eps(r) = 15 eps^4 / (8 pi (r^2 + eps^2)^(7/2))``:

.. math::
    S^{eps}_{ij}(x) = \frac{1}{8 \pi \mu}
        \left[ \delta_{ij} \frac{r^2 + 2 eps^2}{(r^2 + eps^2)^{3/2}}
             + \frac{x_i x_j}{(r^2 + eps^2)^{3/2}} \right]

As ``eps -> 0`` this is the Oseen tensor; at ``r = 0`` it is finite, ``delta_ij / (4 pi mu eps)``, which is
the entire reason a regularised kernel is used at all.  Two structural facts follow and are asserted
rather than measured, because they hold for every argument: ``S`` is symmetric in ``ij``, and it is EVEN
in ``x``.  Together they make the assembled mobility ``M_{mn} = S^{eps}(x_m - x_n)`` exactly symmetric --
Lorentz reciprocity as a property of the discrete operator, not a numerical coincidence.  In IEEE
arithmetic the symmetry is bit-exact: ``(-dx_i)(-dx_j)`` and ``dx_i dx_j`` are the same product.

MAGNITUDE-FREE VS EXACT-MAGNITUDE ORACLES.  Most gates in this project must be shape gates because the
magnitude is a PI-GAP.  Exterior Stokes is the rare case that hands over exact magnitudes: translation
resistance ``6 pi mu a`` and rotation resistance ``8 pi mu a^3`` for a sphere.  Their RATIO,
``8 pi mu a^3 / 6 pi mu a = (4/3) a^2``, is free of the viscosity, which matters here because
``mu_medium`` has no KnowledgeClaim and no SourceEvidence in the knowledge base (queried 2026-07-28) --
it is a PI-GAP.  So every gate in this module is evaluable at ``mu = 1`` and nothing in it depends on the
missing parameter; only a production magnitude does.

WHAT THE REGULARISATION LENGTH IS, AND WHY IT IS NOT A TUNED CONSTANT (Magic-Number Block).  ``eps`` is
set to the mesh's own mean nearest-neighbour spacing ``h`` (:func:`mean_nearest_neighbour_spacing`), so it
is a property of the discretisation and not a free parameter: refine the surface and ``eps`` falls with
it, and the answer must converge.  :func:`sphere_drag_convergence` is what makes that a gate rather than
a claim -- it reports the measured error at a sequence of icosphere refinements and the fitted order.
The gate asserts CONVERGENCE, never an absolute error under a chosen tolerance, precisely so that no
value of ``eps`` can be selected to make something pass.

engine units: length um, force pN, velocity um/s, viscosity Pa*s.  Pa*s and pN*s/um^2 are numerically
equal (``ff/units.py``), so no conversion factor appears anywhere below.

Sanity Gate:
    * dimensional: ``[S] = 1/(mu * L)`` so ``[M F] = (pN*s/um^2 * um)^-1 * pN = um/s`` -- a velocity, and
      ``[B^T M^-1 B] = pN*s/um`` for the translation block, matching ``6 pi mu a``.
    * boundary cases: ``eps -> 0`` recovers the Oseen tensor; ``r -> 0`` stays finite; a single point
      reduces to an isolated regularised Stokeslet with self-mobility ``1/(4 pi mu eps)``.
    * conservation/invariant: ``M`` is exactly symmetric and positive definite, so the dissipation
      ``v.M^-1 v`` is strictly positive for every non-zero surface velocity; the 6x6 grand resistance is
      symmetric positive definite, which is the statement that the six rigid modes are LOADED.
    * numerical: assembled in float64; the icosphere's icosahedral symmetry forces the translation block
      to be exactly isotropic, giving a round-off-level check that needs no reference value.
    * sign-sense: the force the medium exerts on the surface OPPOSES its motion, so the hydrodynamic
      power ``v . f_medium`` is strictly negative for every non-zero ``v`` -- dissipation, never a source.
    * measurement-protocol: the resistance is reported about a stated reduction centre; translation and
      rotation blocks decouple only for a centre at the sphere's centre, which is what is used.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.laws.surface_manifold import icosphere

__all__ = [
    "grand_resistance",
    "hydrodynamic_power",
    "mean_nearest_neighbour_spacing",
    "mobility_matrix",
    "resistance_spectral_upper_bound",
    "regularized_stokeslet",
    "rigid_body_basis",
    "sphere_drag_convergence",
    "sphere_rotation_resistance",
    "sphere_translation_resistance",
    "surface_quadrature_sphere",
]


def regularized_stokeslet(
    dx: npt.ArrayLike, *, epsilon: float, viscosity: float
) -> np.ndarray:
    """Cortez regularised Stokeslet tensor evaluated on a stack of separation vectors.

    Args:
        dx: Separation vectors, shape ``(..., 3)``, in um.
        epsilon: Blob regularisation length in um; must be strictly positive.
        viscosity: Medium dynamic viscosity in Pa*s.

    Returns:
        Array of shape ``(..., 3, 3)`` holding ``S^eps_ij(dx)`` in ``um / (pN * s)``.

    Raises:
        ValueError: If ``epsilon`` or ``viscosity`` is not strictly positive and finite.
    """
    if not (np.isfinite(epsilon) and epsilon > 0.0):
        raise ValueError(f"epsilon must be positive-finite, got {epsilon!r}")
    if not (np.isfinite(viscosity) and viscosity > 0.0):
        raise ValueError(f"viscosity must be positive-finite, got {viscosity!r}")
    x = np.asarray(dx, dtype=np.float64)
    if x.shape[-1] != 3:
        raise ValueError("dx must have a trailing axis of length 3")
    r2 = np.einsum("...i,...i->...", x, x)
    denom = np.power(r2 + epsilon * epsilon, 1.5)
    identity = np.eye(3, dtype=np.float64)
    isotropic = ((r2 + 2.0 * epsilon * epsilon) / denom)[..., None, None] * identity
    dyadic = (x[..., :, None] * x[..., None, :]) / denom[..., None, None]
    return (isotropic + dyadic) / (8.0 * np.pi * viscosity)


def mobility_matrix(
    positions: npt.ArrayLike, *, epsilon: float, viscosity: float
) -> np.ndarray:
    """Assemble the dense ``(3N, 3N)`` regularised-Stokeslet mobility of a surface quadrature.

    The unknowns are point FORCES, not tractions, which is what makes the result exactly symmetric:
    no quadrature weight enters the operator, so ``M_{mn}`` and ``M_{nm}`` are the same tensor.

    Args:
        positions: Surface quadrature points, shape ``(N, 3)``, in um.
        epsilon: Blob regularisation length in um.
        viscosity: Medium dynamic viscosity in Pa*s.

    Returns:
        Symmetric positive-definite mobility ``M`` of shape ``(3N, 3N)`` mapping point forces in pN to
        surface velocities in um/s.
    """
    pts = np.asarray(positions, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    n = pts.shape[0]
    if n == 0:
        raise ValueError("an exterior medium needs at least one surface quadrature point")
    blocks = regularized_stokeslet(
        pts[:, None, :] - pts[None, :, :], epsilon=epsilon, viscosity=viscosity
    )
    return blocks.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)


def resistance_spectral_upper_bound(
    positions: npt.ArrayLike, *, epsilon: float, viscosity: float
) -> float:
    """Return a roundoff-enclosed upper bound on ``lambda_max(M^-1)``.

    This is a numerical stability bound for coupling the mechanistic exterior resistance to an explicit
    inner solve, not a material observable.  The dense matrix is the exact same discrete Cortez operator as
    the CUDA matvec.  For symmetric ``M``, Weyl's inequality bounds the true smallest eigenvalue below by the
    computed value minus a conservative float64 backward-error envelope ``gamma_n ||M||_inf``; reciprocating
    that positive lower endpoint encloses the resistance spectral radius from above.
    """
    mobility = mobility_matrix(positions, epsilon=epsilon, viscosity=viscosity)
    smallest = float(np.linalg.eigvalsh(mobility)[0])
    dimension = int(mobility.shape[0])
    eps64 = float(np.finfo(np.float64).eps)
    gamma_n = dimension * eps64 / (1.0 - dimension * eps64)
    lower = smallest - gamma_n * float(np.linalg.norm(mobility, ord=np.inf))
    if not np.isfinite(lower) or lower <= 0.0:
        raise RuntimeError("could not establish a positive roundoff-enclosed mobility eigenvalue floor")
    return 1.0 / lower


def rigid_body_basis(positions: npt.ArrayLike, *, centre: npt.ArrayLike) -> np.ndarray:
    """Return the ``(3N, 6)`` basis of the six whole-body rigid velocity fields.

    Columns 0-2 are unit translations; columns 3-5 are unit rotations about the reduction centre,
    ``v_n = e_k x (x_n - centre)``.  These are exactly the modes an exterior medium must load and a
    numerical regulariser only masks.

    Args:
        positions: Surface quadrature points, shape ``(N, 3)``, in um.
        centre: Reduction centre, shape ``(3,)``, in um.

    Returns:
        Basis of shape ``(3N, 6)``; translation columns are dimensionless, rotation columns are in um.
    """
    pts = np.asarray(positions, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    origin = np.asarray(centre, dtype=np.float64).reshape(3)
    offset = pts - origin
    n = pts.shape[0]
    basis = np.zeros((n, 3, 6), dtype=np.float64)
    for axis in range(3):
        basis[:, axis, axis] = 1.0
        unit = np.zeros(3, dtype=np.float64)
        unit[axis] = 1.0
        basis[:, :, 3 + axis] = np.cross(unit, offset)
    return basis.reshape(3 * n, 6)


def grand_resistance(
    positions: npt.ArrayLike,
    *,
    epsilon: float,
    viscosity: float,
    centre: npt.ArrayLike,
) -> np.ndarray:
    """Return the symmetric ``6x6`` grand resistance ``R = B^T M^-1 B`` of the surface.

    Prescribing a rigid surface velocity ``B U`` and solving ``M F = B U`` for the point forces gives
    the resultant force/torque ``B^T F``, hence ``R = B^T M^-1 B``.  ``R`` having six strictly positive
    eigenvalues is the precise, falsifiable form of "the rigid modes are held by physics": a numerical
    regulariser produces a multiple of the identity with no relation to the geometry, whereas ``R``'s
    translation block is ``6 pi mu a`` and its rotation block ``8 pi mu a^3`` for a sphere.

    Args:
        positions: Surface quadrature points, shape ``(N, 3)``, in um.
        epsilon: Blob regularisation length in um.
        viscosity: Medium dynamic viscosity in Pa*s.
        centre: Reduction centre for the torque balance, shape ``(3,)``, in um.

    Returns:
        ``(6, 6)`` resistance; the translation block is in pN*s/um and the rotation block in pN*s*um.
    """
    mobility = mobility_matrix(positions, epsilon=epsilon, viscosity=viscosity)
    basis = rigid_body_basis(positions, centre=centre)
    forces = np.linalg.solve(mobility, basis)
    resistance = basis.T @ forces
    return 0.5 * (resistance + resistance.T)


def hydrodynamic_power(
    positions: npt.ArrayLike,
    velocities: npt.ArrayLike,
    *,
    epsilon: float,
    viscosity: float,
) -> float:
    """Return the power the medium removes, ``-v . M^-1 v``, which must be negative for any motion.

    Args:
        positions: Surface quadrature points, shape ``(N, 3)``, in um.
        velocities: Surface velocities, shape ``(N, 3)``, in um/s.
        epsilon: Blob regularisation length in um.
        viscosity: Medium dynamic viscosity in Pa*s.

    Returns:
        Power in pN*um/s; strictly negative unless the surface velocity field is identically zero.
    """
    mobility = mobility_matrix(positions, epsilon=epsilon, viscosity=viscosity)
    flat = np.asarray(velocities, dtype=np.float64).reshape(-1)
    if flat.shape[0] != mobility.shape[0]:
        raise ValueError("velocities must carry one vector per quadrature point")
    return float(-flat @ np.linalg.solve(mobility, flat))


def sphere_translation_resistance(radius: float, viscosity: float) -> float:
    """Stokes' law translation resistance ``6 pi mu a`` in pN*s/um.

    Args:
        radius: Sphere radius in um.
        viscosity: Medium dynamic viscosity in Pa*s.

    Returns:
        Resistance coefficient such that drag force ``= coefficient * velocity``.
    """
    if not (np.isfinite(radius) and radius > 0.0):
        raise ValueError(f"radius must be positive-finite, got {radius!r}")
    if not (np.isfinite(viscosity) and viscosity > 0.0):
        raise ValueError(f"viscosity must be positive-finite, got {viscosity!r}")
    return 6.0 * np.pi * viscosity * radius


def sphere_rotation_resistance(radius: float, viscosity: float) -> float:
    """Stokes rotation resistance ``8 pi mu a^3`` in pN*s*um.

    Args:
        radius: Sphere radius in um.
        viscosity: Medium dynamic viscosity in Pa*s.

    Returns:
        Resistance coefficient such that torque ``= coefficient * angular velocity``.
    """
    if not (np.isfinite(radius) and radius > 0.0):
        raise ValueError(f"radius must be positive-finite, got {radius!r}")
    if not (np.isfinite(viscosity) and viscosity > 0.0):
        raise ValueError(f"viscosity must be positive-finite, got {viscosity!r}")
    return 8.0 * np.pi * viscosity * radius**3


def mean_nearest_neighbour_spacing(positions: npt.ArrayLike) -> float:
    """Return the mean nearest-neighbour distance of a quadrature, the DERIVED blob length.

    This is the whole Magic-Number argument for ``epsilon``: it is read off the discretisation the
    caller already built, it falls under refinement, and it is never selected against an outcome.

    Args:
        positions: Surface quadrature points, shape ``(N, 3)``, in um; at least two points.

    Returns:
        Mean nearest-neighbour spacing in um.
    """
    pts = np.asarray(positions, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] < 2:
        raise ValueError("positions must have shape (N, 3) with N >= 2")
    separation = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
    np.fill_diagonal(separation, np.inf)
    return float(np.mean(np.min(separation, axis=1)))


def surface_quadrature_sphere(
    subdivisions: int, radius: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return icosphere vertices and triangles for a spherical exterior-medium quadrature.

    The icosphere is used rather than a Fibonacci spiral because its icosahedral symmetry forces every
    invariant rank-2 tensor to be isotropic.  The translation block of :func:`grand_resistance` is then
    exactly a multiple of the identity, giving a round-off-level structural check that requires no
    reference value at all.

    Args:
        subdivisions: Icosphere subdivision level (0 gives the base icosahedron).
        radius: Sphere radius in um.

    Returns:
        Tuple ``(vertices, triangles)`` with vertices in um.
    """
    verts, tris = icosphere(subdivisions=subdivisions, radius=radius)
    return np.asarray(verts, dtype=np.float64), np.asarray(tris, dtype=np.int32)


def sphere_drag_convergence(
    *,
    radius: float,
    viscosity: float,
    subdivisions: tuple[int, ...] = (1, 2, 3),
    epsilon_ratio: float = 1.0,
) -> dict[str, object]:
    """Measure how the discrete sphere resistance approaches Stokes' law under refinement.

    This is the module's grid-invariance instrument.  It does NOT choose ``epsilon`` to minimise the
    error: ``epsilon = epsilon_ratio * h`` with ``h`` the mesh's own spacing and the ratio held fixed
    across levels, so the only thing that changes between rows is the discretisation.  A gate reading
    this reports the sequence and the fitted order; it must not silently replace the sequence with a
    single tolerance, because that is how a refinable approximation gets frozen into a magic number.

    Args:
        radius: Sphere radius in um.
        viscosity: Medium dynamic viscosity in Pa*s.
        subdivisions: Icosphere levels to evaluate, in increasing order.
        epsilon_ratio: Fixed ``epsilon / h``, identical at every level.

    Returns:
        Mapping with ``n_points``, ``spacing_um``, ``epsilon_um``, ``translation_resistance``,
        ``rotation_resistance``, ``translation_relative_error``, ``rotation_relative_error``,
        ``translation_isotropy_residual`` and the fitted ``translation_order`` in ``h``.
    """
    if len(subdivisions) < 2 or list(subdivisions) != sorted(set(subdivisions)):
        raise ValueError("subdivisions must be at least two strictly increasing levels")
    if not (np.isfinite(epsilon_ratio) and epsilon_ratio > 0.0):
        raise ValueError(f"epsilon_ratio must be positive-finite, got {epsilon_ratio!r}")
    exact_translation = sphere_translation_resistance(radius, viscosity)
    exact_rotation = sphere_rotation_resistance(radius, viscosity)
    centre = np.zeros(3, dtype=np.float64)
    rows: list[dict[str, float]] = []
    for level in subdivisions:
        verts, _ = surface_quadrature_sphere(level, radius)
        spacing = mean_nearest_neighbour_spacing(verts)
        epsilon = epsilon_ratio * spacing
        resistance = grand_resistance(
            verts, epsilon=epsilon, viscosity=viscosity, centre=centre
        )
        translation_block = resistance[:3, :3]
        rotation_block = resistance[3:, 3:]
        translation = float(np.trace(translation_block) / 3.0)
        rotation = float(np.trace(rotation_block) / 3.0)
        isotropy = float(
            np.linalg.norm(translation_block - translation * np.eye(3))
            / np.linalg.norm(translation_block)
        )
        rows.append(
            {
                "n_points": float(verts.shape[0]),
                "spacing_um": spacing,
                "epsilon_um": epsilon,
                "translation_resistance": translation,
                "rotation_resistance": rotation,
                "translation_relative_error": abs(translation - exact_translation)
                / exact_translation,
                "rotation_relative_error": abs(rotation - exact_rotation) / exact_rotation,
                "translation_isotropy_residual": isotropy,
            }
        )
    spacings = np.array([row["spacing_um"] for row in rows], dtype=np.float64)
    errors = np.array(
        [row["translation_relative_error"] for row in rows], dtype=np.float64
    )
    order = float(np.polyfit(np.log(spacings), np.log(errors), 1)[0])
    return {
        "exact_translation_resistance": exact_translation,
        "exact_rotation_resistance": exact_rotation,
        "epsilon_ratio": float(epsilon_ratio),
        "rows": rows,
        "translation_order": order,
    }
