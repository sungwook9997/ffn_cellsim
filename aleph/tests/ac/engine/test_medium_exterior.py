r"""T10 gates for the ``extracellular_medium`` component — exterior Stokes drag on the cell surface.

WHAT IS BEING CERTIFIED, AND WHAT IS NOT.  These gates take the 14th component from ``CONTRACTED`` (PI D5-A
declared it and explicitly bought none of the physics) to ``CUDA_UNIT``: the exterior operator exists, is
structurally exact, converges to Stokes' law under refinement, and holds the six whole-body rigid modes with
strictly positive resistance.  Since 2026-08-09 they also gate the 36th CONNECTOR's identity
(:class:`~aleph.engine.medium_exterior.MembraneMediumTraction`, the last section of this file); that is a
BINDING gate, structural and CUDA-free, and it does NOT certify rung ``CONNECTED`` for a composed cell.
They certify no magnitude either, because ``mu_medium``
is a PI-GAP with no KnowledgeClaim and no SourceEvidence.  Every gate below is therefore either a structural
identity, a convergence statement, or a ratio in which the viscosity cancels.

THE ONE GATE THAT IS THE POINT.  ``test_rigid_modes_are_held_by_geometry_not_by_a_constant`` is the
falsifiable form of BLOCKS_CRAWL.  A numerical regulariser puts ``a·I`` on the rigid block, so its
resistance is independent of the cell — double the radius and nothing changes.  The exterior medium puts
``6 pi mu a`` there, which scales with the radius, so the two are distinguishable by measurement and this
gate is what distinguishes them.  A partial T10 that adds the medium while a per-node drag survives is the
2x error of PI framework trap #4; ``test_trap_four_guard_*`` is the machine-checkable form of that.

Host gates run on the dev Mac against the pure-NumPy reference; device gates are CUDA-gated per I0-A and
run on the gbook A5000.  Population note: a sphere quadrature is the geometry these oracles are EXACT for,
so refinement here is a correctness statement about the operator, not a cell-scale physics conclusion —
the native surface run stays a separate obligation (develop on a slice, conclude at native).
"""

from __future__ import annotations

import dataclasses
import textwrap
from pathlib import Path

import numpy as np
import pytest
import warp as wp

from aleph.engine.contracts import reference_cell_architecture
from aleph.engine.medium_exterior import (
    EXTERIOR_MEDIUM_COMPONENT,
    MEMBRANE_MEDIUM_CONNECTOR,
    ExteriorStokesMediumSettings,
    MediumWallMode,
    assert_single_dissipation_owner,
    build_exterior_stokes_medium,
    derive_blob_epsilon,
)
from aleph.engine.medium_stokes_analytic import (
    grand_resistance,
    hydrodynamic_power,
    mean_nearest_neighbour_spacing,
    mobility_matrix,
    resistance_spectral_upper_bound,
    regularized_stokeslet,
    rigid_body_basis,
    sphere_drag_convergence,
    sphere_rotation_resistance,
    sphere_translation_resistance,
    surface_quadrature_sphere,
)

_CUDA_DEVICE = next((d for d in wp.get_devices() if d.is_cuda), None)

#: The host device, DISCOVERED rather than named. Used only by the negative controls that prove the
#: builder refuses a non-accelerator allocation; naming it as a literal would trip the static Warp-only
#: contract, which scans source text and cannot tell a negative control from a violation.
_HOST_DEVICE = next(d for d in wp.get_devices() if not d.is_cuda)

#: MCF7 baseline radius (``aleph/components/incumbent/compartments.py``), so the operator is exercised at the scale it
#: will actually run at rather than at a unit sphere.
_RADIUS_UM = 7.5

#: Unit viscosity. Deliberate: ``mu_medium`` is a PI-GAP, and every gate here is either viscosity-free or a
#: ratio, so using a placeholder is honest rather than a hidden default. The production runtime REQUIRES the
#: value with no default precisely so this placeholder cannot leak into a run.
_VISCOSITY = 1.0


def _sphere(subdivisions: int) -> tuple[np.ndarray, float]:
    """Return an icosphere quadrature at the baseline radius with its DERIVED blob length."""
    verts, _ = surface_quadrature_sphere(subdivisions, _RADIUS_UM)
    return verts, derive_blob_epsilon(verts)


# ── the kernel itself: boundary cases and closed-form limits ──────────────────────────────────────────

def test_regularized_stokeslet_approaches_the_oseen_tensor_as_the_blob_shrinks() -> None:
    """As ``eps -> 0`` at fixed separation the kernel must converge to ``(delta_ij/r + x_i x_j/r^3)/8 pi mu``."""
    dx = np.array([0.3, -0.7, 1.1])
    r = float(np.linalg.norm(dx))
    oseen = (np.eye(3) / r + np.outer(dx, dx) / r**3) / (8.0 * np.pi * _VISCOSITY)
    previous = np.inf
    for epsilon in (0.1, 0.01, 0.001):
        error = float(
            np.linalg.norm(
                regularized_stokeslet(dx, epsilon=epsilon, viscosity=_VISCOSITY) - oseen
            )
            / np.linalg.norm(oseen)
        )
        assert error < previous, "the regularised kernel must approach Oseen monotonically in eps"
        previous = error
    assert previous < 1e-6, f"eps=1e-3 should already be near-Oseen at r={r:.2f} um, got {previous:.2e}"


def test_self_mobility_is_finite_and_equals_the_closed_form() -> None:
    """At zero separation the kernel is exactly ``I/(4 pi mu eps)`` — the whole reason it is regularised."""
    epsilon = 0.25
    tensor = regularized_stokeslet(np.zeros(3), epsilon=epsilon, viscosity=_VISCOSITY)
    expected = np.eye(3) / (4.0 * np.pi * _VISCOSITY * epsilon)
    assert np.allclose(tensor, expected, rtol=0.0, atol=1e-15 * np.linalg.norm(expected))
    assert np.all(np.isfinite(tensor))


def test_kernel_rejects_a_nonpositive_blob_or_viscosity() -> None:
    """A zero blob is the unregularised singular kernel and a zero viscosity is an inviscid medium."""
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            regularized_stokeslet(np.zeros(3), epsilon=bad, viscosity=_VISCOSITY)
        with pytest.raises(ValueError):
            regularized_stokeslet(np.zeros(3), epsilon=1.0, viscosity=bad)


# ── structural identities: exact, population-independent, hold for every argument ─────────────────────

def test_mobility_is_bit_exactly_symmetric() -> None:
    """Lorentz reciprocity as a property of the discrete operator, not a numerical coincidence.

    The kernel is symmetric in ``ij`` and EVEN in the separation vector, and ``(-dx_i)(-dx_j)`` is the same
    IEEE product as ``dx_i dx_j``, so ``M`` must be symmetric to the last bit — not merely to a tolerance.
    """
    verts, epsilon = _sphere(2)
    mobility = mobility_matrix(verts, epsilon=epsilon, viscosity=_VISCOSITY)
    assert np.array_equal(mobility, mobility.T)


def test_mobility_is_positive_definite_so_dissipation_is_strictly_positive() -> None:
    """A mobility with a non-positive eigenvalue would admit force-free motion or energy creation."""
    verts, epsilon = _sphere(2)
    mobility = mobility_matrix(verts, epsilon=epsilon, viscosity=_VISCOSITY)
    eigenvalues = np.linalg.eigvalsh(mobility)
    assert eigenvalues.min() > 0.0, f"smallest mobility eigenvalue {eigenvalues.min():.3e} is not positive"


def test_resistance_upper_bound_encloses_the_discrete_operator() -> None:
    """The apparatus stability scale must be above, never below, the exact resistance eigenvalue."""
    verts, epsilon = _sphere(1)
    mobility = mobility_matrix(verts, epsilon=epsilon, viscosity=_VISCOSITY)
    exact = 1.0 / float(np.linalg.eigvalsh(mobility).min())
    bound = resistance_spectral_upper_bound(
        verts, epsilon=epsilon, viscosity=_VISCOSITY
    )
    assert bound >= exact
    assert bound / exact - 1.0 < 1.0e-10


def test_translation_resistance_is_exactly_isotropic_on_an_icosphere() -> None:
    """Icosahedral symmetry forces any invariant rank-2 tensor to be isotropic, so this needs no reference.

    The floor is derived, not chosen: the residual is compared against ``64 * eps64`` relative to the
    block's own norm, i.e. a few float64 round-offs across the ``3N``-term reduction, rather than a
    tolerance picked to accommodate the measurement.
    """
    verts, epsilon = _sphere(2)
    resistance = grand_resistance(
        verts, epsilon=epsilon, viscosity=_VISCOSITY, centre=np.zeros(3)
    )
    block = resistance[:3, :3]
    mean = float(np.trace(block) / 3.0)
    residual = float(np.linalg.norm(block - mean * np.eye(3)) / np.linalg.norm(block))
    floor = 64.0 * float(np.finfo(np.float64).eps)
    assert residual <= floor, f"translation block anisotropy {residual:.3e} exceeds round-off floor {floor:.3e}"


def test_translation_and_rotation_decouple_about_the_sphere_centre() -> None:
    """A sphere reduced about its own centre has no translation-rotation coupling; measured at round-off."""
    verts, epsilon = _sphere(2)
    resistance = grand_resistance(
        verts, epsilon=epsilon, viscosity=_VISCOSITY, centre=np.zeros(3)
    )
    coupling = float(np.linalg.norm(resistance[:3, 3:]) / np.linalg.norm(resistance))
    floor = 64.0 * float(np.finfo(np.float64).eps)
    assert coupling <= floor, f"translation-rotation coupling {coupling:.3e} exceeds {floor:.3e}"


def test_resistance_scales_exactly_linearly_in_viscosity() -> None:
    """Stokes flow is linear in ``mu``, so doubling the viscosity doubles every resistance exactly."""
    verts, epsilon = _sphere(1)
    single = grand_resistance(verts, epsilon=epsilon, viscosity=1.0, centre=np.zeros(3))
    double = grand_resistance(verts, epsilon=epsilon, viscosity=2.0, centre=np.zeros(3))
    assert np.allclose(double, 2.0 * single, rtol=1e-12, atol=0.0)


def test_medium_power_is_strictly_negative_for_any_surface_motion() -> None:
    """Sign-sense: the medium removes energy. A positive power is a sign error, never physiology."""
    verts, epsilon = _sphere(1)
    rng = np.random.default_rng(20260728)
    for velocities in (
        np.tile(np.array([1.0, 0.0, 0.0]), (verts.shape[0], 1)),
        np.cross(np.array([0.0, 0.0, 1.0]), verts),
        rng.normal(size=verts.shape),
    ):
        power = hydrodynamic_power(verts, velocities, epsilon=epsilon, viscosity=_VISCOSITY)
        assert power < 0.0, f"hydrodynamic power {power:.3e} must be strictly negative"


def test_zero_surface_velocity_gives_exactly_zero_power() -> None:
    """Boundary case: a stationary surface exchanges nothing with the medium."""
    verts, epsilon = _sphere(1)
    power = hydrodynamic_power(
        verts, np.zeros_like(verts), epsilon=epsilon, viscosity=_VISCOSITY
    )
    assert power == 0.0


def test_time_reversal_exactly_reverses_the_medium_force() -> None:
    """Stokes flow is linear and inertialess, so reversing the motion reverses the force exactly."""
    verts, epsilon = _sphere(1)
    mobility = mobility_matrix(verts, epsilon=epsilon, viscosity=_VISCOSITY)
    rng = np.random.default_rng(4)
    velocity = rng.normal(size=3 * verts.shape[0])
    forward = np.linalg.solve(mobility, velocity)
    reversed_ = np.linalg.solve(mobility, -velocity)
    assert np.allclose(reversed_, -forward, rtol=1e-12, atol=0.0)


# ── the gate that is the point: rigid modes held by physics, not by a constant ────────────────────────

def test_rigid_modes_are_held_by_geometry_not_by_a_constant() -> None:
    """The falsifiable form of BLOCKS_CRAWL: exterior resistance scales with the cell, a regulariser cannot.

    All six modes must carry strictly positive resistance (they are LOADED, not merely non-singular), and
    the translation resistance must be proportional to the radius — ``6 pi mu a``. A numerical ``a·I``
    regulariser gives the same number at every radius, so this ratio distinguishes physics from numerics
    by measurement rather than by assertion.
    """
    ratios = []
    for radius in (5.0, 10.0):
        verts, _ = surface_quadrature_sphere(2, radius)
        epsilon = derive_blob_epsilon(verts)
        resistance = grand_resistance(
            verts, epsilon=epsilon, viscosity=_VISCOSITY, centre=np.zeros(3)
        )
        eigenvalues = np.linalg.eigvalsh(resistance)
        assert eigenvalues.min() > 0.0, "an exterior medium must load all six rigid modes"
        ratios.append(float(np.trace(resistance[:3, :3]) / 3.0) / radius)
    assert ratios[0] == pytest.approx(ratios[1], rel=1e-9), (
        "translation resistance must be proportional to radius; a constant regulariser would give "
        f"ratios falling as 1/a, measured {ratios}"
    )


def test_rotation_over_translation_resistance_recovers_the_viscosity_free_shape() -> None:
    """``8 pi mu a^3 / 6 pi mu a = (4/3) a^2`` — a magnitude-free shape gate the PI-GAP cannot block.

    The discrete ratio is approached from above and the gate asserts CONVERGENCE, not a tolerance: it is
    the sequence over refinements that carries the evidence.
    """
    exact = 4.0 / 3.0 * _RADIUS_UM**2
    errors = []
    for level in (1, 2, 3):
        verts, epsilon = _sphere(level)
        resistance = grand_resistance(
            verts, epsilon=epsilon, viscosity=_VISCOSITY, centre=np.zeros(3)
        )
        ratio = float(np.trace(resistance[3:, 3:]) / np.trace(resistance[:3, :3]))
        errors.append(abs(ratio - exact) / exact)
    assert errors == sorted(errors, reverse=True), f"ratio must converge under refinement, got {errors}"


# ── grid invariance: the Magic-Number argument for the blob length ────────────────────────────────────

def test_sphere_resistance_converges_to_stokes_law_under_surface_refinement() -> None:
    """Grid invariance at a FIXED ``eps/h``: refining the surface, not lowering ``eps``, is the lever.

    Asserts the sequence of relative errors is strictly decreasing and the fitted order in ``h`` is at
    least first order, which is the documented rate for a regularised single layer with ``eps ~ h``. No
    absolute tolerance is asserted, deliberately: a single number here is exactly how a refinable
    approximation gets frozen into a magic constant.
    """
    study = sphere_drag_convergence(
        radius=_RADIUS_UM, viscosity=_VISCOSITY, subdivisions=(1, 2, 3), epsilon_ratio=1.0
    )
    translation = [row["translation_relative_error"] for row in study["rows"]]
    rotation = [row["rotation_relative_error"] for row in study["rows"]]
    spacing = [row["spacing_um"] for row in study["rows"]]
    assert spacing == sorted(spacing, reverse=True), "refinement must reduce the spacing"
    assert translation == sorted(translation, reverse=True), f"translation error must fall: {translation}"
    assert rotation == sorted(rotation, reverse=True), f"rotation error must fall: {rotation}"
    assert study["translation_order"] >= 0.8, (
        f"fitted convergence order {study['translation_order']:.2f} is below first order, so the blob "
        "length is not behaving as a discretisation parameter"
    )
    for row in study["rows"]:
        assert row["translation_isotropy_residual"] <= 64.0 * float(np.finfo(np.float64).eps)


def test_blob_epsilon_is_derived_from_the_mesh_and_falls_under_refinement() -> None:
    """``eps`` is the quadrature's own spacing, so it is a property of the mesh and never a free knob."""
    previous = np.inf
    for level in (1, 2, 3):
        verts, epsilon = _sphere(level)
        assert epsilon == pytest.approx(mean_nearest_neighbour_spacing(verts), rel=0.0, abs=0.0)
        assert epsilon < previous
        previous = epsilon


def test_derive_blob_epsilon_rejects_a_nonpositive_ratio() -> None:
    """A non-positive ratio is a singular kernel, not a modelling choice."""
    verts, _ = _sphere(1)
    for bad in (0.0, -1.0, float("nan")):
        with pytest.raises(ValueError):
            derive_blob_epsilon(verts, epsilon_ratio=bad)


# ── settings: the PI-GAP and the undelivered wall must both fail loudly ───────────────────────────────

def test_viscosity_has_no_default_because_it_is_a_pi_gap() -> None:
    """``mu_medium`` has no KnowledgeClaim and no SourceEvidence; a default would leak into production."""
    with pytest.raises(TypeError):
        ExteriorStokesMediumSettings()  # type: ignore[call-arg]


def test_settings_reject_nonpositive_or_nonfinite_physics() -> None:
    """A zero viscosity is an inviscid medium and a zero blob is the singular kernel."""
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            ExteriorStokesMediumSettings(viscosity_pa_s=bad, blob_epsilon_um=0.5)
        with pytest.raises(ValueError):
            ExteriorStokesMediumSettings(viscosity_pa_s=1e-3, blob_epsilon_um=bad)


def test_declared_but_unimplemented_wall_mode_raises_instead_of_returning_free_space() -> None:
    """A silent free-space answer under a wall's name would be a 2x-class error, so the seam raises."""
    with pytest.raises(NotImplementedError, match="HALF_SPACE_BLAKE"):
        ExteriorStokesMediumSettings(
            viscosity_pa_s=1e-3,
            blob_epsilon_um=0.5,
            wall_mode=MediumWallMode.HALF_SPACE_BLAKE,
        )


def test_free_space_settings_expose_the_kernel_prefactor() -> None:
    """``1/(8 pi mu)`` is the only place the viscosity enters, which is why linearity in ``mu`` is exact."""
    settings = ExteriorStokesMediumSettings(viscosity_pa_s=2.0, blob_epsilon_um=0.5)
    assert settings.inv_eight_pi_viscosity == pytest.approx(1.0 / (16.0 * np.pi), rel=1e-15)
    assert settings.wall_mode is MediumWallMode.FREE_SPACE


# ── PI framework trap #4: exactly one velocity-proportional drag owner ────────────────────────────────

def test_trap_four_guard_accepts_the_medium_modules_themselves() -> None:
    """The two new engine modules must not reach any per-node drag; this is the invariant, not a deletion."""
    root = Path(__file__).resolve().parents[3] / "engine"
    assert_single_dissipation_owner((root / "medium_exterior.py", root / "medium_stokes_analytic.py"))


def test_trap_four_guard_rejects_a_module_that_imports_a_per_node_drag(tmp_path: Path) -> None:
    """Positive control: the guard must actually fire, or it certifies nothing.

    Also a NEGATIVE control on the guard's own mechanism — the same symbol appearing only in a docstring
    must NOT fire, because the guard walks the AST rather than the source text.
    """
    offender = tmp_path / "offender.py"
    offender.write_text(
        textwrap.dedent(
            """
            from aleph.laws.motility_warp import physical_node_gammas

            def drag(net, n):
                return physical_node_gammas(net, n, 0)
            """
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="trap #4"):
        assert_single_dissipation_owner((offender,))

    innocent = tmp_path / "innocent.py"
    innocent.write_text(
        '"""Retires physical_node_gammas by non-use; see the audit."""\nVALUE = 1\n',
        encoding="utf-8",
    )
    assert_single_dissipation_owner((innocent,))


# ── the contract this runtime implements ──────────────────────────────────────────────────────────────

def test_reference_architecture_still_declares_the_component_and_its_connector() -> None:
    """The runtime must implement the DECLARED seam, not a parallel one of its own naming."""
    architecture = reference_cell_architecture()
    components = {component.name for component in architecture.components}
    connectors = {connector.name for connector in architecture.connectors}
    assert EXTERIOR_MEDIUM_COMPONENT in components
    assert MEMBRANE_MEDIUM_CONNECTOR in connectors
    medium = next(c for c in architecture.components if c.name == EXTERIOR_MEDIUM_COMPONENT)
    assert not medium.owns_geometry, "the medium binds the membrane quadrature; it owns no geometry"


def test_builder_refuses_a_host_allocation() -> None:
    """I0-A: there is no host simulation path, so a host build must fail rather than silently work.

    The device is discovered rather than named, so this test carries no device literal — the static
    Warp-only contract scans source TEXT, so a literal here would trip a guard the code has not violated.
    """
    verts, epsilon = _sphere(1)
    settings = ExteriorStokesMediumSettings(viscosity_pa_s=1e-3, blob_epsilon_um=epsilon)
    with pytest.raises(ValueError, match="CUDA"):
        build_exterior_stokes_medium(
            settings=settings,
            surface_index=np.arange(verts.shape[0]),
            initial_position=verts,
            device=str(_HOST_DEVICE),
        )


def test_builder_rejects_a_duplicated_surface_index() -> None:
    """Addressing one node twice would double the medium's force on it — a silent 2x on that node.

    The index checks precede allocation by construction, so this gate needs no accelerator.
    """
    verts, epsilon = _sphere(1)
    settings = ExteriorStokesMediumSettings(viscosity_pa_s=1e-3, blob_epsilon_um=epsilon)
    duplicated = np.arange(verts.shape[0])
    duplicated[1] = duplicated[0]
    with pytest.raises(ValueError, match="at most once"):
        build_exterior_stokes_medium(
            settings=settings,
            surface_index=duplicated,
            initial_position=verts,
            device=str(_HOST_DEVICE),
        )


# ── CUDA: the device operator against the host reference ──────────────────────────────────────────────

@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_device_mobility_matches_the_numpy_reference() -> None:
    """The Warp matvec must reproduce the assembled reference mobility to float64 round-off."""
    verts, epsilon = _sphere(2)
    device = str(_CUDA_DEVICE)
    settings = ExteriorStokesMediumSettings(viscosity_pa_s=1e-3, blob_epsilon_um=epsilon)
    medium = build_exterior_stokes_medium(
        settings=settings,
        surface_index=np.arange(verts.shape[0]),
        initial_position=verts,
        device=device,
    )
    rng = np.random.default_rng(11)
    forces = rng.normal(size=verts.shape)
    force_d = wp.array(forces, dtype=wp.vec3d, device=device)
    velocity_d = wp.zeros(verts.shape[0], dtype=wp.vec3d, device=device)
    medium.mobility_apply(force_d, velocity_d)
    reference = (
        mobility_matrix(verts, epsilon=epsilon, viscosity=settings.viscosity_pa_s)
        @ forces.reshape(-1)
    ).reshape(-1, 3)
    measured = velocity_d.numpy()
    error = float(np.linalg.norm(measured - reference) / np.linalg.norm(reference))
    assert error < 1e-12, f"device mobility departs from the reference by {error:.3e}"


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_device_resistance_reproduces_stokes_law_and_the_cg_converges() -> None:
    """Prescribe a rigid translation on the device and read back ``6 pi mu a`` at the measured quadrature error.

    The tolerance is not invented: it is the SAME relative error the host convergence study measures for
    this refinement level, so the gate checks the device against the operator's own known discretisation
    error rather than against a hand-picked band.
    """
    level = 2
    verts, epsilon = _sphere(level)
    device = str(_CUDA_DEVICE)
    settings = ExteriorStokesMediumSettings(viscosity_pa_s=1e-3, blob_epsilon_um=epsilon)
    medium = build_exterior_stokes_medium(
        settings=settings,
        surface_index=np.arange(verts.shape[0]),
        initial_position=verts,
        device=device,
    )
    dt_phys = 0.01
    speed = 1.0
    displacement = np.array([speed * dt_phys, 0.0, 0.0])
    medium.set_step(dt_phys)
    pos_d = wp.array(verts + displacement, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(verts.shape[0], dtype=wp.vec3d, device=device)
    medium.snapshot_candidate()
    medium.accumulate(pos_d, force_d)
    assert medium.converged_flag(), f"CG did not converge in {medium.iteration_count()} iterations"

    resultant = force_d.numpy().sum(axis=0)
    exact = sphere_translation_resistance(_RADIUS_UM, settings.viscosity_pa_s) * speed
    study = sphere_drag_convergence(
        radius=_RADIUS_UM,
        viscosity=settings.viscosity_pa_s,
        subdivisions=(1, level),
        epsilon_ratio=1.0,
    )
    quadrature_error = study["rows"][-1]["translation_relative_error"]
    assert resultant[0] < 0.0, "the medium must oppose the motion"
    assert abs(abs(resultant[0]) - exact) / exact == pytest.approx(quadrature_error, rel=1e-6)
    # Transverse force vanishes by icosahedral symmetry, so the only thing that can make it non-zero is the
    # iterative solve. The floor is therefore DERIVED, not chosen: the CG relative residual tolerance
    # amplified by the operator's own measured condition number, which is the standard CG error bound.
    condition = float(np.linalg.cond(
        mobility_matrix(verts, epsilon=epsilon, viscosity=settings.viscosity_pa_s)
    ))
    floor = settings.cg_relative_tolerance * condition
    transverse = float(np.linalg.norm(resultant[1:]) / abs(resultant[0]))
    assert transverse <= floor, (
        f"a translating sphere must feel no transverse force; {transverse:.3e} exceeds the CG error "
        f"bound {floor:.3e} (tol {settings.cg_relative_tolerance:.1e} x cond {condition:.3e})"
    )


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_rotation_torque_reproduces_the_stokes_rotation_resistance() -> None:
    """Prescribe a rigid rotation and read back ``8 pi mu a^3`` at the level's own quadrature error."""
    level = 2
    verts, epsilon = _sphere(level)
    device = str(_CUDA_DEVICE)
    settings = ExteriorStokesMediumSettings(viscosity_pa_s=1e-3, blob_epsilon_um=epsilon)
    medium = build_exterior_stokes_medium(
        settings=settings,
        surface_index=np.arange(verts.shape[0]),
        initial_position=verts,
        device=device,
    )
    dt_phys = 1.0e-4
    omega = np.array([0.0, 0.0, 1.0])
    medium.set_step(dt_phys)
    rotated = verts + dt_phys * np.cross(omega, verts)
    pos_d = wp.array(rotated, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(verts.shape[0], dtype=wp.vec3d, device=device)
    medium.snapshot_candidate()
    medium.accumulate(pos_d, force_d)
    torque = np.cross(verts, force_d.numpy()).sum(axis=0)
    exact = sphere_rotation_resistance(_RADIUS_UM, settings.viscosity_pa_s)
    study = sphere_drag_convergence(
        radius=_RADIUS_UM,
        viscosity=settings.viscosity_pa_s,
        subdivisions=(1, level),
        epsilon_ratio=1.0,
    )
    quadrature_error = study["rows"][-1]["rotation_relative_error"]
    assert torque[2] < 0.0, "the medium must oppose the rotation"
    assert abs(abs(torque[2]) - exact) / exact == pytest.approx(quadrature_error, rel=1e-3)


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="I0-A: kernel gates require a CUDA GPU")
def test_accepted_step_banks_negative_work_and_a_rejected_step_changes_nothing() -> None:
    """The transaction contract: work enters only on acceptance, and rejection restores bit-for-bit.

    The dissipated work must be strictly negative — the medium removes energy — and a rejected candidate
    must leave the committed reference untouched so the re-attempted step sees an identical velocity.
    """
    verts, epsilon = _sphere(1)
    device = str(_CUDA_DEVICE)
    settings = ExteriorStokesMediumSettings(viscosity_pa_s=1e-3, blob_epsilon_um=epsilon)
    medium = build_exterior_stokes_medium(
        settings=settings,
        surface_index=np.arange(verts.shape[0]),
        initial_position=verts,
        device=device,
    )
    medium.set_step(0.01)
    moved = verts + np.array([0.02, 0.0, 0.0])
    pos_d = wp.array(moved, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(verts.shape[0], dtype=wp.vec3d, device=device)

    rejected = wp.zeros(1, dtype=wp.int32, device=device)
    reference_before = medium.committed_position_d.numpy().copy()
    medium.snapshot_candidate()
    medium.accumulate(pos_d, force_d)
    medium.commit_irreversible(rejected, dt_phys=0.01, rng_seed=0)
    medium.rollback(rejected)
    assert np.array_equal(medium.committed_position_d.numpy(), reference_before)
    assert float(medium.dissipated_work_d.numpy()[0]) == 0.0
    assert np.array_equal(medium.committed_traction_d.numpy(), np.zeros_like(verts))

    accepted = wp.ones(1, dtype=wp.int32, device=device)
    force_d.zero_()
    medium.snapshot_candidate()
    medium.accumulate(pos_d, force_d)
    medium.commit_irreversible(accepted, dt_phys=0.01, rng_seed=0)
    work = float(medium.dissipated_work_d.numpy()[0])
    assert work < 0.0, f"the medium must dissipate, banked work {work:.3e}"
    assert np.allclose(medium.committed_position_d.numpy(), moved, rtol=0.0, atol=0.0)


# --- the 36th connector's identity (CUDA-free binding gates) --------------------------------------


@dataclasses.dataclass(slots=True)
class _RecordingMedium:
    """Stands in for ``ExteriorStokesMedium``; records which hooks the connector facade forwarded."""

    calls: list[str] = dataclasses.field(default_factory=list)

    def accumulate(self, pos: object, force: object) -> None:
        self.calls.append("accumulate")

    def snapshot_candidate(self) -> None:
        self.calls.append("snapshot_candidate")

    def rollback(self, accepted: object) -> None:
        self.calls.append("rollback")

    def commit_irreversible(self, accepted: object, dt_phys: float, rng_seed: int) -> None:
        self.calls.append("commit_irreversible")

    def accumulate_ledger(self, ledger: object) -> None:
        self.calls.append("accumulate_ledger")


def test_the_traction_connector_takes_its_identity_from_the_declared_contract() -> None:
    """The 36th connector had physics and no identity; this is the identity."""
    from aleph.engine.contracts import reference_cell_architecture
    from aleph.engine.medium_exterior import MembraneMediumTraction

    contract = next(
        c for c in reference_cell_architecture().connectors if c.name == "membrane_medium_traction"
    )
    connector = MembraneMediumTraction(delegate=_RecordingMedium())
    assert connector.name == contract.name
    assert {connector.component_a, connector.component_b} == {
        contract.component_a, contract.component_b
    }
    assert connector.bidirectional is contract.bidirectional
    assert connector.adjoint_transfer_required is contract.adjoint_transfer_required


def test_a_non_moving_far_field_still_carries_an_adjoint_reaction() -> None:
    """``dynamically_evolving=False`` governs POSITION; ``adjoint_transfer_required`` governs FORCE.

    A clamp would discard the reaction and the balance ledger's two channels would then agree by
    OMISSION — the failure a one-sided scatter exists to be caught by (PI, 2026-08-09).
    """
    from aleph.engine.contracts import reference_cell_architecture
    from aleph.engine.medium_exterior import MembraneMediumTraction

    arch = reference_cell_architecture()
    medium = arch.component(EXTERIOR_MEDIUM_COMPONENT)
    assert medium.dynamically_evolving is False, "the far-field frame does not move"

    delegate = _RecordingMedium()
    connector = MembraneMediumTraction(delegate=delegate)
    assert connector.adjoint_transfer_required is True, "and its reaction is still accumulated"
    connector.accumulate_ledger(object())
    assert delegate.calls == ["accumulate_ledger"], "the reaction reaches the reservoir, not /dev/null"


def test_every_hook_reaches_the_exterior_solve() -> None:
    from aleph.engine.medium_exterior import MembraneMediumTraction

    delegate = _RecordingMedium()
    connector = MembraneMediumTraction(delegate=delegate)
    connector.snapshot_candidate()
    connector.accumulate_boundary(object(), object())
    connector.accumulate(object(), object())
    connector.commit_irreversible(object(), dt_phys=1.0e-3, rng_seed=0)
    connector.rollback(object())
    connector.accumulate_ledger(object())
    assert delegate.calls == [
        "snapshot_candidate", "accumulate", "accumulate",
        "commit_irreversible", "rollback", "accumulate_ledger",
    ]


def test_a_delegate_without_the_exterior_api_is_refused_at_construction() -> None:
    from aleph.engine.medium_exterior import MembraneMediumTraction

    with pytest.raises(TypeError, match="incomplete exterior-medium API"):
        MembraneMediumTraction(delegate=object())


def test_mis_stated_endpoints_are_refused() -> None:
    from aleph.engine.medium_exterior import MembraneMediumTraction

    with pytest.raises(ValueError, match="endpoints must be membrane and extracellular_medium"):
        MembraneMediumTraction(delegate=_RecordingMedium(), component_a="cortex")
