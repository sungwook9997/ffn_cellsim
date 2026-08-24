"""Pixel-integrated optics, extended emitters, z projection, and the evidence provenance DAG.

The optics assertions are written as MEASUREMENTS with stated tolerances rather than as smoke tests:
the point of a pixel-integrated PSF is a number (photon conservation under sub-pixel translation) that
the incumbent point-sampled renderer does not achieve, so the test has to show both sides of it.

The registry assertions centre on one thing: the master plan section 10.2 A6 failure, where synthetic
renderer output launders itself into the biological-truth corpus one derivation at a time.  The
transitive-exclusion test is the most important test in this file.
"""

from __future__ import annotations

import math
import subprocess
import sys

import numpy as np
import pytest
from scipy.special import erfc

from aleph.virtual_cell.optics import (
    UNIMPLEMENTED_AXIAL_MODELS,
    AxialModel,
    LineEmitter,
    OpticsConfig,
    PointEmitter,
    PSFSampling,
    TriangleEmitter,
    TruncationPolicy,
    defocused_sigma_um,
    gaussian_line_field,
    gaussian_pixel_weights_1d,
    gaussian_triangle_mass,
    psf_truncation_tail_fraction,
    render_emitters,
    render_widefield_stack,
)
from aleph.virtual_cell.source_registry import (
    BLOCKED_PUBLIC_MICROSCOPY_DATASETS,
    DECLARED_EVIDENCE_GAPS,
    REQUIRED_DATASET_REQUIREMENTS,
    DatasetIntegrationStatus,
    EvidenceGapVerdict,
    ExclusionReason,
    HeldOutValidationRefused,
    MicroscopyDatasetContract,
    ProvenanceKind,
    SourceCycleError,
    SourceOrigin,
    SourceRecord,
    SourceRegistry,
    SourceStatus,
    UnknownParentError,
    register_public_dataset_contracts,
    worked_fixture_dataset,
)
from aleph.virtual_cell.synthetic_microscopy import MicroscopyConfig, render_point_channel

# --- helpers ---------------------------------------------------------------------------------------


def _microscope(
    *,
    sigma_um: float,
    pixel_um: float,
    size_px: int = 96,
    radius_sigma: float = 8.0,
) -> MicroscopyConfig:
    """A microscope with a DELIBERATELY wide truncation radius; see ``_default_microscope``.

    8 sigma is not the renderer's default (4 is).  It is used by the photon-conservation and
    aliasing measurements in this file, which need the analytic truncation loss pushed below the
    round-off they assert at.  It must never be used to state a property of the SHIPPED renderer:
    the subdivision-invariance tests below did exactly that and hid a six-order-of-magnitude
    regression at the default radius.
    """
    return MicroscopyConfig(
        height_px=size_px,
        width_px=size_px,
        pixel_size_um=pixel_um,
        psf_sigma_um=sigma_um,
        background=0.0,
        psf_radius_sigma=radius_sigma,
    )


def _default_microscope(
    *,
    sigma_um: float,
    pixel_um: float,
    size_px: int = 96,
) -> MicroscopyConfig:
    """A microscope at the renderer's own defaults: ``psf_radius_sigma`` is not passed at all."""
    return MicroscopyConfig(
        height_px=size_px,
        width_px=size_px,
        pixel_size_um=pixel_um,
        psf_sigma_um=sigma_um,
        background=0.0,
    )


def _subdivision_error(whole: np.ndarray, split: np.ndarray) -> tuple[float, float]:
    """Return ``(max absolute, peak-relative)`` difference between two renders of one emitter."""
    absolute = float(np.max(np.abs(split - whole)))
    return absolute, absolute / float(np.max(whole))


def _optics(
    microscope: MicroscopyConfig,
    *,
    sampling: PSFSampling = PSFSampling.PIXEL_INTEGRATED,
    truncation: TruncationPolicy = TruncationPolicy.ANALYTIC,
    nodes: int = 4,
    exposure_s: float = 1.0,
    rayleigh_range_um: float = 0.5,
    focal_plane_z_um: float = 0.0,
) -> OpticsConfig:
    return OpticsConfig(
        microscopy=microscope,
        sampling=sampling,
        truncation=truncation,
        rayleigh_range_um=rayleigh_range_um,
        focal_plane_z_um=focal_plane_z_um,
        exposure_s=exposure_s,
        pixel_quadrature_nodes=nodes,
    )


def _subpixel_totals(
    microscope: MicroscopyConfig, sampling: PSFSampling, *, steps: int = 11
) -> np.ndarray:
    optics = _optics(microscope, sampling=sampling)
    pixel = microscope.pixel_size_um
    return np.array(
        [
            float(
                np.sum(
                    render_emitters(
                        [PointEmitter((12.0 + fraction * pixel, 12.0), 1.0)],
                        optics,
                    )
                )
            )
            for fraction in np.linspace(0.0, 1.0, steps)
        ]
    )


def _second_moment_px(image: np.ndarray) -> float:
    """Isotropic second moment about the intensity centroid, in pixels squared."""
    rows = np.arange(image.shape[0], dtype=np.float64)[:, None]
    columns = np.arange(image.shape[1], dtype=np.float64)[None, :]
    total = float(np.sum(image))
    row_mean = float(np.sum(image * rows) / total)
    column_mean = float(np.sum(image * columns) / total)
    spread = np.square(rows - row_mean) + np.square(columns - column_mean)
    return float(np.sum(image * spread) / total)


# --- pixel-integrated PSF: parity with the incumbent -------------------------------------------------


def test_point_sampled_renormalized_reproduces_the_incumbent_point_renderer() -> None:
    microscope = MicroscopyConfig(height_px=32, width_px=32, pixel_size_um=0.25, psf_sigma_um=0.4)
    incumbent = render_point_channel([[4.0, 4.0], [2.3, 5.7]], [2.0, 1.5], microscope)
    rebuilt = render_emitters(
        [PointEmitter((4.0, 4.0), 2.0), PointEmitter((2.3, 5.7), 1.5)],
        _optics(
            microscope,
            sampling=PSFSampling.POINT_SAMPLED,
            truncation=TruncationPolicy.RENORMALIZED,
        ),
    )
    np.testing.assert_allclose(rebuilt, incumbent, rtol=1e-12, atol=1e-15)


def test_frame_edge_still_loses_photons_instead_of_creating_a_boundary_ghost() -> None:
    microscope = MicroscopyConfig(height_px=32, width_px=32, pixel_size_um=0.25, psf_sigma_um=0.4)
    optics = _optics(microscope)
    edge = render_emitters([PointEmitter((0.05, 4.0), 1.0)], optics)
    outside = render_emitters([PointEmitter((-1.75, 4.0), 1.0)], optics)
    assert 0.0 < float(np.sum(edge)) < 1.0
    assert 0.0 < float(np.sum(outside)) < 1e-4
    assert render_emitters([PointEmitter((1_000.0, 1_000.0), 2.0)], optics).sum() == 0.0


def test_world_y_up_still_maps_larger_y_to_a_smaller_image_row() -> None:
    microscope = _microscope(sigma_um=0.4, pixel_um=0.25, size_px=64)
    optics = _optics(microscope)
    rows = np.arange(microscope.height_px, dtype=np.float64)[:, None]
    lower = render_emitters([PointEmitter((4.0, 3.0), 1.0)], optics)
    upper = render_emitters([PointEmitter((4.0, 4.0), 1.0)], optics)
    lower_row = float(np.sum(lower * rows) / np.sum(lower))
    upper_row = float(np.sum(upper * rows) / np.sum(upper))
    assert lower_row - upper_row == pytest.approx(1.0 / microscope.pixel_size_um, abs=1e-9)


# --- pixel-integrated PSF: the measurements that justify it ------------------------------------------


def test_pixel_weights_are_the_closed_form_erf_difference_and_sum_to_one() -> None:
    weights = gaussian_pixel_weights_1d(-40, 81, 0.37, 2.5)
    assert weights.shape == (81,)
    expected_first = 0.5 * (
        math.erf((-39.5 - 0.37) / (2.5 * math.sqrt(2)))
        - math.erf((-40.5 - 0.37) / (2.5 * math.sqrt(2)))
    )
    assert float(weights[0]) == pytest.approx(expected_first, rel=1e-15, abs=1e-18)
    assert float(np.sum(weights)) == pytest.approx(1.0, abs=1e-15)


def test_integrated_and_point_sampled_agree_when_the_psf_is_much_wider_than_a_pixel() -> None:
    """sigma = 4 px: total agrees to round-off, per-pixel disagreement stays under 1 percent."""
    microscope = _microscope(sigma_um=1.0, pixel_um=0.25, size_px=160)
    emitters = [PointEmitter((20.0, 20.0), 1.0)]
    integrated = render_emitters(
        emitters, _optics(microscope, sampling=PSFSampling.PIXEL_INTEGRATED)
    )
    sampled = render_emitters(emitters, _optics(microscope, sampling=PSFSampling.POINT_SAMPLED))
    assert float(np.sum(integrated)) == pytest.approx(1.0, abs=1e-14)
    assert float(np.sum(sampled)) == pytest.approx(1.0, abs=1e-14)
    relative_peak_difference = float(np.max(np.abs(integrated - sampled)) / np.max(integrated))
    # MEASURED 5.13e-3; the midpoint rule behaves like a Gaussian of variance sigma^2 + 1/12, which
    # predicts a peak deficit of 1 - sigma^2/(sigma^2 + 1/12) = 5.2e-3 at sigma = 4 px.
    assert 1e-3 < relative_peak_difference < 1e-2


def test_narrow_psf_breaks_point_sampling_and_not_the_integrated_form() -> None:
    """sigma = 0.3 px: point sampling loses and gains photons with sub-pixel position."""
    microscope = _microscope(sigma_um=0.075, pixel_um=0.25)
    integrated = _subpixel_totals(microscope, PSFSampling.PIXEL_INTEGRATED)
    sampled = _subpixel_totals(microscope, PSFSampling.POINT_SAMPLED)
    # MEASURED: integrated peak-to-peak 1.1e-16, point-sampled 0.4489 on a declared total of 1.0.
    assert float(np.ptp(integrated)) < 1e-12
    np.testing.assert_allclose(integrated, 1.0, atol=1e-12)
    assert float(np.ptp(sampled)) > 0.4
    assert float(np.min(sampled)) < 0.5
    assert float(np.max(sampled)) < 0.9


def test_sub_pixel_translation_does_not_change_the_integrated_total_at_any_width() -> None:
    """The sharpest test: total signal must be translation invariant to round-off."""
    for sigma_px in (0.3, 0.5, 1.0, 2.0):
        microscope = _microscope(sigma_um=sigma_px * 0.25, pixel_um=0.25)
        totals = _subpixel_totals(microscope, PSFSampling.PIXEL_INTEGRATED, steps=17)
        assert float(np.ptp(totals)) < 1e-12, f"sigma_px={sigma_px}"


def test_point_sampling_error_shrinks_as_the_psf_outgrows_the_pixel() -> None:
    """MEASURED aliasing amplitude: 4.5e-1 at 0.3 px, 2.8e-2 at 0.5 px, 1.1e-8 at 1 px."""
    swings = {}
    for sigma_px in (0.3, 0.5, 1.0):
        microscope = _microscope(sigma_um=sigma_px * 0.25, pixel_um=0.25)
        swings[sigma_px] = float(np.ptp(_subpixel_totals(microscope, PSFSampling.POINT_SAMPLED)))
    assert swings[0.3] > 0.4
    assert 1e-2 < swings[0.5] < 1e-1
    assert swings[1.0] < 1e-6
    assert swings[0.3] > swings[0.5] > swings[1.0]


def test_analytic_truncation_loss_stays_under_the_declared_radius_bound() -> None:
    radius_sigma = 4.0
    microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=128, radius_sigma=radius_sigma)
    total = float(np.sum(render_emitters([PointEmitter((12.0, 12.0), 1.0)], _optics(microscope))))
    bound = 2.0 * float(erfc(radius_sigma / math.sqrt(2.0)))
    # MEASURED loss 1.36e-5 against a 1.27e-4 bound; the square window is wider than R sigma.
    assert 0.0 < 1.0 - total < bound


def test_renormalized_truncation_conserves_the_declared_photon_count_exactly() -> None:
    microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=128, radius_sigma=4.0)
    total = float(
        np.sum(
            render_emitters(
                [PointEmitter((12.0, 12.0), 7.5)],
                _optics(microscope, truncation=TruncationPolicy.RENORMALIZED),
            )
        )
    )
    assert total == pytest.approx(7.5, rel=1e-14)


def test_rendered_total_is_photons_and_scales_with_exposure() -> None:
    microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=128)
    one_second = render_emitters([PointEmitter((12.0, 12.0), 4.0)], _optics(microscope))
    three_seconds = render_emitters(
        [PointEmitter((12.0, 12.0), 4.0)], _optics(microscope, exposure_s=3.0)
    )
    assert float(np.sum(one_second)) == pytest.approx(4.0, abs=1e-12)
    assert float(np.sum(three_seconds)) == pytest.approx(12.0, abs=1e-12)


# --- extended emitters -------------------------------------------------------------------------------


def test_line_field_matches_the_closed_form_and_telescopes_under_splitting() -> None:
    start, end = np.array([2.0, 3.0]), np.array([11.0, 8.0])
    points = np.array([[3.0, 4.0], [7.0, 6.0], [12.0, 9.0], [0.0, 0.0]])
    whole = gaussian_line_field(start, end, 1.7, points)
    for pieces in (2, 3, 7):
        fractions = np.linspace(0.0, 1.0, pieces + 1)
        summed = np.zeros_like(whole)
        for index in range(pieces):
            summed += gaussian_line_field(
                start + fractions[index] * (end - start),
                start + fractions[index + 1] * (end - start),
                1.7,
                points,
            )
        np.testing.assert_allclose(summed, whole, rtol=0.0, atol=1e-15)


def test_line_emitter_total_is_rate_times_length_times_exposure() -> None:
    microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=128)
    emitter = LineEmitter((6.0, 6.0), (14.0, 11.0), 3.0)
    assert emitter.length_um == pytest.approx(math.hypot(8.0, 5.0))
    total = float(np.sum(render_emitters([emitter], _optics(microscope))))
    assert total == pytest.approx(3.0 * emitter.length_um, rel=1e-12)


def test_line_emitter_signal_scales_with_length_not_with_primitive_count() -> None:
    microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=128)
    optics = _optics(microscope)
    short = render_emitters([LineEmitter((6.0, 6.0), (10.0, 8.5), 3.0)], optics)
    long = render_emitters([LineEmitter((6.0, 6.0), (14.0, 11.0), 3.0)], optics)
    assert float(np.sum(long)) / float(np.sum(short)) == pytest.approx(2.0, rel=1e-12)


def _line_renders(microscope: MicroscopyConfig, pieces: int) -> tuple[np.ndarray, np.ndarray]:
    """Render one line emitter whole, and as ``pieces`` collinear sub-segments."""
    optics = _optics(microscope)
    start, end = (6.0, 6.0), (14.0, 11.0)
    whole = render_emitters([LineEmitter(start, end, 3.0)], optics)
    fractions = np.linspace(0.0, 1.0, pieces + 1)
    split = render_emitters(
        [
            LineEmitter(
                (
                    start[0] + fractions[index] * (end[0] - start[0]),
                    start[1] + fractions[index] * (end[1] - start[1]),
                ),
                (
                    start[0] + fractions[index + 1] * (end[0] - start[0]),
                    start[1] + fractions[index + 1] * (end[1] - start[1]),
                ),
                3.0,
            )
            for index in range(pieces)
        ],
        optics,
    )
    return whole, split


def test_the_shipped_truncation_radius_is_four_sigma_and_the_wide_helper_is_not_the_default() -> (
    None
):
    """The premise of every measurement below, asserted rather than assumed.

    The subdivision tests used to run through ``_microscope``, whose 8 sigma is not what
    ``MicroscopyConfig`` ships.  That single substitution turned a 7.51e-6 peak-relative error into
    3.16e-15 and let the file claim round-off invariance for a renderer that does not have it.
    """
    shipped = MicroscopyConfig(height_px=8, width_px=8, pixel_size_um=0.2, psf_sigma_um=0.4)
    assert shipped.psf_radius_sigma == 4.0
    assert _default_microscope(sigma_um=0.4, pixel_um=0.2).psf_radius_sigma == 4.0
    assert _microscope(sigma_um=0.4, pixel_um=0.2).psf_radius_sigma == 8.0


@pytest.mark.parametrize("pieces", [2, 3, 5, 17])
def test_subdividing_a_line_emitter_at_the_shipped_default_stays_under_the_truncated_tail(
    pieces: int,
) -> None:
    """The DEFAULT configuration, which is where the claim has to hold.

    MEASURED at ``psf_radius_sigma = 4``, ``sigma = 2 px``: max absolute difference 8.94e-7 (2
    pieces), 6.05e-7 (3), 2.75e-7 (5), 8.90e-7 (17); peak-relative 7.55e-6, 5.11e-6, 2.32e-6,
    7.51e-6.  These are NOT round-off — each piece clips a different PSF tail — and they sit below
    the truncated-tail scale ``exp(-4^2/2) = 3.35e-4`` the module declares.

    The FRAME TOTAL moves too, by 3.61e-6 relative for the 17-piece split: the pieces' windows do
    not tile the whole emitter's window, so the split render truncates slightly more photons.  The
    superseded test asserted that total to ``rel=1e-13``, which is true only at the 8 sigma its
    helper forced.  Both discrepancies are the same truncated tail and are bounded by the same
    declared scale.
    """
    tail = psf_truncation_tail_fraction(4.0)
    microscope = _default_microscope(sigma_um=0.4, pixel_um=0.2, size_px=128)
    assert microscope.psf_radius_sigma == 4.0
    whole, split = _line_renders(microscope, pieces)
    _, peak_relative = _subdivision_error(whole, split)
    assert peak_relative < tail
    whole_total, split_total = float(np.sum(whole)), float(np.sum(split))
    assert abs(split_total - whole_total) / whole_total < tail


@pytest.mark.parametrize("radius_sigma", [5.0, 6.0, 7.0, 8.0])
def test_widening_the_truncation_radius_buys_back_line_subdivision_invariance(
    radius_sigma: float,
) -> None:
    """The scope of the claim, measured across the radius rather than asserted at one value.

    MEASURED peak-relative difference for the 17-piece split: 7.51e-6 at R = 4 (the default),
    4.25e-8 at 5, 9.91e-11 at 6, 8.78e-14 at 7, 3.16e-15 at 8.  Every value is below
    ``exp(-R^2/2)``, and by R = 7 the discarded tail has fallen under double precision, so what
    remains is round-off and the original round-off claim becomes true.
    """
    microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=128, radius_sigma=radius_sigma)
    whole, split = _line_renders(microscope, 17)
    absolute, peak_relative = _subdivision_error(whole, split)
    assert peak_relative < psf_truncation_tail_fraction(radius_sigma)
    if radius_sigma >= 7.0:
        assert absolute < 1e-13


def test_triangle_mass_is_a_probability_and_is_orientation_independent() -> None:
    vertices = np.array([[0.0, 0.0], [6.0, 1.0], [2.0, 5.0]])
    points = np.array([[1.0, 1.0], [2.0, 2.0], [-8.0, -8.0], [3.0, 1.5]])
    mass = gaussian_triangle_mass(vertices, 1.7, points)
    assert np.all(mass >= 0.0)
    assert np.all(mass <= 1.0)
    reversed_mass = gaussian_triangle_mass(vertices[::-1], 1.7, points)
    np.testing.assert_allclose(reversed_mass, mass, rtol=0.0, atol=1e-15)
    huge = np.array([[-1e4, -1e4], [1e4, -1e4], [0.0, 1e4]])
    assert float(gaussian_triangle_mass(huge, 1.7, np.array([[0.0, 0.0]]))[0]) == pytest.approx(
        1.0, abs=1e-12
    )
    assert float(mass[2]) < 1e-6


def test_triangle_emitter_total_is_rate_times_area_times_exposure() -> None:
    microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=128)
    emitter = TriangleEmitter(((7.0, 7.0), (16.0, 8.5), (10.0, 15.0)), 2.0)
    assert emitter.area_um2 == pytest.approx(33.75)
    total = float(np.sum(render_emitters([emitter], _optics(microscope))))
    assert total == pytest.approx(2.0 * emitter.area_um2, rel=1e-12)


def test_triangle_emitter_signal_scales_with_area() -> None:
    microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=160)
    optics = _optics(microscope)
    small = TriangleEmitter(((8.0, 8.0), (12.0, 8.0), (8.0, 11.0)), 2.0)
    doubled = TriangleEmitter(((8.0, 8.0), (16.0, 8.0), (8.0, 11.0)), 2.0)
    assert doubled.area_um2 == pytest.approx(2.0 * small.area_um2)
    ratio = float(np.sum(render_emitters([doubled], optics))) / float(
        np.sum(render_emitters([small], optics))
    )
    assert ratio == pytest.approx(2.0, rel=1e-12)


def _midpoint_split(triangle: np.ndarray) -> list[np.ndarray]:
    """Split a triangle into four orientation-preserving sub-triangles at its edge midpoints."""
    first, second, third = triangle
    return [
        np.array([first, (first + second) / 2, (third + first) / 2]),
        np.array([(first + second) / 2, second, (second + third) / 2]),
        np.array([(third + first) / 2, (second + third) / 2, third]),
        np.array([(first + second) / 2, (second + third) / 2, (third + first) / 2]),
    ]


def _triangle_renders(microscope: MicroscopyConfig, levels: int) -> tuple[np.ndarray, np.ndarray]:
    """Render one triangle emitter whole, and as ``4 ** levels`` sub-triangles."""
    optics = _optics(microscope)
    vertices = ((7.0, 7.0), (16.0, 8.5), (10.0, 15.0))
    whole = render_emitters([TriangleEmitter(vertices, 2.0)], optics)
    triangles = [np.array(vertices, dtype=np.float64)]
    for _ in range(levels):
        triangles = [piece for triangle in triangles for piece in _midpoint_split(triangle)]
    split = render_emitters(
        [
            TriangleEmitter(tuple(tuple(float(value) for value in v) for v in triangle), 2.0)
            for triangle in triangles
        ],
        optics,
    )
    return whole, split


@pytest.mark.parametrize("levels", [1, 2, 3])
def test_subdividing_a_triangle_emitter_at_the_shipped_default_stays_under_the_truncated_tail(
    levels: int,
) -> None:
    """The DEFAULT configuration again; the triangle side of the same limit.

    MEASURED at ``psf_radius_sigma = 4``: max absolute difference 9.17e-8 (4 sub-triangles),
    1.02e-7 (16), 1.20e-7 (64); peak-relative 1.15e-6, 1.27e-6, 1.50e-6 — below the declared
    ``exp(-4^2/2) = 3.35e-4`` tail scale and four to five orders above round-off.  The frame total
    moves by 4.87e-7 relative at 64 sub-triangles, for the same reason and under the same bound.
    """
    tail = psf_truncation_tail_fraction(4.0)
    microscope = _default_microscope(sigma_um=0.4, pixel_um=0.2, size_px=128)
    assert microscope.psf_radius_sigma == 4.0
    whole, split = _triangle_renders(microscope, levels)
    _, peak_relative = _subdivision_error(whole, split)
    assert peak_relative < tail
    whole_total, split_total = float(np.sum(whole)), float(np.sum(split))
    assert abs(split_total - whole_total) / whole_total < tail


@pytest.mark.parametrize("radius_sigma", [5.0, 6.0, 7.0, 8.0])
def test_widening_the_truncation_radius_buys_back_triangle_subdivision_invariance(
    radius_sigma: float,
) -> None:
    """MEASURED peak-relative difference for the 64-triangle split.

    1.50e-6 at R = 4 (the default), 9.11e-9 at 5, 1.69e-11 at 6, 1.28e-14 at 7, 8.67e-16 at 8.
    """
    microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=128, radius_sigma=radius_sigma)
    whole, split = _triangle_renders(microscope, 3)
    absolute, peak_relative = _subdivision_error(whole, split)
    assert peak_relative < psf_truncation_tail_fraction(radius_sigma)
    if radius_sigma >= 7.0:
        assert absolute < 1e-13


def test_the_subdivision_error_tracks_the_truncated_tail_rather_than_the_piece_count() -> None:
    """Attribution, not just magnitude: the error is the discarded tail and nothing else.

    If the discrepancy were a quadrature or normalisation defect it would follow the number of
    pieces.  It does not: 17 pieces and 64 sub-triangles both fall by roughly the ratio of
    ``exp(-R^2/2)`` between successive radii, and both reach the double-precision floor at the same
    radius.  MEASURED line ratios between consecutive radii 4 -> 5 -> 6: 177x, 429x, against tail
    ratios 90x and 245x.
    """
    radii = (4.0, 5.0, 6.0)
    line_errors = []
    triangle_errors = []
    for radius in radii:
        microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=128, radius_sigma=radius)
        line_errors.append(_subdivision_error(*_line_renders(microscope, 17))[1])
        triangle_errors.append(_subdivision_error(*_triangle_renders(microscope, 3))[1])
    for errors in (line_errors, triangle_errors):
        assert errors[0] > errors[1] > errors[2]
        for error, radius in zip(errors, radii, strict=True):
            assert error < psf_truncation_tail_fraction(radius)
        # The fall must be at least as steep as the tail's own fall, which is what identifies it.
        for index in range(len(radii) - 1):
            tail_ratio = psf_truncation_tail_fraction(radii[index]) / psf_truncation_tail_fraction(
                radii[index + 1]
            )
            assert errors[index] / errors[index + 1] > tail_ratio


def test_the_truncated_tail_fraction_is_the_gaussian_value_at_the_radius() -> None:
    assert psf_truncation_tail_fraction(4.0) == pytest.approx(math.exp(-8.0), rel=1e-15)
    assert psf_truncation_tail_fraction(1.0) == pytest.approx(0.6065306597126334, rel=1e-15)
    assert 0.0 < psf_truncation_tail_fraction(8.0) < psf_truncation_tail_fraction(4.0) < 1.0
    for bad in (0.0, -1.0, float("nan"), float("inf"), True, "4"):
        with pytest.raises(ValueError):
            psf_truncation_tail_fraction(bad)  # type: ignore[arg-type]


def test_extended_emitter_pixel_quadrature_converges_while_the_total_is_exact_throughout() -> None:
    """The Gauss-Legendre pixel rule sets per-pixel accuracy; it never moves the photon budget."""
    microscope = _microscope(sigma_um=0.4, pixel_um=0.2, size_px=128)
    line = LineEmitter((6.0, 6.0), (14.0, 11.0), 3.0)
    triangle = TriangleEmitter(((7.0, 7.0), (16.0, 8.5), (10.0, 15.0)), 2.0)
    for emitter, expected in (
        (line, 3.0 * line.length_um),
        (triangle, 2.0 * triangle.area_um2),
    ):
        reference = render_emitters([emitter], _optics(microscope, nodes=16))
        errors = []
        for nodes in (1, 2, 3, 4):
            image = render_emitters([emitter], _optics(microscope, nodes=nodes))
            errors.append(float(np.max(np.abs(image - reference)) / np.max(reference)))
            assert float(np.sum(image)) == pytest.approx(expected, rel=1e-12)
        # MEASURED for the line at sigma = 2 px: 1.17e-2, 3.22e-5, 6.78e-8, 1.06e-10.
        assert errors == sorted(errors, reverse=True)
        assert errors[0] > 1e-3
        assert errors[3] < 1e-8


def test_extended_emitters_reject_degenerate_geometry() -> None:
    with pytest.raises(ValueError, match="distinct"):
        LineEmitter((1.0, 1.0), (1.0, 1.0), 1.0)
    with pytest.raises(ValueError, match="positive area"):
        TriangleEmitter(((0.0, 0.0), (1.0, 1.0), (2.0, 2.0)), 1.0)
    with pytest.raises(ValueError, match="three"):
        TriangleEmitter(((0.0, 0.0), (1.0, 1.0)), 1.0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="nonnegative"):
        PointEmitter((1.0, 1.0), -1.0)


# --- z projection ------------------------------------------------------------------------------------


def test_defocus_widens_the_lateral_psf_by_the_declared_gaussian_beam_law() -> None:
    optics = _optics(_microscope(sigma_um=0.2, pixel_um=0.1), rayleigh_range_um=0.5)
    assert defocused_sigma_um(0.0, optics) == pytest.approx(0.2)
    assert defocused_sigma_um(0.5, optics) == pytest.approx(0.2 * math.sqrt(2.0))
    assert defocused_sigma_um(-0.5, optics) == pytest.approx(0.2 * math.sqrt(2.0))
    assert defocused_sigma_um(1.0, optics) == pytest.approx(0.2 * math.sqrt(5.0))


def test_an_in_focus_emitter_is_sharper_than_a_defocused_one() -> None:
    microscope = _microscope(sigma_um=0.2, pixel_um=0.1, size_px=192)
    optics = _optics(microscope)
    focused = render_emitters([PointEmitter((9.6, 9.6), 1.0, z_um=0.0)], optics)
    defocused = render_emitters([PointEmitter((9.6, 9.6), 1.0, z_um=1.0)], optics)
    assert float(np.max(focused)) > float(np.max(defocused))
    assert _second_moment_px(focused) < _second_moment_px(defocused)
    # A pixel-INTEGRATED kernel is the Gaussian convolved with the pixel box, so the recovered
    # variance per axis is sigma^2 + 1/12 exactly.  Point sampling would report sigma^2 and thereby
    # under-report the width the camera actually records.
    assert _second_moment_px(focused) == pytest.approx(2.0 * (2.0**2 + 1.0 / 12.0), rel=1e-6)
    assert _second_moment_px(defocused) == pytest.approx(
        2.0 * (defocused_sigma_um(1.0, optics) / 0.1) ** 2 + 2.0 / 12.0, rel=1e-6
    )
    sampled = render_emitters(
        [PointEmitter((9.6, 9.6), 1.0, z_um=0.0)],
        _optics(microscope, sampling=PSFSampling.POINT_SAMPLED),
    )
    assert _second_moment_px(sampled) == pytest.approx(2.0 * 2.0**2, rel=1e-6)


def test_widefield_defocus_spreads_light_without_destroying_it() -> None:
    microscope = _microscope(sigma_um=0.2, pixel_um=0.1, size_px=192)
    optics = _optics(microscope)
    totals = []
    peaks = []
    widths = []
    for depth in (0.0, 0.25, 0.5, 1.0, 2.0):
        image = render_emitters([PointEmitter((9.6, 9.6), 1.0, z_um=depth)], optics)
        totals.append(float(np.sum(image)))
        peaks.append(float(np.max(image)))
        widths.append(defocused_sigma_um(depth, optics) / microscope.pixel_size_um)
    # MEASURED: frame total spread 8.9e-16 across a 4x lateral broadening; peak falls 15.7x.
    assert float(np.ptp(np.array(totals))) < 1e-12
    np.testing.assert_allclose(np.array(totals), 1.0, atol=1e-12)
    assert peaks == sorted(peaks, reverse=True)
    assert peaks[0] / peaks[-1] > 10.0
    scaled = np.array(peaks) * np.square(np.array(widths))
    assert float(np.ptp(scaled) / np.mean(scaled)) < 0.1


def test_widefield_stack_collects_out_of_focus_light_in_every_plane() -> None:
    microscope = _microscope(sigma_um=0.2, pixel_um=0.1, size_px=192)
    optics = _optics(microscope)
    planes = (-1.0, -0.5, 0.0, 0.5, 1.0)
    stack = render_widefield_stack([PointEmitter((9.6, 9.6), 1.0, z_um=0.5)], optics, planes)
    assert stack.shape == (len(planes), microscope.height_px, microscope.width_px)
    assert stack.flags.writeable is False
    per_plane_total = stack.sum(axis=(1, 2))
    np.testing.assert_allclose(per_plane_total, 1.0, atol=1e-12)
    sharpest = int(np.argmax(stack.max(axis=(1, 2))))
    assert planes[sharpest] == 0.5
    with pytest.raises(ValueError, match="at least one plane"):
        render_widefield_stack([PointEmitter((9.6, 9.6), 1.0)], optics, [])


def test_extended_emitters_defocus_too_and_keep_their_photon_budget() -> None:
    microscope = _microscope(sigma_um=0.2, pixel_um=0.1, size_px=256)
    optics = _optics(microscope)
    emitter = LineEmitter((8.0, 8.0), (14.0, 12.0), 5.0, z_um=0.75)
    focused = LineEmitter((8.0, 8.0), (14.0, 12.0), 5.0, z_um=0.0)
    blurred_image = render_emitters([emitter], optics)
    focused_image = render_emitters([focused], optics)
    expected = 5.0 * emitter.length_um
    assert float(np.sum(blurred_image)) == pytest.approx(expected, rel=1e-10)
    assert float(np.sum(focused_image)) == pytest.approx(expected, rel=1e-10)
    assert float(np.max(focused_image)) > float(np.max(blurred_image))


def test_confocal_and_light_sheet_are_refused_rather_than_approximated() -> None:
    microscope = _microscope(sigma_um=0.2, pixel_um=0.1)
    assert set(UNIMPLEMENTED_AXIAL_MODELS) == {AxialModel.CONFOCAL, AxialModel.LIGHT_SHEET}
    for model in (AxialModel.CONFOCAL, AxialModel.LIGHT_SHEET):
        verdict = UNIMPLEMENTED_AXIAL_MODELS[model]
        assert verdict.status == "blocked"
        assert verdict.routes_to
        with pytest.raises(NotImplementedError, match="blocked"):
            OpticsConfig(microscopy=microscope, axial_model=model)


# --- optics configuration ----------------------------------------------------------------------------


def test_optics_config_hash_binds_every_choice_that_changes_the_image() -> None:
    microscope = _microscope(sigma_um=0.2, pixel_um=0.1)
    baseline = _optics(microscope)
    assert len(baseline.config_hash) == 64
    assert baseline.config_hash == _optics(microscope).config_hash
    assert (
        baseline.config_hash != _optics(microscope, sampling=PSFSampling.POINT_SAMPLED).config_hash
    )
    assert (
        baseline.config_hash
        != _optics(microscope, truncation=TruncationPolicy.RENORMALIZED).config_hash
    )
    assert baseline.config_hash != _optics(microscope, exposure_s=2.0).config_hash
    assert baseline.renderer_id == "optics-pixel-integrated-analytic-widefield@1"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"rayleigh_range_um": 0.0}, "rayleigh_range_um"),
        ({"exposure_s": -1.0}, "exposure_s"),
        ({"pixel_quadrature_nodes": 0}, "pixel_quadrature_nodes"),
        ({"focal_plane_z_um": float("nan")}, "focal_plane_z_um"),
    ],
)
def test_optics_config_rejects_unphysical_settings(kwargs: dict[str, object], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        OpticsConfig(microscopy=_microscope(sigma_um=0.2, pixel_um=0.1), **kwargs)  # type: ignore[arg-type]


def test_render_emitters_rejects_the_wrong_shapes_of_input() -> None:
    optics = _optics(_microscope(sigma_um=0.2, pixel_um=0.1))
    with pytest.raises(TypeError, match="sequence of emitters"):
        render_emitters(PointEmitter((1.0, 1.0), 1.0), optics)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="unsupported emitter"):
        render_emitters([object()], optics)  # type: ignore[list-item]
    with pytest.raises(TypeError, match="OpticsConfig"):
        render_emitters([], _microscope(sigma_um=0.2, pixel_um=0.1))  # type: ignore[arg-type]


def test_empty_emitter_sequence_returns_the_background_frame() -> None:
    microscope = MicroscopyConfig(
        height_px=8, width_px=9, pixel_size_um=0.2, psf_sigma_um=0.3, background=1.5
    )
    image = render_emitters([], _optics(microscope))
    assert image.shape == (8, 9)
    assert np.all(image == 1.5)
    assert image.flags.writeable is False


# --- provenance DAG: cycles --------------------------------------------------------------------------


def _literature(source_id: str, **overrides: object) -> SourceRecord:
    payload: dict[str, object] = {
        "source_id": source_id,
        "provenance": ProvenanceKind.SOURCED,
        "origin": SourceOrigin.LITERATURE,
        "status": SourceStatus.ACTIVE,
        "description": f"published measurement {source_id}",
        "citation": f"doi:10.0000/{source_id}",
    }
    payload.update(overrides)
    return SourceRecord(**payload)  # type: ignore[arg-type]


def _derived(source_id: str, parents: tuple[str, ...], **overrides: object) -> SourceRecord:
    payload: dict[str, object] = {
        "source_id": source_id,
        "provenance": ProvenanceKind.DERIVED,
        "origin": SourceOrigin.LITERATURE,
        "status": SourceStatus.ACTIVE,
        "description": f"derived quantity {source_id}",
        "citation": f"derivation note {source_id}",
        "parents": parents,
    }
    payload.update(overrides)
    return SourceRecord(**payload)  # type: ignore[arg-type]


def test_a_self_loop_is_refused_and_reported_as_a_path() -> None:
    registry = SourceRegistry().register(_derived("self-citing", ("self-citing",)))
    with pytest.raises(SourceCycleError) as error:
        registry.freeze()
    assert error.value.cycle == ("self-citing", "self-citing")
    assert "self-citing -> self-citing" in str(error.value)


def test_a_two_node_cycle_is_refused_and_reported_as_a_path() -> None:
    registry = SourceRegistry().register_all(
        [_derived("alpha", ("beta",)), _derived("beta", ("alpha",))]
    )
    with pytest.raises(SourceCycleError) as error:
        registry.freeze()
    assert error.value.cycle[0] == error.value.cycle[-1] == "alpha"
    assert set(error.value.cycle) == {"alpha", "beta"}
    assert len(error.value.cycle) == 3


def test_a_long_cycle_is_refused_with_the_whole_offending_path() -> None:
    names = ("node-a", "node-b", "node-c", "node-d", "node-e")
    parent_of = {name: names[(index + 1) % len(names)] for index, name in enumerate(names)}
    registry = SourceRegistry().register_all([_derived(name, (parent_of[name],)) for name in names])
    with pytest.raises(SourceCycleError) as error:
        registry.freeze()
    cycle = error.value.cycle
    assert len(cycle) == len(names) + 1
    assert cycle[0] == cycle[-1] == "node-a"
    assert set(cycle) == set(names)
    # Every consecutive pair is a real child -> parent edge, so the report is a PATH, not a set.
    for child, parent in zip(cycle[:-1], cycle[1:], strict=True):
        assert parent_of[child] == parent
    assert " -> ".join(cycle) in str(error.value)


def test_a_dangling_parent_reference_is_refused_by_name() -> None:
    registry = SourceRegistry().register(_derived("child", ("never-registered",)))
    with pytest.raises(UnknownParentError) as error:
        registry.freeze()
    assert error.value.source_id == "child"
    assert error.value.parent_id == "never-registered"


def test_an_acyclic_diamond_is_accepted_and_topologically_ordered() -> None:
    graph = (
        SourceRegistry()
        .register_all(
            [
                _literature("root"),
                _derived("left", ("root",)),
                _derived("right", ("root",)),
                _derived("join", ("left", "right")),
            ]
        )
        .freeze()
    )
    order = graph.topological_order
    assert set(order) == {"root", "left", "right", "join"}
    assert order.index("root") < order.index("left") < order.index("join")
    assert order.index("root") < order.index("right") < order.index("join")
    assert graph.ancestors("join") == frozenset({"left", "right", "root"})
    assert graph.ancestors("root") == frozenset()
    assert len(graph) == 4
    assert "join" in graph


def test_duplicate_registration_is_refused() -> None:
    registry = SourceRegistry().register(_literature("only-once"))
    with pytest.raises(ValueError, match="already registered"):
        registry.register(_literature("only-once"))


# --- provenance DAG: the A6 exclusion ----------------------------------------------------------------


def _a6_graph():
    """Renderer output -> encoder corpus -> morphology statistic, plus a clean literature control."""
    return (
        SourceRegistry()
        .register_all(
            [
                SourceRecord(
                    source_id="native-state-trace",
                    provenance=ProvenanceKind.SOURCED,
                    origin=SourceOrigin.SIMULATION,
                    status=SourceStatus.ACTIVE,
                    description="accepted-step simulation state trace",
                ),
                SourceRecord(
                    source_id="synthetic-renderer-output",
                    provenance=ProvenanceKind.DERIVED,
                    origin=SourceOrigin.SYNTHETIC_RENDER,
                    status=SourceStatus.ACTIVE,
                    description="optics.render_emitters output for the state trace",
                    parents=("native-state-trace",),
                ),
                _derived("encoder-training-corpus", ("synthetic-renderer-output",)),
                _derived("morphology-statistic", ("encoder-training-corpus",)),
                _literature("published-cell-morphology"),
                _derived("literature-derived-shape-prior", ("published-cell-morphology",)),
            ]
        )
        .freeze()
    )


def test_a_source_two_hops_downstream_of_a_renderer_is_excluded_transitively() -> None:
    """Master plan 10.2 A6: the laundering happens one derivation at a time, so exclusion is closed."""
    graph = _a6_graph()
    exclusion = graph.exclusion("morphology-statistic")
    assert exclusion is not None
    assert exclusion.reason is ExclusionReason.DERIVED_FROM_EXCLUDED
    assert exclusion.blocking_source_id == "synthetic-renderer-output"
    assert exclusion.path == (
        "morphology-statistic",
        "encoder-training-corpus",
        "synthetic-renderer-output",
    )
    assert "morphology-statistic" not in graph.citable_as_biological_evidence()


def test_the_direct_child_of_a_renderer_is_excluded_as_well() -> None:
    graph = _a6_graph()
    direct = graph.exclusion("encoder-training-corpus")
    assert direct is not None
    assert direct.blocking_source_id == "synthetic-renderer-output"
    assert direct.path == ("encoder-training-corpus", "synthetic-renderer-output")


def test_the_renderer_and_the_simulation_are_excluded_on_their_own_account() -> None:
    graph = _a6_graph()
    for source_id in ("native-state-trace", "synthetic-renderer-output"):
        exclusion = graph.exclusion(source_id)
        assert exclusion is not None
        assert exclusion.reason is ExclusionReason.NON_EMPIRICAL_ORIGIN
        assert exclusion.blocking_source_id == source_id
        assert exclusion.path == (source_id,)


def test_only_the_empirical_chain_survives_the_biological_evidence_query() -> None:
    graph = _a6_graph()
    assert graph.citable_as_biological_evidence() == (
        "literature-derived-shape-prior",
        "published-cell-morphology",
    )
    assert set(graph.excluded_from_biological_evidence()) == {
        "native-state-trace",
        "synthetic-renderer-output",
        "encoder-training-corpus",
        "morphology-statistic",
    }


@pytest.mark.parametrize(
    ("status", "reason"),
    [
        (SourceStatus.REJECTED, ExclusionReason.REJECTED),
        (SourceStatus.RETRACTED, ExclusionReason.RETRACTED),
        (SourceStatus.BLOCKED, ExclusionReason.BLOCKED),
        (SourceStatus.FIXTURE, ExclusionReason.FIXTURE),
    ],
)
def test_bad_standing_propagates_to_every_descendant(
    status: SourceStatus, reason: ExclusionReason
) -> None:
    graph = (
        SourceRegistry()
        .register_all(
            [
                _literature("tainted-root", status=status),
                _derived("first-hop", ("tainted-root",)),
                _derived("second-hop", ("first-hop",)),
            ]
        )
        .freeze()
    )
    assert graph.exclusion("tainted-root").reason is reason  # type: ignore[union-attr]
    deep = graph.exclusion("second-hop")
    assert deep is not None
    assert deep.reason is ExclusionReason.DERIVED_FROM_EXCLUDED
    assert deep.blocking_source_id == "tainted-root"
    assert deep.path == ("second-hop", "first-hop", "tainted-root")
    assert graph.citable_as_biological_evidence() == ()


@pytest.mark.parametrize("provenance", [ProvenanceKind.CONVENIENCE, ProvenanceKind.PI_GAP])
def test_placeholder_provenance_is_not_biological_evidence(provenance: ProvenanceKind) -> None:
    graph = (
        SourceRegistry()
        .register(
            _literature("placeholder", provenance=provenance),
        )
        .freeze()
    )
    exclusion = graph.exclusion("placeholder")
    assert exclusion is not None
    assert exclusion.reason is ExclusionReason.PLACEHOLDER_PROVENANCE


def test_computed_qm_md_is_a_prior_and_not_biological_evidence() -> None:
    graph = (
        SourceRegistry()
        .register(
            SourceRecord(
                source_id="qm-md-binding-free-energy",
                provenance=ProvenanceKind.COMPUTED_QM_MD,
                origin=SourceOrigin.COMPUTATION,
                status=SourceStatus.ACTIVE,
                description="computed prior for a bond rate",
            )
        )
        .freeze()
    )
    exclusion = graph.exclusion("qm-md-binding-free-energy")
    assert exclusion is not None
    assert exclusion.reason is ExclusionReason.NON_EMPIRICAL_ORIGIN


def test_unknown_source_queries_raise_rather_than_answering_optimistically() -> None:
    graph = SourceRegistry().register(_literature("known")).freeze()
    for call in (graph.exclusion, graph.ancestors, graph.calibration_ancestry, graph.record):
        with pytest.raises(KeyError):
            call("absent")
    with pytest.raises(KeyError):
        graph.assert_admissible_as_held_out_validation("absent")


# --- provenance DAG: calibration and held-out validation ---------------------------------------------


def _calibration_graph():
    return (
        SourceRegistry()
        .register_all(
            [
                _literature(
                    "afm-calibration-set",
                    provenance=ProvenanceKind.CALIBRATION_ONLY,
                    held_out_split_id="split-cal",
                ),
                _derived(
                    "calibration-derived-split",
                    ("afm-calibration-set",),
                    held_out_split_id="split-a",
                ),
                _literature("independent-traction-study", held_out_split_id="split-b"),
                _literature("study-without-a-split"),
            ]
        )
        .freeze()
    )


def test_calibration_only_may_not_serve_as_held_out_validation() -> None:
    graph = _calibration_graph()
    with pytest.raises(HeldOutValidationRefused, match="CALIBRATION_ONLY"):
        graph.assert_admissible_as_held_out_validation("afm-calibration-set")


def test_calibration_leaks_transitively_into_anything_derived_from_it() -> None:
    graph = _calibration_graph()
    assert graph.calibration_ancestry("calibration-derived-split") == ("afm-calibration-set",)
    with pytest.raises(HeldOutValidationRefused, match="CALIBRATION_ONLY"):
        graph.assert_admissible_as_held_out_validation("calibration-derived-split")
    # ...while still being citable as biological evidence, which is a different question.
    assert "calibration-derived-split" in graph.citable_as_biological_evidence()


def test_held_out_validation_requires_an_auditable_declared_split() -> None:
    graph = _calibration_graph()
    with pytest.raises(HeldOutValidationRefused, match="held_out_split_id"):
        graph.assert_admissible_as_held_out_validation("study-without-a-split")
    graph.assert_admissible_as_held_out_validation("independent-traction-study")
    assert graph.admissible_held_out_validation_ids() == ("independent-traction-study",)


def test_a_synthetic_ancestry_also_blocks_held_out_validation() -> None:
    graph = _a6_graph()
    with pytest.raises(HeldOutValidationRefused, match="derived-from-excluded"):
        graph.assert_admissible_as_held_out_validation("morphology-statistic")


# --- source record validation -------------------------------------------------------------------------


def test_source_records_refuse_incoherent_provenance() -> None:
    with pytest.raises(ValueError, match="must not name parents"):
        _literature("primary", parents=("other",))
    with pytest.raises(ValueError, match="at least one parent"):
        SourceRecord(
            source_id="derived-from-nothing",
            provenance=ProvenanceKind.DERIVED,
            origin=SourceOrigin.LITERATURE,
            status=SourceStatus.ACTIVE,
            description="no ancestry",
            citation="doi:10.0000/x",
        )
    with pytest.raises(ValueError, match="citation"):
        SourceRecord(
            source_id="uncited",
            provenance=ProvenanceKind.SOURCED,
            origin=SourceOrigin.LITERATURE,
            status=SourceStatus.ACTIVE,
            description="claims literature with no citation",
        )
    with pytest.raises(ValueError, match="PI_GAP"):
        SourceRecord(
            source_id="gap-with-parents",
            provenance=ProvenanceKind.PI_GAP,
            origin=SourceOrigin.EXPERT_JUDGEMENT,
            status=SourceStatus.ACTIVE,
            description="declared absence of evidence",
            parents=("something",),
        )
    with pytest.raises(ValueError, match="lowercase"):
        _literature("Has Spaces And Caps")


def test_graph_hash_is_content_addressed_and_registration_order_independent() -> None:
    records = [_literature("one"), _literature("two"), _derived("three", ("one", "two"))]
    forward = SourceRegistry().register_all(records).freeze()
    backward = SourceRegistry().register_all(list(reversed(records))).freeze()
    assert forward.graph_hash == backward.graph_hash
    assert len(forward.graph_hash) == 64
    changed = SourceRegistry().register_all([*records[:2], _derived("three", ("one",))]).freeze()
    assert changed.graph_hash != forward.graph_hash


# --- public microscopy dataset adapter contract --------------------------------------------------------


def test_no_real_public_dataset_claims_to_be_integrated() -> None:
    assert BLOCKED_PUBLIC_MICROSCOPY_DATASETS
    for contract in BLOCKED_PUBLIC_MICROSCOPY_DATASETS.values():
        assert contract.integration_status is DatasetIntegrationStatus.BLOCKED_NOT_INTEGRATED
        assert not contract.satisfies_contract()
        assert contract.missing_requirements()
        assert contract.to_source_record().status is SourceStatus.BLOCKED


def test_the_worked_fixture_satisfies_every_declared_contract_requirement() -> None:
    contract, payload = worked_fixture_dataset()
    assert contract.missing_requirements() == ()
    assert contract.satisfies_contract()
    assert contract.integration_status is DatasetIntegrationStatus.IN_MEMORY_FIXTURE
    assert set(contract.channel_marker_map) == set(contract.channels) == {"dna", "actin"}
    assert contract.pixel_size_um == pytest.approx(0.325)
    assert contract.bit_depth == 16
    frames = payload["frames"]
    assert set(frames) == {"dna", "actin"}  # type: ignore[arg-type]
    assert payload["shape"] == (4, 4)
    with pytest.raises(TypeError):
        contract.channel_marker_map["dna"] = "tampered"  # type: ignore[index]


def test_a_fixture_is_never_biological_evidence_and_neither_is_anything_derived_from_it() -> None:
    contract, _ = worked_fixture_dataset()
    graph = (
        SourceRegistry()
        .register_all(
            [contract.to_source_record(), _derived("fixture-fit", (contract.dataset_id,))]
        )
        .freeze()
    )
    assert graph.exclusion(contract.dataset_id).reason is ExclusionReason.FIXTURE  # type: ignore[union-attr]
    downstream = graph.exclusion("fixture-fit")
    assert downstream is not None
    assert downstream.blocking_source_id == contract.dataset_id
    assert graph.citable_as_biological_evidence() == ()


def test_registering_the_declared_datasets_yields_an_empty_biological_corpus() -> None:
    graph = register_public_dataset_contracts(SourceRegistry()).freeze()
    assert len(graph) == len(BLOCKED_PUBLIC_MICROSCOPY_DATASETS) + 1
    assert graph.citable_as_biological_evidence() == ()
    assert graph.admissible_held_out_validation_ids() == ()
    with pytest.raises(TypeError, match="SourceRegistry"):
        register_public_dataset_contracts(object())  # type: ignore[arg-type]


def test_a_dataset_cannot_claim_integration_while_missing_a_requirement() -> None:
    with pytest.raises(ValueError, match="claims INTEGRATED but is missing"):
        MicroscopyDatasetContract(
            dataset_id="overclaiming",
            integration_status=DatasetIntegrationStatus.INTEGRATED,
            modality="widefield-fluorescence",
        )


def test_dataset_requirements_are_reported_in_declared_order() -> None:
    contract = MicroscopyDatasetContract(
        dataset_id="partly-declared",
        integration_status=DatasetIntegrationStatus.BLOCKED_NOT_INTEGRATED,
        modality="widefield-fluorescence",
        citation="doi:10.0000/example",
    )
    missing = contract.missing_requirements()
    assert "citation" not in missing
    assert "modality" not in missing
    assert set(missing) < set(REQUIRED_DATASET_REQUIREMENTS)
    assert list(missing) == [name for name in REQUIRED_DATASET_REQUIREMENTS if name in set(missing)]


def test_a_channel_marker_map_must_cover_exactly_the_declared_channels() -> None:
    with pytest.raises(ValueError, match="cover exactly"):
        MicroscopyDatasetContract(
            dataset_id="mismatched-channels",
            integration_status=DatasetIntegrationStatus.BLOCKED_NOT_INTEGRATED,
            channels=("dna", "actin"),
            channel_marker_map={"dna": "Hoechst"},
        )


def test_declared_evidence_gaps_route_to_an_expansion() -> None:
    assert DECLARED_EVIDENCE_GAPS
    for gap in DECLARED_EVIDENCE_GAPS:
        assert isinstance(gap, EvidenceGapVerdict)
        assert gap.status in {"blocked", "oracle-only", "pi-decision"}
        assert "->" in gap.routes_to or "PI" in gap.routes_to
    assert {gap.gap_id for gap in DECLARED_EVIDENCE_GAPS} == {
        "public-microscopy-not-integrated",
        "single-renderer-family",
        "calibration-holdout-partition",
    }


# --- CPU-only ------------------------------------------------------------------------------------------


def test_importing_the_optics_and_registry_modules_does_not_import_warp() -> None:
    probe = (
        "import sys;"
        "import aleph.virtual_cell.optics as o;"
        "import aleph.virtual_cell.source_registry as r;"
        "assert 'warp' not in sys.modules, sorted(k for k in sys.modules if 'warp' in k);"
        "print(o.__name__, r.__name__)"
    )
    completed = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "aleph.virtual_cell.optics" in completed.stdout
    assert "aleph.virtual_cell.source_registry" in completed.stdout
