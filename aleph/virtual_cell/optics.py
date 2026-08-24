"""Pixel-integrated optics: extended emitters and a widefield through-focus projection.

WHY this module exists.  :func:`aleph.virtual_cell.synthetic_microscopy.render_point_channel` is a
correct oracle for one regime and wrong outside it.  It evaluates the Gaussian at the pixel CENTRE,
which is a midpoint rule with a single node per pixel.  When the PSF standard deviation is comparable
to or smaller than the pixel pitch — the regime real cameras are deliberately built for — the sampled
kernel's mass depends on where inside the pixel the emitter sits, so a sub-pixel translation changes
the recorded photon count.  For a separable Gaussian the pixel integral is available in closed form
as a product of ``erf`` differences over the pixel edges, so there is no reason to accept the aliased
midpoint value.  ``PSFSampling.PIXEL_INTEGRATED`` uses that closed form; ``POINT_SAMPLED`` is retained
only so the aliasing can be MEASURED rather than asserted.

WHY extended emitters.  A cortical filament is a line segment and a membrane patch is a surface, not
a bag of points.  Rendering them as points forces the caller to choose a discretisation, and a
per-primitive normalisation then silently makes the image depend on that choice: subdividing a
filament into more pieces would brighten or dim it.  The emitters here carry emission per unit LENGTH
and per unit AREA, and their convolved FIELDS are closed-form and telescoping exactly:

* line segment: the isotropic Gaussian convolved with a straight segment factorises in the segment
  frame into ``exp(-d^2/2s^2)/(s sqrt(2 pi))`` times ``0.5 (erf(t/(s sqrt 2)) - erf((t-L)/(s sqrt 2)))``.
  Splitting ``[0, L]`` telescopes the ``erf`` pair exactly.
* triangle: applying the divergence theorem to the isotropic Gaussian turns the area integral into a
  sum of signed wedge integrals, each of which is an arctangent minus an Owen ``T`` function.  The
  wedge is exactly antisymmetric under edge reversal, so interior edges of any orientation-preserving
  subdivision cancel.

WHERE THAT TELESCOPING STOPS — the rendered IMAGE does not inherit it exactly, and this file used to
say it did.  Every primitive is rendered over its OWN window, ``bounding box + psf_radius_sigma``
sigma, because rendering each one over the whole frame is not affordable at a native filament count.
A subdivided emitter's pieces therefore have smaller windows than the whole, so after the analytic
interior cancellation the two renders discard DIFFERENT tails, and the images differ by the size of
the discarded field rather than by round-off.  Measured at ``sigma = 2 px`` with the default
``psf_radius_sigma = 4`` (128 px frame, four quadrature nodes): a 17-piece line differs from the whole
segment by ``8.90e-7`` photons, ``7.51e-6`` of the image peak, and a 64-triangle subdivision by
``1.20e-7`` / ``1.50e-6``.  At ``psf_radius_sigma = 8`` the same measurements are ``3.75e-16`` /
``3.16e-15`` and ``6.94e-17`` / ``8.67e-16``, i.e. round-off — which is why a test that forced 8 sigma
could report exact invariance while the shipped default was six orders of magnitude worse.

So the claim this module makes is scoped, not universal: **subdivision invariance holds to the
truncated PSF tail at the configured radius**, whose scale is
:func:`psf_truncation_tail_fraction` = ``exp(-R^2/2)``, and reaches round-off only once that tail
underflows double precision (``R >= 7``).  The measured peak-relative discrepancy sits below that
scale at every radius tested (``R = 4, 5, 6, 7, 8``: ``7.55e-6``, ``4.70e-8``, ``1.12e-10``,
``1.01e-13``, ``3.16e-15``).  The FRAME TOTAL is scoped the same way and by the same mechanism: at
the default radius a 17-piece line differs from the whole segment by ``3.61e-6`` of the photon
budget and a 64-triangle patch by ``4.87e-7``, because the pieces' windows do not tile the whole
emitter's window.  A caller who needs bit-level subdivision invariance must raise
``psf_radius_sigma`` and pay for it in window area.  What is NOT scoped is the emission law: the
budget still scales with length and area and never with the primitive count, which is the property
the per-unit-length and per-unit-area normalisation exists to guarantee.

WHY only widefield.  ``AxialModel.WIDEFIELD`` models defocus as a Gaussian-beam lateral broadening
``s(z) = s0 sqrt(1 + ((z - z_focus)/z_R)^2)``.  It spreads light; it never destroys it, so the frame
total is conserved while the peak drops.  A confocal pinhole and a light-sheet excitation profile are
axial GATES: they attenuate out-of-plane emitters and are NOT obtainable by broadening a widefield
kernel.  They are declared in :data:`UNIMPLEMENTED_AXIAL_MODELS` as blocked, not approximated.  A
scalar-diffraction PSF (Gibson-Lanni, Richards-Wolf) is likewise not implemented; the Gaussian-beam
width is a stated approximation to it, not a claim to be it.

Units.  Emitter rates are photons per second, per micrometre per second, and per square micrometre per
second respectively.  Rendered images are EXPECTED PHOTONS per pixel at unit detection gain, which is
the quantity :func:`~aleph.virtual_cell.synthetic_microscopy.apply_camera_model` multiplies by
``photons_per_intensity`` before Poisson sampling.  Geometry is micrometres; internally everything is
converted to pixel units, where the pixel is the unit square, so a rendered patch entry is a photon
count and never a density.

This module renders SYNTHETIC images.  Nothing it produces is biological evidence, and an encoder
trained and scored inside one renderer family is self-scoring (master plan section 10.2, A6).
:mod:`aleph.virtual_cell.source_registry` is what structurally enforces that.

Sanity Gate (self-tested in tests/ac/virtual_cell/test_optics_sources.py):
  * dimensional: emitter rate x exposure x (length | area) gives photons; the rendered frame total is
    photons and is independent of ``pixel_size_um`` for an interior emitter.
  * conservation: the pixel-integrated kernel of a point emitter sums to the emitted photon count up
    to a radius-truncation loss bounded above by ``2 erfc(R/sqrt 2)`` for ``psf_radius_sigma = R``
    (the window is a square at least ``R`` sigma wide, so the realised loss is smaller than the
    bound), and to the count exactly under ``TruncationPolicy.RENORMALIZED``.
  * boundary: PSF mass is normalised BEFORE the frame crop, so photons leaving the field of view are
    lost rather than renormalised into a bright edge (the incumbent behaviour, preserved).
  * sub-pixel invariance: translating an interior emitter by a fraction of a pixel must not change the
    frame total for the pixel-integrated kernel; it does for the point-sampled kernel once the PSF is
    narrower than a pixel.  This is the measurement protocol that separates the two samplings.
  * subdivision invariance: splitting a segment or an orientation-preserving triangle must reproduce
    the same image to within the truncated PSF tail :func:`psf_truncation_tail_fraction` at the
    configured ``psf_radius_sigma``, and to round-off once that tail underflows double precision.
    The measurement is made at the DEFAULT ``psf_radius_sigma``, because the failure this gate exists
    to catch — a per-primitive rather than per-length/per-area normalisation — is invisible at any
    radius, while the truncation limit is visible only at the shipped one.
  * numerical: the closed forms were checked against adaptive quadrature (``scipy.integrate.quad`` and
    ``dblquad``) to ~1e-16 absolute before being written down; the guard
    :data:`_DEGENERATE_EDGE_DISTANCE` only removes a removable 0/0 when the observation point lies on
    an edge line, where the exact wedge value is zero.
  * sign sense: world y handedness and the pixel-EDGE origin follow ``MicroscopyConfig`` exactly; the
    world-to-image map is a reflection, which is an isometry, so isotropic PSF geometry is unaffected.

References: Owen (1956) Ann. Math. Statist. 27, 1075 for ``T(h, a)``; Zhang & Zhu (2007) Appl. Opt. 46,
1819 for the Gaussian approximation to a widefield PSF and its axial broadening.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from types import MappingProxyType

import numpy as np
from scipy.special import erf, owens_t

from aleph.virtual_cell.synthetic_microscopy import MicroscopyConfig

__all__ = [
    "UNIMPLEMENTED_AXIAL_MODELS",
    "AxialModel",
    "Emitter",
    "LineEmitter",
    "OpticsCapabilityVerdict",
    "OpticsConfig",
    "PSFSampling",
    "PointEmitter",
    "TriangleEmitter",
    "TruncationPolicy",
    "defocused_sigma_um",
    "gaussian_line_field",
    "gaussian_pixel_weights_1d",
    "gaussian_triangle_mass",
    "psf_truncation_tail_fraction",
    "render_emitters",
    "render_widefield_stack",
]

_SQRT2 = math.sqrt(2.0)
_SQRT_2PI = math.sqrt(2.0 * math.pi)

# Relative distance below which an observation point is treated as lying ON a triangle edge line.
# The exact wedge value there is identically zero (the arctangent and Owen terms cancel), so this
# removes a removable 0/0 and introduces an error of order the square of this constant.  It is a
# numerical degeneracy guard, not a tuned threshold.
_DEGENERATE_EDGE_DISTANCE = 1e-12


class PSFSampling(StrEnum):
    """How the PSF is converted from a continuous field to per-pixel photon counts."""

    PIXEL_INTEGRATED = "pixel-integrated"
    POINT_SAMPLED = "point-sampled"


class TruncationPolicy(StrEnum):
    """What happens to the PSF mass outside the evaluated support radius."""

    ANALYTIC = "analytic"
    RENORMALIZED = "renormalized"


class AxialModel(StrEnum):
    """Axial/depth model relating an emitter's ``z`` to its contribution."""

    WIDEFIELD = "widefield"
    CONFOCAL = "confocal"
    LIGHT_SHEET = "light-sheet"


@dataclass(frozen=True, slots=True)
class OpticsCapabilityVerdict:
    """A declared-but-unimplemented optical capability and the expansion it routes to."""

    capability: str
    status: str
    why_not_approximated: str
    routes_to: str

    def __post_init__(self) -> None:
        for name in ("capability", "status", "why_not_approximated", "routes_to"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty string")
            object.__setattr__(self, name, value.strip())
        if self.status not in {"blocked", "oracle-only", "pi-decision"}:
            raise ValueError("status must be 'blocked', 'oracle-only', or 'pi-decision'")


UNIMPLEMENTED_AXIAL_MODELS: Mapping[AxialModel, OpticsCapabilityVerdict] = MappingProxyType(
    {
        AxialModel.CONFOCAL: OpticsCapabilityVerdict(
            capability="confocal pinhole axial gate",
            status="blocked",
            why_not_approximated=(
                "a pinhole ATTENUATES out-of-focus emitters; broadening a widefield kernel conserves "
                "their photons, so the widefield model cannot stand in for it"
            ),
            routes_to="observation-model extension: detection-path pinhole transmission operator",
        ),
        AxialModel.LIGHT_SHEET: OpticsCapabilityVerdict(
            capability="light-sheet selective-plane excitation",
            status="blocked",
            why_not_approximated=(
                "the axial selectivity is in the EXCITATION profile, so it multiplies the emitter "
                "rate before the PSF and cannot be folded into a detection-side kernel width"
            ),
            routes_to="observation-model extension: z-dependent excitation field times emitter rate",
        ),
    }
)


def _positive_finite(value: object, *, what: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0.0
    ):
        raise ValueError(f"{what} must be positive-finite")
    return float(value)


def _finite(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{what} must be finite")
    return float(value)


def _nonnegative_finite(value: object, *, what: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0.0
    ):
        raise ValueError(f"{what} must be finite and nonnegative")
    return float(value)


def _xy_pair(value: object, *, what: str) -> tuple[float, float]:
    pair = tuple(value)  # type: ignore[call-overload]
    if len(pair) != 2:
        raise ValueError(f"{what} must be an (x, y) pair")
    return (_finite(pair[0], what=f"{what}.x"), _finite(pair[1], what=f"{what}.y"))


@dataclass(frozen=True, slots=True)
class OpticsConfig:
    """Sampling, truncation, and axial model wrapped around an existing ``MicroscopyConfig``.

    ``microscopy.psf_sigma_um`` is the IN-FOCUS lateral standard deviation; ``rayleigh_range_um`` is
    the defocus at which the widefield lateral width has grown by ``sqrt(2)``.

    ``pixel_quadrature_nodes`` applies to LINE and TRIANGLE emitters only.  Their convolved fields are
    closed form, but an oblique segment or edge is not separable in the pixel frame, so the pixel
    integral is a tensor Gauss-Legendre rule of this order per axis rather than a closed form.  The
    point emitter never uses it.  The rule affects per-pixel values, not the frame total: summing a
    Gauss-Legendre rule over the whole pixel grid is a composite rule over the plane whose error is
    exponentially small in PSF width over pixel pitch, so the emitted photon budget is recovered to
    round-off even at one node.  Per-pixel convergence is spectral (measured against a 16-node
    reference at ``sigma = 2 px``: 1.2e-2 at one node, 3.2e-5 at two, 1.1e-10 at the default four).
    """

    microscopy: MicroscopyConfig
    sampling: PSFSampling = PSFSampling.PIXEL_INTEGRATED
    truncation: TruncationPolicy = TruncationPolicy.ANALYTIC
    axial_model: AxialModel = AxialModel.WIDEFIELD
    rayleigh_range_um: float = 0.5
    focal_plane_z_um: float = 0.0
    exposure_s: float = 1.0
    pixel_quadrature_nodes: int = 4

    def __post_init__(self) -> None:
        if not isinstance(self.microscopy, MicroscopyConfig):
            raise TypeError("microscopy must be a MicroscopyConfig")
        for name, enum_type in (
            ("sampling", PSFSampling),
            ("truncation", TruncationPolicy),
            ("axial_model", AxialModel),
        ):
            if not isinstance(getattr(self, name), enum_type):
                raise TypeError(f"{name} must be a {enum_type.__name__}")
        if self.axial_model in UNIMPLEMENTED_AXIAL_MODELS:
            verdict = UNIMPLEMENTED_AXIAL_MODELS[self.axial_model]
            raise NotImplementedError(
                f"{verdict.capability} is {verdict.status}: {verdict.why_not_approximated}; "
                f"routes to {verdict.routes_to}"
            )
        object.__setattr__(
            self,
            "rayleigh_range_um",
            _positive_finite(self.rayleigh_range_um, what="rayleigh_range_um"),
        )
        object.__setattr__(
            self, "focal_plane_z_um", _finite(self.focal_plane_z_um, what="focal_plane_z_um")
        )
        object.__setattr__(self, "exposure_s", _positive_finite(self.exposure_s, what="exposure_s"))
        if (
            isinstance(self.pixel_quadrature_nodes, bool)
            or not isinstance(self.pixel_quadrature_nodes, int)
            or not 1 <= self.pixel_quadrature_nodes <= 32
        ):
            raise ValueError("pixel_quadrature_nodes must be an integer in [1, 32]")

    @property
    def renderer_id(self) -> str:
        """Stable identifier of the rendering family, for observation provenance."""
        return f"optics-{self.sampling.value}-{self.truncation.value}-{self.axial_model.value}@1"

    @property
    def config_hash(self) -> str:
        """SHA-256 over the optics configuration, including the wrapped microscope."""
        payload = json.dumps(
            {
                "microscopy_config_sha256": self.microscopy.config_hash,
                "sampling": self.sampling.value,
                "truncation": self.truncation.value,
                "axial_model": self.axial_model.value,
                "rayleigh_range_um": self.rayleigh_range_um,
                "focal_plane_z_um": self.focal_plane_z_um,
                "exposure_s": self.exposure_s,
                "pixel_quadrature_nodes": self.pixel_quadrature_nodes,
                "renderer_id": self.renderer_id,
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
        return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True, slots=True)
class PointEmitter:
    """A point label emitting ``photon_rate`` photons per second."""

    xy_um: tuple[float, float]
    photon_rate: float
    z_um: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "xy_um", _xy_pair(self.xy_um, what="xy_um"))
        object.__setattr__(
            self, "photon_rate", _nonnegative_finite(self.photon_rate, what="photon_rate")
        )
        object.__setattr__(self, "z_um", _finite(self.z_um, what="z_um"))

    def total_photons(self, exposure_s: float) -> float:
        """Photons emitted over ``exposure_s`` seconds."""
        return self.photon_rate * exposure_s


@dataclass(frozen=True, slots=True)
class LineEmitter:
    """A straight filament segment emitting per unit LENGTH.

    ``photon_rate_per_length`` is photons per micrometre per second, so the emitted total scales with
    segment length and is invariant to how the segment is subdivided.
    """

    start_xy_um: tuple[float, float]
    end_xy_um: tuple[float, float]
    photon_rate_per_length: float
    z_um: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "start_xy_um", _xy_pair(self.start_xy_um, what="start_xy_um"))
        object.__setattr__(self, "end_xy_um", _xy_pair(self.end_xy_um, what="end_xy_um"))
        object.__setattr__(
            self,
            "photon_rate_per_length",
            _nonnegative_finite(self.photon_rate_per_length, what="photon_rate_per_length"),
        )
        object.__setattr__(self, "z_um", _finite(self.z_um, what="z_um"))
        if self.length_um <= 0.0:
            raise ValueError("line emitter endpoints must be distinct")

    @property
    def length_um(self) -> float:
        """Segment length in micrometres."""
        return math.hypot(
            self.end_xy_um[0] - self.start_xy_um[0],
            self.end_xy_um[1] - self.start_xy_um[1],
        )

    def total_photons(self, exposure_s: float) -> float:
        """Photons emitted over ``exposure_s`` seconds: rate x length x exposure."""
        return self.photon_rate_per_length * self.length_um * exposure_s


@dataclass(frozen=True, slots=True)
class TriangleEmitter:
    """A membrane/surface patch emitting per unit AREA.

    ``photon_rate_per_area`` is photons per square micrometre per second, so the emitted total scales
    with patch area and is invariant to orientation-preserving subdivision of the triangle.
    """

    vertices_xy_um: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    photon_rate_per_area: float
    z_um: float = 0.0

    def __post_init__(self) -> None:
        vertices = tuple(self.vertices_xy_um)
        if len(vertices) != 3:
            raise ValueError("vertices_xy_um must contain exactly three (x, y) pairs")
        object.__setattr__(
            self,
            "vertices_xy_um",
            tuple(_xy_pair(vertex, what="vertices_xy_um entry") for vertex in vertices),
        )
        object.__setattr__(
            self,
            "photon_rate_per_area",
            _nonnegative_finite(self.photon_rate_per_area, what="photon_rate_per_area"),
        )
        object.__setattr__(self, "z_um", _finite(self.z_um, what="z_um"))
        if self.area_um2 <= 0.0:
            raise ValueError("triangle emitter must have positive area")

    @property
    def area_um2(self) -> float:
        """Unsigned triangle area in square micrometres."""
        (ax, ay), (bx, by), (cx, cy) = self.vertices_xy_um
        return 0.5 * abs((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))

    def total_photons(self, exposure_s: float) -> float:
        """Photons emitted over ``exposure_s`` seconds: rate x area x exposure."""
        return self.photon_rate_per_area * self.area_um2 * exposure_s


type Emitter = PointEmitter | LineEmitter | TriangleEmitter


def defocused_sigma_um(z_um: float, config: OpticsConfig) -> float:
    """Widefield lateral PSF width at depth ``z_um``.

    Args:
        z_um: Emitter axial position in micrometres.
        config: Optics configuration supplying the in-focus width, focal plane, and Rayleigh range.

    Returns:
        Lateral Gaussian standard deviation in micrometres, never below the in-focus width.
    """
    if not isinstance(config, OpticsConfig):
        raise TypeError("config must be an OpticsConfig")
    defocus = _finite(z_um, what="z_um") - config.focal_plane_z_um
    broadening = math.sqrt(1.0 + (defocus / config.rayleigh_range_um) ** 2)
    return config.microscopy.psf_sigma_um * broadening


def psf_truncation_tail_fraction(psf_radius_sigma: float) -> float:
    """Scale of the PSF field a per-primitive render window discards, relative to the peak.

    This is the value of a unit-peak 1-D Gaussian at ``psf_radius_sigma`` sigma, ``exp(-R^2/2)``.  It
    is not a tuned tolerance: it is the size of the field that lies just outside the window
    :func:`render_emitters` evaluates, and therefore the scale at which subdivision invariance of the
    rendered IMAGE stops holding, even though the underlying analytic field telescopes exactly.

    The counterexample that made this function necessary: at the default ``psf_radius_sigma = 4`` a
    17-piece subdivision of one line emitter differs from the whole segment by ``7.51e-6`` of the
    image peak, and a 64-triangle subdivision by ``1.50e-6`` — because each piece's window clips a
    different tail.  Raising the radius to 8 drops both to round-off, which is exactly how the
    original test hid the effect.  Every measured peak-relative discrepancy sits below the value
    returned here (``R = 4, 5, 6, 7, 8``: measured ``7.55e-6``, ``4.70e-8``, ``1.12e-10``, ``1.01e-13``,
    ``3.16e-15`` against ``3.35e-4``, ``3.73e-6``, ``1.52e-8``, ``2.29e-11``, ``1.27e-14``).

    Args:
        psf_radius_sigma: Truncation radius in PSF standard deviations, as carried by
            ``MicroscopyConfig.psf_radius_sigma``.

    Returns:
        ``exp(-psf_radius_sigma^2 / 2)``, dimensionless and in ``(0, 1]``.

    Raises:
        ValueError: If ``psf_radius_sigma`` is not positive-finite.
    """
    radius = _positive_finite(psf_radius_sigma, what="psf_radius_sigma")
    return math.exp(-0.5 * radius * radius)


def gaussian_pixel_weights_1d(
    first_index: int,
    count: int,
    centre_px: float,
    sigma_px: float,
) -> np.ndarray:
    """Exact integral of a unit-mass 1-D Gaussian over ``count`` consecutive unit pixels.

    Pixel ``i`` spans ``[i - 0.5, i + 0.5]`` in the centre-index coordinate used by
    ``MicroscopyConfig``.  The integral is the closed-form ``erf`` difference over the pixel EDGES;
    no quadrature is used.

    Args:
        first_index: Index of the first pixel in the window.
        count: Number of consecutive pixels.
        centre_px: Gaussian centre in centre-index pixel coordinates.
        sigma_px: Gaussian standard deviation in pixels.

    Returns:
        Array of shape ``(count,)`` whose entries sum to the Gaussian mass inside the window.
    """
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("count must be a positive integer")
    sigma = _positive_finite(sigma_px, what="sigma_px")
    centre = _finite(centre_px, what="centre_px")
    edges = np.arange(first_index, first_index + count + 1, dtype=np.float64) - 0.5
    cumulative = 0.5 * erf((edges - centre) / (sigma * _SQRT2))
    return np.diff(cumulative)


def gaussian_line_field(
    start_px: np.ndarray,
    end_px: np.ndarray,
    sigma_px: float,
    points_px: np.ndarray,
) -> np.ndarray:
    """Closed-form isotropic-Gaussian convolution of a unit-density straight segment.

    The value is photons per unit pixel area for a segment carrying one photon per unit pixel length.
    In the segment frame the field factorises into a lateral Gaussian and an axial ``erf`` pair, and
    the ``erf`` pair telescopes exactly when the segment is split, which is what makes the rendered
    image subdivision invariant.

    Args:
        start_px: Segment start in centre-index pixel coordinates, shape ``(2,)``.
        end_px: Segment end in centre-index pixel coordinates, shape ``(2,)``.
        sigma_px: Lateral PSF standard deviation in pixels.
        points_px: Evaluation points, shape ``(m, 2)``.

    Returns:
        Array of shape ``(m,)``.
    """
    start = np.asarray(start_px, dtype=np.float64).reshape(2)
    end = np.asarray(end_px, dtype=np.float64).reshape(2)
    points = np.asarray(points_px, dtype=np.float64).reshape(-1, 2)
    sigma = _positive_finite(sigma_px, what="sigma_px")
    edge = end - start
    length = float(np.hypot(edge[0], edge[1]))
    if length <= 0.0:
        raise ValueError("line field needs distinct endpoints")
    axis = edge / length
    normal = np.array([axis[1], -axis[0]], dtype=np.float64)
    relative = points - start
    along = relative @ axis
    across = relative @ normal
    axial = 0.5 * (erf(along / (sigma * _SQRT2)) - erf((along - length) / (sigma * _SQRT2)))
    lateral = np.exp(-0.5 * np.square(across / sigma)) / (sigma * _SQRT_2PI)
    return lateral * axial


def _signed_wedge(
    from_px: np.ndarray,
    to_px: np.ndarray,
    sigma_px: float,
) -> np.ndarray:
    """Signed Gaussian mass of the triangle (observation point, ``from_px``, ``to_px``).

    Both vectors are given RELATIVE to the observation point.  The value is exactly antisymmetric
    under swapping the two vectors, which is what cancels interior edges under subdivision.
    """
    edge = to_px - from_px
    length = np.hypot(edge[:, 0], edge[:, 1])
    result = np.zeros(from_px.shape[0], dtype=np.float64)
    live = length > 0.0
    if not np.any(live):
        return result
    axis_x = edge[live, 0] / length[live]
    axis_y = edge[live, 1] / length[live]
    normal_x, normal_y = axis_y, -axis_x
    distance = normal_x * from_px[live, 0] + normal_y * from_px[live, 1]
    along_from = axis_x * from_px[live, 0] + axis_y * from_px[live, 1]
    along_to = axis_x * to_px[live, 0] + axis_y * to_px[live, 1]
    regular = np.abs(distance) > _DEGENERATE_EDGE_DISTANCE * sigma_px
    ratio_from = np.zeros_like(distance)
    ratio_to = np.zeros_like(distance)
    ratio_from[regular] = along_from[regular] / distance[regular]
    ratio_to[regular] = along_to[regular] / distance[regular]
    height = np.abs(distance) / sigma_px
    wedge = (np.arctan(ratio_to) - np.arctan(ratio_from)) / (2.0 * math.pi) - (
        owens_t(height, ratio_to) - owens_t(height, ratio_from)
    )
    wedge[~regular] = 0.0
    result[live] = wedge
    return result


def gaussian_triangle_mass(
    vertices_px: np.ndarray,
    sigma_px: float,
    points_px: np.ndarray,
) -> np.ndarray:
    """Closed-form isotropic-Gaussian mass inside a triangle, seen from each evaluation point.

    Applying the divergence theorem to the radial Gaussian turns the area integral into three signed
    wedge integrals, each an arctangent minus an Owen ``T`` function.  The result is in ``[0, 1]``:
    it is the fraction of a unit-mass PSF centred at the evaluation point that falls inside the
    triangle, so multiplying by the patch surface density gives photons per unit pixel area.

    Args:
        vertices_px: Triangle vertices in centre-index pixel coordinates, shape ``(3, 2)``.
        sigma_px: PSF standard deviation in pixels.
        points_px: Evaluation points, shape ``(m, 2)``.

    Returns:
        Array of shape ``(m,)``.
    """
    vertices = np.asarray(vertices_px, dtype=np.float64).reshape(3, 2)
    points = np.asarray(points_px, dtype=np.float64).reshape(-1, 2)
    sigma = _positive_finite(sigma_px, what="sigma_px")
    first = vertices[1] - vertices[0]
    second = vertices[2] - vertices[0]
    signed_area2 = float(first[0] * second[1] - first[1] * second[0])
    if signed_area2 == 0.0:
        raise ValueError("triangle mass needs a non-degenerate triangle")
    total = np.zeros(points.shape[0], dtype=np.float64)
    for index in range(3):
        total += _signed_wedge(
            vertices[index][None, :] - points,
            vertices[(index + 1) % 3][None, :] - points,
            sigma,
        )
    return total if signed_area2 > 0.0 else -total


def _pixel_quadrature(nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre nodes/weights on the unit pixel ``[-0.5, 0.5]``; weights sum to one."""
    positions, weights = np.polynomial.legendre.leggauss(nodes)
    return positions * 0.5, weights * 0.5


def _to_pixel_frame(
    points_xy_um: np.ndarray,
    config: MicroscopyConfig,
    origin_xy_um: tuple[float, float],
) -> np.ndarray:
    """Map world micrometres to centre-index pixel coordinates, honouring y handedness."""
    points = np.asarray(points_xy_um, dtype=np.float64).reshape(-1, 2)
    origin = np.asarray(origin_xy_um, dtype=np.float64).reshape(2)
    physical = (points - origin) / config.pixel_size_um
    column = physical[:, 0] - 0.5
    if config.world_y_direction == "up":
        row = config.height_px - physical[:, 1] - 0.5
    else:
        row = physical[:, 1] - 0.5
    return np.stack([column, row], axis=1)


def _window(
    min_px: float,
    max_px: float,
    radius_px: int,
    extent_px: int,
) -> tuple[int, int, int, int] | None:
    """Return ``(full_start, full_stop, clipped_start, clipped_stop)`` or ``None`` if fully outside."""
    full_start = int(math.floor(min_px)) - radius_px
    full_stop = int(math.floor(max_px)) + radius_px + 2
    start = max(0, full_start)
    stop = min(extent_px, full_stop)
    if start >= stop:
        return None
    return full_start, full_stop, start, stop


def _point_patch(
    centre_px: np.ndarray,
    sigma_px: float,
    sampling: PSFSampling,
    columns: tuple[int, int],
    rows: tuple[int, int],
) -> np.ndarray:
    """Unit-mass PSF patch over a rectangular window of pixels."""
    first_column, column_count = columns
    first_row, row_count = rows
    if sampling is PSFSampling.PIXEL_INTEGRATED:
        weights_x = gaussian_pixel_weights_1d(first_column, column_count, centre_px[0], sigma_px)
        weights_y = gaussian_pixel_weights_1d(first_row, row_count, centre_px[1], sigma_px)
        return weights_y[:, None] * weights_x[None, :]
    xs = np.arange(first_column, first_column + column_count, dtype=np.float64)
    ys = np.arange(first_row, first_row + row_count, dtype=np.float64)
    exponent = np.square(xs - centre_px[0])[None, :] + np.square(ys - centre_px[1])[:, None]
    return np.exp(-0.5 * exponent / (sigma_px * sigma_px)) / (2.0 * math.pi * sigma_px * sigma_px)


def _quadrature_nodes(
    first_column: int,
    column_count: int,
    first_row: int,
    row_count: int,
    nodes: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Flattened evaluation points inside a pixel window, plus the 1-D quadrature weights."""
    offsets, weights = _pixel_quadrature(nodes)
    columns = np.arange(first_column, first_column + column_count, dtype=np.float64)
    rows = np.arange(first_row, first_row + row_count, dtype=np.float64)
    node_x = (columns[:, None] + offsets[None, :]).ravel()
    node_y = (rows[:, None] + offsets[None, :]).ravel()
    grid_x, grid_y = np.meshgrid(node_x, node_y)
    points = np.stack([grid_x.ravel(), grid_y.ravel()], axis=1)
    return points, weights


def _contract_quadrature(
    values: np.ndarray,
    row_count: int,
    column_count: int,
    nodes: int,
    weights: np.ndarray,
) -> np.ndarray:
    """Contract per-node field values into per-pixel integrals."""
    shaped = values.reshape(row_count, nodes, column_count, nodes)
    return np.einsum("aibj,i,j->ab", shaped, weights, weights)


def _render_point(
    emitter: PointEmitter,
    optics: OpticsConfig,
    origin_xy_um: tuple[float, float],
    image: np.ndarray,
) -> None:
    config = optics.microscopy
    total = emitter.total_photons(optics.exposure_s)
    if total == 0.0:
        return
    sigma_px = defocused_sigma_um(emitter.z_um, optics) / config.pixel_size_um
    radius_px = max(1, int(math.ceil(config.psf_radius_sigma * sigma_px)))
    centre = _to_pixel_frame(np.array([emitter.xy_um]), config, origin_xy_um)[0]
    horizontal = _window(centre[0], centre[0], radius_px, config.width_px)
    vertical = _window(centre[1], centre[1], radius_px, config.height_px)
    if horizontal is None or vertical is None:
        return
    full_x0, full_x1, x0, x1 = horizontal
    full_y0, full_y1, y0, y1 = vertical
    patch = _point_patch(
        centre,
        sigma_px,
        optics.sampling,
        (full_x0, full_x1 - full_x0),
        (full_y0, full_y1 - full_y0),
    )
    patch = patch * total
    if optics.truncation is TruncationPolicy.RENORMALIZED:
        patch_sum = float(np.sum(patch))
        if patch_sum == 0.0:
            return
        patch = patch * (total / patch_sum)
    image[y0:y1, x0:x1] += patch[y0 - full_y0 : y1 - full_y0, x0 - full_x0 : x1 - full_x0]


def _render_extended(
    emitter: LineEmitter | TriangleEmitter,
    optics: OpticsConfig,
    origin_xy_um: tuple[float, float],
    image: np.ndarray,
) -> None:
    config = optics.microscopy
    total = emitter.total_photons(optics.exposure_s)
    if total == 0.0:
        return
    sigma_px = defocused_sigma_um(emitter.z_um, optics) / config.pixel_size_um
    radius_px = max(1, int(math.ceil(config.psf_radius_sigma * sigma_px)))
    if isinstance(emitter, LineEmitter):
        geometry_um = np.array([emitter.start_xy_um, emitter.end_xy_um], dtype=np.float64)
    else:
        geometry_um = np.array(emitter.vertices_xy_um, dtype=np.float64)
    geometry = _to_pixel_frame(geometry_um, config, origin_xy_um)
    # This window is PER PRIMITIVE, and that is what breaks exact subdivision invariance of the
    # rendered image: the pieces of a split emitter have smaller windows than the whole, so the two
    # renders discard different tails after the analytic interior cancellation.  Counterexample at the
    # default psf_radius_sigma = 4, sigma = 2 px: a 17-piece line differs by 8.90e-7 photons
    # (7.51e-6 of the peak) and a 64-triangle subdivision by 1.20e-7 (1.50e-6); at radius 8 both are
    # round-off.  The frame total moves with it: 3.61e-6 and 4.87e-7 relative for the same two
    # subdivisions.  A window shared across the primitives of one emitter would restore exactness, but
    # render_emitters receives a FLAT sequence with no grouping, and widening every window to the
    # union bounding box makes the cost of a native filament population quadratic in the frame.  The
    # claim is therefore scoped to psf_truncation_tail_fraction(psf_radius_sigma) and measured at the
    # default, rather than being restated at a radius where it happens to hold.
    horizontal = _window(
        float(np.min(geometry[:, 0])), float(np.max(geometry[:, 0])), radius_px, config.width_px
    )
    vertical = _window(
        float(np.min(geometry[:, 1])), float(np.max(geometry[:, 1])), radius_px, config.height_px
    )
    if horizontal is None or vertical is None:
        return
    full_x0, full_x1, x0, x1 = horizontal
    full_y0, full_y1, y0, y1 = vertical
    column_count = full_x1 - full_x0
    row_count = full_y1 - full_y0
    nodes = optics.pixel_quadrature_nodes
    points, weights = _quadrature_nodes(full_x0, column_count, full_y0, row_count, nodes)
    if isinstance(emitter, LineEmitter):
        length_px = float(np.hypot(*(geometry[1] - geometry[0])))
        density = total / length_px
        field = density * gaussian_line_field(geometry[0], geometry[1], sigma_px, points)
    else:
        first = geometry[1] - geometry[0]
        second = geometry[2] - geometry[0]
        area_px = 0.5 * abs(float(first[0] * second[1] - first[1] * second[0]))
        density = total / area_px
        field = density * gaussian_triangle_mass(geometry, sigma_px, points)
    patch = _contract_quadrature(field, row_count, column_count, nodes, weights)
    if optics.truncation is TruncationPolicy.RENORMALIZED:
        patch_sum = float(np.sum(patch))
        if patch_sum == 0.0:
            return
        patch = patch * (total / patch_sum)
    image[y0:y1, x0:x1] += patch[y0 - full_y0 : y1 - full_y0, x0 - full_x0 : x1 - full_x0]


def render_emitters(
    emitters: Sequence[Emitter],
    optics: OpticsConfig,
    *,
    origin_xy_um: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    """Render point, line, and triangle emitters into one expected-photon frame.

    Each emitter contributes at the lateral PSF width implied by its own ``z_um``, so the returned
    frame is the widefield z-projection of the emitter distribution at the configured focal plane.
    The PSF is normalised BEFORE the frame crop, so photons leaving the field of view are lost rather
    than renormalised into a boundary ghost.

    Args:
        emitters: Emitters to render; an empty sequence returns the background frame.
        optics: Sampling, truncation, axial model, exposure, and wrapped microscope geometry.
        origin_xy_um: World coordinate of the lower-left image EDGE when world y points up, and of
            the upper-left edge when it points down.

    Returns:
        Read-only ``(height_px, width_px)`` array of expected photons per pixel.

    Raises:
        TypeError: If ``optics`` is not an ``OpticsConfig`` or an emitter has an unsupported type.
        ValueError: If ``origin_xy_um`` is not a finite ``(x, y)`` pair.
    """
    if not isinstance(optics, OpticsConfig):
        raise TypeError("optics must be an OpticsConfig")
    if isinstance(emitters, (PointEmitter, LineEmitter, TriangleEmitter)):
        raise TypeError("emitters must be a sequence of emitters, not a single emitter")
    origin = _xy_pair(origin_xy_um, what="origin_xy_um")
    config = optics.microscopy
    image = np.full((config.height_px, config.width_px), float(config.background), dtype=np.float64)
    for emitter in emitters:
        if isinstance(emitter, PointEmitter):
            _render_point(emitter, optics, origin, image)
        elif isinstance(emitter, (LineEmitter, TriangleEmitter)):
            _render_extended(emitter, optics, origin, image)
        else:
            raise TypeError(f"unsupported emitter type: {type(emitter).__name__}")
    image.setflags(write=False)
    return image


def render_widefield_stack(
    emitters: Sequence[Emitter],
    optics: OpticsConfig,
    focal_planes_um: Sequence[float],
    *,
    origin_xy_um: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    """Render one widefield frame per focal plane, giving a through-focus z stack.

    Every plane sees every emitter: widefield collects out-of-focus light rather than rejecting it,
    which is exactly why a widefield stack is not a confocal stack.

    Args:
        emitters: Emitters to render.
        optics: Optics configuration; its ``focal_plane_z_um`` is overridden per plane.
        focal_planes_um: Focal-plane positions in micrometres, at least one.
        origin_xy_um: As in :func:`render_emitters`.

    Returns:
        Read-only ``(len(focal_planes_um), height_px, width_px)`` array of expected photons.

    Raises:
        ValueError: If ``focal_planes_um`` is empty or contains a non-finite entry.
    """
    planes = [_finite(plane, what="focal_planes_um entry") for plane in focal_planes_um]
    if not planes:
        raise ValueError("focal_planes_um must contain at least one plane")
    frames = []
    for plane in planes:
        plane_optics = OpticsConfig(
            microscopy=optics.microscopy,
            sampling=optics.sampling,
            truncation=optics.truncation,
            axial_model=optics.axial_model,
            rayleigh_range_um=optics.rayleigh_range_um,
            focal_plane_z_um=plane,
            exposure_s=optics.exposure_s,
            pixel_quadrature_nodes=optics.pixel_quadrature_nodes,
        )
        frames.append(render_emitters(emitters, plane_optics, origin_xy_um=origin_xy_um))
    stack = np.stack(frames, axis=0)
    stack.setflags(write=False)
    return stack
