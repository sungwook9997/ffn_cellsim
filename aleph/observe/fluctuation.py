r"""The passive membrane thermal fluctuation spectrum — Aleph's first observation operator.

`ALEPH-DQ-102` chose this modality first because it needs no probe geometry, no contact model and no
inverse-model regularization, which removes the three largest sources of protocol ambiguity from the
very first measurement the project makes. What it produces is a mode power spectrum ``P(q)``; turning
that into ``(kappa, sigma)`` is inference (L8) and judging it against the closed form is the analytic
oracle's job (L4). Neither is imported here, and the firewall test enforces the second.

Units, matching L2's agreed convention: ``q`` in 1/µm, areas in µm², mode amplitudes in µm.

1. Projection
-------------
Write the normal displacement field on ``V`` vertices as ``h_v(t)`` [µm]. A mode basis supplies
``Phi`` of shape ``(V, M)`` whose columns ``phi_n`` each carry a wavevector magnitude ``q_n`` [1/µm],
plus per-vertex areas ``A_v`` [µm²]. The basis is orthonormal under the **area-weighted** inner
product — the discretisation of the continuum ``L²`` product on the surface::

    <f, g> = sum_v A_v f_v g_v / A_total,    <phi_m, phi_n> = delta_mn

so that::

    a_n(t) = sum_v A_v phi_n(v) h_v(t) / A_total          [µm]

The area weighting is not cosmetic. An unweighted dot product weights densely meshed regions more
heavily, so the recovered spectrum picks up whatever the mesher happened to refine. That is a mesh
artefact and it would be read as physics. ``test_unweighted_projection_biases_a_nonuniform_mesh``
demonstrates the failure directly, so the justification is executable rather than asserted.

2. Power spectrum
-----------------
``P(q_n) = < |a_n(t)|^2 >_t`` [µm²], with the per-mode mean subtracted first and ``ddof = 1``.

Not subtracting the mean measures the *static* shape of the configuration rather than its
fluctuation; for a mode with a non-zero equilibrium offset that offset dominates the variance
outright. ``ddof = 1`` because the mean was estimated from the same samples.

Per-mode sample counts are reported twice, and the distinction carries weight:

``n_samples``
    Frames that entered the estimate.

``n_effective``
    The independent-sample count from :mod:`aleph.observe.stationarity`. Consecutive frames of a
    relaxing membrane are strongly correlated, so ``n_samples`` overstates the precision of
    ``P(q_n)`` by ``sqrt(2 * tau_samples)`` — the ALEPH-PORT-501 correction, applied where it
    actually bites. Long-wavelength modes relax slowly, so the inflation is *worst at small ``q``*,
    which is exactly the end a bending-rigidity fit leans on.

The reported standard error uses ``n_effective``. Never ``n_samples``.

3. The trusted wavevector range
-------------------------------
A discrete mesh cannot report a spectrum at arbitrary ``q``. Two bounds, both measured from the mesh
rather than supplied:

**Upper — the mesh Nyquist limit.** With mean vertex spacing ``a`` [µm] the shortest resolvable
wavelength is ``2a``, so ``q_max = pi / a``. Above it a mode is aliased: the sampled field cannot
distinguish it from a longer-wavelength one, and any power reported there is some other mode's power
wearing this mode's label. This is the bound most easily violated by accident, because a mode basis
is usually generated from a formula that will happily produce as many modes as it is asked for.

**Lower — the finite patch.** The longest wavelength that fits in a patch of area ``A_total`` is
``L = sqrt(A_total)``, so ``q_min = 2*pi / L``. Below it no full period fits inside the domain; the
"mode" is a gradient across the patch, and its variance is set by the boundary condition rather than
by the fluctuation spectrum.

Modes outside ``[q_min, q_max]`` are **excluded from the reported spectrum and listed by index with
a reason**. They are not silently dropped, and they are not reported with a caveat attached to a
number that a reader will use anyway. ``trusted_wavevector_range`` travels in the result, because a
spectrum without a stated range is a spectrum somebody will extrapolate.

4. Coordination with ``aleph.state.spectral``
---------------------------------------------
L2 owns the mode basis on the real manifold, and at the time of writing ``aleph/state/spectral.py``
does not exist. Rather than block or guess at its API, this module defines the :class:`ModeBasis`
**Protocol** — the concept and the units, not the implementation. Any L2 object exposing those
members satisfies it structurally, so neither lane imports the other and neither can break the
other.

:class:`PlaneWaveModeBasis` ships here as the concrete implementation for a flat periodic patch.
That is not a stand-in: a flat periodic patch *is* the Helfrich geometry of `ALEPH-DQ-102`, the
geometry the closed-form oracle is written for. When L2's curved basis lands it should satisfy the
same Protocol and reduce to this one in the flat limit — a cross-lane test worth writing once both
exist.

5. Why the apparent quantity is the equipartition ratio
-------------------------------------------------------
:meth:`FluctuationSpectrumOperator.raw_observable` returns the spectrum — what the instrument
produces. :meth:`~FluctuationSpectrumOperator.apparent_quantity` reports one scalar from it, and the
default is the mean mode energy in units of ``kB*T/2``.

That is deliberately the quantity that can be *unambiguously wrong*. `ALEPH-DQ-102`'s argument for
this modality is that on the exact taut-string case equipartition is an identity, not a limit: an
estimator that cannot recover ``<E_n> = kB*T/2`` there is simply broken, with no fitting and no
threshold involved. The ``(kappa, sigma)`` fit is a harder question with a known degeneracy
structure and it lives elsewhere.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final, Protocol, runtime_checkable

import numpy as np

from aleph.observe.manifest import ProtocolManifest
from aleph.observe.operator import (
    AcceptedState,
    ApparentQuantity,
    ApplicabilityDomain,
    ObservationContext,
    RawObservable,
    Refusal,
    RefusalCode,
)
from aleph.observe.stationarity import integrated_autocorrelation_time

__all__ = [
    "FluctuationSpectrumOperator",
    "ModeBasis",
    "ModeExclusion",
    "PlaneWaveModeBasis",
    "SpectrumResult",
    "TrustedRange",
    "mesh_nyquist_wavevector",
    "patch_fundamental_wavevector",
]

_OPERATOR_ID: Final[str] = "aleph.observe.fluctuation.FluctuationSpectrumOperator/v1"

#: Minimum frames before a variance means anything. Two frames give one degree of freedom, which is
#: a number but not an estimate; four is the same floor the correlation-time estimator imposes.
MIN_FRAMES: Final[int] = 4

EPS64: Final[float] = float(np.finfo(np.float64).eps)

#: Round-off allowance on the *lower* band edge, in epsilons per vertex.
#:
#: ``q_min = 2*pi/sqrt(sum_v A_v)`` and a basis fundamental of ``2*pi/L`` are the same number
#: reached by different arithmetic, so a bare ``q < q_min`` has no margin at all and drops the
#: longest wavelength in the patch whenever the last bit lands the wrong way. The area sum
#: accumulates round-off across the vertices and the square root carries half of it through, so the
#: allowance is scaled by the vertex count at the call site. It is deliberately far too small to
#: admit anything physical: a mode genuinely longer than the patch sits a *factor* below the
#: fundamental, not an ULP, and `test_a_mode_genuinely_longer_than_the_patch_is_still_excluded`
#: holds that line.
BAND_EDGE_EPS_PER_VERTEX: Final[float] = 8.0


# ---------------------------------------------------------------------------------------------
# The mode basis seam
# ---------------------------------------------------------------------------------------------


@runtime_checkable
class ModeBasis(Protocol):
    """A set of surface modes with wavevectors, areas, and a spacing — units fixed, source free.

    This is the seam with L2's ``aleph.state.spectral``. It is a Protocol so that neither lane has
    to import the other: any object carrying these five members satisfies it structurally.

    ``basis_matrix``
        ``(V, M)`` float array. Column ``n`` is mode ``n`` sampled at the ``V`` vertices.
        Orthonormal under the area-weighted inner product of §1.

    ``wavevectors_inv_um``
        ``(M,)`` float array of wavevector magnitudes ``q_n`` [1/µm], one per column.

    ``vertex_areas_um2``
        ``(V,)`` float array of per-vertex areas [µm²] — the quadrature weights of the inner
        product, not a decoration.

    ``total_area_um2``
        ``sum_v A_v`` [µm²]. Sets the lower trusted bound.

    ``mean_vertex_spacing_um``
        Mean spacing between neighbouring vertices [µm]. Sets the Nyquist limit, so it must be the
        *sampling* spacing, not a bounding-box size.
    """

    @property
    def basis_matrix(self) -> np.ndarray: ...

    @property
    def wavevectors_inv_um(self) -> np.ndarray: ...

    @property
    def vertex_areas_um2(self) -> np.ndarray: ...

    @property
    def total_area_um2(self) -> float: ...

    @property
    def mean_vertex_spacing_um(self) -> float: ...


def mesh_nyquist_wavevector(mean_vertex_spacing_um: float) -> float:
    """``q_max = pi / a`` — the largest wavevector a mesh of spacing ``a`` can represent.

    Two samples per period is the minimum that distinguishes a wave from its aliases, so the
    shortest resolvable wavelength is ``2a`` and ``q_max = 2*pi/(2a) = pi/a``.
    """
    spacing = float(mean_vertex_spacing_um)
    if not (math.isfinite(spacing) and spacing > 0.0):
        raise ValueError(f"mean vertex spacing must be positive and finite [µm]; got {spacing!r}")
    return math.pi / spacing


def patch_fundamental_wavevector(total_area_um2: float) -> float:
    """``q_min = 2*pi / sqrt(A)`` — the longest wavelength that fits in a patch of area ``A``."""
    area = float(total_area_um2)
    if not (math.isfinite(area) and area > 0.0):
        raise ValueError(f"total area must be positive and finite [µm²]; got {area!r}")
    return 2.0 * math.pi / math.sqrt(area)


@dataclass(frozen=True, slots=True)
class PlaneWaveModeBasis:
    """Real plane-wave modes on a square periodic patch — the flat Helfrich geometry.

    Built on a ``side x side`` grid of ``n_side²`` vertices with periodic wrapping, so vertex areas
    are uniform and the modes are exactly orthonormal under the area-weighted product. Modes come in
    cosine/sine pairs at each admissible integer wavevector of the **half** plane ``kx > 0`` plus the
    ``kx == 0, ky > 0`` axis, ordered by ``|q|``. The half plane rather than a quadrant: ``(kx, -ky)``
    is the other diagonal orientation and an independent shape, so a quadrant would leave the
    instrument blind to one orientation at every off-axis ``|q|``.

    The normalisation is ``sqrt(2)`` on the trigonometric functions, which is what makes
    ``<phi_n, phi_n> = 1`` under ``sum_v A_v phi_v^2 / A_total`` — the mean of ``2*cos^2`` over a
    full period is 1.

    Args:
        side_um: Patch edge length [µm].
        n_side: Vertices per edge. The grid is ``n_side x n_side``.
        n_modes: How many modes to build, taken in ascending ``|q|``.
        include_beyond_nyquist: Build modes past ``pi / a`` as well. Defaults to ``False``, because
            a basis that cannot express an aliased mode cannot accidentally report one. Set it
            ``True`` only to exercise the exclusion machinery — which is exactly what
            ``test_modes_beyond_mesh_nyquist_are_excluded_with_a_reason`` does. Note that such a
            basis is **not orthonormal**: a mode labelled with ``q`` above the Nyquist limit samples
            to the same values as a lower-``q`` one, so the extra columns are duplicates wearing
            different labels. That is the aliasing being demonstrated, not a construction bug.
    """

    side_um: float
    n_side: int
    n_modes: int
    include_beyond_nyquist: bool = False

    _matrix: np.ndarray = field(init=False, repr=False)
    _q: np.ndarray = field(init=False, repr=False)
    _areas: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.side_um <= 0.0 or not math.isfinite(self.side_um):
            raise ValueError(f"side_um must be positive and finite; got {self.side_um!r}")
        if self.n_side < 2:
            raise ValueError(f"need at least a 2x2 grid; got n_side={self.n_side!r}")
        if self.n_modes < 1:
            raise ValueError(f"need at least one mode; got n_modes={self.n_modes!r}")

        n = int(self.n_side)
        side = float(self.side_um)
        spacing = side / n
        axis = (np.arange(n, dtype=np.float64) + 0.5) * spacing
        gx, gy = np.meshgrid(axis, axis, indexing="ij")
        x = gx.ravel()
        y = gy.ravel()

        q_unit = 2.0 * math.pi / side
        nyquist = math.pi / spacing
        # Integer wavevector limit. `n//2` is the grid's Nyquist index, and on a cell-centred grid
        # that index is degenerate — cos(pi*(j+1/2)) samples to exactly zero — so the honest default
        # stops one short of it. Building past it produces columns that are ALIASES of lower-k ones
        # rather than new modes, which is exactly why `include_beyond_nyquist` is off by default and
        # why a basis built with it on is deliberately not orthonormal.
        k_limit = (n - 1) if self.include_beyond_nyquist else max(1, n // 2 - 1)

        # Real modes are indexed by a HALF plane, not a quadrant. `-k` reproduces the same cos/sin
        # pair as `+k`, so one of each opposite pair is redundant — but `(kx, -ky)` is not the
        # opposite of `(kx, ky)`, it is the other diagonal orientation, and it is an independent
        # shape. Iterating the quadrant drops it at every off-axis |q|, which halves the mode count
        # on those shells and makes the missing orientation invisible rather than merely absent.
        # The half plane taken here is `kx > 0` at any `ky`, plus the `kx == 0` axis at `ky > 0`.
        candidates: list[tuple[float, int, int, str]] = []
        for kx in range(0, k_limit + 1):
            ky_low = 1 if kx == 0 else -k_limit
            for ky in range(ky_low, k_limit + 1):
                q = q_unit * math.hypot(kx, ky)
                if not self.include_beyond_nyquist and q > nyquist:
                    continue
                candidates.append((q, kx, ky, "cos"))
                candidates.append((q, kx, ky, "sin"))
        candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3]))

        columns: list[np.ndarray] = []
        q_values: list[float] = []
        for q, kx, ky, kind in candidates:
            if len(columns) == int(self.n_modes):
                break
            phase = q_unit * (kx * x + ky * y)
            wave = np.cos(phase) if kind == "cos" else np.sin(phase)
            # Normalise to exactly unit area-weighted norm. The analytic factor is sqrt(2) for a
            # generic mode, but a mode at the grid's Nyquist index samples to zero, so the norm is
            # measured. A degenerate column is skipped rather than rescaled: it is not representable
            # on this grid, which is a fact about the grid and not a failure.
            norm = math.sqrt(float(np.mean(wave * wave)))
            if norm <= 1e-12:
                continue
            columns.append(wave / norm)
            q_values.append(q)

        if len(columns) < int(self.n_modes):
            raise ValueError(
                f"requested {self.n_modes} modes but only {len(columns)} are representable on a "
                f"{n}x{n} grid under the current Nyquist policy; refuse to pad with duplicates"
            )

        object.__setattr__(self, "_matrix", np.column_stack(columns))
        object.__setattr__(self, "_q", np.asarray(q_values, dtype=np.float64))
        object.__setattr__(
            self, "_areas", np.full(n * n, spacing * spacing, dtype=np.float64)
        )

    @property
    def basis_matrix(self) -> np.ndarray:
        return self._matrix

    @property
    def wavevectors_inv_um(self) -> np.ndarray:
        return self._q

    @property
    def vertex_areas_um2(self) -> np.ndarray:
        return self._areas

    @property
    def total_area_um2(self) -> float:
        return float(np.sum(self._areas))

    @property
    def mean_vertex_spacing_um(self) -> float:
        return float(self.side_um) / float(self.n_side)


# ---------------------------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class TrustedRange:
    """The wavevector band a spectrum may be reported over, with both bounds' provenance.

    A range object rather than a pair of floats, because the *reasons* have to travel too: a reader
    deciding whether to fit over this band needs to know that the upper bound is aliasing and the
    lower bound is the box, since those fail differently.
    """

    q_min_inv_um: float
    q_max_inv_um: float
    q_min_basis: str
    q_max_basis: str

    def contains(self, q: float) -> bool:
        return self.q_min_inv_um <= float(q) <= self.q_max_inv_um

    def as_dict(self) -> dict[str, Any]:
        return {
            "q_min_inv_um": self.q_min_inv_um,
            "q_max_inv_um": self.q_max_inv_um,
            "q_min_basis": self.q_min_basis,
            "q_max_basis": self.q_max_basis,
        }


@dataclass(frozen=True, slots=True)
class ModeExclusion:
    """One mode that was measured and then not reported, with the reason it was not.

    Excluded modes are listed rather than dropped. A spectrum that quietly contains fewer modes than
    the basis it was computed from is a spectrum whose ``q`` axis nobody can reconstruct.
    """

    index: int
    q_inv_um: float
    reason: str


@dataclass(frozen=True, slots=True)
class SpectrumResult:
    """A mode power spectrum with its sample accounting, its trusted band, and its manifest hash.

    Attributes:
        q_inv_um: Wavevector magnitudes of the reported modes [1/µm], ascending.
        power_um2: ``P(q) = <|a_n|^2>`` for those modes [µm²].
        n_samples: Frames entering each mode's estimate.
        n_effective: Independent-sample count per mode, from ``tau_int``. Always at or below
            ``n_samples``, and much below it at small ``q``.
        standard_error_um2: Error on each ``P``, computed from ``n_effective``. For a variance
            estimated from ``m`` independent Gaussian samples the relative error is
            ``sqrt(2/(m-1))``, so ``se = P * sqrt(2 / (n_effective - 1))``.
        tau_int_s: Correlation time per mode [s]. Reported because its ``q``-dependence is a
            physical signal in its own right, not merely a correction factor.
        trusted_wavevector_range: The band the spectrum is claimed over.
        excluded_modes: Modes measured and not reported, with reasons.
        manifest_hash: Digest of the protocol manifest in force.
        n_frames: Frames in the input series.
        detail: Anything else on the record.
    """

    q_inv_um: np.ndarray
    power_um2: np.ndarray
    n_samples: np.ndarray
    n_effective: np.ndarray
    standard_error_um2: np.ndarray
    tau_int_s: np.ndarray
    trusted_wavevector_range: TrustedRange
    excluded_modes: tuple[ModeExclusion, ...]
    manifest_hash: str
    n_frames: int
    detail: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.trusted_wavevector_range, TrustedRange):
            raise TypeError(
                "a spectrum must carry its trusted wavevector range; a spectrum without one is a "
                "spectrum somebody will extrapolate past the mesh Nyquist limit"
            )
        if not self.manifest_hash:
            raise ValueError(
                "a spectrum must carry the hash of the protocol manifest that produced it"
            )
        n = self.q_inv_um.size
        widths = {
            "power_um2": self.power_um2.size,
            "n_samples": self.n_samples.size,
            "n_effective": self.n_effective.size,
            "standard_error_um2": self.standard_error_um2.size,
            "tau_int_s": self.tau_int_s.size,
        }
        bad = {name: size for name, size in widths.items() if size != n}
        if bad:
            raise ValueError(f"spectrum arrays disagree with {n} wavevectors: {bad}")

    def as_dict(self) -> dict[str, Any]:
        """JSON-able record. The spectrum is a field and is kept in full."""
        return {
            "q_inv_um": self.q_inv_um.tolist(),
            "power_um2": self.power_um2.tolist(),
            "n_samples": self.n_samples.tolist(),
            "n_effective": self.n_effective.tolist(),
            "standard_error_um2": self.standard_error_um2.tolist(),
            "tau_int_s": self.tau_int_s.tolist(),
            "trusted_wavevector_range": self.trusted_wavevector_range.as_dict(),
            "excluded_modes": [
                {"index": e.index, "q_inv_um": e.q_inv_um, "reason": e.reason}
                for e in self.excluded_modes
            ],
            "manifest_hash": self.manifest_hash,
            "n_frames": self.n_frames,
            "detail": dict(self.detail),
        }


# ---------------------------------------------------------------------------------------------
# The operator
# ---------------------------------------------------------------------------------------------


def project_onto_modes(
    displacements_um: np.ndarray, basis: ModeBasis, *, area_weighted: bool = True
) -> np.ndarray:
    """Project a ``(T, V)`` displacement series onto the mode basis. Returns ``(T, M)`` in µm.

    Args:
        displacements_um: Normal displacements per frame per vertex [µm].
        basis: The mode basis.
        area_weighted: Use the area-weighted inner product of §1. ``False`` exists solely so that
            ``test_unweighted_projection_biases_a_nonuniform_mesh`` can demonstrate the failure the
            weighting prevents. It is never the right choice for a measurement.

    Returns:
        Mode amplitudes ``a_n(t)`` [µm].
    """
    h = np.asarray(displacements_um, dtype=np.float64)
    phi = np.asarray(basis.basis_matrix, dtype=np.float64)
    if area_weighted:
        areas = np.asarray(basis.vertex_areas_um2, dtype=np.float64)
        total = float(basis.total_area_um2)
        return (h * areas) @ phi / total
    return h @ phi / float(phi.shape[0])


@dataclass(frozen=True, slots=True)
class FluctuationSpectrumOperator:
    """Passive thermal fluctuation spectrum of a membrane. Satisfies ``ObservationOperator``.

    Args:
        manifest: The protocol in force. Build it with
            :func:`~aleph.observe.manifest.passive_fluctuation_manifest`.
        basis: The mode basis. Any object satisfying :class:`ModeBasis`.
        owner: State owner holding the displacement series. Default ``"membrane"``.
        array_name: Array within that owner. Default ``"normal_displacement_um"``, shape ``(T, V)``.
        min_temperature_k: Lower applicability bound. Below the freezing point of the medium the
            fluctuation picture does not describe a fluid membrane at all.
        max_temperature_k: Upper applicability bound.
    """

    manifest: ProtocolManifest
    basis: ModeBasis
    owner: str = "membrane"
    array_name: str = "normal_displacement_um"
    min_temperature_k: float = 250.0
    max_temperature_k: float = 330.0

    @property
    def operator_id(self) -> str:
        return _OPERATOR_ID

    @property
    def required_state_owners(self) -> tuple[str, ...]:
        return (self.owner,)

    @property
    def identifiability_directions(self) -> tuple[str, ...]:
        """What this observation can actually constrain — including where it cannot.

        The honest statement has three parts, and the third is the one that gets forgotten:

        * ``sigma`` is constrained where ``q << q* = sqrt(sigma/kappa)``, the band in which the
          spectrum goes as ``q^-2``.
        * ``kappa`` is constrained where ``q >> q*``, where it goes as ``q^-4``.
        * **Nothing** is constrained along an overall amplitude rescale: multiplying every ``P(q)``
          by ``f`` is absorbed exactly by ``(kappa, sigma) -> (kappa/f, sigma/f)``. A calibration
          error in the displacement scale is therefore invisible to this observation, and it is
          invisible in a way that moves both parameters together.

        The list is a declaration, not a result. It is here so that a degeneracy is known before a
        posterior fails to contract rather than discovered afterwards.
        """
        return (
            "sigma: constrained on the q^-2 branch, q << sqrt(sigma/kappa)",
            "kappa: constrained on the q^-4 branch, q >> sqrt(sigma/kappa)",
            "NOT constrained: overall amplitude rescale, exactly degenerate with "
            "(kappa, sigma) -> (kappa/f, sigma/f)",
        )

    @property
    def applicability(self) -> ApplicabilityDomain:
        return ApplicabilityDomain(
            description=(
                "passive, unforced membrane in a fluid medium, sampled long enough for the slowest "
                "reported mode to have relaxed; no probe in contact"
            ),
            bounds={
                "temperature_k": (self.min_temperature_k, self.max_temperature_k),
                "sample_interval_s": (0.0, math.inf),
            },
            required_state_owners=(self.owner,),
        )

    def trusted_range(self) -> TrustedRange:
        """The band this basis and mesh can support, derived in §3. No caller input."""
        q_max = mesh_nyquist_wavevector(self.basis.mean_vertex_spacing_um)
        q_min = patch_fundamental_wavevector(self.basis.total_area_um2)
        return TrustedRange(
            q_min_inv_um=q_min,
            q_max_inv_um=q_max,
            q_min_basis=(
                "2*pi/sqrt(A_total): the longest wavelength that fits in the patch. Below it no "
                "full period is contained in the domain and the variance is set by the boundary "
                "condition rather than by the fluctuation spectrum"
            ),
            q_max_basis=(
                "pi/a with a the mean vertex spacing: the mesh Nyquist limit. Above it a mode is "
                "aliased and any power reported there belongs to a different mode"
            ),
        )

    # -- the two halves, kept apart ------------------------------------------------------------

    def raw_observable(
        self, state: AcceptedState, context: ObservationContext
    ) -> RawObservable | Refusal:
        """Project the displacement series and form the mode power spectrum.

        Refuses, in this order and before any arithmetic that could produce a number: missing state
        owner, out-of-domain context, missing array, wrong shape, too few frames, mismatched vertex
        count, non-finite samples, and finally an empty trusted band.
        """
        domain = self.applicability

        owner_refusal = domain.check_owners(self.operator_id, state.owners)
        if owner_refusal is not None:
            return owner_refusal

        bounds_refusal = domain.check(
            self.operator_id,
            temperature_k=context.temperature_k,
            sample_interval_s=context.sample_interval_s,
        )
        if bounds_refusal is not None:
            return bounds_refusal

        if not (context.sample_interval_s > 0.0):
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"sample_interval_s = {context.sample_interval_s!r}; a correlation time cannot "
                    "be expressed without a positive sampling interval, and defaulting it would "
                    "silently return tau in frames"
                ),
                operator_id=self.operator_id,
                detail={"sample_interval_s": context.sample_interval_s},
            )

        try:
            raw_series = state.array(self.owner, self.array_name)
        except KeyError as exc:
            return Refusal(
                code=RefusalCode.MISSING_STATE_OWNER,
                reason=str(exc),
                operator_id=self.operator_id,
                detail={"owner": self.owner, "array": self.array_name},
            )

        h = np.asarray(raw_series, dtype=np.float64)
        if h.ndim != 2:
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"displacement series must be (frames, vertices); got shape {h.shape}. A 1-D "
                    "array is one frame and has no ensemble to average over"
                ),
                operator_id=self.operator_id,
                detail={"shape": list(h.shape)},
            )

        n_frames, n_vertices = h.shape
        expected_vertices = int(np.asarray(self.basis.basis_matrix).shape[0])
        if n_vertices != expected_vertices:
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"displacement series has {n_vertices} vertices but the mode basis is built on "
                    f"{expected_vertices}; the projection would be meaningless"
                ),
                operator_id=self.operator_id,
                detail={"series_vertices": n_vertices, "basis_vertices": expected_vertices},
            )

        if n_frames < MIN_FRAMES:
            return Refusal(
                code=RefusalCode.INSUFFICIENT_SAMPLES,
                reason=(
                    f"{n_frames} frame(s); a power spectrum is an ensemble average and needs at "
                    f"least {MIN_FRAMES}. One snapshot has a shape, not a spectrum"
                ),
                operator_id=self.operator_id,
                detail={"n_frames": int(n_frames), "minimum": MIN_FRAMES},
            )

        if not np.all(np.isfinite(h)):
            n_bad = int(np.count_nonzero(~np.isfinite(h)))
            return Refusal(
                code=RefusalCode.NON_FINITE_INPUT,
                reason=(
                    f"{n_bad} non-finite displacement sample(s). Refused rather than masked: "
                    "dropping frames changes the sampling interval that every correlation time "
                    "here is expressed in"
                ),
                operator_id=self.operator_id,
                detail={"n_non_finite": n_bad},
            )

        band = self.trusted_range()
        q_all = np.asarray(self.basis.wavevectors_inv_um, dtype=np.float64)

        # The lower edge is compared against a round-off-widened copy of itself; see
        # BAND_EDGE_EPS_PER_VERTEX. `band` keeps the un-widened bound, because the reported band is
        # the physical one and the slack is an artefact of how the two sides were computed.
        n_vertices = int(np.asarray(self.basis.vertex_areas_um2).size)
        q_min_floor = band.q_min_inv_um * (
            1.0 - BAND_EDGE_EPS_PER_VERTEX * float(n_vertices) * EPS64
        )

        exclusions: list[ModeExclusion] = []
        keep: list[int] = []
        for index, q in enumerate(q_all):
            if q > band.q_max_inv_um:
                exclusions.append(
                    ModeExclusion(
                        index=index,
                        q_inv_um=float(q),
                        reason=(
                            f"q = {q:.6g} 1/µm is above the mesh Nyquist limit "
                            f"{band.q_max_inv_um:.6g} 1/µm; the power here belongs to an alias, "
                            "not to this mode"
                        ),
                    )
                )
            elif q < q_min_floor:
                exclusions.append(
                    ModeExclusion(
                        index=index,
                        q_inv_um=float(q),
                        reason=(
                            f"q = {q:.6g} 1/µm is below the patch fundamental "
                            f"{band.q_min_inv_um:.6g} 1/µm; no full period fits in the domain"
                        ),
                    )
                )
            else:
                keep.append(index)

        if not keep:
            return Refusal(
                code=RefusalCode.BEYOND_RESOLUTION_LIMIT,
                reason=(
                    f"every one of the {q_all.size} basis modes falls outside the trusted band "
                    f"[{band.q_min_inv_um:.6g}, {band.q_max_inv_um:.6g}] 1/µm. There is nothing "
                    "this mesh can report; reporting the modes anyway would be reporting aliases"
                ),
                operator_id=self.operator_id,
                detail={
                    "n_modes": int(q_all.size),
                    "trusted_range": band.as_dict(),
                    "excluded": [e.index for e in exclusions],
                },
            )

        amplitudes = project_onto_modes(h, self.basis)  # (T, M) in µm
        kept = np.asarray(keep, dtype=np.intp)
        a = amplitudes[:, kept]

        # Mean subtracted per mode: a static offset is the equilibrium shape, not a fluctuation.
        power = np.var(a, axis=0, ddof=1)

        n_samples = np.full(kept.size, int(n_frames), dtype=np.int64)
        n_effective = np.empty(kept.size, dtype=np.float64)
        tau_int = np.empty(kept.size, dtype=np.float64)
        for column in range(kept.size):
            # tau_int of the SQUARED amplitude: the spectrum is an average of a_n^2, so the
            # correlation time that governs the error on that average is the one of a_n^2, not of
            # a_n. For an Ornstein-Uhlenbeck a_n the two differ by a factor of two, and using the
            # amplitude's own tau would overstate the independent-sample count by that factor.
            estimate = integrated_autocorrelation_time(
                a[:, column] ** 2, context.sample_interval_s
            )
            n_effective[column] = estimate.n_effective
            tau_int[column] = estimate.tau_int_s

        # Relative error of a variance from m independent Gaussian samples is sqrt(2/(m-1)).
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.sqrt(2.0 / np.maximum(n_effective - 1.0, 1e-12))
        standard_error = power * rel

        order = np.argsort(np.asarray(self.basis.wavevectors_inv_um)[kept], kind="stable")

        result = SpectrumResult(
            q_inv_um=np.asarray(self.basis.wavevectors_inv_um)[kept][order],
            power_um2=power[order],
            n_samples=n_samples[order],
            n_effective=n_effective[order],
            standard_error_um2=standard_error[order],
            tau_int_s=tau_int[order],
            trusted_wavevector_range=band,
            excluded_modes=tuple(exclusions),
            manifest_hash=self.manifest.manifest_hash(),
            n_frames=int(n_frames),
            detail={
                "projection": "area-weighted, sum_v A_v phi_n(v) h_v / A_total",
                "mean_subtracted": True,
                "ddof": 1,
                "n_effective_basis": (
                    "tau_int of the squared amplitude series; the spectrum averages a_n^2, so that "
                    "is the correlation time governing its error"
                ),
                "temperature_k": float(context.temperature_k),
                "sample_interval_s": float(context.sample_interval_s),
            },
        )

        return RawObservable(
            name="membrane_mode_power_spectrum",
            operator_id=self.operator_id,
            values=result.power_um2,
            units="um^2",
            n_samples=result.n_effective,
            manifest_hash=result.manifest_hash,
            support={"q_inv_um": result.q_inv_um},
            detail={"spectrum": result},
        )

    def apparent_quantity(self, raw: RawObservable) -> ApparentQuantity | Refusal:
        """Report the mean mode energy in units of ``kB*T/2`` — the equipartition ratio.

        On the exact taut-string case this is an identity: every mode carries ``kB*T/2``, so the
        ratio is 1 with no fit, no threshold and no free parameter. That makes it the one apparent
        quantity from this operator that can be *unambiguously* wrong, which is why it is the
        default. See §5.

        The energy of mode ``n`` is ``E_n = (1/2) * K_n * <|a_n|^2>`` with ``K_n`` the mode
        stiffness, so the ratio is ``K_n * <|a_n|^2> / (kB*T)``. The stiffnesses come from the raw
        observable's detail, placed there by whoever built the state; this operator does not invent
        them, and refuses if they are absent.
        """
        spectrum = raw.detail.get("spectrum")
        if not isinstance(spectrum, SpectrumResult):
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason="raw observable does not carry a SpectrumResult; it was not produced here",
                operator_id=self.operator_id,
                detail={"raw_name": raw.name},
            )

        stiffness = raw.detail.get("mode_stiffness_pn_per_um")
        kbt = raw.detail.get("kbt_pn_um")
        if stiffness is None or kbt is None:
            return Refusal(
                code=RefusalCode.UNSUPPORTED_QUANTITY,
                reason=(
                    "the equipartition ratio needs the mode stiffnesses K_n [pN/µm] and kB*T "
                    "[pN·µm], and neither is derivable from the spectrum alone. Supply them on the "
                    "raw observable's detail, or ask a different operator"
                ),
                operator_id=self.operator_id,
                detail={"available": sorted(raw.detail)},
            )

        k = np.asarray(stiffness, dtype=np.float64)
        if k.shape != spectrum.power_um2.shape:
            return Refusal(
                code=RefusalCode.INCONSISTENT_INPUT,
                reason=(
                    f"mode stiffness array has shape {k.shape} against {spectrum.power_um2.shape} "
                    "reported modes"
                ),
                operator_id=self.operator_id,
                detail={},
            )

        ratios = k * spectrum.power_um2 / float(kbt)
        value = float(np.mean(ratios))

        # The error on the mean ratio uses the per-mode effective sample counts, not the frame
        # count and not the mode count. Modes are independent by construction (the basis is
        # orthonormal), so their errors combine in quadrature.
        per_mode_rel = np.sqrt(2.0 / np.maximum(spectrum.n_effective - 1.0, 1e-12))
        per_mode_err = ratios * per_mode_rel
        standard_error = float(np.sqrt(np.sum(per_mode_err**2)) / ratios.size)

        return ApparentQuantity(
            name="equipartition_ratio",
            operator_id=self.operator_id,
            value=value,
            units="dimensionless (mode energy in units of kB*T/2, so 1.0 is exact equipartition)",
            standard_error=standard_error,
            n_effective=float(np.sum(spectrum.n_effective)),
            manifest_hash=spectrum.manifest_hash,
            analysis_chain=(
                "area-weighted projection of vertex displacements onto the mode basis",
                "per-mode mean subtraction",
                "per-mode variance with ddof=1",
                "restriction to the trusted wavevector band [2*pi/sqrt(A), pi/a]",
                "per-mode effective sample count from tau_int of the squared amplitude",
                "energy E_n = (1/2) K_n <|a_n|^2>, reported as E_n / (kB*T/2)",
                "mean over reported modes, errors combined in quadrature",
            ),
            identifiability=self.identifiability_directions,
            detail={
                "n_modes_reported": int(ratios.size),
                "n_modes_excluded": len(spectrum.excluded_modes),
                "trusted_wavevector_range": spectrum.trusted_wavevector_range.as_dict(),
                "per_mode_ratio": ratios.tolist(),
            },
        )


def attach_equipartition_inputs(
    raw: RawObservable, *, mode_stiffness_pn_per_um: Sequence[float] | np.ndarray, kbt_pn_um: float
) -> RawObservable:
    """Return a copy of ``raw`` carrying the stiffnesses and ``kB*T`` the apparent quantity needs.

    A separate function on purpose. These are properties of the *model*, not of the measurement, and
    letting the operator reach into the state for them would let it choose a stiffness that makes
    the ratio come out at 1.
    """
    detail = dict(raw.detail)
    detail["mode_stiffness_pn_per_um"] = np.asarray(mode_stiffness_pn_per_um, dtype=np.float64)
    detail["kbt_pn_um"] = float(kbt_pn_um)
    return RawObservable(
        name=raw.name,
        operator_id=raw.operator_id,
        values=raw.values,
        units=raw.units,
        n_samples=raw.n_samples,
        manifest_hash=raw.manifest_hash,
        support=raw.support,
        detail=detail,
    )
