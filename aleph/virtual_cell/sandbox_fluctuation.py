"""S1 — the passive thermal-fluctuation canon, and the machinery that judges a spectrum against it.

**This module makes no claim about the engine's membrane.**  It never imports Warp, never reads a run
record, never touches ``aleph.ac``, and it has not measured a single amplitude produced by the
Warp/CUDA membrane shell.  It is the *oracle the engine will be judged against* once an accepted
native fluctuation trace exists; everything it emits is
:data:`~aleph.virtual_cell.contracts.EvidenceSource.ANALYTIC_ORACLE`, and nothing here is native
evidence or may close a physics gate.

The sandbox question (master plan §8, row S1) is::

    does the membrane/FEM thermal fluctuation match the analytic spectrum, MODE BY MODE?

and its declared falsifier is a per-mode spectrum mismatch, which opens a **stochastic boundary law**
project.  The hard part is not writing ``kT / (A (kappa q**4 + sigma q**2))``.  The hard part is that a
mode-by-mode check is *defeatable*, in four specific ways this module measures rather than asserts:

1. **A two-parameter fit hides a wrong exponent.**  Over a narrow band of ``q`` a pure ``q**-2``
   spectrum is fitted by the Helfrich form with ``kappa -> 0`` at excellent chi-square.  Fit quality is
   therefore not the test; the *exponent*, estimated over a stated decade span, is.
   :func:`decades_required_for_exponent_separation` states how many decades of ``q`` are needed before
   ``q**-2`` and ``q**-4`` are separable at a declared confidence, and
   :func:`measure_exponent_separation_ladder` measures it by sampling.  That number is a requirement
   levied on a future native run: it fixes how large a patch must be simulated.
2. **Equipartition is only falsifiable when ``kappa`` and ``sigma`` are known independently.**  With
   both free, scaling every amplitude by ``sqrt(f)`` is *exactly* absorbed by
   ``(kappa, sigma) -> (kappa/f, sigma/f)``; the residuals do not move by one part in ``1e-12``.
   :data:`EQUIPARTITION_SCALE_DEGENERACY` names this and :func:`check_equipartition` refuses to run
   without an independently supplied spectrum.
3. **A per-mode variance check passes correlated modes.**  Independence is a separate measurement
   (:func:`check_mode_independence`), not a corollary of the variances being right.
4. **A mesh cutoff is not an absence of physics.**  A spectrum truncated at high ``q`` by the grid
   looks like a stiffer membrane if it is fitted.  :func:`trusted_wavevector_range` finds the largest
   low-``q`` prefix that passes a goodness-of-fit declared *before* the data, and the fit is taken
   there.  How MUCH stiffer is not a constant: it is a monotone, unbounded function of where the
   cutoff sits inside the fitted band, measured by :func:`measure_cutoff_damage_ladder` and scoped by
   :data:`CUTOFF_DAMAGE_IS_NOT_A_MULTIPLE`.  Any single multiple quoted from that curve is a property
   of one cutoff position and must be reported with the design that produced it.

The 1-D taut string is included because it is the one case with no approximation at all: with fixed
ends the sine coefficients are exactly independent quadratic degrees of freedom, so
``<E_n> = kB*T/2`` per mode is an identity, not a limit.  If the estimator cannot recover that, it is
broken, and no amount of agreement on the 2-D membrane would mean anything.

Provenance of the physical constants (honest, and inherited rather than verified):

* ``kB`` is the exact SI defining constant; ``T = 310 K`` is the physiological operating point.
* :data:`ENGINE_KAPPA_J` and :data:`ENGINE_SIGMA_N_PER_M` are read off this repository's own engine
  constants (``aleph/components/incumbent/compartments.py`` ``MEMBRANE_DEFAULTS``), tagged there ``KB-3.B1.2``
  (attributed to Rawicz 2000) and ``KB-3.B1.1`` (attributed to Diz-Munoz 2013).  **This lane has read
  neither paper.**  The provenance recorded is "repository constant with an inherited citation", which
  is why :data:`ENGINE_MEMBRANE_PROVENANCE` states it that way and why every public entry point takes
  ``kappa`` and ``sigma`` as arguments instead of defaulting to them.

Sanity Gate (recorded before first execution):

* **Dimensional** — enforced, not asserted in prose.  :class:`Dimension` is a four-exponent
  (length, mass, time, temperature) algebra and :func:`helfrich_variance_dimension` composes
  ``kB*T / (A * (kappa*q**4 + sigma*q**2))`` symbolically; it must equal ``length**2`` and it *raises*
  if the two terms of the denominator are not dimensionally commensurate.  A ``q**3`` typo therefore
  fails at the dimension layer rather than surviving as a plausible-looking fit.  SI throughout:
  ``kappa`` in J, ``sigma`` in N/m, ``q`` in 1/m, ``A`` in m**2, variance in m**2.  Engine units
  (pN, um) are converted at the boundary by :func:`kappa_from_pn_um` / :func:`sigma_from_pn_per_um`,
  never mixed inside a formula.
* **Boundary** — ``q = 0`` has infinite variance (the translation mode is not a fluctuation mode) and
  is refused, not silently regularised.  ``sigma = 0`` is legal (pure bending) and ``kappa = 0`` is
  legal (pure tension), but not both; a negative ``kappa`` or ``sigma`` raises.  The crossover
  ``q* = sqrt(sigma/kappa)`` is undefined for ``kappa = 0`` and returns ``inf``, which is the correct
  statement that the spectrum is tension-dominated everywhere.
* **Conservation / invariant** — equipartition itself is the invariant: mean energy per independent
  quadratic mode is ``kB*T/2`` exactly.  :func:`mean_energy_per_quadratic_dof` measures it in units of
  ``kB*T/2`` so the target is the literal number 1.  A fitted variance is a second moment and is
  non-negative by construction; the fit is constrained to ``kappa >= 0, sigma >= 0`` because a negative
  bending modulus is not a smaller answer, it is an unstable membrane.
* **Measurement protocol** — every estimator is told ``n_samples`` and uses the *known* sampling law of
  a variance estimator rather than an empirical scatter: for ``N`` zero-mean real Gaussian draws,
  ``S_hat = mean(a**2)`` is ``V * chi2_N / N``, so ``log S_hat`` has variance ``2/N`` independent of
  ``q``.  That is what makes an unweighted log-space fit correct here, and it is what makes the
  chi-square goodness-of-fit meaningful.  Confidence levels, the ``n_sigma`` separation requirement,
  and the Bonferroni family size are all arguments fixed by the caller before the data, never
  re-thresholded afterwards.
* **Numerical** — the fit is Gauss-Newton on log residuals (the model is linear in ``(kappa, sigma)``
  *inside* the log), seeded from a weighted linear solve; it is declared non-converged rather than
  returning the last iterate.  Parameter covariance uses the a-priori ``2/N`` log variance, not the
  residual scatter, so a *wrong model* inflates chi-square instead of quietly inflating the error bars
  until it fits.  All spectra are strictly positive, so no logarithm is taken of a non-positive number
  without raising first.

Known and measured limitation, stated rather than hidden: ``kappa`` and ``sigma`` are constrained to
be non-negative, so when the truth sits ON that boundary (a pure ``q**-2`` spectrum has ``kappa = 0``)
the fit consumes fewer than two effective parameters while :attr:`SpectrumFit.dof` still reports
``M - 2``.  The goodness-of-fit is then CONSERVATIVE — it over-rejects.  Measured: for a ``q**-2``
spectrum at ``N = 1000`` the chi-square rejection rate at ``alpha = 0.01`` is 0.130 with 7 modes but
falls to the nominal 0.010 by 21 modes.  The error is in the safe direction and it disappears with
mode count, but a 7-mode fit's p-value is not a calibrated p-value and must not be quoted as one.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from statistics import NormalDist
from typing import Any

import numpy as np
from scipy.special import gammaincc

from aleph.virtual_cell.contracts import EvidenceSource, SandboxExperimentCard

__all__ = [
    "BOLTZMANN_J_PER_K",
    "CUTOFF_DAMAGE_IS_NOT_A_MULTIPLE",
    "EQUIPARTITION_SCALE_DEGENERACY",
    "EVIDENCE_AUTHORITY",
    "EVIDENCE_SOURCE",
    "ENGINE_KAPPA_J",
    "ENGINE_MEMBRANE_PROVENANCE",
    "ENGINE_SIGMA_N_PER_M",
    "HELFRICH_FAMILY_ABSORPTION",
    "PATCH_BOUNDARY_CONDITION_IS_A_DECLARATION",
    "PATCH_PERIODIC_RATIONALE",
    "PHYSIOLOGICAL_TEMPERATURE_K",
    "CutoffDamageLadder",
    "Dimension",
    "EquipartitionReport",
    "ExponentFit",
    "ExponentSeparationLadder",
    "FluctuationFailure",
    "FluctuationRegime",
    "FluctuationVerdict",
    "HelfrichSpectrum",
    "ModeIndependenceReport",
    "PatchBoundaryCondition",
    "PatchDesign",
    "SpectrumFit",
    "TautStringSpectrum",
    "TrustedRangeReport",
    "adjudicate_fluctuation_spectrum",
    "chi_square_survival",
    "correlated_mode_amplitudes",
    "decades_required_for_exponent_separation",
    "fit_helfrich_spectrum",
    "fit_power_law_exponent",
    "helfrich_variance_dimension",
    "kappa_from_pn_um",
    "log_spaced_wavevectors",
    "mean_energy_per_quadratic_dof",
    "measure_cutoff_damage_ladder",
    "measure_exponent_separation_ladder",
    "measure_mode_variance",
    "mesh_nyquist_wavevector",
    "mesh_truncated_variances",
    "native_patch_requirement",
    "power_law_variances",
    "s1_experiment_card",
    "sample_gaussian_modes",
    "sigma_from_pn_per_um",
    "check_equipartition",
    "check_mode_independence",
    "thermal_energy",
    "trusted_wavevector_range",
]

# --------------------------------------------------------------------------------------------
# constants and provenance
# --------------------------------------------------------------------------------------------

#: Boltzmann constant, J/K.  Exact by the 2019 SI redefinition, so this is a definition and not a
#: measurement with provenance.
BOLTZMANN_J_PER_K = 1.380649e-23

#: Physiological operating point.  Every spectrum in this module is evaluated at a real temperature;
#: a "unit temperature" convenience would be a non-physical baseline.
PHYSIOLOGICAL_TEMPERATURE_K = 310.0

#: 1 pN*um in joules: 1e-12 N * 1e-6 m.
_JOULE_PER_PN_UM = 1.0e-18
#: 1 pN/um in N/m: 1e-12 N / 1e-6 m.
_N_PER_M_PER_PN_PER_UM = 1.0e-6

#: Helfrich bending modulus taken from this repository's engine constant
#: ``MEMBRANE_DEFAULTS["kappa_m"] = 0.0828 pN*um``.  See :data:`ENGINE_MEMBRANE_PROVENANCE` —
#: the citation attached to it in the tree was NOT read by this lane.
ENGINE_KAPPA_J = 0.0828 * _JOULE_PER_PN_UM
#: Bilayer in-plane tension from ``MEMBRANE_DEFAULTS["gamma_mem"] = 10.0 pN/um``.
ENGINE_SIGMA_N_PER_M = 10.0 * _N_PER_M_PER_PN_PER_UM

#: How the two numbers above are allowed to be described in a deliverable.
ENGINE_MEMBRANE_PROVENANCE: dict[str, str] = {
    "kappa_source": "repository constant aleph/components/incumbent/compartments.py MEMBRANE_DEFAULTS['kappa_m']",
    "kappa_inherited_citation": "KB-3.B1.2 (attributed in-tree to Rawicz 2000) — NOT READ BY THIS LANE",
    "kappa_engine_value": "0.0828 pN*um",
    "kappa_temperature_note": (
        "the in-tree comment calls 0.0828 pN*um '20 kBT'; 20*kB*T is 0.0828 pN*um at T = 300.0 K and "
        "0.0856 pN*um at the physiological 310 K, so the engine constant carries a 300 K convention"
    ),
    "sigma_source": "repository constant aleph/components/incumbent/compartments.py MEMBRANE_DEFAULTS['gamma_mem']",
    "sigma_inherited_citation": (
        "KB-3.B1.1 (attributed in-tree to Diz-Munoz 2013) — NOT READ BY THIS LANE"
    ),
    "sigma_engine_value": "10.0 pN/um",
    "status": (
        "INHERITED-UNVERIFIED: usable as a parameterisation for oracle work, not quotable as a "
        "sourced literature value by this lane; treat as PI_GAP for any claim about a real membrane"
    ),
}

#: Measured structural fact, not an opinion: a pure ``q**-2`` spectrum is a MEMBER of the Helfrich
#: family (it is the ``kappa = 0`` case), so the two-parameter fit's goodness-of-fit can never reject
#: it — measured chi-square rejection rate stays at the nominal alpha out to a four-decade span.  Only
#: the exponent test, or the fit's own report that ``kappa`` is unresolved, catches "tension mistaken
#: for bending".  A spectrum OUTSIDE the family (``q**-3``) is rejected by chi-square, but needs a
#: much wider span than the exponent test does.
HELFRICH_FAMILY_ABSORPTION = (
    "q**-2 and q**-4 are both members of the Helfrich two-parameter family, so no chi-square "
    "goodness-of-fit on that family can discriminate them at any span.  Report the EXPONENT interval "
    "and the kappa-resolved flag; never report fit quality as evidence of bending."
)

#: The structural reason an equipartition violation is invisible to a free two-parameter fit.
EQUIPARTITION_SCALE_DEGENERACY = (
    "scaling every mode amplitude by sqrt(f) multiplies every mode variance by f, which the Helfrich "
    "form absorbs exactly via (kappa, sigma) -> (kappa/f, sigma/f).  The log-space residuals are "
    "invariant, so chi-square cannot see it.  Equipartition is falsifiable ONLY against an "
    "independently determined kappa and sigma (or an independently determined temperature)."
)

#: What this layer may claim.  Mirrors the vocabulary used by the rest of the virtual-cell lanes.
EVIDENCE_AUTHORITY = "unverified-observer"
#: Everything this module produces is an analytic oracle number.
EVIDENCE_SOURCE = EvidenceSource.ANALYTIC_ORACLE


# --------------------------------------------------------------------------------------------
# dimensional algebra — the wrong-power-of-q trap, caught structurally
# --------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Dimension:
    """A physical dimension as exponents of (length, mass, time, temperature).

    This exists for one reason: ``<|h_q|^2> = kB*T / (A * (kappa*q**4 + sigma*q**2))`` is a formula in
    which a wrong power of ``q`` produces a curve that still looks like a spectrum and still admits a
    two-parameter fit.  Composing the formula in this algebra makes that mistake a raised exception.
    """

    length: float = 0.0
    mass: float = 0.0
    time: float = 0.0
    temperature: float = 0.0

    def __mul__(self, other: Dimension) -> Dimension:
        return Dimension(
            self.length + other.length,
            self.mass + other.mass,
            self.time + other.time,
            self.temperature + other.temperature,
        )

    def __truediv__(self, other: Dimension) -> Dimension:
        return Dimension(
            self.length - other.length,
            self.mass - other.mass,
            self.time - other.time,
            self.temperature - other.temperature,
        )

    def __pow__(self, exponent: float) -> Dimension:
        return Dimension(
            self.length * exponent,
            self.mass * exponent,
            self.time * exponent,
            self.temperature * exponent,
        )

    def __str__(self) -> str:
        parts = []
        for name, value in (
            ("L", self.length),
            ("M", self.mass),
            ("T", self.time),
            ("K", self.temperature),
        ):
            if value:
                parts.append(f"{name}^{value:g}")
        return "*".join(parts) if parts else "1"


_DIMENSIONLESS = Dimension()
_LENGTH = Dimension(length=1.0)
_MASS = Dimension(mass=1.0)
_TIME = Dimension(time=1.0)
_TEMPERATURE = Dimension(temperature=1.0)
_ENERGY = _MASS * _LENGTH**2 / _TIME**2
_FORCE = _MASS * _LENGTH / _TIME**2
_BOLTZMANN_DIM = _ENERGY / _TEMPERATURE
_WAVEVECTOR = _LENGTH**-1
_AREA = _LENGTH**2
_KAPPA_DIM = _ENERGY  # J
_SIGMA_DIM = _FORCE / _LENGTH  # N/m


def helfrich_variance_dimension(*, bending_power: int = 4, tension_power: int = 2) -> Dimension:
    """Compose the Helfrich mode-variance formula symbolically and return its dimension.

    Args:
        bending_power: exponent of ``q`` multiplying ``kappa``.  The physical value is ``4``.
        tension_power: exponent of ``q`` multiplying ``sigma``.  The physical value is ``2``.

    Returns:
        The dimension of ``kB*T / (A * (kappa*q**b + sigma*q**t))``.

    Raises:
        ValueError: if the two denominator terms are not dimensionally commensurate, which is exactly
            what happens for any ``(bending_power, tension_power)`` other than ``(4, 2)`` up to a
            common shift, or if the composed result is not an area.
    """
    bending_term = _KAPPA_DIM * _WAVEVECTOR**bending_power
    tension_term = _SIGMA_DIM * _WAVEVECTOR**tension_power
    if bending_term != tension_term:
        raise ValueError(
            f"Helfrich denominator terms are not commensurate: kappa*q^{bending_power} is "
            f"[{bending_term}] but sigma*q^{tension_power} is [{tension_term}]"
        )
    result = (_BOLTZMANN_DIM * _TEMPERATURE) / (_AREA * bending_term)
    if result != _AREA:
        raise ValueError(f"Helfrich mode variance must be an area; got [{result}]")
    return result


# --------------------------------------------------------------------------------------------
# small numerical helpers
# --------------------------------------------------------------------------------------------


def _positive_float(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{what} must be finite and strictly positive; got {value!r}")
    return number


def _nonnegative_float(value: object, *, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a real number")
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{what} must be finite and non-negative; got {value!r}")
    return number


def _positive_int(value: object, *, what: str, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{what} must be an int")
    if value < minimum:
        raise ValueError(f"{what} must be >= {minimum}; got {value}")
    return value


def _confidence(value: object) -> float:
    number = _positive_float(value, what="confidence")
    if not 0.0 < number < 1.0:
        raise ValueError(f"confidence must lie strictly in (0, 1); got {number}")
    return number


def _positive_array(values: Sequence[float] | np.ndarray, *, what: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"{what} must be a non-empty 1-D array")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{what} must be finite")
    if np.any(array <= 0.0):
        raise ValueError(f"{what} must be strictly positive (q = 0 is the translation mode)")
    return array


def chi_square_survival(statistic: float, dof: int) -> float:
    """Return ``P(chi2_dof > statistic)``.

    Args:
        statistic: observed chi-square, ``>= 0``.
        dof: degrees of freedom, ``>= 1``.

    Returns:
        Upper-tail probability.
    """
    value = _nonnegative_float(statistic, what="statistic")
    k = _positive_int(dof, what="dof")
    return float(gammaincc(0.5 * k, 0.5 * value))


def thermal_energy(temperature_K: float) -> float:
    """Return ``kB*T`` in joules.

    Args:
        temperature_K: absolute temperature in kelvin.

    Returns:
        Thermal energy scale in J.
    """
    return BOLTZMANN_J_PER_K * _positive_float(temperature_K, what="temperature_K")


def kappa_from_pn_um(kappa_pn_um: float) -> float:
    """Convert a bending modulus from the engine's ``pN*um`` to SI joules."""
    return _nonnegative_float(kappa_pn_um, what="kappa_pn_um") * _JOULE_PER_PN_UM


def sigma_from_pn_per_um(sigma_pn_per_um: float) -> float:
    """Convert a surface tension from the engine's ``pN/um`` to SI N/m."""
    return _nonnegative_float(sigma_pn_per_um, what="sigma_pn_per_um") * _N_PER_M_PER_PN_PER_UM


def mesh_nyquist_wavevector(spacing_m: float) -> float:
    """Largest wavevector a mesh of the given node spacing can represent, ``pi / a`` in 1/m.

    Args:
        spacing_m: node-to-node spacing in metres.

    Returns:
        Nyquist wavevector in 1/m.  Power reported above this is a discretisation artifact, never
        physics.
    """
    return math.pi / _positive_float(spacing_m, what="spacing_m")


def log_spaced_wavevectors(q_min: float, q_max: float, n_modes: int) -> np.ndarray:
    """Return ``n_modes`` logarithmically spaced wavevectors on ``[q_min, q_max]``.

    Log spacing is the right protocol for an exponent measurement: it makes the regressor
    ``log q`` uniform, so every mode contributes the same leverage.
    """
    low = _positive_float(q_min, what="q_min")
    high = _positive_float(q_max, what="q_max")
    count = _positive_int(n_modes, what="n_modes", minimum=2)
    if high <= low:
        raise ValueError("q_max must exceed q_min")
    return np.logspace(math.log10(low), math.log10(high), count)


# --------------------------------------------------------------------------------------------
# the analytic canon
# --------------------------------------------------------------------------------------------


class FluctuationRegime(StrEnum):
    """Which term of the Helfrich denominator dominates at a given wavevector."""

    BENDING_DOMINATED = "bending-dominated"
    TENSION_DOMINATED = "tension-dominated"
    CROSSOVER = "crossover"


@dataclass(frozen=True, slots=True)
class HelfrichSpectrum:
    """Equilibrium mode spectrum of a nearly-flat thermally fluctuating Helfrich patch.

    The mean square amplitude of the mode at wavevector ``q`` is

    ``<|h_q|^2> = kB*T / (A * (kappa*q**4 + sigma*q**2))``

    which is equipartition applied to the quadratic Helfrich energy
    ``E = (A/2) * sum_q (kappa*q**4 + sigma*q**2) |h_q|^2``.  It is an *equilibrium* statement: it
    constrains no dynamics, no relaxation rate and no viscosity, and a simulation can reproduce it
    exactly while getting the boundary law completely wrong — which is why S1's failure route is a
    stochastic boundary law project and not a "tune kappa" project.

    Attributes:
        kappa_J: bending modulus in joules.  May be ``0`` (pure tension).
        sigma_N_per_m: surface tension in N/m.  May be ``0`` (pure bending).
        area_m2: projected patch area in m**2.
        temperature_K: absolute temperature in kelvin.
    """

    kappa_J: float
    sigma_N_per_m: float
    area_m2: float
    temperature_K: float = PHYSIOLOGICAL_TEMPERATURE_K

    def __post_init__(self) -> None:
        object.__setattr__(self, "kappa_J", _nonnegative_float(self.kappa_J, what="kappa_J"))
        object.__setattr__(
            self, "sigma_N_per_m", _nonnegative_float(self.sigma_N_per_m, what="sigma_N_per_m")
        )
        object.__setattr__(self, "area_m2", _positive_float(self.area_m2, what="area_m2"))
        object.__setattr__(
            self, "temperature_K", _positive_float(self.temperature_K, what="temperature_K")
        )
        if self.kappa_J == 0.0 and self.sigma_N_per_m == 0.0:
            raise ValueError(
                "kappa and sigma cannot both be zero: a surface with no restoring term has no "
                "equilibrium spectrum, it has an unbounded one"
            )
        # Dimensional check runs on construction, not only in a test.
        helfrich_variance_dimension()

    @property
    def thermal_energy_J(self) -> float:
        """``kB*T`` in joules."""
        return thermal_energy(self.temperature_K)

    @property
    def crossover_wavevector_per_m(self) -> float:
        """``q* = sqrt(sigma/kappa)`` in 1/m — a checkable prediction, not a fitted quantity.

        Below ``q*`` the spectrum is tension-dominated and falls as ``q**-2``; above it, bending
        dominates and it falls as ``q**-4``.  Returns ``inf`` when ``kappa == 0`` (tension-dominated
        at every wavevector) and ``0.0`` when ``sigma == 0`` (bending-dominated everywhere).
        """
        if self.kappa_J == 0.0:
            return math.inf
        return math.sqrt(self.sigma_N_per_m / self.kappa_J)

    @property
    def crossover_wavelength_m(self) -> float:
        """``2*pi/q*`` in metres, i.e. the real-space length at which the two terms are equal."""
        q_star = self.crossover_wavevector_per_m
        if q_star == 0.0:
            return math.inf
        if math.isinf(q_star):
            return 0.0
        return 2.0 * math.pi / q_star

    def restoring_density(self, q: Sequence[float] | np.ndarray) -> np.ndarray:
        """Return ``kappa*q**4 + sigma*q**2`` in N/m**3 (energy per unit area per unit amplitude**2)."""
        wavevectors = _positive_array(q, what="q")
        return self.kappa_J * wavevectors**4 + self.sigma_N_per_m * wavevectors**2

    def mode_variance(self, q: Sequence[float] | np.ndarray) -> np.ndarray:
        """Return ``<|h_q|^2>`` in m**2 for each wavevector.

        Args:
            q: wavevectors in 1/m, strictly positive.  ``q = 0`` is the rigid translation of the
                patch, not a fluctuation mode, and is refused.

        Returns:
            Mean square amplitudes in m**2.
        """
        return self.thermal_energy_J / (self.area_m2 * self.restoring_density(q))

    def bending_limit_variance(self, q: Sequence[float] | np.ndarray) -> np.ndarray:
        """Pure-bending limit ``kB*T / (A*kappa*q**4)``; valid for ``q >> q*``."""
        if self.kappa_J == 0.0:
            raise ValueError("the bending limit is undefined for kappa = 0")
        wavevectors = _positive_array(q, what="q")
        return self.thermal_energy_J / (self.area_m2 * self.kappa_J * wavevectors**4)

    def tension_limit_variance(self, q: Sequence[float] | np.ndarray) -> np.ndarray:
        """Pure-tension limit ``kB*T / (A*sigma*q**2)``; valid for ``q << q*``."""
        if self.sigma_N_per_m == 0.0:
            raise ValueError("the tension limit is undefined for sigma = 0")
        wavevectors = _positive_array(q, what="q")
        return self.thermal_energy_J / (self.area_m2 * self.sigma_N_per_m * wavevectors**2)

    def regime_at(self, q: float, *, dominance_ratio: float = 10.0) -> FluctuationRegime:
        """Classify a single wavevector by which denominator term dominates.

        Args:
            q: wavevector in 1/m.
            dominance_ratio: how many times larger one term must be before the label is committed.
                The default ``10`` means a decade; anything less is honestly ``CROSSOVER``.

        Returns:
            The regime label.
        """
        wavevector = _positive_float(q, what="q")
        ratio_floor = _positive_float(dominance_ratio, what="dominance_ratio")
        if ratio_floor <= 1.0:
            raise ValueError("dominance_ratio must exceed 1")
        bending = self.kappa_J * wavevector**4
        tension = self.sigma_N_per_m * wavevector**2
        if tension == 0.0 or bending > ratio_floor * tension:
            return FluctuationRegime.BENDING_DOMINATED
        if bending == 0.0 or tension > ratio_floor * bending:
            return FluctuationRegime.TENSION_DOMINATED
        return FluctuationRegime.CROSSOVER

    def expected_exponent(self, q: float, *, dominance_ratio: float = 10.0) -> float | None:
        """Return the local log-log slope expected at ``q``, or ``None`` inside the crossover."""
        regime = self.regime_at(q, dominance_ratio=dominance_ratio)
        if regime is FluctuationRegime.BENDING_DOMINATED:
            return -4.0
        if regime is FluctuationRegime.TENSION_DOMINATED:
            return -2.0
        return None


@dataclass(frozen=True, slots=True)
class TautStringSpectrum:
    """1-D analogue: a taut string of length ``L`` with fixed ends at temperature ``T``.

    Written in the real sine basis ``h(x) = sum_n a_n sin(n*pi*x/L)``, the elastic energy
    ``(tau/2) * integral (dh/dx)^2 dx`` is exactly ``sum_n (tau*L/4) * q_n**2 * a_n**2``.  Each ``a_n``
    is therefore an independent quadratic degree of freedom with stiffness ``k_n = tau*L*q_n**2/2``,
    and equipartition gives ``<E_n> = kB*T/2`` and ``<a_n**2> = 2*kB*T/(tau*L*q_n**2)``.

    There is no approximation anywhere in that chain — no small-gradient expansion beyond the quadratic
    energy that defines the model, no continuum limit, no mode truncation.  It is the case where the
    estimator can be validated against an identity.
    """

    tension_N: float
    length_m: float
    temperature_K: float = PHYSIOLOGICAL_TEMPERATURE_K

    def __post_init__(self) -> None:
        object.__setattr__(self, "tension_N", _positive_float(self.tension_N, what="tension_N"))
        object.__setattr__(self, "length_m", _positive_float(self.length_m, what="length_m"))
        object.__setattr__(
            self, "temperature_K", _positive_float(self.temperature_K, what="temperature_K")
        )

    @property
    def thermal_energy_J(self) -> float:
        """``kB*T`` in joules."""
        return thermal_energy(self.temperature_K)

    def mode_wavevector(self, n: Sequence[int] | np.ndarray) -> np.ndarray:
        """Return ``q_n = n*pi/L`` in 1/m for mode indices ``n >= 1``."""
        indices = np.asarray(n, dtype=float)
        if indices.ndim != 1 or indices.size == 0:
            raise ValueError("mode indices must be a non-empty 1-D array")
        if np.any(indices < 1) or np.any(indices != np.floor(indices)):
            raise ValueError("mode indices must be integers >= 1")
        return indices * math.pi / self.length_m

    def mode_stiffness(self, n: Sequence[int] | np.ndarray) -> np.ndarray:
        """Return ``k_n = tau*L*q_n**2/2`` in N/m so that ``E_n = k_n*a_n**2/2``."""
        q = self.mode_wavevector(n)
        return 0.5 * self.tension_N * self.length_m * q**2

    def mode_variance(self, n: Sequence[int] | np.ndarray) -> np.ndarray:
        """Return ``<a_n**2> = kB*T / k_n = 2*kB*T/(tau*L*q_n**2)`` in m**2."""
        return self.thermal_energy_J / self.mode_stiffness(n)


def mean_energy_per_quadratic_dof(
    amplitudes: np.ndarray,
    stiffnesses: Sequence[float] | np.ndarray,
    *,
    temperature_K: float,
) -> np.ndarray:
    """Measure ``<E_n>`` per mode in units of ``kB*T/2``.

    The equipartition invariant is the literal number ``1`` for every mode.  Expressing the result in
    units of ``kB*T/2`` rather than in joules is deliberate: it removes any chance of a plausible-
    looking joule value passing inspection while being wrong by a factor of two.

    Args:
        amplitudes: array of shape ``(n_samples, n_modes)`` of real mode amplitudes in m.
        stiffnesses: per-mode stiffness ``k_n`` in N/m, so that ``E_n = k_n*a_n**2/2``.
        temperature_K: absolute temperature the amplitudes were drawn at.

    Returns:
        Array of shape ``(n_modes,)`` giving ``<k_n*a_n**2/2> / (kB*T/2)``.
    """
    samples = np.asarray(amplitudes, dtype=float)
    if samples.ndim != 2 or samples.size == 0:
        raise ValueError("amplitudes must be a non-empty (n_samples, n_modes) array")
    if not np.all(np.isfinite(samples)):
        raise ValueError("amplitudes must be finite")
    k = _positive_array(stiffnesses, what="stiffnesses")
    if k.size != samples.shape[1]:
        raise ValueError("stiffnesses must have one entry per mode column")
    mean_energy = 0.5 * k * np.mean(samples**2, axis=0)
    return mean_energy / (0.5 * thermal_energy(temperature_K))


# --------------------------------------------------------------------------------------------
# samplers: the positive control and the four negative controls
# --------------------------------------------------------------------------------------------


def sample_gaussian_modes(
    variances: Sequence[float] | np.ndarray,
    *,
    n_samples: int,
    rng: np.random.Generator,
    correlation: np.ndarray | None = None,
) -> np.ndarray:
    """Draw zero-mean Gaussian mode amplitudes with the given per-mode variances.

    This is the *positive control generator*: the truth is exactly the analytic variance, so anything
    the estimator fails to recover here is an estimator defect, not a physics result.

    Args:
        variances: per-mode variance in m**2, strictly positive.
        n_samples: number of independent configurations to draw.
        rng: seeded NumPy generator; the seed is the measurement protocol and must be recorded.
        correlation: optional ``(n_modes, n_modes)`` correlation matrix.  ``None`` means independent
            modes, which is what the physics says.  Supplying one builds the *correlated-mode*
            negative control: the per-mode variances are unchanged, so a variance-only check passes.

    Returns:
        Array of shape ``(n_samples, n_modes)``.
    """
    var = _positive_array(variances, what="variances")
    count = _positive_int(n_samples, what="n_samples", minimum=2)
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy.random.Generator so the seed is explicit")
    n_modes = var.size
    if correlation is None:
        return rng.normal(0.0, np.sqrt(var), size=(count, n_modes))

    corr = np.asarray(correlation, dtype=float)
    if corr.shape != (n_modes, n_modes):
        raise ValueError("correlation must be square with one row per mode")
    if not np.allclose(corr, corr.T, rtol=0.0, atol=1e-12):
        raise ValueError("correlation must be symmetric")
    if not np.allclose(np.diag(corr), 1.0, rtol=0.0, atol=1e-12):
        raise ValueError("correlation must have unit diagonal")
    eigenvalues = np.linalg.eigvalsh(corr)
    if eigenvalues.min() <= 0.0:
        raise ValueError("correlation must be positive definite")
    scale = np.sqrt(var)
    covariance = corr * np.outer(scale, scale)
    factor = np.linalg.cholesky(covariance)
    white = rng.normal(0.0, 1.0, size=(count, n_modes))
    return white @ factor.T


def correlated_mode_amplitudes(
    variances: Sequence[float] | np.ndarray,
    *,
    n_samples: int,
    rng: np.random.Generator,
    rho: float,
) -> np.ndarray:
    """Negative control 3: modes with the correct variances but an equicorrelated covariance.

    Args:
        variances: per-mode variance in m**2.
        n_samples: number of configurations.
        rng: seeded generator.
        rho: off-diagonal correlation, ``|rho| < 1``.

    Returns:
        Array of shape ``(n_samples, n_modes)`` whose per-mode variances match ``variances`` in
        expectation while the modes are not independent.
    """
    var = _positive_array(variances, what="variances")
    if isinstance(rho, bool) or not isinstance(rho, (int, float)) or not -1.0 < float(rho) < 1.0:
        raise ValueError("rho must be a real number strictly inside (-1, 1)")
    n_modes = var.size
    corr = np.full((n_modes, n_modes), float(rho))
    np.fill_diagonal(corr, 1.0)
    return sample_gaussian_modes(var, n_samples=n_samples, rng=rng, correlation=corr)


def power_law_variances(
    q: Sequence[float] | np.ndarray,
    *,
    exponent: float,
    reference_q: float,
    reference_variance: float,
) -> np.ndarray:
    """Negative control 1: a single pure power law ``V(q) = V0 * (q/q0)**exponent``.

    Setting ``exponent = -2`` where the physics says ``-4`` is "tension mistaken for bending".  It is
    pinned to pass through ``(reference_q, reference_variance)`` so that it is *indistinguishable from
    the truth at one point* — the failure only becomes visible across a span of ``q``.
    """
    wavevectors = _positive_array(q, what="q")
    if isinstance(exponent, bool) or not isinstance(exponent, (int, float)):
        raise TypeError("exponent must be a real number")
    q0 = _positive_float(reference_q, what="reference_q")
    v0 = _positive_float(reference_variance, what="reference_variance")
    return v0 * (wavevectors / q0) ** float(exponent)


def mesh_truncated_variances(
    q: Sequence[float] | np.ndarray,
    variances: Sequence[float] | np.ndarray,
    *,
    q_cutoff: float,
    rolloff_exponent: float = 4.0,
) -> np.ndarray:
    """Negative control 4: the true spectrum multiplied by a discretisation roll-off.

    A finite mesh cannot represent amplitude near its Nyquist wavevector, so the measured spectrum
    bends *down* there.  Fitted, that looks like a stiffer membrane.  The roll-off used here is
    ``exp(-(q/q_cutoff)**p)``, which is smooth and leaves the low-``q`` decades untouched, so the
    detector has to find the boundary rather than being handed a discontinuity.
    """
    wavevectors = _positive_array(q, what="q")
    var = _positive_array(variances, what="variances")
    if var.size != wavevectors.size:
        raise ValueError("variances must have one entry per wavevector")
    cutoff = _positive_float(q_cutoff, what="q_cutoff")
    power = _positive_float(rolloff_exponent, what="rolloff_exponent")
    return var * np.exp(-((wavevectors / cutoff) ** power))


def measure_mode_variance(amplitudes: np.ndarray) -> np.ndarray:
    """Estimate the per-mode variance from amplitudes whose mean is known to be zero.

    ``mean(a**2)`` is used rather than a sample variance about an estimated mean: the mean of a
    fluctuation mode is zero by construction, so subtracting an estimated mean would throw away one
    degree of freedom and bias the estimate low.
    """
    samples = np.asarray(amplitudes, dtype=float)
    if samples.ndim != 2 or samples.size == 0:
        raise ValueError("amplitudes must be a non-empty (n_samples, n_modes) array")
    if not np.all(np.isfinite(samples)):
        raise ValueError("amplitudes must be finite")
    return np.mean(samples**2, axis=0)


def _log_variance_of_variance_estimator(n_samples: int) -> float:
    """Variance of ``log(S_hat)`` where ``S_hat`` averages ``N`` squared zero-mean Gaussian draws.

    ``S_hat`` is distributed as ``V * chi2_N / N``, so ``Var[log S_hat] = psi'(N/2)``, whose
    asymptotic form ``2/N`` is used here.  It is independent of ``q`` and of ``V``, which is the
    single fact that makes an unweighted log-space least squares the correct estimator.
    """
    return 2.0 / _positive_int(n_samples, what="n_samples", minimum=2)


# --------------------------------------------------------------------------------------------
# estimators
# --------------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SpectrumFit:
    """Result of fitting ``(kappa, sigma)`` to a measured mode spectrum."""

    kappa_J: float
    sigma_N_per_m: float
    kappa_ci: tuple[float, float]
    sigma_ci: tuple[float, float]
    confidence: float
    chi_square: float
    dof: int
    p_value: float
    n_modes: int
    n_samples: int
    q_min: float
    q_max: float
    decades: float
    converged: bool
    iterations: int

    @property
    def kappa_resolved(self) -> bool:
        """True when the bending term is distinguishable from zero at the declared confidence."""
        return self.kappa_ci[0] > 0.0

    @property
    def sigma_resolved(self) -> bool:
        """True when the tension term is distinguishable from zero at the declared confidence."""
        return self.sigma_ci[0] > 0.0

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of the fit."""
        return {
            "kappa_J": self.kappa_J,
            "sigma_N_per_m": self.sigma_N_per_m,
            "kappa_ci": list(self.kappa_ci),
            "sigma_ci": list(self.sigma_ci),
            "confidence": self.confidence,
            "chi_square": self.chi_square,
            "dof": self.dof,
            "p_value": self.p_value,
            "n_modes": self.n_modes,
            "n_samples": self.n_samples,
            "q_min": self.q_min,
            "q_max": self.q_max,
            "decades": self.decades,
            "converged": self.converged,
            "iterations": self.iterations,
            "kappa_resolved": self.kappa_resolved,
            "sigma_resolved": self.sigma_resolved,
            "evidence_source": EVIDENCE_SOURCE.value,
        }


def fit_helfrich_spectrum(
    q: Sequence[float] | np.ndarray,
    measured_variance: Sequence[float] | np.ndarray,
    *,
    area_m2: float,
    temperature_K: float,
    n_samples: int,
    confidence: float = 0.95,
    max_iterations: int = 200,
    tolerance: float = 1e-12,
) -> SpectrumFit:
    """Recover ``kappa`` and ``sigma`` from a measured per-mode variance spectrum.

    The fit is a Gauss-Newton minimisation of the *log-space* residual

    ``r_i = log(kB*T / (A * S_i)) - log(kappa*q_i**4 + sigma*q_i**2)``

    with the a-priori residual standard deviation ``sqrt(2/n_samples)``.  Log space is not cosmetic:
    the sampling error of a variance estimator is multiplicative and ``q``-independent there, which is
    what makes an unweighted least squares the correct estimator and what makes the returned
    chi-square an honest goodness-of-fit rather than a rescaled residual norm.

    Args:
        q: wavevectors in 1/m.
        measured_variance: measured ``<|h_q|^2>`` in m**2, one per wavevector.
        area_m2: projected patch area in m**2.
        temperature_K: temperature the amplitudes were sampled at.
        n_samples: number of configurations behind each variance estimate.
        confidence: two-sided confidence level for the returned intervals.
        max_iterations: Gauss-Newton iteration cap.
        tolerance: relative parameter change at which iteration stops.

    Returns:
        A :class:`SpectrumFit`.  ``converged`` is reported, never silently assumed.

    Raises:
        ValueError: on non-positive inputs, mismatched lengths, fewer than three modes (two parameters
            need at least one degree of freedom for the goodness-of-fit to mean anything), or a
            singular design.
    """
    wavevectors = _positive_array(q, what="q")
    variance = _positive_array(measured_variance, what="measured_variance")
    if variance.size != wavevectors.size:
        raise ValueError("measured_variance must have one entry per wavevector")
    if wavevectors.size < 3:
        raise ValueError("at least three modes are needed to fit two parameters and test the fit")
    area = _positive_float(area_m2, what="area_m2")
    kt = thermal_energy(temperature_K)
    count = _positive_int(n_samples, what="n_samples", minimum=2)
    level = _confidence(confidence)
    log_sigma = math.sqrt(_log_variance_of_variance_estimator(count))

    # Target: y_i = kB*T/(A*S_i) = kappa*q^4 + sigma*q^2 in the noiseless case.
    y = kt / (area * variance)
    log_y = np.log(y)

    # Non-dimensionalise before solving.  In SI a membrane band spans q^4 ~ 1e36 against q^2 ~ 1e18,
    # so the design matrix is numerically rank deficient in raw units and the solve reports "sigma is
    # not constrained" for data that constrain it perfectly well.  Working in u = q/q_ref and
    # y/y_ref makes both columns O(1); the parameters are transformed back exactly afterwards.
    q_ref = math.sqrt(float(wavevectors.min()) * float(wavevectors.max()))
    y_ref = float(np.exp(np.mean(log_y)))
    u = wavevectors / q_ref
    basis = np.column_stack((u**4, u**2))
    log_y_scaled = log_y - math.log(y_ref)
    y_scaled = y / y_ref

    # Seed with a weighted linear solve in y-space (weights 1/y approximate log space).
    weights = 1.0 / y_scaled
    seed, *_ = np.linalg.lstsq(basis * weights[:, None], y_scaled * weights, rcond=None)
    params = np.clip(np.asarray(seed, dtype=float), 0.0, None)
    if not np.any(params > 0.0):
        params = np.array([float(np.mean(y_scaled) / np.mean(u**4)), 0.0])

    converged = False
    iterations = 0
    step_tolerance = _positive_float(tolerance, what="tolerance")
    for attempt in range(1, _positive_int(max_iterations, what="max_iterations") + 1):
        iterations = attempt
        density = basis @ params
        if np.any(density <= 0.0):
            raise ValueError(
                "fit produced a non-positive restoring density; inputs are inconsistent"
            )
        residual = log_y_scaled - np.log(density)
        jacobian = basis / density[:, None]
        step, *_ = np.linalg.lstsq(jacobian, residual, rcond=None)
        trial = np.clip(params + step, 0.0, None)
        scale = max(float(np.max(np.abs(params))), 1e-300)
        moved = float(np.max(np.abs(trial - params)) / scale)
        params = trial
        if moved < step_tolerance:
            converged = True
            break

    density = basis @ params
    if np.any(density <= 0.0):
        raise ValueError("fit produced a non-positive restoring density; inputs are inconsistent")
    residual = log_y_scaled - np.log(density)
    chi_square = float(np.sum((residual / log_sigma) ** 2))
    dof = int(wavevectors.size - 2)
    p_value = chi_square_survival(chi_square, dof) if dof >= 1 else float("nan")

    jacobian = basis / density[:, None]
    gram = jacobian.T @ jacobian
    if np.linalg.matrix_rank(gram) < 2:
        raise ValueError("the wavevector range does not constrain both kappa and sigma")
    covariance = np.linalg.inv(gram) * log_sigma**2
    z = NormalDist().inv_cdf(0.5 + 0.5 * level)
    errors = np.sqrt(np.clip(np.diag(covariance), 0.0, None))

    # Undo the non-dimensionalisation exactly: kappa = a * y_ref / q_ref**4, sigma = b * y_ref/q_ref**2.
    kappa_scale = y_ref / q_ref**4
    sigma_scale = y_ref / q_ref**2
    kappa = float(params[0]) * kappa_scale
    sigma = float(params[1]) * sigma_scale
    kappa_err = float(errors[0]) * kappa_scale
    sigma_err = float(errors[1]) * sigma_scale
    kappa_ci = (max(0.0, kappa - z * kappa_err), kappa + z * kappa_err)
    sigma_ci = (max(0.0, sigma - z * sigma_err), sigma + z * sigma_err)

    q_low = float(wavevectors.min())
    q_high = float(wavevectors.max())
    return SpectrumFit(
        kappa_J=kappa,
        sigma_N_per_m=sigma,
        kappa_ci=kappa_ci,
        sigma_ci=sigma_ci,
        confidence=level,
        chi_square=chi_square,
        dof=dof,
        p_value=float(p_value),
        n_modes=int(wavevectors.size),
        n_samples=count,
        q_min=q_low,
        q_max=q_high,
        decades=math.log10(q_high / q_low),
        converged=converged,
        iterations=iterations,
    )


@dataclass(frozen=True, slots=True)
class ExponentFit:
    """Result of a log-log power-law exponent measurement over a stated ``q`` span."""

    exponent: float
    standard_error: float
    ci: tuple[float, float]
    confidence: float
    decades: float
    n_modes: int
    n_samples: int
    q_min: float
    q_max: float

    def excludes(self, candidate: float) -> bool:
        """True when the confidence interval does not contain ``candidate``.

        This is the test that a two-parameter fit cannot fake: it asks about the *slope*, not about
        whether some ``(kappa, sigma)`` reproduces the data.
        """
        return not (self.ci[0] <= float(candidate) <= self.ci[1])

    def separation_sigma(self, candidate: float) -> float:
        """Distance from ``candidate`` in units of the exponent's own standard error."""
        if self.standard_error == 0.0:
            return math.inf
        return abs(self.exponent - float(candidate)) / self.standard_error


def fit_power_law_exponent(
    q: Sequence[float] | np.ndarray,
    measured_variance: Sequence[float] | np.ndarray,
    *,
    n_samples: int,
    confidence: float = 0.95,
) -> ExponentFit:
    """Estimate the log-log slope of a measured spectrum with an a-priori error bar.

    The standard error uses the known ``2/n_samples`` log-variance of the variance estimator rather
    than the residual scatter.  That matters: a *curved* spectrum fitted with a straight line has
    inflated residuals, and using them would inflate the error bar until the wrong exponent became
    "consistent".  Here the error bar is fixed by the sampling protocol, so curvature shows up as a
    rejected exponent, which is the honest outcome.

    Args:
        q: wavevectors in 1/m, at least three, not all equal.
        measured_variance: measured ``<|h_q|^2>`` in m**2.
        n_samples: configurations behind each variance estimate.
        confidence: two-sided confidence level.

    Returns:
        An :class:`ExponentFit`.
    """
    wavevectors = _positive_array(q, what="q")
    variance = _positive_array(measured_variance, what="measured_variance")
    if variance.size != wavevectors.size:
        raise ValueError("measured_variance must have one entry per wavevector")
    if wavevectors.size < 3:
        raise ValueError("at least three modes are needed to estimate a slope and its error")
    count = _positive_int(n_samples, what="n_samples", minimum=2)
    level = _confidence(confidence)

    x = np.log(wavevectors)
    y = np.log(variance)
    centred = x - x.mean()
    sxx = float(np.sum(centred**2))
    if sxx <= 0.0:
        raise ValueError("wavevectors span no range; an exponent is not identifiable")
    slope = float(np.sum(centred * y) / sxx)
    log_sigma = math.sqrt(_log_variance_of_variance_estimator(count))
    stderr = log_sigma / math.sqrt(sxx)
    z = NormalDist().inv_cdf(0.5 + 0.5 * level)
    return ExponentFit(
        exponent=slope,
        standard_error=stderr,
        ci=(slope - z * stderr, slope + z * stderr),
        confidence=level,
        decades=math.log10(float(wavevectors.max()) / float(wavevectors.min())),
        n_modes=int(wavevectors.size),
        n_samples=count,
        q_min=float(wavevectors.min()),
        q_max=float(wavevectors.max()),
    )


def decades_required_for_exponent_separation(
    *,
    delta_exponent: float = 2.0,
    n_samples: int,
    modes_per_decade: float,
    n_sigma: float = 5.0,
) -> float:
    """Predict the ``q`` span needed to separate two exponents, in decades.

    For ``M`` log-uniformly spaced modes over ``D`` decades the regressor spread is
    ``Sxx = M*(M+1)/(12*(M-1)) * (D*ln10)**2``, so the slope standard error is
    ``sqrt(2/N)/sqrt(Sxx)``.  Requiring ``n_sigma * stderr <= delta_exponent`` and solving for ``D``
    (``M`` itself depends on ``D``) gives the answer; the root is found by bisection because of that
    coupling.

    This is a *design* number: it tells a future native run how large a membrane patch it must
    simulate before a mode-by-mode spectrum check can distinguish ``q**-2`` from ``q**-4`` at all.

    Args:
        delta_exponent: the exponent difference to resolve.  ``2`` for tension-vs-bending.
        n_samples: independent configurations per mode.
        modes_per_decade: mode density of the planned measurement.
        n_sigma: required separation in standard errors.

    Returns:
        The required span in decades of ``q``.

    Raises:
        ValueError: if no span below 100 decades suffices, which means the sample count is the binding
            constraint rather than the patch size.
    """
    delta = _positive_float(delta_exponent, what="delta_exponent")
    count = _positive_int(n_samples, what="n_samples", minimum=2)
    density = _positive_float(modes_per_decade, what="modes_per_decade")
    sigmas = _positive_float(n_sigma, what="n_sigma")
    log_sigma = math.sqrt(_log_variance_of_variance_estimator(count))

    def shortfall(decades: float) -> float:
        modes = max(3.0, math.floor(density * decades) + 1.0)
        sxx = modes * (modes + 1.0) / (12.0 * (modes - 1.0)) * (decades * math.log(10.0)) ** 2
        return sigmas * log_sigma / math.sqrt(sxx) - delta

    low, high = 1e-6, 100.0
    if shortfall(low) < 0.0:
        return low
    if shortfall(high) > 0.0:
        raise ValueError(
            "no span below 100 decades separates these exponents; increase n_samples or "
            "modes_per_decade instead of the patch size"
        )
    for _ in range(200):
        mid = 0.5 * (low + high)
        if shortfall(mid) > 0.0:
            low = mid
        else:
            high = mid
    return high


@dataclass(frozen=True, slots=True)
class ExponentSeparationLadder:
    """Measured rejection rate of a wrong exponent as a function of the ``q`` span."""

    decades: tuple[float, ...]
    rejection_rate: tuple[float, ...]
    mean_separation_sigma: tuple[float, ...]
    n_trials: int
    n_samples: int
    modes_per_decade: float
    truth_exponent: float
    expected_exponent: float
    confidence: float
    required_power: float
    seed: int
    predicted_decades: float

    @property
    def measured_decades(self) -> float | None:
        """Smallest tested span whose rejection rate reaches ``required_power``, or ``None``."""
        for span, rate in zip(self.decades, self.rejection_rate, strict=True):
            if rate >= self.required_power:
                return span
        return None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of the ladder."""
        return {
            "decades": list(self.decades),
            "rejection_rate": list(self.rejection_rate),
            "mean_separation_sigma": list(self.mean_separation_sigma),
            "n_trials": self.n_trials,
            "n_samples": self.n_samples,
            "modes_per_decade": self.modes_per_decade,
            "truth_exponent": self.truth_exponent,
            "expected_exponent": self.expected_exponent,
            "confidence": self.confidence,
            "required_power": self.required_power,
            "seed": self.seed,
            "predicted_decades": self.predicted_decades,
            "measured_decades": self.measured_decades,
            "evidence_source": EVIDENCE_SOURCE.value,
        }


def measure_exponent_separation_ladder(
    *,
    decades: Sequence[float],
    modes_per_decade: float,
    n_samples: int,
    n_trials: int,
    seed: int,
    truth_exponent: float = -2.0,
    expected_exponent: float = -4.0,
    confidence: float = 0.95,
    required_power: float = 1.0,
    q_min: float = 1.0e5,
) -> ExponentSeparationLadder:
    """Measure, by sampling, how many decades of ``q`` are needed to reject a wrong exponent.

    For each span, ``n_trials`` synthetic spectra are drawn from the *wrong* power law and the
    exponent is fitted; a trial counts as a rejection when the confidence interval excludes the
    expected exponent.  Nothing here is a prediction: the returned rates are counted outcomes for the
    stated seed, and :meth:`ExponentSeparationLadder.measured_decades` reports the smallest tested
    span that reached the required power.

    Args:
        decades: ascending spans to test, in decades of ``q``.
        modes_per_decade: mode density.
        n_samples: configurations behind each per-mode variance.
        n_trials: independent repetitions per span.
        seed: RNG seed; part of the measurement record.
        truth_exponent: exponent the data are actually drawn from.
        expected_exponent: exponent the checker expects, and which must be rejected.
        confidence: confidence level of the exponent interval.
        required_power: rejection rate that counts as "separated"; ``1.0`` means every trial.
        q_min: lower end of the tested band in 1/m.  Only the span matters for the slope error, so
            this is a placement choice, not a result.

    Returns:
        An :class:`ExponentSeparationLadder`.
    """
    spans = tuple(float(_positive_float(value, what="decades entry")) for value in decades)
    if not spans:
        raise ValueError("decades must not be empty")
    if any(later <= earlier for earlier, later in zip(spans, spans[1:], strict=False)):
        raise ValueError("decades must be strictly ascending")
    density = _positive_float(modes_per_decade, what="modes_per_decade")
    trials = _positive_int(n_trials, what="n_trials")
    level = _confidence(confidence)
    power = _confidence(required_power) if required_power < 1.0 else 1.0
    low_q = _positive_float(q_min, what="q_min")
    rng = np.random.default_rng(_positive_int(seed, what="seed", minimum=0))

    rates: list[float] = []
    separations: list[float] = []
    for span in spans:
        n_modes = max(3, int(math.floor(density * span)) + 1)
        q = log_spaced_wavevectors(low_q, low_q * 10.0**span, n_modes)
        truth = power_law_variances(
            q, exponent=truth_exponent, reference_q=low_q, reference_variance=1.0e-18
        )
        rejected = 0
        sigma_sum = 0.0
        for _ in range(trials):
            amplitudes = sample_gaussian_modes(truth, n_samples=n_samples, rng=rng)
            measured = measure_mode_variance(amplitudes)
            fit = fit_power_law_exponent(q, measured, n_samples=n_samples, confidence=level)
            if fit.excludes(expected_exponent):
                rejected += 1
            sigma_sum += fit.separation_sigma(expected_exponent)
        rates.append(rejected / trials)
        separations.append(sigma_sum / trials)

    z = NormalDist().inv_cdf(0.5 + 0.5 * level)
    predicted = decades_required_for_exponent_separation(
        delta_exponent=abs(float(expected_exponent) - float(truth_exponent)),
        n_samples=n_samples,
        modes_per_decade=density,
        n_sigma=z,
    )
    return ExponentSeparationLadder(
        decades=spans,
        rejection_rate=tuple(rates),
        mean_separation_sigma=tuple(separations),
        n_trials=trials,
        n_samples=_positive_int(n_samples, what="n_samples", minimum=2),
        modes_per_decade=density,
        truth_exponent=float(truth_exponent),
        expected_exponent=float(expected_exponent),
        confidence=level,
        required_power=power,
        seed=int(seed),
        predicted_decades=predicted,
    )


class PatchBoundaryCondition(StrEnum):
    """Which boundary condition a native membrane patch is declared to carry.

    This is not a formatting choice.  The relation between the patch side ``L`` and the smallest
    representable wavevector ``q_min`` — the *fundamental* — is a property of the boundary condition,
    and so is whether the grid's far edge is a distinct node or the same node as the near edge.  Mixing
    the two conventions produces a node count that is right under neither, which is exactly what an
    external review found in the version this replaces (see :data:`PATCH_BOUNDARY_CONDITION_IS_A_DECLARATION`).
    """

    #: Periodic (torus) patch.  Fundamental ``q_min = 2*pi/L``; the grid has ``L/a`` unique nodes per
    #: side with NO duplicated endpoint, because node ``L`` *is* node ``0``.
    PERIODIC = "periodic"
    #: Fixed (Dirichlet) edges, the 2-D analogue of this module's :class:`TautStringSpectrum`.  The
    #: mode functions are sines, so the fundamental is ``q_min = pi/L``, and the grid is
    #: endpoint-inclusive with ``L/a + 1`` nodes per side of which the two boundary rows are pinned.
    FIXED_EDGE = "fixed-edge"


#: Why :func:`native_patch_requirement` refuses to guess a boundary condition.
#:
#: Correction from the external review of 2026-07-29 (finding 10).  The superseded implementation
#: derived the patch side from the PERIODIC fundamental ``q_min = 2*pi/L`` — giving ``L/a =
#: 2*q_max/q_min = 200`` for the engine constants — and then added an endpoint for 201 nodes per side,
#: i.e. ``201**2 = 40,401``.  A periodic grid has no duplicated endpoint, so that number is right under
#: neither convention: periodic gives ``200**2 = 40,000`` and fixed-edge gives ``101**2 = 10,201``
#: (the fixed-edge patch is also only half as wide, because its fundamental is ``pi/L``, not
#: ``2*pi/L``).  The counterexample is arithmetic, not statistical, and the number was a candidate
#: design input for a native run.
PATCH_BOUNDARY_CONDITION_IS_A_DECLARATION = (
    "the fundamental relation between patch side and smallest wavevector is a property of the "
    "boundary condition (periodic: q_min = 2*pi/L; fixed edges: q_min = pi/L), and so is whether the "
    "far edge is a distinct node.  A patch design that does not state its boundary condition is not a "
    "design.  This module will not choose one for the caller: which one applies is a statement about "
    "what the engine's membrane shell implements, and this lane has not read that code."
)

#: The reason to PREFER periodic for an S1 patch — stated as a recommendation with its argument, not
#: as a default.  Ratifying it is a PI decision, and it is only correct if the engine agrees.
PATCH_PERIODIC_RATIONALE = (
    "S1's estimator assumes the modes are INDEPENDENT with variance kB*T/(A*(kappa*q^4 + sigma*q^2)), "
    "which is the periodic Fourier result for a translation-invariant patch; a fixed-edge patch has "
    "sine modes, pinned boundary rows and edge-localised corrections that the variance law above does "
    "not carry.  A representative patch cut from a much larger membrane is therefore normally "
    "periodic.  That is an argument, not a measurement: this lane did not read the engine's membrane "
    "shell, so the caller must declare, and ratifying 'periodic' for the native run is a PI decision."
)


@dataclass(frozen=True, slots=True)
class PatchDesign:
    """What a native membrane patch must be, for an S1 spectrum check to be able to conclude.

    Every geometric field is conditioned on :attr:`boundary_condition`; none of them is a
    boundary-condition-free fact about the ``q`` band.
    """

    boundary_condition: PatchBoundaryCondition
    fundamental_relation: str
    q_min_per_m: float
    q_max_per_m: float
    decades: float
    crossover_q_per_m: float
    patch_side_m: float
    mesh_spacing_m: float
    nodes_per_side: int
    total_nodes: int
    free_nodes_per_side: int
    total_free_nodes: int

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able design record."""
        return {
            "boundary_condition": self.boundary_condition.value,
            "fundamental_relation": self.fundamental_relation,
            "q_min_per_m": self.q_min_per_m,
            "q_max_per_m": self.q_max_per_m,
            "decades": self.decades,
            "crossover_q_per_m": self.crossover_q_per_m,
            "patch_side_m": self.patch_side_m,
            "mesh_spacing_m": self.mesh_spacing_m,
            "nodes_per_side": self.nodes_per_side,
            "total_nodes": self.total_nodes,
            "free_nodes_per_side": self.free_nodes_per_side,
            "total_free_nodes": self.total_free_nodes,
            "evidence_source": EVIDENCE_SOURCE.value,
        }


def native_patch_requirement(
    spectrum: HelfrichSpectrum,
    *,
    boundary_condition: PatchBoundaryCondition | str,
    decades_below_crossover: float = 1.0,
    decades_above_crossover: float = 1.0,
) -> PatchDesign:
    """Turn a required ``q`` span into a patch size and a mesh spacing a native run must meet.

    A decade requirement is abstract until it is converted into geometry, and the conversion is what
    makes S1 a constraint on a future native run rather than a comment about it.  The largest
    wavevector a mesh of spacing ``a`` can carry is ``pi/a`` under either boundary condition, so the
    mesh spacing is common; the patch SIDE and the NODE COUNT are not, because both follow from the
    fundamental, which is a property of the boundary condition:

    ==============  ===================  ==============================  ==================
    boundary        fundamental          side for a given ``q_min``      nodes per side
    ==============  ===================  ==============================  ==================
    periodic        ``q_min = 2*pi/L``   ``L = 2*pi/q_min``              ``L/a`` (no endpoint)
    fixed edges     ``q_min = pi/L``     ``L = pi/q_min``                ``L/a + 1``
    ==============  ===================  ==============================  ==================

    ``boundary_condition`` is therefore REQUIRED and has no default: see
    :data:`PATCH_BOUNDARY_CONDITION_IS_A_DECLARATION` for the counterexample that made it required, and
    :data:`PATCH_PERIODIC_RATIONALE` for the argument — not the decision — in favour of periodic.

    Anchoring the band on the crossover ``q* = sqrt(sigma/kappa)`` is the demanding requirement,
    because resolving BOTH moduli needs both regimes present — an exponent measurement alone needs far
    less span.

    Args:
        spectrum: the spectrum whose crossover anchors the band.
        boundary_condition: the declared boundary condition, as a :class:`PatchBoundaryCondition` or
            its string value.  There is deliberately no default.
        decades_below_crossover: how far below ``q*`` the band must reach.
        decades_above_crossover: how far above ``q*`` the band must reach.

    Returns:
        A :class:`PatchDesign`.

    Raises:
        ValueError: if the spectrum has no finite crossover (pure bending or pure tension), in which
            case there is no two-regime band to design for; or if ``boundary_condition`` is not one of
            the declared conventions.
    """
    if not isinstance(spectrum, HelfrichSpectrum):
        raise TypeError("spectrum must be a HelfrichSpectrum")
    try:
        condition = PatchBoundaryCondition(boundary_condition)
    except ValueError as exc:
        allowed = ", ".join(member.value for member in PatchBoundaryCondition)
        raise ValueError(
            f"boundary_condition must be declared as one of ({allowed}); "
            f"{PATCH_BOUNDARY_CONDITION_IS_A_DECLARATION}"
        ) from exc
    q_star = spectrum.crossover_wavevector_per_m
    if not math.isfinite(q_star) or q_star <= 0.0:
        raise ValueError(
            "the spectrum has no finite crossover, so there is no two-regime band to design for"
        )
    below = _positive_float(decades_below_crossover, what="decades_below_crossover")
    above = _positive_float(decades_above_crossover, what="decades_above_crossover")
    q_min = q_star / 10.0**below
    q_max = q_star * 10.0**above
    spacing = math.pi / q_max
    if condition is PatchBoundaryCondition.PERIODIC:
        fundamental = "q_min = 2*pi/L (periodic fundamental)"
        side = 2.0 * math.pi / q_min
        # A periodic grid does NOT duplicate its endpoint: node L is node 0.  ceil, not ceil+1.
        nodes_per_side = int(math.ceil(side / spacing))
        # Every node of a torus is free; there are no pinned boundary rows.
        free_per_side = nodes_per_side
    else:
        fundamental = "q_min = pi/L (fixed-edge / Dirichlet fundamental, as for the taut string)"
        side = math.pi / q_min
        # Endpoint-inclusive: both edges are nodes, and both are pinned.
        nodes_per_side = int(math.ceil(side / spacing)) + 1
        free_per_side = nodes_per_side - 2
    return PatchDesign(
        boundary_condition=condition,
        fundamental_relation=fundamental,
        q_min_per_m=q_min,
        q_max_per_m=q_max,
        decades=below + above,
        crossover_q_per_m=q_star,
        patch_side_m=side,
        mesh_spacing_m=spacing,
        nodes_per_side=nodes_per_side,
        total_nodes=nodes_per_side**2,
        free_nodes_per_side=free_per_side,
        total_free_nodes=free_per_side**2,
    )


@dataclass(frozen=True, slots=True)
class EquipartitionReport:
    """Whether the measured amplitudes carry the thermal energy they are supposed to."""

    recovered_thermal_energy_J: float
    reference_thermal_energy_J: float
    ratio: float
    ratio_ci: tuple[float, float]
    confidence: float
    n_modes: int
    n_samples: int
    consistent: bool
    degeneracy_note: str = EQUIPARTITION_SCALE_DEGENERACY

    @property
    def implied_temperature_K(self) -> float:
        """Temperature the amplitudes would correspond to, in kelvin."""
        return self.recovered_thermal_energy_J / BOLTZMANN_J_PER_K


def check_equipartition(
    q: Sequence[float] | np.ndarray,
    measured_variance: Sequence[float] | np.ndarray,
    *,
    spectrum: HelfrichSpectrum,
    n_samples: int,
    confidence: float = 0.99,
) -> EquipartitionReport:
    """Negative control 2: catch amplitudes scaled so the implied ``kB*T`` is wrong.

    ``kappa`` and ``sigma`` must come from ``spectrum`` and must have been determined
    *independently of these amplitudes*.  That is not a stylistic preference — see
    :data:`EQUIPARTITION_SCALE_DEGENERACY`.  Given them, each mode yields an independent estimate
    ``kB*T_hat = S_i * A * (kappa*q_i**4 + sigma*q_i**2)``; these are combined in log space, where the
    sampling error is ``sqrt(2/(n_samples*n_modes))``.

    Args:
        q: wavevectors in 1/m.
        measured_variance: measured ``<|h_q|^2>`` in m**2.
        spectrum: independently determined spectrum supplying ``kappa``, ``sigma``, ``A`` and the
            reference temperature.
        n_samples: configurations behind each variance.
        confidence: two-sided confidence level.

    Returns:
        An :class:`EquipartitionReport`.  ``consistent`` is true when the interval on the ratio
        contains ``1``.
    """
    if not isinstance(spectrum, HelfrichSpectrum):
        raise TypeError("spectrum must be a HelfrichSpectrum determined independently of the data")
    wavevectors = _positive_array(q, what="q")
    variance = _positive_array(measured_variance, what="measured_variance")
    if variance.size != wavevectors.size:
        raise ValueError("measured_variance must have one entry per wavevector")
    count = _positive_int(n_samples, what="n_samples", minimum=2)
    level = _confidence(confidence)

    per_mode = variance * spectrum.area_m2 * spectrum.restoring_density(wavevectors)
    log_mean = float(np.mean(np.log(per_mode)))
    recovered = math.exp(log_mean)
    reference = spectrum.thermal_energy_J
    stderr = math.sqrt(_log_variance_of_variance_estimator(count) / wavevectors.size)
    z = NormalDist().inv_cdf(0.5 + 0.5 * level)
    ratio = recovered / reference
    ci = (ratio * math.exp(-z * stderr), ratio * math.exp(z * stderr))
    return EquipartitionReport(
        recovered_thermal_energy_J=recovered,
        reference_thermal_energy_J=reference,
        ratio=ratio,
        ratio_ci=ci,
        confidence=level,
        n_modes=int(wavevectors.size),
        n_samples=count,
        consistent=ci[0] <= 1.0 <= ci[1],
    )


@dataclass(frozen=True, slots=True)
class ModeIndependenceReport:
    """Whether modes that must be independent actually are."""

    n_modes: int
    n_samples: int
    max_abs_correlation: float
    max_pair: tuple[int, int]
    min_p_value: float
    family_size: int
    bonferroni_alpha: float
    independent: bool


def check_mode_independence(
    amplitudes: np.ndarray, *, confidence: float = 0.99
) -> ModeIndependenceReport:
    """Negative control 3: a per-mode variance check passes correlated modes, so test separately.

    Under the Helfrich quadratic form distinct modes are independent, which is a statement about the
    *joint* law and not about any marginal.  Correlating two modes while leaving both variances
    correct is therefore invisible to every check in :func:`fit_helfrich_spectrum`.

    The statistic is the Fisher ``z`` transform of each pairwise sample correlation,
    ``z = atanh(r)*sqrt(N-3)``, with a Bonferroni correction over all ``M*(M-1)/2`` pairs.

    Args:
        amplitudes: ``(n_samples, n_modes)`` array.
        confidence: family-wise confidence level.

    Returns:
        A :class:`ModeIndependenceReport`.
    """
    samples = np.asarray(amplitudes, dtype=float)
    if samples.ndim != 2 or samples.shape[1] < 2:
        raise ValueError("amplitudes must be (n_samples, n_modes) with at least two modes")
    if samples.shape[0] < 5:
        raise ValueError("at least five samples are needed for a Fisher z test")
    if not np.all(np.isfinite(samples)):
        raise ValueError("amplitudes must be finite")
    level = _confidence(confidence)

    n_obs, n_modes = samples.shape
    corr = np.corrcoef(samples, rowvar=False)
    if not np.all(np.isfinite(corr)):
        raise ValueError("a mode column is constant; a correlation is undefined")
    upper = np.triu_indices(n_modes, k=1)
    values = corr[upper]
    index = int(np.argmax(np.abs(values)))
    max_r = float(np.abs(values[index]))
    max_pair = (int(upper[0][index]), int(upper[1][index]))

    clipped = min(max_r, 1.0 - 1e-15)
    z_stat = math.atanh(clipped) * math.sqrt(n_obs - 3)
    min_p = 2.0 * (1.0 - NormalDist().cdf(abs(z_stat)))
    family = n_modes * (n_modes - 1) // 2
    alpha = (1.0 - level) / family
    return ModeIndependenceReport(
        n_modes=n_modes,
        n_samples=n_obs,
        max_abs_correlation=max_r,
        max_pair=max_pair,
        min_p_value=min_p,
        family_size=family,
        bonferroni_alpha=alpha,
        independent=min_p >= alpha,
    )


@dataclass(frozen=True, slots=True)
class TrustedRangeReport:
    """Which part of the measured ``q`` range may be fitted, and why the rest may not."""

    q_trust_max: float
    n_modes_trusted: int
    n_modes_rejected: int
    cutoff_detected: bool
    chi_square: float
    dof: int
    p_value: float
    confidence: float
    reason: str
    nyquist_q: float | None = None


def trusted_wavevector_range(
    q: Sequence[float] | np.ndarray,
    measured_variance: Sequence[float] | np.ndarray,
    *,
    area_m2: float,
    temperature_K: float,
    n_samples: int,
    confidence: float = 0.99,
    min_modes: int = 5,
    nyquist_q: float | None = None,
) -> TrustedRangeReport:
    """Negative control 4: separate "the physics has no power there" from "my grid cannot represent it".

    A mesh-truncated spectrum is still a smooth decreasing curve, and fitting it returns a *larger*
    ``kappa`` with a perfectly respectable-looking parameter value.  The detector is deliberately
    one-sided and low-``q``-anchored: the trusted range is the **largest low-``q`` prefix whose
    Helfrich fit passes a chi-square goodness-of-fit at the declared level**.  The level is an
    argument, fixed before the data; the procedure never re-thresholds after seeing a result.

    Args:
        q: wavevectors in 1/m, ascending.
        measured_variance: measured ``<|h_q|^2>`` in m**2.
        area_m2: patch area in m**2.
        temperature_K: temperature.
        n_samples: configurations behind each variance.
        confidence: goodness-of-fit level; the prefix is accepted when ``p > 1 - confidence``.
        min_modes: smallest prefix that may be reported as trusted.
        nyquist_q: optional independently known grid Nyquist wavevector, recorded for comparison.
            It is *not* used to make the decision — the point is to detect the cutoff from the data.

    Returns:
        A :class:`TrustedRangeReport`.
    """
    wavevectors = _positive_array(q, what="q")
    variance = _positive_array(measured_variance, what="measured_variance")
    if variance.size != wavevectors.size:
        raise ValueError("measured_variance must have one entry per wavevector")
    if np.any(np.diff(wavevectors) <= 0.0):
        raise ValueError("q must be strictly ascending")
    floor = _positive_int(min_modes, what="min_modes", minimum=3)
    if wavevectors.size < floor:
        raise ValueError("fewer modes than min_modes")
    level = _confidence(confidence)
    alpha = 1.0 - level

    best: tuple[int, SpectrumFit] | None = None
    for cut in range(wavevectors.size, floor - 1, -1):
        fit = fit_helfrich_spectrum(
            wavevectors[:cut],
            variance[:cut],
            area_m2=area_m2,
            temperature_K=temperature_K,
            n_samples=n_samples,
            confidence=level,
        )
        if fit.p_value > alpha:
            best = (cut, fit)
            break

    if best is None:
        fit = fit_helfrich_spectrum(
            wavevectors[:floor],
            variance[:floor],
            area_m2=area_m2,
            temperature_K=temperature_K,
            n_samples=n_samples,
            confidence=level,
        )
        return TrustedRangeReport(
            q_trust_max=float(wavevectors[floor - 1]),
            n_modes_trusted=0,
            n_modes_rejected=int(wavevectors.size),
            cutoff_detected=True,
            chi_square=fit.chi_square,
            dof=fit.dof,
            p_value=fit.p_value,
            confidence=level,
            reason=(
                "no low-q prefix of at least min_modes passes the declared goodness-of-fit; the "
                "spectrum is not Helfrich anywhere in the measured band and must not be fitted"
            ),
            nyquist_q=nyquist_q,
        )

    cut, fit = best
    rejected = int(wavevectors.size - cut)
    if rejected == 0:
        reason = "the whole measured band passes the declared goodness-of-fit"
    else:
        reason = (
            f"the top {rejected} mode(s) are excluded: including them fails the declared "
            "chi-square goodness-of-fit, which is the signature of a grid cutoff rather than of "
            "physics with no power there"
        )
    return TrustedRangeReport(
        q_trust_max=float(wavevectors[cut - 1]),
        n_modes_trusted=int(cut),
        n_modes_rejected=rejected,
        cutoff_detected=rejected > 0,
        chi_square=fit.chi_square,
        dof=fit.dof,
        p_value=fit.p_value,
        confidence=level,
        reason=reason,
        nyquist_q=nyquist_q,
    )


#: What may and may not be said about the size of the mesh-cutoff damage.
#:
#: Correction from the external review of 2026-07-29 (finding 9).  The S1 lane had reported "fitting
#: kappa through a mesh cutoff inflates it by 33.63x" as if that multiple characterised the failure.
#: It does not: it is the value at ONE cutoff position.  Holding the whole design fixed (band
#: ``1e5..3e8`` 1/m, 40 modes, ``N = 4000``, seed 20260729) and moving only ``q_cutoff``, the
#: reviewer measured — and this lane reproduced to five significant figures — ``8e7 -> 1.9346e5``,
#: ``1.0e8 -> 476.53``, ``1.2e8 -> 33.628``, ``1.5e8 -> 5.358``, ``2e8 -> 1.8145``,
#: ``3e8 -> 1.1304``.  The damage is a monotone function of where the cutoff sits inside the fitted
#: band; it is unbounded as the cutoff descends into that band and tends to 1 as it leaves.
#: :func:`measure_cutoff_damage_ladder` measures that curve, which is the honest deliverable.
CUTOFF_DAMAGE_IS_NOT_A_MULTIPLE = (
    "the inflation of kappa from fitting through a mesh cutoff is a FUNCTION of where the cutoff sits "
    "inside the fitted band, not a transferable multiple.  It is unbounded as the cutoff descends "
    "into the band and tends to 1 as it rises out of it, so a single number such as '33.63x' may only "
    "be quoted together with the band, the mode count, the sample count and the cutoff position that "
    "produced it.  The transferable claim is the SHAPE — monotone and unbounded — plus the measured "
    "fact that the trusted-range guard recovers kappa at every cutoff position on the ladder."
)


@dataclass(frozen=True, slots=True)
class CutoffDamageLadder:
    """Measured kappa inflation as a function of mesh-cutoff position, with and without the guard.

    The point of this object is that :attr:`naive_kappa_ratio` is a curve.  Reporting one entry of it
    as "the" damage multiple is the overgeneralisation :data:`CUTOFF_DAMAGE_IS_NOT_A_MULTIPLE`
    forbids, and the reason this ladder exists at all.
    """

    q_cutoffs_per_m: tuple[float, ...]
    cutoff_over_fitted_q_max: tuple[float, ...]
    decades_of_band_above_cutoff: tuple[float, ...]
    naive_kappa_ratio: tuple[float, ...]
    guarded_kappa_ratio: tuple[float | None, ...]
    n_modes_trusted: tuple[int, ...]
    fitted_q_min_per_m: float
    fitted_q_max_per_m: float
    n_modes: int
    n_samples: int
    confidence: float
    seed: int
    scope_note: str = CUTOFF_DAMAGE_IS_NOT_A_MULTIPLE

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able record of the ladder."""
        return {
            "q_cutoffs_per_m": list(self.q_cutoffs_per_m),
            "cutoff_over_fitted_q_max": list(self.cutoff_over_fitted_q_max),
            "decades_of_band_above_cutoff": list(self.decades_of_band_above_cutoff),
            "naive_kappa_ratio": list(self.naive_kappa_ratio),
            "guarded_kappa_ratio": list(self.guarded_kappa_ratio),
            "n_modes_trusted": list(self.n_modes_trusted),
            "fitted_q_min_per_m": self.fitted_q_min_per_m,
            "fitted_q_max_per_m": self.fitted_q_max_per_m,
            "n_modes": self.n_modes,
            "n_samples": self.n_samples,
            "confidence": self.confidence,
            "seed": self.seed,
            "scope_note": self.scope_note,
            "evidence_source": EVIDENCE_SOURCE.value,
        }


def measure_cutoff_damage_ladder(
    spectrum: HelfrichSpectrum,
    *,
    q: Sequence[float] | np.ndarray,
    q_cutoffs: Sequence[float],
    n_samples: int,
    seed: int,
    confidence: float,
    min_modes: int = 5,
    rolloff_exponent: float = 4.0,
) -> CutoffDamageLadder:
    """Measure the kappa inflation as a function of mesh-cutoff position, not at one position.

    For each cutoff the true spectrum is truncated by the same roll-off, amplitudes are drawn from the
    *same* seed so the cutoff position is the only thing that varies, and two fits are compared: the
    naive fit over the whole band, and the fit taken on the prefix :func:`trusted_wavevector_range`
    declares trusted.  The comparison is the deliverable — see
    :data:`CUTOFF_DAMAGE_IS_NOT_A_MULTIPLE` for what a single entry of it may be said to mean.

    ``guarded_kappa_ratio`` is ``None`` at any cutoff where the guard trusts fewer modes than a
    two-parameter fit can support; that is a refusal, and it is recorded rather than filled in.

    Args:
        spectrum: the true spectrum the amplitudes are drawn from.
        q: the fitted band, ascending, in 1/m.
        q_cutoffs: cutoff wavevectors to test, in 1/m.  The caller chooses these; the function
            supplies no ladder of its own, because a default ladder would smuggle in a preferred
            answer.
        n_samples: configurations behind each per-mode variance.
        seed: RNG seed; part of the measurement record.
        confidence: goodness-of-fit level handed to the guard.  No default: it is a pre-declared
            threshold and the caller must state it.
        min_modes: smallest prefix the guard may report as trusted.
        rolloff_exponent: exponent of the ``exp(-(q/q_cutoff)**p)`` discretisation roll-off.

    Returns:
        A :class:`CutoffDamageLadder`.

    Raises:
        ValueError: if ``q_cutoffs`` is empty.
    """
    if not isinstance(spectrum, HelfrichSpectrum):
        raise TypeError("spectrum must be a HelfrichSpectrum")
    if spectrum.kappa_J <= 0.0:
        raise ValueError(
            "a kappa inflation ratio is undefined for a spectrum whose true kappa is zero"
        )
    wavevectors = _positive_array(q, what="q")
    if np.any(np.diff(wavevectors) <= 0.0):
        raise ValueError("q must be strictly ascending")
    cutoffs = tuple(float(_positive_float(value, what="q_cutoffs entry")) for value in q_cutoffs)
    if not cutoffs:
        raise ValueError("q_cutoffs must not be empty")
    count = _positive_int(n_samples, what="n_samples", minimum=2)
    level = _confidence(confidence)
    truth = spectrum.mode_variance(wavevectors)
    q_top = float(wavevectors[-1])
    common = {
        "area_m2": spectrum.area_m2,
        "temperature_K": spectrum.temperature_K,
        "n_samples": count,
    }

    naive_ratios: list[float] = []
    guarded_ratios: list[float | None] = []
    trusted_counts: list[int] = []
    for cutoff in cutoffs:
        truncated = mesh_truncated_variances(
            wavevectors, truth, q_cutoff=cutoff, rolloff_exponent=rolloff_exponent
        )
        rng = np.random.default_rng(_positive_int(seed, what="seed", minimum=0))
        measured = measure_mode_variance(sample_gaussian_modes(truncated, n_samples=count, rng=rng))
        naive = fit_helfrich_spectrum(wavevectors, measured, confidence=level, **common)
        naive_ratios.append(naive.kappa_J / spectrum.kappa_J)
        report = trusted_wavevector_range(
            wavevectors, measured, confidence=level, min_modes=min_modes, **common
        )
        trusted_counts.append(report.n_modes_trusted)
        if report.n_modes_trusted >= 3:
            guarded = fit_helfrich_spectrum(
                wavevectors[: report.n_modes_trusted],
                measured[: report.n_modes_trusted],
                confidence=level,
                **common,
            )
            guarded_ratios.append(guarded.kappa_J / spectrum.kappa_J)
        else:
            guarded_ratios.append(None)

    return CutoffDamageLadder(
        q_cutoffs_per_m=cutoffs,
        cutoff_over_fitted_q_max=tuple(cutoff / q_top for cutoff in cutoffs),
        decades_of_band_above_cutoff=tuple(math.log10(q_top / cutoff) for cutoff in cutoffs),
        naive_kappa_ratio=tuple(naive_ratios),
        guarded_kappa_ratio=tuple(guarded_ratios),
        n_modes_trusted=tuple(trusted_counts),
        fitted_q_min_per_m=float(wavevectors[0]),
        fitted_q_max_per_m=q_top,
        n_modes=int(wavevectors.size),
        n_samples=count,
        confidence=level,
        seed=int(seed),
    )


# --------------------------------------------------------------------------------------------
# verdict and failure routing
# --------------------------------------------------------------------------------------------


class FluctuationFailure(StrEnum):
    """The distinct ways an S1 comparison can fail; each routes to a different expansion."""

    WRONG_EXPONENT = "wrong-exponent"
    EQUIPARTITION_VIOLATED = "equipartition-violated"
    MODES_CORRELATED = "modes-correlated"
    DISCRETISATION_CUTOFF = "discretisation-cutoff"
    PARAMETER_MISMATCH = "parameter-mismatch"
    SPAN_TOO_NARROW = "span-too-narrow"


#: Where each failure sends the programme.  A failure is never a reason to delete scope.
FAILURE_EXPANSIONS: dict[FluctuationFailure, str] = {
    FluctuationFailure.WRONG_EXPONENT: (
        "stochastic boundary law — the restoring operator, not its coefficients, is wrong; open a "
        "boundary-law project rather than refitting kappa and sigma"
    ),
    FluctuationFailure.EQUIPARTITION_VIOLATED: (
        "stochastic boundary law — the fluctuation-dissipation balance of the thermostat/boundary is "
        "wrong; the noise amplitude and the drag are not the same operator"
    ),
    FluctuationFailure.MODES_CORRELATED: (
        "stochastic boundary law — correlated noise injection or a non-diagonal mobility; open a "
        "correlated-noise / mobility-operator project"
    ),
    FluctuationFailure.DISCRETISATION_CUTOFF: (
        "numerical-invariance battery (S0) — resolve the mesh and re-measure before any physics "
        "claim; the high-q band is a discretisation artifact, not an observable"
    ),
    FluctuationFailure.PARAMETER_MISMATCH: (
        "missing component / missing law — the spectral SHAPE is right but the recovered moduli are "
        "not the input, so something else contributes to the restoring force"
    ),
    FluctuationFailure.SPAN_TOO_NARROW: (
        "adaptive design — the measurement cannot discriminate at this patch size; enlarge the patch "
        "or the sample count before drawing any conclusion"
    ),
}


@dataclass(frozen=True, slots=True)
class FluctuationVerdict:
    """Conjunctive adjudication of a measured spectrum against the analytic canon."""

    passed: bool
    failures: tuple[FluctuationFailure, ...]
    expansions: tuple[str, ...]
    fit: SpectrumFit
    exponent: ExponentFit | None
    equipartition: EquipartitionReport | None
    independence: ModeIndependenceReport | None
    trusted_range: TrustedRangeReport
    evidence_source: EvidenceSource = EVIDENCE_SOURCE
    evidence_authority: str = EVIDENCE_AUTHORITY

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-able adjudication record."""
        return {
            "passed": self.passed,
            "failures": [failure.value for failure in self.failures],
            "expansions": list(self.expansions),
            "fit": self.fit.to_dict(),
            "exponent": None if self.exponent is None else self.exponent.exponent,
            "equipartition_ratio": (
                None if self.equipartition is None else self.equipartition.ratio
            ),
            "modes_independent": (
                None if self.independence is None else self.independence.independent
            ),
            "q_trust_max": self.trusted_range.q_trust_max,
            "cutoff_detected": self.trusted_range.cutoff_detected,
            "evidence_source": self.evidence_source.value,
            "evidence_authority": self.evidence_authority,
        }


def adjudicate_fluctuation_spectrum(
    q: Sequence[float] | np.ndarray,
    amplitudes: np.ndarray,
    *,
    reference: HelfrichSpectrum,
    confidence: float = 0.99,
    exponent_confidence: float = 0.95,
    parameter_tolerance: float = 0.25,
    min_modes: int = 5,
    nyquist_q: float | None = None,
) -> FluctuationVerdict:
    """Run the whole S1 comparison and return a conjunctive verdict with failure routing.

    The order matters and is not cosmetic.  The trusted range is established *first*, because every
    later statistic is meaningless if it is computed over a band the grid cannot represent.  Only then
    are the exponent, the moduli, equipartition and mode independence tested, all on the trusted
    modes.  Any single failure fails the verdict — a spectrum is not "mostly right".

    Args:
        q: wavevectors in 1/m, ascending.
        amplitudes: ``(n_samples, n_modes)`` real mode amplitudes in m.
        reference: the independently determined spectrum the data claim to realise.
        confidence: level used for the goodness-of-fit, equipartition and independence tests.
        exponent_confidence: level used for the exponent interval.
        parameter_tolerance: relative tolerance on the recovered ``kappa`` and ``sigma`` against the
            reference, applied only to a parameter the fit actually resolves.
        min_modes: smallest trusted prefix.
        nyquist_q: optional known grid Nyquist wavevector, recorded not used.

    Returns:
        A :class:`FluctuationVerdict`.
    """
    if not isinstance(reference, HelfrichSpectrum):
        raise TypeError("reference must be a HelfrichSpectrum")
    wavevectors = _positive_array(q, what="q")
    samples = np.asarray(amplitudes, dtype=float)
    if samples.ndim != 2 or samples.shape[1] != wavevectors.size:
        raise ValueError("amplitudes must be (n_samples, n_modes) matching q")
    n_samples = int(samples.shape[0])
    tolerance = _positive_float(parameter_tolerance, what="parameter_tolerance")
    level = _confidence(confidence)

    measured = measure_mode_variance(samples)
    trusted = trusted_wavevector_range(
        wavevectors,
        measured,
        area_m2=reference.area_m2,
        temperature_K=reference.temperature_K,
        n_samples=n_samples,
        confidence=level,
        min_modes=min_modes,
        nyquist_q=nyquist_q,
    )

    failures: list[FluctuationFailure] = []
    if trusted.cutoff_detected:
        failures.append(FluctuationFailure.DISCRETISATION_CUTOFF)

    keep = max(trusted.n_modes_trusted, min_modes)
    keep = min(keep, wavevectors.size)
    q_trusted = wavevectors[:keep]
    v_trusted = measured[:keep]
    a_trusted = samples[:, :keep]

    fit = fit_helfrich_spectrum(
        q_trusted,
        v_trusted,
        area_m2=reference.area_m2,
        temperature_K=reference.temperature_K,
        n_samples=n_samples,
        confidence=level,
    )

    exponent: ExponentFit | None = None
    if keep >= 3:
        exponent = fit_power_law_exponent(
            q_trusted, v_trusted, n_samples=n_samples, confidence=exponent_confidence
        )
        span_regimes = {
            reference.regime_at(float(value)) for value in (q_trusted[0], q_trusted[-1])
        }
        if span_regimes == {FluctuationRegime.BENDING_DOMINATED}:
            if exponent.excludes(-4.0):
                failures.append(FluctuationFailure.WRONG_EXPONENT)
        elif span_regimes == {FluctuationRegime.TENSION_DOMINATED} and exponent.excludes(-2.0):
            failures.append(FluctuationFailure.WRONG_EXPONENT)

    equipartition = check_equipartition(
        q_trusted,
        v_trusted,
        spectrum=reference,
        n_samples=n_samples,
        confidence=level,
    )
    if not equipartition.consistent:
        failures.append(FluctuationFailure.EQUIPARTITION_VIOLATED)

    independence: ModeIndependenceReport | None = None
    if keep >= 2 and n_samples >= 5:
        independence = check_mode_independence(a_trusted, confidence=level)
        if not independence.independent:
            failures.append(FluctuationFailure.MODES_CORRELATED)

    if (
        fit.kappa_resolved
        and reference.kappa_J > 0.0
        and abs(fit.kappa_J - reference.kappa_J) / reference.kappa_J > tolerance
    ):
        failures.append(FluctuationFailure.PARAMETER_MISMATCH)
    if (
        FluctuationFailure.PARAMETER_MISMATCH not in failures
        and fit.sigma_resolved
        and reference.sigma_N_per_m > 0.0
        and abs(fit.sigma_N_per_m - reference.sigma_N_per_m) / reference.sigma_N_per_m > tolerance
    ):
        failures.append(FluctuationFailure.PARAMETER_MISMATCH)

    if not fit.kappa_resolved and not fit.sigma_resolved:
        failures.append(FluctuationFailure.SPAN_TOO_NARROW)

    ordered = tuple(dict.fromkeys(failures))
    return FluctuationVerdict(
        passed=not ordered,
        failures=ordered,
        expansions=tuple(FAILURE_EXPANSIONS[failure] for failure in ordered),
        fit=fit,
        exponent=exponent,
        equipartition=equipartition,
        independence=independence,
        trusted_range=trusted,
    )


# --------------------------------------------------------------------------------------------
# pre-registration
# --------------------------------------------------------------------------------------------


def s1_experiment_card(
    manifest_hashes: Sequence[str], *, budget_class: str = "cpu-oracle"
) -> SandboxExperimentCard:
    """Return the pre-registration card for sandbox S1.

    The card is written before any native amplitude exists, which is the only time a falsifier can be
    written honestly.  Its ``failure_expansions`` are the routes in :data:`FAILURE_EXPANSIONS`, headed
    by the stochastic boundary law project the master plan names for this row.

    Args:
        manifest_hashes: cell-state manifest hashes the experiment applies to.
        budget_class: declared cost class.

    Returns:
        A :class:`~aleph.virtual_cell.contracts.SandboxExperimentCard`.
    """
    return SandboxExperimentCard(
        experiment_id="S1-passive-fluctuation-canon",
        question=(
            "does the simulated membrane's thermal fluctuation reproduce the analytic Helfrich mode "
            "spectrum kB*T/(A*(kappa*q^4 + sigma*q^2)) MODE BY MODE, at the physiological "
            "temperature and with the modes independent?"
        ),
        manifest_hashes=tuple(manifest_hashes),
        intervention=(
            "no active intervention: the membrane is held at thermal equilibrium with no motor, no "
            "pressure step and no adhesion, and only the temperature and the declared moduli are set"
        ),
        observables=(
            "per-mode mean square amplitude <|h_q|^2>",
            "log-log spectral exponent over a stated q decade span",
            "recovered kappa and sigma with confidence intervals",
            "implied kB*T from independently known moduli",
            "pairwise mode correlation matrix",
            "largest q whose spectrum passes the declared goodness-of-fit",
        ),
        positive_controls=(
            "amplitudes drawn from the analytic Gaussian: the estimator must recover kappa and sigma "
            "inside their stated intervals",
            "1-D taut string with fixed ends: mean energy per mode must equal kB*T/2, an identity "
            "with no approximation anywhere",
        ),
        negative_controls=(
            "wrong power law: a q^-2 spectrum must be REJECTED by the exponent test rather than "
            "absorbed by refitting kappa and sigma",
            "equipartition violated by a factor: amplitudes scaled so the implied kB*T is wrong must "
            "be caught against independently determined moduli",
            "correlated modes: modes with correct variances but nonzero covariance must be caught by "
            "an independence test, which a per-mode variance check cannot do",
            "discretisation cutoff: a mesh-truncated spectrum must be excluded from the fit rather "
            "than fitted as a stiffer membrane",
        ),
        adversarial_controls=(
            "narrow-band fit: over a span shorter than the measured separation requirement, a q^-2 "
            "spectrum is fitted by the Helfrich form at excellent chi-square, so fit quality alone "
            "must never be reported as agreement",
            "global amplitude rescale with kappa and sigma free: exactly degenerate, so a free "
            "two-parameter fit must never be presented as an equipartition check",
        ),
        held_out_interventions=(
            "a second temperature not used in any fit: the spectrum must scale linearly in T",
            "a second patch area not used in any fit: <|h_q|^2> must scale as 1/A at fixed q",
            "an imposed tension step: q* = sqrt(sigma/kappa) must move as the square root of it",
        ),
        falsifier=(
            "per-mode mismatch: on the trusted q range, either the fitted exponent interval excludes "
            "the regime's expected value, or the implied kB*T interval excludes the set temperature, "
            "or the mode-independence test rejects at the declared family-wise level"
        ),
        failure_expansions=tuple(dict.fromkeys(FAILURE_EXPANSIONS.values())),
        budget_class=budget_class,
        wet_lab=False,
        information_gain_rationale=None,
    )
