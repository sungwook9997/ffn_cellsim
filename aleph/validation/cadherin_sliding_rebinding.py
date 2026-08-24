"""E-cadherin X-dimer catch-slip bond — Rakshit 2012 sliding-rebinding model (L2.5).

Runtime-importable closed form (pure NumPy), placed at the ``validation/`` level exactly like
``validation/pereverzev.py`` (the catch-slip k_off that ``bridge/integrin_bonds.py`` imports)
— NOT under ``validation/oracles/`` (that dir is runtime-import-forbidden). Faithful
implementation of the sliding-rebinding catch-bond model that Rakshit et al. 2012 fit to
E-cadherin X-dimer single-molecule force spectroscopy — the PI-chosen "option B" (the
two-pathway/Pereverzev form did NOT fit their data; the sliding-rebinding model did). The
runtime ``the retired HOOMD tree (deleted 1a9ded66; `git show 1a9ded66^:`).spheroid.cadherin_bonds.CadherinBondUpdater`` samples bond lifetimes from
this module's force-dependent off-rate, and it is also the closed-form right-hand side of the
emergent-vs-oracle gate (G5).

Model (Rakshit 2012 SI, three transient states + absorbing dissociated state)
-----------------------------------------------------------------------------
States: ``P11`` (X-dimer, both pseudoatom pairs bound), ``P10`` (one pair dissociated,
pre-slide — lumps P10/P01), ``P10p`` (P10', a new interaction formed after sliding), and the
absorbing ``P00`` (dissociated). Survival probability S(t) = P11 + P10 + P10p. Rate equations
(Rakshit 2012 SI Eq. set, Fig. S1):

    dP11 /dt = 2 k+1 P10  + k+2 P10p − k-2 P11
    dP10 /dt = k-2 P11 − 2(k+1 + k-1) P10
    dP10p/dt = 2 Pn k-1 P10 − (k+2 + k-1) P10p

with the force-dependent single-pair off-rate and the shared-force two-pair off-rate

    k-1(f) = k-1^0 · exp(+f · x / kB T)          (Bell SLIP; x = bound→transition distance)
    k-2(f) = 2 · k-1(f/2)                         (two pairs share the load)

and the new-interaction (sliding→rebinding) probability

    Pn(f) = 0                                   (f < 0)
          = { 0.5 [1 + sin(π f / f0 − π/2)] }^n  = sin^{2n}(π f / 2 f0)   (0 ≤ f ≤ f0)
          = 1                                   (f > f0)

Pn rises 0→1 over [0, f0]: under tension the X-dimer interface slides and rebinds, *creating*
interactions faster than they break → the bond LIVES LONGER (catch) up to f0≈29 pN, beyond
which Pn saturates and Bell dissociation (k-1 ↑ with force) wins → it shortens (slip). That
catch regime is the tension proliferation generates, so a catch-bond cohesion should resist
the L2.4 growth-driven fragmentation that the static Morse well could not.

Off-rate sign (documented inference)
------------------------------------
The single-pair off-rate uses the Bell SLIP form ``k-1 = k-1^0 · exp(+f·x/kBT)`` (off-rate
INCREASES with force), matching the source model **Lou & Zhu 2007 Eq. 1**
(``k-1 = k-1^0 exp(a f / n kB T)``, positive). The Rakshit SI renders the exponent with a
minus sign, but that is unphysical here: with a negative exponent k-1 → 0 under load and, with
Pn → 1 and the large rebinding rate k+2, the lifetime DIVERGES (pure catch, no slip) — it
cannot reproduce the measured biphasic peak. The positive (slip) exponent is required for the
catch (from sliding-rebinding Pn) to hand off to slip beyond f0, and it reproduces the measured
catch peak at f ≈ 29 pN (test). Treated as a sign-convention/typo in the rendered SI.

Faithful-transcription note (documented inference)
--------------------------------------------------
The SI's rendered dP10p/dt sliding term reads ambiguously as ``k+1`` vs ``k-1`` (the ±1
subscript is not legible in the rendered PDF). Probability conservation settles it: with the
stated P10 diagonal −2(k+1 + k-1), the outflows of P10 are 2 k+1 (rebind → P11) + 2 k-1
(dissociate), and the dissociation flux 2 k-1 splits Pn → P10p and (1−Pn) → P00. Hence the
slide term is ``2 Pn k-1 P10`` (NOT k+1): only then are the column leaks to P00 non-negative
(P11: 0, P10: 2(1−Pn)k-1, P10p: k-1) and the chain a proper absorbing one. This implementation
is validated to reproduce the measured catch peak at f ≈ f0 ≈ 29 pN (test).

Mean bond lifetime
------------------
For the linear absorbing system dP/dt = M·P (P = [P11, P10, P10p], P00 absorbing), the mean
time to dissociation from the freshly-bound X-dimer (P(0) = [1, 0, 0]) is the exact

    τ(f) = ∫_0^∞ S(t) dt = 1ᵀ ∫_0^∞ P dt = −1ᵀ M⁻¹ P(0)

(no numerical integration). The effective single-bond off-rate the updater samples is
k_off(f) = 1 / τ(f).

Provenance
----------
Rakshit S, Zhang Y, Manibog K, Shafraz O, Sivasankar S (2012) "Ideal, catch, and slip bonds
in cadherin adhesion." PNAS 109(46):18815-18820. DOI 10.1073/pnas.1208349109 (PMC3503169);
model parameters from SI Table S1 (W2A E-cadherin X-dimer fit). Sliding-rebinding model:
Lou J, Zhu C (2007) Biophys J 92(5):1471-1485, DOI 10.1529/biophysj.106.097048 (PMC1796828).
Both retrieved via PubMed; PDFs in aleph/references/.

Sanity Gate
-----------
- Dimensional: f [N], x [m], kB T [J] → k-1 [s⁻¹]; rates [s⁻¹]; τ [s]; k_off [s⁻¹].
- Boundary: f=0 → Pn=0, pure decay, τ(0) finite > 0. f→∞ → Pn=1, k-1→0 (x>0), slip via the
  finite rebinding pool. f<0 rejected (tensile domain only).
- Conservation: M column leaks to P00 are all ≥ 0 (asserted in tests); M is Hurwitz so τ>0.
- Catch signature: τ(f) increases from f=0, peaks near f0≈29 pN, then decreases (test).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

__all__ = [
    "CadherinCatchParams",
    "RAKSHIT_W2A",
    "new_interaction_probability",
    "single_pair_off_rate",
    "mean_lifetime",
    "effective_k_off",
    "catch_peak_force",
]

# Boltzmann constant and body temperature (match the retired HOOMD tree (deleted 1a9ded66; `git show 1a9ded66^:`).bridge.ligand_species: 37 °C).
_KB: float = 1.380649e-23      # J/K (CODATA)
_T_BODY: float = 310.15        # K
_KBT: float = _KB * _T_BODY    # ≈ 4.2816e-21 J


@dataclass(frozen=True, slots=True)
class CadherinCatchParams:
    """Sliding-rebinding parameters for an E-cadherin X-dimer catch-slip bond (SI units).

    Defaults = Rakshit 2012 SI Table S1 (W2A E-cadherin X-dimer fit). These are MEASURED-FIT
    acceptance constants (not tuning knobs).
    """

    k_off0: float = 30.4           # s⁻¹  intrinsic single-pair off-rate k-1^0
    x: float = 0.34e-9             # m    bound→transition distance (0.34 nm)
    k_on1: float = 5.3             # s⁻¹  rebinding rate k+1
    k_on2: float = 1985.9          # s⁻¹  rebinding rate k+2 (P10' → P11)
    f0: float = 29.2e-12           # N    catch→slip transition force (29.2 pN)
    n: float = 4.8                 # –    interfacial-angle fitting exponent
    kT: float = _KBT               # J

    def __post_init__(self) -> None:
        for name, v in (
            ("k_off0", self.k_off0), ("x", self.x), ("k_on1", self.k_on1),
            ("k_on2", self.k_on2), ("f0", self.f0), ("n", self.n), ("kT", self.kT),
        ):
            if not (np.isfinite(v) and v > 0.0):
                raise ValueError(f"CadherinCatchParams.{name} must be finite and > 0; got {v!r}")


RAKSHIT_W2A = CadherinCatchParams()


def new_interaction_probability(
    f: float | npt.NDArray[np.float64], params: CadherinCatchParams = RAKSHIT_W2A
) -> float | npt.NDArray[np.float64]:
    """Pn(f): probability a slide creates a new interaction. 0 below 0, sin^{2n} ramp, 1 above f0."""
    f_arr = np.asarray(f, dtype=np.float64)
    ramp = np.sin(np.pi * np.clip(f_arr, 0.0, params.f0) / (2.0 * params.f0)) ** (2.0 * params.n)
    pn = np.where(f_arr <= 0.0, 0.0, np.where(f_arr >= params.f0, 1.0, ramp))
    return float(pn) if np.ndim(f) == 0 else pn


def single_pair_off_rate(
    f: float | npt.NDArray[np.float64], params: CadherinCatchParams = RAKSHIT_W2A
) -> float | npt.NDArray[np.float64]:
    """Bell single-pair off-rate k-1(f) = k-1^0 · exp(+f·x / kB T)  [s⁻¹] (f ≥ 0, slip form)."""
    f_arr = np.asarray(f, dtype=np.float64)
    if np.any(f_arr < 0.0):
        raise ValueError("single_pair_off_rate: f must be ≥ 0 (tensile loading only).")
    k = params.k_off0 * np.exp(f_arr * params.x / params.kT)
    return float(k) if np.ndim(f) == 0 else k


def _generator(f: float, params: CadherinCatchParams) -> npt.NDArray[np.float64]:
    """Transient-state generator M (dP/dt = M·P, P = [P11, P10, P10p]); P00 absorbing."""
    k1 = float(single_pair_off_rate(f, params))            # k-1(f)
    k1_half = float(single_pair_off_rate(0.5 * f, params))  # k-1(f/2)
    k2 = 2.0 * k1_half                                      # k-2(f)
    pn = float(new_interaction_probability(f, params))
    kp1, kp2 = params.k_on1, params.k_on2
    return np.array(
        [
            [-k2,          2.0 * kp1,            kp2],
            [k2,          -2.0 * (kp1 + k1),     0.0],
            [0.0,          2.0 * pn * k1,       -(kp2 + k1)],
        ],
        dtype=np.float64,
    )


def mean_lifetime(
    f: float, params: CadherinCatchParams = RAKSHIT_W2A
) -> float:
    """Mean X-dimer bond lifetime τ(f) [s] = −1ᵀ M⁻¹ [1,0,0]ᵀ (exact mean time to dissociation)."""
    M = _generator(float(f), params)
    integral = np.linalg.solve(M, np.array([1.0, 0.0, 0.0]))  # ∫P dt = −M⁻¹ P(0)
    tau = -float(np.sum(integral))
    if not (np.isfinite(tau) and tau > 0.0):
        raise FloatingPointError(f"non-physical mean lifetime {tau!r} at f={f!r} (check M stability).")
    return tau


def effective_k_off(
    f: float, params: CadherinCatchParams = RAKSHIT_W2A
) -> float:
    """Effective single-bond off-rate k_off(f) = 1/τ(f) [s⁻¹] (what the runtime updater samples)."""
    return 1.0 / mean_lifetime(f, params)


def catch_peak_force(
    params: CadherinCatchParams = RAKSHIT_W2A, *, f_max_pN: float = 80.0, n_grid: int = 1601
) -> tuple[float, float]:
    """Numerically locate the lifetime catch peak: returns (F* [N], τ(F*) [s]).

    Scans τ(f) on a grid f ∈ [0, f_max_pN] pN; the maximiser is the catch peak (≈ f0 for a
    well-formed catch bond). Used by the G5 sanity gate (peak must sit near f0 ≈ 29 pN).
    """
    fs = np.linspace(0.0, f_max_pN * 1e-12, n_grid)
    taus = np.array([mean_lifetime(float(f), params) for f in fs])
    k = int(np.argmax(taus))
    return float(fs[k]), float(taus[k])
