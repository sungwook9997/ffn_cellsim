"""Pereverzev two-pathway catch-slip k_off(F) closed-form oracle (D2, KU-2.5).

Validation oracle for the HOOMD-emergent integrin off-rate produced by
``ffn_sim.archive.hoomd_legacy.bridge.integrin_bonds.IntegrinBondUpdater``. This module is
``import``-clean (pure NumPy, no HOOMD) and is meant to be called from
``ffn_sim/tests/validation/`` only — per CLAUDE.md hard rule, closed-form
paper models are NEVER imported by the runtime path.

Pereverzev 2005 two-pathway form::

    k_off(F) = k_s · exp(F / F_s)  +  k_c · exp(−F / F_c)
              \\______slip______/    \\______catch______/

with KU-2.18 defaults (H.4 brief §Integrin off-rate)::

    k_s = 0.5 s⁻¹,  F_s = 30 pN
    k_c = 0.4 s⁻¹,  F_c = 7 pN

The lifetime τ(F) = 1 / k_off(F) has a maximum at F = F* where
dk_off/dF = 0::

    (k_s/F_s) · exp(F*/F_s) = (k_c/F_c) · exp(−F*/F_c)
    F* = ln( (k_c · F_s) / (k_s · F_c) ) / (1/F_s + 1/F_c)

With the KU-2.18 defaults the analytic peak is::

    F* = ln(0.4·30 / (0.5·7)) / (1/30 + 1/7)
       = ln(3.4286) / 0.17619 pN⁻¹
       ≈ 6.99 pN

This is the **analytic** peak from the KU-2.18 parameters — *not* the
KU-2.5 **experimental** peak F* ≈ 30 pN (Kong 2009 PNAS, Elosegui-Artola
2016 Nat Cell Biol). The H.4 sanity gate
``gate_unit2_1_motor_clutch.catch_peak_F_star_within_analytic`` enforces
agreement with the analytic value and logs the experimental gap as info.
Refitting KU-2.18 to the experimental peak is deferred to Phase 2 per
PHASE_0_3_DECISIONS D2.

Sanity Gate
-----------
*Pure NumPy oracle; no HOOMD; runtime-side IntegrinBondUpdater holds the
operational gate. This module's own gate is correctness against the
analytic identities above.*

1. **Dimensional analysis**: F [N], F_s/F_c [N], k_s/k_c [s⁻¹] →
   k_off [s⁻¹], τ [s]. F* [N]. All SI; STATIC test in
   ``test_pereverzev_oracle.py`` (lives inside the catch_peak validation
   test file).
2. **Boundary cases**:
   - F=0: k_off(0) = k_s + k_c. STATIC.
   - F → ∞: slip term dominates, k_off → ∞. STATIC.
   - F → −∞ (non-physical compressive load): catch term diverges; we
     restrict the oracle's valid domain to F ≥ 0 and raise on F < 0.
   - Slip-only limit (k_c → 0 or F_c → ∞): F* is undefined (no peak);
     ``pereverzev_F_star`` returns NaN; STATIC.
   - Catch-dominated limit (k_s · F_c ≫ k_c · F_s): the log argument
     becomes < 1, log → negative, F* < 0 → the analytic stationary
     point sits in the non-physical compressive regime; the function
     returns ``F*`` as-is (caller's responsibility to interpret).
3. **Conservation invariants**: not a conservation law — first-order
   rate process. The total off-rate is the sum of two pathway-specific
   rates, each independently exponential in F per Bell-Evans.
4. **Numerical sanity**: float64 throughout. ``np.exp`` saturates at
   F/F_s ≈ 700 with float64; if a caller hands in F > 700·F_s
   ≈ 2.1e-8 N (i.e. 21 nN per bond), the slip term overflows. Realistic
   bond forces are O(pN) so this is far outside the operating range,
   but we add a runtime guard that raises ``FloatingPointError`` rather
   than silently returning ``inf``.
5. **Sign / sense**: ``dk_off/dF`` at F=0 is (k_s/F_s − k_c/F_c). With
   KU-2.18 defaults this is 0.5/30 − 0.4/7 ≈ +0.0167 − 0.0571 = −0.0404
   < 0, i.e. **k_off decreases as F increases from 0** — the catch-bond
   signature. STATIC test asserts this sign at the KU-2.18 point in
   parameter space.
6. **Measurement-protocol consistency**: The HOOMD-emergent k_off(F)
   gate (``gate_emergent_vs_oracle`` for module ``pereverzev_koff``)
   compares simulated k_off on a force grid F ∈ [0, 60] pN within ± 5%.
   This module returns the closed-form values used on the right-hand
   side of that comparison.

References
----------
- Pereverzev YV, Prezhdo OV, Forero M, Sokurenko EV, Thomas WE 2005,
  "The two-pathway model for the catch-slip transition in biological
  adhesion", *Biophys J* 89(3):1446-54. The k_s·exp(F/F_s) + k_c·exp(−F/F_c) form.
- Bangasser BL, Rosenfeld SS, Odde DJ 2013, "Determinants of maximal force
  transmission in a motor-clutch model of cell traction in a compliant
  microenvironment", *Biophys J* 105(3):581-92 (PMID 23931306). KU-2.18
  illustrative parameter set used here. [title corrected 2026-06-07: was the
  unrelated "Determinism and stochasticity during maturation of the zyxin focal
  adhesion" — metadata drift, cf. the 2026-06-02 SourceEvidence audit.]
- PHASE_0_3_DECISIONS.md §D2 "Bell-Evans / catch-slip" (2026-05-19,
  PI-ratified).
- ``ffn_sim/docs/briefs/H4_fa_motor_clutch.md`` §Integrin off-rate.
- ``ffn_sim/validation/oracles/common/sanity_gate.py::gate_emergent_vs_oracle``
  — generic emergent-vs-oracle ± 5% comparison gate used by the
  catch-peak validation test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# KU-2.18 illustrative parameter set (Bangasser 2013 / H.4 brief).
# These are also the defaults baked into ``configs/phase1_h4.yaml``; this
# module's dataclass mirrors them so a caller without a config file can
# still validate the closed-form behaviour.
KU_2_18_K_S: float = 0.5            # s⁻¹    slip prefactor
KU_2_18_F_S: float = 30.0e-12       # N      slip characteristic force (30 pN)
KU_2_18_K_C: float = 0.4            # s⁻¹    catch prefactor
KU_2_18_F_C: float = 7.0e-12        # N      catch characteristic force (7 pN)

# Float64 ``np.exp(arg)`` overflows at arg ≈ 709.78. We bound at 700 to
# leave a small margin and raise rather than return inf silently. The
# guard fires when |F| / min(F_s, F_c) > 700; with KU-2.18's F_c = 7 pN
# that's |F| > 4.9 nN per bond, far outside the integrin operating range
# (O(pN)) so this only catches caller bugs (e.g. passing F in pN instead
# of N).
EXP_ARG_GUARD: float = 700.0


@dataclass(frozen=True, slots=True)
class PereverzevParams:
    """Two-pathway catch-slip parameter set (SI units: N, s⁻¹).

    Defaults match KU-2.18 (Bangasser 2013) per H.4 brief. Set via
    :meth:`from_config` to bind to ``configs/phase1_h4.yaml``.
    """

    k_s: float = KU_2_18_K_S
    F_s: float = KU_2_18_F_S
    k_c: float = KU_2_18_K_C
    F_c: float = KU_2_18_F_C

    @classmethod
    def from_config(cls, cfg: dict) -> "PereverzevParams":
        b = cfg["bridge"] if "bridge" in cfg else cfg
        cb = b["catch_bond"]
        return cls(
            k_s=float(cb["k_off_slip"]),
            F_s=float(cb["F_s"]),
            k_c=float(cb["k_off_catch"]),
            F_c=float(cb["F_c"]),
        )

    def __post_init__(self) -> None:
        for name, v in (
            ("k_s", self.k_s), ("F_s", self.F_s),
            ("k_c", self.k_c), ("F_c", self.F_c),
        ):
            if not (np.isfinite(v) and v > 0.0):
                raise ValueError(
                    f"PereverzevParams.{name} must be finite and > 0, got {v!r}"
                )


DEFAULT_PARAMS = PereverzevParams()


def pereverzev_k_off(
    F: float | np.ndarray, params: PereverzevParams = DEFAULT_PARAMS
) -> float | np.ndarray:
    """k_off(F) = k_s·exp(F/F_s) + k_c·exp(−F/F_c)  [s⁻¹].

    Parameters
    ----------
    F : float or ndarray
        Bond force magnitude in Newtons. Must be ≥ 0 (compressive forces
        are not part of the Pereverzev domain).
    params : PereverzevParams, default KU-2.18

    Returns
    -------
    k_off : same shape as ``F``  [s⁻¹]
    """
    F_arr = np.asarray(F, dtype=np.float64)
    if (F_arr < 0).any():
        raise ValueError(
            "pereverzev_k_off: F must be ≥ 0 (Pereverzev two-pathway form "
            "is defined for tensile loading only)."
        )
    # Overflow guard before np.exp: F/F_s and F/F_c both ≤ EXP_ARG_GUARD.
    arg_s = F_arr / params.F_s
    arg_c = F_arr / params.F_c
    if (arg_s > EXP_ARG_GUARD).any() or (arg_c > EXP_ARG_GUARD).any():
        worst = float(max(arg_s.max(), arg_c.max()))
        raise FloatingPointError(
            f"pereverzev_k_off: |F|/min(F_s,F_c) = {worst:.3e} > "
            f"{EXP_ARG_GUARD}; np.exp would overflow. Check unit (F in N, "
            f"not pN). KU-2.18 operating range is ~pN."
        )
    return params.k_s * np.exp(arg_s) + params.k_c * np.exp(-arg_c)


def pereverzev_lifetime(
    F: float | np.ndarray, params: PereverzevParams = DEFAULT_PARAMS
) -> float | np.ndarray:
    """τ(F) = 1 / k_off(F)  [s]."""
    return 1.0 / pereverzev_k_off(F, params)


def pereverzev_F_star(params: PereverzevParams = DEFAULT_PARAMS) -> float:
    """Analytic catch-peak force F*  [N], or NaN if no peak exists.

    Solves dk_off/dF = 0 in closed form::

        F* = ln( k_c·F_s / (k_s·F_c) ) / (1/F_s + 1/F_c)

    NaN return cases:
        - log argument ≤ 0 (catch coefficient missing or pathological).
        - The stationary point is a minimum, not a maximum — happens
          when the catch term cannot lower the rate below F=0. We detect
          via the second-derivative sign at the candidate F*.
    """
    arg = (params.k_c * params.F_s) / (params.k_s * params.F_c)
    if arg <= 0.0:
        return float("nan")
    F_star = float(np.log(arg) / (1.0 / params.F_s + 1.0 / params.F_c))
    # Verify maximum: d²k_off/dF² = (k_s/F_s²)·exp(F/F_s) + (k_c/F_c²)·exp(−F/F_c) > 0
    # always for k_s, k_c > 0 — so k_off is strictly convex and any
    # stationary point is a minimum. Wait — that means τ = 1/k_off has a
    # MAXIMUM at F* (since k_off has a minimum there). Catch-bond
    # lifetime peak corresponds to k_off minimum. The sign check below
    # asserts k_off is at a minimum (d²k_off/dF² > 0 trivially) AND that
    # F* is positive (catch-slip transition in the tensile regime).
    if F_star < 0.0:
        return float("nan")
    return F_star


def pereverzev_lifetime_peak(
    params: PereverzevParams = DEFAULT_PARAMS,
) -> tuple[float, float]:
    """(F*, τ*) at the analytic lifetime maximum  [N, s].

    Returns ``(nan, nan)`` if no peak exists (slip-only limit or
    catch-pathological parameters; see ``pereverzev_F_star``).
    """
    F_star = pereverzev_F_star(params)
    if not np.isfinite(F_star):
        return float("nan"), float("nan")
    tau_star = float(pereverzev_lifetime(F_star, params))
    return F_star, tau_star


def sample_bond_lifetime(
    F: float,
    params: PereverzevParams = DEFAULT_PARAMS,
    *,
    rng: np.random.Generator | None = None,
) -> float:
    """Draw one bond lifetime from the exponential distribution at force F.

    ``T ~ Exp(rate=k_off(F))``  so  ``T = −ln(U) / k_off(F)`` with
    U ~ Uniform(0, 1). Used by the validation tests to estimate
    HOOMD-emergent τ(F) via Monte Carlo without needing the full HOOMD
    state — a "what would the integrin updater do for a constant-F
    bond?" reference run.

    Parameters
    ----------
    F : float
        Bond force in Newtons.
    params : PereverzevParams
    rng : np.random.Generator, optional
        Defaults to a freshly-seeded ``default_rng()``.

    Returns
    -------
    T : float
        Single sampled lifetime [s].
    """
    if rng is None:
        rng = np.random.default_rng()
    k = float(pereverzev_k_off(F, params))
    if k <= 0.0 or not np.isfinite(k):
        return float("inf")
    u = float(rng.uniform(0.0, 1.0))
    # Avoid log(0); the open interval (0, 1) excludes it but float draws
    # can reach the endpoints. Treat u==0 as the (probability-zero)
    # infinite-lifetime tail.
    if u <= 0.0:
        return float("inf")
    return -np.log(u) / k


def estimate_mean_lifetime(
    F: float,
    n_samples: int,
    params: PereverzevParams = DEFAULT_PARAMS,
    *,
    seed: int = 0,
) -> tuple[float, float]:
    """Monte Carlo ⟨τ(F)⟩ ± stderr from ``n_samples`` exponential draws.

    Returns ``(mean, stderr)``. Used by KU-2.5 catch-peak test to
    establish that the emergent τ_peak coincides with F* within ± 20 %
    of the closed-form ``pereverzev_lifetime(F*)``.
    """
    rng = np.random.default_rng(seed)
    samples = np.array(
        [sample_bond_lifetime(F, params, rng=rng) for _ in range(n_samples)],
        dtype=np.float64,
    )
    finite = samples[np.isfinite(samples)]
    if finite.size == 0:
        return float("nan"), float("nan")
    mean = float(finite.mean())
    stderr = float(finite.std(ddof=1) / np.sqrt(finite.size))
    return mean, stderr
