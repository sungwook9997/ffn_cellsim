"""Osmotic / volume-response prediction-envelope oracle (dimensionless, literature-range).

ACCEPTANCE ORACLE — closed-form only, runtime-import-forbidden (CLAUDE.md). This is the
Fluid/fields component's **prediction-envelope oracle + acceptance target** for the osmotic
axis, transcribed from the PI-specified design
`docs/v2_audit/OSMOTIC_RESPONSE_MODEL_DESIGN_2026-07-23.md` (Option A: built AFTER GATE A
banks the converged in-medium resting cell as the `r = 1`, `V₀ = E₀ = O₀ = 1` reference).

It is NOT the runtime mechanism. The runtime stays the fine-grained Fluid/fields vertical
(per-solute conserved pools Nᵢ with reflection σᵢ, explicit water/solute flux, RVD/RVI as
transporter events). The mechanistic vertical must reproduce the envelopes this module emits.

Governing discipline (PI, hard): with NO experimental data the deliverable is a
**literature-range uncertainty band**, never a single deterministic value. Every scalar law
below is therefore wrapped by a parameter sweep (b, τ_reg/τ_w, R, τ_s/τ_w, α, m) whose min /
median / max form the reported envelope. Callers must present the band, not the median
(see `Envelope` and `condition_envelope`). Absolute second/minute times and a real Young's
modulus need later calibration — this is a possibility-region model, not a personalised one.

Dimensionless model (design §"Dimensionless model")
---------------------------------------------------
Reference (current medium):  V₀ = E₀ = O₀ = 1 ;  medium ratio  r = O_out / O₀.

Passive equilibrium volume — Ponder–Boyle–van 't Hoff:
    v_eq = b + (1 − b)/r                         b = osmotically-inactive fraction ∈ {0.1,0.2,0.3}

Time course (two timescales):
    v_passive(t) = v_eq + (1 − v_eq)·e^(−t/τ_w)                                (water flux only)
    v(t)         = v_passive(t) + R·(1 − v_eq)·(1 − e^(−t/τ_reg))             (+ active RVD/RVI)
    R ∈ {0, 0.5, 1}   (no / partial / full recovery)

Permeant solute (mannitol vs glycerol):
    σ(t)     = σ₀·e^(−t/τ_s)          r_eff(t) = 1 + σ(t)·(r − 1)
    mannitol/sucrose σ≈1 (sustained shrink); glycerol/urea σ₀≈1→0 (shrink then re-swell).

Shape (separate "volume ↓" from "flattened / distorted"):
    λx = λy = v^α ,  λz = v^(1−2α)      A/A₀ = v^(2α) ,  h/h₀ = v^(1−2α)
    α = 1/3 isotropic; α = 0 fixed footprint (height-only); 0<α<1/3 adherent partial.

Stiffness (sensitivity, NOT a law):
    E/E₀ = v^(−m)      m ∈ {0, 1, 2}   (m = 2 ≈ Guo et al. 2017 V^−2 in one condition; NOT universal)
    Guo et al. 2017, PMID 28973866.

Sanity Gate (verified in `_sanity()` and the pytest suite)
- Dimensional: r, v, α, m, b, R, σ all dimensionless; τ_w, τ_reg, τ_s share one time unit
  (or use t/τ_w). v_eq, v, A/A₀, h/h₀, E/E₀ pure ratios. OK.
- Boundary: r=1 → v_eq=1 (∀ b); r→∞ → v_eq→b; r→0⁺ → v_eq→∞ (swell). t=0 → v=1;
  t→∞ passive → v_eq; t→∞ with R=1 → 1 (full recovery to reference); R=0 → passive.
- Conservation: λx·λy·λz = v exactly (shape factors preserve the volume ratio).
- Sign-sense: hyper-osmotic (r>1) shrinks (v<1) and — for m>0 — stiffens (E/E₀>1). Hypo swells.
- Permeant: σ(0)=σ₀, σ(∞)=0; r_eff(0)=r, r_eff(∞)=1 → glycerol trajectory re-swells to ~1.
- Measurement-protocol: v is V/V₀ of the SAME cell relative to its in-medium reference; there
  is deliberately no "osmotic pressure = 0" control (design §Governing principle 2).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import Callable, Iterable, Mapping, NamedTuple, Sequence

import numpy as np

ArrayLike = float | Sequence[float] | np.ndarray

# ---------------------------------------------------------------------------
# Literature-range sweep grids (design §"Conditions", §"Dimensionless model").
# These are the defaults; callers may override any grid.
# ---------------------------------------------------------------------------
B_GRID: tuple[float, ...] = (0.1, 0.2, 0.3)              # osmotically-inactive fraction
R_GRID: tuple[float, ...] = (0.0, 0.5, 1.0)              # RVD/RVI recovery fraction
TAU_REG_OVER_TAU_W_GRID: tuple[float, ...] = (1.0, 3.0, 10.0)   # regulatory vs water timescale
TAU_S_OVER_TAU_W_GRID: tuple[float, ...] = (1.0, 10.0, 100.0)   # permeant vs water timescale
ALPHA_GRID: tuple[float, ...] = (0.0, 1.0 / 6.0, 1.0 / 3.0)     # fixed-footprint → isotropic
M_GRID: tuple[int, ...] = (0, 1, 2)                     # stiffness sensitivity exponent

# Minimum condition set (design §"Conditions"). r is the medium osmolality ratio.
CONDITIONS: dict[str, float] = {
    "reference": 1.0,
    "weak_hyper": 1.15,
    "strong_hyper": 1.30,
    "weak_hypo": 0.85,
    "strong_hypo": 0.70,
}


# ===========================================================================
# Core scalar closed forms (all vectorised over NumPy arrays; GPU-agnostic).
# ===========================================================================
def pbvh_v_eq(r: ArrayLike, b: float) -> np.ndarray:
    """Ponder–Boyle–van 't Hoff passive equilibrium volume ratio v_eq = V_eq/V₀.

    v_eq = b + (1 − b)/r.

    Args:
        r: medium osmolality ratio O_out/O₀ (r>1 hyper-, r<1 hypo-osmotic). Must be > 0.
        b: osmotically-inactive volume fraction ∈ [0, 1).

    Returns:
        v_eq [dimensionless]. r=1 → 1; r→∞ → b; r→0⁺ → ∞ (swelling).
    """
    r_arr = np.asarray(r, dtype=float)
    if np.any(r_arr <= 0.0):
        raise ValueError("r must be > 0 (osmolality ratio)")
    if not (0.0 <= b < 1.0):
        raise ValueError(f"b must be in [0, 1), got {b}")
    return b + (1.0 - b) / r_arr


def v_passive(t: ArrayLike, v_eq: float, tau_w: float = 1.0) -> np.ndarray:
    """Water-flux-only volume trajectory v(t) = v_eq + (1 − v_eq)·e^(−t/τ_w).

    Starts at v(0)=1 and relaxes to v_eq with the single water timescale τ_w. Use τ_w=1 for
    the dimensionless time t/τ_w (design default when no literature τ_w is swept).
    """
    if tau_w <= 0.0:
        raise ValueError("tau_w must be > 0")
    t_arr = np.asarray(t, dtype=float)
    return v_eq + (1.0 - v_eq) * np.exp(-t_arr / tau_w)


def v_active(t: ArrayLike, v_eq: float, tau_w: float = 1.0,
             recovery: float = 0.0, tau_reg: float = 3.0) -> np.ndarray:
    """Volume trajectory with active RVD/RVI regulatory recovery.

    v(t) = v_eq + (1 − v_eq)·e^(−t/τ_w) + R·(1 − v_eq)·(1 − e^(−t/τ_reg)).

    R=0 reproduces `v_passive`; R=1 recovers fully to the reference v=1 as t→∞; 0<R<1 is a
    partial regulatory volume decrease/increase. Never assert "recovers"/"doesn't" — sweep R.

    Args:
        t: time (same unit as τ_w, τ_reg) or dimensionless t/τ_w with τ_w=1.
        v_eq: passive equilibrium volume ratio (from `pbvh_v_eq`).
        tau_w: water-flux timescale.
        recovery: R ∈ [0, 1], the RVD/RVI recovery fraction.
        tau_reg: regulatory timescale (typically ≥ τ_w).
    """
    if not (0.0 <= recovery <= 1.0):
        raise ValueError(f"recovery R must be in [0, 1], got {recovery}")
    if tau_reg <= 0.0:
        raise ValueError("tau_reg must be > 0")
    t_arr = np.asarray(t, dtype=float)
    passive = v_passive(t_arr, v_eq, tau_w)
    regulatory = recovery * (1.0 - v_eq) * (1.0 - np.exp(-t_arr / tau_reg))
    return passive + regulatory


def sigma_of_t(t: ArrayLike, sigma0: float = 1.0, tau_s: float = 10.0) -> np.ndarray:
    """Effective reflection coefficient decay σ(t) = σ₀·e^(−t/τ_s) for a permeant solute.

    σ₀≈1 (impermeant limit, mannitol/sucrose) with τ_s → ∞ gives a sustained shrink; a finite
    τ_s (glycerol/urea) lets the solute equilibrate so σ→0.
    """
    if tau_s <= 0.0:
        raise ValueError("tau_s must be > 0")
    if not (0.0 <= sigma0 <= 1.0):
        raise ValueError(f"sigma0 must be in [0, 1], got {sigma0}")
    return sigma0 * np.exp(-np.asarray(t, dtype=float) / tau_s)


def r_eff_of_t(t: ArrayLike, r: float, sigma0: float = 1.0, tau_s: float = 10.0) -> np.ndarray:
    """Time-varying effective osmolality ratio r_eff(t) = 1 + σ(t)·(r − 1).

    r_eff(0)=r (full osmotic drive), r_eff(∞)=1 (solute equilibrated → drive gone).
    """
    return 1.0 + sigma_of_t(t, sigma0, tau_s) * (r - 1.0)


def v_permeant(t: ArrayLike, r: float, b: float, tau_w: float = 1.0,
               sigma0: float = 1.0, tau_s: float = 10.0,
               recovery: float = 0.0, tau_reg: float = 3.0,
               n_substep: int = 4096) -> np.ndarray:
    """Water-lagged volume trajectory for a permeant solute (shrink-then-re-swell).

    The permeant target is the quasi-static PBvH volume evaluated at the *time-varying* ratio
    r_eff(t): v_eq(t) = b + (1 − b)/r_eff(t). The cell tracks this moving target through the
    water timescale τ_w by the linear relaxation dv/dt = (v_eq(t) − v)/τ_w, integrated here
    with an exact exponential (piecewise-constant-target) integrator on a dense internal grid —
    unconditionally stable and closed-form per substep, so this stays an oracle, not a solver.

    An optional active recovery term (R, τ_reg) is added exactly as in `v_active`, driven by the
    instantaneous (1 − v_eq(t)) deficit. With σ₀=1, τ_s→∞ this collapses onto `v_active`.

    Returns v(t) sampled at the requested `t` (which need not be uniform; t≥0 assumed, v(0)=1).
    """
    t_arr = np.asarray(t, dtype=float)
    if t_arr.ndim == 0:
        t_arr = t_arr.reshape(1)
    if np.any(t_arr < 0.0):
        raise ValueError("t must be >= 0")
    t_max = float(t_arr.max())
    if t_max == 0.0:
        return np.ones_like(t_arr)

    dense = np.linspace(0.0, t_max, n_substep + 1)
    dt = dense[1] - dense[0]
    decay = np.exp(-dt / tau_w)
    # target evaluated at substep midpoints (piecewise-constant over each substep)
    mid = 0.5 * (dense[:-1] + dense[1:])
    v_target = pbvh_v_eq(r_eff_of_t(mid, r, sigma0, tau_s), b)

    v_dense = np.empty_like(dense)
    v_dense[0] = 1.0
    for k in range(n_substep):
        v_dense[k + 1] = v_target[k] + (v_dense[k] - v_target[k]) * decay

    if recovery > 0.0:
        if not (0.0 <= recovery <= 1.0):
            raise ValueError(f"recovery R must be in [0, 1], got {recovery}")
        if tau_reg <= 0.0:
            raise ValueError("tau_reg must be > 0")
        v_eq_inst = pbvh_v_eq(r_eff_of_t(dense, r, sigma0, tau_s), b)
        v_dense = v_dense + recovery * (1.0 - v_eq_inst) * (1.0 - np.exp(-dense / tau_reg))

    out = np.interp(t_arr, dense, v_dense)
    return out.reshape(np.asarray(t).shape) if np.asarray(t).ndim else float(out[0])


# ===========================================================================
# Shape + stiffness (map a scalar/array volume ratio to observables).
# ===========================================================================
class ShapeFactors(NamedTuple):
    """Anisotropic stretch factors and observables for a volume ratio v at exponent α."""
    lam_x: np.ndarray       # lateral stretch (= v^α)
    lam_y: np.ndarray       # lateral stretch (= v^α)
    lam_z: np.ndarray       # axial/height stretch (= v^(1−2α))
    area_ratio: np.ndarray  # A/A₀ = v^(2α)  (top-down footprint)
    height_ratio: np.ndarray  # h/h₀ = v^(1−2α)


def shape_factors(v: ArrayLike, alpha: float = 1.0 / 3.0) -> ShapeFactors:
    """Decompose a volume ratio v into lateral/axial stretch + footprint/height observables.

    λx = λy = v^α, λz = v^(1−2α); A/A₀ = v^(2α), h/h₀ = v^(1−2α). By construction
    λx·λy·λz = v exactly (volume-ratio conservation). α=1/3 isotropic; α=0 fixed footprint
    (height-only response); 0<α<1/3 adherent partial lateral shrink.
    """
    if not (0.0 <= alpha <= 0.5):
        raise ValueError(f"alpha must be in [0, 0.5], got {alpha}")
    v_arr = np.asarray(v, dtype=float)
    if np.any(v_arr <= 0.0):
        raise ValueError("v must be > 0")
    lam_lat = v_arr ** alpha
    lam_ax = v_arr ** (1.0 - 2.0 * alpha)
    return ShapeFactors(
        lam_x=lam_lat, lam_y=lam_lat, lam_z=lam_ax,
        area_ratio=v_arr ** (2.0 * alpha), height_ratio=lam_ax,
    )


def stiffness_ratio(v: ArrayLike, m: float = 0.0) -> np.ndarray:
    """Apparent-stiffness sensitivity E/E₀ = v^(−m). Sensitivity, NOT a validated law.

    m=0 → no stiffness change; m=2 ≈ Guo et al. 2017 (PMID 28973866) V^−2 in one condition,
    explicitly not universal. Shrinkage (v<1) with m>0 stiffens (E/E₀>1).
    """
    if m < 0.0:
        raise ValueError(f"m must be >= 0, got {m}")
    v_arr = np.asarray(v, dtype=float)
    if np.any(v_arr <= 0.0):
        raise ValueError("v must be > 0")
    return v_arr ** (-float(m))


# ===========================================================================
# Uncertainty-band machinery (the point of the whole module).
# ===========================================================================
class Envelope(NamedTuple):
    """Min / median / max band across a parameter sweep, aligned to a sample axis.

    `curves` keeps every realisation (shape (n_combos, n_samples)) so callers can plot the
    per-realisation spread with the band overlaid (design + CLAUDE.md viz-integrity rule).
    `labels` names each realisation by its parameter combo.
    """
    x: np.ndarray            # sample axis (time, r, or v) the band is a function of
    lo: np.ndarray           # per-sample minimum across the sweep
    med: np.ndarray          # per-sample median
    hi: np.ndarray           # per-sample maximum
    curves: np.ndarray       # (n_combos, n_samples) every swept realisation
    labels: list[str]        # parameter-combo label per row of `curves`


def sweep(func: Callable[..., np.ndarray], x: np.ndarray,
          grids: Mapping[str, Iterable[float]]) -> Envelope:
    """Evaluate `func(x, **combo)` over the Cartesian product of `grids` → an `Envelope`.

    Args:
        func: closed form taking the sample axis `x` as first positional arg plus one keyword
            per grid key; returns an array shaped like `x`.
        x: sample axis (time / r / v).
        grids: {param_name: iterable_of_values}. The full Cartesian product is swept — this is
            the "report the full envelope, not the median" discipline made mechanical.

    Returns:
        `Envelope` with lo/med/hi = min/median/max across all combos at each sample.
    """
    x = np.asarray(x, dtype=float)
    keys = list(grids.keys())
    value_lists = [list(grids[k]) for k in keys]
    curves: list[np.ndarray] = []
    labels: list[str] = []
    for combo in product(*value_lists):
        kwargs = dict(zip(keys, combo))
        curves.append(np.asarray(func(x, **kwargs), dtype=float).reshape(x.shape))
        labels.append(", ".join(f"{k}={v:g}" for k, v in kwargs.items()))
    stack = np.vstack(curves) if curves else np.empty((0, x.size))
    return Envelope(
        x=x,
        lo=stack.min(axis=0), med=np.median(stack, axis=0), hi=stack.max(axis=0),
        curves=stack, labels=labels,
    )


def passive_shrink_envelope(r_values: ArrayLike = (1.15, 1.30, 0.85, 0.70),
                            b_grid: Iterable[float] = B_GRID) -> Envelope:
    """Initial passive equilibrium-volume band vs medium ratio r, swept over b.

    x-axis = r; band = PBvH v_eq over b ∈ b_grid. E.g. r=1.15 → v_eq band for the "weak hyper"
    condition. This is the "for assumed b=0.1–0.3, initial volume decrease is ~X–Y%" envelope.
    """
    r_arr = np.asarray(r_values, dtype=float)
    return sweep(lambda r, b: pbvh_v_eq(r, b), r_arr, {"b": b_grid})


def timecourse_envelope(t: ArrayLike, r: float, tau_w: float = 1.0,
                        b_grid: Iterable[float] = B_GRID,
                        recovery_grid: Iterable[float] = R_GRID,
                        tau_reg_ratio_grid: Iterable[float] = TAU_REG_OVER_TAU_W_GRID) -> Envelope:
    """Volume-vs-time band for an impermeant medium ratio r, swept over (b, R, τ_reg/τ_w).

    The band spans no-recovery to full-recovery scenarios — deliberately, so the reader sees the
    R=0/0.5/1 fan rather than a single "recovers"/"doesn't" claim.
    """
    def _curve(tt: np.ndarray, b: float, R: float, tau_reg_ratio: float) -> np.ndarray:
        v_eq = float(pbvh_v_eq(r, b))
        return v_active(tt, v_eq, tau_w=tau_w, recovery=R, tau_reg=tau_reg_ratio * tau_w)

    return sweep(_curve, np.asarray(t, dtype=float),
                 {"b": b_grid, "R": recovery_grid, "tau_reg_ratio": tau_reg_ratio_grid})


def permeant_envelope(t: ArrayLike, r: float, tau_w: float = 1.0,
                      b_grid: Iterable[float] = B_GRID,
                      tau_s_ratio_grid: Iterable[float] = TAU_S_OVER_TAU_W_GRID,
                      sigma0: float = 1.0) -> Envelope:
    """Permeant-solute volume-vs-time band, swept over (b, τ_s/τ_w).

    Small τ_s/τ_w → fast solute equilibration → shrink-then-re-swell toward v≈1 (glycerol);
    large τ_s/τ_w → sustained shrink (mannitol-like). The band brackets both regimes.
    """
    def _curve(tt: np.ndarray, b: float, tau_s_ratio: float) -> np.ndarray:
        return v_permeant(tt, r=r, b=b, tau_w=tau_w, sigma0=sigma0, tau_s=tau_s_ratio * tau_w)

    return sweep(_curve, np.asarray(t, dtype=float),
                 {"b": b_grid, "tau_s_ratio": tau_s_ratio_grid})


def shape_envelope(v_values: ArrayLike, alpha_grid: Iterable[float] = ALPHA_GRID
                   ) -> dict[str, Envelope]:
    """Footprint (A/A₀) and height (h/h₀) bands vs volume ratio v, swept over α.

    Returns {"area_ratio": Envelope, "height_ratio": Envelope}. Separates "volume ↓" from
    "flattened/distorted": at fixed v the α sweep sets how the loss splits between footprint
    and height.
    """
    v_arr = np.asarray(v_values, dtype=float)
    area = sweep(lambda vv, alpha: shape_factors(vv, alpha).area_ratio, v_arr, {"alpha": alpha_grid})
    height = sweep(lambda vv, alpha: shape_factors(vv, alpha).height_ratio, v_arr, {"alpha": alpha_grid})
    return {"area_ratio": area, "height_ratio": height}


def stiffness_envelope(v_values: ArrayLike, m_grid: Iterable[float] = M_GRID) -> Envelope:
    """Apparent-stiffness band E/E₀ vs volume ratio v, swept over m ∈ {0,1,2}.

    The m=0 floor is "no change"; the band's top (m=2) is the Guo-2017-like maximal sensitivity.
    Report as "for m=0–2, stiffness change ranges from none to ~N×", never "stiffness increases".
    """
    v_arr = np.asarray(v_values, dtype=float)
    return sweep(lambda vv, m: stiffness_ratio(vv, m), v_arr, {"m": m_grid})


# ===========================================================================
# Condition report (design §"Conditions" minimum set → structured envelopes).
# ===========================================================================
@dataclass
class ConditionEnvelope:
    """Full band package for one medium condition (design §Conditions row)."""
    name: str
    r: float
    v_eq_lo: float                       # PBvH v_eq band over b (initial equilibrium)
    v_eq_hi: float
    timecourse: Envelope = field(repr=False)   # v(t) band over (b, R, τ_reg/τ_w)
    stiffness: Envelope = field(repr=False)     # E/E₀ band over m at the v_eq band edges

    @property
    def initial_volume_change_pct(self) -> tuple[float, float]:
        """(min%, max%) volume change at passive equilibrium; negative = shrink."""
        return (100.0 * (self.v_eq_lo - 1.0), 100.0 * (self.v_eq_hi - 1.0))


def condition_envelope(name: str, r: float, t: ArrayLike | None = None,
                       tau_w: float = 1.0) -> ConditionEnvelope:
    """Assemble the PBvH + time-course + stiffness band package for one condition."""
    t_arr = np.linspace(0.0, 8.0 * tau_w, 200) if t is None else np.asarray(t, dtype=float)
    v_eqs = np.array([float(pbvh_v_eq(r, b)) for b in B_GRID])
    tc = timecourse_envelope(t_arr, r, tau_w=tau_w)
    stiff = stiffness_envelope(np.array([v_eqs.min(), v_eqs.max()]))
    return ConditionEnvelope(
        name=name, r=r,
        v_eq_lo=float(v_eqs.min()), v_eq_hi=float(v_eqs.max()),
        timecourse=tc, stiffness=stiff,
    )


def all_condition_envelopes(tau_w: float = 1.0) -> list[ConditionEnvelope]:
    """Band package for the full minimum condition set (design §Conditions)."""
    return [condition_envelope(name, r, tau_w=tau_w) for name, r in CONDITIONS.items()]


# ===========================================================================
# Sanity gate (run: python -m aleph.validation.oracles.osmotic_response).
# ===========================================================================
def _sanity() -> None:
    import math

    # --- PBvH boundaries + Boyle-van't Hoff linearity ---
    for b in B_GRID:
        assert abs(float(pbvh_v_eq(1.0, b)) - 1.0) < 1e-12, "r=1 → v_eq=1 ∀ b"
        assert abs(float(pbvh_v_eq(1e12, b)) - b) < 1e-6, "r→∞ → v_eq→b"
    # linear in 1/r with intercept b, slope (1-b)
    inv_r = np.array([1 / 0.7, 1 / 0.85, 1.0, 1 / 1.15, 1 / 1.30])
    v = pbvh_v_eq(1.0 / inv_r, 0.2)
    slope, intercept = np.polyfit(inv_r, v, 1)
    assert abs(slope - 0.8) < 1e-9 and abs(intercept - 0.2) < 1e-9, "Boyle-van't Hoff linearity"
    # design worked example: r=1.2 → v_eq ∈ [0.850, 0.883]
    lo = float(pbvh_v_eq(1.2, 0.1)); hi = float(pbvh_v_eq(1.2, 0.3))
    assert abs(lo - 0.85) < 1e-3 and abs(hi - 0.8833) < 1e-3, f"r=1.2 band {lo:.3f},{hi:.3f}"

    # --- time course boundaries ---
    v_eq = float(pbvh_v_eq(1.2, 0.2))
    assert abs(float(v_passive(0.0, v_eq)) - 1.0) < 1e-12, "t=0 → v=1"
    assert abs(float(v_passive(1e6, v_eq)) - v_eq) < 1e-9, "t→∞ passive → v_eq"
    assert abs(float(v_active(1e6, v_eq, recovery=1.0)) - 1.0) < 1e-9, "R=1, t→∞ → 1 (full recovery)"
    assert np.allclose(v_active(np.linspace(0, 5, 20), v_eq, recovery=0.0),
                       v_passive(np.linspace(0, 5, 20), v_eq)), "R=0 == passive"
    # partial recovery lands strictly between passive and reference
    v_end_half = float(v_active(1e6, v_eq, recovery=0.5))
    assert v_eq < v_end_half < 1.0, "0<R<1 → partial recovery"

    # --- permeant σ / r_eff + shrink-then-re-swell ---
    assert abs(float(sigma_of_t(0.0, 1.0, 10.0)) - 1.0) < 1e-12, "σ(0)=σ₀"
    assert float(sigma_of_t(1e6, 1.0, 10.0)) < 1e-9, "σ(∞)=0"
    assert abs(float(r_eff_of_t(0.0, 1.3)) - 1.3) < 1e-12, "r_eff(0)=r"
    assert abs(float(r_eff_of_t(1e6, 1.3)) - 1.0) < 1e-6, "r_eff(∞)=1"
    tt = np.linspace(0.0, 200.0, 400)
    vp = v_permeant(tt, r=1.3, b=0.2, tau_w=1.0, sigma0=1.0, tau_s=10.0)
    assert abs(vp[0] - 1.0) < 1e-6, "permeant v(0)=1"
    assert vp.min() < 0.98, "permeant shrinks"
    assert vp[-1] > vp.min() and vp[-1] > 0.98, "permeant re-swells toward ~1 (glycerol)"
    # mannitol limit (tau_s→∞) collapses onto impermeant active with R=0
    vp_mann = v_permeant(tt, r=1.3, b=0.2, tau_w=1.0, sigma0=1.0, tau_s=1e9)
    va = v_active(tt, float(pbvh_v_eq(1.3, 0.2)), tau_w=1.0)
    assert np.allclose(vp_mann, va, atol=2e-3), "τ_s→∞ permeant == impermeant passive"

    # --- shape: volume-ratio conservation + boundaries ---
    for alpha in ALPHA_GRID:
        sf = shape_factors(0.85, alpha)
        assert abs(float(sf.lam_x * sf.lam_y * sf.lam_z) - 0.85) < 1e-12, "λx·λy·λz = v"
    sf13 = shape_factors(0.85, 1.0 / 3.0)   # isotropic
    assert abs(float(sf13.height_ratio) - 0.85 ** (1 / 3)) < 1e-12, "isotropic height = v^(1/3)"
    assert abs(float(sf13.height_ratio) - 0.947) < 1e-3, "v=0.85 isotropic → height −5.3%"
    sf0 = shape_factors(0.85, 0.0)          # fixed footprint
    assert abs(float(sf0.area_ratio) - 1.0) < 1e-12, "α=0 → A/A₀=1 (fixed footprint)"
    assert abs(float(sf0.height_ratio) - 0.85) < 1e-12, "α=0 → h/h₀=v (height-only, −15%)"

    # --- stiffness ---
    assert np.allclose(stiffness_ratio(np.array([0.7, 0.85, 1.2]), 0), 1.0), "m=0 → no change"
    assert abs(float(stiffness_ratio(0.85, 2)) - 0.85 ** -2) < 1e-12, "m=2 → v^-2"
    assert float(stiffness_ratio(0.85, 2)) > 1.0, "shrink + m>0 → stiffer"

    # --- envelope machinery ---
    env = passive_shrink_envelope((1.15, 1.30))
    assert env.lo.shape == (2,) and np.all(env.lo <= env.hi), "band ordered"
    assert env.curves.shape[0] == len(B_GRID), "one curve per b"
    tenv = timecourse_envelope(np.linspace(0, 8, 50), 1.30)
    assert np.all(tenv.lo <= tenv.med) and np.all(tenv.med <= tenv.hi), "lo≤med≤hi"

    print("osmotic_response sanity gate: PASS "
          "(PBvH linearity+boundaries, two-timescale RVD/RVI, permeant re-swell, "
          "shape volume-conservation, stiffness sensitivity, envelope ordering)")


if __name__ == "__main__":
    _sanity()
