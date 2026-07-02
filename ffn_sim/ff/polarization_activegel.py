"""FF active movement piece-2/5 — cell polarization by actomyosin contractile-flow instability.

PI /goal 2026-07-02: symmetry-breaking = the front-rear axis a cell must pick before it can migrate directionally.
Mechanism A (chosen in FF_POLARIZATION_LITERATURE_2026-07-02 over Rho wave-pinning, per the mechanistic-not-lumped
rule): a uniform contractile actomyosin cortex is linearly UNSTABLE. Myosin generates active stress → stress
gradients drive cortical flow v → flow ADVECTS myosin → local density ↑ → active stress ↑ (positive feedback);
myosin diffusion + turnover oppose it. Above a contractility threshold ζ_c a single high-myosin cap (= the rear)
forms with steady cortical flow → polarity. (Bois, Jülicher & Grill 2011, PRL 106:028103; Mayer et al. 2010,
Nature 467:617.)

This is the active-gel reduction (1-D periodic cortex ring, arclength x∈[0,L=2πR)) — the mechanistic
hydrodynamic theory of the same explicit myosin the engine already carries. The coupled fields:

    myosin density :  ∂_t c = -∂_x(c v) + D ∂²_x c - k_off (c - c0)        (advection + diffusion + turnover)
    force balance  :  γ v - η ∂²_x v = ∂_x σ_a ,   σ_a = ζ f(c) ,  f(c)=c/(1+c/c*)   ;  ℓ ≡ √(η/γ)

solved spectrally (periodic). Validated against its own ANALYTIC linear dispersion (oracle-is-crosscheck rule):

    λ(k) = [c0 ζ f'(c0)/γ] · k²/(1 + ℓ²k²) - D k² - k_off .

Constants (µm, pN, s — FF): D, k_off from Bois 2011 (order-of-magnitude estimates); the hydrodynamic length
ℓ≈14 µm and cortical-flow-speed regime are the MEASURED anchors (Mayer 2010, C. elegans — ⚠️ no MCF7 datum,
PI-gated). ζ (contractility) is the SWEPT control variable (like the γ-floor myosin sweep — NOT tuned to an
outcome). γ is a normalization (absolute cortical friction is unmeasured — Mayer gives only η/γ=ℓ²); so absolute
flow SPEED is ζ/γ-scaled and not claimed against Mayer — only the length scale ℓ and the instability onset are.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Lit-anchored (DOIs in FF_POLARIZATION_LITERATURE_2026-07-02; flagged for PI KB-registration)
D_MYO_UM2_S = 1.0            # myosin diffusion (Bois 2011 estimate)
K_OFF_S = 0.1               # myosin turnover (Bois 2011 estimate)
ELL_HYDRO_UM = 14.0         # hydrodynamic length √(η/γ) — Mayer 2010 MEASURED (C. elegans; PI-gated for MCF7)


@dataclass(frozen=True, slots=True)
class ActiveGelParams:
    """1-D active-gel cortex-ring parameters (µm, pN, s)."""

    L: float                # domain length = 2πR [µm]
    D: float                # myosin diffusion [µm²/s]
    k_off: float            # turnover rate [s⁻¹]
    ell: float              # hydrodynamic length √(η/γ) [µm]
    gamma: float            # friction (normalization) [pN·s/µm³ · µm = pN·s/µm²]
    zeta: float             # contractility / active-stress amplitude [pN] — SWEPT control variable
    c0: float               # baseline myosin density [normalized]
    cstar: float            # active-stress saturation density [normalized]

    @property
    def eta(self) -> float:
        """Cortical viscosity η = γ·ℓ² (from the definition ℓ=√(η/γ))."""
        return self.gamma * self.ell ** 2


def resolve_activegel(*, R_um: float = 7.5, zeta: float = 1.0, gamma: float = 1.0,
                      c0: float = 1.0, cstar: float = 1.0, ell_um: float = ELL_HYDRO_UM,
                      D: float = D_MYO_UM2_S, k_off: float = K_OFF_S) -> ActiveGelParams:
    """Build the active-gel params for a cell of radius ``R_um`` (ring L=2πR). ``zeta`` is swept."""
    if R_um <= 0 or zeta < 0 or gamma <= 0 or cstar <= 0:
        raise ValueError("R_um,gamma,cstar>0 and zeta>=0 required")
    return ActiveGelParams(L=2.0 * np.pi * R_um, D=D, k_off=k_off, ell=ell_um,
                           gamma=gamma, zeta=zeta, c0=c0, cstar=cstar)


def f_active(c, cstar):
    """Saturating active-stress response σ_a/ζ = c/(1+c/c*) (monotone, saturates at high density)."""
    return c / (1.0 + c / cstar)


def f_active_prime(c, cstar):
    """d f/dc = 1/(1+c/c*)²."""
    return 1.0 / (1.0 + c / cstar) ** 2


def dispersion(k, p: ActiveGelParams):
    """Analytic linear growth rate λ(k) about the homogeneous state (the always-on validation ground truth)."""
    fp = f_active_prime(p.c0, p.cstar)
    return (p.c0 * p.zeta * fp / p.gamma) * k ** 2 / (1.0 + p.ell ** 2 * k ** 2) - p.D * k ** 2 - p.k_off


def zeta_critical(p: ActiveGelParams, k):
    """Contractility ζ at which λ(k)=0 for wavenumber k (marginal stability)."""
    fp = f_active_prime(p.c0, p.cstar)
    return (p.D * k ** 2 + p.k_off) * p.gamma * (1.0 + p.ell ** 2 * k ** 2) / (p.c0 * fp * k ** 2)


def _velocity(c, p: ActiveGelParams, k_grid):
    """Cortical flow v from overdamped force balance (γ - η∂²)v = ∂_x σ_a, solved spectrally (periodic)."""
    sa = p.zeta * f_active(c, p.cstar)
    sa_hat = np.fft.rfft(sa)
    v_hat = (1j * k_grid * sa_hat) / (p.gamma + p.eta * k_grid ** 2)
    return np.fft.irfft(v_hat, n=c.shape[0])


def step(c, p: ActiveGelParams, dx, dt, k_grid):
    """One explicit time step of ∂_t c = -∂_x(c v) + D ∂²_x c - k_off(c-c0) (spectral space derivatives)."""
    v = _velocity(c, p, k_grid)
    cv_hat = np.fft.rfft(c * v)
    c_hat = np.fft.rfft(c)
    dflux = np.fft.irfft(1j * k_grid * cv_hat, n=c.shape[0])          # ∂_x(c v)
    diff = np.fft.irfft(-(k_grid ** 2) * c_hat, n=c.shape[0])         # ∂²_x c
    dc = -dflux + p.D * diff - p.k_off * (c - p.c0)
    return c + dt * dc, v


def make_grid(p: ActiveGelParams, N=256):
    """Return (x, dx, k_grid) for the periodic ring of length L with N points (rfft wavenumbers)."""
    x = np.linspace(0.0, p.L, N, endpoint=False)
    dx = p.L / N
    k_grid = 2.0 * np.pi * np.fft.rfftfreq(N, d=dx)
    return x, dx, k_grid


def measure_growth_rate(p: ActiveGelParams, mode: int, *, N=256, eps=1e-4, n_steps=200, dt=None):
    """Seed c = c0 + eps·cos(mode·2πx/L), evolve, measure the numerical growth rate λ_num of that Fourier mode.

    Returns (lambda_numeric, lambda_analytic) for cross-check against :func:`dispersion`."""
    x, dx, k_grid = make_grid(p, N)
    k = mode * 2.0 * np.pi / p.L
    if dt is None:
        # CFL-ish: stable for diffusion + the active term; conservative
        dt = 0.05 * min(dx ** 2 / max(p.D, 1e-9), 1.0 / max(p.k_off, 1e-9))
    c = p.c0 + eps * np.cos(mode * 2.0 * np.pi * x / p.L)
    amps, ts = [], []
    for s in range(n_steps):
        ch = np.fft.rfft(c - p.c0)
        amps.append(np.abs(ch[mode])); ts.append(s * dt)
        c, _ = step(c, p, dx, dt, k_grid)
    amps = np.array(amps); ts = np.array(ts)
    good = amps > 0
    lam_num = np.polyfit(ts[good], np.log(amps[good]), 1)[0]
    return lam_num, dispersion(k, p)
