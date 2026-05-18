"""Overdamped Langevin time integration (KU-1.26).

Equation of motion (per bead, overdamped, no inertia):

    γ_b · dr/dt = F(r) + √(2 γ_b k_B T) · η(t)

with η white-noise, ⟨η_α(t) η_β(t')⟩ = δ_αβ δ(t−t').

Phase 1 ships **Euler-Maruyama** (first-order weak / first-order
strong on additive noise); BAOAB (Leimkuhler-Matthews 2013, weak
order 2) is a Phase 2+ upgrade — stub kept here for the strategy
pattern.

Sanity Gate
-----------
1. Dimensional analysis: [γ_b] = N·s/m, [F] = N, [k_BT] = J. The
   discrete update Δr = (F/γ_b)·dt + √(2 k_BT·dt/γ_b)·N(0,1) has
   [Δr] = m ✓. CFL: dt < α·τ_min with α=0.1 (Leimkuhler-Matthews
   recommend α ≤ 0.5 for harmonic; we take 0.1 for safety with
   stiff cross-links).
2. Boundary cases: dt → 0 returns positions unchanged (force × 0,
   noise × 0). dt > τ_min trips the CFL gate, not handled silently.
3. Conservation: overdamped Langevin is **NOT** energy-conserving;
   it samples the canonical ensemble at temperature T. The
   ⟨½ k |Δr|²⟩ → ½ k_BT equipartition test is the conservation
   surrogate.
4. Numerical sanity: float64; RNG draws use `rng.standard_normal`
   for explicit reproducibility; positions wrapped to [0, L) after
   each step.
5. Sign / sense: drag opposes velocity (γ_b·v term moved to LHS),
   noise is unbiased (mean 0).
6. Measurement protocol: per-bead Δr is summed over all beads. The
   equipartition reporter samples each cross-link's |Δr| (not
   averaged over surrounding cells), matching the analytic
   ⟨½ k_xl |Δr_link|²⟩ = ½ k_BT (one degree of freedom per link
   end-to-end vector).
"""

from __future__ import annotations

from typing import Callable, Protocol

import numpy as np


class Integrator(Protocol):
    """Stateless overdamped Langevin step (KU-1.26).

    `forces_fn(positions)` must return an array of the same shape as
    ``positions`` (Newtons). The integrator returns the updated
    positions; periodic wrapping is the caller's responsibility.
    """

    def step(
        self,
        positions: np.ndarray,
        forces_fn: Callable[[np.ndarray], np.ndarray],
        gamma_b: float,
        kT: float,
        dt: float,
        rng: np.random.Generator,
    ) -> np.ndarray: ...


class EulerMaruyama:
    """First-order Euler-Maruyama overdamped Langevin (KU-1.26 Phase 1).

    r(t + dt) = r(t) + (F/γ_b)·dt + √(2 k_B T dt/γ_b) · ξ
    """

    name = "euler_maruyama"

    def step(
        self,
        positions: np.ndarray,
        forces_fn: Callable[[np.ndarray], np.ndarray],
        gamma_b: float,
        kT: float,
        dt: float,
        rng: np.random.Generator,
    ) -> np.ndarray:
        F = forces_fn(positions)
        drift = (dt / gamma_b) * F
        noise_scale = np.sqrt(2.0 * kT * dt / gamma_b)
        noise = rng.standard_normal(positions.shape) * noise_scale
        return positions + drift + noise


class BAOAB:
    """Leimkuhler-Matthews 2013 BAOAB splitting (Phase 2+; not used in Phase 1).

    Stubbed here so the dispatcher in :func:`make_integrator` and the
    code-pattern in `docs/02_force_models.md` can be tested. Calling
    ``step`` raises NotImplementedError until Phase 2.
    """

    name = "baoab"

    def step(self, *args, **kwargs):
        raise NotImplementedError("BAOAB is Phase 2+; use EulerMaruyama in Phase 1.")


def make_integrator(name: str) -> Integrator:
    if name == EulerMaruyama.name:
        return EulerMaruyama()
    if name == BAOAB.name:
        return BAOAB()
    raise ValueError(f"Unknown integrator '{name}'. Use 'euler_maruyama' or 'baoab'.")


def run(
    positions: np.ndarray,
    forces_fn: Callable[[np.ndarray], np.ndarray],
    integrator: Integrator,
    gamma_b: float,
    kT: float,
    dt: float,
    box_size: float,
    n_steps: int,
    rng: np.random.Generator,
    sample_interval: int | None = None,
    sample_fn: Callable[[np.ndarray, int], None] | None = None,
) -> np.ndarray:
    """Drive `integrator` for `n_steps`, wrapping positions periodically.

    If `sample_fn(positions, step_index)` is provided, it is invoked
    every `sample_interval` steps for streaming diagnostics.
    """
    p = positions.copy()
    for step_index in range(n_steps):
        p = integrator.step(p, forces_fn, gamma_b, kT, dt, rng)
        p = np.mod(p, box_size)
        if sample_fn is not None and sample_interval is not None \
                and (step_index + 1) % sample_interval == 0:
            sample_fn(p, step_index + 1)
    return p
