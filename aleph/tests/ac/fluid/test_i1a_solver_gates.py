"""I1a discrete-solver gates — the conservative FV reference verified against the closed-form oracles.

Pure NumPy (no Warp/CUDA): these gate the DISCRETE stencil that the Warp-CUDA ``biot_substrate.py``
kernel ports to the device — the conservation identities a closed form cannot express, plus refinement
convergence back to Terzaghi/Green. The seven gates named in ``ac/fluid/__init__``'s Sanity Gate:

  1. constant-state preservation (static and cut/masked domain);
  2. impermeable-limit mass conservation to round-off;
  3. global content balance == integrated boundary flux + interior source;
  4. pressure-work sign (no-flux source-free energy is non-increasing; discrete Laplacian is SPD);
  5. CFL (stable at the registered limit, unstable above it);
  6. convergence to Terzaghi U(T_v) at ~2nd order in dx;
  7. convergence to the Green kernel spread rate d<r^2>/dt = 2 d c_v.

Uses the I0-B1 ledger closure c_v = mobility / S = 50 um^2/s.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.fluid.consolidation_analytic import (
    terzaghi_degree_of_consolidation,
    terzaghi_excess_pressure,
    time_factor,
)
from aleph.components.fluid.fv_reference import (
    FLUID,
    NUCLEUS,
    OUTSIDE,
    BiotFVReference,
    boundary_flux_integral,
    darcy_divergence,
)

# I0-B1 ledger closure (um-pN-s-Pa).
C_V = 50.0          # um^2/s
MOBILITY = 5.0e-3   # k/mu, um^2/(Pa*s)
STORAGE_S = 1.0e-4  # 1/M, 1/Pa


def test_ledger_closure_is_consistent() -> None:
    """c_v = mobility / S must hold (shared closure with the analytic oracles)."""
    assert C_V == pytest.approx(MOBILITY / STORAGE_S, rel=1e-12)


# ---------------------------------------------------------------------------------------------------
# Gate 1 — constant-state preservation
# ---------------------------------------------------------------------------------------------------
def test_constant_state_preserved_uniform_box() -> None:
    """A uniform field with zero source and no solid motion does not drift (all-fluid box)."""
    ref = BiotFVReference((16, 16, 16), dx=0.5, mobility=MOBILITY, storage_S=STORAGE_S)
    p = np.full(ref.shape, 37.0)  # Pi_0-like resting mean
    for _ in range(50):
        p = ref.step(p, ref.cfl_dt())
    assert np.allclose(p, 37.0, atol=1e-12)


def test_constant_state_preserved_cut_domain() -> None:
    """Uniform state is preserved even with a nucleus inclusion + outside cells (cut-cell mask)."""
    n = 21
    mask = np.full((n, n, n), FLUID, dtype=np.int_)
    xs = (np.arange(n) - n // 2) * 0.5
    xx, yy, zz = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(xx**2 + yy**2 + zz**2)
    mask[r > 4.5] = OUTSIDE      # membrane cut
    mask[r < 1.5] = NUCLEUS      # nuclear inclusion
    ref = BiotFVReference((n, n, n), dx=0.5, mobility=MOBILITY, storage_S=STORAGE_S, mask=mask)
    p = np.full(ref.shape, 40.0)
    p0 = p.copy()
    for _ in range(30):
        p = ref.step(p, ref.cfl_dt())  # membrane_flux defaults to 0 (no drive)
    assert np.allclose(p[mask == FLUID], p0[mask == FLUID], atol=1e-12)


# ---------------------------------------------------------------------------------------------------
# Gate 2 — impermeable-limit mass conservation
# ---------------------------------------------------------------------------------------------------
def test_impermeable_limit_conserves_mass() -> None:
    """No-flux box, zero source: total fluid content sum_i S p_i V_i is conserved to round-off."""
    rng = np.random.default_rng(1)
    ref = BiotFVReference((24, 24, 24), dx=0.4, mobility=MOBILITY, storage_S=STORAGE_S)
    p = 40.0 + 5.0 * rng.standard_normal(ref.shape)  # arbitrary non-uniform initial field
    c0 = ref.total_content(p)
    for _ in range(400):
        p = ref.step(p, ref.cfl_dt())
    assert ref.total_content(p) == pytest.approx(c0, rel=1e-11)


def test_impermeable_field_relaxes_to_the_mean() -> None:
    """The no-flux box relaxes toward the (conserved) spatial mean, not some other constant."""
    rng = np.random.default_rng(2)
    ref = BiotFVReference((20, 20, 20), dx=0.4, mobility=MOBILITY, storage_S=STORAGE_S)
    p = 40.0 + 8.0 * rng.standard_normal(ref.shape)
    mean0 = p.mean()
    for _ in range(6000):
        p = ref.step(p, ref.cfl_dt())
    assert np.allclose(p, mean0, atol=1e-3)


# ---------------------------------------------------------------------------------------------------
# Gate 3 — global content balance == integrated boundary flux + interior source
# ---------------------------------------------------------------------------------------------------
def test_content_balance_equals_membrane_flux() -> None:
    """d/dt(total content) == integrated membrane boundary flux (divergence theorem), to round-off.

    A spherical fluid blob (OUTSIDE beyond it) is driven by a uniform outward membrane flux; the
    content change rate must equal the integrated boundary flux exactly (discrete Gauss).
    """
    n = 25
    xs = (np.arange(n) - n // 2) * 0.5
    xx, yy, zz = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(xx**2 + yy**2 + zz**2)
    mask = np.where(r <= 5.0, FLUID, OUTSIDE).astype(np.int_)
    ref = BiotFVReference((n, n, n), dx=0.5, mobility=MOBILITY, storage_S=STORAGE_S, mask=mask)
    p = np.full(ref.shape, 40.0)
    memb = -3.0e-4  # inward flux (negative outward) -> content should RISE

    dcontent_dt = float(
        STORAGE_S
        * np.sum(ref.rhs(p, membrane_flux=memb)[mask == FLUID])
        * ref.dx**ref.dim
    )
    integrated_flux = boundary_flux_integral(p, MOBILITY, ref.dx, mask, membrane_flux=memb)
    # d(content)/dt = -(outward flux)   (outflux drains content)
    assert dcontent_dt == pytest.approx(-integrated_flux, rel=1e-10, abs=1e-14)


def test_interior_source_changes_content_at_its_integral() -> None:
    """With a no-flux box and an interior source s_water, d(content)/dt == integral of s_water."""
    ref = BiotFVReference((18, 18, 18), dx=0.5, mobility=MOBILITY, storage_S=STORAGE_S)
    rng = np.random.default_rng(3)
    p = 40.0 + rng.standard_normal(ref.shape)
    s = 1.0e-3 * rng.standard_normal(ref.shape)
    dcontent_dt = float(STORAGE_S * np.sum(ref.rhs(p, s_water=s)) * ref.dx**ref.dim)
    assert dcontent_dt == pytest.approx(float(np.sum(s) * ref.dx**ref.dim), rel=1e-10)


def test_solid_dilatation_is_a_content_source() -> None:
    """alpha div(v_s) enters as a source: uniform compression div(v_s)<0 raises fluid content."""
    ref = BiotFVReference((16, 16, 16), dx=0.5, mobility=MOBILITY, storage_S=STORAGE_S)
    p = np.full(ref.shape, 40.0)
    dil = np.full(ref.shape, -2.0e-3)  # skeleton compressing -> squeezes fluid in
    rate = ref.rhs(p, solid_dilatation_rate=dil, alpha=1.0)
    assert np.all(rate > 0.0)  # p rises everywhere


# ---------------------------------------------------------------------------------------------------
# Gate 4 — pressure-work sign (energy non-increasing; discrete Laplacian SPD)
# ---------------------------------------------------------------------------------------------------
def test_diffusion_energy_monotone_nonincreasing() -> None:
    """No-flux, source-free: the Lyapunov energy 0.5 S sum p^2 V never increases (Darcy dissipation)."""
    rng = np.random.default_rng(4)
    ref = BiotFVReference((22, 22, 22), dx=0.4, mobility=MOBILITY, storage_S=STORAGE_S)
    p = 40.0 + 6.0 * rng.standard_normal(ref.shape)
    e_prev = ref.diffusion_energy(p)
    for _ in range(300):
        p = ref.step(p, ref.cfl_dt())
        e = ref.diffusion_energy(p)
        assert e <= e_prev + 1e-12
        e_prev = e


def test_discrete_laplacian_is_symmetric_negative_semidefinite() -> None:
    """Assemble -div(q)/mobility as a matrix on a small no-flux grid: symmetric, eigenvalues <= 0."""
    n = 6
    ref = BiotFVReference((n, n), dx=0.5, mobility=MOBILITY, storage_S=STORAGE_S)
    size = n * n
    lap = np.zeros((size, size))
    for j in range(size):
        e = np.zeros((n, n))
        e.flat[j] = 1.0
        # -div(q) = mobility * laplacian; divide out mobility to get the pure operator
        lap[:, j] = (-darcy_divergence(e, MOBILITY, ref.dx) / MOBILITY).ravel()
    assert np.allclose(lap, lap.T, atol=1e-12)              # symmetric
    evals = np.linalg.eigvalsh(lap)
    assert np.all(evals <= 1e-9)                            # negative semidefinite
    assert abs(evals[-1]) < 1e-9                            # a constant is the zero mode (no-flux)


# ---------------------------------------------------------------------------------------------------
# Gate 5 — CFL
# ---------------------------------------------------------------------------------------------------
def test_cfl_stable_below_unstable_above() -> None:
    """Stable at the registered CFL limit; forward Euler blows up above it (a genuine limit)."""
    rng = np.random.default_rng(5)
    ref = BiotFVReference((16, 16, 16), dx=0.4, mobility=MOBILITY, storage_S=STORAGE_S)
    p0 = 40.0 + rng.standard_normal(ref.shape)

    p = p0.copy()
    for _ in range(200):
        p = ref.step(p, ref.cfl_dt(safety=0.9))
    assert np.isfinite(p).all() and p.std() < p0.std()  # decays, stays finite

    dt_unstable = ref.cfl_dt(safety=1.0) / 0.9 * 1.3     # 1.3x the true limit
    p = p0.copy()
    for _ in range(200):
        p = ref.step(p, dt_unstable)
    assert not np.isfinite(p).all() or p.std() > 10.0 * p0.std()  # amplifies


# ---------------------------------------------------------------------------------------------------
# Gate 6 — convergence to Terzaghi
# ---------------------------------------------------------------------------------------------------
def _run_terzaghi(n_cells: int, t_target: float, u0: float = 1.0) -> tuple[np.ndarray, float]:
    """Explicit 1-D doubly-drained consolidation; return (final pressure profile, T_v reached)."""
    half_h = 5.0                       # drainage path H (um)
    length = 2.0 * half_h
    dx = length / n_cells
    ref = BiotFVReference((n_cells,), dx=dx, mobility=MOBILITY, storage_S=STORAGE_S,
                          dirichlet={(0, 0): 0.0, (0, 1): 0.0})
    p = np.full(n_cells, u0)
    dt = ref.cfl_dt(safety=0.5)
    n_steps = int(np.ceil(t_target / dt))
    dt = t_target / n_steps            # land exactly on t_target
    for _ in range(n_steps):
        p = ref.step(p, dt)
    tv = float(time_factor(C_V, t_target, half_h))
    return p, tv


def test_terzaghi_degree_of_consolidation_matches_oracle() -> None:
    """Mean-pressure degree of consolidation U reproduces the Terzaghi series."""
    t_target = 0.5  # s  -> T_v = c_v t / H^2 = 50*0.5/25 = 1.0
    p, tv = _run_terzaghi(160, t_target)
    u_numeric = 1.0 - p.mean() / 1.0
    u_oracle = float(terzaghi_degree_of_consolidation(np.array([tv]))[0])
    assert u_numeric == pytest.approx(u_oracle, abs=3e-3)


def test_terzaghi_profile_second_order_convergence() -> None:
    """The pressure profile converges to the Terzaghi series at ~2nd order in dx."""
    t_target = 0.15  # s -> T_v = 0.3 (mid-dissipation, series well-resolved)
    half_h = 5.0
    errors = []
    grids = (40, 80, 160)
    for n in grids:
        p, tv = _run_terzaghi(n, t_target)
        z_over_h = (np.arange(n) + 0.5) * (2.0 * half_h / n) / half_h   # cell centres in z/H
        p_oracle = terzaghi_excess_pressure(z_over_h, tv, u0=1.0)
        errors.append(float(np.sqrt(np.mean((p - p_oracle) ** 2))))
    # error should drop ~4x per halving of dx (2nd order)
    order_lo = np.log2(errors[0] / errors[1])
    order_hi = np.log2(errors[1] / errors[2])
    assert order_lo > 1.7 and order_hi > 1.7
    assert errors[-1] < 2e-3


# ---------------------------------------------------------------------------------------------------
# Gate 7 — convergence to the Green kernel spread
# ---------------------------------------------------------------------------------------------------
def test_greens_spread_rate_matches_2dcv() -> None:
    """A compact packet in a no-flux box spreads at d<r^2>/dt = 2 d c_v (diffusion law), reading c_v."""
    n = 61
    dx = 0.6
    ref = BiotFVReference((n, n, n), dx=dx, mobility=MOBILITY, storage_S=STORAGE_S)
    xs = (np.arange(n) - n // 2) * dx
    xx, yy, zz = np.meshgrid(xs, xs, xs, indexing="ij")
    r2 = xx**2 + yy**2 + zz**2
    sig0 = 2.0 * dx
    p = np.exp(-r2 / (2.0 * sig0**2))  # compact Gaussian packet at the centre

    def mean_sq_radius(field: np.ndarray) -> float:
        w = np.clip(field, 0.0, None)
        return float(np.sum(r2 * w) / np.sum(w))

    dt = ref.cfl_dt(safety=0.5)
    n_steps = 60
    r2_0 = mean_sq_radius(p)
    for _ in range(n_steps):
        p = ref.step(p, dt)
    r2_1 = mean_sq_radius(p)
    rate = (r2_1 - r2_0) / (n_steps * dt)
    assert rate == pytest.approx(2.0 * 3.0 * C_V, rel=0.05)  # d=3
