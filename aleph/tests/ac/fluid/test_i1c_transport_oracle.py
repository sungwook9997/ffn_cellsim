"""I1c conservative RAD transport gates — total-actin conservation, FRAP, advection-follows-v_f.

Pure NumPy (no Warp/CUDA). Gates the monomer economy d(phi c)/dt + div(phi c v_f - phi D_c grad c) = R:
  * total actin integral(phi c) + N_polymer conserved to MACHINE PRECISION (the hard non-vacuous invariant);
  * positivity under the CFL;
  * FRAP: a cosine bleach recovers at D_c|k|^2 (reads D_c);
  * advection follows the absolute pore-fluid velocity v_f, NOT the discharge q;
  * uniform-state preserved (no spurious interior source);
  * treadmill balance: barbed consumption == pointed release => stationary free pool.

D_c is the I0-B1c draft ~3 um^2/s (crowded cytoplasm); phi is the I0-B1b GAP (swept as an oracle param).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.fluid.transport_analytic import frap_cosine_recovery_rate, treadmill_net_source
from aleph.components.fluid.transport_reference import RADTransportReference, rad_cfl_dt

D_C = 3.0    # um^2/s (I0-B1c draft, crowded cytoplasm)
PHI = 0.7    # I0-B1b GAP — swept; used here only to relate u=phi*c


def test_total_actin_conserved_to_machine_precision() -> None:
    """Closed domain, arbitrary v_f + reaction: integral(phi c) + bound is invariant to round-off."""
    rng = np.random.default_rng(0)
    shape = (16, 16, 16)
    ref = RADTransportReference(shape, dx=0.4, phi=PHI, d_c=D_C)
    u = np.abs(rng.standard_normal(shape)) * PHI  # non-negative density phi*c
    v_f = 0.3 * rng.standard_normal((*shape, 3))  # um/s
    reaction = 1e-2 * rng.standard_normal(shape)  # exchange with the bound pool
    ref.bound_monomer = 5.0
    total0 = ref.total_actin(u)
    dt = rad_cfl_dt(ref.dx, D_C, float(np.max(np.abs(v_f))), ref.dim)
    for _ in range(200):
        u = ref.step(u, dt, v_f=v_f, reaction=reaction)
    assert ref.total_actin(u) == pytest.approx(total0, rel=1e-11)


def test_positivity_preserved() -> None:
    """A non-negative field stays non-negative under upwind advection + diffusion at the CFL."""
    rng = np.random.default_rng(1)
    shape = (20, 20)
    ref = RADTransportReference(shape, dx=0.5, phi=PHI, d_c=D_C)
    u = np.zeros(shape)
    u[8:12, 8:12] = 2.0  # a positive patch
    v_f = np.zeros((*shape, 2))
    v_f[..., 0] = 0.8
    dt = rad_cfl_dt(ref.dx, D_C, 0.8, ref.dim)
    for _ in range(150):
        u = ref.step(u, dt, v_f=v_f)
    assert np.all(u >= -1e-15)


def test_frap_cosine_recovers_at_dc_k_squared() -> None:
    """Pure diffusion (v_f=0, R=0): a cosine bleach pattern recovers at rate ~ D_c |k|^2 (reads D_c)."""
    n = 200
    L = 20.0
    dx = L / n
    ref = RADTransportReference((n,), dx=dx, phi=PHI, d_c=D_C)
    x = (np.arange(n) + 0.5) * dx
    k = np.pi / L  # zero-gradient (no-flux) at both walls
    u = PHI * (1.0 + 0.3 * np.cos(k * x))  # bleach = cosine perturbation on a uniform background
    amp0 = 0.3 * PHI
    dt = rad_cfl_dt(dx, D_C, 0.0, 1)
    t_total = 0.0
    for _ in range(400):
        u = ref.step(u, dt)
        t_total += dt
    # perturbation amplitude = (max-min)/2 of the deviation from the mean
    dev = u - u.mean()
    amp = 0.5 * (dev.max() - dev.min())
    rate_measured = -np.log(amp / amp0) / t_total
    assert rate_measured == pytest.approx(frap_cosine_recovery_rate([k], D_C), rel=0.03)


def test_advection_front_follows_v_f_not_q() -> None:
    """Pure advection (D_c=0): a blob centroid travels at v_f, distinguishable from q = phi*v_f."""
    n = 400
    dx = 0.25
    ref = RADTransportReference((n,), dx=dx, phi=PHI, d_c=0.0)
    x = (np.arange(n) + 0.5) * dx
    v_f = 1.0  # um/s
    u = np.exp(-((x - 20.0) ** 2) / (2.0 * (2.0 * dx) ** 2))  # blob near the inflow quarter
    vf = np.full((n, 1), v_f)
    dt = rad_cfl_dt(dx, 0.0, v_f, 1)
    n_steps = 200

    def centroid(f: np.ndarray) -> float:
        return float(np.sum(x * f) / np.sum(f))

    c0 = centroid(u)
    for _ in range(n_steps):
        u = ref.step(u, dt, v_f=vf)
    speed = (centroid(u) - c0) / (n_steps * dt)
    assert speed == pytest.approx(v_f, rel=1e-6)          # rides v_f exactly
    assert not np.isclose(speed, PHI * v_f, rtol=0.05)    # NOT the discharge q = phi*v_f


def test_uniform_state_has_no_spurious_interior_source() -> None:
    """Uniform u advected by a uniform v_f has zero divergence in the strict interior (no spurious source)."""
    shape = (16, 16, 16)
    ref = RADTransportReference(shape, dx=0.5, phi=PHI, d_c=D_C)
    u = np.full(shape, 0.55)
    v_f = np.zeros((*shape, 3))
    v_f[..., 0] = 0.4
    v_f[..., 1] = -0.2
    div = ref._flux_divergence(u, v_f)
    interior = div[2:-2, 2:-2, 2:-2]
    assert np.allclose(interior, 0.0, atol=1e-12)


def test_uniform_state_preserved_quiescent() -> None:
    """A uniform field is exactly preserved with no flow and no reaction."""
    shape = (12, 12, 12)
    ref = RADTransportReference(shape, dx=0.5, phi=PHI, d_c=D_C)
    u = np.full(shape, 0.7)
    for _ in range(50):
        u = ref.step(u, rad_cfl_dt(0.5, D_C, 0.0, 3))
    assert np.allclose(u, 0.7, atol=1e-13)


def test_treadmill_balance_stationary_pool() -> None:
    """Barbed consumption balanced by pointed release => integral R = 0 and the free pool is stationary."""
    assert treadmill_net_source(barbed_consumption=4.2, pointed_release=4.2) == pytest.approx(0.0)
    shape = (14, 14, 14)
    ref = RADTransportReference(shape, dx=0.5, phi=PHI, d_c=D_C)
    u = np.full(shape, 0.5)
    reaction = np.zeros(shape)
    reaction[3, 3, 3] = -0.1   # barbed consumption (sink)
    reaction[10, 10, 10] = 0.1  # pointed release (source) — equal and opposite
    content0 = ref.field_content(u)
    bound0 = ref.bound_monomer
    dt = rad_cfl_dt(0.5, D_C, 0.0, 3)
    for _ in range(80):
        u = ref.step(u, dt, reaction=reaction)
    assert ref.field_content(u) == pytest.approx(content0, rel=1e-10)   # total free pool unchanged
    assert ref.bound_monomer == pytest.approx(bound0, abs=1e-12)        # net exchange zero
