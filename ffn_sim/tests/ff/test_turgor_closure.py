"""Osmotic turgor closure ΔP(V) — analytic ground truth + the DERIVED (no-magic) bulk modulus.

The resting turgor and its stiffness are NOT free parameters: `gamma_floor.turgor_pressure` implements
the Guo-2017 entropic enclosed-volume closure

    Π_in(V) = N·kB·T/(V − Vmin),   Vmin = vmin_frac·V0,   N·kB·T = Π_in0·(V0 − Vmin)
    ΔP(V)   = Π_in0·(V0 − Vmin)/(V − Vmin) − (Π_in0 − dP0),     ΔP(V0) = dP0

from LIT-anchored inputs only (Π_in0 = c_osm·R·T at c_osm≈200 mM; vmin_frac≈0.30 Venkova/Adar;
dP0 = 40 Pa resting). The bulk modulus is then DERIVED, not chosen — directly supporting the
"no empirical magic number" hard rule (the magic `K_vol` is removed):

    K_vol = −V·dΔP/dV|_{V0} = Π_in0/(1 − vmin_frac)

This pins the closure against its analytic form and the derived K_vol against its identity. (Tested on
the compression side, where ΔP > 0 and the physical max(ΔP,0) self-relief floor is inactive.)
"""

import numpy as np
import pytest

from ffn_sim.ff.gamma_floor import (
    TURGOR_DP0,
    TURGOR_PI_IN0,
    VMIN_FRAC,
    CortexParams,
    build_crosslinked_cortex,
    turgor_pressure,
)


@pytest.fixture(scope="module")
def shell():
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=200, n_xl=200, n_myo=20,
                                  rng=np.random.default_rng(0))
    ctr = cx.net.pos.mean(axis=0)
    R0 = float(np.linalg.norm(cx.net.pos - ctr, axis=1).mean())
    cx.R0_mean = R0                                    # self-consistent rest reference
    return cx, ctr, R0


def _dP(cx, ctr, scale):
    """ΔP at the shell scaled to radius·`scale` (uniform radial breathing about the centroid)."""
    return turgor_pressure(cx, ctr + (cx.net.pos - ctr) * scale)[0]


def _closure(R0, scale):
    """The analytic Guo closure ΔP(V) (unfloored) at radius·`scale`."""
    V0 = (4.0 / 3.0) * np.pi * R0**3
    vmin = VMIN_FRAC * V0
    V = (4.0 / 3.0) * np.pi * (R0 * scale) ** 3
    return TURGOR_PI_IN0 * (V0 - vmin) / (V - vmin) - (TURGOR_PI_IN0 - TURGOR_DP0)


def test_resting_pressure_is_dP0(shell):
    """At the rest volume V0 the net turgor is exactly the resting dP0 (40 Pa) — the physiological
    baseline the cortex is pre-tensioned against (not a convenient zero)."""
    cx, ctr, _ = shell
    assert _dP(cx, ctr, 1.0) == pytest.approx(TURGOR_DP0, rel=1e-9)


def test_closure_matches_guo_formula(shell):
    """ΔP(V) equals the Guo entropic closure to machine precision (compression side, unfloored)."""
    cx, ctr, R0 = shell
    for s in (0.995, 0.99, 0.98):
        assert _dP(cx, ctr, s) == pytest.approx(_closure(R0, s), rel=1e-9)


def test_derived_bulk_modulus_no_magic_Kvol(shell):
    """K_vol = −V·dΔP/dV|_{V0} = Π_in0/(1−vmin_frac) — the bulk modulus is DERIVED from the osmotic
    inputs, not a magic constant. Compression-side finite-difference (avoids the self-relief floor)."""
    cx, ctr, R0 = shell
    V0 = (4.0 / 3.0) * np.pi * R0**3
    s = 1e-4
    Vm = (4.0 / 3.0) * np.pi * (R0 * (1 - s)) ** 3
    K_fd = -V0 * (_dP(cx, ctr, 1.0) - _dP(cx, ctr, 1 - s)) / (V0 - Vm)
    K_analytic = TURGOR_PI_IN0 / (1.0 - VMIN_FRAC)
    assert K_fd == pytest.approx(K_analytic, rel=2e-3)


def test_compression_rises_expansion_self_relieves(shell):
    """ΔP is the physiological counter-force: compression (V↓) steeply RAISES it; a small expansion
    DROPS it below dP0 and the floor self-relieves the resting excess to ΔP=0 (Guo pump-leak)."""
    cx, ctr, _ = shell
    assert _dP(cx, ctr, 0.99) > 10 * TURGOR_DP0       # compression → steep counter-pressure
    assert _dP(cx, ctr, 1.01) == pytest.approx(0.0, abs=1e-9)   # small expansion → self-relieved (floored)
