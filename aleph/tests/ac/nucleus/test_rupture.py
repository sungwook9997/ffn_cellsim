"""Self-test of the EMERGENT envelope-rupture oracle (pure NumPy — no Warp/CUDA).

Rupture is a FALSIFIABLE gate: it appears past a sourced envelope-strain threshold and is ABSENT below
it. The threshold is GAP (I0-B2) — the ON/OFF STRUCTURE is gated here, never the magnitude (do NOT tune).
"""
from __future__ import annotations

import numpy as np

from aleph.components.nucleus.lamina_analytic import LaminaParams, is_ruptured, lamina_tension

P = LaminaParams(k_chrom=2.0, k_lamin_b=1.0, k_lamin_ac=9.0, knee_strain=0.10, eps_rupture=0.50,
                 kappa_ne=0.0828)


def test_rupture_absent_below_threshold() -> None:
    """No face ruptures while strain < threshold; tension stays finite (intact envelope)."""
    eps = np.linspace(0.0, P.eps_rupture - 1e-3, 100)
    assert not np.any(is_ruptured(eps, P.eps_rupture))
    assert np.all(lamina_tension(eps, P) > 0.0 - 1e-12)


def test_rupture_emerges_above_threshold() -> None:
    """Faces strained past the threshold rupture, and their tension drops to exactly 0 (the tear)."""
    eps = np.array([P.eps_rupture + 1e-3, 0.7, 1.2])
    assert np.all(is_ruptured(eps, P.eps_rupture))
    assert np.allclose(lamina_tension(eps, P), 0.0)


def test_rupture_is_a_step_at_the_threshold() -> None:
    """The rupture mask flips exactly at eps_rupture (falsifiable on/off, not a gradual fade)."""
    assert not bool(is_ruptured(P.eps_rupture - 1e-9, P.eps_rupture))
    assert bool(is_ruptured(P.eps_rupture + 1e-9, P.eps_rupture))


def test_threshold_is_not_tuned_to_pass() -> None:
    """Sanity: the gate is threshold-parameterized, so a DIFFERENT sourced threshold still gives a
    clean on/off (the gate does not bake in a specific magnitude — report-not-tune)."""
    for thr in (0.3, 0.5, 0.8):
        assert not bool(is_ruptured(thr - 1e-6, thr))
        assert bool(is_ruptured(thr + 1e-6, thr))
