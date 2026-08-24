"""Self-test of the framework-#6 lamin-split areal-tension oracle (pure NumPy — no Warp/CUDA).

STRUCTURAL gates (magnitude-independent — they pass before the I0-B2 moduli close): the tension is
continuous at the knee (C0), the tangent modulus jumps UP by K_laminAC/K_soft (strain-STIFFENING, not
softening), tension is monotone up to rupture, and Young-Laplace ΔP=2σ/R closes the pressurised shell.
"""
from __future__ import annotations

import numpy as np
import pytest

from aleph.components.nucleus.lamina_analytic import (
    LaminaParams,
    lamina_tangent_modulus,
    lamina_tension,
    young_laplace_pressure,
)

# structural test params (NOT the sourced physiological values — those live in params_i0b2.yaml)
P = LaminaParams(k_chrom=2.0, k_lamin_b=1.0, k_lamin_ac=9.0, knee_strain=0.10, eps_rupture=0.50,
                 kappa_ne=0.0828)


def test_tension_continuous_at_knee() -> None:
    """σ is C0-continuous across the lamin-A/C engagement knee (no jump in tension, only in slope)."""
    lo = float(lamina_tension(P.knee_strain - 1e-9, P))
    hi = float(lamina_tension(P.knee_strain + 1e-9, P))
    assert lo == pytest.approx(hi, abs=1e-6)
    assert float(lamina_tension(0.0, P)) == pytest.approx(0.0, abs=1e-12)


def test_tangent_modulus_jumps_up_at_knee() -> None:
    """The tangent areal modulus is K_soft below the knee and K_laminAC above — a jump UP
    (strain-stiffening). The ratio is a REPORTED finding (K_laminAC/K_soft), not tuned."""
    below = float(lamina_tangent_modulus(0.05, P))
    above = float(lamina_tangent_modulus(0.20, P))
    assert below == pytest.approx(P.k_soft)
    assert above == pytest.approx(P.k_lamin_ac)
    assert above > below                                        # stiffening, not softening
    assert P.stiffening_ratio == pytest.approx(P.k_lamin_ac / P.k_soft)


def test_tension_monotone_increasing_below_rupture() -> None:
    """σ(ε) strictly increases with strain up to the rupture threshold."""
    e = np.linspace(0.0, P.eps_rupture - 1e-6, 200)
    s = lamina_tension(e, P)
    assert np.all(np.diff(s) > 0.0)


def test_young_laplace() -> None:
    """ΔP = 2σ/R for a spherical pressurised shell (the nucleoplasm-volume balance target)."""
    assert young_laplace_pressure(5.0, 3.0) == pytest.approx(2 * 5.0 / 3.0)
    assert young_laplace_pressure(1.0, 1.0) == pytest.approx(2.0)
    with pytest.raises(ValueError):
        young_laplace_pressure(1.0, 0.0)


def test_soft_regime_is_chromatin_plus_laminb() -> None:
    """The sub-knee tangent = chromatin + lamin-B (framework #6: both soft networks share small strain)."""
    assert P.k_soft == pytest.approx(P.k_chrom + P.k_lamin_b)
