r"""Host gate for the assembled-cell cortical-tension harness — the actin-axial-tension extraction physics.

Validates :func:`actin_axial_tension` (the external-force running sum that recovers the actin inextensibility
constraint tension) and :func:`crosslink_tension`. The device orchestration is a runtime path; here we gate the
extraction MATH (pure NumPy, out-of-hot-loop).
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.incumbent.assembled_cortical_stress import actin_axial_tension, crosslink_tension


def test_uniform_tension_from_end_loads() -> None:
    """A filament pulled apart only at its two ends carries a UNIFORM axial tension = the end force."""
    pos = np.array([[0.0, 0, 0], [0.2, 0, 0], [0.4, 0, 0], [0.6, 0, 0]])   # straight, +x
    f_ext = np.array([[-1.0, 0, 0], [0, 0, 0], [0, 0, 0], [1.0, 0, 0]])    # ends pulled APART ⇒ tension
    off = np.array([0, 4])
    sa, sb, T = actin_axial_tension(f_ext, pos, off)
    assert list(sa) == [0, 1, 2] and list(sb) == [1, 2, 3]
    assert T == pytest.approx([1.0, 1.0, 1.0])                             # uniform +1 (tension) everywhere


def test_compression_sign() -> None:
    """Ends pushed together (contractile inward pull at the ends) → negative tension (compression)."""
    pos = np.array([[0.0, 0, 0], [0.3, 0, 0], [0.6, 0, 0]])
    f_ext = np.array([[1.0, 0, 0], [0, 0, 0], [-1.0, 0, 0]])               # ends pushed IN ⇒ compression
    _, _, T = actin_axial_tension(f_ext, pos, np.array([0, 3]))
    assert T == pytest.approx([-1.0, -1.0])


def test_distributed_load_running_sum() -> None:
    """An interior load relieves the tension beyond it (the running sum steps down)."""
    pos = np.array([[0.0, 0, 0], [0.2, 0, 0], [0.4, 0, 0], [0.6, 0, 0]])
    f_ext = np.array([[-2.0, 0, 0], [0, 0, 0], [1.0, 0, 0], [1.0, 0, 0]])  # Σ = 0 (equilibrium)
    _, _, T = actin_axial_tension(f_ext, pos, np.array([0, 4]))
    assert T == pytest.approx([2.0, 2.0, 1.0])                             # 2 → 2 → (relieved by +1 at node 2) → 1


def test_two_fibers_and_short_fiber_skipped() -> None:
    """Multiple fibers are handled independently; a <2-node fiber contributes no segment."""
    pos = np.array([[0, 0, 0], [0.2, 0, 0],           # fiber 0: 2 nodes → 1 segment
                    [0, 1, 0],                          # fiber 1: 1 node → skipped
                    [0, 2, 0], [0, 2.2, 0], [0, 2.4, 0]], float)   # fiber 2: 3 nodes → 2 segments
    f_ext = np.zeros((6, 3)); f_ext[0] = [-1, 0, 0]; f_ext[1] = [1, 0, 0]
    off = np.array([0, 2, 3, 6])
    sa, sb, T = actin_axial_tension(f_ext, pos, off)
    assert list(sa) == [0, 3, 4]                        # fiber1's lone node produced no segment
    assert T[0] == pytest.approx(1.0)                   # fiber0: ends pulled apart ⇒ +1 tension
    assert T[1:] == pytest.approx([0.0, 0.0])           # fiber2 unloaded


def test_transverse_load_does_not_add_axial_tension() -> None:
    """A force perpendicular to the filament adds bending, not axial tension (projection is zero)."""
    pos = np.array([[0.0, 0, 0], [0.2, 0, 0], [0.4, 0, 0]])
    f_ext = np.array([[0.0, 1.0, 0], [0, -2.0, 0], [0, 1.0, 0]])           # pure transverse (Σ=0)
    _, _, T = actin_axial_tension(f_ext, pos, np.array([0, 3]))
    assert T == pytest.approx([0.0, 0.0], abs=1e-12)


def test_crosslink_tension_hookean() -> None:
    """Crosslink tension is k(|r|−r0); stretched → +, compressed → −."""
    pos = np.array([[0, 0, 0], [0.1, 0, 0], [0, 0, 0], [0.02, 0, 0]], float)
    xl = np.array([[0, 1], [2, 3]])
    a, b, T = crosslink_tension(pos, xl, kxl=np.array([1000.0, 500.0]), r0xl=np.array([0.05, 0.05]))
    assert T == pytest.approx([1000.0 * (0.1 - 0.05), 500.0 * (0.02 - 0.05)])   # +50, −15
