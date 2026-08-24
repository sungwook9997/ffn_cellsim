"""Pure-NumPy reference gates for the fiber-arclength multigrid building blocks (2:1 R/P + O(L) line smoother).

GPU-only contract (CLAUDE.md I0-A, PI 2026-07-16): ``ac`` kernels are Warp-CUDA and are **never launched on the
CPU device**, not even in tests. So — following the sanctioned dev-box pattern used by the fluid-substrate /
cytosol tests — this module validates the MATH of ``ac/cell/fiber_arclength_mg.py`` with a **pure-NumPy
reference reimplementation** of the restriction ``R``, prolongation ``P`` and the tridiagonal line smoother,
and asserts the Warp **kernel source** carries the matching stencil for the one property (the device Thomas
finite-latch) that only the kernel can produce. The Warp-kernel launch is native-validated by the Lead on the
A5000 GPU (§4.4 ladder of ``aleph/docs/v2_audit/FINE_MESH_MULTIGRID_DESIGN_2026-07-24.md``).

``build_arclength_restriction`` is a pure host (NumPy) function, so its 2:1 coarsening table + coarse ``foff``
are exercised DIRECTLY. ``R``/``P``/the line smoother are matrix-free device kernels — each is mirrored here by
a faithful NumPy function read off the kernel body (``arclength_restrict_kernel`` = ``np.add.at`` full-weighting;
``arclength_prolong_add_kernel`` = weighted gather; ``fiber_line_smooth_kernel`` = a direct solve of the
row-sum-lumped SPD tridiagonal ``M_f`` the Thomas sweep factors).

Covered (task deliverable (a)-(d)):
    (a) linear-interpolation exactness of the NumPy ``P`` on arclength-linear data + ``R∘P`` full-weighting scale;
    (b) transpose consistency ``<R x, y>_c == <x, P y>_f`` to round-off (the Galerkin requirement ``R == P^T``);
    (c) the NumPy line smoother = a solve of the analytic row-sum-lumped tridiagonal ``M_f``; it damps every
        high-arclength-frequency mode of the exact pentadiagonal fiber operator by a solid factor and is a
        near no-op on the smooth (near-null) component;
    (d) odd-length + 2-node fibers coarsen + smooth without crashing (finite predicate stays set in NumPy; the
        DEVICE non-finite-pivot latch is asserted at the kernel-source level, native-validated by the Lead).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

# The module defines @wp.kernel bodies, so importing the host helpers requires warp importable — but NOTHING
# here launches a kernel or allocates a device array (GPU-only contract: no CPU-Warp launches in ac tests).
pytest.importorskip("warp")

import aleph.components.incumbent.fiber_arclength_mg as _mg  # noqa: E402
from aleph.components.incumbent.fiber_arclength_mg import build_arclength_restriction  # noqa: E402

_MG_SRC = Path(_mg.__file__).read_text(encoding="utf-8")


# ----------------------------------------------------------------------------------------------------
# Pure-NumPy references (faithful to the kernel bodies; native launch validated by the Lead).
# ----------------------------------------------------------------------------------------------------


def _restrict_np(level: dict, fine: np.ndarray) -> np.ndarray:
    """NumPy mirror of ``arclength_restrict_kernel``: ``coarse = R fine`` with ``R = P^T`` (full-weighting).

    The kernel scatters ``w0*fine`` to ``c0`` and ``w1*fine`` to ``c1`` with atomic adds; ``np.add.at`` is the
    exact serial equivalent (and reads the SAME ``(c, w)`` table as :func:`_prolong_np`, so ``R == P^T``)."""
    fine = np.ascontiguousarray(fine, dtype=np.float64)
    n_coarse = int(level["n_coarse"])
    coarse = np.zeros((max(n_coarse, 1), 3), dtype=np.float64)
    c0, w0, c1, w1 = level["c0"], level["w0"], level["c1"], level["w1"]
    m0 = c0 >= 0
    np.add.at(coarse, c0[m0], w0[m0, None] * fine[m0])
    m1 = c1 >= 0
    np.add.at(coarse, c1[m1], w1[m1, None] * fine[m1])
    return coarse[:n_coarse]


def _prolong_np(level: dict, coarse: np.ndarray, n_fine: int) -> np.ndarray:
    """NumPy mirror of ``arclength_prolong_add_kernel``: ``fine = P coarse`` (1-D linear arclength interp)."""
    coarse = np.ascontiguousarray(coarse, dtype=np.float64)
    fine = np.zeros((n_fine, 3), dtype=np.float64)
    c0, w0, c1, w1 = level["c0"], level["w0"], level["c1"], level["w1"]
    m0 = c0 >= 0
    fine[m0] += w0[m0, None] * coarse[c0[m0]]
    m1 = c1 >= 0
    fine[m1] += w1[m1, None] * coarse[c1[m1]]
    return fine


def _lumped_bending(length: int, alpha: np.ndarray) -> np.ndarray:
    """The row-sum-lumped tridiagonal bending band ``fiber_line_smooth_kernel`` assembles (one component, L x L).

    Per triple ``t`` on nodes ``(t, t+1, t+2)``: diag ``+2a, +4a, +2a`` and sub/super ``-2a, -2a`` (the two
    NF2007 pentadiagonal corners lumped onto the diagonal) — exactly the kernel's ``two_a``/``four_a`` updates."""
    diag = np.zeros(length, dtype=np.float64)
    off = np.zeros(max(length - 1, 0), dtype=np.float64)  # off[k] between nodes k and k+1
    for t in range(length - 2):
        a = float(alpha[t])
        diag[t] += 2.0 * a
        diag[t + 1] += 4.0 * a
        diag[t + 2] += 2.0 * a
        off[t] -= 2.0 * a
        off[t + 1] -= 2.0 * a
    matrix = np.diag(diag)
    for k in range(length - 1):
        matrix[k + 1, k] = off[k]
        matrix[k, k + 1] = off[k]
    return matrix


def _lumped_tridiagonal(ext_scalar: np.ndarray, alpha: np.ndarray) -> np.ndarray:
    """Full per-fiber ``M_f = D_external + lumped(D2^T diag(alpha) D2)`` for one Cartesian component (L x L)."""
    return _lumped_bending(ext_scalar.shape[0], alpha) + np.diag(ext_scalar.astype(np.float64))


def _pentadiagonal_operator(ext: np.ndarray, alpha: np.ndarray) -> np.ndarray:
    """Exact fiber sub-operator ``D_external + D2^T diag(alpha) D2`` (pentadiagonal, one component)."""
    length = ext.shape[0]
    d2 = np.zeros((length - 2, length), dtype=np.float64)
    for t in range(length - 2):
        d2[t, t] = 1.0
        d2[t, t + 1] = -2.0
        d2[t, t + 2] = 1.0
    return np.diag(ext.astype(np.float64)) + d2.T @ np.diag(alpha) @ d2


def _line_smooth_np(residual, x, foff, toff, alpha, ext, omega):
    """NumPy mirror of ``fiber_line_smooth_kernel``: per fiber ``x += omega * M_f^-1 residual`` (direct solve).

    The device kernel Thomas-solves the SPD tridiagonal ``M_f`` = ``ext ⊙ I + lumped bending``; a direct
    ``np.linalg.solve`` of the SAME matrix is the exact reference (three componentwise systems that share the
    scalar bending band but keep their own ``ext`` component). Returns ``(x, finite)``; the well-conditioned
    NumPy path never trips the pivot latch (that DEVICE guard is checked at the kernel-source level)."""
    x = np.array(x, dtype=np.float64, copy=True)
    residual = np.asarray(residual, dtype=np.float64)
    ext = np.asarray(ext, dtype=np.float64)
    for f in range(len(foff) - 1):
        begin, end = int(foff[f]), int(foff[f + 1])
        length = end - begin
        if length <= 0:
            continue
        tb = int(toff[f])
        alpha_f = np.asarray(alpha[tb:tb + max(length - 2, 0)], dtype=np.float64)
        bend = _lumped_bending(length, alpha_f)
        for c in range(3):
            matrix = bend + np.diag(ext[begin:end, c])
            correction = np.linalg.solve(matrix, residual[begin:end, c])
            x[begin:end, c] += omega * correction
    return x, 1


# ----------------------------------------------------------------------------------------------------
# (a) linear-interpolation exactness + R∘P full-weighting scale.
# ----------------------------------------------------------------------------------------------------


def test_prolong_is_linear_exact_on_arclength_linear_data() -> None:
    """``P`` reproduces an arclength-linear coarse field exactly (odd fibers: no trailing boundary node)."""
    foff = np.array([0, 41, 82, 121], dtype=np.int32)  # three odd-length fibers (41, 41, 39)
    level = build_arclength_restriction(foff)
    n_fine = int(foff[-1])
    slopes = np.array([0.37, -1.1, 2.4])
    intercepts = np.array([5.0, -2.0, 0.5])
    coarse = np.zeros((int(level["n_coarse"]), 3))
    coff = level["foff_coarse"].astype(np.int64)
    for f in range(foff.shape[0] - 1):
        cbase = int(coff[f])
        length_c = int(coff[f + 1]) - cbase
        for j in range(length_c):
            coarse[cbase + j] = intercepts + slopes * (2 * j)  # coarse node j sits at fine index 2j
    fine = _prolong_np(level, coarse, n_fine)
    expected = np.zeros((n_fine, 3))
    for f in range(foff.shape[0] - 1):
        begin = int(foff[f])
        length = int(foff[f + 1]) - begin
        for k in range(length):
            expected[begin + k] = intercepts + slopes * k
    assert np.allclose(fine, expected, atol=1e-12), float(np.abs(fine - expected).max())


def test_restrict_prolong_constant_recovers_full_weighting_scale() -> None:
    """``P`` maps a constant to a constant; ``R∘P`` on a constant = the (column-sum) full-weighting scale."""
    foff = np.array([0, 7, 20, 22], dtype=np.int32)  # odd (7), odd (13), 2-node
    level = build_arclength_restriction(foff)
    n_fine = int(foff[-1])
    n_coarse = int(level["n_coarse"])
    const = np.tile([2.0, -3.0, 4.0], (n_coarse, 1))
    fine = _prolong_np(level, const, n_fine)
    assert np.allclose(fine, np.tile([2.0, -3.0, 4.0], (n_fine, 1)), atol=1e-12)
    rp = _restrict_np(level, fine)
    ones_fine = np.ones((n_fine, 3))
    colsum = _restrict_np(level, ones_fine)  # column sums of P == R applied to all-ones fine
    assert np.allclose(rp, colsum * np.array([2.0, -3.0, 4.0]), atol=1e-12)
    # Interior coarse nodes of an odd-length fiber see the classic full-weighting scale 1 + 1/2 + 1/2 == 2
    # (endpoints get 1 + 1/2 == 1.5); no coarse node exceeds the full-weighting scale here.
    assert np.isclose(colsum[:, 0].max(), 2.0)
    assert np.isclose(colsum[:, 0].min(), 1.5)


# ----------------------------------------------------------------------------------------------------
# (b) transpose consistency  <R x, y>_c == <x, P y>_f   (Galerkin requirement R == P^T).
# ----------------------------------------------------------------------------------------------------


def test_restrict_and_prolong_are_transpose_consistent() -> None:
    """``<R x, y>_c == <x, P y>_f`` to round-off — R and P share one weight table so R == P^T exactly."""
    rng = np.random.default_rng(20260724)
    foff = np.array([0, 41, 45, 47, 84, 91], dtype=np.int32)  # 41 (odd), 4 (arp23), 2, 37 (odd), 7
    level = build_arclength_restriction(foff)
    n_fine = int(foff[-1])
    n_coarse = int(level["n_coarse"])
    x_fine = rng.standard_normal((n_fine, 3))
    y_coarse = rng.standard_normal((n_coarse, 3))
    lhs = float(np.sum(_restrict_np(level, x_fine) * y_coarse))   # <R x, y>_c
    rhs = float(np.sum(x_fine * _prolong_np(level, y_coarse, n_fine)))  # <x, P y>_f
    residual = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0)
    assert residual < 1e-12, residual


# ----------------------------------------------------------------------------------------------------
# (c) line smoother: it solves the analytic lumped tridiagonal M_f; damps HF, near no-op on the smooth mode.
# ----------------------------------------------------------------------------------------------------


def test_line_smoother_solves_the_lumped_tridiagonal() -> None:
    """The NumPy line smoother (omega=1, x=0) inverts EXACTLY the analytic row-sum-lumped tridiagonal ``M_f``,
    and the kernel SOURCE assembles the identical ``[2a, 4a, 2a] / -2a`` band (native-validated by the Lead)."""
    length = 9
    foff = np.array([0, length], dtype=np.int32)
    toff = np.array([0, length - 2], dtype=np.int32)
    rng = np.random.default_rng(11)
    alpha = rng.uniform(0.5, 3.0, size=length - 2)
    ext_scalar = rng.uniform(0.01, 0.1, size=length)
    ext = np.stack([ext_scalar, ext_scalar, ext_scalar], axis=1)  # identical components => M same per comp
    # Reconstruct M^-1 column-by-column from the smoother; column j = solution for residual = e_j (component 0).
    m_inv = np.zeros((length, length))
    for j in range(length):
        residual = np.zeros((length, 3))
        residual[j, 0] = 1.0
        x, ok = _line_smooth_np(residual, np.zeros((length, 3)), foff, toff, alpha, ext, omega=1.0)
        assert ok == 1
        m_inv[:, j] = x[:, 0]
    m_smoother = np.linalg.inv(m_inv)
    m_analytic = _lumped_tridiagonal(ext_scalar, alpha)
    assert np.allclose(m_smoother, m_analytic, atol=1e-9), float(np.abs(m_smoother - m_analytic).max())
    # Row sums of the BENDING part are zero => the constant near-null is annihilated (smooth-mode no-op).
    bending = m_analytic - np.diag(ext_scalar)
    assert np.allclose(bending.sum(axis=1), 0.0, atol=1e-9)
    # SOURCE: the Warp kernel builds the SAME lumped band (the two pentadiagonal corners moved to the diagonal).
    # native-validated by the Lead (device Thomas sweep on the A5000).
    assert "two_a = wp.float64(2.0) * a" in _MG_SRC
    assert "four_a = wp.float64(4.0) * a" in _MG_SRC
    assert "diag_bend[begin + t + 1] = diag_bend[begin + t + 1] + four_a" in _MG_SRC
    assert "off_sub[begin + t + 1] = off_sub[begin + t + 1] - two_a" in _MG_SRC


def _one_sweep_reduction(profile, operator, foff, toff, alpha, ext, omega):
    """One NumPy smoother sweep on the pure mode ``profile`` (b = 0 => residual = -A x); ``||x_1|| / ||x_0||``."""
    x_np = np.stack([profile, profile, profile], axis=1).astype(np.float64)
    norm0 = np.linalg.norm(x_np)
    residual = -(operator @ x_np)
    x_np, ok = _line_smooth_np(residual, x_np, foff, toff, alpha, ext, omega)
    assert ok == 1
    return float(np.linalg.norm(x_np) / norm0)


def test_line_smoother_damps_high_frequency_and_preserves_smooth() -> None:
    """Classic smoothing factor: one sweep crushes every high-arclength-frequency mode of the exact fiber
    operator while a bending-dominated soft mode is a near no-op; a few sweeps then kill an oscillatory error."""
    length = 33
    foff = np.array([0, length], dtype=np.int32)
    toff = np.array([0, length - 2], dtype=np.int32)
    alpha = np.ones(length - 2)                       # stiff, bending-dominated single fiber
    ext_scalar = np.full(length, 1.0e-3)              # small regularizer floor (a << bending)
    ext = np.stack([ext_scalar, ext_scalar, ext_scalar], axis=1)
    operator = _pentadiagonal_operator(ext_scalar, alpha)  # exact pentadiagonal fiber sub-operator
    omega = 4.0 / 3.0                                 # classic 1-D line-relaxation smoothing weight

    eigvals, eigvecs = np.linalg.eigh(operator)       # ascending: soft -> stiff arclength modes

    # (i) The smoothing factor: EVERY mode in the stiff half of the spectrum is damped by a solid factor.
    stiff_reductions = [
        _one_sweep_reduction(eigvecs[:, i], operator, foff, toff, alpha, ext, omega)
        for i in range(length // 2, length)
    ]
    assert max(stiff_reductions) < 0.4, max(stiff_reductions)

    # (ii) A bending-dominated soft mode (well above the ext floor) is a near no-op — coarse grid's job.
    soft_reduction = _one_sweep_reduction(eigvecs[:, 4], operator, foff, toff, alpha, ext, omega)
    assert soft_reduction > 0.9, soft_reduction

    # (iii) A few sweeps on the pure theta=pi oscillatory error crush it (compounded smoothing factor).
    x_np = np.tile(((-1.0) ** np.arange(length))[:, None], (1, 3)).astype(np.float64)
    norm0 = np.linalg.norm(x_np)
    for _ in range(4):
        x_np, ok = _line_smooth_np(-(operator @ x_np), x_np, foff, toff, alpha, ext, omega)
        assert ok == 1
    assert float(np.linalg.norm(x_np) / norm0) < 0.05, float(np.linalg.norm(x_np) / norm0)


# ----------------------------------------------------------------------------------------------------
# (d) heterogeneous odd + 2-node fibers coarsen + smooth without crashing.
# ----------------------------------------------------------------------------------------------------


def test_mixed_and_short_fibers_build_and_smooth() -> None:
    """Odd (5, 3), even (4), and degenerate 2-node fibers coarsen + smooth without crashing (finite stays 1)."""
    foff = np.array([0, 5, 8, 12, 14], dtype=np.int32)   # lengths 5, 3, 4, 2
    toff = np.array([0, 3, 4, 6, 6], dtype=np.int32)     # per-fiber triple counts: 3, 1, 2, 0
    n_fine = int(foff[-1])
    n_triples = int(toff[-1])

    level = build_arclength_restriction(foff)
    assert list(np.diff(level["foff_coarse"])) == [3, 2, 2, 1]  # ceil(L/2): 3, 2, 2, 1
    assert int(level["n_coarse"]) == 8

    rng = np.random.default_rng(3)
    # Round-trip R/P must stay finite and transpose-consistent on this heterogeneous set.
    x_fine = rng.standard_normal((n_fine, 3))
    y_coarse = rng.standard_normal((int(level["n_coarse"]), 3))
    lhs = float(np.sum(_restrict_np(level, x_fine) * y_coarse))
    rhs = float(np.sum(x_fine * _prolong_np(level, y_coarse, n_fine)))
    assert abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0) < 1e-12

    alpha = rng.uniform(0.5, 2.0, size=max(n_triples, 1))
    ext = np.tile([0.05, 0.05, 0.05], (n_fine, 1))
    residual = rng.standard_normal((n_fine, 3))
    x, ok = _line_smooth_np(residual, np.zeros((n_fine, 3)), foff, toff, alpha, ext, omega=1.0)
    assert ok == 1
    assert np.all(np.isfinite(x))
    # The 2-node fiber (no triples) reduces to a pure diagonal solve: x = residual / ext exactly.
    b2 = int(foff[3])
    assert np.allclose(x[b2:b2 + 2], residual[b2:b2 + 2] / 0.05, atol=1e-12)
    # SOURCE: the DEVICE Thomas sweep latches the finite predicate to 0 on a non-positive pivot (non-finite
    # input) rather than emitting a spurious correction — a kernel-only guard, native-validated by the Lead.
    assert "finite[0] = wp.int32(0)" in _MG_SRC
    assert "if b0x <= wp.float64(1.0e-300)" in _MG_SRC
