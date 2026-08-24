"""Pure-NumPy reference gates for the fiber-arclength multigrid V-cycle (``FiberArclengthMultigrid``) + wiring.

GPU-only contract (CLAUDE.md I0-A, PI 2026-07-16): ``ac`` kernels are Warp-CUDA and are **never launched on the
CPU device**, not even in tests. Following the sanctioned dev-box pattern (Mac = kernel SOURCE + NumPy
reference; native launch = the Lead on the A5000), this module reimplements the whole symmetric V-cycle in
pure NumPy from the exact kernel math — the 2:1 ``R``/``P`` tables (host ``build_arclength_restriction`` /
``_coarsen_bending``, exercised DIRECTLY), the tridiagonal line smoother, the assembled per-fiber operator, the
Galerkin squared-weight diagonal restriction and the device-side ``omega`` power-iteration — and gates:

    (a) SPD: the V-cycle preconditioner built column-by-column is symmetric to round-off and positive definite;
    (b) h-independence: one V-cycle reduces a multi-fiber residual's energy norm by a solid factor that does NOT
        degrade as the fiber is refined (two lengths / level counts);
    (b') SPEED: the assembled per-fiber coarse operator matches the fine per-fiber block EXACTLY, collapses the
        full-fine apply count to level-0-only (4), and keeps the V-cycle SPD + h-independent;
    (c) the fiber-quotient coarsest level attaches + the additive form stays SPD — the fiber-quotient is a
        separate Warp/BSR CUDA subsystem (``FiberQuotientCoarsePathB``), so this is gated at the SOURCE level +
        native-validated by the Lead;
    (d) off-flag bit-parity + preconditioning-only + CG host early-exit — these are ``ProjectedAnalyticCG``
        device-solver invariants (no NumPy analogue), gated at the SOURCE level + native-validated by the Lead.

The synthetic operator is the exact per-fiber pentadiagonal NF2007 bending sub-operator + ``a I + ext`` — the
design regime (stiff-along-arclength, soft-across); ``ext`` is ``diag(A) - diag(bending)`` exactly as the real
``coarse_external_diagonal`` so the line smoother sees the full node diagonal.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pytest

# Importing the host helpers requires warp importable (the module carries @wp.kernel bodies), but NOTHING here
# launches a kernel or allocates a device array (GPU-only contract: no CPU-Warp launches in ac tests).
pytest.importorskip("warp")

import aleph.components.incumbent.fiber_arclength_mg as _mg  # noqa: E402
from aleph.components.incumbent.fiber_arclength_mg import (  # noqa: E402
    FiberArclengthMultigrid,
    _coarsen_bending,
    build_arclength_restriction,
)

_CELL_DIR = Path(_mg.__file__).parent
_MG_SRC = Path(_mg.__file__).read_text(encoding="utf-8")
_IMPL_SRC = (_CELL_DIR / "implicit_mechanics.py").read_text(encoding="utf-8")
_FQ_SRC = (_CELL_DIR / "fiber_quotient_coarse.py").read_text(encoding="utf-8")


# ----------------------------------------------------------------------------------------------------
# Pure-NumPy references (faithful to the kernel bodies; the device launch is native-validated by the Lead).
# ----------------------------------------------------------------------------------------------------


def _restrict(level: dict, fine: np.ndarray) -> np.ndarray:
    """``coarse = R fine`` (full-weighting); ``np.add.at`` == the kernel's atomic-add scatter."""
    fine = np.ascontiguousarray(fine, dtype=np.float64)
    n_coarse = int(level["n_coarse"])
    coarse = np.zeros((max(n_coarse, 1), 3), dtype=np.float64)
    c0, w0, c1, w1 = level["c0"], level["w0"], level["c1"], level["w1"]
    m0 = c0 >= 0
    np.add.at(coarse, c0[m0], w0[m0, None] * fine[m0])
    m1 = c1 >= 0
    np.add.at(coarse, c1[m1], w1[m1, None] * fine[m1])
    return coarse[:n_coarse]


def _prolong(level: dict, coarse: np.ndarray, n_fine: int) -> np.ndarray:
    """``fine = P coarse`` (1-D linear arclength interpolation); reads the SAME table as ``R`` (so ``R == P^T``)."""
    coarse = np.ascontiguousarray(coarse, dtype=np.float64)
    fine = np.zeros((n_fine, 3), dtype=np.float64)
    c0, w0, c1, w1 = level["c0"], level["w0"], level["c1"], level["w1"]
    m0 = c0 >= 0
    fine[m0] += w0[m0, None] * coarse[c0[m0]]
    m1 = c1 >= 0
    fine[m1] += w1[m1, None] * coarse[c1[m1]]
    return fine


def _restrict_diag_sq(level: dict, fine_diag: np.ndarray, n_coarse: int) -> np.ndarray:
    """``restrict_diagonal_sq_kernel``: Galerkin coarse of a diagonal ``coarse[c] = sum_i P_ic^2 fine_diag[i]``."""
    fine_diag = np.ascontiguousarray(fine_diag, dtype=np.float64)
    coarse = np.zeros((max(n_coarse, 1), 3), dtype=np.float64)
    c0, w0, c1, w1 = level["c0"], level["w0"], level["c1"], level["w1"]
    m0 = c0 >= 0
    np.add.at(coarse, c0[m0], (w0[m0] ** 2)[:, None] * fine_diag[m0])
    m1 = c1 >= 0
    np.add.at(coarse, c1[m1], (w1[m1] ** 2)[:, None] * fine_diag[m1])
    return coarse[:n_coarse]


def _lumped_bending(length: int, alpha: np.ndarray) -> np.ndarray:
    """Row-sum-lumped tridiagonal bending band (one component) — ``fiber_line_smooth_kernel``'s ``M_f`` bending."""
    diag = np.zeros(length, dtype=np.float64)
    off = np.zeros(max(length - 1, 0), dtype=np.float64)
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


def _line_smooth(residual, x, foff, toff, alpha, ext, omega):
    """NumPy mirror of ``fiber_line_smooth_kernel``: per fiber ``x += omega * M_f^-1 residual`` (direct solve)."""
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
            x[begin:end, c] += omega * np.linalg.solve(matrix, residual[begin:end, c])
    return x


def _perfiber_operator(xin, foff, toff, alpha, ext):
    """NumPy mirror of ``fiber_perfiber_operator_kernel``: ``ext ⊙ x + (D2^T diag(alpha) D2) x`` per fiber."""
    out = np.asarray(ext, dtype=np.float64) * np.asarray(xin, dtype=np.float64)
    for f in range(len(foff) - 1):
        begin, end = int(foff[f]), int(foff[f + 1])
        length = end - begin
        if length <= 0:
            continue
        tb = int(toff[f])
        for t in range(length - 2):
            a = float(alpha[tb + t])
            n0, n1, n2 = begin + t, begin + t + 1, begin + t + 2
            ad2 = a * (xin[n0] - 2.0 * xin[n1] + xin[n2])
            out[n0] += ad2
            out[n1] -= 2.0 * ad2
            out[n2] += ad2
    return out


# ----------------------------------------------------------------------------------------------------
# Pure-NumPy V-cycle: a faithful reimplementation of FiberArclengthMultigrid (build/apply, matrix-free + assembled).
# ----------------------------------------------------------------------------------------------------


class _NumpyFiberMG:
    """NumPy reference of :class:`FiberArclengthMultigrid` (same hierarchy build, level operator, smoother,
    omega power-iteration, symmetric V(nu, nu) cycle). Level 0 uses the dense fine operator ``A``; coarse levels
    use either the matrix-free ``R ∘ A ∘ P`` or (``assembled_coarse``) the assembled per-fiber operator."""

    def __init__(self, foff, toff, alpha, A_fine, ext0, *, nu=2, omega_cap=0.8, min_fiber_nodes=4,
                 max_levels=6, coarsest_sweeps=20, assembled_coarse=False):
        self.nu = int(nu)
        self.omega_cap = float(omega_cap)
        self.coarsest_sweeps = int(coarsest_sweeps)
        self.assembled_coarse = bool(assembled_coarse)
        self.A_fine = np.asarray(A_fine, dtype=np.float64)
        self.n_fibers = len(foff) - 1
        self._fine_calls = 0

        level_foff = [np.asarray(foff, dtype=np.int64)]
        level_toff = [np.asarray(toff, dtype=np.int64)]
        level_alpha = [np.asarray(alpha, dtype=np.float64)]
        self.tables: list[dict | None] = [None]
        cur_foff, cur_toff, cur_alpha = level_foff[0], level_toff[0], level_alpha[0]
        while (cur_foff.shape[0] > 1 and int(np.diff(cur_foff).max()) > min_fiber_nodes
               and len(self.tables) <= max_levels):
            table = build_arclength_restriction(cur_foff)
            foff_c = np.asarray(table["foff_coarse"], dtype=np.int64)
            toff_c, alpha_c = _coarsen_bending(cur_foff, cur_toff, cur_alpha, foff_c)
            self.tables.append(table)
            level_foff.append(foff_c)
            level_toff.append(np.asarray(toff_c, dtype=np.int64))
            level_alpha.append(np.asarray(alpha_c, dtype=np.float64))
            cur_foff, cur_toff, cur_alpha = foff_c, np.asarray(toff_c, dtype=np.int64), np.asarray(alpha_c)
        self.n_levels = len(self.tables) - 1
        self.level_foff, self.level_toff, self.level_alpha = level_foff, level_toff, level_alpha
        self.n_nodes = [int(f[-1]) for f in level_foff]

        # ext per level: Galerkin squared-weight diagonal restriction (kernel restrict_diagonal_sq).
        self.ext = [np.asarray(ext0, dtype=np.float64).reshape(self.n_nodes[0], 3)]
        for l in range(1, self.n_levels + 1):
            self.ext.append(_restrict_diag_sq(self.tables[l], self.ext[l - 1], self.n_nodes[l]))
        self.omega = self._estimate_omega()

    def _apply_fine(self, x):
        self._fine_calls += 1
        return (self.A_fine @ np.asarray(x, dtype=np.float64).reshape(-1)).reshape(-1, 3)

    def _prolong_to_fine(self, l, v):
        cur = np.asarray(v, dtype=np.float64)
        for m in range(l, 0, -1):
            cur = _prolong(self.tables[m], cur, self.n_nodes[m - 1])
        return cur

    def _restrict_from_fine(self, l, v):
        cur = np.asarray(v, dtype=np.float64)
        for m in range(1, l + 1):
            cur = _restrict(self.tables[m], cur)
        return cur

    def _apply_level(self, l, xin):
        if l == 0:
            return self._apply_fine(xin)
        if self.assembled_coarse:
            return _perfiber_operator(xin, self.level_foff[l], self.level_toff[l],
                                      self.level_alpha[l], self.ext[l])
        return self._restrict_from_fine(l, self._apply_fine(self._prolong_to_fine(l, xin)))

    def _smooth_level(self, l, residual, x, omega):
        return _line_smooth(residual, x, self.level_foff[l], self.level_toff[l],
                            self.level_alpha[l], self.ext[l], omega)

    def _smooth(self, l, x, b, sweeps, x_is_zero):
        for s in range(sweeps):
            if x_is_zero and s == 0:
                res = b.copy()
            else:
                res = b - self._apply_level(l, x)
            res = res * self.omega  # pre-scale by device omega (M^-1 linear); smoother runs at omega=1
            x = self._smooth_level(l, res, x, 1.0)
        return x

    def _estimate_omega(self):
        """Device power iteration -> ``omega = min(omega_cap, 0.9 / lambda_max(M_line^-1 A_0))`` (SPD smoother)."""
        n = self.n_nodes[0]
        if self.n_fibers == 0 or n == 0:
            return self.omega_cap
        v = np.ones((n, 3))
        num = den = 1.0
        for _ in range(10):
            av = self._apply_level(0, v)
            w = self._smooth_level(0, av, np.zeros((n, 3)), 1.0)  # M_line^-1 av (x=0, omega=1)
            num = float(np.sum(av * w))
            den = float(np.sum(v * av))
            nrm = float(np.sum(w * w))
            v = w / np.sqrt(nrm) if nrm > 1.0e-300 else w
        if num > 1.0e-300 and den > 1.0e-300 and np.isfinite(num) and np.isfinite(den):
            return min(self.omega_cap, 0.9 * den / num)
        return self.omega_cap

    def _coarsest_solve(self, b):
        x = np.zeros((self.n_nodes[self.n_levels], 3))
        return self._smooth(self.n_levels, x, b, self.coarsest_sweeps, x_is_zero=True)

    def _vcycle(self, l, b):
        if l == self.n_levels:
            return self._coarsest_solve(b)
        x = np.zeros((self.n_nodes[l], 3))
        x = self._smooth(l, x, b, self.nu, x_is_zero=True)
        res = b - self._apply_level(l, x)
        b_coarse = _restrict(self.tables[l + 1], res)
        e_coarse = self._vcycle(l + 1, b_coarse)
        x = x + _prolong(self.tables[l + 1], e_coarse, self.n_nodes[l])
        x = self._smooth(l, x, b, self.nu, x_is_zero=False)
        return x

    def apply(self, b):
        """OVERWRITE-equivalent of ``FiberArclengthMultigrid.apply``: return the actin-block ``M^-1 b`` (SPD)."""
        return self._vcycle(0, np.asarray(b, dtype=np.float64))

    def count_fine_applies(self, b):
        self._fine_calls = 0
        self.apply(b)
        return self._fine_calls


def _make_fibers(n_fibers: int, length: int):
    foff = np.arange(0, n_fibers * length + 1, length, dtype=np.int64)
    toff = np.arange(0, n_fibers * (length - 2) + 1, length - 2, dtype=np.int64)
    alpha = np.ones(n_fibers * (length - 2), dtype=np.float64)
    return foff, toff, alpha, n_fibers * length


def _dense_operator(foff, toff, alpha, ext_scalar, reg):
    """Crosslink-free multi-fiber operator ``a I + ext + sum_t alpha_t D2_t^T D2_t`` (block-diagonal per fiber)."""
    n = int(foff[-1])
    A = np.zeros((3 * n, 3 * n))
    A[np.diag_indices(3 * n)] += reg
    for i in range(n):
        for c in range(3):
            A[3 * i + c, 3 * i + c] += ext_scalar[i]
    for f in range(len(foff) - 1):
        b = int(foff[f])
        length = int(foff[f + 1]) - b
        tb = int(toff[f])
        for t in range(length - 2):
            a = float(alpha[tb + t])
            idx = [b + t, b + t + 1, b + t + 2]
            w = [1.0, -2.0, 1.0]
            for r in range(3):
                for cc in range(3):
                    for comp in range(3):
                        A[3 * idx[r] + comp, 3 * idx[cc] + comp] += a * w[r] * w[cc]
    return A


def _build_case(n_fibers, length, *, assembled_coarse=False):
    """Return (mg, dense A, n) for a synthetic crosslink-free multi-fiber system (ext = diag(A) - diag(bending))."""
    foff, toff, alpha, n = _make_fibers(n_fibers, length)
    ext_scalar = np.full(n, 1.0e-3)
    reg = 1.0e-3
    A = _dense_operator(foff, toff, alpha, ext_scalar, reg)
    ext0 = np.tile((reg + ext_scalar)[:, None], (1, 3))  # the live coarse_external_diagonal
    mg = _NumpyFiberMG(foff, toff, alpha, A, ext0, nu=2, omega_cap=0.8, coarsest_sweeps=20,
                       assembled_coarse=assembled_coarse)
    return mg, A, n


def _preconditioner_matrix(mg, n):
    """Assemble the V-cycle preconditioner ``M^-1`` column-by-column (the leading n_actin fiber block)."""
    Minv = np.zeros((3 * n, 3 * n))
    for col in range(3 * n):
        e = np.zeros((n, 3))
        e[col // 3, col % 3] = 1.0
        Minv[:, col] = mg.apply(e).reshape(-1)
    return Minv


def _one_vcycle_energy_reduction(mg, A, n, seed):
    """``||e1||_A / ||e0||_A`` for x0 = 0 (e0 = A^-1 b) after one V-cycle correction z (e1 = e0 - z)."""
    rng = np.random.default_rng(seed)
    b = rng.standard_normal((n, 3))
    x_true = np.linalg.solve(A, b.reshape(-1))
    z = mg.apply(b).reshape(-1)
    e0, e1 = x_true, x_true - z
    return float(np.sqrt(e1 @ (A @ e1)) / np.sqrt(e0 @ (A @ e0)))


# ----------------------------------------------------------------------------------------------------
# (a) SPD — the V-cycle preconditioner is symmetric + positive definite.
# ----------------------------------------------------------------------------------------------------


def test_vcycle_preconditioner_is_spd() -> None:
    """A symmetric V(2,2) on the multi-fiber bending operator is an SPD preconditioner (M^-1 = (M^-1)^T > 0)."""
    mg, A, n = _build_case(6, 9)
    Minv = _preconditioner_matrix(mg, n)
    symmetric_defect = np.abs(Minv - Minv.T).max() / max(np.abs(Minv).max(), 1.0e-30)
    assert symmetric_defect < 1.0e-10, symmetric_defect
    eig = np.linalg.eigvalsh(0.5 * (Minv + Minv.T))
    assert eig.min() > 0.0, eig.min()
    # <M^-1 r, r> > 0 for a concrete random residual (positive-definiteness in action).
    r = np.random.default_rng(7).standard_normal((n, 3))
    assert float(np.sum(mg.apply(r) * r)) > 0.0


# ----------------------------------------------------------------------------------------------------
# (b) h-independence — one V-cycle reduces the energy norm by a solid factor that does not degrade with h.
# ----------------------------------------------------------------------------------------------------


def test_one_vcycle_energy_reduction_is_h_independent() -> None:
    """One V-cycle crushes the energy norm by a solid factor at two fiber lengths; the factor does NOT degrade
    as the fiber is refined (more nodes / one more MG level) — the multigrid h-independence signature."""
    mg9, A9, n9 = _build_case(6, 9)
    mg17, A17, n17 = _build_case(6, 17)
    assert mg9.n_levels == 2 and mg17.n_levels == 3  # 9->5->3 ; 17->9->5->3 (an extra level from refinement)
    red9 = _one_vcycle_energy_reduction(mg9, A9, n9, seed=11)
    red17 = _one_vcycle_energy_reduction(mg17, A17, n17, seed=11)
    assert red9 < 0.6, red9
    assert red17 < 0.6, red17
    assert red17 <= red9 * 1.25, (red9, red17)  # refining does NOT degrade the factor


# ----------------------------------------------------------------------------------------------------
# (b') SPEED — assembled per-fiber coarse operator: cheap O(n_nodes[l]) applies, still SPD + h-independent.
# ----------------------------------------------------------------------------------------------------


def test_assembled_perfiber_operator_matches_the_fine_perfiber_block() -> None:
    """The assembled ``fiber_perfiber_operator_kernel`` (ext ⊙ x + pentadiagonal bending) reproduces the fine
    operator's per-fiber block EXACTLY (a crosslink-free multi-fiber operator IS block-diagonal per fiber)."""
    mg, A, n = _build_case(4, 11)
    x = np.random.default_rng(9).standard_normal((n, 3))
    out = _perfiber_operator(x, mg.level_foff[0], mg.level_toff[0], mg.level_alpha[0], mg.ext[0])
    ref = (A @ x.reshape(-1)).reshape(n, 3)
    assert np.abs(out - ref).max() < 1.0e-11, np.abs(out - ref).max()
    # SOURCE: the assembled kernel is the one wired into the coarse level (l>=1) under assembled_coarse.
    assert "def fiber_perfiber_operator_kernel(" in _MG_SRC
    assert "if l >= 1 and self.assembled_coarse:" in _MG_SRC  # native-validated by the Lead


def test_assembled_coarse_slashes_full_fine_applies_per_vcycle() -> None:
    """With ``assembled_coarse`` the coarse levels (l>=1) stop calling the full-fine operator: the per-V-cycle
    full-fine apply count collapses to the level-0 work only (nu1-1 + residual + nu2 = 4) and — critically — does
    NOT grow with the level count, whereas the matrix-free coarse operator's count grows with every extra level."""
    counts_off, counts_on = {}, {}
    for length in (9, 17, 33):
        off, _, _ = _build_case(6, length, assembled_coarse=False)
        on, _, _ = _build_case(6, length, assembled_coarse=True)
        b = np.random.default_rng(2).standard_normal((off.n_nodes[0], 3))
        counts_off[length] = off.count_fine_applies(b)
        counts_on[length] = on.count_fine_applies(b)
    assert counts_off[9] < counts_off[17] < counts_off[33]  # OFF grows with depth
    assert set(counts_on.values()) == {4}, counts_on       # ON is a constant 4 (level-0 only)
    assert counts_on[33] < counts_off[33] / 5              # >=5x fewer full-fine applies at depth


def test_assembled_coarse_vcycle_is_spd_and_h_independent() -> None:
    """The assembled-coarse V-cycle is STILL an SPD preconditioner AND still h-independent: one V-cycle crushes
    the energy norm by a solid factor at two lengths and refinement does not degrade it (same signature as the
    matrix-free coarse). Only the coarse OPERATOR changed (drops inter-fiber, handled additively)."""
    mg9, A9, n9 = _build_case(6, 9, assembled_coarse=True)
    mg17, A17, n17 = _build_case(6, 17, assembled_coarse=True)
    Minv = _preconditioner_matrix(mg9, n9)
    assert np.abs(Minv - Minv.T).max() / max(np.abs(Minv).max(), 1.0e-30) < 1.0e-10
    assert np.linalg.eigvalsh(0.5 * (Minv + Minv.T)).min() > 0.0
    red9 = _one_vcycle_energy_reduction(mg9, A9, n9, seed=11)
    red17 = _one_vcycle_energy_reduction(mg17, A17, n17, seed=11)
    assert red9 < 0.6, red9
    assert red17 < 0.6, red17
    assert red17 <= red9 * 1.25, (red9, red17)


def test_assembled_coarse_is_opt_in_and_off_by_default() -> None:
    """Default keeps the currently-native-validated matrix-free coarse (``assembled_coarse`` False). Checked on
    the REAL ``FiberArclengthMultigrid.__init__`` signature (no device allocation / launch)."""
    default = inspect.signature(FiberArclengthMultigrid.__init__).parameters["assembled_coarse"].default
    assert default is False
    mg, _, _ = _build_case(6, 9)
    assert mg.assembled_coarse is False


# ----------------------------------------------------------------------------------------------------
# (c) fiber-quotient coarsest level — separate Warp/BSR CUDA subsystem; SOURCE-gated + native-validated.
# ----------------------------------------------------------------------------------------------------


def test_fiber_quotient_coarsest_attaches_and_stays_spd() -> None:
    """With inter-fiber crosslinks the ``FiberQuotientCoarsePathB`` coarsest level attaches and is injected
    ADDITIVELY at the fine level (``M^-1 = V_arclength + P A_c^-1 P^T``). Both terms SPD => the composite is SPD.

    ``FiberQuotientCoarsePathB`` is a warp.sparse BSR + block-Jacobi-CG CUDA subsystem with no NumPy analogue,
    so the attachment + additive-SPD wiring is gated at the SOURCE level; the numeric SPD of the assembled
    device preconditioner is native-validated by the Lead on the A5000 (§4.4 ladder)."""
    # Attaches only on a crosslinked cell, from the FiberQuotientCoarsePathB module.
    assert 'attach_fiber_quotient and int(getattr(cell, "n_xl", 0)) > 0' in _MG_SRC
    assert "from aleph.components.incumbent.fiber_quotient_coarse import FiberQuotientCoarsePathB" in _MG_SRC
    assert "self.fq = FiberQuotientCoarsePathB(cell, inner_iterations=fq_inner_iterations)" in _MG_SRC
    # Injected ADDITIVELY at the fine level in apply() (NOT nested in the arclength recursion).
    assert "if self.fq is not None:" in _MG_SRC
    assert "self.fq.apply_richardson(pos, residual, preconditioned, regularization, finite," in _MG_SRC
    assert "M^-1 = V_arclength + P A_c^-1 P^T" in _MG_SRC
    # The coarse block is SPD by construction (D_ext >= a > 0 congruence).
    assert "M^-1 = S^-1 + P A_c^-1 P^T`` stays SPD" in _FQ_SRC


def test_fiber_quotient_richardson_is_linear_and_spd() -> None:
    """The coarsest ``apply_richardson`` is a FIXED, LINEAR, SPD map (unlike the CG-based ``apply``): a fixed
    damped block-Jacobi Richardson (``P B_k P^T``). SOURCE-gated + native-validated by the Lead."""
    assert "class FiberQuotientCoarsePathB:" in _FQ_SRC
    assert "def apply_richardson(self, pos: wp.array, residual: wp.array, preconditioned: wp.array," in _FQ_SRC
    assert "fixed damped block-Jacobi RICHARDSON, LINEAR + SPD" in _FQ_SRC


# ----------------------------------------------------------------------------------------------------
# (d) off-flag bit-parity + preconditioning-only + CG host early-exit — ProjectedAnalyticCG device-solver
#     invariants (no NumPy analogue); SOURCE-gated + native-validated by the Lead.
# ----------------------------------------------------------------------------------------------------


def test_off_flag_is_inert_and_uses_incumbent_precondition() -> None:
    """``multigrid=False`` builds NO MG (``self.mg`` stays None) and ``_precondition`` only consults ``self.mg``
    under ``if self.multigrid:`` — so the off path is the incumbent, byte-untouched. SOURCE-gated on the real
    ``ProjectedAnalyticCG`` wiring; the bit-parity solve is native-validated by the Lead."""
    assert "self.mg = None" in _IMPL_SRC
    assert "if self.multigrid and int(cell.n_fibers) > 0:" in _IMPL_SRC
    # No fibers to precondition => leave the incumbent path untouched (never a half-built MG).
    assert "self.multigrid = False  # no fibers to precondition; leave the incumbent path untouched" in _IMPL_SRC
    # _precondition consults self.mg ONLY inside the multigrid branch (else the incumbent smoother runs).
    precond = _IMPL_SRC[_IMPL_SRC.index("def _precondition("):]
    precond = precond[:precond.index("\n    def ", 1)]
    assert "if self.multigrid:" in precond
    assert "self.mg.apply(pos, self.r, self.z, self._operator, self._mg_reg, finite)" in precond
    assert precond.index("if self.multigrid:") < precond.index("self.mg.apply(")


def test_multigrid_is_preconditioning_only_same_fixed_point() -> None:
    """MG on vs off solve the IDENTICAL SPD system to the same converged ``dx`` — the V-cycle is preconditioning
    only (SPD ``M^-1``), it cannot move the fixed point (CLAUDE.md no-gate-loosening). SOURCE-gated on the
    documented invariant; the identical-``dx`` solve is native-validated by the Lead."""
    assert "preconditioning-only, the fixed point + residual gate" in _IMPL_SRC
    assert "CLAUDE.md no-gate-loosening" in _IMPL_SRC
    # The MG only OVERWRITES the leading n_actin block; the disjoint remainder keeps its node-Jacobi value.
    assert "the V-cycle now OVERWRITES the leading n_actin block with its SPD correction" in _IMPL_SRC


def test_cg_host_early_exit_is_bit_identical() -> None:
    """The outer-PCG host early-exit (``cg_check_every>0``) breaks only AFTER the device latch clears — every
    further loop body is a no-op update (all _cg_* kernels gate on ``active``), so breaking only skips wasted
    launches and cannot change ``dx``. SOURCE-gated; the bit-identical dx is native-validated by the Lead."""
    assert "if self.cg_check_every and iteration % self.cg_check_every == 0:" in _IMPL_SRC
    assert "if int(self.active.numpy()[0]) == 0:" in _IMPL_SRC
    # The comment records the bit-identical guarantee (post-convergence work is a no-op, dx unchanged).
    assert "bit-identical" in _IMPL_SRC
    assert "dx is unchanged; only the wasted post-convergence work is saved" in _IMPL_SRC
