"""Constrained Leimkuhler-Matthews BAOAB-limit overdamped Langevin (rigid bonds).

Sibling of the FROZEN ``baoab.py`` (CLAUDE.md §integrator-freeze — this is a
NEW module, baoab.py is untouched). Replaces the stiff backbone *stretch*
harmonic bond with a RIGID DISTANCE CONSTRAINT so the CFL timestep is set by
the bending mode (τ_bend ≈ 698 μs) instead of the stretch mode
(τ_stretch ≈ 130 ns) — a ~5357× headroom that lets the multi-second KU-3.x
emergent-physics gates (rounding / tension / blebbistatin) become
minutes-to-hours runs instead of months.

This is *mechanistically more correct*, not a shortcut (CLAUDE.md hard rule
+ ``feedback_acs_no_abstractions``): an actin segment's axial stiffness
≫ its bending stiffness, so it is effectively inextensible. The 130 ns
stretch vibration is a numerical proxy with zero contribution to the
μs–ms rounding/tension physics; a rigid bond is the fine-grained truth.

Algorithm (per step), positions ``r``, drag γ_i, mobility M_i = 1/γ_i,
holonomic constraints ``σ_a(r) = |r_{i(a)} − r_{j(a)}|² − ℓ_a² = 0``:

1. **Fixman pseudo-force** ``F_F = −∇ U_F``, ``U_F = ½ kT ln det G`` where
   ``G_ab = Σ_i M_i ∇_i σ_a · ∇_i σ_b`` is the constraint mobility metric.
   For a chain of rigid bonds with a *soft* bending angle, ``det G`` is
   configuration-dependent **through the inter-bond angles** (G is
   tridiagonal per filament, off-diagonal ∝ −cos ψ), so this correction
   is NOT negligible — without it the L_p / equipartition gates are
   biased (Fixman 1974; Hinch 1994). PI-ratified route 1 (2026-05-28):
   apply the explicit pseudo-force.
2. **L-M predictor** (identical formula to ``baoab.py``; here the force is
   ``net_force + F_F`` and carries NO stretch-bond term — the stretch DOF
   is the constraint):
   ``r̃_i = r_i + (F_i/γ_i)·Δt + √(kT/(2 γ_i Δt))·(W_n + W_{n-1})·Δt``.
3. **SHAKE projection** ``r̃ → r⁺`` satisfying every ``σ_a(r⁺)=0`` by
   Gauss-Seidel iteration of mobility-weighted Lagrange corrections.
   (Overdamped BD has no velocities → the RATTLE velocity half collapses;
   only the position projection remains.)
4. Minimum-image wrap into the box (reuses ``baoab._wrap_into_box``).

Sanity Gate
-----------
*Per CLAUDE.md "Sanity Gate Protocol mandatory before first execution."
STATIC + analytic checks in ``tests/test_constrained_baoab.py``;
RUNTIME checks asserted in ``act`` / ``attach``.*

1. **Dimensional analysis** — SHAKE correction
   ``Δr_i = M_i·g·d0`` with ``g = (|s|²−ℓ²)/(2 (s·d0)(M_i+M_j))``:
   ``[m²]/([m²]·[m/N·s]⁻¹... )`` reduces so ``Δr_i`` is [m]. The Fixman
   force ``F_F = −½kT ∇ ln det G`` is [J]·[1/m] = [N]. STATIC (sympy +
   finite-difference) checks both.
2. **Boundary cases** — empty constraint set ⇒ trajectory identical to
   ``baoab.make_baoab_updater`` with the same seed (RUNTIME + STATIC).
   Dimer (1 bond, 0 angles) ⇒ ``det G`` const ⇒ zero Fixman force.
   N=2 chain ⇒ no angle.
3. **Conservation** — constraint corrections are internal (Newton 3rd
   law): for equal mobility the mobility-weighted COM is unchanged by
   SHAKE. STATIC.
4. **Numerical sanity** — constraint drift ``max_a ||b_a|−ℓ_a| ≤ tol``
   asserted every step (RUNTIME, mirrors baoab NaN guard). SHAKE
   non-convergence raises.
5. **Sign-sense** — a stretched rigid bond (|s|>ℓ) is pulled in by SHAKE;
   compressed is pushed out. STATIC.
6. **Measurement-protocol consistency** — the decisive Fixman test: a
   rigid-bond + soft-angle trimer must reproduce the SAME bending-angle
   PDF as the stiff-harmonic-bond + soft-angle trimer (and a free rigid
   dimer must diffuse with D_com = kT/(γ_i+γ_j) and hold |b| fixed).
   These prove the Fixman sign + magnitude. The H.2 L_p strict bands and
   H.3 equipartition ⟨E_bend⟩ vs 0.9898 kT are the production re-validation
   (Milestone 2/3) — those touch the oracle ``tau_min`` contract and are
   gated on separate PI sign-off.

References
----------
- Fixman 1974, PNAS 71(8):3050 (metric tensor / pseudo-potential).
- Hinch 1994, JFM 271:219 (Brownian dynamics with constraints).
- ``ffn_sim/integrator/baoab.py`` (the unconstrained L-M contract).
- PHASE_0_3_DECISIONS §D3; CLAUDE.md §integrator-freeze.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

import hoomd
import hoomd.custom

from ffn_sim.integrator.baoab import _wrap_into_box

# Sign of the Fixman pseudo-potential U_F = FIXMAN_SIGN · ½ kT ln det G.
# +1 is the textbook convention (Fixman 1978; Hinch 1994): naive constrained
# Brownian dynamics over-samples by det(G)^{+1/2}, and U_F = +½kT ln det G
# cancels it to recover the unconstrained (stiff-spring) marginal —
# exp(-(U+U_F)/kT)·det(G)^{1/2} ∝ exp(-U/kT). Confirmed analytically for the
# trimer: naive-rigid ⟨E_bend⟩=0.8691 (det T^{+1/2}) → +1 Fixman → 0.8390 =
# flexible target. NOTE: the trimer CANNOT arbitrate the sign empirically —
# the metric effect there is only ±3–4%, below the MD seed noise (~6–8% at
# 2e6 steps), so it is a consistency check only. The DECISIVE empirical sign
# confirmation is Milestone 2 (single-filament L_p, 19 cumulative angles +
# the tight H.2 strict band). (A prior −1 flip on 2026-05-28 was a
# noise-driven over-conclusion, reverted.)
FIXMAN_SIGN: float = +1.0


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def _min_image_orthorhombic(dr: np.ndarray, box) -> np.ndarray:
    """Minimum-image displacement for an orthorhombic (cube) box.

    Accepts EITHER a HOOMD box object (legacy) or a precomputed (3,)
    ``np.ndarray`` of ``[Lx, Ly, Lz]`` (preferred — avoids repeated
    HOOMD-property getattr that dominated profile traces). Cortex /
    single-filament boxes are cubes with no tilt; rigid bonds are never
    used with Lees-Edwards shear, so orthorhombic is exact. ℓ₀ ≪ L
    guarantees the nearest image is the physical bond.
    """
    if isinstance(box, np.ndarray):
        L = box
    else:
        L = np.array([box.Lx, box.Ly, box.Lz], dtype=np.float64)
    return dr - L * np.round(dr / L)


def _box_L(box) -> np.ndarray:
    """Cache HOOMD box → (3,) numpy array once per act() to avoid
    per-call getattr overhead (the dominant cost in the profile)."""
    return np.array([box.Lx, box.Ly, box.Lz], dtype=np.float64)


# ---------------------------------------------------------------------------
# SHAKE projection (pure function — unit-testable without HOOMD)
# ---------------------------------------------------------------------------
def shake_project(
    pred_pos: np.ndarray,
    ref_pos: np.ndarray,
    pairs: np.ndarray,
    lengths: np.ndarray,
    inv_mass: np.ndarray,
    box: hoomd.box.Box,
    *,
    tol: float = 1.0e-10,
    max_iter: int = 500,
) -> np.ndarray:
    """Project ``pred_pos`` onto the rigid-bond constraint manifold.

    Gauss-Seidel SHAKE: for each constraint a=(i,j) with target length
    ℓ_a, using the reference (start-of-step) bond vector ``d0`` as the
    gradient direction (standard SHAKE), iterate mobility-weighted
    Lagrange corrections until ``max_a ||b_a| − ℓ_a| / ℓ_a ≤ tol``.

    Parameters
    ----------
    pred_pos : (N, 3) post-predictor positions (modified copy returned).
    ref_pos  : (N, 3) start-of-step positions (gradient reference d0).
    pairs    : (C, 2) row-indices of the two beads of each constraint.
    lengths  : (C,)  target bond length ℓ_a.
    inv_mass : (N,)  mobility M_i = 1/γ_i (constraint correction weight).
    box      : HOOMD box (for minimum-image bond vectors).

    Returns
    -------
    (N, 3) projected positions satisfying all constraints to ``tol``.
    """
    pos = pred_pos.copy()
    i = pairs[:, 0]
    j = pairs[:, 1]
    L2 = lengths * lengths
    # Reference bond vectors (minimum-image) — fixed across iterations.
    d0 = _min_image_orthorhombic(ref_pos[i] - ref_pos[j], box)
    Mi = inv_mass[i]
    Mj = inv_mass[j]
    for _ in range(max_iter):
        s = _min_image_orthorhombic(pos[i] - pos[j], box)
        diff = np.einsum("ab,ab->a", s, s) - L2
        if np.max(np.abs(diff) / L2) <= tol:
            return pos
        # g_a = diff / (2 (s·d0)(M_i + M_j)); Gauss-Seidel applies each
        # constraint in turn so shared beads see updated positions.
        sd0 = np.einsum("ab,ab->a", s, d0)
        denom = 2.0 * sd0 * (Mi + Mj)
        g = diff / denom
        # Sequential (Gauss-Seidel) update — a chain's tridiagonal coupling
        # converges in O(few) sweeps; vectorised Jacobi would be slower.
        for a in range(pairs.shape[0]):
            corr = g[a] * d0[a]
            pos[i[a]] -= Mi[a] * corr
            pos[j[a]] += Mj[a] * corr
    raise RuntimeError(
        f"SHAKE failed to converge in {max_iter} iters; "
        f"max relative drift = {float(np.max(np.abs(diff) / L2)):.3e} > tol={tol}."
    )


def _thomas(sub: np.ndarray, diag: np.ndarray, sup: np.ndarray,
            rhs: np.ndarray) -> np.ndarray:
    """Solve a tridiagonal system (Thomas algorithm). sub[0], sup[-1] unused."""
    n = diag.shape[0]
    cp = np.empty(n); dp = np.empty(n)
    cp[0] = sup[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for k in range(1, n):
        m = diag[k] - sub[k] * cp[k - 1]
        cp[k] = sup[k] / m
        dp[k] = (rhs[k] - sub[k] * dp[k - 1]) / m
    x = np.empty(n)
    x[-1] = dp[-1]
    for k in range(n - 2, -1, -1):
        x[k] = dp[k] - cp[k] * x[k + 1]
    return x


def _thomas_batched(sub: np.ndarray, diag: np.ndarray, sup: np.ndarray,
                    rhs: np.ndarray) -> np.ndarray:
    """Batched tridiagonal solve over the leading axis. All (F, m)."""
    F, m = diag.shape
    cp = np.empty((F, m)); dp = np.empty((F, m))
    cp[:, 0] = sup[:, 0] / diag[:, 0]
    dp[:, 0] = rhs[:, 0] / diag[:, 0]
    for k in range(1, m):
        den = diag[:, k] - sub[:, k] * cp[:, k - 1]
        cp[:, k] = sup[:, k] / den
        dp[:, k] = (rhs[:, k] - sub[:, k] * dp[:, k - 1]) / den
    x = np.empty((F, m))
    x[:, -1] = dp[:, -1]
    for k in range(m - 2, -1, -1):
        x[:, k] = dp[:, k] - cp[:, k] * x[:, k + 1]
    return x


def shake_project_chains(
    pred_pos: np.ndarray,
    ref_pos: np.ndarray,
    chains: "list[np.ndarray]",
    rest_length: float,
    inv_mass: np.ndarray,
    box: hoomd.box.Box,
    *,
    tol: float = 1.0e-10,
    max_iter: int = 100,
    return_lambdas: bool = False,
):
    """Matrix-SHAKE (tridiagonal Newton) for linear-chain bond constraints.

    Same convention as :func:`shake_project` (mobility-weighted corrections
    along the reference-configuration bond gradient ``d0``), but for each
    chain it solves the tridiagonal constraint Jacobian ``J Δλ = −g`` with
    the Thomas algorithm and iterates Newton steps — converging in a handful
    of iterations instead of the O(N²) Gauss-Seidel sweeps. ``chains`` are
    ordered bead row-index arrays; all bonds use length ``rest_length``.

    Per-bead correction with multipliers λ (bond a between beads p_a, p_{a+1},
    reference bond vector ``d0_a = r_ref[p_a] − r_ref[p_{a+1}]``):
        Δr[p_k] = M_{p_k} (λ_{k-1} d0_{k-1} − λ_k d0_k).
    Newton linearisation of ``|s_a + Δ|² = ℓ²`` gives the tridiagonal system
        2 s_a·[M_{p_a} d0_{a-1} λ_{a-1}
               − (M_{p_a}+M_{p_{a+1}}) d0_a λ_a
               + M_{p_{a+1}} d0_{a+1} λ_{a+1}] = −g_a .

    Lagrange-multiplier exposure (PI 2026-05-28 verbal, R1 in
    RIGID_LAGRANGE_TENSION_DESIGN.md): when ``return_lambdas=True`` the
    accumulated per-bond Lagrange multiplier ``λ_total = Σ_iter Δλ`` is
    returned alongside the projected positions. ``λ_total`` carries the
    full constraint impulse applied this step (units: ``[M] · [r] = [r]/[F]``
    in our scaled mobility convention — the consumer converts to per-bond
    scalar tension via ``T_bond = λ_total · r₀ / Δt`` and projects onto
    the cut normal for method-of-planes γ).

    Backward-compat: ``return_lambdas`` defaults to False; existing callers
    that only expect a position array get bit-for-bit prior behaviour. The
    only in-tree caller (``ConstrainedLeimkuhlerMatthewsBAOAB.act``) opts
    in via the flag and unpacks the tuple.
    """
    pos = pred_pos.copy()
    L2 = rest_length * rest_length

    # ---- Vectorised fast path: all chains the same length (cortex N=7) ----
    # Accept EITHER a pre-stacked (F, N) ndarray (fast path — saves ~10%
    # per profile) OR a list of per-chain arrays (legacy).
    if isinstance(chains, np.ndarray) and chains.ndim == 2 and chains.shape[1] >= 3:
        P = chains; F, Np1 = P.shape; m = Np1 - 1
        _uniform_chains = True
    else:
        lens = [np.asarray(c).shape[0] for c in chains]
        _uniform_chains = bool(chains) and len(set(lens)) == 1 and lens[0] >= 3
        if _uniform_chains:
            P = np.stack([np.asarray(c) for c in chains], axis=0)
            F, Np1 = P.shape; m = Np1 - 1
    if _uniform_chains:
        M = inv_mass[P]                                          # (F, m+1)
        d0 = _min_image_orthorhombic(ref_pos[P[:, :-1]] - ref_pos[P[:, 1:]], box)
        # Accumulated per-bond Lagrange multiplier across Newton iterations.
        # Allocated only when the caller actually needs it so the default
        # path stays zero-allocation extra.
        lambda_total = np.zeros((F, m)) if return_lambdas else None
        for _ in range(max_iter):
            s = _min_image_orthorhombic(pos[P[:, :-1]] - pos[P[:, 1:]], box)
            g = np.einsum("fab,fab->fa", s, s) - L2             # (F, m)
            if np.max(np.abs(g)) / L2 <= tol:
                break
            sd = np.einsum("fab,fab->fa", s, d0)
            diag = -2.0 * (M[:, :-1] + M[:, 1:]) * sd           # (F, m)
            sub = np.zeros((F, m)); sup = np.zeros((F, m))
            if m > 1:
                sub[:, 1:] = 2.0 * M[:, 1:-1] * np.einsum(
                    "fab,fab->fa", s[:, 1:], d0[:, :-1])
                sup[:, :-1] = 2.0 * M[:, 1:-1] * np.einsum(
                    "fab,fab->fa", s[:, :-1], d0[:, 1:])
            lam = _thomas_batched(sub, diag, sup, -g)            # (F, m)
            if lambda_total is not None:
                lambda_total += lam
            disp = np.zeros((F, m + 1, 3))
            disp[:, :-1] -= (M[:, :-1] * lam)[:, :, None] * d0
            disp[:, 1:] += (M[:, 1:] * lam)[:, :, None] * d0
            pos[P] += disp                                       # chains disjoint
        else:
            drift = float(np.max(np.abs(g)) / L2)
            raise RuntimeError(
                f"M-SHAKE (vectorised, {F} chains len {m + 1}) failed in "
                f"{max_iter} iters; max relative drift = {drift:.3e} > tol={tol}.")
        if return_lambdas:
            return pos, lambda_total
        return pos

    # ---- Ragged fallback: per-chain (mixed lengths) ----
    # Per-chain lambda accumulators (list of (m,) arrays) when requested.
    lambda_per_chain = [] if return_lambdas else None
    for chain in chains:
        p = np.asarray(chain)
        m = p.shape[0] - 1
        if m == 0:
            if return_lambdas:
                lambda_per_chain.append(np.zeros(0))
            continue
        M = inv_mass[p]                                   # (m+1,)
        d0 = _min_image_orthorhombic(ref_pos[p[:-1]] - ref_pos[p[1:]], box)  # (m,3)
        lambda_total = np.zeros(m) if return_lambdas else None
        ok = False
        for _ in range(max_iter):
            s = _min_image_orthorhombic(pos[p[:-1]] - pos[p[1:]], box)       # (m,3)
            g = np.einsum("ab,ab->a", s, s) - L2
            if np.max(np.abs(g)) / L2 <= tol:
                ok = True
                break
            sd_diag = np.einsum("ab,ab->a", s, d0)                          # s_a·d0_a
            # tridiagonal coefficients
            diag = -2.0 * (M[:-1] + M[1:]) * sd_diag                        # (m,)
            sub = np.zeros(m); sup = np.zeros(m)
            if m > 1:
                sd_lower = np.einsum("ab,ab->a", s[1:], d0[:-1])            # s_a·d0_{a-1}
                sub[1:] = 2.0 * M[1:-1] * sd_lower
                sd_upper = np.einsum("ab,ab->a", s[:-1], d0[1:])           # s_a·d0_{a+1}
                sup[:-1] = 2.0 * M[1:-1] * sd_upper
            lam = _thomas(sub, diag, sup, -g)
            if lambda_total is not None:
                lambda_total += lam
            # apply Δr[p_k] = M_k (λ_{k-1} d0_{k-1} − λ_k d0_k)
            disp = np.zeros((m + 1, 3))
            disp[:-1] -= (M[:-1] * lam)[:, None] * d0      # −λ_a d0_a on bead p_a
            disp[1:] += (M[1:] * lam)[:, None] * d0        # +λ_a d0_a on bead p_{a+1}
            pos[p] += disp
        if not ok:
            s = _min_image_orthorhombic(pos[p[:-1]] - pos[p[1:]], box)
            drift = float(np.max(np.abs(np.einsum("ab,ab->a", s, s) - L2)) / L2)
            raise RuntimeError(
                f"M-SHAKE chain (len {m + 1}) failed in {max_iter} iters; "
                f"max relative drift = {drift:.3e} > tol={tol}."
            )
        if return_lambdas:
            lambda_per_chain.append(lambda_total)
    if return_lambdas:
        return pos, lambda_per_chain
    return pos


# ---------------------------------------------------------------------------
# Fixman metric pseudo-force (pure function — unit-testable without HOOMD)
# ---------------------------------------------------------------------------
def fixman_logdet_and_force(
    pos: np.ndarray,
    chains: Sequence[np.ndarray],
    kT: float,
    inv_gamma: np.ndarray,
    box: hoomd.box.Box,
) -> tuple[float, np.ndarray]:
    """Fixman pseudo-potential ``U_F`` and force ``F_F = −∇U_F``.

    For each chain (ordered bead row-indices ``[p0..pm]``, m rigid bonds)
    builds the mobility-metric Gram matrix ``G`` (m×m, tridiagonal):

        G_aa     =  4 |b_a|² (M_{p_a} + M_{p_{a+1}})
        G_{a,a+1}= −4 M_{p_{a+1}} (b_a · b_{a+1})

    and accumulates ``U_F = FIXMAN_SIGN · ½ kT Σ_chains ln det G`` and its
    analytic gradient. Chains with < 2 bonds have a configuration-
    independent ``det G`` (bond lengths are fixed) → zero force.

    Returns
    -------
    (U_F, force) with force shape (N, 3).
    """
    N = pos.shape[0]
    force = np.zeros((N, 3), dtype=np.float64)
    U_F = 0.0
    half_kT = 0.5 * kT

    # ---- Vectorised fast path: all chains the same length (cortex N=7) ----
    # Accept pre-stacked (F, N) ndarray (avoids per-step np.stack — was the
    # #2 profile hotspot) OR a list of per-chain arrays.
    _is_stacked = isinstance(chains, np.ndarray) and chains.ndim == 2 and chains.shape[1] >= 3
    if _is_stacked:
        P = chains; F = P.shape[0]; m = P.shape[1] - 1
        _uniform_chains = True
    else:
        lens = [np.asarray(c).shape[0] for c in chains]
        _uniform_chains = bool(chains) and len(set(lens)) == 1 and lens[0] >= 3
        if _uniform_chains:
            P = np.stack([np.asarray(c) for c in chains], axis=0)
            F = P.shape[0]; m = P.shape[1] - 1
    if _uniform_chains:
        b = _min_image_orthorhombic(pos[P[:, 1:]] - pos[P[:, :-1]], box)  # (F,m,3)
        M = inv_gamma[P]                                         # (F, m+1)
        b2 = np.einsum("fab,fab->fa", b, b)
        G = np.zeros((F, m, m))
        for a in range(m):
            G[:, a, a] = 4.0 * b2[:, a] * (M[:, a] + M[:, a + 1])
        bdot = np.einsum("fab,fab->fa", b[:, :-1], b[:, 1:])     # (F, m-1)
        for a in range(m - 1):
            off = -4.0 * M[:, a + 1] * bdot[:, a]
            G[:, a, a + 1] = off; G[:, a + 1, a] = off
        sign, logdet = np.linalg.slogdet(G)                     # (F,)
        if np.any(sign <= 0):
            raise FloatingPointError("Fixman metric det G non-positive (vectorised).")
        U_F = FIXMAN_SIGN * half_kT * float(logdet.sum())
        Ginv = np.linalg.inv(G)                                 # (F, m, m)
        dlogdet_db = np.zeros((F, m, 3))
        for a in range(m):
            dlogdet_db[:, a, :] += (Ginv[:, a, a, None] * 8.0
                                    * (M[:, a] + M[:, a + 1])[:, None] * b[:, a, :])
            if a >= 1:
                dlogdet_db[:, a, :] += (2.0 * Ginv[:, a - 1, a, None]
                                        * (-4.0 * M[:, a, None]) * b[:, a - 1, :])
            if a <= m - 2:
                dlogdet_db[:, a, :] += (2.0 * Ginv[:, a, a + 1, None]
                                        * (-4.0 * M[:, a + 1, None]) * b[:, a + 1, :])
        grad = FIXMAN_SIGN * half_kT * dlogdet_db               # (F, m, 3)
        bf = np.zeros((F, m + 1, 3))
        bf[:, :-1, :] += grad
        bf[:, 1:, :] -= grad
        force[P] += bf                                          # chains disjoint
        return U_F, force

    # ---- Ragged fallback: per-chain (mixed lengths) ----
    for chain in chains:
        p = np.asarray(chain)
        m = p.shape[0] - 1            # number of bonds
        if m < 2:
            continue                  # det G config-independent → no force
        # Bond vectors (minimum-image).
        b = _min_image_orthorhombic(pos[p[1:]] - pos[p[:-1]], box)   # (m, 3)
        M = inv_gamma[p]                                             # (m+1,)
        # Build G (m, m) tridiagonal.
        G = np.zeros((m, m), dtype=np.float64)
        b2 = np.einsum("ab,ab->a", b, b)
        for a in range(m):
            G[a, a] = 4.0 * b2[a] * (M[a] + M[a + 1])
        bdot = np.einsum("ab,ab->a", b[:-1], b[1:])                  # (m-1,)
        for a in range(m - 1):
            off = -4.0 * M[a + 1] * bdot[a]
            G[a, a + 1] = off
            G[a + 1, a] = off
        sign, logdet = np.linalg.slogdet(G)
        if sign <= 0:
            raise FloatingPointError(
                f"Fixman metric det G non-positive (sign={sign}); chain {p}."
            )
        U_F += FIXMAN_SIGN * half_kT * logdet
        Ginv = np.linalg.inv(G)
        # d ln det G / d b_a  (3-vector), assembled from entries touching b_a.
        dlogdet_db = np.zeros((m, 3), dtype=np.float64)
        for a in range(m):
            # Diagonal G[a,a] = 4 b2[a] (M_a+M_{a+1}) → ∂/∂b_a = 8(M_a+M_{a+1}) b_a
            dlogdet_db[a] += Ginv[a, a] * 8.0 * (M[a] + M[a + 1]) * b[a]
            # Coupling to previous bond: entries (a-1,a),(a,a-1) = −4 M_a (b_{a-1}·b_a)
            if a >= 1:
                dlogdet_db[a] += 2.0 * Ginv[a - 1, a] * (-4.0 * M[a]) * b[a - 1]
            # Coupling to next bond: entries (a,a+1),(a+1,a) = −4 M_{a+1}(b_a·b_{a+1})
            if a <= m - 2:
                dlogdet_db[a] += 2.0 * Ginv[a, a + 1] * (-4.0 * M[a + 1]) * b[a + 1]
        # Chain rule b_a = r_{p_{a+1}} − r_{p_a}:  ∂/∂r_{p_a} = −, ∂/∂r_{p_{a+1}} = +.
        # F_{p_i} = −FIXMAN_SIGN·½kT · d ln det G / d r_{p_i}
        for a in range(m):
            grad = FIXMAN_SIGN * half_kT * dlogdet_db[a]
            force[p[a]] += grad        # from −(∂/∂r_{p_a}) = +dlogdet_db[a]
            force[p[a + 1]] -= grad     # from −(∂/∂r_{p_{a+1}}) = −dlogdet_db[a]
    return U_F, force


# ---------------------------------------------------------------------------
# The constrained L-M Action
# ---------------------------------------------------------------------------
class ConstrainedLeimkuhlerMatthewsBAOAB(hoomd.custom.Action):
    """L-M BAOAB-limit step with rigid backbone bonds (SHAKE) + Fixman.

    Mirrors :class:`ffn_sim.integrator.baoab.LeimkuhlerMatthewsBAOAB` (same
    predictor, prv_rnds memory, tag-indexed buffers, RUNTIME guards) and
    adds the Fixman pseudo-force + SHAKE projection. The stretch-bond
    Harmonic force MUST be omitted from the ``md.Integrator.forces`` when
    this Action is used (the bond is now the constraint); angle / LJ /
    ERM / xlink / myosin forces stay and are read from ``net_force``.

    Parameters
    ----------
    kT, gamma, dt, seed : as in the unconstrained Action.
    constraint_pairs : (C, 2) int — bead TAGS of each rigid bond.
    constraint_lengths : (C,) float — rigid bond length ℓ_a (SI [m]).
    chains : sequence of int arrays — ordered bead TAGS per filament,
        used to build the per-filament Fixman metric. Pass ``None`` or
        empty for a pure-constraint run with no Fixman correction (e.g.
        an all-dimer system).
    shake_tol, shake_max_iter : SHAKE convergence controls.
    """

    def __init__(
        self,
        *,
        kT: float,
        gamma: Mapping[str, float],
        dt: float,
        constraint_pairs: np.ndarray,
        constraint_lengths: np.ndarray,
        chains: Sequence[np.ndarray] | None = None,
        seed: int = 0,
        shake_tol: float = 1.0e-10,
        shake_max_iter: int = 500,
        record_lambda: bool = False,
    ) -> None:
        super().__init__()
        if not (np.isfinite(kT) and kT >= 0.0):
            raise ValueError(f"kT must be finite and ≥ 0, got {kT!r}")
        if not (np.isfinite(dt) and dt > 0.0):
            raise ValueError(f"dt must be finite and > 0, got {dt!r}")
        for typ, g in gamma.items():
            if not (np.isfinite(g) and g > 0.0):
                raise ValueError(f"gamma['{typ}'] must be finite and > 0, got {g!r}")
        cp = np.asarray(constraint_pairs, dtype=np.int64)
        cl = np.asarray(constraint_lengths, dtype=np.float64)
        if cp.ndim != 2 or cp.shape[1] != 2:
            raise ValueError(f"constraint_pairs must be (C, 2); got {cp.shape}")
        if cl.shape != (cp.shape[0],):
            raise ValueError(
                f"constraint_lengths must be (C,)={cp.shape[0]}; got {cl.shape}"
            )
        if cp.shape[0] and not np.all(np.isfinite(cl)) or np.any(cl <= 0.0):
            raise ValueError("constraint_lengths must all be finite and > 0.")

        self.kT = float(kT)
        self.dt = float(dt)
        self.gamma_map = dict(gamma)
        self._rng = np.random.default_rng(seed)
        self._constraint_pairs_tag = cp
        self._constraint_lengths = cl
        self._chains_tag = (
            [np.asarray(c, dtype=np.int64) for c in chains] if chains else []
        )
        # Pre-stack chain TAGs to a (F, N) array (chains don't change shape
        # during a run). Each step's chains_row is then ONE numpy index op,
        # not a Python list comprehension + stack — saves ~10% per the profile
        # (np.stack was the #2 hotspot at 142 ms / 1000 steps).
        self._chains_tag_stacked: np.ndarray | None = None
        if self._chains_tag:
            lens = {c.shape[0] for c in self._chains_tag}
            if len(lens) == 1:
                self._chains_tag_stacked = np.stack(self._chains_tag, axis=0)
        # Uniform bond length enables the fast per-chain tridiagonal M-SHAKE;
        # mixed lengths fall back to Gauss-Seidel shake_project.
        self._chain_rest_length: float | None = (
            float(cl[0]) if cl.shape[0] and np.allclose(cl, cl[0], rtol=1e-12, atol=0)
            else None
        )
        self.shake_tol = float(shake_tol)
        self.shake_max_iter = int(shake_max_iter)

        self._gamma_by_tag: np.ndarray | None = None
        self._inv_gamma_by_tag: np.ndarray | None = None
        self._bd_prefactor_by_tag: np.ndarray | None = None
        self._prv_rnds: np.ndarray | None = None
        self._sim_ref: hoomd.Simulation | None = None
        self._steps_run = 0
        self._max_drift = 0.0
        # Rigid-bond Lagrange-multiplier exposure (R1 in
        # RIGID_LAGRANGE_TENSION_DESIGN.md, PI verbal 2026-05-28). When
        # ``record_lambda=True`` the M-SHAKE-converged per-bond Lagrange
        # multiplier vector is captured into ``self._lambda_buf`` each step
        # so consumers (KU-3.5 tension method-of-planes) can include the
        # rigid-bond shell tension that is otherwise invisible to soft-bond
        # summation. ``record_lambda=False`` → zero overhead (the
        # accumulator isn't allocated inside ``shake_project_chains``).
        self.record_lambda = bool(record_lambda)
        self._lambda_buf = None  # most-recent step's λ; shape depends on chains

    # ------------------------------------------------------------------
    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation
        ig = simulation.operations.integrator
        if ig is None:
            raise RuntimeError(
                "ConstrainedLeimkuhlerMatthewsBAOAB requires an md.Integrator "
                "(forces attached, methods=[]). None is set."
            )
        if not np.isclose(float(ig.dt), self.dt, rtol=0, atol=0):
            raise RuntimeError(
                f"Integrator dt={float(ig.dt)!r} must equal Action dt={self.dt!r}."
            )
        if len(ig.methods) != 0:
            raise RuntimeError(
                "md.Integrator.methods must be empty; got "
                f"{list(ig.methods)!r} (a Method would double-step)."
            )

        with simulation.state.cpu_local_snapshot as snap:
            type_ids = np.asarray(snap.particles.typeid)
            tags = np.asarray(snap.particles.tag)
            type_names = list(simulation.state.particle_types)
            N = type_ids.shape[0]
            if N == 0:
                raise RuntimeError("Constrained Action attached to empty state.")
            missing = [t for t in type_names if t not in self.gamma_map]
            if missing:
                raise RuntimeError(f"gamma missing types {missing}; have {type_names}.")
            # Single-rank dense-tag assumption (same as baoab.py; constraints
            # spanning a domain boundary would break SHAKE under MPI).
            if int(tags.max()) >= N or int(tags.min()) < 0:
                raise RuntimeError(
                    "Constrained Action requires dense tags in [0, N); got "
                    f"[{int(tags.min())}, {int(tags.max())}] for N={N}."
                )
            gamma_by_typeid = np.array(
                [self.gamma_map[t] for t in type_names], dtype=np.float64
            )
            gamma_by_tag = np.empty(N, dtype=np.float64)
            gamma_by_tag[tags] = gamma_by_typeid[type_ids]

        self._gamma_by_tag = gamma_by_tag
        self._inv_gamma_by_tag = 1.0 / gamma_by_tag
        self._bd_prefactor_by_tag = np.sqrt(
            self.kT / (2.0 * gamma_by_tag * self.dt)
        ).reshape(-1, 1)
        self._prv_rnds = np.zeros((N, 3), dtype=np.float64)
        self._steps_run = 0
        # Constrained beads must share a common mobility for the simplified
        # equal-γ SHAKE/Fixman metric (backbone actin beads do). Verify.
        if self._constraint_pairs_tag.shape[0]:
            cg = gamma_by_tag[self._constraint_pairs_tag.ravel()]
            if not np.allclose(cg, cg[0], rtol=1e-12, atol=0):
                raise RuntimeError(
                    "Constrained beads must share one γ (backbone actin); "
                    "mixed-γ rigid bonds are out of scope for this Action."
                )

    # ------------------------------------------------------------------
    def act(self, timestep: int) -> None:
        if self._prv_rnds is None or self._inv_gamma_by_tag is None:
            raise RuntimeError("act() called before attach().")
        sim = self._sim_ref
        assert sim is not None
        with sim.state.cpu_local_snapshot as snap:
            pos = np.asarray(snap.particles.position)
            F = np.asarray(snap.particles.net_force)
            image = np.asarray(snap.particles.image)
            tag = np.asarray(snap.particles.tag)
            N = pos.shape[0]
            if N != self._prv_rnds.shape[0]:
                raise RuntimeError(
                    f"Particle count changed ({self._prv_rnds.shape[0]}→{N})."
                )
            if not np.all(np.isfinite(F)):
                bad = np.argwhere(~np.isfinite(F))
                raise FloatingPointError(
                    f"Non-finite net_force at timestep={timestep}; "
                    f"first: {bad[:5].tolist()}."
                )

            box = sim.state.box
            # Cache box dimensions ONCE per act() — avoids the dominant
            # per-call getattr overhead seen in cProfile (box.Lx/Ly/Lz +
            # box.L + _vec3_to_array totalled ~30% of step time).
            box_L = _box_L(box)
            inv_gamma_row = self._inv_gamma_by_tag[tag].reshape(-1, 1)
            bd_prefactor_row = self._bd_prefactor_by_tag[tag]
            prv_W_row = self._prv_rnds[tag]
            W_row = self._rng.standard_normal(size=(N, 3))

            # tag → row map so constraint/chain TAG lists index current rows.
            row_of_tag = np.empty(N, dtype=np.int64)
            row_of_tag[tag] = np.arange(N, dtype=np.int64)

            # 1. Fixman pseudo-force (added to net_force before predictor).
            F_total = F
            if self._chains_tag:
                # Vectorised: stacked tag array → one indexing op (avoids
                # Python list comprehension + per-step np.stack hotspot).
                if self._chains_tag_stacked is not None:
                    chains_row_arg = row_of_tag[self._chains_tag_stacked]
                else:
                    chains_row_arg = [row_of_tag[c] for c in self._chains_tag]
                _, F_fixman = fixman_logdet_and_force(
                    pos, chains_row_arg, self.kT, self._inv_gamma_by_tag[tag], box_L
                )
                F_total = F + F_fixman

            # 2. L-M predictor.
            dr = F_total * inv_gamma_row * self.dt + (
                bd_prefactor_row * (W_row + prv_W_row) * self.dt
            )
            pred = pos + dr

            # 3. SHAKE projection onto the rigid-bond manifold.
            if self._constraint_pairs_tag.shape[0]:
                pairs_row = row_of_tag[self._constraint_pairs_tag]
                inv_g = self._inv_gamma_by_tag[tag]
                if self._chains_tag and self._chain_rest_length is not None:
                    # Fast tridiagonal M-SHAKE per chain (uniform bond length).
                    if self._chains_tag_stacked is not None:
                        chains_row_arg = row_of_tag[self._chains_tag_stacked]
                    else:
                        chains_row_arg = [row_of_tag[c] for c in self._chains_tag]
                    if self.record_lambda:
                        projected, self._lambda_buf = shake_project_chains(
                            pred, pos, chains_row_arg,
                            self._chain_rest_length, inv_g, box_L,
                            tol=self.shake_tol, max_iter=self.shake_max_iter,
                            return_lambdas=True,
                        )
                    else:
                        projected = shake_project_chains(
                            pred, pos, chains_row_arg,
                            self._chain_rest_length, inv_g, box_L,
                            tol=self.shake_tol, max_iter=self.shake_max_iter,
                        )
                else:
                    projected = shake_project(
                        pred, pos, pairs_row, self._constraint_lengths,
                        inv_g, box_L,
                        tol=self.shake_tol, max_iter=self.shake_max_iter,
                    )
                # §4 drift guard.
                s = _min_image_orthorhombic(
                    projected[pairs_row[:, 0]] - projected[pairs_row[:, 1]], box_L
                )
                drift = np.abs(
                    np.linalg.norm(s, axis=1) - self._constraint_lengths
                ) / self._constraint_lengths
                self._max_drift = float(drift.max())
            else:
                projected = pred

            if not np.all(np.isfinite(projected)):
                bad = np.argwhere(~np.isfinite(projected))
                raise FloatingPointError(
                    f"Non-finite position after SHAKE at timestep={timestep}; "
                    f"first: {bad[:5].tolist()}."
                )

            # 4. Wrap.
            wrapped, img_delta = _wrap_into_box(projected, box)
            pos[:] = wrapped
            image[:] = image + img_delta
            self._prv_rnds[tag] = W_row
        self._steps_run += 1

    # ------------------------------------------------------------------
    @property
    def prv_rnds(self) -> np.ndarray | None:
        if self._prv_rnds is None:
            return None
        v = self._prv_rnds.view()
        v.flags.writeable = False
        return v

    @property
    def steps_run(self) -> int:
        return self._steps_run

    @property
    def max_constraint_drift(self) -> float:
        """Largest relative bond-length drift after the last SHAKE."""
        return self._max_drift

    @property
    def lambda_buf(self):
        """Most-recent step's accumulated per-bond Lagrange-multiplier vector.

        Shape (F, m) ndarray for the uniform-chain fast path (cortex case);
        list of (m_chain,) ndarrays for the ragged fallback. ``None`` if
        ``record_lambda=False`` or before any ``act()`` call. Per R1 in
        RIGID_LAGRANGE_TENSION_DESIGN.md: consumers convert each
        ``λ`` element to physical scalar bond tension via
        ``T_bond = λ · r₀ / Δt`` and project onto the cut normal for
        method-of-planes γ.
        """
        return self._lambda_buf


def make_constrained_baoab_updater(
    *,
    kT: float,
    gamma: Mapping[str, float],
    dt: float,
    constraint_pairs: np.ndarray,
    constraint_lengths: np.ndarray,
    chains: Sequence[np.ndarray] | None = None,
    seed: int = 0,
    shake_tol: float = 1.0e-10,
    shake_max_iter: int = 500,
    record_lambda: bool = False,
) -> tuple[ConstrainedLeimkuhlerMatthewsBAOAB, hoomd.update.CustomUpdater]:
    """Build the constrained L-M Action wrapped in a per-step CustomUpdater."""
    action = ConstrainedLeimkuhlerMatthewsBAOAB(
        kT=kT, gamma=gamma, dt=dt,
        constraint_pairs=constraint_pairs, constraint_lengths=constraint_lengths,
        chains=chains, seed=seed, shake_tol=shake_tol, shake_max_iter=shake_max_iter,
        record_lambda=record_lambda,
    )
    updater = hoomd.update.CustomUpdater(
        action=action, trigger=hoomd.trigger.Periodic(1)
    )
    return action, updater
