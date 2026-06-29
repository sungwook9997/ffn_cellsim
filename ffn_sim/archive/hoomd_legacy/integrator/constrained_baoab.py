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

from ffn_sim.archive.hoomd_legacy.integrator.baoab import _wrap_into_box

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
# Device backend dispatch (GPU-main port 2026-06-01)
# ---------------------------------------------------------------------------
def array_backend(use_gpu: bool):
    """Return the array module for the constrained Action (cupy or numpy).

    Centralises the ``xp = cp if gpu else np`` dispatch. The cupy import is
    deferred so the CPU/dev path (and every existing caller, which defaults
    to ``xp=np``) never requires cupy to be installed. On GPU hosts (gbook,
    A5000) ``use_gpu=True`` returns the cupy module; the pure constrained
    functions then run device-resident over ``gpu_local_snapshot`` arrays.
    """
    if use_gpu:
        import cupy as cp  # deferred: only needed on GPU hosts
        return cp
    return np


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def _min_image_orthorhombic(dr, box, xp=np):
    """Minimum-image displacement for an orthorhombic (cube) box.

    Accepts EITHER a HOOMD box object (legacy) or a precomputed (3,)
    array of ``[Lx, Ly, Lz]`` (preferred — avoids repeated HOOMD-property
    getattr that dominated profile traces). Cortex / single-filament boxes
    are cubes with no tilt; rigid bonds are never used with Lees-Edwards
    shear, so orthorhombic is exact. ℓ₀ ≪ L guarantees the nearest image
    is the physical bond.

    Device dispatch (GPU-main port 2026-06-01): ``xp`` selects the array
    backend (``np`` default → bit-identical to the pre-port CPU path; ``cp``
    for cupy/GPU). When ``box`` is a precomputed array the caller is
    responsible for placing it on the same device as ``dr`` (the GPU ``act``
    builds ``box_L`` as a cupy array once per step). The HOOMD-box legacy
    path stays host/numpy.
    """
    if hasattr(box, "Lx"):
        L = np.array([box.Lx, box.Ly, box.Lz], dtype=np.float64)
    else:
        L = box  # precomputed (3,) array, already on the xp device
    return dr - L * xp.round(dr / L)


def _box_L(box) -> np.ndarray:
    """Cache HOOMD box → (3,) numpy array once per act() to avoid
    per-call getattr overhead (the dominant cost in the profile)."""
    return np.array([box.Lx, box.Ly, box.Lz], dtype=np.float64)


def _wrap_into_box_xp(pos, box, xp):
    """Device-aware minimum-image box wrap (GPU-main port 2026-06-01).

    Sibling of the FROZEN ``baoab._wrap_into_box`` (CLAUDE.md
    §integrator-freeze — that module is untouched; the CPU constrained path
    still calls it for bit-identical behaviour). Same upper-triangular
    fractional-coordinate algorithm + the int32-image overflow guard, but
    every array op dispatches through ``xp`` so it runs on cupy
    (``gpu_local_snapshot`` positions stay device-resident). Box dimensions
    are HOST scalars (reading ``box.Lx`` etc. does not sync device arrays);
    ``pos`` is the only device array.

    Returns ``(wrapped_pos, image_delta)`` as ``xp`` arrays (float64 / int32),
    matching the frozen reference's contract.
    """
    Lx, Ly, Lz = box.Lx, box.Ly, box.Lz
    xy, xz, yz = box.xy, box.xz, box.yz

    rx = pos[:, 0]; ry = pos[:, 1]; rz = pos[:, 2]
    fz = rz / Lz
    fy = (ry - yz * Lz * fz) / Ly
    fx = (rx - xy * Ly * fy - xz * Lz * fz) / Lx

    nx = xp.round(fx); ny = xp.round(fy); nz = xp.round(fz)

    # §4 int32-image overflow guard (mirrors baoab._wrap_into_box). The three
    # ``.any()`` reductions are 0-d device→host scalars (a tiny, necessary
    # safety sync — not a full-array transfer).
    INT32_GUARD = 1.0e8
    if (
        bool((xp.abs(nx) > INT32_GUARD).any())
        or bool((xp.abs(ny) > INT32_GUARD).any())
        or bool((xp.abs(nz) > INT32_GUARD).any())
    ):
        worst = float(
            max(float(xp.abs(nx).max()), float(xp.abs(ny).max()),
                float(xp.abs(nz).max()))
        )
        raise FloatingPointError(
            "Constrained BAOAB _wrap_into_box_xp: |fractional coord| exceeded "
            f"the int32-image guard (worst |round(f)|={worst:.3e} > "
            f"{INT32_GUARD:.0e}); see baoab._wrap_into_box for the rationale."
        )

    fx = fx - nx; fy = fy - ny; fz = fz - nz

    out = xp.empty_like(pos)
    out[:, 0] = Lx * fx + xy * Ly * fy + xz * Lz * fz
    out[:, 1] = Ly * fy + yz * Lz * fz
    out[:, 2] = Lz * fz

    img_delta = xp.empty_like(pos, dtype=np.int32)
    img_delta[:, 0] = nx.astype(np.int32)
    img_delta[:, 1] = ny.astype(np.int32)
    img_delta[:, 2] = nz.astype(np.int32)
    return out, img_delta


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
            rhs: np.ndarray, xp=np) -> np.ndarray:
    """Solve a tridiagonal system (Thomas algorithm). sub[0], sup[-1] unused.

    ``xp`` is the array backend (np default → bit-identical CPU path; cp for
    cupy). The forward/backward recurrence is sequential over the chain
    length so each ``k`` launches one small ``xp`` op — cheap (m≈6).
    """
    n = diag.shape[0]
    cc = xp.empty(n); dp = xp.empty(n)
    cc[0] = sup[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for k in range(1, n):
        m = diag[k] - sub[k] * cc[k - 1]
        cc[k] = sup[k] / m
        dp[k] = (rhs[k] - sub[k] * dp[k - 1]) / m
    x = xp.empty(n)
    x[-1] = dp[-1]
    for k in range(n - 2, -1, -1):
        x[k] = dp[k] - cc[k] * x[k + 1]
    return x


def _thomas_batched(sub: np.ndarray, diag: np.ndarray, sup: np.ndarray,
                    rhs: np.ndarray, xp=np) -> np.ndarray:
    """Batched tridiagonal solve over the leading axis. All (F, m).

    ``xp`` is the array backend (np default → bit-identical CPU path; cp for
    cupy). The ``k``-loop over chain length m (~6) is a sequential recurrence
    that stays a Python loop; each iteration is one batched ``xp`` op over the
    F chains, so it launches m cupy kernels (cheap vs the per-step host sync
    it removes).
    """
    F, m = diag.shape
    cc = xp.empty((F, m)); dp = xp.empty((F, m))
    cc[:, 0] = sup[:, 0] / diag[:, 0]
    dp[:, 0] = rhs[:, 0] / diag[:, 0]
    for k in range(1, m):
        den = diag[:, k] - sub[:, k] * cc[:, k - 1]
        cc[:, k] = sup[:, k] / den
        dp[:, k] = (rhs[:, k] - sub[:, k] * dp[:, k - 1]) / den
    x = xp.empty((F, m))
    x[:, -1] = dp[:, -1]
    for k in range(m - 2, -1, -1):
        x[:, k] = dp[:, k] - cc[:, k] * x[:, k + 1]
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
    compression_release: bool = False,
    release_load_crit: float | None = None,
    dt: float | None = None,
    eligible_mask=None,
    xp=np,
):
    """Matrix-SHAKE (tridiagonal Newton) for linear-chain bond constraints.

    Gate-B relaxed mode (``compression_release=True``, 2026-06-08): the rigid
    equality bond constraint ``|s_a| = ℓ₀`` becomes the **unilateral** (one-sided)
    constraint ``|s_a| ≤ ℓ₀`` — the projection onto the closed ball. A bond is
    enforced only while it is in **tension** (stretched, ``g_a = |s|²−ℓ² ≥ 0``,
    pulled back to ℓ₀ — the physiological inextensible actin axis, k_axial≈154
    pN/nm); a bond in **compression** (``g_a < 0``) is **released** (λ_a = 0, free
    to shorten), permitting the filament to buckle / the network to condense
    (the symmetry-breaking a fixed-length rod suppresses). This is the H.7 Gate-B
    transmission lever (contract §3): tension side inextensible, compression side
    yields. ``compression_release=False`` (default) is the exact rigid M-SHAKE —
    bit-identical to the pre-Gate-B path (control §5.1 superset parity).

    Active-set numerics: a released bond becomes an identity row of the
    tridiagonal Jacobian (λ_a = 0) AND its off-diagonal couplings to active
    neighbours are zeroed, so the chain cleanly decomposes into independent
    constrained runs at each released joint (a tridiagonal solve cannot otherwise
    drop a mid-chain bond — the shared-bead coupling would leak through). The
    active set is re-evaluated each Newton iteration from the current ``g``;
    convergence is ``max_a max(g_a, 0)/ℓ² ≤ tol`` (released bonds, g<0, are
    already feasible and excluded from the violation norm).

    Euler-thresholded release (``release_load_crit`` set, the physiological
    Gate-B mode, contract §8): a compressed bond is released ONLY if its
    *compressive constraint force* exceeds the critical buckling load
    ``F_crit`` — NOT at every sub-ℓ₀ thermal fluctuation. A first rigid pass
    yields the per-bond Lagrange multiplier ``λ_rigid`` (bond tension
    ``T = λ·r₀/Δt``; λ>0 tension, λ<0 compression); a bond is *buckle-eligible*
    iff ``λ_rigid < −F_crit·Δt/r₀`` (compressive load past the Euler threshold).
    The unilateral solve then releases a bond iff it is BOTH compressed (g<0)
    AND eligible — sub-threshold compression stays rigid (the rod holds), and
    the tension side is always enforced. This restores the passive control
    (few unloaded bonds exceed F_crit ⇒ relaxed-OFF ≈ rigid-OFF) and isolates
    myosin-loaded buckling. ``release_load_crit=None`` ⇒ the pure geometric
    unilateral release (every g<0 bond; used by the unit tests), a zero-
    threshold limit that over-condenses and is NOT the production mode.

    Eligibility timescale (PI 2026-06-08): the per-step rigid λ is thermal-noise
    dominated (~6–45 pN, ≫ F_crit even passively), so the INTERNAL
    ``release_load_crit`` path (instantaneous λ) over-releases and is for unit
    tests only. The PRODUCTION gate passes a precomputed ``eligible_mask`` (F,m
    bool) built by the Action from a τ_bend-windowed EMA of the per-bond
    *sustained* compressive load. When ``eligible_mask`` is given it overrides
    the internal computation.

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

    Device dispatch (GPU-main port 2026-06-01): ``xp`` is the array backend
    (np default → bit-identical CPU reference; cp for cupy/GPU). The GPU path
    requires the uniform fast path — pass a pre-stacked ``(F, m+1)`` cupy index
    array as ``chains`` with cupy ``pred_pos``/``ref_pos``/``inv_mass``; the
    ragged per-chain fallback stays numpy-only (mixed-length chains never
    occur on the GPU cortex path). The ``if drift <= tol`` break reads one
    0-d scalar back each iteration (a tiny, unavoidable sync for the
    convergence test).
    """
    pos = pred_pos.copy()
    L2 = rest_length * rest_length

    # ---- Vectorised fast path: all chains the same length (cortex N=7) ----
    # Accept EITHER a pre-stacked (F, N) array — numpy OR cupy (GPU path) — or
    # a list of per-chain arrays (legacy ragged path, numpy-only fallback
    # below). Detection is duck-typed (``ndim``) so a stacked cupy array takes
    # the fast path too; ``xp`` dispatches every array op (np default →
    # bit-identical CPU reference).
    _stacked = getattr(chains, "ndim", 0) == 2 and chains.shape[1] >= 3
    if _stacked:
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
        d0 = _min_image_orthorhombic(
            ref_pos[P[:, :-1]] - ref_pos[P[:, 1:]], box, xp=xp)
        # Euler-thresholded eligibility (contract §8): a first RIGID pass gives
        # λ_rigid; a bond may buckle only if its compressive constraint load
        # exceeds F_crit (λ_rigid < −F_crit·Δt/r₀). None ⇒ geometric (all
        # bonds eligible — the zero-threshold test/limit mode).
        eligible = None
        if compression_release and eligible_mask is not None:
            eligible = eligible_mask                             # (F, m) precomputed
        elif compression_release and release_load_crit is not None:
            if dt is None or dt <= 0.0:
                raise ValueError("release_load_crit requires dt > 0 (T=λ·r₀/Δt).")
            _, lam_rigid = shake_project_chains(
                pred_pos, ref_pos, P, rest_length, inv_mass, box,
                tol=tol, max_iter=max_iter, return_lambdas=True,
                compression_release=False, xp=xp)
            lam_crit = release_load_crit * dt / rest_length
            eligible = lam_rigid < -lam_crit                     # (F, m) bool
        # Accumulated per-bond Lagrange multiplier across Newton iterations.
        # Allocated only when the caller actually needs it so the default
        # path stays zero-allocation extra.
        lambda_total = xp.zeros((F, m)) if return_lambdas else None
        for _ in range(max_iter):
            s = _min_image_orthorhombic(
                pos[P[:, :-1]] - pos[P[:, 1:]], box, xp=xp)
            g = xp.einsum("fab,fab->fa", s, s) - L2             # (F, m)
            if compression_release:
                # Released set: compressed (g<0) AND buckle-eligible (Euler gate;
                # all g<0 in geometric mode). Enforced = everything else (tension
                # OR sub-threshold compression) and must reach g=0 → only those
                # count toward convergence.
                rel = g < 0.0                                   # (F, m) bool
                if eligible is not None:
                    rel = rel & eligible
                if xp.max(xp.where(rel, 0.0, xp.abs(g))) / L2 <= tol:
                    break
            elif xp.max(xp.abs(g)) / L2 <= tol:                 # 0-d sync (1 scalar)
                break
            sd = xp.einsum("fab,fab->fa", s, d0)
            diag = -2.0 * (M[:, :-1] + M[:, 1:]) * sd           # (F, m)
            sub = xp.zeros((F, m)); sup = xp.zeros((F, m))
            if m > 1:
                sub[:, 1:] = 2.0 * M[:, 1:-1] * xp.einsum(
                    "fab,fab->fa", s[:, 1:], d0[:, :-1])
                sup[:, :-1] = 2.0 * M[:, 1:-1] * xp.einsum(
                    "fab,fab->fa", s[:, :-1], d0[:, 1:])
            rhs = -g
            if compression_release:
                # Released bonds → identity row (λ=0); decouple active neighbours
                # from the released joint so the tridiagonal solve splits into
                # independent active runs (see docstring).
                diag = xp.where(rel, 1.0, diag)
                sub = xp.where(rel, 0.0, sub)
                sup = xp.where(rel, 0.0, sup)
                rhs = xp.where(rel, 0.0, rhs)
                if m > 1:
                    rel_next = xp.zeros((F, m), dtype=bool)
                    rel_prev = xp.zeros((F, m), dtype=bool)
                    rel_next[:, :-1] = rel[:, 1:]   # bond a's a+1 neighbour released
                    rel_prev[:, 1:] = rel[:, :-1]   # bond a's a-1 neighbour released
                    sup = xp.where(rel_next, 0.0, sup)
                    sub = xp.where(rel_prev, 0.0, sub)
            lam = _thomas_batched(sub, diag, sup, rhs, xp=xp)    # (F, m)
            if compression_release:
                lam = xp.where(rel, 0.0, lam)                   # released λ ≡ 0
            if lambda_total is not None:
                lambda_total += lam
            disp = xp.zeros((F, m + 1, 3))
            disp[:, :-1] -= (M[:, :-1] * lam)[:, :, None] * d0
            disp[:, 1:] += (M[:, 1:] * lam)[:, :, None] * d0
            pos[P] += disp                                       # chains disjoint
        else:
            _gd = xp.where(rel, 0.0, xp.abs(g)) if compression_release else xp.abs(g)
            drift = float(xp.max(_gd) / L2)
            raise RuntimeError(
                f"M-SHAKE (vectorised, {F} chains len {m + 1}) failed in "
                f"{max_iter} iters; max relative drift = {drift:.3e} > tol={tol}.")
        if return_lambdas:
            return pos, lambda_total
        return pos

    # ---- Ragged fallback: per-chain (mixed lengths) ----
    if compression_release and release_load_crit is not None:
        raise NotImplementedError(
            "Euler-thresholded release (release_load_crit) is implemented only on "
            "the uniform stacked fast path; the ragged mixed-length fallback "
            "supports geometric release (release_load_crit=None) only. The cortex "
            "backbone chains are uniform (7 beads), so production takes the fast path."
        )
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
            if compression_release:
                if np.max(np.where(g > 0.0, g, 0.0)) / L2 <= tol:
                    ok = True
                    break
            elif np.max(np.abs(g)) / L2 <= tol:
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
            rhs = -g
            if compression_release:
                # Unilateral: release compressed bonds (λ=0) + decouple their
                # active neighbours (identity-row chain decomposition).
                rel = g < 0.0
                diag = np.where(rel, 1.0, diag)
                sub = np.where(rel, 0.0, sub)
                sup = np.where(rel, 0.0, sup)
                rhs = np.where(rel, 0.0, rhs)
                if m > 1:
                    rel_next = np.zeros(m, dtype=bool); rel_next[:-1] = rel[1:]
                    rel_prev = np.zeros(m, dtype=bool); rel_prev[1:] = rel[:-1]
                    sup = np.where(rel_next, 0.0, sup)
                    sub = np.where(rel_prev, 0.0, sub)
            lam = _thomas(sub, diag, sup, rhs)
            if compression_release:
                lam = np.where(g < 0.0, 0.0, lam)
            if lambda_total is not None:
                lambda_total += lam
            # apply Δr[p_k] = M_k (λ_{k-1} d0_{k-1} − λ_k d0_k)
            disp = np.zeros((m + 1, 3))
            disp[:-1] -= (M[:-1] * lam)[:, None] * d0      # −λ_a d0_a on bead p_a
            disp[1:] += (M[1:] * lam)[:, None] * d0        # +λ_a d0_a on bead p_{a+1}
            pos[p] += disp
        if not ok:
            s = _min_image_orthorhombic(pos[p[:-1]] - pos[p[1:]], box)
            _g = np.einsum("ab,ab->a", s, s) - L2
            _gd = np.where(_g > 0.0, _g, 0.0) if compression_release else np.abs(_g)
            drift = float(np.max(_gd) / L2)
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
    xp=np,
) -> tuple[float, np.ndarray]:
    """Fixman pseudo-potential ``U_F`` and force ``F_F = −∇U_F``.

    For each chain (ordered bead row-indices ``[p0..pm]``, m rigid bonds)
    builds the mobility-metric Gram matrix ``G`` (m×m, tridiagonal):

        G_aa     =  4 |b_a|² (M_{p_a} + M_{p_{a+1}})
        G_{a,a+1}= −4 M_{p_{a+1}} (b_a · b_{a+1})

    and accumulates ``U_F = FIXMAN_SIGN · ½ kT Σ_chains ln det G`` and its
    analytic gradient. Chains with < 2 bonds have a configuration-
    independent ``det G`` (bond lengths are fixed) → zero force.

    Device dispatch (GPU-main port 2026-06-01): ``xp`` is the array backend
    (np default → bit-identical CPU reference; cp for cupy/GPU). The batched
    ``slogdet``/``inv`` use ``xp.linalg`` (present in cupy 14). The GPU path
    needs the uniform fast path (pre-stacked ``(F, m+1)`` cupy index array +
    cupy ``pos``/``inv_gamma``); the ragged per-chain fallback stays
    numpy-only. ``U_F`` is returned as a python float (one 0-d sync); ``force``
    is an ``xp`` array.

    Returns
    -------
    (U_F, force) with force shape (N, 3).
    """
    N = pos.shape[0]
    force = xp.zeros((N, 3), dtype=np.float64)
    U_F = 0.0
    half_kT = 0.5 * kT

    # ---- Vectorised fast path: all chains the same length (cortex N=7) ----
    # Accept pre-stacked (F, N) array — numpy OR cupy (GPU path) — (avoids
    # per-step np.stack, the #2 profile hotspot) OR a list of per-chain arrays
    # (legacy ragged, numpy-only). Detection is duck-typed (``ndim``).
    _is_stacked = getattr(chains, "ndim", 0) == 2 and chains.shape[1] >= 3
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
        b = _min_image_orthorhombic(
            pos[P[:, 1:]] - pos[P[:, :-1]], box, xp=xp)  # (F,m,3)
        M = inv_gamma[P]                                         # (F, m+1)
        b2 = xp.einsum("fab,fab->fa", b, b)
        G = xp.zeros((F, m, m))
        for a in range(m):
            G[:, a, a] = 4.0 * b2[:, a] * (M[:, a] + M[:, a + 1])
        bdot = xp.einsum("fab,fab->fa", b[:, :-1], b[:, 1:])     # (F, m-1)
        for a in range(m - 1):
            off = -4.0 * M[:, a + 1] * bdot[:, a]
            G[:, a, a + 1] = off; G[:, a + 1, a] = off
        sign, logdet = xp.linalg.slogdet(G)                     # (F,)
        if bool(xp.any(sign <= 0)):
            raise FloatingPointError("Fixman metric det G non-positive (vectorised).")
        U_F = FIXMAN_SIGN * half_kT * float(logdet.sum())
        Ginv = xp.linalg.inv(G)                                 # (F, m, m)
        dlogdet_db = xp.zeros((F, m, 3))
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
        bf = xp.zeros((F, m + 1, 3))
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

    Mirrors :class:`ffn_sim.archive.hoomd_legacy.integrator.baoab.LeimkuhlerMatthewsBAOAB` (same
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
        compression_release: bool = False,
        release_load_crit: float | None = None,
        load_tau: float | None = None,
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
        self._seed = int(seed)
        # Device backend, resolved in attach() from sim.device. ``xp=np`` /
        # CPU snapshot until then; on a GPU device attach() flips these to
        # cupy + gpu_local_snapshot (GPU-main port 2026-06-01). The CPU path
        # is bit-identical to the pre-port Action.
        self._xp = np
        self._on_gpu = False
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

        # H.7 Gate-B relaxed-constraint mode (2026-06-08). When True the rigid
        # equality bond constraint becomes unilateral |s|≤ℓ₀ (compression side
        # released → buckling/condensation permitted; tension side inextensible).
        # False = exact rigid M-SHAKE (control §5.1 superset parity). The
        # generic Gauss-Seidel ``shake_project`` (pairs path) does NOT support
        # it — chains are required (attach() / act() enforce this).
        self.compression_release = bool(compression_release)
        # Euler buckling load F_crit [N] — a compressed backbone bond is released
        # only if its compressive constraint force exceeds this (contract §8;
        # derived π²κ/L², NOT a free knob). None ⇒ geometric zero-threshold
        # release (test/limit mode, over-condenses — not production).
        self.release_load_crit = (
            float(release_load_crit) if release_load_crit is not None else None
        )
        # τ_bend EMA window for the SUSTAINED compressive load (PI 2026-06-08):
        # the per-step rigid λ is thermal-noise dominated, so eligibility is gated
        # on an EMA of the per-bond load over τ_bend = γ_b·ℓ₀³/κ_B (the segment
        # bending-relaxation time — the timescale a segment actually buckles on).
        # When set with compression_release + release_load_crit, the Action runs a
        # rigid pre-pass each step, EMAs the load, and passes the eligible_mask.
        self.load_tau = float(load_tau) if load_tau is not None else None
        self._load_ema = None   # (F, m) sustained per-bond constraint force [N]
        # Most-recent step's released (compression) bond fraction — the
        # condensation-engagement monitor (contract §5.3). 0.0 in rigid mode.
        self._last_released_frac: float = 0.0
        if self.compression_release and (
            not self._chains_tag or self._chain_rest_length is None
        ):
            raise ValueError(
                "compression_release=True requires uniform-rest-length chain "
                "constraints (the unilateral M-SHAKE lives in the chain path); "
                "the generic Gauss-Seidel pairs fallback does not support it."
            )

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

        # ---- Device backend (GPU-main port 2026-06-01) ----
        # Follow the simulation's device: a GPU device flips the per-step Action
        # onto gpu_local_snapshot + cupy so positions/forces stay device-resident
        # (removing the per-step host sync that throttled the C++ 71%). A CPU
        # device leaves everything on the bit-identical numpy path.
        self._on_gpu = isinstance(simulation.device, hoomd.device.GPU)
        if self._on_gpu:
            xp = array_backend(True)
            self._xp = xp
            # cupy RNG for the per-step Gaussian (own stream — gates are
            # statistical, re-passed not bit-matched vs the numpy stream).
            self._rng = xp.random.default_rng(self._seed)
            # GPU requires the vectorised uniform fast path (stacked chains);
            # the ragged per-chain + Gauss-Seidel shake_project fallbacks are
            # numpy-only. Cortex backbone chains are uniform → stacked.
            if self._chains_tag and self._chains_tag_stacked is None:
                raise RuntimeError(
                    "GPU constrained Action requires uniform-length chains "
                    "(stacked fast path); got ragged chains. Run on CPU or "
                    "pad/uniformise the backbone chains."
                )
            if (self._constraint_pairs_tag.shape[0]
                    and (not self._chains_tag or self._chain_rest_length is None)):
                raise RuntimeError(
                    "GPU constrained Action requires uniform-length chain "
                    "constraints (the cupy M-SHAKE fast path); the generic "
                    "shake_project fallback is CPU-only."
                )
            # Move the per-tag + topology buffers onto the device once.
            self._inv_gamma_by_tag = xp.asarray(self._inv_gamma_by_tag)
            self._bd_prefactor_by_tag = xp.asarray(self._bd_prefactor_by_tag)
            self._prv_rnds = xp.asarray(self._prv_rnds)
            self._constraint_pairs_tag = xp.asarray(self._constraint_pairs_tag)
            self._constraint_lengths = xp.asarray(self._constraint_lengths)
            if self._chains_tag_stacked is not None:
                self._chains_tag_stacked = xp.asarray(self._chains_tag_stacked)

    # ------------------------------------------------------------------
    def act(self, timestep: int) -> None:
        if self._prv_rnds is None or self._inv_gamma_by_tag is None:
            raise RuntimeError("act() called before attach().")
        sim = self._sim_ref
        assert sim is not None
        xp = self._xp
        # GPU-main port 2026-06-01: on a GPU device the whole step stays
        # device-resident via gpu_local_snapshot + cupy (xp=cp), removing the
        # per-step host sync. On CPU xp=np and this is the bit-identical
        # pre-port path (xp.asarray on a numpy view is a no-op).
        snap_ctx = (sim.state.gpu_local_snapshot if self._on_gpu
                    else sim.state.cpu_local_snapshot)
        with snap_ctx as snap:
            pos = xp.asarray(snap.particles.position)
            F = xp.asarray(snap.particles.net_force)
            image = xp.asarray(snap.particles.image)
            tag = xp.asarray(snap.particles.tag)
            N = pos.shape[0]
            if N != self._prv_rnds.shape[0]:
                raise RuntimeError(
                    f"Particle count changed ({self._prv_rnds.shape[0]}→{N})."
                )
            if not bool(xp.all(xp.isfinite(F))):
                bad = xp.argwhere(~xp.isfinite(F))[:5]
                if self._on_gpu:
                    bad = xp.asnumpy(bad)
                raise FloatingPointError(
                    f"Non-finite net_force at timestep={timestep}; "
                    f"first: {bad.tolist()}."
                )

            box = sim.state.box
            # Cache box dimensions ONCE per act() — avoids the dominant
            # per-call getattr overhead seen in cProfile (box.Lx/Ly/Lz +
            # box.L + _vec3_to_array totalled ~30% of step time). On GPU the
            # (3,) L array is moved to the device once (tiny — not a sync).
            box_L = xp.asarray(_box_L(box)) if self._on_gpu else _box_L(box)
            inv_gamma_row = self._inv_gamma_by_tag[tag].reshape(-1, 1)
            bd_prefactor_row = self._bd_prefactor_by_tag[tag]
            prv_W_row = self._prv_rnds[tag]
            W_row = self._rng.standard_normal(size=(N, 3))

            # tag → row map so constraint/chain TAG lists index current rows.
            row_of_tag = xp.empty(N, dtype=xp.int64)
            row_of_tag[tag] = xp.arange(N, dtype=xp.int64)

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
                    pos, chains_row_arg, self.kT,
                    self._inv_gamma_by_tag[tag], box_L, xp=xp,
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
                    # τ_bend-windowed buckle eligibility (PI 2026-06-08): a rigid
                    # pre-pass gives the per-bond constraint force; an EMA over
                    # τ_bend isolates the SUSTAINED compressive load from the
                    # thermal per-step impulse. A bond is buckle-eligible iff its
                    # sustained load is compressive past F_crit. Only on the
                    # stacked fast path (production cortex); requires F_crit+τ.
                    eligible_mask = None
                    if (self.compression_release
                            and self.release_load_crit is not None
                            and self.load_tau is not None
                            and self._chains_tag_stacked is not None):
                        _, lam_rigid = shake_project_chains(
                            pred, pos, chains_row_arg,
                            self._chain_rest_length, inv_g, box_L,
                            tol=self.shake_tol, max_iter=self.shake_max_iter,
                            return_lambdas=True, compression_release=False, xp=xp,
                        )
                        T_rigid = lam_rigid * (self._chain_rest_length / self.dt)
                        if self._load_ema is None:
                            self._load_ema = xp.zeros_like(T_rigid)
                        a_ema = self.dt / self.load_tau
                        self._load_ema = self._load_ema + a_ema * (T_rigid - self._load_ema)
                        eligible_mask = self._load_ema < -self.release_load_crit
                    if self.record_lambda:
                        projected, self._lambda_buf = shake_project_chains(
                            pred, pos, chains_row_arg,
                            self._chain_rest_length, inv_g, box_L,
                            tol=self.shake_tol, max_iter=self.shake_max_iter,
                            return_lambdas=True,
                            compression_release=self.compression_release,
                            release_load_crit=self.release_load_crit, dt=self.dt,
                            eligible_mask=eligible_mask, xp=xp,
                        )
                    else:
                        projected = shake_project_chains(
                            pred, pos, chains_row_arg,
                            self._chain_rest_length, inv_g, box_L,
                            tol=self.shake_tol, max_iter=self.shake_max_iter,
                            compression_release=self.compression_release,
                            release_load_crit=self.release_load_crit, dt=self.dt,
                            eligible_mask=eligible_mask, xp=xp,
                        )
                    if self.compression_release and self._chains_tag_stacked is not None:
                        # Condensation monitor (§5.3): fraction of backbone bonds
                        # left in compression (|s|<ℓ₀) after the unilateral solve.
                        cs = _min_image_orthorhombic(
                            projected[chains_row_arg[:, :-1]]
                            - projected[chains_row_arg[:, 1:]], box_L, xp=xp)
                        clen = xp.sqrt(xp.einsum("fab,fab->fa", cs, cs))
                        self._last_released_frac = float(
                            xp.mean((clen < self._chain_rest_length).astype(np.float64))
                        )
                else:
                    # Generic Gauss-Seidel fallback (CPU-only; attach() forbids
                    # this path on GPU — uniform chains are required there).
                    projected = shake_project(
                        pred, pos, pairs_row, self._constraint_lengths,
                        inv_g, box_L,
                        tol=self.shake_tol, max_iter=self.shake_max_iter,
                    )
                # §4 drift guard. In relaxed mode the constraint is unilateral
                # (|s|≤ℓ₀): released compression bonds legitimately sit below ℓ₀,
                # so the drift metric is the TENSION-side overshoot only
                # (max(‖s‖−ℓ₀, 0)/ℓ₀) — a compressed bond is feasible, not drift.
                s = _min_image_orthorhombic(
                    projected[pairs_row[:, 0]] - projected[pairs_row[:, 1]],
                    box_L, xp=xp,
                )
                _excess = xp.linalg.norm(s, axis=1) - self._constraint_lengths
                if self.compression_release:
                    _excess = xp.where(_excess > 0.0, _excess, 0.0)
                drift = xp.abs(_excess) / self._constraint_lengths
                self._max_drift = float(drift.max())
            else:
                projected = pred

            if not bool(xp.all(xp.isfinite(projected))):
                bad = xp.argwhere(~xp.isfinite(projected))[:5]
                if self._on_gpu:
                    bad = xp.asnumpy(bad)
                raise FloatingPointError(
                    f"Non-finite position after SHAKE at timestep={timestep}; "
                    f"first: {bad.tolist()}."
                )

            # 4. Wrap (device-aware sibling on GPU; frozen baoab helper on CPU).
            if self._on_gpu:
                wrapped, img_delta = _wrap_into_box_xp(projected, box, xp)
            else:
                wrapped, img_delta = _wrap_into_box(projected, box)
            pos[:] = wrapped
            image[:] = image + img_delta
            self._prv_rnds[tag] = W_row
        self._steps_run += 1

    # ------------------------------------------------------------------
    @property
    def prv_rnds(self) -> np.ndarray | None:
        """Read-only host (numpy) view of the previous-step Gaussian buffer.

        On a GPU device the buffer is cupy; it is copied to host here so the
        accessor always returns a numpy array (cupy has no ``.flags.writeable``,
        so the read-only view below would otherwise raise on GPU). Off the hot
        path — only used by tests / inspection."""
        if self._prv_rnds is None:
            return None
        if self._on_gpu:
            return self._xp.asnumpy(self._prv_rnds)
        v = self._prv_rnds.view()
        v.flags.writeable = False
        return v

    @property
    def steps_run(self) -> int:
        return self._steps_run

    @property
    def max_constraint_drift(self) -> float:
        """Largest relative bond-length drift after the last SHAKE.

        Rigid mode: |‖s‖−ℓ₀|/ℓ₀ over all bonds. Relaxed mode
        (``compression_release``): the tension-side overshoot only — released
        compression bonds (‖s‖<ℓ₀) are feasible, not drift."""
        return self._max_drift

    @property
    def released_fraction(self) -> float:
        """Fraction of backbone bonds left in compression (‖s‖<ℓ₀) after the
        last unilateral M-SHAKE — the condensation-engagement monitor
        (contract §5.3). 0.0 in rigid mode (``compression_release=False``)."""
        return self._last_released_frac

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

        On GPU the buffer is moved to host numpy here so downstream consumers
        keep the numpy contract — λ is read at most once per binding batch
        (~1% of steps), so the device→host copy is cheap and not on the hot
        per-step path.
        """
        buf = self._lambda_buf
        if buf is not None and self._on_gpu:
            buf = self._xp.asnumpy(buf)
        return buf

    @property
    def chains_tag_stacked(self):
        """Host (numpy) view of the stacked per-filament bead-TAG array.

        External consumers (e.g. the KU-3.5 method-of-planes rigid-bond tension
        in ``scripts/h3_ku35_tension.py``) index host positions with this, so
        on GPU the cupy hot-path buffer is copied to host here. Read off the
        per-step path (once per tension measurement), so the device→host copy
        is cheap. ``None`` if the Action carries no uniform chains. Use this
        rather than the private ``_chains_tag_stacked`` (which is device-resident
        cupy on a GPU device — the port moved it there for the hot path)."""
        c = self._chains_tag_stacked
        if c is not None and self._on_gpu:
            c = self._xp.asnumpy(c)
        return c


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
    compression_release: bool = False,
    release_load_crit: float | None = None,
    load_tau: float | None = None,
) -> tuple[ConstrainedLeimkuhlerMatthewsBAOAB, hoomd.update.CustomUpdater]:
    """Build the constrained L-M Action wrapped in a per-step CustomUpdater."""
    action = ConstrainedLeimkuhlerMatthewsBAOAB(
        kT=kT, gamma=gamma, dt=dt,
        constraint_pairs=constraint_pairs, constraint_lengths=constraint_lengths,
        chains=chains, seed=seed, shake_tol=shake_tol, shake_max_iter=shake_max_iter,
        record_lambda=record_lambda, compression_release=compression_release,
        release_load_crit=release_load_crit, load_tau=load_tau,
    )
    updater = hoomd.update.CustomUpdater(
        action=action, trigger=hoomd.trigger.Periodic(1)
    )
    return action, updater
