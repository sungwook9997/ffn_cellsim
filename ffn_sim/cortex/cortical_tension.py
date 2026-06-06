"""Unified 3-channel cortical-tension (γ) estimator — the GATE-B measurement tool.

This module is the single, self-contained tool a GATE-B sweep uses to read the
cortical tension γ of a built cell. It reports the tension through **three
physically distinct channels, kept SEPARATE**:

1. ``gamma_soft`` — *active / soft method-of-planes*. The harmonic-bond
   contribution: every soft bond (myosin head springs + head–actin attach,
   crosslinker intra + attach, ERM, any soft cortex component) carries a
   scalar tension ``T = k·(L − r₀)`` which is projected onto a set of
   isotropic diametral cut-planes through the cell centre and divided by the
   cut circumference ``2πR``. This is the channel myosin contraction drives.
2. ``gamma_rigid`` — *rigid / Lagrange-multiplier method-of-planes*. The
   M-SHAKE-constrained actin backbone carries tension as a Lagrange multiplier
   ``λ`` (it is a rigid constraint, not a harmonic bond), recovered as
   ``T_bond = λ·r₀/Δt`` and projected exactly like the soft channel. Without
   this channel the backbone stress (the dominant cortex load path) is invisible
   to a soft-bond-only sum (it under-reports ~200×).
3. ``gamma_passive`` — *passive / Young-Laplace turgor*. The resting osmotic
   turgor pre-tensions the closed cortex shell via Young-Laplace
   ``γ_passive = ΔP·R/2``, with ``ΔP`` the enclosed-volume osmotic pressure
   (``turgor_dP0`` at the construction volume V=V0). This is a passive,
   pressure-balance contribution — NOT actomyosin-generated tension.

GATE-B separation rule (B3 — HARD)
----------------------------------
``gamma_passive`` (turgor) MUST stay a separate channel and is **never** folded
into a single hidden "total". A measurement that happens to land an aggregate
in-band does **NOT** mean the active cortex has closed the band: the passive
turgor term can carry an arbitrary fraction of it. ``gamma_structural`` is the
*mechanical* (actomyosin-generated) tension ``gamma_soft + gamma_rigid`` ONLY;
the turgor channel is reported alongside it, never inside it. Always read the
channels separately so the active-vs-passive split is visible.

Cortical-tension reference band
-------------------------------
The literature cortical-tension band ``[0.35, 0.65] mN/m`` (Salbreux, Charras &
Paluch 2012, *Trends Cell Biol*) is exposed here as a documentation constant
:data:`CORTICAL_TENSION_BAND_N_PER_M`. **Nothing in this module is tuned to it**
— it is provided so callers can overlay the band; the estimator computes γ from
the cell's mechanical state alone.

Provenance
----------
- Soft method-of-planes: factored from
  ``ffn_sim/scripts/h3_ku35_tension.py:_tension_method_of_planes``.
- Rigid Lagrange method-of-planes: factored from
  ``ffn_sim/scripts/h3_ku35_tension.py:_tension_method_of_planes_rigid``
  (``RIGID_LAGRANGE_TENSION_DESIGN.md`` §2, R1).
- Passive Young-Laplace: ``ffn_sim/cortex/enclosed_volume.py`` (turgor_dP0).

The ``|û·n̂|`` absolute value on the cut-plane projection is REQUIRED in both
method-of-planes channels: the raw HOOMD A→B bond direction has an arbitrary
sign relative to a cut plane, so without ``|·|`` the per-plane sum cancels as
√N instead of N, fabricating a ~97× suppressed "isotropy" floor (estimator
audit 2026-06-02; corrected == Irving-Kirkwood to <0.5% on a synthetic shell).
The scalar tension ``T`` keeps its own sign (+ tension / − compression).
"""

from __future__ import annotations

from typing import Any

import numpy as np

# Literature cortical-tension band [N/m] (Salbreux, Charras & Paluch 2012,
# Trends Cell Biol 22(10):536-545). Documentation constant ONLY — the estimator
# is NOT tuned to it; callers may overlay it.
CORTICAL_TENSION_BAND_N_PER_M: tuple[float, float] = (0.35e-3, 0.65e-3)

# Number of isotropic diametral cut-planes for the method-of-planes average.
# Fixed by the estimator design (Fibonacci sphere sampling), not a tuned knob.
_N_PLANES_DEFAULT: int = 12


def _fibonacci_plane_normals(n_planes: int) -> np.ndarray:
    """Return ``n_planes`` near-isotropic unit normals on the sphere.

    Fibonacci-lattice sampling so the diametral cut planes are spread evenly
    over orientation (the same convention used by both method-of-planes
    channels, so soft and rigid contributions are sampled consistently).

    Args:
        n_planes: Number of cut-plane orientations to sample.

    Returns:
        ``(n_planes, 3)`` array of unit normal vectors.
    """
    phi = (1.0 + 5.0**0.5) / 2.0
    i = np.arange(n_planes, dtype=np.float64)
    z = 1.0 - 2.0 * (i + 0.5) / n_planes
    rxy = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    theta = 2.0 * np.pi * i / phi
    return np.stack([rxy * np.cos(theta), rxy * np.sin(theta), z], axis=1)


def _method_of_planes_gamma(
    rA: np.ndarray,
    rB: np.ndarray,
    u: np.ndarray,
    T: np.ndarray,
    R_cell: float,
    n_planes: int,
) -> float:
    """Method-of-planes cortical tension γ [N/m] from per-bond scalar tensions.

    Shared kernel for both the soft and rigid channels. For each of
    ``n_planes`` diametral planes through the cell centre, sum the tensile
    force ``T·|û·n̂|`` of every bond that crosses the plane, divide by the cut
    circumference ``2πR_cell``, and average ``|γ|`` over plane orientations.

    Args:
        rA: ``(M, 3)`` bond endpoint A positions [m] (centroid-relative or
            absolute — only the sign of the dot with each plane normal matters
            for crossing detection, so an isotropic average is centre-agnostic
            for a centred shell).
        rB: ``(M, 3)`` bond endpoint B positions [m].
        u: ``(M, 3)`` unit bond directions (A→B).
        T: ``(M,)`` signed scalar bond tensions [N] (+ tension, − compression).
        R_cell: Cell radius [m] (the great-circle radius for the cut
            circumference).
        n_planes: Number of cut-plane orientations to average over.

    Returns:
        Mean over orientations of ``|γ|`` [N/m]. ``0.0`` if there are no bonds.
    """
    if rA.shape[0] == 0:
        return 0.0
    normals = _fibonacci_plane_normals(n_planes)
    gammas: list[float] = []
    circ = 2.0 * np.pi * R_cell
    for n_hat in normals:
        a_side = rA @ n_hat
        b_side = rB @ n_hat
        crossing = (a_side * b_side) < 0.0
        if not crossing.any():
            gammas.append(0.0)
            continue
        # |û·n̂| (absolute) is REQUIRED — see module docstring (√N vs N floor).
        f_cut = T[crossing] * np.abs(u[crossing] @ n_hat)
        gammas.append(float(np.sum(f_cut)) / circ)
    return float(np.mean(np.abs(np.asarray(gammas))))


def _read_tag_ordered_positions_and_bonds(
    sim: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Read tag-ordered positions + bond topology from a HOOMD simulation.

    Args:
        sim: A built HOOMD ``Simulation`` (its state holds particles + bonds).

    Returns:
        Tuple ``(pos_by_tag, bond_group, bond_typeid, bond_types)``:

        * ``pos_by_tag``: ``(N, 3)`` positions re-indexed into tag order [m].
        * ``bond_group``: ``(B, 2)`` per-bond endpoint tag indices.
        * ``bond_typeid``: ``(B,)`` per-bond type id.
        * ``bond_types``: list of bond-type name strings.
    """
    bond_types = list(sim.state.bond_types)
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        pos_by_tag = pos[inv].copy()
        bond_group = np.asarray(s.bonds.group).copy()
        bond_typeid = np.asarray(s.bonds.typeid).copy()
    return pos_by_tag, bond_group, bond_typeid, bond_types


def _gamma_soft(sim: Any, R_cell: float, n_planes: int) -> dict[str, float]:
    """Soft (harmonic-bond) method-of-planes cortical tension γ [N/m].

    Reuses the formula from
    ``scripts/h3_ku35_tension.py:_tension_method_of_planes``: every harmonic
    bond present in the integrator carries ``T = k·(|r| − r₀)``; these are
    projected onto isotropic cut planes (see :func:`_method_of_planes_gamma`).

    Args:
        sim: Built HOOMD ``Simulation`` with an integrator carrying a harmonic
            bond force.
        R_cell: Cell radius [m].
        n_planes: Number of cut-plane orientations.

    Returns:
        Dict with ``gamma`` [N/m] and ``n_bonds`` (int). ``gamma`` is ``nan``
        if no harmonic bond force is found on the integrator.
    """
    pos_by_tag, bg, bt, bond_types = _read_tag_ordered_positions_and_bonds(sim)

    # Locate the harmonic bond force on the integrator (the force whose params
    # cover the registered bond-type names).
    bond_force = None
    integrator = sim.operations.integrator
    if integrator is not None:
        for f in integrator.forces:
            params = getattr(f, "params", None)
            if params is not None and any(t in bond_types for t in params):
                bond_force = f
                break
    if bond_force is None or bg.shape[0] == 0:
        return {"gamma": float("nan"), "n_bonds": int(bg.shape[0])}

    # Per-bond k, r0 (tag-ordered frame).
    k_arr = np.zeros(bg.shape[0], dtype=np.float64)
    r0_arr = np.zeros(bg.shape[0], dtype=np.float64)
    for i, tname in enumerate(bond_types):
        try:
            kp = bond_force.params[tname]
            mask = bt == i
            k_arr[mask] = float(kp["k"])
            r0_arr[mask] = float(kp["r0"])
        except Exception:  # noqa: BLE001  (a type with no harmonic params → 0)
            pass

    rA = pos_by_tag[bg[:, 0]]
    rB = pos_by_tag[bg[:, 1]]
    d = rB - rA
    L = np.linalg.norm(d, axis=1)
    L_safe = np.where(L > 0.0, L, 1.0)
    u = d / L_safe[:, None]
    T = k_arr * (L - r0_arr)  # signed scalar tension [N]
    gamma = _method_of_planes_gamma(rA, rB, u, T, R_cell, n_planes)
    return {"gamma": float(gamma), "n_bonds": int(bg.shape[0])}


def _gamma_rigid(
    sim: Any,
    lambda_accumulator: Any,
    R_cell: float,
    dt: float | None,
    n_planes: int,
) -> dict[str, float]:
    """Rigid (Lagrange-multiplier) method-of-planes cortical tension γ [N/m].

    Reuses the formula from
    ``scripts/h3_ku35_tension.py:_tension_method_of_planes_rigid``
    (``RIGID_LAGRANGE_TENSION_DESIGN.md`` §2, R1): the M-SHAKE constraint
    carries the rigid backbone's bond tension as a Lagrange multiplier ``λ``,
    converted via ``T_bond = λ·r₀/Δt`` and projected like the soft channel.

    Args:
        sim: Built HOOMD ``Simulation``.
        lambda_accumulator: The constrained-BAOAB Action that captures λ
            (must expose ``.lambda_buf``, ``.chains_tag_stacked`` and
            ``._chain_rest_length``; e.g. the cell's ``baoab_action`` built
            with ``record_lambda=True``). ``None`` for an unconstrained cell.
        R_cell: Cell radius [m].
        dt: Constrained integrator timestep [s] (required to convert λ → T).
        n_planes: Number of cut-plane orientations.

    Returns:
        Dict with ``gamma`` [N/m] and ``available`` (bool). ``gamma`` is
        ``0.0`` and ``available`` False when there is no Lagrange buffer
        (unconstrained run, ``record_lambda=False``, or no step taken yet) —
        the rigid channel is simply absent, not an error.
    """
    if lambda_accumulator is None:
        return {"gamma": 0.0, "available": False}
    lam = getattr(lambda_accumulator, "lambda_buf", None)
    chains = getattr(lambda_accumulator, "chains_tag_stacked", None)
    r0_attr = getattr(lambda_accumulator, "_chain_rest_length", None)
    if (
        lam is None
        or not isinstance(lam, np.ndarray)
        or chains is None
        or getattr(chains, "ndim", 0) != 2
        or r0_attr is None
        or dt is None
        or dt <= 0.0
    ):
        return {"gamma": 0.0, "available": False}
    r0 = float(r0_attr)

    pos_by_tag, _bg, _bt, _bond_types = _read_tag_ordered_positions_and_bonds(sim)
    rA = pos_by_tag[chains[:, :-1]]  # (F, m, 3)
    rB = pos_by_tag[chains[:, 1:]]   # (F, m, 3)
    d = rB - rA
    L = np.linalg.norm(d, axis=-1)
    L_safe = np.where(L > 0.0, L, 1.0)
    u = d / L_safe[..., None]
    T = lam * (r0 / dt)  # signed scalar tension [N] = λ·r₀/Δt

    gamma = _method_of_planes_gamma(
        rA.reshape(-1, 3),
        rB.reshape(-1, 3),
        u.reshape(-1, 3),
        T.reshape(-1),
        R_cell,
        n_planes,
    )
    return {"gamma": float(gamma), "available": True}


def _gamma_passive(p_enclosed_volume: Any, R_cell: float) -> dict[str, float]:
    """Passive Young-Laplace cortical tension from osmotic turgor [N/m].

    ``γ_passive = ΔP·R/2`` (Young-Laplace for a sphere), with ``ΔP`` the
    enclosed-volume osmotic pressure. At the construction volume (V = V0) the
    elastic bulk term vanishes, so ``ΔP = turgor_dP0`` (the resting turgor
    ``Π₀``). This is the pressure that pre-tensions the closed cortex shell
    BEFORE any actomyosin tension — a passive contribution, reported separately.

    Args:
        p_enclosed_volume: A ``ResolvedEnclosedVolume`` (exposes ``turgor_dP0``
            [Pa]; ``dP_ref`` is the fallback documentation anchor). ``None``
            when the cell has no enclosed-volume compartment.
        R_cell: Cell radius [m].

    Returns:
        Dict with ``gamma`` [N/m] and ``dP`` [Pa]. ``gamma`` is ``0.0`` when
        there is no enclosed-volume compartment (no turgor → no passive
        pre-tension).
    """
    if p_enclosed_volume is None:
        return {"gamma": 0.0, "dP": 0.0}
    # At V = V0 the elastic bulk term is zero, so ΔP = turgor_dP0.
    dP = float(getattr(p_enclosed_volume, "turgor_dP0", 0.0))
    gamma = dP * R_cell / 2.0
    return {"gamma": float(gamma), "dP": dP}


def measure_cortical_tension(
    sim: Any,
    *,
    R_cell: float,
    p_enclosed_volume: Any | None = None,
    lambda_accumulator: Any | None = None,
    dt: float | None = None,
    n_planes: int = _N_PLANES_DEFAULT,
) -> dict[str, Any]:
    """Measure cortical tension γ through three separate channels.

    The single GATE-B measurement entry point. Reports the soft (active
    harmonic-bond), rigid (M-SHAKE Lagrange), and passive (Young-Laplace
    turgor) channels SEPARATELY (B3 / GATE-B rule), plus the *mechanical*
    structural total ``gamma_structural = gamma_soft + gamma_rigid``.

    GATE-B separation rule (HARD): ``gamma_passive`` (turgor) is NEVER folded
    into ``gamma_structural`` or any single hidden "total". An in-band aggregate
    must NOT be read as active-cortex closure — the passive turgor term can
    carry an arbitrary share of it. Always read the channels separately so the
    active-vs-passive(turgor) split stays visible.

    Args:
        sim: A built HOOMD ``Simulation`` (e.g. ``cell.simulation``).
        R_cell: Cell radius [m] (e.g. ``cell.p_cortex.R_cell``).
        p_enclosed_volume: ``ResolvedEnclosedVolume`` for the passive turgor
            channel (e.g. ``cell.p_enclosed_volume``). ``None`` → no passive
            contribution.
        lambda_accumulator: Constrained-BAOAB Action that captures the SHAKE
            Lagrange multipliers (e.g.
            ``cell.extras["handles"]["baoab_action"]`` built with
            ``record_lambda=True``). ``None`` → no rigid contribution.
        dt: Constrained integrator timestep [s], needed to convert λ → tension.
            Required for a non-zero rigid channel.
        n_planes: Number of method-of-planes cut-plane orientations.

    Returns:
        A dict with the channel summary and detail::

            {
              "gamma_soft": float,        # active harmonic-bond MOP [N/m]
              "gamma_rigid": float,       # rigid M-SHAKE Lagrange MOP [N/m]
              "gamma_passive": float,     # passive Young-Laplace turgor [N/m]
              "gamma_structural": float,  # = gamma_soft + gamma_rigid [N/m]
                                          #   (mechanical / actomyosin ONLY —
                                          #    turgor is NOT included here)
              "R_cell": float,            # [m]
              "band_N_per_m": (lo, hi),   # literature band overlay (NOT a fit)
              "channels": {
                  "soft":    {"gamma": float, "n_bonds": int},
                  "rigid":   {"gamma": float, "available": bool},
                  "passive": {"gamma": float, "dP": float},  # dP in Pa
              },
            }
    """
    soft = _gamma_soft(sim, R_cell, n_planes)
    rigid = _gamma_rigid(sim, lambda_accumulator, R_cell, dt, n_planes)
    passive = _gamma_passive(p_enclosed_volume, R_cell)

    g_soft = soft["gamma"]
    g_rigid = rigid["gamma"]
    g_passive = passive["gamma"]

    # Structural = mechanical (actomyosin-generated) tension ONLY. The passive
    # turgor channel is reported alongside, NEVER summed in (GATE-B B3 rule).
    # nan-safe: a missing soft channel (no harmonic bond force) must not poison
    # the structural total when the rigid channel is present.
    g_struct = float(np.nansum([g_soft, g_rigid]))

    return {
        "gamma_soft": float(g_soft),
        "gamma_rigid": float(g_rigid),
        "gamma_passive": float(g_passive),
        "gamma_structural": g_struct,
        "R_cell": float(R_cell),
        "band_N_per_m": CORTICAL_TENSION_BAND_N_PER_M,
        "channels": {
            "soft": soft,
            "rigid": rigid,
            "passive": passive,
        },
    }
