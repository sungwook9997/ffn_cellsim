"""Slater 2021 radial stress(r) estimator — method-of-planes on radial shells.

A SECOND, independent stress estimator for the contractile cortex/cell, distinct
from both the scalar cortical-tension γ (:mod:`ffn_sim.cortex.cortical_tension`)
and the coherent focal-adhesion traction differential. Where the γ soft channel
projects soft-bond tensions onto *diametral cut-planes through the centre* and
divides by the great-circle circumference (a surface tension γ [N/m]), this
estimator projects the SAME bond tensions onto *concentric shells of radius r*
and divides by the shell area, giving a radial *stress* profile σ(r) [Pa] — the
quantity Slater et al. 2021 (*Soft Matter* 17:10274, DOI 10.1039/d0sm01911a;
same Taeyoon-Kim discrete actomyosin+ECM family) use to characterise how a
contracting network transmits stress outward.

Why a second estimator
----------------------
The magnitude line's active-γ readout is the per-SF coherent FA-traction
differential, which loop24b found seed-unstable. σ(r) reads the SAME realised
bond-tension state through an INDEPENDENT geometric projection (radial shells,
not the FA differential), so disagreement between the two flags an estimator
artefact rather than physics, and agreement raises confidence. It is also the
readout the Slater cross-model transmission oracle (O1: σ(r) ∝ r⁻² in 3D) needs
— implemented here as a MEASUREMENT only; wiring it as an acceptance GATE is a
gate-contract change reserved for the PI (see
``docs/SLATER2021_TRANSMISSION_ORACLE_DRAFT_2026-06-11.md``).

Method (radial method-of-planes)
--------------------------------
For a shell of radius ``r`` centred on the cell centre, a soft bond (endpoints
``rA``, ``rB`` relative to the centre) CROSSES the shell iff
``(|rA| − r)·(|rB| − r) < 0``. Each crossing bond transmits its scalar tension
``T = k·(L − r₀)`` across the shell; the radial component is ``T·|û·r̂|`` with
``r̂`` the outward radial unit vector at the bond midpoint. The stress is the
crossing sum over the shell area::

    σ(r) = Σ_crossing  T·|û·r̂|  /  A(r)
    A(r) = 4πr²        (sphere, 3D single cell  → O1 expects σ ∝ r⁻²)
         = 2πr·height  (cylinder, 2D-like slab  → O1 expects σ ∝ r⁻¹)

The ``|û·r̂|`` absolute value is REQUIRED for the same reason the cut-plane
channel needs ``|û·n̂|`` (:mod:`ffn_sim.cortex.cortical_tension` module
docstring): the raw HOOMD A→B bond direction has an arbitrary sign relative to
the shell normal, so an unsigned per-bond projection would cancel as √N instead
of N. The scalar tension ``T`` keeps its own sign (+ tension / − compression),
so σ(r) is signed.

Sanity Gate
-----------
- *Dimensional*: ``T`` [N] · dimensionless / ``A`` [m²] = [Pa]. ✔
- *Boundary*: no bonds, or no bonds crossing a shell → ``σ(r) = 0`` (not nan). ✔
- *Conservation / scaling*: on a force-balanced isotropic radial field (every
  shell crossed by the same total radial force ``F``), ``σ(r)·A(r) = F`` is
  constant → ``σ ∝ r⁻²`` (sphere) / ``σ ∝ r⁻¹`` (cylinder). Verified in
  ``scripts/h7_radial_stress_validate.py`` and ``tests/test_radial_stress.py``.
- *Sign-sense*: tensile field (T>0) → σ>0; compressive (T<0) → σ<0. ✔
- *Labeling invariance*: swapping a bond's A/B endpoints leaves σ unchanged (the
  ``|û·r̂|`` removes the HOOMD storage-order sign). ✔
- *Protocol consistency*: reuses the SAME soft-bond tension state and cortical
  bond-type filter as the γ soft channel (no second hard-coded param map). ✔

Provenance: the soft-bond tension state + cortical filter match
``cortical_tension.py:_gamma_soft``; the crossing+projection pattern generalises
``cortical_tension.py:_method_of_planes_gamma`` from a flat plane to a shell.
This is a diagnostic/measurement tool, NOT a runtime hot-path module.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from ffn_sim.cortex.cortical_tension import (
    _read_tag_ordered_positions_and_bonds,
    cortical_bond_typeid_mask,
)


def radial_stress_profile(
    rA: np.ndarray,
    rB: np.ndarray,
    u: np.ndarray,
    T: np.ndarray,
    radii: np.ndarray,
    *,
    geometry: str = "sphere",
    height: float | None = None,
) -> np.ndarray:
    """Radial stress σ(r) [Pa] by the radial method-of-planes.

    Args:
        rA: ``(M, 3)`` bond endpoint A positions RELATIVE TO THE CELL CENTRE [m].
        rB: ``(M, 3)`` bond endpoint B positions relative to the cell centre [m].
        u: ``(M, 3)`` unit bond directions (A→B); the sign is arbitrary and is
            removed by the ``|û·r̂|`` projection.
        T: ``(M,)`` signed scalar bond tensions [N] (+ tension, − compression).
        radii: ``(K,)`` shell radii to evaluate [m], strictly > 0.
        geometry: ``"sphere"`` (area ``4πr²``) or ``"cylinder"`` (area
            ``2πr·height``).
        height: Cylinder height [m]; required iff ``geometry == "cylinder"``.

    Returns:
        ``(K,)`` signed radial stress ``σ(r)`` [Pa]; ``0.0`` at radii with no
        crossing bond.
    """
    radii = np.asarray(radii, dtype=np.float64)
    out = np.zeros(radii.shape[0], dtype=np.float64)
    if rA.shape[0] == 0:
        return out

    aR = np.linalg.norm(rA, axis=1)
    bR = np.linalg.norm(rB, axis=1)
    rmid = 0.5 * (rA + rB)
    rmid_n = np.linalg.norm(rmid, axis=1)
    rhat = rmid / np.where(rmid_n > 0.0, rmid_n, 1.0)[:, None]
    # |û·r̂| (absolute) — labeling-invariant; T carries the physical sign.
    Tproj = T * np.abs(np.sum(u * rhat, axis=1))

    if geometry == "sphere":
        area = 4.0 * np.pi * radii ** 2
    elif geometry == "cylinder":
        if height is None or height <= 0.0:
            raise ValueError("geometry='cylinder' requires height > 0 [m]")
        area = 2.0 * np.pi * radii * float(height)
    else:
        raise ValueError(
            f"geometry must be 'sphere' or 'cylinder', got {geometry!r}"
        )

    for i, r in enumerate(radii):
        crossing = (aR - r) * (bR - r) < 0.0
        if not crossing.any():
            continue
        out[i] = float(np.sum(Tproj[crossing])) / float(area[i])
    return out


def radial_crossing_counts(
    rA: np.ndarray, rB: np.ndarray, radii: np.ndarray
) -> np.ndarray:
    """Number of bonds crossing each shell radius (sampling-support diagnostic).

    Args:
        rA: ``(M, 3)`` endpoint A positions relative to the cell centre [m].
        rB: ``(M, 3)`` endpoint B positions relative to the cell centre [m].
        radii: ``(K,)`` shell radii [m].

    Returns:
        ``(K,)`` int array of crossing-bond counts; read it alongside
        :func:`radial_stress_profile` so a low-support shell is not mistaken for
        a physical σ(r) feature.
    """
    radii = np.asarray(radii, dtype=np.float64)
    if rA.shape[0] == 0:
        return np.zeros(radii.shape[0], dtype=np.int64)
    aR = np.linalg.norm(rA, axis=1)
    bR = np.linalg.norm(rB, axis=1)
    return np.array(
        [int(np.count_nonzero((aR - r) * (bR - r) < 0.0)) for r in radii],
        dtype=np.int64,
    )


def fit_power_law(
    radii: np.ndarray,
    stress: np.ndarray,
    *,
    r_min: float | None = None,
    r_max: float | None = None,
) -> dict[str, float]:
    """Least-squares slope of ``log|σ|`` vs ``log r`` over ``[r_min, r_max]``.

    Fits ``σ ≈ a·rᵖ``. The Slater O1 law expects ``p ≈ −2`` (sphere / 3D) or
    ``p ≈ −1`` (cylinder / 2D) in the strong-transmission regime. This is a
    DIAGNOSTIC scaling readout, NOT an acceptance gate — a gate is PI-authored
    (``docs/SLATER2021_TRANSMISSION_ORACLE_DRAFT_2026-06-11.md``).

    Args:
        radii: ``(K,)`` shell radii [m].
        stress: ``(K,)`` σ(r) [Pa] from :func:`radial_stress_profile`.
        r_min: Lower radius bound for the fit window [m] (``None`` → min radius
            with non-zero σ).
        r_max: Upper radius bound for the fit window [m] (``None`` → max radius
            with non-zero σ).

    Returns:
        Dict ``{"exponent", "prefactor", "r2", "n_points"}``. ``exponent`` is
        ``nan`` if fewer than two usable (positive-radius, non-zero-σ) points lie
        in the window.
    """
    radii = np.asarray(radii, dtype=np.float64)
    stress = np.asarray(stress, dtype=np.float64)
    good = (radii > 0.0) & (stress != 0.0) & np.isfinite(stress)
    if r_min is not None:
        good &= radii >= r_min
    if r_max is not None:
        good &= radii <= r_max
    nan = {"exponent": float("nan"), "prefactor": float("nan"),
           "r2": float("nan"), "n_points": int(np.count_nonzero(good))}
    if np.count_nonzero(good) < 2:
        return nan
    x = np.log(radii[good])
    y = np.log(np.abs(stress[good]))
    p, c = np.polyfit(x, y, 1)
    resid = y - (p * x + c)
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
    return {
        "exponent": float(p),
        "prefactor": float(np.exp(c)),
        "r2": float(r2),
        "n_points": int(np.count_nonzero(good)),
    }


def radial_stress_from_sim(
    sim: Any,
    *,
    radii: np.ndarray | None = None,
    n_radii: int = 24,
    r_min: float | None = None,
    r_max: float | None = None,
    cortical_bond_types: set[str] | frozenset[str] | None = None,
    geometry: str = "sphere",
    height: float | None = None,
    center: np.ndarray | None = None,
) -> dict[str, Any]:
    """Radial stress profile σ(r) of a built cell's soft cortical bonds.

    Reads the live harmonic-bond state from a HOOMD simulation, computes each
    cortical bond's scalar tension ``T = k·(L − r₀)`` (the SAME state and
    cortical/adhesion filter as the γ soft channel), centres positions on the
    cortical-bond COM (or a supplied centre), and evaluates
    :func:`radial_stress_profile` on a set of shells.

    Args:
        sim: A built HOOMD ``Simulation`` with an integrator carrying a harmonic
            bond force (e.g. ``cell.simulation``).
        radii: Explicit ``(K,)`` shell radii [m]. ``None`` → log-spaced between
            ``r_min`` and ``r_max``.
        n_radii: Number of log-spaced shells when ``radii`` is ``None``.
        r_min: Inner shell radius [m] when ``radii`` is ``None`` (default
            ``0.05·r_max``).
        r_max: Outer shell radius [m] when ``radii`` is ``None`` (default = the
            largest cortical-bond endpoint radius from the centre).
        cortical_bond_types: Cortical bond-type allowlist override, or ``None``
            for the robust adhesion-denylist default
            (:func:`ffn_sim.cortex.cortical_tension.cortical_bond_typeid_mask`).
        geometry: ``"sphere"`` (3D single cell) or ``"cylinder"``.
        height: Cylinder height [m]; required iff ``geometry == "cylinder"``.
        center: ``(3,)`` cell centre [m]; ``None`` → mean of cortical-bond
            endpoint positions.

    Returns:
        Dict with::

            {
              "radii":      (K,) shell radii [m],
              "stress":     (K,) signed σ(r) [Pa],
              "n_crossing": (K,) crossing-bond count per shell,
              "geometry":   "sphere" | "cylinder",
              "center":     (3,) cell centre used [m],
              "n_bonds":    cortical bonds counted,
              "n_bonds_total": all bonds present,
              "n_bonds_excluded": adhesion/substrate bonds filtered out,
              "fit":        fit_power_law(...) over the non-zero shells,
            }

        ``stress`` is all-nan if no harmonic bond force is found on the
        integrator.
    """
    pos_by_tag, bg, bt, bond_types = _read_tag_ordered_positions_and_bonds(sim)

    # Locate the harmonic bond force on the integrator — same selection rule as
    # cortical_tension.py:_gamma_soft (the force whose params cover the bond
    # types). Kept local rather than refactored into the γ path to leave the
    # GATE-B measurement byte-identical.
    bond_force = None
    integrator = sim.operations.integrator
    if integrator is not None:
        for f in integrator.forces:
            params = getattr(f, "params", None)
            if params is not None and any(t in bond_types for t in params):
                bond_force = f
                break

    n_total = int(bg.shape[0])
    empty = {
        "radii": np.zeros(0), "stress": np.full(0, np.nan),
        "n_crossing": np.zeros(0, dtype=np.int64), "geometry": geometry,
        "center": np.zeros(3), "n_bonds": 0, "n_bonds_total": n_total,
        "n_bonds_excluded": 0,
        "fit": {"exponent": float("nan"), "prefactor": float("nan"),
                "r2": float("nan"), "n_points": 0},
    }
    if bond_force is None or n_total == 0:
        return empty

    # Cortical bond selection (denylist default / allowlist override), then the
    # per-bond (k, r₀) → T = k·(L − r₀) state — provenance: cortical_tension.py
    # :_gamma_soft (lines 354-380).
    type_is_cortical = cortical_bond_typeid_mask(bond_types, cortical_bond_types)
    cortical_mask = type_is_cortical[bt]
    n_cortical = int(np.count_nonzero(cortical_mask))
    bg_c = bg[cortical_mask]
    bt_c = bt[cortical_mask]
    if bg_c.shape[0] == 0:
        return {**empty, "n_bonds": 0, "n_bonds_excluded": n_total}

    k_arr = np.zeros(bg_c.shape[0], dtype=np.float64)
    r0_arr = np.zeros(bg_c.shape[0], dtype=np.float64)
    for i, tname in enumerate(bond_types):
        try:
            kp = bond_force.params[tname]
            mask = bt_c == i
            k_arr[mask] = float(kp["k"])
            r0_arr[mask] = float(kp["r0"])
        except Exception:  # noqa: BLE001  (a type with no harmonic params → 0)
            pass

    rA_abs = pos_by_tag[bg_c[:, 0]]
    rB_abs = pos_by_tag[bg_c[:, 1]]
    if center is None:
        center = 0.5 * (rA_abs + rB_abs).mean(axis=0)
    center = np.asarray(center, dtype=np.float64).reshape(3)
    rA = rA_abs - center
    rB = rB_abs - center

    d = rB - rA
    L = np.linalg.norm(d, axis=1)
    L_safe = np.where(L > 0.0, L, 1.0)
    u = d / L_safe[:, None]
    T = k_arr * (L - r0_arr)  # signed scalar tension [N]

    if radii is None:
        endpoint_r = np.concatenate(
            [np.linalg.norm(rA, axis=1), np.linalg.norm(rB, axis=1)]
        )
        rmax = float(r_max) if r_max is not None else float(endpoint_r.max())
        rmin = float(r_min) if r_min is not None else 0.05 * rmax
        rmin = max(rmin, 1e-12)
        radii = np.logspace(np.log10(rmin), np.log10(rmax), int(n_radii))
    radii = np.asarray(radii, dtype=np.float64)

    stress = radial_stress_profile(
        rA, rB, u, T, radii, geometry=geometry, height=height
    )
    n_crossing = radial_crossing_counts(rA, rB, radii)
    return {
        "radii": radii,
        "stress": stress,
        "n_crossing": n_crossing,
        "geometry": geometry,
        "center": center,
        "n_bonds": n_cortical,
        "n_bonds_total": n_total,
        "n_bonds_excluded": n_total - n_cortical,
        "fit": fit_power_law(radii, stress),
    }
