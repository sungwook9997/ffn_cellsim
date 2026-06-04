"""Layer-2 §E: coherent collective traction — SPP plithotaxis (the DRIVING/COORDINATION fix).

§B2/§C/§D eliminated scale, measurement-definition, and cohesion as the ~5-9x magnitude-gap
suspect; the `daxis_yield` table is decisive — a radially-symmetric outward crawl, uniform on
every basal cell and recomputed each epoch, does **net ~0** to spread the sheet (A/A0 1.35 at
f=0 vs 1.34 at f=9.4 nN). That force is an *isotropic internal pressure*: cohesion (surface
tension) balances it exactly, so the cluster sits at a static force-balance equilibrium and
raising f only moves the equilibrium, it does not spread. The missing ingredient is
**coordination, not force magnitude** — persistence, spatial correlation, and free-edge polarity.

This module is the literature-faithful active-matter rendering of *plithotaxis* / collective
epithelial migration (Trepat 2009; Tambe 2011) on the center-based CBM: each cell carries an
in-plane **polarization** `p_i = (cos θ_i, sin θ_i)` and is self-propelled along it, the traction
reacted by the **SUBSTRATE** (drag = clutch γ_cell, NOT cohesion-competing — the `substrate_crawl`
lineage). The polarity evolves by the **Smeets-2016 CIL-SPP rule** (`PNAS 113:14621`, the closest
breast-epithelial-anchored model — MCF10A):

    θ̇_i = − f_cil · w_edge,i · (θ_i − θ_i^free)  +  √(2 D_r) · η_i(t)        [ + J·sin(θ̄_i − θ_i) ]

* **Persistence** (rotational diffusion `D_r`, persistence time τ = 1/D_r ≈ 20 min, MCF10A
  Smeets 2016): a cell keeps crawling in one direction over τ instead of being re-randomised —
  the ingredient the §D radial field wholly lacked.
* **CIL free-edge polarity** (`f_cil` ≈ 0.1 min⁻¹, Smeets 2016): a cell at the free boundary
  repolarises toward open space (away from the weighted mean position of its contacting
  neighbours), `θ_i^free`; the rate is gated by an edge-ness weight `w_edge,i` (geometry) so bulk
  cells with no free space are not spuriously driven. This is the directed, continuously-advancing
  front — not an isotropic pressure.
* **Emergent correlation, NOT imposed alignment.** Garcia 2015 (PNAS 112:15314) + Henkes/
  Marchetti 2020 (Nat Commun) show the ~200 µm velocity-correlation length of collective
  epithelia **emerges from per-cell persistence + the cohesive elastic coupling, with NO explicit
  Vicsek alignment**. So the optional Vicsek term `J` defaults **OFF** (imposing it would add an
  unanchored tunable — ξ emerges without it — and double-count the cohesion's coordination). The
  catch+turnover cohesion (`cadherin_bonds`, `yield_remodel=True`) supplies the elastic coupling;
  the emergent ξ is *measured* and overlay-validated against ~200 µm, never fitted.

The result is a continuously migrating, ductile-cohesion sheet (not a static equilibrium) whose
footprint grows over the 60 h window → A/A0 can rise toward the PI magnitude. Every INPUT
(`D_r`, `f_cil`, `f_active`) is a measured value; ξ and intercellular stress are emergent overlay
targets; nothing is fitted to the PI A/A0 (overlay-only HARD rule). The per-step force is a fixed
`SettableForce` vector (GPU-resident); the polarity update is host-side numpy + one cKDTree per
epoch (the same cadence/pattern as the existing crawl recompute), so the BAOAB integrator is
untouched. See `outputs/layer2/DESIGN_plithotaxis.md` and
`outputs/tag_kb/SE_REGISTRATION_CANDIDATES_2026-06-04_collective-migration.md`.

Sanity Gate
-----------
- Dimensional: θ [rad]; D_r, f_cil, J [s⁻¹]; Δt [s]; f_active [N]; force [N]. `√(2 D_r Δt)` [rad];
  `f_cil·Δt`, `J·Δt` [dimensionless·rad]. ✓
- Boundary: f_active = 0 ⇒ zero force (reduces to the D-axis substrate baseline). f_cil=J=0 ⇒
  pure rotational diffusion ⇒ isotropic, COM diffusive not ballistic, net spread → 0 (recovers
  the §D null as a *limit*, confirming the new spread comes from coordination). A cell far above
  the substrate (z ≫ z_substrate + basal_band) feels ~0 force (only basal cells crawl).
- Sign/sense: larger f_cil ⇒ stronger outward front ⇒ more spread; larger D_r (smaller τ) ⇒
  shorter persistence ⇒ less spread; larger J ⇒ longer ξ. CIL turns θ toward θ^free (outward).
- Conservation: z-force ≡ 0 (the substrate reservoir reacts z); count conserved. In-plane net
  force is NOT required zero — the front is genuinely propelled against substrate drag (unlike a
  zero-mean random self-propulsion), which is the point.
- Numerical: the discrete angle update sub-steps adaptively so the per-sub-step deterministic
  increment stays < 0.2 rad and the noise std < 0.25 rad — the discrete map tracks the SDE; this
  is host-side and does NOT touch the MD integrator/timestep.
- Measurement consistency: the emergent velocity-correlation length and front velocity are
  reported against the measured ξ≈200 µm (Petitjean 2010) and v_front (Poujade 2007) as an
  overlay self-consistency check — never fitted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy.spatial import cKDTree

from ffn_sim.spheroid.params import ResolvedL2

__all__ = ["ResolvedPlithotaxis", "resolve_plithotaxis", "PolarizationField"]

# --- Anchored MCF10A breast-epithelial CIL-SPP constants (Smeets et al. 2016, PNAS 113:14621;
# arXiv:1612.06901). MSD-fit direct values; the closest available breast-epithelial proxy for
# MCF7 (no MCF7-specific data exists — deep-research wf_e0df20b3-470). NOT fitted to PI A/A0. ---
_D_R_PER_MIN: float = 0.05      # min⁻¹  rotational diffusion (⇒ persistence τ = 1/D_r = 20 min)
_F_CIL_PER_MIN: float = 0.1     # min⁻¹  CIL repolarisation rate toward free space
# n_bulk = first-shell coordination of a cell fully embedded in a 2D sheet (hex packing) — the
# geometric reference for edge-ness; grid-invariant (a packing number, not a tuned constant).
_N_BULK_2D: int = 6


@dataclass(frozen=True)
class ResolvedPlithotaxis:
    """Anchored SPP-plithotaxis parameters (SI). See module docstring + SE candidates file."""

    D_r: float                   # s⁻¹   rotational diffusion (persistence τ = 1/D_r)
    f_cil: float                 # s⁻¹   CIL repolarisation rate toward free space
    J: float                     # s⁻¹   OPTIONAL Vicsek neighbour-alignment rate (default 0)
    interaction_radius: float    # m     in-plane neighbour/contact radius (free-edge detection)
    n_bulk: int                  # –     bulk first-shell coordination (edge-ness reference)
    tau: float                   # s     persistence time = 1/D_r (provenance/readout)
    psi: float                   # –     f_cil/(2 D_r) — Smeets MCF10A consistency check (≈1)


def resolve_plithotaxis(resolved: ResolvedL2, *, J: float = 0.0) -> ResolvedPlithotaxis:
    """Resolve the anchored CIL-SPP parameters (no sim; literature constants → SI).

    Args:
        resolved: resolved Layer-2 params (provides ``morse_r0`` for the geometric neighbour radius).
        J: optional Vicsek alignment rate (s⁻¹). Default 0 (OFF) — the ~200 µm correlation emerges
            from persistence + cohesion (Garcia 2015 / Henkes 2020); imposing J adds an unanchored
            tunable. Set > 0 only to test the alignment-augmented variant.

    Returns:
        ``ResolvedPlithotaxis`` with D_r, f_cil in SI (s⁻¹), the geometric neighbour radius, and the
        ψ = f_cil/(2 D_r) consistency check (≈ 1 reproduces Smeets' MCF10A estimate).
    """
    if J < 0.0:
        raise ValueError("J (Vicsek alignment rate) must be >= 0.")
    d_r = _D_R_PER_MIN / 60.0    # min⁻¹ → s⁻¹
    f_cil = _F_CIL_PER_MIN / 60.0
    r0 = float(resolved.morse_r0)
    # In-plane contact radius: two cohesive cells sit at ~r0; 1.3·r0 captures the first contact
    # shell (geometry, scales with r0 — grid-invariant, not a tuned constant).
    interaction_radius = 1.3 * r0
    return ResolvedPlithotaxis(
        D_r=d_r,
        f_cil=f_cil,
        J=float(J),
        interaction_radius=interaction_radius,
        n_bulk=_N_BULK_2D,
        tau=1.0 / d_r,
        psi=f_cil / (2.0 * d_r),
    )


class PolarizationField:
    """Stateful per-cell in-plane polarity θ evolving by the Smeets-2016 CIL-SPP rule.

    Holds one angle per pool slot (``n_max``) so a void→cell division daughter already carries a
    random initial polarity. ``update`` evolves θ over an epoch (host-side; one cKDTree); ``forces``
    maps θ → an in-plane, substrate-reacted, basal-gated force in the order of ``positions[active]``
    so it drops into the existing proliferation epoch loop unchanged.
    """

    def __init__(self, n_max: int, plith: ResolvedPlithotaxis, *, seed: int) -> None:
        self._p = plith
        self._rng = np.random.default_rng(seed)
        # random initial polarity for every slot (new daughters inherit a defined random angle)
        self._theta = self._rng.uniform(-np.pi, np.pi, int(n_max))

    @property
    def theta(self) -> npt.NDArray[np.float64]:
        """Current polarity angles (n_max,) [rad] — read-only view for diagnostics/tests."""
        return self._theta

    def polarization(self, active: npt.NDArray[np.bool_]) -> npt.NDArray[np.float64]:
        """Unit polarity vectors p̂ = (cosθ, sinθ) for the active cells, in active order (N,2)."""
        th = self._theta[active]
        return np.column_stack([np.cos(th), np.sin(th)])

    def update(
        self,
        positions: npt.NDArray[np.float64],
        active: npt.NDArray[np.bool_],
        dt: float,
        *,
        z_substrate: float,
        basal_band: float,
    ) -> None:
        """Evolve θ over one epoch of duration ``dt`` (s) for the basal active cells.

        CIL + rotational diffusion (+ optional Vicsek) act on cells touching the substrate
        (``z ≤ z_substrate + basal_band``); the in-plane neighbour graph and free-edge direction
        are computed among those basal cells. Non-basal active cells get pure rotational diffusion
        (their force is 0 anyway). Adaptive host-side sub-stepping keeps the discrete map faithful
        to the SDE without touching the MD integrator.
        """
        if dt <= 0.0:
            raise ValueError("dt must be > 0.")
        if basal_band <= 0.0:
            raise ValueError("basal_band must be > 0.")
        p = self._p
        idx = np.where(active)[0]
        if idx.size == 0:
            return
        pos = np.asarray(positions, dtype=np.float64)[idx]
        dz = pos[:, 2] - z_substrate
        basal = (dz <= basal_band) & (dz >= -basal_band)

        # --- adaptive sub-stepping: keep per-sub-step deterministic increment < 0.2 rad and the
        # noise std < 0.25 rad so the discrete update tracks θ̇ = −f_cil(θ−θ^free)+√(2D_r)η. ---
        det_rate = max(p.f_cil, p.J)
        n_sub = 1
        if det_rate * dt > 0.2:
            n_sub = max(n_sub, int(np.ceil(det_rate * dt / 0.2)))
        if 2.0 * p.D_r * dt > 0.25**2:
            n_sub = max(n_sub, int(np.ceil(2.0 * p.D_r * dt / 0.25**2)))
        dt_sub = dt / n_sub
        noise_std = np.sqrt(2.0 * p.D_r * dt_sub)

        # free-edge geometry (computed once per epoch from the start-of-epoch basal positions —
        # positions are fixed during the host-side angle relaxation, so θ^free / w_edge are constant)
        theta_free = np.full(idx.size, np.nan)        # NaN where undefined (no neighbours)
        w_edge = np.ones(idx.size)                    # isolated cell ⇒ fully "edge"
        theta_bar = np.full(idx.size, np.nan)         # mean-neighbour polarity (Vicsek, if used)
        b_idx = np.where(basal)[0]
        if b_idx.size >= 2:
            bxy = pos[b_idx, :2]
            tree = cKDTree(bxy)
            neigh = tree.query_ball_point(bxy, r=p.interaction_radius)
            th_b = self._theta[idx[b_idx]]
            for li, nb in enumerate(neigh):
                others = [k for k in nb if k != li]
                ncount = len(others)
                w_edge[b_idx[li]] = np.clip(1.0 - ncount / p.n_bulk, 0.0, 1.0)
                if ncount > 0:
                    cen = bxy[others].mean(axis=0)          # neighbour centroid
                    d = bxy[li] - cen                        # away-from-neighbours = toward free space
                    if d[0] != 0.0 or d[1] != 0.0:
                        theta_free[b_idx[li]] = np.arctan2(d[1], d[0])
                    if p.J > 0.0:
                        cx = np.cos(th_b[others]).sum(); cy = np.sin(th_b[others]).sum()
                        if cx != 0.0 or cy != 0.0:
                            theta_bar[b_idx[li]] = np.arctan2(cy, cx)

        th = self._theta[idx].copy()
        has_free = ~np.isnan(theta_free)
        has_bar = ~np.isnan(theta_bar)
        for _ in range(n_sub):
            dth = np.zeros(idx.size)
            # CIL: turn toward free space, gated by edge-ness (basal cells with a defined θ^free)
            cil = basal & has_free
            dth[cil] += p.f_cil * w_edge[cil] * dt_sub * np.sin(theta_free[cil] - th[cil])
            # optional Vicsek alignment (default OFF, J=0 ⇒ skipped)
            if p.J > 0.0:
                al = basal & has_bar
                dth[al] += p.J * dt_sub * np.sin(theta_bar[al] - th[al])
            # rotational diffusion (all active cells — persistence)
            dth += noise_std * self._rng.standard_normal(idx.size)
            th = th + dth
        # wrap to (-pi, pi]
        self._theta[idx] = (th + np.pi) % (2.0 * np.pi) - np.pi

    def forces(
        self,
        positions: npt.NDArray[np.float64],
        active: npt.NDArray[np.bool_],
        *,
        f_active: float,
        z_substrate: float,
        basal_band: float,
    ) -> npt.NDArray[np.float64]:
        """In-plane, substrate-reacted self-propulsion force per ACTIVE cell (N_active, 3).

        Force = ``f_active · w_basal · p̂`` with z-component ≡ 0 (the substrate reacts z). The basal
        membership weight ``w_basal`` is the same half-cosine taper as ``substrate_crawl`` (geometry).
        Returned in the order of ``positions[active]`` so it scatters straight into the epoch loop.
        """
        idx = np.where(active)[0]
        out = np.zeros((idx.size, 3), dtype=np.float64)
        if f_active == 0.0 or idx.size == 0:
            return out
        if basal_band <= 0.0:
            raise ValueError("basal_band must be > 0.")
        pos = np.asarray(positions, dtype=np.float64)[idx]
        dz = pos[:, 2] - z_substrate
        w = np.zeros(idx.size)
        in_band = (dz <= basal_band) & (dz >= -basal_band)
        w[in_band] = 0.5 * (1.0 + np.cos(np.pi * np.clip(np.abs(dz[in_band]) / basal_band, 0.0, 1.0)))
        th = self._theta[idx]
        out[:, 0] = f_active * w * np.cos(th)
        out[:, 1] = f_active * w * np.sin(th)
        # out[:, 2] stays 0 — crawl is tangent to the dish (substrate reacts z).
        return out
