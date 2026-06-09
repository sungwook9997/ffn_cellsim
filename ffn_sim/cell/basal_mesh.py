"""Connected basal contractile actin mesh — geometry layer (B1; KU-3.5).

PI directive (a), 2026-06-09 (``docs/v2_audit/BASAL_MESH_DESIGN_2026-06-09.md``):
the ventral stress-fiber apparatus is a CONNECTED actin mesh laid on the basal
plane, with the long-axis SF cables woven into it, FA-anchored to the substrate
and contracted by NMII. A connected mesh has MANY filaments, so the Stam-Hocky
bipolar NMII gate (which needs two antiparallel filaments) works — unlike a single
FA→FA chain.

This module is the **geometry layer (increment B1)**: it lays out the basal mesh
as a BIMODAL filament field on a FLAT basal DISK — the exact analogue of the cortex
connected-mesh construction (:func:`ffn_sim.cortex.cortex.generate_bimodal_cortex_layout`)
but projected onto a plane instead of the sphere shell:

* **short Arp2/3-infill filaments** (≈ one segment) → the connecting basal meshwork;
* **long formin filaments** (≈ 10× longer) → the ventral stress-fiber CABLES,
  oriented along the cell long axis (the basal footprint's principal axis);
* DISORDERED-ISOTROPIC short filaments (uniform in-plane azimuth) so the infill
  percolates in every direction (Li-Gao-Xu 2022 cortex-rheology rationale, applied
  in-plane).

The connectivity (bridge-different-filament crosslinkers → giant ≥ 0.9, z ∈ [3,3.5]),
the FA anchoring, the ``sf_myosin_`` NMII placement and the equilibrated Kumar /
Balaban active gate are the FOLLOWING increments (B2–B5). This file adds NO force and
NO HOOMD state — it returns a pure geometry layout (reusing the cortex
:class:`VariableLengthCortexLayout`, which is prefix/force-agnostic).

Why reuse ``VariableLengthCortexLayout``
----------------------------------------
That dataclass is a pure flat-indexed filament topology (positions, per-filament
bead counts/starts, tangents, chain bond/angle groups, ``is_formin`` mask). It
carries no shell-specific assumption and no particle/bond TYPE name, so it equally
describes a planar basal mesh. The downstream snapshot-extender (B3/B4) mints the
``sf_``-prefixed types from this geometry (so every basal bond is γ-denylisted).

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. Checks in ``ffn_sim/tests/test_basal_mesh.py``.*

1. **Dimensional analysis** — ``ell0`` [m] bead spacing; ``footprint_radius`` [m];
   ``z_basal`` [m]; ``band_thickness`` [m]. Lengths only; no force here.
2. **Boundary cases** — ``n_filaments ≤ 0`` → ValueError; ``footprint_radius ≤ 0``
   → ValueError; ``ell0 ≤ 0`` → ValueError; ``formin_fraction`` ∈ [0,1].
3. **Conservation / count** — particle count = Σ N_beads_i; bonds = Σ (N_i − 1);
   angles = Σ (N_i − 2); ``is_formin`` sums to ≈ ``formin_fraction · F``.
4. **Numerical** — all positions finite float64; no NaN.
5. **Sign / sense (PLANAR + force-free)** — every bead's z lies in the basal band
   ``[z_basal − band_thickness, z_basal]`` (the mesh is FLAT, not bulging); every
   CoM lies inside the footprint disk (radial ≤ footprint_radius); every chain bond
   length equals ``ell0`` exactly (born force-free under a Harmonic r0 = ell0).
6. **Measurement-protocol** — long (formin = cable) filaments align with the long
   axis (mean |cos| ≥ ``align_floor``, default 0.8 — the SF-gate convention); short
   (infill) filaments are isotropic in-plane (mean |cos| ≈ 2/π ≈ 0.637).

References
----------
- Bimodal cortex construction (the shell analogue this mirrors):
  :func:`ffn_sim.cortex.cortex.generate_bimodal_cortex_layout`; Fritzsche 2016/2017
  (bimodal short Arp2/3 + long formin); Li-Gao-Xu 2022 (disordered isotropic).
- vSF long-axis alignment (cables): Tojkander 2012; Hotulainen-Lappalainen 2006;
  PI 2026-06-09 long-axis-aligned pairing rule (``select_aligned_fa_pairs``).
- Design: ``docs/v2_audit/BASAL_MESH_DESIGN_2026-06-09.md``.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from ffn_sim.cortex.cortex import VariableLengthCortexLayout


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def _inplane_long_axis(long_axis: np.ndarray | None) -> np.ndarray:
    """Return a unit 2-vector (x, y) for the cable orientation (default +x)."""
    if long_axis is None:
        return np.array([1.0, 0.0], dtype=np.float64)
    v = np.asarray(long_axis, dtype=np.float64).reshape(-1)[:2]
    n = float(np.linalg.norm(v))
    if not (math.isfinite(n) and n > 0.0):
        return np.array([1.0, 0.0], dtype=np.float64)
    return v / n


def generate_basal_mesh_layout(
    *,
    n_filaments: int,
    ell0: float,
    footprint_radius: float,
    z_basal: float,
    formin_fraction: float = 0.12,
    L_short_mean: float | None = None,
    L_long_mean: float | None = None,
    long_axis: np.ndarray | None = None,
    long_axis_jitter_deg: float = 0.0,
    band_thickness: float = 200.0e-9,
    seed: int = 73,
    rng: np.random.Generator | None = None,
    n_beads_min: int = 2,
    n_beads_max: int | None = None,
) -> VariableLengthCortexLayout:
    """Lay out a bimodal basal actin mesh on a FLAT disk (geometry only).

    Mirrors :func:`generate_bimodal_cortex_layout` but the centres of mass are
    sampled in an in-plane DISK (radius ``footprint_radius``, at height
    ``z_basal``) and the tangents lie IN the basal plane: long (formin = cable)
    filaments oriented along ``long_axis`` (± with optional jitter), short (Arp2/3
    = infill) filaments at a uniform in-plane azimuth. No shell projection; each
    filament sits at one depth within the basal band so crossing filaments are
    z-separated (no construction overlap), mirroring the shell-band trick.

    Args:
        n_filaments: number of basal filaments (cables + infill).
        ell0: bead rest length [m] (the cortex ℓ0; the Harmonic r0 downstream).
        footprint_radius: basal contact-disk radius [m] (from the FA footprint).
        z_basal: basal plane height [m] (the substrate-contact z of the cell).
        formin_fraction: fraction of filaments that are LONG cables (Fritzsche
            ~0.10-0.20; default 0.12 = the cortex value).
        L_short_mean: mean infill length [m] (default ``ell0``).
        L_long_mean: mean cable length [m] (default ``10·L_short_mean``).
        long_axis: in-plane cable orientation (3- or 2-vector; default +x).
        long_axis_jitter_deg: Gaussian angular spread of the cables about the
            long axis [deg]. Default 0.0 (perfectly aligned cables — no free
            parameter; raise it only if a measured SF dispersion is supplied).
        band_thickness: basal z-band depth [m] for radial separation (KU-3.17
            cortex thickness analogue; default 200 nm — a discretisation knob).
        seed / rng: RNG control.
        n_beads_min / n_beads_max: per-filament bead-count clamps.

    Returns:
        A :class:`VariableLengthCortexLayout` (planar; ``is_formin`` marks cables).
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    F = int(n_filaments)
    if F <= 0:
        raise ValueError(f"n_filaments must be > 0; got {F}")
    _require_finite_positive("ell0", ell0)
    _require_finite_positive("footprint_radius", footprint_radius)
    _require_finite_positive("band_thickness", band_thickness)
    if not math.isfinite(z_basal):
        raise ValueError(f"z_basal must be finite; got {z_basal!r}")
    if not (0.0 <= formin_fraction <= 1.0):
        raise ValueError(f"formin_fraction must be in [0, 1]; got {formin_fraction}")
    if n_beads_min < 2:
        raise ValueError(f"n_beads_min must be ≥ 2; got {n_beads_min}")

    L0 = float(ell0)
    if L_short_mean is None:
        L_short_mean = L0
    if L_long_mean is None:
        L_long_mean = 10.0 * L_short_mean
    _require_finite_positive("L_short_mean", L_short_mean)
    if not (math.isfinite(L_long_mean) and L_long_mean >= L_short_mean):
        raise ValueError(
            f"L_long_mean must be finite ≥ L_short_mean; got {L_long_mean}"
        )

    if n_beads_max is None:
        # A cable chord cannot exceed the footprint diameter.
        L_cap = min(2.0 * footprint_radius, max(L_long_mean * 6.0, 12.0 * L0))
        n_beads_max = max(n_beads_min, int(round(L_cap / L0)) + 1)

    # ---- per-filament class (cable=formin) + exponential length ----
    is_formin = rng.random(F) < formin_fraction
    means = np.where(is_formin, L_long_mean, L_short_mean)
    L_samples = rng.exponential(means)
    n_beads_per = np.clip(
        np.round(L_samples / L0).astype(np.int64) + 1, n_beads_min, n_beads_max
    )
    L_realised = (n_beads_per - 1).astype(np.float64) * L0

    # ---- in-plane disk CoM (uniform-area) + per-filament band depth ----
    u = rng.uniform(0.0, 1.0, F)
    r_disk = footprint_radius * np.sqrt(u)            # uniform-area disk sampling
    theta = rng.uniform(0.0, 2.0 * math.pi, F)
    depth = rng.uniform(0.0, band_thickness, F)        # toward the substrate
    centers = np.empty((F, 3), dtype=np.float64)
    centers[:, 0] = r_disk * np.cos(theta)
    centers[:, 1] = r_disk * np.sin(theta)
    centers[:, 2] = z_basal - depth

    # ---- in-plane tangents: cables along the long axis, infill isotropic ----
    ax2 = _inplane_long_axis(long_axis)
    base_ang = math.atan2(ax2[1], ax2[0])
    angles = np.empty(F, dtype=np.float64)
    short_mask = ~is_formin
    n_short = int(short_mask.sum())
    if n_short > 0:
        angles[short_mask] = rng.uniform(0.0, 2.0 * math.pi, n_short)
    n_long = int(is_formin.sum())
    if n_long > 0:
        jit = math.radians(long_axis_jitter_deg)
        # ± along the axis (a cable has no head/tail polarity here), with optional
        # Gaussian jitter. jitter 0 → perfectly aligned cables (no free param).
        sign = rng.choice([0.0, math.pi], size=n_long)
        noise = rng.normal(0.0, jit, n_long) if jit > 0.0 else np.zeros(n_long)
        angles[is_formin] = base_ang + sign + noise
    tangents = np.stack(
        [np.cos(angles), np.sin(angles), np.zeros(F)], axis=1
    )
    tangents = tangents / np.linalg.norm(
        tangents, axis=1, keepdims=True
    ).clip(min=1e-30)

    # ---- flat layout (planar; z constant within a filament) ----
    n_total_beads = int(n_beads_per.sum())
    filament_starts = np.zeros(F, dtype=np.int64)
    filament_starts[1:] = np.cumsum(n_beads_per[:-1])
    positions_flat = np.empty((n_total_beads, 3), dtype=np.float64)
    bond_pairs: list[tuple[int, int]] = []
    angle_triplets: list[tuple[int, int, int]] = []
    for f in range(F):
        N_f = int(n_beads_per[f])
        start = int(filament_starts[f])
        offsets = (np.arange(N_f, dtype=np.float64) - 0.5 * (N_f - 1)) * L0
        positions_flat[start:start + N_f] = (
            centers[f] + offsets[:, None] * tangents[f]
        )
        for j in range(N_f - 1):
            bond_pairs.append((start + j, start + j + 1))
        for j in range(N_f - 2):
            angle_triplets.append((start + j, start + j + 1, start + j + 2))

    bond_groups = np.array(bond_pairs, dtype=np.int64).reshape(-1, 2)
    angle_groups = (
        np.array(angle_triplets, dtype=np.int64).reshape(-1, 3)
        if angle_triplets else np.empty((0, 3), dtype=np.int64)
    )

    return VariableLengthCortexLayout(
        positions_flat=positions_flat,
        n_beads_per_filament=n_beads_per,
        filament_starts=filament_starts,
        L_per_filament=L_realised,
        centers_of_mass=centers,
        tangents=tangents,
        bond_groups=bond_groups,
        angle_groups=angle_groups,
        is_formin=is_formin,
    )


def basal_mesh_build_report(
    layout: VariableLengthCortexLayout,
    *,
    ell0: float,
    footprint_radius: float,
    z_basal: float,
    band_thickness: float,
    long_axis: np.ndarray | None = None,
    align_floor: float = 0.8,
) -> dict[str, Any]:
    """Build-time B1 gate metrics (PLANAR + force-free + long-axis cables).

    Returns a dict with the controls and a ``verdict`` ("PASS"/"REVIEW"):
      * planar_ok          — every bead z ∈ [z_basal − band_thickness, z_basal].
      * within_footprint   — every CoM radial (x,y) ≤ footprint_radius (+ε).
      * force_free         — every chain bond length == ell0 (max strain < 1e-9).
      * cables_aligned     — mean |cos| of LONG (cable) tangents vs long axis ≥ floor.
      * infill_isotropic   — mean |cos| of SHORT (infill) tangents ≈ 2/π (in [0.5,0.75]).
    """
    pos = np.asarray(layout.positions_flat, dtype=np.float64)
    z = pos[:, 2]
    planar_ok = bool(
        np.all(z <= z_basal + 1e-12)
        and np.all(z >= z_basal - band_thickness - 1e-12)
    )

    com_xy = np.asarray(layout.centers_of_mass, dtype=np.float64)[:, :2]
    radial = np.linalg.norm(com_xy, axis=1)
    within_footprint = bool(np.all(radial <= footprint_radius + 1e-12))

    # force-free: chain bond lengths == ell0.
    bg = np.asarray(layout.bond_groups, dtype=np.int64).reshape(-1, 2)
    if bg.shape[0] > 0:
        seg = pos[bg[:, 0]] - pos[bg[:, 1]]
        ln = np.linalg.norm(seg, axis=1)
        max_strain = float(np.max(np.abs(ln - ell0) / ell0))
    else:
        max_strain = 0.0
    force_free = bool(max_strain < 1e-9)

    ax2 = _inplane_long_axis(long_axis)
    tang = np.asarray(layout.tangents, dtype=np.float64)[:, :2]
    tnrm = tang / np.linalg.norm(tang, axis=1, keepdims=True).clip(min=1e-30)
    absdot = np.abs(tnrm @ ax2)
    is_formin = np.asarray(layout.is_formin, dtype=bool)
    cable_align = float(np.mean(absdot[is_formin])) if is_formin.any() else float("nan")
    infill_align = (
        float(np.mean(absdot[~is_formin])) if (~is_formin).any() else float("nan")
    )
    cables_aligned = bool(np.isfinite(cable_align) and cable_align >= align_floor)
    infill_isotropic = bool(
        np.isfinite(infill_align) and 0.5 <= infill_align <= 0.75
    )

    n_fil = int(layout.n_beads_per_filament.shape[0])
    controls = {
        "n_filaments": n_fil,
        "n_cables": int(is_formin.sum()),
        "n_infill": int((~is_formin).sum()),
        "n_beads_total": int(pos.shape[0]),
        "planar_ok": planar_ok,
        "within_footprint": within_footprint,
        "force_free": {"ok": force_free, "max_strain": max_strain},
        "cables_aligned": {"ok": cables_aligned, "mean_abs_cos": cable_align},
        "infill_isotropic": {"ok": infill_isotropic, "mean_abs_cos": infill_align},
    }
    ok = (
        planar_ok and within_footprint and force_free
        and cables_aligned and infill_isotropic
    )
    return {"verdict": "PASS" if ok else "REVIEW", "controls": controls}
