"""H.5 Lamellipodium (Bieling/Funk greenfield) — dendritic actin meshwork
at the leading-edge WAVE membrane plane.

Phase 1 H.5 brief `ffn_sim/docs/briefs/H5_lamellipodium.md`,
PHASE_0_3_DECISIONS.md §D1, AFINES_ALGORITHM_NOTES.md §5.

Implements three batched D2 Updaters per the D1 spec:

* :class:`BarbedEndElongationUpdater` — per-barbed-end Bell-Evans
  slip elongation ``k_elong(F) = k_elong⁰ · exp(−F δ_elong / kT)``
  with ``k_elong⁰ = 11.6 s⁻¹``, ``δ_elong = 2.7 nm`` (Bieling 2016).
  On fire, appends an actin bead at rest length ``ℓ_0`` along the
  filament's local tangent direction (mother grows away from WAVE,
  daughters grow at 72° off their mother).

* :class:`ArpBranchingUpdater` — per-WAVE Arp2/3 daughter nucleation
  ``k_b(F) = k_b⁰ · (1 − 0.2 · F / F_stall)`` (Bieling 2016) with
  ``k_b⁰ = 0.037 s⁻¹`` per WAVE molecule.  On fire, inserts a
  DAUGHTER actin bead bonded to the mother's barbed-end via a
  branch-angle harmonic ``md.angle.Harmonic`` at ``t0 = 72°``
  (Arp2/3 crystal-structure angle).  Funk 2021 Nat Commun abortive-failure
  branch emerges naturally: above the per-WAVE force threshold
  ``F > 500 Pa · WAVE_area`` the ``(1 − 0.2 F/F_stall)`` factor
  becomes negative and is clamped to zero — no separate scalar
  override.

* :class:`CappingUpdater` — per-barbed-end Bell-Evans
  ``k_cap(F) = k_cap⁰ · exp(−F · δ_cap · sin θ / kT)`` with
  ``k_cap⁰ = 3 s⁻¹``, ``δ_cap = 0.3 nm`` (Li/Bieling 2022; value
  PI-pending per audit: may be 2.7 nm = δ_elong, Li/Bieling 2022).  On fire,
  marks the barbed end CAPPED — no further elongation or
  branching; the actin chain is preserved.

Topology
--------
* WAVE plane at ``y = +Y_max`` (top of box; opposite the cortex shell
  which sits at the box origin).  ``n_WAVE = σ_NPF · WAVE_area`` particles.
* WAVE particles pinned to the plane via a radial-harmonic
  ``md.force.Custom`` field (same pattern as ``cortex/erm.py``,
  with CFL gate against the host integrator's ``dt``).
* One mother actin seed per WAVE at construction (chosen per H.5
  brief open question: "(a) seeded near WAVE plane at construction").
  Mother grows from WAVE INTO cytosol (barbed end pointed away from
  membrane, toward ``−y``).

Particle types
--------------
* ``wave_particle`` — membrane WAVE/NPF reservoir.  Pinned at plane.
* ``actin_lamel`` — lamellipodial actin bead.  Per-bead state tracks
  whether it is a BARBED-END (eligible for elongation / branching /
  capping) or INTERIOR (frozen mid-chain) bead — same Python-level
  bookkeeping pattern as `crosslinkers.XlinkBondUpdater._head_bound_to_actin`.

Bond types
----------
* ``lamel_actin_bond`` — actin-actin harmonic backbone.
* ``lamel_branch_bond`` — Arp2/3 mother-to-daughter branch bond.
* ``lamel_wave_anchor`` — WAVE-to-mother-bead initial anchor (optional;
  pinned via the WAVE membrane plane force compute and the harmonic
  bond keeps mother attached to its WAVE seed).

Angle types
-----------
* ``lamel_branch_angle`` — Arp2/3 branch angle harmonic at ``t0 = 72°``,
  ``k_angle = 100 pN·μm/rad²`` (literature TBD, brief default).

Sanity Gate
-----------
*Per CLAUDE.md hard rule.  STATIC checks in
``ffn_sim/tests/test_lamellipodium.py``.*

1. **Dimensional analysis**
   - All rates [1/s]; forces [N]; ``δ_elong``, ``δ_cap`` [m]; ``kT`` [J].
   - Bell-Evans exponent ``F · δ / kT`` dimensionless. ✓
   - ``k_b(F) = k_b⁰ · (1 − 0.2 F/F_stall)`` linear; ``F_stall`` [N];
     dimensionless ratio in parens. ✓
   - ``WAVE_area`` [m²] · ``500 Pa`` → [N] = abortive-branch threshold per WAVE.
   - Branch angle ``t0 = 72°`` = ``5π/12`` rad.  ``k_angle`` [J/rad²].
   - ``dt_batch · k_max`` dimensionless → D2 batch CFL.

2. **Boundary cases**
   - ``n_WAVE == 0``: empty membrane plane, all Updaters idle.
   - ``F > F_stall / 0.2`` in branching rate: ``k_b`` would go negative;
     clamp to zero (Funk 2021 Nat Commun abortive-failure regime emerges naturally).
   - Membrane plane outside box: validated via ``Y_max ≤ L_box / 2``.
   - All barbed ends capped: elongation + branching idle; capping idle
     (no eligible ends).

3. **Conservation invariants**
   - Total bond count grows monotonically (no bond removal in
     lamellipodium — capping mutates per-bead state, not topology).
   - Per-WAVE mother count is invariant at 1 (one mother per WAVE seed,
     daughters branch off in the actin tree without consuming the mother).
   - Capped-bead count monotonically increases over time.

4. **Numerical sanity**
   - All ``np.isfinite`` checks at resolve time.
   - D2 batch CFL ``batch_steps · dt · max(k_elong, k_b, k_cap) ≤ 1e-3``
     enforced; auto-shrinks ``batch_steps`` per the cortex `XlinkBondUpdater` precedent.
   - WAVE plane pinning ``k_wave_pin`` must satisfy ERM-style CFL
     ``dt ≤ 0.1 · γ_b / k_wave_pin``; gate raises otherwise.

5. **Sign / sense**
   - Elongation: barbed end advances along its LOCAL TANGENT direction
     (mother: ``−ŷ`` initial; daughter: 72° off its mother).  STATIC test.
   - Branching: daughter inserted at mother's barbed-end + ``ℓ_0`` along
     a tangent rotated 72° from the mother's tangent (random azimuthal
     angle in the plane perpendicular to the mother).
   - Capping: pure state mutation; no force-direction sign concern.

6. **Measurement protocol**
   - Dendritic density: per-WAVE barbed-end count averaged over the
     lamellipodium volume.  KU-5.1 target ≈ 100 / μm².
   - Force-velocity: membrane-recession rate under imposed cortex-side
     load.  KU-5.2 Bieling 2016 reference curve.
   - Abortive-branching: branching-event count per WAVE per second,
     plotted vs F.  KU-5.3 Funk 2021 Nat Commun reference: 50 % drop at 500 Pa
     · WAVE_area.

References
----------
- Brief: `ffn_sim/docs/briefs/H5_lamellipodium.md`.
- PHASE_0_3_DECISIONS.md §D1 (greenfield Arp2/3 Bieling/Funk).
- Bieling 2016 Cell (force-velocity + dendritic density).
- Funk 2021 Nat Commun (abortive branching/NPF); Li/Bieling 2022 eLife (capping Bell-Evans).
- Mullins 1998 Nature (Arp2/3 72° crystal structure).
- ``ffn_sim/cortex/crosslinkers.py`` (D2 batched Updater pattern).
- ``ffn_sim/cortex/erm.py`` (md.force.Custom membrane confinement pattern).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedH5:
    """KU-5.x lamellipodium parameters + Bieling/Funk D1 rates."""

    # WAVE membrane plane
    Y_max: float                  # m   membrane plane y-coordinate (top of box)
    n_WAVE: int                   # WAVE particle count on the plane
    wave_area: float              # m²   WAVE plane area (= L_box² if full plane)
    sigma_NPF: float              # 1/m² emergent areal density of WAVE
    k_wave_pin: float             # N/m harmonic pin of WAVE to plane

    # Actin per-filament force constants (mirrors H.2 single-actin-bundle)
    bead_radius: float            # m   actin bead radius (30 nm, ×40 bundle)
    rest_length: float            # m   ℓ_0 backbone bond rest length (0.5 μm)
    bond_k: float                 # N/m actin-actin harmonic spring (= μ/ℓ_0)
    angle_branch_t0: float        # rad branch angle (72° = 5π/12)
    angle_branch_k: float         # J/rad² Arp2/3 angle stiffness (≈ 100 pN·μm/rad²)

    # D1 elongation (Bieling 2016 slip Bell-Evans)
    k_elong_0: float              # 1/s 11.6 unloaded per barbed end
    delta_elong: float            # m   2.7 nm G-actin attachment length
    F_stall_elong: float          # N   per-barbed-end stall force (Bieling)

    # D1 Arp2/3 branching (Bieling 2016 + Funk 2021 mechanistic coupling, γ-Phase 1)
    # PI-ratified 2026-05-29: branching rate is gated by the CP-NPF release
    # mechanism (rate ∝ n_capped / n_total = "free-NPF fraction"), and the
    # geometry searches ANY actin_lamel bead within r_branch_eff (Arp2/3
    # binds along F-actin side, not only at the barbed tip).
    # See docs/briefs/H5_GAMMA_BRANCHING_DESIGN.md (commit 62d6fbc) §3.
    k_b_0: float                  # 1/s 0.037 per WAVE unloaded (Bieling 2016)
    F_stall_branch: float         # N   per-WAVE branching stall (Bieling)
    r_branch_eff: float           # m   WAVE→any actin_lamel reach (PI 2026-05-29: 100 nm,
                                  #     ≈ Arp2/3 physical reach: Arp2/3 ⌀22nm + F-actin ⌀7nm
                                  #     + access window. Funk 2021 + Bieling 2023 review.)
    abortive_pressure_Pa: float   # 500 Pa Funk 2021 Nat Commun threshold

    # D1 capping (Li/Bieling 2022, eLife 11:e73145)
    k_cap_0: float                # 1/s 3 at 100 nM CP
    delta_cap: float              # m   Li/Bieling 2022 ≈ 0.3 nm equiv (we store as length)

    kT: float                     # J
    seed: int

    # Box / integration
    L_box: float                  # m   host box edge
    dt: float                     # s   host integrator dt

    # D2 batch
    batch_steps: int              # BAOAB steps between Updater ticks

    # Acceptance (first-principles bands)
    dendritic_density_target: float  # 1/m² KU-5.1 Bieling ~100 /μm²
    force_velocity_rel_tol: float    # ±30 % KU-5.2 Bieling
    abortive_drop_min: float         # 0.5 KU-5.3 Funk 50 % drop at 500 Pa

    # Derived
    batch_dt: float = 0.0
    k_max_for_cfl: float = 0.0
    F_abortive_per_wave: float = 0.0
    extras: dict[str, Any] = field(default_factory=dict)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def resolve_h5_lamellipodium(
    cfg: dict, *, L_box: float, dt: float, kT: float = 4.28e-21,
) -> ResolvedH5:
    """Resolve cortex.h5 / lamellipodium config block with D1 + D2 gates.

    ``cfg`` is the YAML root or the ``lamellipodium`` sub-dict.  ``L_box``
    and ``dt`` come from the host cortex resolved params (lamellipodium
    inherits the integration dt and the box geometry).
    """
    if "lamellipodium" in cfg:
        cfg = cfg["lamellipodium"]

    # Y_max default = L_box / 2 − margin (top of box, leaving a sliver).
    Y_max = float(cfg.get("Y_max", 0.45 * L_box))
    wave_area = float(cfg.get("wave_area", L_box * L_box))
    sigma_NPF = float(cfg["sigma_NPF"])
    n_WAVE = int(cfg.get("n_WAVE", int(round(sigma_NPF * wave_area))))

    p = ResolvedH5(
        Y_max=Y_max,
        n_WAVE=n_WAVE,
        wave_area=wave_area,
        sigma_NPF=sigma_NPF,
        k_wave_pin=float(cfg["k_wave_pin"]),
        bead_radius=float(cfg["bead_radius"]),
        rest_length=float(cfg["rest_length"]),
        bond_k=float(cfg["bond_k"]),
        angle_branch_t0=float(cfg.get("angle_branch_t0", 5.0 * math.pi / 12.0)),
        angle_branch_k=float(cfg["angle_branch_k"]),
        k_elong_0=float(cfg["k_elong_0"]),
        delta_elong=float(cfg["delta_elong"]),
        F_stall_elong=float(cfg["F_stall_elong"]),
        k_b_0=float(cfg["k_b_0"]),
        F_stall_branch=float(cfg["F_stall_branch"]),
        r_branch_eff=float(cfg["r_branch_eff"]),
        abortive_pressure_Pa=float(cfg.get("abortive_pressure_Pa", 500.0)),
        k_cap_0=float(cfg["k_cap_0"]),
        delta_cap=float(cfg["delta_cap"]),
        kT=float(kT),
        seed=int(cfg.get("seed", 45)),
        L_box=float(L_box),
        dt=float(dt),
        batch_steps=int(cfg.get("batch_steps", 100)),
        dendritic_density_target=float(cfg["acceptance"]["dendritic_density_target"]),
        force_velocity_rel_tol=float(cfg["acceptance"]["force_velocity_rel_tol"]),
        abortive_drop_min=float(cfg["acceptance"]["abortive_drop_min"]),
    )

    # ---- §2 boundary checks ----
    for name, x in [
        ("Y_max", abs(p.Y_max)), ("wave_area", p.wave_area),
        ("sigma_NPF", p.sigma_NPF), ("k_wave_pin", p.k_wave_pin),
        ("bead_radius", p.bead_radius), ("rest_length", p.rest_length),
        ("bond_k", p.bond_k), ("angle_branch_k", p.angle_branch_k),
        ("k_elong_0", p.k_elong_0), ("delta_elong", p.delta_elong),
        ("F_stall_elong", p.F_stall_elong), ("k_b_0", p.k_b_0),
        ("F_stall_branch", p.F_stall_branch), ("r_branch_eff", p.r_branch_eff),
        ("k_cap_0", p.k_cap_0), ("delta_cap", p.delta_cap),
        ("kT", p.kT), ("L_box", p.L_box), ("dt", p.dt),
    ]:
        _require_finite_positive(name, x)
    if p.n_WAVE < 0:
        raise ValueError(f"n_WAVE must be ≥ 0; got {p.n_WAVE}")
    if p.batch_steps < 1:
        raise ValueError(f"batch_steps must be ≥ 1; got {p.batch_steps}")
    if abs(p.Y_max) >= 0.5 * p.L_box:
        raise ValueError(
            f"Y_max = {p.Y_max:.3e} m must be inside box half-edge "
            f"L_box/2 = {0.5 * p.L_box:.3e} m."
        )

    # ---- §1 derived: D2 batch CFL + Funk abortive threshold ----
    p.batch_dt = p.batch_steps * p.dt
    p.k_max_for_cfl = max(p.k_elong_0, p.k_b_0, p.k_cap_0)
    cfl_product = p.batch_dt * p.k_max_for_cfl
    if cfl_product > 1.0e-3:
        target_batch_dt = 1.0e-3 / p.k_max_for_cfl
        new_batch_steps = max(1, int(math.floor(target_batch_dt / p.dt)))
        p.batch_steps = new_batch_steps
        p.batch_dt = p.batch_steps * p.dt
        p.extras["batch_steps_shrunk_from"] = int(cfg.get("batch_steps", 100))

    # Per-WAVE Funk abortive force threshold: F = pressure · area_per_WAVE.
    area_per_wave = p.wave_area / max(p.n_WAVE, 1)
    p.F_abortive_per_wave = p.abortive_pressure_Pa * area_per_wave

    return p


# ---------------------------------------------------------------------------
# WAVE membrane plane + seed actin topology
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class LamellipodiumLayout:
    """Per-WAVE membrane + initial mother actin layout (construction state).

    Attributes
    ----------
    wave_positions : (n_WAVE, 3)
    wave_tag_start : int — global tag offset for the first WAVE particle.
    mother_seed_positions : (n_WAVE, 3) — initial mother actin bead per WAVE
        (placed at ``Y_max − ℓ_0`` directly below its WAVE).
    mother_tag_start : int — global tag offset for the first mother actin.
    wave_to_mother_bond_pairs : (n_WAVE, 2) — (WAVE_tag, mother_actin_tag)
        anchor-bond pairs (kept at rest length ℓ_0 with k_wave_pin/10 — soft
        so it does not break the WAVE-pin CFL).
    """

    wave_positions: np.ndarray
    wave_tag_start: int
    mother_seed_positions: np.ndarray
    mother_tag_start: int
    wave_to_mother_bond_pairs: np.ndarray
    # Per-WAVE initial mother-actin barbed-end unit tangent (n_WAVE, 3). None for
    # the flat_plane geometry (seed loops then default to the legacy -ŷ). The
    # basal_ring / polarized_patch geometries set it to outward-radial / forward
    # directions so the seed loops grow the dendritic network in-geometry.
    mother_tangents: np.ndarray | None = None


def generate_lamellipodium_layout(
    p: ResolvedH5,
    *,
    wave_tag_start: int,
    rng: np.random.Generator | None = None,
    geometry: str = "flat_plane",
    R_cell: float | None = None,
    cap_depth: float | None = None,
    contact_radius_frac: float | None = None,
    basal_z_offset: float | None = None,
    polarization: "np.ndarray | tuple[float, float, float] | None" = None,
    half_angle_azimuth: float = math.pi / 6.0,
    half_angle_linear: float | None = None,
    fan_spread: float | None = None,
) -> LamellipodiumLayout:
    """Place WAVE particles + one mother-actin seed per WAVE.

    ``geometry`` selects the leading-edge placement (H.7 active-spreading):

    * ``"flat_plane"`` (default) — the original in-vitro reconstitution geometry:
      WAVE uniform on the ``y = Y_max`` plane, mothers one ℓ₀ below, growing −ŷ.
      Bit-for-bit unchanged from the pre-H.7 builder.
    * ``"basal_ring"`` — single-cell ISOTROPIC spreading: WAVE on the basal
      contact ring, mothers growing outward-radially in the basal plane
      (:func:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium_basal_ring.generate_basal_ring_lamellipodium_layout`).
    * ``"polarized_patch"`` — single-cell MIGRATING: a localized leading-edge
      patch about a polarization ``p̂``, mothers growing forward
      (:func:`ffn_sim.archive.hoomd_legacy.cell.lamellipodium_polarized_patch.generate_polarized_patch_layout`).

    The two single-cell geometries require ``R_cell`` (the cortex shell radius)
    and return a drop-in :class:`LamellipodiumLayout` with ``mother_tangents``
    populated. The Bieling/Funk dendritic mechanism is unchanged — only the
    leading-edge geometry differs (membrane load is OFF for step-1).
    """
    if geometry != "flat_plane":
        if R_cell is None:
            raise ValueError(
                f"geometry={geometry!r} requires R_cell (the cortex shell radius)."
            )
        if geometry == "basal_ring":
            from ffn_sim.archive.hoomd_legacy.cell.lamellipodium_basal_ring import (
                generate_basal_ring_lamellipodium_layout,
            )
            wrapped = generate_basal_ring_lamellipodium_layout(
                p, wave_tag_start=wave_tag_start, R_cell=R_cell,
                cap_depth=cap_depth, contact_radius_frac=contact_radius_frac,
                basal_z_offset=basal_z_offset, rng=rng,
            )
        elif geometry == "polarized_patch":
            from ffn_sim.archive.hoomd_legacy.cell.lamellipodium_polarized_patch import (
                generate_polarized_patch_layout,
            )
            wrapped = generate_polarized_patch_layout(
                p, wave_tag_start=wave_tag_start, R_cell=R_cell,
                polarization=polarization, half_angle_azimuth=half_angle_azimuth,
                half_angle_linear=half_angle_linear, fan_spread=fan_spread, rng=rng,
            )
        else:
            raise ValueError(
                f"unknown geometry={geometry!r}; expected one of "
                "'flat_plane', 'basal_ring', 'polarized_patch'."
            )
        layout = wrapped.layout
        layout.mother_tangents = wrapped.mother_tangents
        return layout

    if rng is None:
        rng = np.random.default_rng(p.seed)

    # WAVE particle grid: uniform random (x, z) on the membrane plane.
    half = 0.5 * p.L_box
    n = p.n_WAVE
    if n == 0:
        return LamellipodiumLayout(
            wave_positions=np.empty((0, 3), dtype=np.float64),
            wave_tag_start=wave_tag_start,
            mother_seed_positions=np.empty((0, 3), dtype=np.float64),
            mother_tag_start=wave_tag_start,
            wave_to_mother_bond_pairs=np.empty((0, 2), dtype=np.int64),
        )

    xs = rng.uniform(-half * 0.95, half * 0.95, n)
    zs = rng.uniform(-half * 0.95, half * 0.95, n)
    wave_positions = np.column_stack([xs, np.full(n, p.Y_max), zs])

    # Mother seed: one actin bead at y = Y_max - ℓ_0 directly below WAVE.
    mother_positions = wave_positions.copy()
    mother_positions[:, 1] -= p.rest_length

    # Anchor bonds: (WAVE_tag, mother_actin_tag).  WAVE tags go first, then
    # mother actins.
    wave_tags = np.arange(n, dtype=np.int64) + wave_tag_start
    mother_tag_start = wave_tag_start + n
    mother_tags = np.arange(n, dtype=np.int64) + mother_tag_start
    anchor_pairs = np.column_stack([wave_tags, mother_tags])

    return LamellipodiumLayout(
        wave_positions=wave_positions,
        wave_tag_start=wave_tag_start,
        mother_seed_positions=mother_positions,
        mother_tag_start=mother_tag_start,
        wave_to_mother_bond_pairs=anchor_pairs,
    )


# ---------------------------------------------------------------------------
# WAVE membrane plane pin (md.force.Custom, mirrors cortex/erm.py)
# ---------------------------------------------------------------------------
class WaveMembranePin(md.force.Custom):
    """Harmonic pin on WAVE particles to the membrane plane ``y = Y_max``.

    ``U_i = ½ k_wave_pin (y_i − Y_max)²`` per WAVE particle.  Force only on
    particles with tag in ``[wave_tag_start, wave_tag_start + n_WAVE)``.
    """

    def __init__(
        self, p: ResolvedH5, *, wave_tag_start: int, n_WAVE: int,
        aniso: bool = False,
    ) -> None:
        super().__init__(aniso=aniso)
        self.p = p
        self.tag_start = int(wave_tag_start)
        self.tag_end = int(wave_tag_start + n_WAVE)

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        mask = (tag >= self.tag_start) & (tag < self.tag_end)
        dy = pos[:, 1] - self.p.Y_max
        F_vec = np.zeros_like(pos)
        F_vec[:, 1] = -self.p.k_wave_pin * dy
        U_per = 0.5 * self.p.k_wave_pin * dy * dy
        F_vec[~mask] = 0.0
        U_per[~mask] = 0.0
        with self.cpu_local_force_arrays as arrays:
            arrays.force[:] = F_vec
            arrays.potential_energy[:] = U_per


# ---------------------------------------------------------------------------
# Per-barbed-end / per-WAVE state — extends snapshot via Python dicts
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class LamellipodiumState:
    """Mutable runtime state — barbed-end status + branch tree bookkeeping.

    Attributes
    ----------
    barbed_end_tags : list[int]
        Global tags of CURRENTLY-eligible barbed-end actin beads
        (one per filament until capped).
    capped_tags : set[int]
        Global tags of capped barbed-end actin beads.
    parent_of : dict[int, int]
        For each daughter actin bead, the parent (mother barbed-end) tag.
    tangent_of : dict[int, np.ndarray]
        Per-barbed-end unit tangent vector (direction of next elongation).
    actin_next_tag : int
        Next available global tag for newly-elongated actin beads.
    """

    barbed_end_tags: list[int] = field(default_factory=list)
    capped_tags: set[int] = field(default_factory=set)
    parent_of: dict[int, int] = field(default_factory=dict)
    tangent_of: dict[int, np.ndarray] = field(default_factory=dict)
    actin_next_tag: int = 0


# ---------------------------------------------------------------------------
# Updaters — D2 batched (mirrors crosslinkers.XlinkBondUpdater pattern)
# ---------------------------------------------------------------------------
class _BatchedLamelUpdater(hoomd.custom.Action):
    """Shared bookkeeping for the three D2-batched lamellipodium Updaters."""

    def __init__(
        self, *, p: ResolvedH5, lamel_state: LamellipodiumState,
        seed_offset: int,
    ) -> None:
        super().__init__()
        self.p = p
        self.state = lamel_state
        self._rng = np.random.default_rng(p.seed + seed_offset)
        self._sim_ref: hoomd.Simulation | None = None
        self._steps_run = 0

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    @property
    def steps_run(self) -> int:
        return self._steps_run


class BarbedEndElongationUpdater(_BatchedLamelUpdater):
    """D1 elongation: ``k_elong(F) = k_elong⁰ · exp(−F · δ_elong / kT)``.

    Per batch tick, for every eligible (non-capped) barbed-end tag:

    1. Compute per-barbed-end load ``F = F_partitioned`` (force partition:
       network load divided by N_free_barbed_ends per the D1 spec).
    2. Sample ``p_elong = 1 − exp(−k_elong(F) · Δt_batch)``.
    3. On fire: append one actin bead at rest length ``ℓ_0`` along the
       barbed end's local tangent, demote the old barbed end to interior,
       promote the new bead to barbed-end status.

    Force partition uses a constant network-load handle ``F_network_total``
    set via ``set_network_load(F)`` from the test harness or Cell-level
    diagnostic.  Default = 0 (unloaded).
    """

    def __init__(
        self, *, p: ResolvedH5, lamel_state: LamellipodiumState,
        seed_offset: int = 4,
    ) -> None:
        super().__init__(p=p, lamel_state=lamel_state, seed_offset=seed_offset)
        self._F_network_total = 0.0
        self._n_elongation_events = 0

    def set_network_load(self, F_total: float) -> None:
        """Set the network-level load distributed across free barbed ends."""
        self._F_network_total = float(F_total)

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        assert sim is not None
        if not self.state.barbed_end_tags:
            self._steps_run += 1
            return

        read_snap = sim.state.get_snapshot()
        if read_snap.communicator.rank != 0:
            return

        pos = np.asarray(read_snap.particles.position, dtype=np.float64).copy()
        n_part = int(read_snap.particles.N)

        free_be = [
            t for t in self.state.barbed_end_tags if t not in self.state.capped_tags
        ]
        if not free_be:
            self._steps_run += 1
            return

        # D1 force partition: F per barbed end = F_total / N_free.
        F_per = self._F_network_total / max(len(free_be), 1)
        bell_exponent = -F_per * self.p.delta_elong / self.p.kT
        k_elong = self.p.k_elong_0 * math.exp(bell_exponent)
        p_elong = 1.0 - math.exp(-k_elong * self.p.batch_dt)

        # Sample which barbed ends fire.
        u = self._rng.uniform(0.0, 1.0, size=len(free_be))
        fires = [be for be, ui in zip(free_be, u) if ui < p_elong]
        if not fires:
            self._steps_run += 1
            return

        # Build new actin beads on snapshot.  Each firing barbed end
        # appends ONE new bead at rest_length along its tangent.
        new_positions = []
        new_bonds = []
        new_be_tags = []
        for be_tag in fires:
            if be_tag not in self.state.tangent_of:
                continue
            tangent = self.state.tangent_of[be_tag]
            r_be = pos[be_tag]
            r_new = r_be + self.p.rest_length * tangent
            new_tag = self.state.actin_next_tag
            self.state.actin_next_tag += 1
            new_positions.append(r_new)
            new_bonds.append((be_tag, new_tag))
            new_be_tags.append(new_tag)
            # Move barbed-end status from old → new bead.
            self.state.parent_of[new_tag] = be_tag
            self.state.tangent_of[new_tag] = tangent.copy()
            # Remove old barbed end from eligibility (becomes interior).
            try:
                self.state.barbed_end_tags.remove(be_tag)
            except ValueError:
                pass
            self.state.barbed_end_tags.append(new_tag)

        if not new_positions:
            self._steps_run += 1
            return

        # Append new particles + bonds to snapshot.
        write_snap = _extend_snapshot_with_new_actins(
            read_snap, new_positions, new_bonds,
            actin_type_name="actin_lamel",
            bond_type_name="lamel_actin_bond",
        )
        sim.state.set_snapshot(write_snap)
        self._n_elongation_events += len(fires)
        self._steps_run += 1

    @property
    def n_elongation_events(self) -> int:
        return self._n_elongation_events


class ArpBranchingUpdater(_BatchedLamelUpdater):
    """D1 Arp2/3 daughter nucleation (γ-Phase 1, PI-ratified 2026-05-29).

    Funk 2021 mechanistic coupling — branching rate is gated by the
    "free-NPF fraction" (CP-mediated release of NPF sequestered on
    barbed ends). Brief implementation pre-γ treated branching and
    capping as independent and located branching only at the barbed
    tip; both were incorrect per Funk 2021 + Bieling 2016 + Abraham
    1999 evidence (KU-5.1 v1 production sweep, 2026-05-29 Lead session).

    Rate (per batch tick, per WAVE)::

        free_npf_fraction = n_capped / max(n_barbed + n_capped, 1)
        k_b_eff = k_b⁰ · (1 − 0.2 · F / F_stall_branch) · free_npf_fraction
        p_branch = 1 − exp(−k_b_eff · Δt_batch)

    - When ``n_capped = 0`` (initial): ``k_b_eff = 0`` → no branching
      (matches Funk 2021 sequestered-NPF regime).
    - When ``n_capped / n_total → 1`` (all ends capped):
      ``k_b_eff → k_b⁰ · force_factor`` (asymptotic free-NPF regime).

    Geometry (per WAVE that fires):

    1. Enumerate ALL ``actin_lamel`` particles in the snapshot
       (not only ``barbed_end_tags``) — Arp2/3 binds along F-actin
       side, not only at the barbed tip.
    2. Filter to particles within ``r_branch_eff`` of the WAVE
       (PI-ratified 2026-05-29: 100 nm, ≈ Arp2/3 physical reach).
    3. Pick mother from candidates weighted by 1/r² (diffusion-limited
       kinetics, PI-ratified 2026-05-29).
    4. Insert DAUGHTER actin bead at ``r_mother + ℓ₀ · t_daughter``
       where ``t_daughter`` is the mother's local tangent rotated 72°
       around a random perpendicular axis (Arp2/3 crystal angle).
    5. Add ``lamel_branch_bond`` (mother→daughter) + (if mother has
       a parent) ``lamel_branch_angle`` harmonic.
    6. Daughter starts as a new barbed end.

    Funk 2021 Nat Commun abortive emerges: when ``F > F_stall / 0.2 = 5 · F_stall``,
    the force factor ``(1 − 0.2 F / F_stall)`` clamps to 0 (no separate
    scalar override).

    Design: ``docs/briefs/H5_GAMMA_BRANCHING_DESIGN.md`` (commit 62d6fbc).
    """

    def __init__(
        self, *, p: ResolvedH5, lamel_state: LamellipodiumState,
        wave_tag_start: int, n_WAVE: int, seed_offset: int = 5,
    ) -> None:
        super().__init__(p=p, lamel_state=lamel_state, seed_offset=seed_offset)
        self.wave_tag_start = int(wave_tag_start)
        self.n_WAVE = int(n_WAVE)
        self._F_network_total = 0.0
        self._n_branch_events = 0

    def set_network_load(self, F_total: float) -> None:
        self._F_network_total = float(F_total)

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        assert sim is not None
        if self.n_WAVE == 0:
            self._steps_run += 1
            return

        # γ-Phase 1 rate factor (Funk 2021): free-NPF fraction. With
        # n_capped = 0 (start of simulation), branching is OFF until
        # capping has released NPFs from sequestered barbed ends.
        n_barbed = len(self.state.barbed_end_tags)
        n_capped = len(self.state.capped_tags)
        n_total = n_barbed + n_capped
        if n_total == 0:
            self._steps_run += 1
            return
        free_npf_fraction = n_capped / n_total
        if free_npf_fraction == 0.0:
            self._steps_run += 1
            return

        # Per-WAVE force allotment (Bieling 2016, unchanged from pre-γ).
        F_per_wave = self._F_network_total / max(self.n_WAVE, 1)
        ratio = 0.2 * F_per_wave / max(self.p.F_stall_branch, 1.0e-30)
        k_b_eff = max(0.0, self.p.k_b_0 * (1.0 - ratio)) * free_npf_fraction
        if k_b_eff == 0.0:
            self._steps_run += 1
            return
        p_branch = 1.0 - math.exp(-k_b_eff * self.p.batch_dt)

        read_snap = sim.state.get_snapshot()
        if read_snap.communicator.rank != 0:
            return
        pos = np.asarray(read_snap.particles.position, dtype=np.float64).copy()
        typeid = np.asarray(read_snap.particles.typeid)
        type_names = list(read_snap.particles.types)
        if "actin_lamel" not in type_names:
            self._steps_run += 1
            return
        actin_lamel_tid = type_names.index("actin_lamel")

        # The global (rank-0 aggregated) snapshot is tag-ordered:
        # `position[i]`, `typeid[i]` correspond to particle with tag=i. So
        # row index IS the tag here — no separate tag lookup needed. This
        # is different from the local snapshot (cpu_local_snapshot) which
        # is row-permuted by HOOMD's ParticleSorter.

        # γ-Phase 1: enumerate ALL actin_lamel beads (mother sites can be
        # anywhere along F-actin per Funk 2021 + Bieling 2023 review).
        actin_rows = np.where(typeid == actin_lamel_tid)[0]
        if actin_rows.size == 0:
            self._steps_run += 1
            return
        actin_pos = pos[actin_rows]
        actin_tags = actin_rows  # row == tag in global snapshot

        # Pairwise WAVE→actin distance. WAVE tags are contiguous starting
        # at wave_tag_start (set at Cell.build time); row == tag here.
        wave_rows = np.arange(
            self.wave_tag_start, self.wave_tag_start + self.n_WAVE,
        )
        wave_pos = pos[wave_rows]
        d = np.linalg.norm(
            wave_pos[:, None, :] - actin_pos[None, :, :], axis=-1
        )  # (n_WAVE, n_actin)
        eligible = d < self.p.r_branch_eff
        any_eligible = eligible.any(axis=1)

        u = self._rng.uniform(0.0, 1.0, size=self.n_WAVE)
        fires = (u < p_branch) & any_eligible
        if not fires.any():
            self._steps_run += 1
            return

        new_positions = []
        new_actin_bonds: list = []  # γ-Phase 1 has no backbone bonds, branch-only
        new_branch_bonds = []
        new_branch_angles = []
        for w in np.nonzero(fires)[0]:
            # γ-Phase 1: pick mother from all eligible actin_lamel within
            # r_branch_eff, weighted by 1/r² (diffusion-limited kinetics,
            # PI-ratified 2026-05-29). Floor on r² to avoid div-zero when
            # WAVE sits on top of an actin bead (numerical edge).
            elig_idx = np.where(eligible[w])[0]
            d_w = d[w, elig_idx]
            r2_floor = (0.1 * self.p.r_branch_eff) ** 2  # 10% of reach
            weights = 1.0 / np.maximum(d_w * d_w, r2_floor)
            weights /= weights.sum()
            mother_local = int(self._rng.choice(elig_idx.size, p=weights))
            mother_row = int(elig_idx[mother_local])
            mother_tag = int(actin_tags[mother_row])

            # Mother tangent lookup. Priority: tangent_of (set at construction
            # / elongation / branching) → derive from parent_of → fallback
            # random (only for orphan beads that should not exist).
            if mother_tag in self.state.tangent_of:
                t_mother = self.state.tangent_of[mother_tag]
            elif mother_tag in self.state.parent_of:
                # Global snapshot row == tag, so parent_tag indexes pos directly.
                parent_tag = int(self.state.parent_of[mother_tag])
                t_mother = pos[parent_tag] - actin_pos[mother_row]  # parent→mother
                norm = float(np.linalg.norm(t_mother))
                if norm < 1.0e-12:
                    continue  # degenerate; skip this WAVE this tick
                t_mother = t_mother / norm
            else:
                # Orphan — should not happen for well-formed assemblies.
                # Skip rather than fabricate a direction.
                continue

            # Daughter tangent: rotate mother tangent 72° around a random
            # perpendicular axis (Arp2/3 crystal, Mullins 1998).
            perp = _random_perpendicular(t_mother, self._rng)
            theta = self.p.angle_branch_t0
            t_daughter = math.cos(theta) * t_mother + math.sin(theta) * perp
            t_daughter /= np.linalg.norm(t_daughter)
            r_mother = actin_pos[mother_row]
            r_daughter = r_mother + self.p.rest_length * t_daughter

            new_tag = self.state.actin_next_tag
            self.state.actin_next_tag += 1
            new_positions.append(r_daughter)
            new_branch_bonds.append((mother_tag, new_tag))
            if mother_tag in self.state.parent_of:
                grand_parent = int(self.state.parent_of[mother_tag])
                new_branch_angles.append((grand_parent, mother_tag, new_tag))

            self.state.parent_of[new_tag] = mother_tag
            self.state.tangent_of[new_tag] = t_daughter.copy()
            # γ-Phase 1: daughter is a new barbed end. Mother stays whatever
            # it was (barbed, capped, or interior — Funk 2021 + Bieling 2023
            # both report branch site is NOT bound to mother's barbed-vs-other
            # state).
            self.state.barbed_end_tags.append(new_tag)

        if not new_positions:
            self._steps_run += 1
            return

        write_snap = _extend_snapshot_with_new_actins(
            read_snap, new_positions, new_actin_bonds,
            actin_type_name="actin_lamel",
            bond_type_name="lamel_actin_bond",
            extra_bond_groups=new_branch_bonds,
            extra_bond_type_name="lamel_branch_bond",
            extra_angle_groups=new_branch_angles,
            extra_angle_type_name="lamel_branch_angle",
        )
        sim.state.set_snapshot(write_snap)
        self._n_branch_events += len(new_positions)
        self._steps_run += 1

    @property
    def n_branch_events(self) -> int:
        return self._n_branch_events


class CappingUpdater(_BatchedLamelUpdater):
    """D1 capping: ``k_cap(F) = k_cap⁰ · exp(−F · δ_cap · sin θ / kT)``.

    Per batch tick, for every eligible barbed end:

    1. Compute per-barbed-end load ``F``.
    2. Compute ``sin θ`` between barbed-end tangent and membrane normal (+ŷ).
    3. Sample ``p_cap = 1 − exp(−k_cap(F) · Δt_batch)``.
    4. On fire: mark the barbed end CAPPED (move tag from
       ``barbed_end_tags`` to ``capped_tags``).  No topology mutation.
    """

    def __init__(
        self, *, p: ResolvedH5, lamel_state: LamellipodiumState,
        seed_offset: int = 6,
    ) -> None:
        super().__init__(p=p, lamel_state=lamel_state, seed_offset=seed_offset)
        self._F_network_total = 0.0
        self._n_cap_events = 0

    def set_network_load(self, F_total: float) -> None:
        self._F_network_total = float(F_total)

    def act(self, timestep: int) -> None:  # noqa: D401
        if not self.state.barbed_end_tags:
            self._steps_run += 1
            return

        free_be = [
            t for t in self.state.barbed_end_tags if t not in self.state.capped_tags
        ]
        if not free_be:
            self._steps_run += 1
            return

        F_per = self._F_network_total / max(len(free_be), 1)
        # Per-bead sin θ — assume membrane normal = +ŷ.  Tangent's y-
        # component gives cos θ; sin θ = √(1 − cos²θ).
        cap_events = []
        for be in free_be:
            tangent = self.state.tangent_of.get(be)
            if tangent is None:
                continue
            # STEP-2 TODO (H.7 single-cell geometries): the membrane normal is
            # hard-coded +ŷ here. For basal_ring / polarized_patch the leading-edge
            # normal is per-WAVE radial/forward, so this cos_theta is wrong for
            # those geometries. It is DORMANT while membrane-load is OFF (F_per=0
            # zeroes bell_exponent -> k_cap=k_cap_0, geometry-independent), so the
            # membrane-OFF step-1 A/A0 comparison is unaffected; it MUST be fixed
            # (pass the per-WAVE membrane normal) before membrane-load is wired on
            # for the curved geometries. (direction-review 2026-06-07.)
            cos_theta = float(tangent[1])  # tangent · ŷ  (flat-plane normal)
            sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
            bell_exponent = -F_per * self.p.delta_cap * sin_theta / self.p.kT
            k_cap = self.p.k_cap_0 * math.exp(bell_exponent)
            p_cap = 1.0 - math.exp(-k_cap * self.p.batch_dt)
            if self._rng.uniform() < p_cap:
                cap_events.append(be)

        for be in cap_events:
            self.state.capped_tags.add(be)
            try:
                self.state.barbed_end_tags.remove(be)
            except ValueError:
                pass

        self._n_cap_events += len(cap_events)
        self._steps_run += 1

    @property
    def n_cap_events(self) -> int:
        return self._n_cap_events


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _random_perpendicular(
    v: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Return a unit vector perpendicular to ``v``, uniform in the
    perpendicular plane (random azimuth)."""
    v = v / max(np.linalg.norm(v), 1e-30)
    # Build any vector not parallel to v.
    if abs(v[0]) < 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    else:
        ref = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(ref, v)
    e1 /= max(np.linalg.norm(e1), 1e-30)
    e2 = np.cross(v, e1)
    phi = rng.uniform(0.0, 2.0 * math.pi)
    return math.cos(phi) * e1 + math.sin(phi) * e2


def _extend_snapshot_with_new_actins(
    read_snap,
    new_positions: list,
    new_actin_bonds: list,
    *,
    actin_type_name: str,
    bond_type_name: str,
    extra_bond_groups: list | None = None,
    extra_bond_type_name: str | None = None,
    extra_angle_groups: list | None = None,
    extra_angle_type_name: str | None = None,
):
    """Build a fresh Snapshot extending the read snapshot with new actin
    beads + optional extra bonds (branch) + optional extra angles.

    Mirrors the cortex.crosslinkers `act()` snapshot-rebuild pattern so
    the runtime Updater dynamically grows the lamellipodium.
    """
    n_old = int(read_snap.particles.N)
    n_new = len(new_positions)
    n_total = n_old + n_new

    write_snap = hoomd.Snapshot()
    write_snap.particles.N = n_total
    # Particle types: ensure actin_type_name registered.
    types = list(read_snap.particles.types)
    if actin_type_name not in types:
        types.append(actin_type_name)
    actin_typeid = types.index(actin_type_name)
    write_snap.particles.types = types

    old_typeid = np.asarray(read_snap.particles.typeid)
    new_typeid = np.full(n_new, actin_typeid, dtype=np.uint32)
    write_snap.particles.typeid[:] = np.concatenate([old_typeid, new_typeid])

    old_pos = np.asarray(read_snap.particles.position)
    new_pos_arr = np.array(new_positions, dtype=np.float64).reshape(-1, 3)
    write_snap.particles.position[:] = np.concatenate([old_pos, new_pos_arr], axis=0)

    old_vel = np.asarray(read_snap.particles.velocity)
    new_vel = np.zeros((n_new, 3), dtype=np.float64)
    write_snap.particles.velocity[:] = np.concatenate([old_vel, new_vel], axis=0)

    old_mass = np.asarray(read_snap.particles.mass)
    new_mass = np.ones(n_new, dtype=np.float64)
    write_snap.particles.mass[:] = np.concatenate([old_mass, new_mass])

    old_image = np.asarray(read_snap.particles.image)
    new_image = np.zeros((n_new, 3), dtype=np.int32)
    write_snap.particles.image[:] = np.concatenate([old_image, new_image], axis=0)

    write_snap.configuration.box = list(read_snap.configuration.box)

    # Bonds: existing + actin backbone + optional branch.
    bond_types = list(read_snap.bonds.types)
    if bond_type_name not in bond_types:
        bond_types.append(bond_type_name)
    bond_typeid = bond_types.index(bond_type_name)
    if extra_bond_type_name is not None and extra_bond_type_name not in bond_types:
        bond_types.append(extra_bond_type_name)

    old_bg = np.asarray(read_snap.bonds.group, dtype=np.int64)
    old_bt = np.asarray(read_snap.bonds.typeid, dtype=np.uint32)
    new_actin_bg = (
        np.array(new_actin_bonds, dtype=np.int64).reshape(-1, 2)
        if new_actin_bonds else np.empty((0, 2), dtype=np.int64)
    )
    new_actin_bt = np.full(new_actin_bg.shape[0], bond_typeid, dtype=np.uint32)
    bg_pieces = [old_bg, new_actin_bg]
    bt_pieces = [old_bt, new_actin_bt]
    if extra_bond_groups and extra_bond_type_name is not None:
        extra_bg = np.array(extra_bond_groups, dtype=np.int64).reshape(-1, 2)
        extra_typeid = bond_types.index(extra_bond_type_name)
        extra_bt = np.full(extra_bg.shape[0], extra_typeid, dtype=np.uint32)
        bg_pieces.append(extra_bg)
        bt_pieces.append(extra_bt)
    merged_bg = np.concatenate(bg_pieces, axis=0).astype(np.uint32)
    merged_bt = np.concatenate(bt_pieces).astype(np.uint32)
    write_snap.bonds.N = int(merged_bg.shape[0])
    write_snap.bonds.types = bond_types
    if merged_bg.shape[0] > 0:
        write_snap.bonds.group[:] = merged_bg
        write_snap.bonds.typeid[:] = merged_bt

    # Angles: existing + optional branch.
    angle_types = list(read_snap.angles.types)
    old_ag = np.asarray(read_snap.angles.group, dtype=np.int64)
    old_at = np.asarray(read_snap.angles.typeid, dtype=np.uint32)
    if extra_angle_groups and extra_angle_type_name is not None:
        if extra_angle_type_name not in angle_types:
            angle_types.append(extra_angle_type_name)
        extra_ag = np.array(extra_angle_groups, dtype=np.int64).reshape(-1, 3)
        extra_typeid = angle_types.index(extra_angle_type_name)
        extra_at = np.full(extra_ag.shape[0], extra_typeid, dtype=np.uint32)
        merged_ag = np.concatenate([old_ag, extra_ag], axis=0).astype(np.uint32)
        merged_at = np.concatenate([old_at, extra_at]).astype(np.uint32)
    else:
        merged_ag = old_ag.astype(np.uint32)
        merged_at = old_at.astype(np.uint32)
    if merged_ag.shape[0] > 0:
        write_snap.angles.N = int(merged_ag.shape[0])
        write_snap.angles.types = angle_types
        write_snap.angles.group[:] = merged_ag
        write_snap.angles.typeid[:] = merged_at
    elif angle_types:
        write_snap.angles.N = 0
        write_snap.angles.types = angle_types

    # Pass-through dihedrals/impropers if any.
    for grp_name in ("dihedrals", "impropers"):
        src = getattr(read_snap, grp_name)
        dst = getattr(write_snap, grp_name)
        if int(src.N) > 0:
            dst.N = int(src.N)
            dst.types = list(src.types)
            dst.group[:] = np.asarray(src.group)
            dst.typeid[:] = np.asarray(src.typeid)

    return write_snap


# ---------------------------------------------------------------------------
# Public factory: build complete lamellipodium-only sim
# ---------------------------------------------------------------------------
def build_lamellipodium_simulation(
    p: ResolvedH5,
    *,
    device: hoomd.device.Device | None = None,
    rng: np.random.Generator | None = None,
):
    """End-to-end lamellipodium-only HOOMD sim builder (for KU-5.x test
    harness).  Returns dict of handles: sim, layout, state,
    elong_updater, branch_updater, cap_updater, wave_pin_force.
    """
    import gsd.hoomd

    layout = generate_lamellipodium_layout(p, wave_tag_start=0, rng=rng)

    # State setup.
    state = LamellipodiumState()
    # Mothers start as barbed ends, tangent = -ŷ (away from membrane).
    default_tangent = np.array([0.0, -1.0, 0.0], dtype=np.float64)
    tangents = layout.mother_tangents  # per-WAVE for non-flat geometries, else None
    n_seed = p.n_WAVE
    for i in range(n_seed):
        mother_tag = layout.mother_tag_start + i
        state.barbed_end_tags.append(mother_tag)
        t = default_tangent if tangents is None else np.asarray(tangents[i], dtype=np.float64)
        state.tangent_of[mother_tag] = t.copy()
    state.actin_next_tag = layout.mother_tag_start + n_seed

    # Build snapshot: wave particles + mother seeds + initial WAVE-mother anchor bonds.
    n_total = layout.wave_positions.shape[0] + layout.mother_seed_positions.shape[0]
    snap = gsd.hoomd.Frame()
    snap.particles.N = n_total
    snap.particles.types = ["wave_particle", "actin_lamel"]
    typeids = np.empty(n_total, dtype=np.uint32)
    typeids[:layout.wave_positions.shape[0]] = 0   # wave
    typeids[layout.wave_positions.shape[0]:] = 1   # actin_lamel
    snap.particles.typeid = typeids
    snap.particles.position = np.concatenate(
        [layout.wave_positions, layout.mother_seed_positions], axis=0
    )
    snap.particles.mass = np.ones(n_total, dtype=np.float64)

    # Bonds: WAVE-mother anchor (lamel_wave_anchor) + (no internal yet).
    snap.bonds.N = layout.wave_to_mother_bond_pairs.shape[0]
    snap.bonds.types = ["lamel_wave_anchor", "lamel_actin_bond", "lamel_branch_bond"]
    snap.bonds.typeid = np.zeros(snap.bonds.N, dtype=np.uint32)
    snap.bonds.group = layout.wave_to_mother_bond_pairs.astype(np.uint32)

    # Angles: branch angle type (no instances at construction).
    snap.angles.N = 0
    snap.angles.types = ["lamel_branch_angle"]

    snap.configuration.box = [p.L_box, p.L_box, p.L_box, 0.0, 0.0, 0.0]

    sim = hoomd.Simulation(
        device=device or hoomd.device.CPU(), seed=p.seed
    )
    sim.create_state_from_snapshot(snap)

    # Forces: WAVE-anchor + actin backbone + branch + branch-angle.
    bond = md.bond.Harmonic()
    bond.params["lamel_wave_anchor"] = dict(k=p.k_wave_pin * 0.1, r0=p.rest_length)
    bond.params["lamel_actin_bond"] = dict(k=p.bond_k, r0=p.rest_length)
    bond.params["lamel_branch_bond"] = dict(k=p.bond_k, r0=p.rest_length)

    angle = md.angle.Harmonic()
    angle.params["lamel_branch_angle"] = dict(
        k=p.angle_branch_k, t0=p.angle_branch_t0
    )

    # WAVE plane pin custom force.
    wave_pin = WaveMembranePin(
        p, wave_tag_start=layout.wave_tag_start, n_WAVE=p.n_WAVE,
    )

    ig = md.Integrator(dt=p.dt)
    ig.forces.append(bond)
    ig.forces.append(angle)
    ig.forces.append(wave_pin)
    sim.operations.integrator = ig

    # Updaters: BAOAB (D3) + three D1 batched updaters.
    from ffn_sim.archive.hoomd_legacy.integrator.baoab import make_baoab_updater
    gamma_map = {
        "actin_lamel": 6.0 * math.pi * 6.913e-4 * p.bead_radius,
        "wave_particle": 6.0 * math.pi * 6.913e-4 * p.bead_radius,
    }
    _, baoab_updater = make_baoab_updater(
        kT=p.kT, gamma=gamma_map, dt=p.dt, seed=p.seed,
    )
    sim.operations.updaters.append(baoab_updater)

    elong_action = BarbedEndElongationUpdater(p=p, lamel_state=state)
    elong_updater = hoomd.update.CustomUpdater(
        action=elong_action, trigger=hoomd.trigger.Periodic(p.batch_steps),
    )
    branch_action = ArpBranchingUpdater(
        p=p, lamel_state=state,
        wave_tag_start=layout.wave_tag_start, n_WAVE=p.n_WAVE,
    )
    branch_updater = hoomd.update.CustomUpdater(
        action=branch_action, trigger=hoomd.trigger.Periodic(p.batch_steps),
    )
    cap_action = CappingUpdater(p=p, lamel_state=state)
    cap_updater = hoomd.update.CustomUpdater(
        action=cap_action, trigger=hoomd.trigger.Periodic(p.batch_steps),
    )
    sim.operations.updaters.append(elong_updater)
    sim.operations.updaters.append(branch_updater)
    sim.operations.updaters.append(cap_updater)

    return {
        "sim": sim,
        "layout": layout,
        "state": state,
        "elong_updater": elong_updater,
        "elong_action": elong_action,
        "branch_updater": branch_updater,
        "branch_action": branch_action,
        "cap_updater": cap_updater,
        "cap_action": cap_action,
        "wave_pin_force": wave_pin,
        "n_wave": p.n_WAVE,
    }


# ---------------------------------------------------------------------------
# H.5 단계 2 — Cell composition helpers
# ---------------------------------------------------------------------------
# These two helpers expose the lamellipodium subsystem as an additive layer
# on top of an existing cortex (+ optional xlinks / myosin / ERM) HOOMD
# Simulation, mirroring the cortex.crosslinkers.extend_cortex_state_with_xlinks
# + cortex.myosin.extend_state_with_cortex_myosin pattern.
#
# The Cell.build() factory in `cell/cell.py` calls
# extend_cortex_snapshot_with_lamellipodium() to grow the cortex snapshot
# with WAVE + mother-actin particles + the WAVE-anchor bonds, then calls
# attach_lamellipodium_to_simulation() to wire bond / angle / membrane-pin
# forces + the three D2-batched Updaters into the live simulation.
#
# Both helpers assume the host simulation has NOT yet been ``sim.run(...)``
# (the BAOAB Action attaches at first run, so mutating its ``gamma_map`` is
# safe up until that point).


def extend_cortex_snapshot_with_lamellipodium(
    base_snap,
    p: ResolvedH5,
    *,
    wave_tag_start: int,
    rng: np.random.Generator | None = None,
    geometry: str = "flat_plane",
    R_cell: float | None = None,
    cap_depth: float | None = None,
    contact_radius_frac: float | None = None,
    basal_z_offset: float | None = None,
    polarization: "np.ndarray | tuple[float, float, float] | None" = None,
    half_angle_azimuth: float = math.pi / 6.0,
    half_angle_linear: float | None = None,
    fan_spread: float | None = None,
) -> tuple["LamellipodiumLayout", "LamellipodiumState", Any]:
    """Extend an existing cortex (+ optional xlinks / myosin) snapshot.

    Mirrors :func:`ffn_sim.archive.hoomd_legacy.cortex.crosslinkers.extend_cortex_state_with_xlinks`
    and :func:`ffn_sim.archive.hoomd_legacy.cortex.myosin.extend_state_with_cortex_myosin`.

    Appends ``n_WAVE`` ``wave_particle`` beads at the WAVE membrane plane,
    ``n_WAVE`` ``actin_lamel`` mother-actin seed beads one rest-length
    below each WAVE, plus ``n_WAVE`` ``lamel_wave_anchor`` bonds joining
    each WAVE to its mother seed.  The bond types ``lamel_actin_bond``
    and ``lamel_branch_bond`` and the angle type ``lamel_branch_angle``
    are registered (with zero initial instances) so the runtime Updaters'
    snapshot-rebuild calls do not have to grow the type registry.

    Parameters
    ----------
    base_snap : gsd.hoomd.Frame | hoomd.Snapshot
        Snapshot from the cortex (+ xlinks / myosin) construction step,
        prior to ``sim.create_state_from_snapshot`` OR after via
        ``sim.state.get_snapshot()`` — both share the .particles/.bonds/
        .angles attribute layout.
    p : ResolvedH5
        Resolved lamellipodium parameters (provides ``n_WAVE``,
        ``Y_max``, ``rest_length``, etc.).
    wave_tag_start : int
        Global tag offset for the first WAVE particle (= current snapshot
        ``particles.N`` — the lamellipodium block is appended at the end
        of the tag space, after cortex / xlinks / myosin).
    rng : np.random.Generator, optional
        WAVE-placement RNG. Defaults to ``np.random.default_rng(p.seed)``.

    Returns
    -------
    snap_new : gsd.hoomd.Frame
        New snapshot extending ``base_snap`` with the lamellipodium beads
        and bonds.  Particle / bond / angle TYPES are registered for the
        downstream D2-batched Updaters even when the per-type instance
        count is zero at construction.
    layout : LamellipodiumLayout
        WAVE + mother-seed layout (positions, tag ranges, anchor pairs).
    state : LamellipodiumState
        Pre-seeded runtime state: each mother is registered as a
        barbed-end with unit tangent ``-ŷ`` (grows away from the WAVE
        plane at ``+Y_max`` into cytosol).
    """
    import gsd.hoomd

    if rng is None:
        rng = np.random.default_rng(p.seed)

    layout = generate_lamellipodium_layout(
        p, wave_tag_start=wave_tag_start, rng=rng,
        geometry=geometry, R_cell=R_cell, cap_depth=cap_depth,
        contact_radius_frac=contact_radius_frac, basal_z_offset=basal_z_offset,
        polarization=polarization, half_angle_azimuth=half_angle_azimuth,
        half_angle_linear=half_angle_linear, fan_spread=fan_spread,
    )

    snap_old = base_snap
    n_part_old = int(snap_old.particles.N)
    n_WAVE = p.n_WAVE
    n_new_particles = 2 * n_WAVE   # WAVE + mother actin (one each)

    snap = gsd.hoomd.Frame()
    snap.particles.N = n_part_old + n_new_particles

    # ---- Particle types: existing + (wave_particle, actin_lamel) ----
    old_ptypes = list(snap_old.particles.types)
    new_ptypes = list(old_ptypes)
    if "wave_particle" not in new_ptypes:
        new_ptypes.append("wave_particle")
    if "actin_lamel" not in new_ptypes:
        new_ptypes.append("actin_lamel")
    wave_typeid = new_ptypes.index("wave_particle")
    actin_lamel_typeid = new_ptypes.index("actin_lamel")
    snap.particles.types = new_ptypes

    typeids = np.empty(snap.particles.N, dtype=np.uint32)
    typeids[:n_part_old] = np.asarray(snap_old.particles.typeid)
    if n_WAVE > 0:
        typeids[n_part_old:n_part_old + n_WAVE] = wave_typeid
        typeids[n_part_old + n_WAVE:] = actin_lamel_typeid
    snap.particles.typeid = typeids

    # ---- Positions ----
    pos_new = np.empty((snap.particles.N, 3), dtype=np.float64)
    pos_new[:n_part_old] = np.asarray(snap_old.particles.position)
    if n_WAVE > 0:
        pos_new[n_part_old:n_part_old + n_WAVE] = layout.wave_positions
        pos_new[n_part_old + n_WAVE:] = layout.mother_seed_positions
    snap.particles.position = pos_new

    # ---- Masses ----
    mass_new = np.empty(snap.particles.N, dtype=np.float64)
    mass_new[:n_part_old] = np.asarray(snap_old.particles.mass)
    mass_new[n_part_old:] = 1.0
    snap.particles.mass = mass_new

    # ---- Bond types: existing + lamel_wave_anchor + lamel_actin_bond +
    #      lamel_branch_bond.  Only lamel_wave_anchor has instances at
    #      construction; the two filament-growth types are pre-registered
    #      so the runtime Updaters' snapshot-rebuilds don't grow the
    #      registry mid-simulation. ----
    old_btypes = list(snap_old.bonds.types)
    new_btypes = list(old_btypes)
    for name in (
        "lamel_wave_anchor", "lamel_actin_bond", "lamel_branch_bond",
    ):
        if name not in new_btypes:
            new_btypes.append(name)
    anchor_btid = new_btypes.index("lamel_wave_anchor")

    old_bg = np.asarray(snap_old.bonds.group, dtype=np.uint32)
    old_bt = np.asarray(snap_old.bonds.typeid, dtype=np.uint32)
    if n_WAVE > 0:
        anchor_bg = layout.wave_to_mother_bond_pairs.astype(np.uint32)
        anchor_bt = np.full(n_WAVE, anchor_btid, dtype=np.uint32)
        merged_bg = np.concatenate([old_bg, anchor_bg], axis=0)
        merged_bt = np.concatenate([old_bt, anchor_bt])
    else:
        merged_bg = old_bg
        merged_bt = old_bt
    snap.bonds.N = int(merged_bg.shape[0])
    snap.bonds.types = new_btypes
    snap.bonds.group = merged_bg
    snap.bonds.typeid = merged_bt

    # ---- Angle types: existing + lamel_branch_angle (no instances yet). ----
    old_atypes = list(snap_old.angles.types) if int(snap_old.angles.N) > 0 \
        else list(snap_old.angles.types)
    new_atypes = list(old_atypes)
    if "lamel_branch_angle" not in new_atypes:
        new_atypes.append("lamel_branch_angle")
    snap.angles.N = int(snap_old.angles.N)
    snap.angles.types = new_atypes
    if int(snap_old.angles.N) > 0:
        snap.angles.group = np.asarray(snap_old.angles.group)
        snap.angles.typeid = np.asarray(snap_old.angles.typeid)

    # ---- Pass-through dihedrals / impropers if present. ----
    for grp_name in ("dihedrals", "impropers"):
        src = getattr(snap_old, grp_name)
        dst = getattr(snap, grp_name)
        if int(src.N) > 0:
            dst.N = int(src.N)
            dst.types = list(src.types)
            dst.group = np.asarray(src.group)
            dst.typeid = np.asarray(src.typeid)

    snap.configuration.box = list(snap_old.configuration.box)

    # ---- Initialise the runtime state: each mother is a barbed end with
    #      tangent -ŷ (away from membrane). ----
    state = LamellipodiumState()
    default_tangent = np.array([0.0, -1.0, 0.0], dtype=np.float64)
    tangents = layout.mother_tangents  # per-WAVE for non-flat geometries, else None
    for i in range(n_WAVE):
        mother_tag = layout.mother_tag_start + i
        state.barbed_end_tags.append(mother_tag)
        t = default_tangent if tangents is None else np.asarray(tangents[i], dtype=np.float64)
        state.tangent_of[mother_tag] = t.copy()
    state.actin_next_tag = layout.mother_tag_start + n_WAVE

    return snap, layout, state


def attach_lamellipodium_to_simulation(
    sim: hoomd.Simulation,
    p: ResolvedH5,
    *,
    layout: "LamellipodiumLayout",
    state: "LamellipodiumState",
    baoab_action: Any | None = None,
) -> dict:
    """Wire lamellipodium bond / angle / membrane-pin forces + 3 Updaters.

    Pre-conditions
    --------------
    * ``sim.state`` already contains the WAVE + mother-actin particles
      AND the ``wave_particle`` / ``actin_lamel`` particle types AND the
      ``lamel_wave_anchor`` / ``lamel_actin_bond`` / ``lamel_branch_bond``
      bond types AND the ``lamel_branch_angle`` angle type.  Caller must
      have invoked :func:`extend_cortex_snapshot_with_lamellipodium` and
      built the sim from that extended snapshot via
      ``sim.create_state_from_snapshot`` (note: HOOMD's ``set_snapshot``
      CANNOT add new particle / bond / angle types post-init; only
      ``create_state_from_snapshot`` registers them).
    * ``sim.operations.integrator`` exists with a ``md.bond.Harmonic``
      force, a ``md.angle.Harmonic`` force, and ``md.Integrator.dt == p.dt``.
    * If ``baoab_action`` is provided, it has NOT yet been
      ``sim.run(...)``-attached (its ``gamma_map`` is still mutable).

    Side-effects on ``sim``:

    1. Registers ``lamel_wave_anchor`` / ``lamel_actin_bond`` /
       ``lamel_branch_bond`` params on the existing ``bond.Harmonic`` force.
    2. Registers ``lamel_branch_angle`` params on the existing
       ``angle.Harmonic`` force.
    3. Appends a fresh :class:`WaveMembranePin` custom force to the
       integrator.
    4. If ``baoab_action`` is passed, extends its ``gamma_map`` with the
       new particle types so the BAOAB attach-time validation passes.
    5. Appends the three D2-batched lamellipodium Updaters
       (``BarbedEndElongationUpdater``, ``ArpBranchingUpdater``,
       ``CappingUpdater``) wrapped in ``hoomd.update.CustomUpdater`` with
       ``hoomd.trigger.Periodic(p.batch_steps)``.

    Returns
    -------
    dict
        Handles: ``layout``, ``state``, ``wave_pin_force``,
        ``elong_action`` / ``elong_updater``, ``branch_action`` /
        ``branch_updater``, ``cap_action`` / ``cap_updater``,
        ``n_wave_particles``, ``n_lamellipodium_actin``.
    """
    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "attach_lamellipodium_to_simulation requires "
            "sim.operations.integrator to be wired."
        )
    if not math.isclose(float(ig.dt), p.dt, rel_tol=0.0, abs_tol=0.0):
        raise ValueError(
            f"Lamellipodium dt={p.dt!r} must equal host integrator "
            f"dt={float(ig.dt)!r}; they share the BAOAB step."
        )

    # 1. Locate existing bond / angle forces and register the new types.
    bond_force: md.bond.Harmonic | None = None
    angle_force: md.angle.Harmonic | None = None
    for f in ig.forces:
        if isinstance(f, md.bond.Harmonic) and bond_force is None:
            bond_force = f
        elif isinstance(f, md.angle.Harmonic) and angle_force is None:
            angle_force = f
    if bond_force is None or angle_force is None:
        raise RuntimeError(
            "attach_lamellipodium_to_simulation requires the host "
            "integrator to expose both md.bond.Harmonic and md.angle.Harmonic "
            "force objects (cortex builders attach both)."
        )

    # WAVE-anchor harmonic is intentionally SOFT (k = 0.1 · k_wave_pin) so
    # it never tightens the integration CFL relative to the WAVE plane
    # pin itself (which is gated against dt at resolve time).
    bond_force.params["lamel_wave_anchor"] = dict(
        k=p.k_wave_pin * 0.1, r0=p.rest_length,
    )
    bond_force.params["lamel_actin_bond"] = dict(k=p.bond_k, r0=p.rest_length)
    bond_force.params["lamel_branch_bond"] = dict(k=p.bond_k, r0=p.rest_length)
    angle_force.params["lamel_branch_angle"] = dict(
        k=p.angle_branch_k, t0=p.angle_branch_t0,
    )

    # 2. WAVE plane membrane pin force compute.
    wave_pin = WaveMembranePin(
        p, wave_tag_start=layout.wave_tag_start, n_WAVE=p.n_WAVE,
    )
    ig.forces.append(wave_pin)

    # 3. Extend the BAOAB gamma_map with the new particle types so attach()
    # doesn't trip the missing-types guard.  WAVE + actin_lamel inherit
    # the cortex bead drag via the standard Stokes formula 6πη·R.
    if baoab_action is not None:
        gamma_lamel = 6.0 * math.pi * 6.913e-4 * p.bead_radius
        baoab_action.gamma_map.setdefault("wave_particle", gamma_lamel)
        baoab_action.gamma_map.setdefault("actin_lamel", gamma_lamel)

    # 4. Three D2-batched lamellipodium Updaters.
    elong_action = BarbedEndElongationUpdater(p=p, lamel_state=state)
    elong_updater = hoomd.update.CustomUpdater(
        action=elong_action, trigger=hoomd.trigger.Periodic(p.batch_steps),
    )
    branch_action = ArpBranchingUpdater(
        p=p, lamel_state=state,
        wave_tag_start=layout.wave_tag_start, n_WAVE=p.n_WAVE,
    )
    branch_updater = hoomd.update.CustomUpdater(
        action=branch_action, trigger=hoomd.trigger.Periodic(p.batch_steps),
    )
    cap_action = CappingUpdater(p=p, lamel_state=state)
    cap_updater = hoomd.update.CustomUpdater(
        action=cap_action, trigger=hoomd.trigger.Periodic(p.batch_steps),
    )
    sim.operations.updaters.append(elong_updater)
    sim.operations.updaters.append(branch_updater)
    sim.operations.updaters.append(cap_updater)

    return {
        "layout": layout,
        "state": state,
        "wave_pin_force": wave_pin,
        "elong_action": elong_action,
        "elong_updater": elong_updater,
        "branch_action": branch_action,
        "branch_updater": branch_updater,
        "cap_action": cap_action,
        "cap_updater": cap_updater,
        "n_wave_particles": int(p.n_WAVE),
        "n_lamellipodium_actin": int(p.n_WAVE),  # 1 mother per WAVE at construction
    }
