"""H.3 cortical myosin minifilament topology + runtime stepping (D5 + D6).

Phase 1 H.3 brief §Myosin minifilaments. Stam-Hocky bipolar
architecture per PHASE_0_3_DECISIONS D5 + Hill force-velocity per D6.

This module is the CORTEX-SHELL placement of the D5 minifilament; the
underlying Stam-Hocky construction kernel + D6 Hill closed-form live in
``ffn_sim/bridge/motor.py`` (H.4 Sub-owned, shared module per H.4
brief). cortex/myosin.py ADDS:

* Cortical placement: ``n_motors_per_cell`` minifilaments randomly
  positioned on the cortex shell ``r = R_cell``, oriented tangentially
  (random in-plane direction).
* State extension: appends minifilament particles + intra-minifilament
  bonds onto an existing cortex (+xlink, +ERM, ...) state frame.
* Optional runtime stepper ``MyosinStepUpdater`` (D2 Bell-Evans slip
  head-to-actin bonds + D6 Hill stepping kernel; batched
  ``hoomd.custom.Action``).

Architecture per minifilament (brief literal):

* **Backbone**: ``n_backbone = 14`` beads, total length ``700 nm``,
  approximated as a stiff-harmonic rigid rod (``k_backbone =
  10 · k_head_spring`` matches H.4 motor.py convention). Particle type
  ``cortex_myosin_backbone``.
* **Heads**: ``n_heads_per_side = 10`` cross-bridge heads on EACH side
  of the backbone, perpendicular offset ``r0 = 200 nm`` via harmonic
  spring ``k_head_spring = 1 pN/μm = 1e-6 N/m``. Particle type
  ``cortex_myosin_head``.
* **Totals**: 14 + 2·10 = ``34 particles per minifilament``;
  ``34 · 100 = 3 400`` motor beads per cell (per brief §Myosin).

Per-head runtime dynamics (when ``MyosinStepUpdater`` attached):

* **D2 Bell-Evans slip** head ↔ actin attach bonds (mirrors
  ``ffn_sim/cortex/crosslinkers.py`` Bell-Evans Updater; reuses the
  same ``xlink_attach_b{i}`` per-r0 binning convention applied to
  the dynamic head-to-actin pair).
* **D6 Hill stepping** while engaged: head's stepping velocity
  ``v_step(F) = v0 · (F_s − F) / (F_s + F/a_over_F_stall)`` (Hill 1938
  with Kovács 2003 ``a/F_s = 0.5`` for NMII per brief literal).
  Implementation: advance the ``head_attach_bond.r0`` per tick by
  ``v_step · batch_dt`` so the rest-length of the dynamic attach bond
  shrinks → equivalent to the head ratcheting along actin toward the
  filament minus end. Hill velocity computation uses
  ``ffn_sim.bridge.motor.hill_velocity_clamped``.

CFL note
--------
Brief literal ``k_head_spring = 1 pN/μm = 1e-6 N/m``. With per-bead
``γ_b ≈ 3.91·10⁻¹⁰ N·s/m``, ``τ_head = γ_b/k_head_spring ≈ 0.4 ms``.
Backbone ``k_backbone = 10·k_head_spring = 1e-5 N/m`` gives
``τ_backbone ≈ 39 μs``. Both timescales are FAR ABOVE the cortex
``dt_CFL = 13 ns``, so attaching myosin does NOT tighten the
integration step (unlike ERM where ``k_ERM = 0.1 N/m`` violates CFL).
No CFL gate raised for myosin.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC checks in
``ffn_sim/tests/test_myosin.py``.*

1. **Dimensional analysis**
   - ``k_head_spring`` [N/m] · ``head_rest_length`` [m] → force [N]. ✓
   - ``k_backbone`` [N/m] · ``backbone_segment_length`` [m] → force [N].
   - Hill: ``v0`` [m/s], ``F_stall`` [N], ``F`` [N] → ``v_step`` [m/s].
     ``v_step · batch_dt`` [m] = bond r0 advance per tick.
   - Bell-Evans slip: ``k_off(F) = k_off⁰ · exp(F · x_β / kT)`` [1/s].

2. **Boundary cases**
   - ``n_motors_per_cell = 0``: empty layout; extend_state is no-op.
   - ``n_backbone < 2``: D5 forbidden (matches H.4 motor.py guard).
   - ``n_heads_per_side ≤ 0``: D5 forbidden.
   - ``backbone_length ≤ 0`` or ``head_rest_length ≤ 0``: ValueError.

3. **Conservation invariants**
   - Particle count: extend adds exactly
     ``n_motors_per_cell · (n_backbone + 2·n_heads_per_side)``
     particles.
   - Bond count: extend adds
     ``n_motors_per_cell · ((n_backbone-1) + 2·n_heads_per_side)``
     static bonds (backbone-backbone + head-backbone).
   - Dynamic head-actin bonds are added/removed by MyosinStepUpdater at
     runtime; bounded by ``2 · n_motors_per_cell · n_heads_per_side``
     per cell.
   - Newton 3rd law inherited from md.bond.Harmonic.

4. **Numerical sanity**
   - All positions / forces float64.
   - dt-CFL for myosin: ``τ_head, τ_backbone ≫ dt_cortex`` (see CFL
     note above). No additional gate required.
   - Hill velocity clamped to zero on super-stall lengthening branch
     (``hill_velocity_clamped`` from H.4 motor.py).

5. **Sign / sense**
   - Hill stepping is forward (positive ``v_step``) when ``F < F_stall``,
     zero in the super-stall regime. STATIC test asserts.
   - Bell-Evans slip head-actin: ``k_off`` monotonically INCREASES with
     ``|F|`` (slip), same as ``cortex/crosslinkers.py`` Bell-Evans gate.

6. **Measurement protocol**
   - Cortex tension γ_cortex emerges from collective myosin head action
     pulling on the cortex actin shell. KU-3.5 brief gate: γ_cortex ≈
     0.5 mN/m ± 30 %. Implemented in
     ``tests/validation/test_ku35_tension.py`` (separate test file,
     wired to Cell.build with myosin enabled).

References
----------
- Brief: ``ffn_sim/docs/briefs/H3_cortex.md`` §Myosin minifilaments.
- KU-3.5 (γ_cortex = 0.5 mN/m), KU-3.18 (blebbistatin).
- ``ffn_sim/bridge/motor.py`` (H.4 shared D5 + D6 module; read-only
  import here per CLAUDE.md cross-module ownership table).
- AFINES algorithm notes §4.3 (myosin stepping kernel).
- Stam-Hocky bipolar minifilament: Stam et al. 2017 PNAS.
- Hill 1938; Kovács 2003 (NMII tuning, a/F_s = 0.5).
- Salbreux 2012 (cortex myosin density 3/μm²).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md

# Read-only import from Sub-owned bridge/ (CLAUDE.md cross-module read).
from ffn_sim.bridge.motor import hill_velocity_clamped


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedCortexMyosin:
    """Cortical myosin (Stam-Hocky bipolar minifilament + D6 Hill) params."""

    n_motors_per_cell: int       # KU-3.x 100 per cell (Salbreux 2012)
    n_backbone: int              # 14 (brief literal)
    n_heads_per_side: int        # 10 (brief literal)
    backbone_length: float       # m   700 nm (brief literal)
    head_rest_length: float      # m   200 nm perpendicular head offset

    # Spring constants
    k_head_spring: float         # N/m head-backbone spring k (1 pN/μm = 1e-6 N/m)
    k_backbone_factor: float     # backbone k = factor · k_head_spring (10 per H.4)
    k_head_actin: float          # N/m head-actin attach bond k (=k_head_spring)

    # D6 Hill
    v0_per_head: float           # m/s 1 μm/s per brief (Kovács 2003)
    F_stall_per_head: float      # N   0.5 pN per brief
    a_over_F_stall: float        # 0.5 per brief

    # D2 Bell-Evans (slip) head-actin
    head_actin_k_off0: float     # 1/s unloaded off-rate
    head_actin_x_beta: float     # m   Bell strength length
    head_actin_k_on: float       # 1/s per-head binding rate when acceptor present
    head_actin_max_bind_dist: float   # m   binding-acceptor search radius

    # D2 batch
    batch_steps: int             # BAOAB steps between MyosinStepUpdater ticks
    dt: float                    # host sim dt [s]
    n_bins: int                  # per-r0 binning count for dynamic attach bonds

    # Seed
    seed: int

    # Derived
    backbone_segment_length: float = 0.0       # backbone_length / (n_backbone-1)
    k_backbone: float = 0.0                    # k_head_spring · k_backbone_factor
    batch_dt: float = 0.0                      # batch_steps · dt
    k_off_max_at_zero_load: float = 0.0
    n_particles_per_motor: int = 0             # n_backbone + 2·n_heads_per_side
    n_static_bonds_per_motor: int = 0          # (n_backbone-1) + 2·n_heads_per_side
    extras: dict[str, Any] = field(default_factory=dict)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def resolve_cortex_myosin(cfg: dict, *, dt: float) -> ResolvedCortexMyosin:
    """Resolve cortical-myosin config block with D5 boundary + CFL gates."""
    if "cortex" in cfg:
        cfg = cfg["cortex"]
    if "myosin" in cfg:
        cfg = cfg["myosin"]

    p = ResolvedCortexMyosin(
        n_motors_per_cell=int(cfg["n_motors_per_cell"]),
        n_backbone=int(cfg["n_backbone"]),
        n_heads_per_side=int(cfg["n_heads_per_side"]),
        backbone_length=float(cfg["backbone_length"]),
        head_rest_length=float(cfg["head_rest_length"]),
        k_head_spring=float(cfg["k_head_spring"]),
        k_backbone_factor=float(cfg.get("k_backbone_factor", 10.0)),
        k_head_actin=float(cfg.get("k_head_actin", cfg["k_head_spring"])),
        v0_per_head=float(cfg["v0_per_head"]),
        F_stall_per_head=float(cfg["F_stall_per_head"]),
        a_over_F_stall=float(cfg.get("a_over_F_stall", 0.5)),
        head_actin_k_off0=float(cfg["head_actin_k_off0"]),
        head_actin_x_beta=float(cfg["head_actin_x_beta"]),
        head_actin_k_on=float(cfg["head_actin_k_on"]),
        head_actin_max_bind_dist=float(cfg["head_actin_max_bind_dist"]),
        batch_steps=int(cfg["batch_steps"]),
        dt=float(dt),
        n_bins=int(cfg.get("n_bins", 10)),
        seed=int(cfg.get("seed", 44)),
    )

    # §2 boundary checks
    if p.n_motors_per_cell < 0:
        raise ValueError(f"n_motors_per_cell must be ≥ 0; got {p.n_motors_per_cell}")
    if p.n_backbone < 2:
        raise ValueError(
            f"n_backbone must be ≥ 2 (D5 brief literal = 14); got {p.n_backbone}"
        )
    if p.n_heads_per_side <= 0:
        raise ValueError(
            f"n_heads_per_side must be ≥ 1 (D5 brief literal = 10); "
            f"got {p.n_heads_per_side}"
        )
    for name, x in [
        ("backbone_length", p.backbone_length),
        ("head_rest_length", p.head_rest_length),
        ("k_head_spring", p.k_head_spring),
        ("v0_per_head", p.v0_per_head),
        ("F_stall_per_head", p.F_stall_per_head),
        ("head_actin_k_off0", p.head_actin_k_off0),
        ("head_actin_x_beta", p.head_actin_x_beta),
        ("head_actin_k_on", p.head_actin_k_on),
        ("head_actin_max_bind_dist", p.head_actin_max_bind_dist),
        ("dt", p.dt),
    ]:
        _require_finite_positive(name, x)
    if p.batch_steps < 1:
        raise ValueError(f"batch_steps must be ≥ 1; got {p.batch_steps}")
    if p.n_bins < 1:
        raise ValueError(f"n_bins must be ≥ 1; got {p.n_bins}")

    # §1 derived
    p.backbone_segment_length = p.backbone_length / (p.n_backbone - 1)
    p.k_backbone = p.k_head_spring * p.k_backbone_factor
    p.batch_dt = p.batch_steps * p.dt
    p.k_off_max_at_zero_load = p.head_actin_k_off0  # F = 0 floor
    p.n_particles_per_motor = p.n_backbone + 2 * p.n_heads_per_side
    p.n_static_bonds_per_motor = (p.n_backbone - 1) + 2 * p.n_heads_per_side

    # §2 D2 batch CFL: batch_dt · k_off_max ≤ 1e-3 (mirrors crosslinkers.py)
    cfl_product = p.batch_dt * p.k_off_max_at_zero_load
    if cfl_product > 1.0e-3:
        target_batch_dt = 1.0e-3 / p.k_off_max_at_zero_load
        new_batch_steps = max(1, int(math.floor(target_batch_dt / p.dt)))
        if new_batch_steps < 1:
            raise RuntimeError(
                f"D2 batch CFL impossible: dt={p.dt:.3e} > 1e-3/k_off_max."
            )
        p.batch_steps = new_batch_steps
        p.batch_dt = p.batch_steps * p.dt
        p.extras["batch_steps_shrunk_from"] = int(cfg["batch_steps"])

    return p


# ---------------------------------------------------------------------------
# Bond type names + per-r0 binning
# ---------------------------------------------------------------------------
BOND_TYPE_MYOSIN_BACKBONE = "cortex_myosin_backbone"
BOND_TYPE_MYOSIN_HEAD_BACKBONE = "cortex_myosin_head_backbone"


def cortex_myosin_attach_bin_names(n_bins: int) -> list[str]:
    """Dynamic head-actin attach bond type names (per-r0 binned)."""
    return [f"cortex_myosin_attach_b{i}" for i in range(n_bins)]


def cortex_myosin_attach_bin_rest_lengths(
    n_bins: int, max_bind_dist: float
) -> np.ndarray:
    edges = np.linspace(0.0, max_bind_dist, n_bins + 1, dtype=np.float64)
    return 0.5 * (edges[:-1] + edges[1:])


# ---------------------------------------------------------------------------
# Cortex shell topology generator
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class CortexMyosinLayout:
    """Per-minifilament tag bookkeeping in the cortex-myosin frame.

    Tag conventions (relative to the start of the myosin block):
    - Tag offset 0..n_backbone-1: backbone beads (in axis order).
    - Tag offset n_backbone..n_backbone+n_heads-1: + polarity heads.
    - Tag offset n_backbone+n_heads..n_backbone+2·n_heads-1: − polarity heads.

    Attributes
    ----------
    centers : (n_motors, 3) float64 — minifilament center positions on shell.
    axes : (n_motors, 3) float64 — random tangent-plane unit axis per motor.
    positions : (n_motors, n_part_per_motor, 3) float64 — bead positions.
    motor_tag_start : int — global tag index where the first myosin particle sits.
    """

    centers: np.ndarray
    axes: np.ndarray
    positions: np.ndarray
    motor_tag_start: int


def _sample_sphere_surface(rng: np.random.Generator, n: int, R: float) -> np.ndarray:
    """Marsaglia 1972 — same routine as cortex.py."""
    pts = np.empty((n, 3), dtype=np.float64)
    drawn = 0
    while drawn < n:
        batch = max(2 * (n - drawn), 64)
        u = rng.uniform(-1.0, 1.0, batch)
        v = rng.uniform(-1.0, 1.0, batch)
        s = u * u + v * v
        mask = s < 1.0
        u, v, s = u[mask], v[mask], s[mask]
        k = min(len(u), n - drawn)
        if k == 0:
            continue
        u, v, s = u[:k], v[:k], s[:k]
        factor = 2.0 * np.sqrt(1.0 - s)
        pts[drawn:drawn + k, 0] = R * u * factor
        pts[drawn:drawn + k, 1] = R * v * factor
        pts[drawn:drawn + k, 2] = R * (1.0 - 2.0 * s)
        drawn += k
    return pts


def _tangent_plane_basis(normals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = normals.shape[0]
    abs_n = np.abs(normals)
    min_axis = np.argmin(abs_n, axis=1)
    ref = np.zeros_like(normals)
    ref[np.arange(n), min_axis] = 1.0
    e1 = np.cross(ref, normals)
    e1 = e1 / np.linalg.norm(e1, axis=1, keepdims=True).clip(min=1e-30)
    e2 = np.cross(normals, e1)
    return e1, e2


def generate_cortex_myosin_layout(
    p: ResolvedCortexMyosin, R_cell: float,
    *, motor_tag_start: int, rng: np.random.Generator | None = None,
) -> CortexMyosinLayout:
    """Place ``n_motors_per_cell`` minifilaments on the cortex shell.

    Each minifilament:

    1. Center sampled uniformly on sphere ``r = R_cell`` (Marsaglia).
    2. Random tangent-plane axis ``u`` (uniform azimuth).
    3. Backbone beads laid along ``u`` at ``backbone_segment_length``
       spacing, centered on the minifilament center.
    4. Heads placed perpendicular to ``u`` along the local outward
       normal ``n``: + polarity heads at ``+head_rest_length · n``,
       − polarity heads at ``−head_rest_length · n``, distributed evenly
       along the backbone span.
    """
    if rng is None:
        rng = np.random.default_rng(p.seed)

    M = p.n_motors_per_cell
    N = p.n_backbone
    H = p.n_heads_per_side
    L = p.backbone_length
    seg = p.backbone_segment_length
    head_off = p.head_rest_length

    if M == 0:
        return CortexMyosinLayout(
            centers=np.empty((0, 3), dtype=np.float64),
            axes=np.empty((0, 3), dtype=np.float64),
            positions=np.empty((0, N + 2 * H, 3), dtype=np.float64),
            motor_tag_start=motor_tag_start,
        )

    centers = _sample_sphere_surface(rng, M, R_cell)
    normals = centers / R_cell
    e1, e2 = _tangent_plane_basis(normals)
    phi = rng.uniform(0.0, 2.0 * math.pi, M)
    axes = (
        np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2
    )
    axes = axes / np.linalg.norm(axes, axis=1, keepdims=True).clip(min=1e-30)

    positions = np.empty((M, N + 2 * H, 3), dtype=np.float64)
    # Backbone bead offsets along axis: i − (N−1)/2, scaled by seg.
    backbone_offsets = (np.arange(N, dtype=np.float64) - 0.5 * (N - 1)) * seg
    # Per-head offset along axis (heads distributed evenly along backbone span).
    if H > 1:
        head_axis_offsets = (np.arange(H, dtype=np.float64) - 0.5 * (H - 1)) * (L / max(H - 1, 1))
    else:
        head_axis_offsets = np.zeros(H, dtype=np.float64)

    for m in range(M):
        cm = centers[m]
        u = axes[m]
        n = normals[m]
        # Backbone
        for i in range(N):
            positions[m, i] = cm + backbone_offsets[i] * u
        # + heads (outward direction)
        for i in range(H):
            positions[m, N + i] = cm + head_axis_offsets[i] * u + head_off * n
        # − heads (inward direction)
        for i in range(H):
            positions[m, N + H + i] = cm + head_axis_offsets[i] * u - head_off * n

    return CortexMyosinLayout(
        centers=centers,
        axes=axes,
        positions=positions,
        motor_tag_start=motor_tag_start,
    )


# ---------------------------------------------------------------------------
# Frame builder: extend an existing cortex (+ xlinks + erm) frame
# ---------------------------------------------------------------------------
def extend_state_with_cortex_myosin(
    base_snap,
    layout: CortexMyosinLayout,
    p_myo: ResolvedCortexMyosin,
):
    """Append cortex_myosin particles + intra-minifilament bonds.

    Adds two new particle types (``cortex_myosin_backbone`` +
    ``cortex_myosin_head``) and three new bond types:
    ``cortex_myosin_backbone``, ``cortex_myosin_head_backbone``, and
    ``cortex_myosin_attach_b{i}`` (i=0..n_bins-1, for the dynamic head-
    to-actin bonds; no bonds of this type are added at construction).
    """
    import gsd.hoomd

    snap_old = base_snap
    n_part_old = int(snap_old.particles.N)
    M = p_myo.n_motors_per_cell
    n_part_per_motor = p_myo.n_particles_per_motor
    n_new_particles = M * n_part_per_motor

    snap = gsd.hoomd.Frame()
    snap.particles.N = n_part_old + n_new_particles

    # Particle types: existing + (backbone, head)
    old_types = list(snap_old.particles.types)
    backbone_type = "cortex_myosin_backbone"
    head_type = "cortex_myosin_head"
    new_particle_types = list(old_types)
    if backbone_type not in new_particle_types:
        new_particle_types.append(backbone_type)
    if head_type not in new_particle_types:
        new_particle_types.append(head_type)
    backbone_typeid = new_particle_types.index(backbone_type)
    head_typeid = new_particle_types.index(head_type)
    snap.particles.types = new_particle_types

    # typeids
    typeids = np.empty(snap.particles.N, dtype=np.uint32)
    typeids[:n_part_old] = np.asarray(snap_old.particles.typeid)
    if M > 0:
        N = p_myo.n_backbone
        H = p_myo.n_heads_per_side
        per_motor_tids = np.empty(n_part_per_motor, dtype=np.uint32)
        per_motor_tids[:N] = backbone_typeid
        per_motor_tids[N:] = head_typeid
        typeids[n_part_old:] = np.tile(per_motor_tids, M)
    snap.particles.typeid = typeids

    # positions
    pos_new = np.empty((snap.particles.N, 3), dtype=np.float64)
    pos_new[:n_part_old] = np.asarray(snap_old.particles.position)
    if M > 0:
        pos_new[n_part_old:] = layout.positions.reshape(n_new_particles, 3)
    snap.particles.position = pos_new

    # masses
    mass_new = np.empty(snap.particles.N, dtype=np.float64)
    mass_new[:n_part_old] = np.asarray(snap_old.particles.mass)
    mass_new[n_part_old:] = 1.0
    snap.particles.mass = mass_new

    # Bond types: existing + 3 new
    old_bond_types = list(snap_old.bonds.types)
    new_bond_types = list(old_bond_types)
    for new_name in (BOND_TYPE_MYOSIN_BACKBONE, BOND_TYPE_MYOSIN_HEAD_BACKBONE):
        if new_name not in new_bond_types:
            new_bond_types.append(new_name)
    backbone_bond_tid = new_bond_types.index(BOND_TYPE_MYOSIN_BACKBONE)
    head_backbone_bond_tid = new_bond_types.index(BOND_TYPE_MYOSIN_HEAD_BACKBONE)
    for name in cortex_myosin_attach_bin_names(p_myo.n_bins):
        if name not in new_bond_types:
            new_bond_types.append(name)

    # Build minifilament bonds
    old_bg = np.asarray(snap_old.bonds.group, dtype=np.uint32)
    old_bt = np.asarray(snap_old.bonds.typeid, dtype=np.uint32)

    if M > 0:
        N = p_myo.n_backbone
        H = p_myo.n_heads_per_side
        per_motor_static_bonds = p_myo.n_static_bonds_per_motor
        new_bonds = np.empty((M * per_motor_static_bonds, 2), dtype=np.uint32)
        new_btids = np.empty(M * per_motor_static_bonds, dtype=np.uint32)
        ptr = 0
        for m in range(M):
            base = n_part_old + m * n_part_per_motor
            # Backbone-backbone (N-1 bonds)
            for i in range(N - 1):
                new_bonds[ptr] = (base + i, base + i + 1)
                new_btids[ptr] = backbone_bond_tid
                ptr += 1
            # Head-backbone (+ heads): head i attaches to backbone bead i'
            # where i' is the nearest backbone bead along the axis. Since
            # head_axis_offsets are uniform along the backbone span, the
            # mapping is i_head → i_backbone = round(i_head · (N-1)/(H-1))
            # when H>1; trivially i_backbone=N//2 when H=1.
            for i in range(H):
                if H > 1:
                    ib = int(round(i * (N - 1) / max(H - 1, 1)))
                else:
                    ib = N // 2
                head_tag = base + N + i
                bb_tag = base + ib
                new_bonds[ptr] = (head_tag, bb_tag)
                new_btids[ptr] = head_backbone_bond_tid
                ptr += 1
            # Head-backbone (− heads)
            for i in range(H):
                if H > 1:
                    ib = int(round(i * (N - 1) / max(H - 1, 1)))
                else:
                    ib = N // 2
                head_tag = base + N + H + i
                bb_tag = base + ib
                new_bonds[ptr] = (head_tag, bb_tag)
                new_btids[ptr] = head_backbone_bond_tid
                ptr += 1
        snap.bonds.N = int(old_bg.shape[0] + new_bonds.shape[0])
        snap.bonds.types = new_bond_types
        snap.bonds.group = np.concatenate([old_bg, new_bonds], axis=0)
        snap.bonds.typeid = np.concatenate([old_bt, new_btids], axis=0)
    else:
        snap.bonds.N = int(old_bg.shape[0])
        snap.bonds.types = new_bond_types
        snap.bonds.group = old_bg
        snap.bonds.typeid = old_bt

    # Angles + others: pass through unchanged.
    if int(snap_old.angles.N) > 0:
        snap.angles.N = int(snap_old.angles.N)
        snap.angles.types = list(snap_old.angles.types)
        snap.angles.typeid = np.asarray(snap_old.angles.typeid)
        snap.angles.group = np.asarray(snap_old.angles.group)

    snap.configuration.box = list(snap_old.configuration.box)
    return snap


def register_cortex_myosin_bond_params(
    bond: md.bond.Harmonic, p_myo: ResolvedCortexMyosin
) -> None:
    """Wire bond.Harmonic params for the three cortex-myosin bond types."""
    bond.params[BOND_TYPE_MYOSIN_BACKBONE] = dict(
        k=p_myo.k_backbone, r0=p_myo.backbone_segment_length
    )
    bond.params[BOND_TYPE_MYOSIN_HEAD_BACKBONE] = dict(
        k=p_myo.k_head_spring, r0=p_myo.head_rest_length
    )
    bin_r0 = cortex_myosin_attach_bin_rest_lengths(
        p_myo.n_bins, p_myo.head_actin_max_bind_dist
    )
    for name, r0 in zip(cortex_myosin_attach_bin_names(p_myo.n_bins), bin_r0):
        bond.params[name] = dict(k=p_myo.k_head_actin, r0=float(r0))


# ---------------------------------------------------------------------------
# Runtime D2 Bell-Evans + D6 Hill stepping Updater
# ---------------------------------------------------------------------------
def _bell_evans_k_off(
    F: np.ndarray, k_off0: float, x_beta: float, kT: float,
) -> np.ndarray:
    return k_off0 * np.exp(F * x_beta / kT)


class MyosinStepUpdater(hoomd.custom.Action):
    """D2 Bell-Evans head-actin bond kinetics + D6 Hill stepping.

    Per batch tick:

    1. Read snapshot. Partition bonds into engaged attach bonds vs others.
    2. **Bell-Evans unbinding**: compute F for each engaged attach bond
       (using the bond's bin r0), sample break probability, remove if
       fires.
    3. **Binding**: for each unbound head, KDTree-query actin_cortex
       beads within ``head_actin_max_bind_dist``; bind to nearest with
       prob ``1 − exp(−k_on · batch_dt)``.
    4. **D6 Hill stepping**: for each surviving (post step 2-3) engaged
       attach bond, compute Hill velocity v(F) and SHRINK the bond's
       bin index by ``v · batch_dt / bin_width`` (rounded to nearest
       integer; clamped to bin 0). This advances the head along actin
       toward the minus-end (rest length goes to zero).

    Mirrors ``crosslinkers.XlinkBondUpdater`` structure exactly except
    for the additional Hill stepping kernel.
    """

    def __init__(
        self,
        *,
        p_myo: ResolvedCortexMyosin,
        layout: CortexMyosinLayout,
        kT: float,
        n_cortex_actin: int,
        seed_offset: int = 3,
    ) -> None:
        super().__init__()
        self.p = p_myo
        self.layout = layout
        self.kT = float(kT)
        self.n_cortex_actin = int(n_cortex_actin)
        self._rng = np.random.default_rng(p_myo.seed + seed_offset)

        # Per-head state — bound to which actin tag (−1 if free).
        n_heads_total = 2 * p_myo.n_heads_per_side * p_myo.n_motors_per_cell
        self._head_bound_to_actin = np.full(n_heads_total, -1, dtype=np.int64)
        # head local index → (motor_idx, head_offset_in_motor)
        self._sim_ref: hoomd.Simulation | None = None
        self._steps_run = 0
        self._n_break_total = 0
        self._n_bind_total = 0
        self._n_step_advances_total = 0

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    def _head_global_tag(self, head_local_idx: int) -> int:
        """Map head_local_idx ∈ [0, 2·H·M) → global particle tag."""
        H = self.p.n_heads_per_side
        N = self.p.n_backbone
        per_motor = self.p.n_particles_per_motor
        motor_idx = head_local_idx // (2 * H)
        within = head_local_idx % (2 * H)
        # within ∈ [0, 2H); first H are + heads, last H are − heads.
        # Layout per motor: [backbone (N), +heads (H), −heads (H)]
        return self.layout.motor_tag_start + motor_idx * per_motor + N + within

    def _head_local_from_tag(self, tag: int) -> int:
        """Inverse of _head_global_tag, returns −1 if tag is not a head."""
        N = self.p.n_backbone
        H = self.p.n_heads_per_side
        per_motor = self.p.n_particles_per_motor
        offset = tag - self.layout.motor_tag_start
        if offset < 0:
            return -1
        motor_idx = offset // per_motor
        within_motor = offset % per_motor
        if within_motor < N:
            return -1   # backbone bead
        # Head local index inside this motor
        head_within = within_motor - N
        return motor_idx * 2 * H + head_within

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        assert sim is not None
        if self.p.n_motors_per_cell == 0:
            return

        read_snap = sim.state.get_snapshot()
        if read_snap.communicator.rank != 0:
            return

        pos = np.asarray(read_snap.particles.position, dtype=np.float64).copy()
        bg = np.asarray(read_snap.bonds.group, dtype=np.int64).copy()
        bt = np.asarray(read_snap.bonds.typeid, dtype=np.uint32).copy()
        bond_type_names = list(read_snap.bonds.types)

        attach_bin_typeids = [
            i for i, name in enumerate(bond_type_names)
            if name.startswith("cortex_myosin_attach_b")
        ]
        if not attach_bin_typeids:
            raise RuntimeError(
                "snap.bonds.types missing cortex_myosin_attach_b* types; "
                "did you forget extend_state_with_cortex_myosin + "
                "register_cortex_myosin_bond_params?"
            )

        is_attach = np.isin(bt, np.asarray(attach_bin_typeids, dtype=bt.dtype))
        attach_bonds = bg[is_attach]
        other_bg = bg[~is_attach]
        other_bt = bt[~is_attach]

        bin_r0 = cortex_myosin_attach_bin_rest_lengths(
            self.p.n_bins, self.p.head_actin_max_bind_dist
        )
        bin_width = self.p.head_actin_max_bind_dist / self.p.n_bins

        # ---- Step 1: Bell-Evans unbinding ----
        if attach_bonds.shape[0] > 0:
            head_tags = attach_bonds[:, 0]
            actin_tags = attach_bonds[:, 1]
            r_head = pos[head_tags]
            r_actin = pos[actin_tags]
            r = np.linalg.norm(r_head - r_actin, axis=1)
            bond_bins = bt[is_attach] - attach_bin_typeids[0]
            r0_per_bond = bin_r0[bond_bins]
            F_mag = self.p.k_head_actin * np.clip(r - r0_per_bond, 0.0, None)
            k_off = _bell_evans_k_off(
                F_mag, self.p.head_actin_k_off0, self.p.head_actin_x_beta, self.kT,
            )
            p_break = 1.0 - np.exp(-k_off * self.p.batch_dt)
            u = self._rng.uniform(0.0, 1.0, size=attach_bonds.shape[0])
            broke = u < p_break
            n_broke = int(broke.sum())
            self._n_break_total += n_broke
            if n_broke > 0:
                broken_head_tags = head_tags[broke]
                for ht in broken_head_tags:
                    h_local = self._head_local_from_tag(int(ht))
                    if h_local >= 0:
                        self._head_bound_to_actin[h_local] = -1
                attach_bonds = attach_bonds[~broke]
                bond_bins = bond_bins[~broke]
        else:
            bond_bins = np.empty((0,), dtype=np.int64)

        # ---- Step 2: binding (KDTree) ----
        from scipy.spatial import cKDTree
        unbound_head_locals = np.flatnonzero(self._head_bound_to_actin < 0)
        if unbound_head_locals.size > 0:
            unbound_head_tags = np.array(
                [self._head_global_tag(int(h)) for h in unbound_head_locals],
                dtype=np.int64,
            )
            r_heads = pos[unbound_head_tags]
            r_actin_all = pos[:self.n_cortex_actin]
            tree = cKDTree(r_actin_all)
            nbr_lists = tree.query_ball_point(
                r_heads, r=self.p.head_actin_max_bind_dist
            )
            p_bind = 1.0 - np.exp(-self.p.head_actin_k_on * self.p.batch_dt)
            u2 = self._rng.uniform(0.0, 1.0, size=unbound_head_locals.size)
            new_bonds_list = []
            new_bins_list = []
            for k, nbrs in enumerate(nbr_lists):
                if len(nbrs) == 0 or u2[k] >= p_bind:
                    continue
                nbrs_arr = np.asarray(nbrs, dtype=np.int64)
                d_nbrs = np.linalg.norm(
                    r_actin_all[nbrs_arr] - r_heads[k], axis=1
                )
                nearest_local = int(np.argmin(d_nbrs))
                actin_tag = int(nbrs_arr[nearest_local])
                d_use = float(d_nbrs[nearest_local])
                idx_bin = int(min(self.p.n_bins - 1, max(0, int(d_use / bin_width))))
                head_tag = int(unbound_head_tags[k])
                new_bonds_list.append((head_tag, actin_tag))
                new_bins_list.append(idx_bin)
                self._head_bound_to_actin[unbound_head_locals[k]] = actin_tag
            if new_bonds_list:
                new_bonds = np.array(new_bonds_list, dtype=np.int64)
                new_bins = np.array(new_bins_list, dtype=np.int64)
                attach_bonds = np.concatenate(
                    [attach_bonds, new_bonds], axis=0
                )
                bond_bins = np.concatenate([bond_bins, new_bins], axis=0)
                self._n_bind_total += int(new_bonds.shape[0])

        # ---- Step 3: D6 Hill stepping (advance bin) ----
        # For each engaged attach bond, advance bin index "inward" by
        # the Hill v(F) · batch_dt / bin_width. Clamped to bin 0.
        if attach_bonds.shape[0] > 0:
            head_tags = attach_bonds[:, 0]
            actin_tags = attach_bonds[:, 1]
            r_head = pos[head_tags]
            r_actin = pos[actin_tags]
            r = np.linalg.norm(r_head - r_actin, axis=1)
            r0_per_bond = bin_r0[bond_bins]
            F_mag = self.p.k_head_actin * np.clip(r - r0_per_bond, 0.0, None)
            v_step = hill_velocity_clamped(
                F_mag, v0=self.p.v0_per_head, F_stall=self.p.F_stall_per_head,
                a_over_F_stall=self.p.a_over_F_stall,
            )
            d_bin = (v_step * self.p.batch_dt) / bin_width
            new_bins = np.clip(
                bond_bins - np.round(d_bin).astype(np.int64),
                0, self.p.n_bins - 1,
            )
            n_advanced = int((new_bins != bond_bins).sum())
            self._n_step_advances_total += n_advanced
            bond_bins = new_bins

        # ---- Rebuild bonds.group + bonds.typeid + write_snap ----
        attach_typeids = (
            np.asarray(bond_bins, dtype=np.uint32) + attach_bin_typeids[0]
        )
        if other_bg.shape[0] > 0:
            new_bg = np.concatenate(
                [attach_bonds.astype(np.uint32), other_bg.astype(np.uint32)],
                axis=0,
            )
            new_bt = np.concatenate(
                [attach_typeids, other_bt.astype(np.uint32)]
            )
        else:
            new_bg = attach_bonds.astype(np.uint32)
            new_bt = attach_typeids

        write_snap = hoomd.Snapshot()
        N_part = int(read_snap.particles.N)
        write_snap.particles.N = N_part
        write_snap.particles.types = list(read_snap.particles.types)
        write_snap.particles.typeid[:] = np.asarray(read_snap.particles.typeid)
        write_snap.particles.position[:] = pos
        write_snap.particles.velocity[:] = np.asarray(read_snap.particles.velocity)
        write_snap.particles.mass[:] = np.asarray(read_snap.particles.mass)
        write_snap.particles.image[:] = np.asarray(read_snap.particles.image)
        box = read_snap.configuration.box
        write_snap.configuration.box = list(box)

        write_snap.bonds.N = int(new_bg.shape[0])
        write_snap.bonds.types = list(bond_type_names)
        if new_bg.shape[0] > 0:
            write_snap.bonds.group[:] = new_bg
            write_snap.bonds.typeid[:] = new_bt

        for grp_name in ("angles", "dihedrals", "impropers"):
            src = getattr(read_snap, grp_name)
            dst = getattr(write_snap, grp_name)
            if int(src.N) > 0:
                dst.N = int(src.N)
                dst.types = list(src.types)
                dst.group[:] = np.asarray(src.group)
                dst.typeid[:] = np.asarray(src.typeid)

        sim.state.set_snapshot(write_snap)
        self._steps_run += 1

    # --- Introspection ---
    @property
    def n_engaged(self) -> int:
        return int((self._head_bound_to_actin >= 0).sum())

    @property
    def n_break_total(self) -> int:
        return self._n_break_total

    @property
    def n_bind_total(self) -> int:
        return self._n_bind_total

    @property
    def n_step_advances_total(self) -> int:
        return self._n_step_advances_total

    @property
    def steps_run(self) -> int:
        return self._steps_run


def make_cortex_myosin_updater(
    *,
    p_myo: ResolvedCortexMyosin,
    layout: CortexMyosinLayout,
    kT: float,
    n_cortex_actin: int,
    seed_offset: int = 3,
) -> tuple[MyosinStepUpdater, hoomd.update.CustomUpdater]:
    action = MyosinStepUpdater(
        p_myo=p_myo, layout=layout, kT=kT,
        n_cortex_actin=n_cortex_actin, seed_offset=seed_offset,
    )
    updater = hoomd.update.CustomUpdater(
        action=action, trigger=hoomd.trigger.Periodic(p_myo.batch_steps)
    )
    return action, updater
