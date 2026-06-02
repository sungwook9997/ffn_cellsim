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
    # head_actin_max_bind_dist: BOND R0 BIN-SCHEME max (largest head-bead bond r0
    # accepted). Segment-derived: √((ℓ₀/2)² + capture_perp²) ≈ 327nm at ℓ₀=500nm,
    # capture_perp=210nm. Larger than the legacy 50nm because heads now sit
    # head_off≈200nm laterally from the actin segment they bind (Option C).
    head_actin_max_bind_dist: float   # m   max head-bead bond r0 (binning ceiling)
    # head_actin_capture_perp: PERPENDICULAR distance from a head to an actin
    # segment LINE that counts as binding-eligible (segment-projection,
    # KU-3.5 fix 2026-05-29 Option C). Set to head_rest_length + small slack
    # (=210nm for head_off=200nm) — the physical head reach via its spring.
    head_actin_capture_perp: float    # m  segment perp capture (derived from head_off+slack)

    # D2 batch
    batch_steps: int             # BAOAB steps between MyosinStepUpdater ticks
    dt: float                    # host sim dt [s]
    n_bins: int                  # per-r0 binning count for dynamic attach bonds

    # Contractile stepping mechanism (KU-3.5 grip-walk redesign 2026-05-31).
    # "binned_r0"  — legacy lumped proxy: a "step" relabels the head-actin bond
    #                to a lower-r0 bin on the SAME actin bead (transports no
    #                material → the diagnosed KU-3.5 floor). DEFAULT (A/B base).
    # "grip_walk"  — AFINES pos_a_end: the head's grip point WALKS toward the
    #                filament minus end (re-targets the bond to downstream beads),
    #                a continuous commanded sub-bead stretch s_grip carries force
    #                (r0_eff = max(r − s_grip, 0)); sustained contraction.
    # PI-ratified 2026-05-31 (decisions 1.4 continuous sub-bead / 2 Option-A
    # polarity / 3 bipolar antiparallel gate / 4 no param change). Opt-in.
    stepping_mode: str           # "binned_r0" | "grip_walk"

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


def resolve_cortex_myosin(
    cfg: dict, *, dt: float, R_cell: float | None = None
) -> ResolvedCortexMyosin:
    """Resolve cortical-myosin config block with D5 boundary + CFL gates.

    ``R_cell`` is required only when ``mesoscale_force_scaling`` is enabled
    (KU-3.5 Route B) — it sets the native minifilament count from areal density.
    """
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
        head_actin_capture_perp=float(
            cfg.get("head_actin_capture_perp", float(cfg["head_rest_length"]) + 10.0e-9)
        ),
        batch_steps=int(cfg["batch_steps"]),
        dt=float(dt),
        n_bins=int(cfg.get("n_bins", 10)),
        stepping_mode=str(cfg.get("stepping_mode", "binned_r0")),
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
    if p.stepping_mode not in ("binned_r0", "grip_walk"):
        raise ValueError(
            f"stepping_mode must be 'binned_r0' or 'grip_walk'; "
            f"got {p.stepping_mode!r}"
        )

    # Mesoscale myosin FORCE scaling (KU-3.5 Route B, PI-ratified 2026-05-31).
    # The ×40 mesoscopic coarse-graining reduces the MOTOR count (native
    # ~areal_density·4πR² minifilaments → n_motors_per_cell), so each EFFECTIVE
    # minifilament stands in for `factor = native/effective` native ones. To
    # carry the aggregate contractile force — the KU-3.5 floor's SECOND cause,
    # a ~63× force-budget gap at n_motors=100 (impl doc §4b) — model the
    # effective minifilament as `factor` native ones IN PARALLEL:
    #   • springs add in parallel  → k_head_spring, k_head_actin ×= factor
    #     (k_backbone follows, since k_backbone = k_head_spring·k_backbone_factor)
    #   • stall force adds          → F_stall_per_head ×= factor
    #     (so s_grip can still grow to ℓ₀: s_grip_max = F_stall/k is INVARIANT
    #      under the scaling — the stretch stays physical)
    #   • Bell-Evans load-sensitivity per NATIVE head is preserved → x_β ÷= factor
    #     (a parallel bundle shares load: each native bond bears F_eff/factor, so
    #      k_off(F_eff; x_β/factor) = k_off(F_native; x_β) — without this the
    #      ×factor force would strip the effective head ~exp(factor)× faster).
    # `factor` is DERIVED (density·area/n_motors), grid-invariant, not tuned —
    # satisfies the Magic-Number Block. Opt-in + grip_walk only (binned_r0 legacy
    # proxy stays byte-identical). v0 (per-motor kinetics) is unchanged.
    p.extras["mesoscale_force_scaling"] = False
    if bool(cfg.get("mesoscale_force_scaling", False)) and p.stepping_mode == "grip_walk":
        if R_cell is None or not (math.isfinite(R_cell) and R_cell > 0.0):
            raise ValueError(
                "mesoscale_force_scaling requires a finite R_cell > 0 (to derive "
                "the native minifilament count from areal density × surface area)."
            )
        density_per_um2 = float(cfg.get("areal_density_per_um2", 3.0))  # Salbreux 2012
        native_n_motors = density_per_um2 * 1.0e12 * 4.0 * math.pi * R_cell ** 2
        factor = native_n_motors / max(p.n_motors_per_cell, 1)
        p.k_head_spring *= factor
        p.k_head_actin *= factor
        p.F_stall_per_head *= factor
        p.head_actin_x_beta /= factor
        p.extras.update(
            mesoscale_force_scaling=True,
            mesoscale_force_factor=float(factor),
            native_n_motors=float(native_n_motors),
        )

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
    cortex_positions: np.ndarray | None = None,
    cortex_tangents: np.ndarray | None = None,
    beads_per_filament: int | None = None,
) -> CortexMyosinLayout:
    """Place ``n_motors_per_cell`` minifilaments on the cortex shell.

    Two placement modes:

    **Actin-aware (preferred)** — when ``cortex_positions`` (F·N, 3) +
    ``cortex_tangents`` (F, 3) + ``beads_per_filament`` are provided:

    1. Pick ``M`` random cortex actin BEADS as minifilament centres
       (without replacement, so no two minifilaments share a center bead
       — gives overlap-avoidance at construction).
    2. Backbone axis ``u`` = the local actin filament tangent (the
       minifilament lies along the actin it sits on, biologically
       meaningful — myosin minifilaments crosslink cortical actin).
    3. Heads offset in the **tangent plane** perpendicular to the
       backbone: lateral direction ``w = normalize(cross(n, u))`` where
       ``n`` is the local outward shell normal. + polarity heads at
       ``+head_rest_length · w``, − heads at ``−head_rest_length · w``.
       Heads stay IN the cortex shell (radius ≈ R_cell), within reach of
       the actin filament they sit on (and adjacent ones).

    **Legacy sphere-random** — when cortex info is omitted: keeps the
    original (buggy for binding — heads radially off shell, 824 nm median
    from actin) placement for back-compat with tests that don't have
    a cortex topology in hand. KU-3.5/3.1/3.18 production needs the
    actin-aware mode (PI ratified 2026-05-29 — Option C).
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

    actin_aware = (
        cortex_positions is not None
        and cortex_tangents is not None
        and beads_per_filament is not None
    )

    if actin_aware:
        # Pick M actin beads as minifilament centers with PAIRWISE SPACING
        # ≥ backbone_length + slack, so two minifilaments don't overlap
        # (intra-myosin LJ would blow up the integrator otherwise — the
        # backbones span ~700 nm + heads ±head_off laterally).
        n_actin = cortex_positions.shape[0]
        if M > n_actin:
            raise ValueError(
                f"n_motors_per_cell={M} > n_cortex_actin={n_actin}; "
                "cannot place each minifilament at a unique actin bead.")
        from scipy.spatial import cKDTree
        # Greedy: shuffle, accept if far enough from already-picked.
        min_sep = float(L) + 100.0e-9   # backbone_length + 100 nm slack
        order = rng.permutation(n_actin)
        picked: list[int] = []
        picked_pos: list[np.ndarray] = []
        for idx in order:
            if not picked_pos:
                picked.append(int(idx)); picked_pos.append(cortex_positions[idx])
                if len(picked) == M: break
                continue
            d = np.linalg.norm(
                np.stack(picked_pos) - cortex_positions[idx], axis=1
            ).min()
            if d >= min_sep:
                picked.append(int(idx)); picked_pos.append(cortex_positions[idx])
                if len(picked) == M: break
        if len(picked) < M:
            raise RuntimeError(
                f"Could not place {M} minifilaments with ≥{min_sep*1e9:.0f}nm "
                f"pairwise spacing on {n_actin} actin beads (got {len(picked)}). "
                "Reduce n_motors_per_cell or relax spacing.")
        bead_choice = np.asarray(picked, dtype=np.int64)
        centers = cortex_positions[bead_choice].copy()
        # Filament index for each chosen bead (contiguous bead blocks).
        fil_idx = bead_choice // beads_per_filament
        axes = cortex_tangents[fil_idx].copy()
        axes = axes / np.linalg.norm(axes, axis=1, keepdims=True).clip(min=1e-30)
        normals = centers / np.linalg.norm(centers, axis=1, keepdims=True).clip(min=1e-30)
        # In-plane lateral direction perpendicular to backbone, in tangent plane.
        lateral = np.cross(normals, axes)
        lateral = lateral / np.linalg.norm(lateral, axis=1, keepdims=True).clip(min=1e-30)
    else:
        centers = _sample_sphere_surface(rng, M, R_cell)
        normals = centers / R_cell
        e1, e2 = _tangent_plane_basis(normals)
        phi = rng.uniform(0.0, 2.0 * math.pi, M)
        axes = (
            np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2
        )
        axes = axes / np.linalg.norm(axes, axis=1, keepdims=True).clip(min=1e-30)
        lateral = normals  # legacy: heads offset RADIALLY (the binding bug)

    # DIAGNOSTIC (γ-flat root-cause, 2026-06-02; env-gated, default OFF so the
    # physical random/filament-tangent orientation is unchanged): override the
    # minifilament axes with a COHERENT global meridional director field (the
    # projection of +z onto each tangent plane) instead of the
    # isotropic/filament-tangent orientation. Tests whether force ISOTROPY is the
    # cortical-tension wall — aligned bipolar motors sum into net great-circle
    # tension, isotropic ones cancel. NOT a physical claim that cortex myosin is
    # globally aligned; a probe to separate isotropy from transmission loss.
    import os
    if os.environ.get("FFN_MYOSIN_ALIGN") == "meridional":
        zg = np.array([0.0, 0.0, 1.0])
        proj = zg[None, :] - (normals @ zg)[:, None] * normals
        nrm = np.linalg.norm(proj, axis=1, keepdims=True)
        e1_pole = _tangent_plane_basis(normals)[0]
        bad = nrm[:, 0] < 1e-6                              # near the poles
        axes = np.where(bad[:, None], e1_pole, proj / nrm.clip(min=1e-30))
        axes = axes / np.linalg.norm(axes, axis=1, keepdims=True).clip(min=1e-30)
        lateral = np.cross(normals, axes)
        lateral = lateral / np.linalg.norm(lateral, axis=1, keepdims=True).clip(min=1e-30)

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
        w = lateral[m]    # tangent-plane lateral (actin-aware) OR n (legacy)
        # Backbone
        for i in range(N):
            positions[m, i] = cm + backbone_offsets[i] * u
        # + heads (one lateral direction)
        for i in range(H):
            positions[m, N + i] = cm + head_axis_offsets[i] * u + head_off * w
        # − heads (opposite lateral)
        for i in range(H):
            positions[m, N + H + i] = cm + head_axis_offsets[i] * u - head_off * w

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
        cortex_bond_groups: np.ndarray | None = None,
        ell0_cortex: float | None = None,
        cortex_beads_per_filament: int | None = None,
        seed_offset: int = 3,
    ) -> None:
        super().__init__()
        self.p = p_myo
        self.layout = layout
        self.kT = float(kT)
        self.n_cortex_actin = int(n_cortex_actin)
        self._rng = np.random.default_rng(p_myo.seed + seed_offset)
        # Grip-walk geometry (KU-3.5 redesign). Required iff stepping_mode is
        # "grip_walk": the bead-tag ↔ (filament, position) map + the bead
        # spacing ℓ₀ that one walked sub-bead step is measured in.
        self.stepping_mode = str(p_myo.stepping_mode)
        self._ell0_cortex = (
            float(ell0_cortex) if ell0_cortex is not None else None
        )
        self._cortex_beads_per_filament = (
            int(cortex_beads_per_filament)
            if cortex_beads_per_filament is not None else None
        )
        if self.stepping_mode == "grip_walk":
            if self._ell0_cortex is None or self._cortex_beads_per_filament is None:
                raise ValueError(
                    "grip_walk stepping_mode requires ell0_cortex + "
                    "cortex_beads_per_filament (the fixed-N bead-tag ↔ "
                    "(filament, pos) map); got None."
                )
        # Segment-projection binding precompute (KU-3.5 Option C). When
        # cortex_bond_groups provided: for each actin bead, list the segment
        # indices it participates in (1 or 2 segments per bead in a chain).
        # Segments use the bond-group ordering; each segment is the pair
        # (bg[s,0], bg[s,1]) which are two consecutive cortex actin beads.
        self._cortex_bond_groups = (
            np.asarray(cortex_bond_groups, dtype=np.int64)
            if cortex_bond_groups is not None else None
        )
        if self._cortex_bond_groups is not None:
            # Per-bead adjacency: which segments touch each bead.
            adj: list[list[int]] = [[] for _ in range(int(n_cortex_actin))]
            for s, (a, b) in enumerate(self._cortex_bond_groups):
                if 0 <= a < n_cortex_actin: adj[int(a)].append(s)
                if 0 <= b < n_cortex_actin: adj[int(b)].append(s)
            self._bead_to_segs = [np.asarray(v, dtype=np.int64) for v in adj]
        else:
            self._bead_to_segs = None

        # Per-head state — bound to which actin tag (−1 if free).
        n_heads_total = 2 * p_myo.n_heads_per_side * p_myo.n_motors_per_cell
        self._head_bound_to_actin = np.full(n_heads_total, -1, dtype=np.int64)
        # Per-head Hill-stepping FRACTIONAL ACCUMULATOR (KU-3.5 단계-4 fix
        # 2026-05-29). Without this the original `round(d_bin)` was always 0
        # at the model's timescale (v·batch_dt = 70 pm ≪ bin_width = 33 nm)
        # → no stepping → no contraction. Accumulate d_bin (float) per tick;
        # advance an integer bin only when accum ≥ 1, retain the remainder.
        # Reset on unbind so re-binding starts fresh.
        self._head_step_accum = np.zeros(n_heads_total, dtype=np.float64)
        # --- Grip-walk per-head state (KU-3.5 redesign, used iff grip_walk) ---
        # The decomposition of _head_bound_to_actin (global bead tag) into the
        # (filament, position-within-filament) that the walk increments toward
        # the minus end (Option A: minus = bead 0), plus the continuous
        # commanded sub-bead stretch s_grip ∈ [0, ℓ₀) (AFINES pos_a_end).
        self._head_bound_bead_pos = np.full(n_heads_total, -1, dtype=np.int64)
        self._head_bound_filament = np.full(n_heads_total, -1, dtype=np.int64)
        self._head_grip_s = np.zeros(n_heads_total, dtype=np.float64)
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

    # --- Grip-walk geometry helpers (fixed-N cortex; Option-A polarity) ---
    def _tag_to_fil_pos(self, bead_tag: int) -> tuple[int, int]:
        """Decompose a cortex actin bead tag → (filament, pos-in-filament).

        Fixed-N contiguous block ``[0, n_cortex_actin)``:
        ``filament = tag // N``, ``pos = tag % N`` (exact O(1)).
        """
        nb = self._cortex_beads_per_filament
        return bead_tag // nb, bead_tag % nb

    def _bead_tag(self, fil: int, pos: int) -> int:
        """Inverse of :meth:`_tag_to_fil_pos` (fixed-N)."""
        return fil * self._cortex_beads_per_filament + pos

    def _walk_toward_minus(self, pos_j: int, n: int) -> int:
        """Advance a grip position ``n`` beads toward the minus end.

        Option A (PI-ratified 2026-05-31): minus end = bead ``j=0``, so
        walking decrements the position and clamps at 0 (AFINES at-minus-end
        latch — the head dwells, keeps pulling, never walks off the filament).
        """
        return max(0, pos_j - n)

    def _bipolar_accepts(
        self, head_local: int, bead_tag: int, pos: np.ndarray
    ) -> bool:
        """Bipolar sidedness gate (PI decision 3); grip_walk only.

        Accept a candidate head→bead bind iff (a) the candidate filament's
        minus-end-ward tangent ``m̂`` is oriented with the rod axis ``û`` per
        the head's side (+ side: ``m̂·û > 0``; − side: ``< 0``), AND (b) the
        motor's OTHER side does not already grip this filament (the
        different-filament clause that removes the zero-dipole degeneracy).
        """
        H = self.p.n_heads_per_side
        nb = self._cortex_beads_per_filament
        motor_idx = head_local // (2 * H)
        head_within = head_local % (2 * H)
        side_plus = head_within < H
        fil, pos_j = self._tag_to_fil_pos(int(bead_tag))
        # Minus-end-ward tangent m̂ from bead positions (Option A: minus=bead 0).
        if pos_j > 0:
            m_vec = pos[self._bead_tag(fil, pos_j - 1)] - pos[bead_tag]
        elif nb > 1:
            # At the minus end: extrapolate the minus direction from j=1→j=0.
            m_vec = pos[bead_tag] - pos[self._bead_tag(fil, 1)]
        else:
            return True  # single-bead filament: no defined polarity
        norm = float(np.linalg.norm(m_vec))
        if norm <= 0.0:
            return True
        dot = float((m_vec / norm) @ self.layout.axes[motor_idx])
        if side_plus and not (dot > 0.0):
            return False
        if (not side_plus) and not (dot < 0.0):
            return False
        base = motor_idx * 2 * H
        other = (
            slice(base + H, base + 2 * H) if side_plus
            else slice(base, base + H)
        )
        if np.any(self._head_bound_filament[other] == fil):
            return False
        return True

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

        # Per-cortex-bead total bonded-degree budget (HOOMD nlist exclusion cap).
        # MAX_HEADS_PER_BEAD below caps only MYOSIN attach bonds and is blind to
        # backbone + FA clutch + xlink, so 3 heads on an already-loaded bead can
        # push its raw degree to 8 -> HOOMD's "Too many bonds to process
        # exclusions" crash (hard limit 7; the pre-existing long-run v4 crash).
        # Also enforce the SHARED budget over ALL current bonds on the bead.
        # Margin of 1 under 7 covers sequential-binder timing. (bg is reassigned
        # to the cortex backbone groups later, so compute this from the full
        # bonds.group HERE.)
        _flat_bg = bg.reshape(-1)
        _cortex_bead_degree = np.bincount(
            _flat_bg[_flat_bg < self.n_cortex_actin],
            minlength=self.n_cortex_actin,
        ).astype(np.int64)
        _MAX_CORTEX_BEAD_DEGREE = 6

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
                        self._head_step_accum[h_local] = 0.0  # fresh on re-bind
                        self._head_bound_bead_pos[h_local] = -1
                        self._head_bound_filament[h_local] = -1
                        self._head_grip_s[h_local] = 0.0
                attach_bonds = attach_bonds[~broke]
                bond_bins = bond_bins[~broke]
        else:
            bond_bins = np.empty((0,), dtype=np.int64)

        # ---- Step 2: binding ----
        # Two modes:
        #  - segment-projection (KU-3.5 Option C, when cortex_bond_groups
        #    is set at construction): a head is eligible if its perpendicular
        #    distance to an actin SEGMENT line is ≤ capture_perp AND its
        #    projection lies within the segment (or the endpoint is within
        #    capture_perp). Bind to the segment's NEARER endpoint bead.
        #  - legacy bead-center (when no bond_groups): KDTree within
        #    max_bind_dist. Preserved for back-compat with old tests.
        from scipy.spatial import cKDTree
        unbound_head_locals = np.flatnonzero(self._head_bound_to_actin < 0)
        if unbound_head_locals.size > 0:
            unbound_head_tags = np.array(
                [self._head_global_tag(int(h)) for h in unbound_head_locals],
                dtype=np.int64,
            )
            r_heads = pos[unbound_head_tags]
            r_actin_all = pos[:self.n_cortex_actin]
            p_bind = 1.0 - np.exp(-self.p.head_actin_k_on * self.p.batch_dt)
            u2 = self._rng.uniform(0.0, 1.0, size=unbound_head_locals.size)
            new_bonds_list = []
            new_bins_list = []
            if self._cortex_bond_groups is not None and self._bead_to_segs is not None:
                # Segment-projection mode.
                bg = self._cortex_bond_groups
                tree = cKDTree(r_actin_all)
                # Search radius covers any bead whose adjacent segment could
                # be perpendicular-eligible: head can be √(capture_perp² +
                # (ℓ₀/2)²) from a bead's center yet still be perp-close to
                # the segment midpoint. Use head_actin_max_bind_dist as that
                # geometric ceiling (already sized).
                nbr_lists = tree.query_ball_point(
                    r_heads, r=self.p.head_actin_max_bind_dist
                )
                cap = self.p.head_actin_capture_perp
                # Per-actin-bead binding cap. Each actin bead already carries
                # ~2 backbone bonds + ~2 1-3 angle exclusions; HOOMD's nlist
                # has a compile-time max of 7 exclusions per particle. Cap
                # head-actin attach bonds at 3 per bead → ≤7 total. Also
                # biologically defensible (steric: one cortical actin bead
                # ≈ 500nm × ~90 monomers cannot host arbitrarily many myosin
                # heads simultaneously).
                MAX_HEADS_PER_BEAD = 3
                bead_attach_count = np.zeros(self.n_cortex_actin, dtype=np.int32)
                if attach_bonds.shape[0] > 0:
                    for ab in attach_bonds[:, 1]:
                        ab = int(ab)
                        if 0 <= ab < self.n_cortex_actin:
                            bead_attach_count[ab] += 1
                for k, nbrs in enumerate(nbr_lists):
                    if len(nbrs) == 0 or u2[k] >= p_bind:
                        continue
                    h = r_heads[k]
                    # Collect candidate segments from each nearby bead's adjacency.
                    cand_segs = np.unique(np.concatenate(
                        [self._bead_to_segs[int(b)] for b in nbrs]
                        + [np.empty(0, dtype=np.int64)]
                    ))
                    best_perp = np.inf; best_bead = -1; best_d_use = 0.0
                    for s in cand_segs:
                        a_idx = int(bg[s, 0]); b_idx = int(bg[s, 1])
                        A = r_actin_all[a_idx]; B = r_actin_all[b_idx]
                        seg = B - A; L2 = float(seg @ seg)
                        if L2 <= 0.0: continue
                        t = float((h - A) @ seg) / L2          # axial parameter
                        t_cl = max(0.0, min(1.0, t))           # clamp to segment
                        closest = A + t_cl * seg
                        perp = float(np.linalg.norm(h - closest))
                        if perp <= cap and perp < best_perp:
                            # Bind to the NEARER endpoint bead.
                            dA = float(np.linalg.norm(h - A))
                            dB = float(np.linalg.norm(h - B))
                            if dA <= dB:
                                best_bead = a_idx; best_d_use = dA
                            else:
                                best_bead = b_idx; best_d_use = dB
                            best_perp = perp
                    if best_bead < 0:
                        continue
                    # Enforce per-bead exclusion-safe cap (myosin-specific count
                    # AND the shared total-degree budget across all bond sources).
                    if bead_attach_count[best_bead] >= MAX_HEADS_PER_BEAD:
                        continue
                    if _cortex_bead_degree[best_bead] >= _MAX_CORTEX_BEAD_DEGREE:
                        continue
                    head_local = int(unbound_head_locals[k])
                    bw_fil = bw_pos = -1
                    if self.stepping_mode == "grip_walk":
                        # ---- Bipolar sidedness gate (PI decision 3) ----
                        # Organize the +/− head sets into a net contractile
                        # dipole (Stam–Hocky): + side binds filaments whose
                        # minus-end-ward tangent m̂ points along the rod axis û,
                        # − side binds the antiparallel ones; and the two sides
                        # must grip DIFFERENT filaments (no zero-dipole
                        # degeneracy). m̂ is computed from bead positions under
                        # Option-A polarity (minus = bead 0).
                        if not self._bipolar_accepts(
                            head_local, best_bead, pos
                        ):
                            continue
                        bw_fil, bw_pos = self._tag_to_fil_pos(best_bead)
                    # Bin r0 from actual head-bead distance (clamp to bin range).
                    d_use = min(best_d_use, self.p.head_actin_max_bind_dist - 1e-12)
                    idx_bin = int(min(self.p.n_bins - 1,
                                      max(0, int(d_use / bin_width))))
                    head_tag = int(unbound_head_tags[k])
                    new_bonds_list.append((head_tag, best_bead))
                    new_bins_list.append(idx_bin)
                    self._head_bound_to_actin[head_local] = best_bead
                    if self.stepping_mode == "grip_walk":
                        # Initialise the grip at the bound bead, zero stretch
                        # (r0_eff = r at bind → force-free construction, §6.2).
                        self._head_bound_filament[head_local] = bw_fil
                        self._head_bound_bead_pos[head_local] = bw_pos
                        self._head_grip_s[head_local] = 0.0
                    bead_attach_count[best_bead] += 1
                    _cortex_bead_degree[best_bead] += 1
            else:
                # Legacy bead-center mode.
                tree = cKDTree(r_actin_all)
                nbr_lists = tree.query_ball_point(
                    r_heads, r=self.p.head_actin_max_bind_dist
                )
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
                    idx_bin = int(min(self.p.n_bins - 1,
                                      max(0, int(d_use / bin_width))))
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

        # ---- Step 3: D6 Hill stepping ----
        if attach_bonds.shape[0] > 0:
            head_tags = attach_bonds[:, 0]
            actin_tags = attach_bonds[:, 1]
            r = np.linalg.norm(pos[head_tags] - pos[actin_tags], axis=1)
            # Map head_tags → head local indices (vectorised arithmetic
            # inverse of _head_global_tag).
            H = self.p.n_heads_per_side
            N = self.p.n_backbone
            per_motor = self.p.n_particles_per_motor
            offset = head_tags - self.layout.motor_tag_start
            head_locals = (offset // per_motor) * (2 * H) + (
                (offset % per_motor) - N
            )

            if self.stepping_mode == "binned_r0":
                # Legacy lumped proxy (FRACTIONAL ACCUMULATOR — 단계-4 fix):
                # accumulate Hill v(F)·batch_dt/bin_width into a per-head
                # accumulator and advance integer bins when accum ≥ 1, clamped
                # to bin 0. This relabels the bond's r0 on the SAME bead — the
                # diagnosed KU-3.5 lumped proxy (transports no material).
                r0_per_bond = bin_r0[bond_bins]
                F_mag = self.p.k_head_actin * np.clip(r - r0_per_bond, 0.0, None)
                v_step = hill_velocity_clamped(
                    F_mag, v0=self.p.v0_per_head,
                    F_stall=self.p.F_stall_per_head,
                    a_over_F_stall=self.p.a_over_F_stall,
                )
                d_bin = (v_step * self.p.batch_dt) / bin_width  # float, per-bond
                self._head_step_accum[head_locals] += d_bin
                adv = np.floor(
                    self._head_step_accum[head_locals]
                ).astype(np.int64)
                self._head_step_accum[head_locals] -= adv.astype(np.float64)
                new_bins = np.clip(bond_bins - adv, 0, self.p.n_bins - 1)
                n_advanced = int((new_bins != bond_bins).sum())
                self._n_step_advances_total += n_advanced
                bond_bins = new_bins
            else:
                # ---- grip_walk (AFINES pos_a_end; PI-ratified 2026-05-31) ----
                # The commanded sub-bead stretch s_grip carries the force via
                # r0_eff = max(r − s_grip, 0)  ⇒  F = k·min(s_grip, r). Crucially
                # r0_eff is RE-DERIVED from the current r every tick, so the
                # force does NOT relax to 0 as the bead approaches (the proxy's
                # failure); s_grip GROWS by Hill v(F)·batch_dt (pos_a_end), and
                # when it overflows one bead spacing ℓ₀ the grip RE-TARGETS to
                # the next minus-ward bead (material transport).
                ell0 = self._ell0_cortex
                s = self._head_grip_s[head_locals]            # commanded stretch
                F_mag = self.p.k_head_actin * np.clip(
                    np.minimum(s, r), 0.0, None
                )
                v_step = hill_velocity_clamped(
                    F_mag, v0=self.p.v0_per_head,
                    F_stall=self.p.F_stall_per_head,
                    a_over_F_stall=self.p.a_over_F_stall,
                )
                s_new = s + v_step * self.p.batch_dt          # pos_a_end advance
                # Per-head overflow → minus-ward re-target (with the shared
                # degree budget re-checked AT WALK TIME, not just bind time).
                for row in range(attach_bonds.shape[0]):
                    h = int(head_locals[row])
                    s_h = float(s_new[row])
                    if s_h < ell0:
                        self._head_grip_s[h] = s_h
                        continue
                    fil = int(self._head_bound_filament[h])
                    pos_j = int(self._head_bound_bead_pos[h])
                    while s_h >= ell0:
                        if pos_j <= 0:
                            # AFINES minus-end latch: dwell, hold max sub-bead
                            # stretch (just under ℓ₀), keep pulling.
                            s_h = ell0 * (1.0 - 1e-9)
                            break
                        new_j = self._walk_toward_minus(pos_j, 1)
                        new_tag = self._bead_tag(fil, new_j)
                        if (_cortex_bead_degree[new_tag]
                                >= _MAX_CORTEX_BEAD_DEGREE):
                            # Downstream bead full → defer the walk (clean
                            # stall), hold just under overflow, retry next tick.
                            s_h = ell0 * (1.0 - 1e-9)
                            break
                        old_tag = self._bead_tag(fil, pos_j)
                        _cortex_bead_degree[old_tag] = max(
                            0, _cortex_bead_degree[old_tag] - 1
                        )
                        _cortex_bead_degree[new_tag] += 1
                        pos_j = new_j
                        s_h -= ell0
                        attach_bonds[row, 1] = new_tag
                        self._head_bound_to_actin[h] = new_tag
                        self._head_bound_bead_pos[h] = pos_j
                        self._n_step_advances_total += 1
                    self._head_grip_s[h] = s_h
                # Requantize r0_eff into the bin coordinate (HOOMD Harmonic
                # carries r0 per bond TYPE → pick the bin whose r0 ≈ the
                # commanded rest length). Recompute r at the (possibly
                # re-targeted) bead.
                r2 = np.linalg.norm(
                    pos[head_tags] - pos[attach_bonds[:, 1]], axis=1
                )
                r0_eff = np.clip(r2 - self._head_grip_s[head_locals], 0.0, None)
                bond_bins = np.clip(
                    np.round(r0_eff / bin_width).astype(np.int64),
                    0, self.p.n_bins - 1,
                )

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
    cortex_bond_groups: np.ndarray | None = None,
    ell0_cortex: float | None = None,
    cortex_beads_per_filament: int | None = None,
    seed_offset: int = 3,
) -> tuple[MyosinStepUpdater, hoomd.update.CustomUpdater]:
    action = MyosinStepUpdater(
        p_myo=p_myo, layout=layout, kT=kT,
        n_cortex_actin=n_cortex_actin,
        cortex_bond_groups=cortex_bond_groups,
        ell0_cortex=ell0_cortex,
        cortex_beads_per_filament=cortex_beads_per_filament,
        seed_offset=seed_offset,
    )
    updater = hoomd.update.CustomUpdater(
        action=action, trigger=hoomd.trigger.Periodic(p_myo.batch_steps)
    )
    return action, updater
