"""D5 Stam-Hocky bipolar myosin minifilament + D6 Hill force-velocity.

Phase 1 H.4 brief §Motor stepping. Shared module — imported by both
``ffn_sim.archive.hoomd_legacy.bridge.fa`` (FA-attached motor population) and (when the H.5
cortex unit lands) ``ffn_sim.archive.hoomd_legacy.cortex.myosin`` (cortical myosin pool). The
common implementation lives here per the H.4 brief deliverable table.

Mechanistic choices (PHASE_0_3_DECISIONS):

- **D5 Stam-Hocky bipolar minifilament**: rigid-rod backbone
  (n_backbone_beads ≈ 14, ~700 nm), n_heads_per_side ≈ 10 cross-bridge
  heads on EACH polarity side (~34 particles per minifilament). NOT the
  AFINES single 2-head spring (which would fail
  ``gate_minifilament_topology.bipolar_two_sides``).
- **D6 Hill force-velocity**: ``v(F) = v0 · (F_s − F) / (F_s + F/a)``,
  Hill 1938 with Kovács 2003 NMII tuning (a/F_s = 0.5). NOT AFINES
  piecewise-linear stall.

H.4 week 1–2 scope ships the **construction + closed-form** layer:

  - ``build_minifilaments_for_fas`` : populate a HOOMD frame with
    bipolar-minifilament particles + (backbone, head-backbone) bond
    topology. Head ↔ actin bonds are dynamic (added at runtime when the
    actin filament couples in week 3) so this builder leaves them out.
  - ``register_motor_bond_params`` : wire bond.Harmonic params for the
    three motor bond types — backbone-backbone, head-backbone, and
    head-actin (the last carries k=k_xb, r0=0 so the runtime stepper
    can append head-actin bonds without re-registering).
  - ``hill_velocity(F)`` : D6 closed-form, used as oracle in
    ``tests/validation/test_ku24_motor_clutch.py`` and (later) by the
    runtime motor-stepping Action.
  - ``motor_gammas(p)`` : per-type Stokes drags for the BAOAB Updater.

The runtime per-head stepping Action (``MyosinStepUpdater``, advances
``pos_a_end`` each step per Hill v(F)) is deferred to H.5 — the H.4
brief explicitly notes the motor module is shared, so we land the
closed-form + topology here and let H.5 add the runtime advance on top.
The KU-2.4 biphasic gate uses the Hill form directly as the oracle, so
the runtime advance is not on the H.4 critical path.

Sanity Gate
-----------
*STATIC checks in ``ffn_sim/tests/test_h4_topology.py``;
``gate_minifilament_topology`` (D5) + ``gate_emergent_vs_oracle`` for
``hill_velocity`` (D6).*

1. Dimensional analysis: backbone bond r0 = backbone_length/(N-1) [m];
   head ⊥ backbone r0 = head_rest_length [m]; head ↔ actin r0 = 0
   (idealised attachment, F = k_xb·|r|). Hill: v0 [m/s], F_s [N],
   a [N], F [N] → v [m/s]. ✓
2. Boundary cases: n_backbone < 2 raises (AFINES 1-bead regression);
   n_heads_per_side ≤ 0 raises (degenerate 1-spring motor); F=0 ⇒ v=v0;
   F=F_s ⇒ v=0; F > F_s ⇒ v < 0 (clamped to 0 by
   ``hill_velocity_clamped``).
3. Conservation: bipolar symmetry — exactly n_heads_per_side on each
   side (``gate_minifilament_topology.bipolar_two_sides``). Bond count
   per minifilament = (n_backbone − 1) backbone + 2·n_heads_per_side
   head-backbone bonds.
4. Numerical sanity: positions float64, bond groups int64 (caller casts
   to uint32 when stacking into the GSD frame). γ_motor_* > 0.
5. Sign / sense: backbone bonds keep the rod straight; head-backbone
   bond restores the head to its perpendicular equilibrium offset. Hill
   v(F): positive F (resisting load) → smaller positive v; stall at
   F = F_s.
6. Measurement-protocol consistency: KU-2.4 biphasic uses the Hill v(F)
   as the actin-bundle velocity in the quasi-static FA force balance;
   the runtime per-head advance arrives at H.5.

References
----------
- Stam et al. 2024 (Stam-Hocky bipolar minifilament).
- Hill AV 1938, *Proc R Soc B* 126:136-95.
- Kovács et al. 2003, *J Biol Chem* 278(40):38132-40 (NMII a/F_s=0.5).
- PHASE_0_3_DECISIONS.md §D5 / §D6.
- ``ffn_sim/docs/briefs/H4_fa_motor_clutch.md`` §Motor stepping.
- ``ffn_sim/validation/oracles/common/sanity_gate.py::
  gate_minifilament_topology``.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

import hoomd.md as md


# Three motor bond types — backbone-backbone (rigid rod), head-backbone
# (perpendicular offset), head-actin (dynamic cross-bridge). Externalised
# so fa.py + tests agree on the strings.
BOND_TYPE_MOTOR_BACKBONE: str = "motor_backbone"
BOND_TYPE_MOTOR_HEAD_BACKBONE: str = "motor_head_backbone"
BOND_TYPE_MOTOR_HEAD_ACTIN: str = "motor_head_actin"

# Geometric reference values; overridden via ResolvedH4.motor.
BACKBONE_R0_M: float = 700.0e-9 / 13.0
HEAD_REST_LENGTH_M: float = 200.0e-9

if TYPE_CHECKING:
    from ffn_sim.archive.hoomd_legacy.bridge.fa import FALayout, ResolvedH4


# ---------------------------------------------------------------------------
# D6 Hill force-velocity (closed-form oracle + runtime kernel reference)
# ---------------------------------------------------------------------------
def hill_velocity(
    F: float | np.ndarray,
    *,
    v0: float,
    F_stall: float,
    a_over_F_stall: float = 0.5,
) -> float | np.ndarray:
    """D6 Hill force-velocity v(F) = v0·(F_s − F)/(F_s + F/a)  [m/s].

    Sign convention: positive F = load resisting head advance (i.e. F
    pulls actin retrograde while the head walks anterograde). Positive v
    = head still advancing under the load. v < 0 corresponds to
    super-stall lengthening; ``hill_velocity_clamped`` zeros that branch.

    Implementation note: Hill's canonical (P + a)(v + b) = (P_0 + a)·b
    has ``a`` as a force-dimension constant. We parametrise with the
    dimensionless ratio ``a_over_F_stall = a / F_stall``, giving the
    equivalent form ``v(F) = v0·(F_s − F)/(F_s + F/a_over_F_stall)``:
    F/a_over_F_stall has units of force (F over dimensionless), so the
    denominator is force-dimensioned and the ratio is dimensionless. The
    1938 NMII reference (Kovács 2003) gives a_over_F_stall ≈ 0.5.
    """
    F_arr = np.asarray(F, dtype=np.float64)
    return v0 * (F_stall - F_arr) / (F_stall + F_arr / a_over_F_stall)


def hill_velocity_clamped(
    F: float | np.ndarray,
    *,
    v0: float,
    F_stall: float,
    a_over_F_stall: float = 0.5,
) -> float | np.ndarray:
    """v(F) with the super-stall lengthening branch clamped to zero."""
    v = hill_velocity(F, v0=v0, F_stall=F_stall, a_over_F_stall=a_over_F_stall)
    return np.clip(v, 0.0, None)


def hill_per_minifilament_velocity(
    F_total: float | np.ndarray, *, p: "ResolvedH4",
) -> float | np.ndarray:
    """Per-minifilament v(F_total) given the per-FA brief defaults.

    ``F_total`` is the load on the whole minifilament; divide by the
    per-side head count (head force-sharing within the rigid-bipolar
    limit) before evaluating Hill.
    """
    motor = p.motor
    return hill_velocity_clamped(
        F_total / max(int(motor["n_heads_per_side"]), 1),
        v0=float(motor["v0_per_head"]),
        F_stall=float(motor["F_stall_per_head"]),
        a_over_F_stall=float(motor["a_over_F_stall"]),
    )


# ---------------------------------------------------------------------------
# BAOAB per-type Stokes drag
# ---------------------------------------------------------------------------
def motor_gammas(p: "ResolvedH4") -> dict[str, float]:
    """Per-type Stokes drags for the BAOAB Updater.

    Phase 1 default: motor backbone + motor head use the same bead radius
    as integrin (≈ 5 nm). Refinement to per-component radii deferred to
    H.5 when the cortex myosin measurements ratify a separate choice.
    """
    R_motor = float(p.integrin_radius)
    eta = float(p.water_viscosity)
    gamma = 6.0 * math.pi * eta * R_motor
    return {"motor_backbone": gamma, "motor_head": gamma}


# ---------------------------------------------------------------------------
# bond.Harmonic parameter wiring
# ---------------------------------------------------------------------------
def register_motor_bond_params(
    bond: md.bond.Harmonic, p: "ResolvedH4"
) -> None:
    """Wire bond.Harmonic params for the three motor bond types.

    Called by ``build_h4_simulation(with_motors=True)`` after registering
    the integrin_ligand bond type.

    - motor_backbone: backbone-backbone, k = 10·head_k (rigid-rod), r0 =
      backbone_length/(N-1).
    - motor_head_backbone: head-backbone perpendicular offset, k =
      head_k, r0 = head_rest_length.
    - motor_head_actin: head-actin cross-bridge, k = head_k, r0 = 0
      (idealised attachment). Pre-registered here so the runtime
      stepper (H.5) can dynamically add head-actin bonds without a
      param re-registration round trip.
    """
    motor = p.motor
    n_backbone = int(motor["n_backbone_beads"])
    backbone_length = float(motor["backbone_length"])
    head_k = float(motor["head_spring_k"])
    head_r0 = float(motor["head_rest_length"])
    backbone_r0 = backbone_length / max(n_backbone - 1, 1)

    bond.params[BOND_TYPE_MOTOR_BACKBONE] = dict(k=head_k * 10.0, r0=backbone_r0)
    bond.params[BOND_TYPE_MOTOR_HEAD_BACKBONE] = dict(k=head_k, r0=head_r0)
    bond.params[BOND_TYPE_MOTOR_HEAD_ACTIN] = dict(k=head_k, r0=0.0)


# ---------------------------------------------------------------------------
# Minifilament topology builder
# ---------------------------------------------------------------------------
def build_minifilaments_for_fas(
    p: "ResolvedH4",
    layouts: "list[FALayout]",
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], dict]:
    """Construct ``n_motors_per_fa`` bipolar minifilaments per FA.

    Layout per minifilament: ``[backbone_0..backbone_{N-1},
    head_+_0..head_+_{H-1}, head_-_0..head_-_{H-1}]``, total
    ``n_backbone + 2·n_heads_per_side`` particles. Backbone beads lie
    along a randomly-oriented in-plane axis at z = 2·h_integrin (above
    the integrin plane). Heads sit perpendicular to the backbone at the
    head_rest_length offset, polarity-separated onto the two sides of
    the rod.

    Returns
    -------
    pos : (n_motor_particles, 3) float64
    typeids : (n_motor_particles,) uint32
        0 = motor_backbone (relative), 1 = motor_head (relative). The
        caller offsets these to the global ``particles.types`` index
        when merging into the H.4 frame.
    bonds : (n_bonds, 2) int64
        Tag pairs, **local to this motor block** (caller adds the
        integrin+ligand tag offset when merging into the global frame).
    bond_type_names : list[str]
        Ordered: ``[BOND_TYPE_MOTOR_BACKBONE, BOND_TYPE_MOTOR_HEAD_BACKBONE]``.
        ``BOND_TYPE_MOTOR_HEAD_ACTIN`` is registered in the bond.Harmonic
        params but no bonds of that type exist at construction (added at
        runtime).
    meta : dict
        ``{'n_motor_particles', 'n_motor_bonds', 'motor_bond_typeids'}``.
        ``motor_bond_typeids[i]`` is the index into ``bond_type_names``
        of the i-th bond.
    """
    motor = p.motor
    n_backbone = int(motor["n_backbone_beads"])
    n_heads_per_side = int(motor["n_heads_per_side"])
    backbone_length = float(motor["backbone_length"])
    head_offset = float(motor["head_rest_length"])
    n_motors_per_fa = int(p.n_motors_per_fa)

    if n_backbone < 2:
        raise ValueError(
            f"motor.n_backbone_beads={n_backbone}; D5 requires ≥ 2 "
            "(1-bead backbone = AFINES regression)."
        )
    if n_heads_per_side <= 0:
        raise ValueError(
            f"motor.n_heads_per_side={n_heads_per_side}; D5 requires "
            "≥ 1 head per side (zero = degenerate single-spring motor)."
        )

    n_part_per_mini = n_backbone + 2 * n_heads_per_side
    z_motor = 2.0 * p.h_integrin_above_substrate
    n_FAs = len(layouts)
    n_total_motors = n_FAs * n_motors_per_fa
    n_total_particles = n_total_motors * n_part_per_mini

    pos = np.zeros((n_total_particles, 3), dtype=np.float64)
    typeids = np.zeros(n_total_particles, dtype=np.uint32)
    bonds_all: list[np.ndarray] = []
    bond_typeids_all: list[np.ndarray] = []

    # bond_type_names ordering: backbone (0), head-backbone (1).
    # head-actin is intentionally absent at construction (registered in
    # bond.Harmonic but no instances; added at runtime when motor steps).
    bond_type_names = [BOND_TYPE_MOTOR_BACKBONE, BOND_TYPE_MOTOR_HEAD_BACKBONE]
    BACKBONE_TID = 0
    HEAD_BACKBONE_TID = 1

    motor_idx_global = 0
    for fa_idx, L in enumerate(layouts):
        cxy = L.centre_xy
        A = p.A_mature if L.is_mature else p.A_nascent
        r_disk = math.sqrt(A / math.pi)
        for _ in range(n_motors_per_fa):
            theta_axis = rng.uniform(0.0, 2.0 * math.pi)
            ux = math.cos(theta_axis)
            uy = math.sin(theta_axis)
            offset_r = r_disk * math.sqrt(rng.uniform(0.0, 1.0))
            offset_theta = rng.uniform(0.0, 2.0 * math.pi)
            cx = cxy[0] + offset_r * math.cos(offset_theta)
            cy = cxy[1] + offset_r * math.sin(offset_theta)

            tag_local0 = motor_idx_global * n_part_per_mini

            # Backbone beads along the in-plane axis.
            s_backbone = np.linspace(
                -0.5 * backbone_length, 0.5 * backbone_length, n_backbone
            )
            for i in range(n_backbone):
                tag = tag_local0 + i
                pos[tag, 0] = cx + s_backbone[i] * ux
                pos[tag, 1] = cy + s_backbone[i] * uy
                pos[tag, 2] = z_motor
                typeids[tag] = 0  # backbone

            perp_x = -uy
            perp_y = ux

            # +polarity heads attached to backbone beads on the first
            # half of the rod, perpendicular offset +head_offset.
            # Distribute heads evenly across the available backbone
            # beads on each polarity side (max one head per bead, mod
            # backbone count if heads > beads per side).
            backbone_per_side = max(n_backbone // 2, 1)
            for h in range(n_heads_per_side):
                tag = tag_local0 + n_backbone + h
                back_idx = (h * backbone_per_side) // max(n_heads_per_side, 1)
                back_idx = min(back_idx, n_backbone - 1)
                bb_x = cx + s_backbone[back_idx] * ux
                bb_y = cy + s_backbone[back_idx] * uy
                pos[tag, 0] = bb_x + head_offset * perp_x
                pos[tag, 1] = bb_y + head_offset * perp_y
                pos[tag, 2] = z_motor
                typeids[tag] = 1  # head
                bonds_all.append(np.array([[tag, tag_local0 + back_idx]]))
                bond_typeids_all.append(np.array([HEAD_BACKBONE_TID], dtype=np.uint32))

            # -polarity heads on the other side of the rod, perpendicular
            # offset −head_offset.
            for h in range(n_heads_per_side):
                tag = tag_local0 + n_backbone + n_heads_per_side + h
                back_idx_from_end = (h * backbone_per_side) // max(n_heads_per_side, 1)
                back_idx = n_backbone - 1 - back_idx_from_end
                back_idx = max(back_idx, 0)
                bb_x = cx + s_backbone[back_idx] * ux
                bb_y = cy + s_backbone[back_idx] * uy
                pos[tag, 0] = bb_x - head_offset * perp_x
                pos[tag, 1] = bb_y - head_offset * perp_y
                pos[tag, 2] = z_motor
                typeids[tag] = 1
                bonds_all.append(np.array([[tag, tag_local0 + back_idx]]))
                bond_typeids_all.append(np.array([HEAD_BACKBONE_TID], dtype=np.uint32))

            # Backbone↔backbone bonds (rigid rod approximation).
            for i in range(n_backbone - 1):
                bonds_all.append(
                    np.array([[tag_local0 + i, tag_local0 + i + 1]])
                )
                bond_typeids_all.append(np.array([BACKBONE_TID], dtype=np.uint32))

            motor_idx_global += 1

    if bonds_all:
        bonds = np.concatenate(bonds_all, axis=0)
        motor_bond_typeids = np.concatenate(bond_typeids_all)
    else:
        bonds = np.zeros((0, 2), dtype=np.int64)
        motor_bond_typeids = np.zeros(0, dtype=np.uint32)

    meta = {
        "n_motor_particles": int(n_total_particles),
        "n_motor_bonds": int(bonds.shape[0]),
        "motor_bond_typeids": motor_bond_typeids,
        "n_motors_total": int(n_total_motors),
        "n_part_per_minifilament": int(n_part_per_mini),
    }
    return pos, typeids, bonds, bond_type_names, meta
