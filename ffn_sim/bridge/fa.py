"""H.4 focal-adhesion topology + HOOMD wiring (KU-2.4 / KU-2.18).

Phase 1 Unit H.4 (PHASE_0_3_DECISIONS D2/D5/D6/D7).

Builds a HOOMD-blue simulation of focal adhesions as **dynamic bond groups**
between ``integrin`` particles (clustered at FA centres) and ``ligand``
particles (a synthetic substrate plane at z=0 for week 1–2 standalone
testing; replaced by H.1 ECM particles in week 3).

Each FA is a HOOMD bond group + Python-side metadata dict (NOT the v1
archived ``acs_kb/bridge/types.py::FocalAdhesion`` dataclass, which lived
on top of the v1 closed-form Chan-Odde force balance and is preserved
read-only at ``ffn_sim/validation/oracles/bridge/types.py``). H.4 v2's
runtime is HOOMD-native:

    H = ½ k_int Σ (|r_integrin − r_ligand| − r0)²    (md.bond.Harmonic)
      + ½ k_xb  Σ (|r_head − r_actin|     − r0_xb)²  (D5 minifilament)
      + Σ LJ_WCA(r_ij; ε, σ, r_cut)                  (D7, repulsive only)

Dynamic binding/unbinding is delegated to
``ffn_sim.bridge.integrin_bonds.IntegrinBondUpdater`` (a hoomd.custom.Action
that mutates the bond topology via snapshot get/set every ``D2 batch_steps``).

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule "Sanity Gate Protocol mandatory before first
execution of any physics/numerics module." STATIC checks live in
``ffn_sim/tests/test_h4_topology.py``; RUNTIME checks are enforced in
``resolve_h4`` / ``build_h4_state``.*

1. **Dimensional analysis**
   - ``k_int`` [N/m] · m² → [J] (bond energy). HOOMD's bond.Harmonic
     ``k`` is [E/length²] = J/m² = N/m. ✓
   - ``γ_int = 6π η R_int`` [Pa·s · m] = [N·s/m]. ✓
   - ``ε_LJ = epsilon_kT · kT`` [J]; ``σ_LJ = R_int + R_lig`` [m];
     ``r_cut = 2^(1/6) σ_LJ`` [m]. ✓
   - ``dt_CFL = α · 1/max(k_on, k_off(F_high))`` [s] (D2 batch CFL).

2. **Boundary cases**
   - ``n_nascent + n_mature == 0``: forbidden (no FAs); RUNTIME raises.
   - ``capture_radius_R_FA ≤ h_integrin_above_substrate``: integrins
     cannot reach any ligand at construction; RUNTIME raises (binding
     would never fire).
   - All FAs placed at the same point: physically degenerate but not
     a numerical failure; warned only.

3. **Conservation invariants**
   - At construction (all integrins disengaged), bond force = 0
     everywhere → ``potential_energy ≈ 0`` to numerical precision.
     STATIC test asserts.
   - Newton 3rd law per bond: HOOMD's ``md.bond.Harmonic`` is symmetric.
   - Particle count: ``n_FAs · n_total_per_fa`` integrins + N_lig
     ligands + N_motors · particles_per_minifilament motor particles.
     STATIC.

4. **Numerical sanity**
   - Bond-group indices in ``[0, N)`` — enforced when building snapshot.
   - Per-particle γ > 0 for every type (integrin, ligand, motor heads).
   - Float64 positions; uint32 bond groups (HOOMD GSD requirement).

5. **Sign / sense**
   - At construction, integrin-ligand bond r0 = 0 with integrin sitting
     ``h_integrin_above_substrate`` above ligand → bond force pulls
     them together (attractive, like a Hookean spring). When the
     IntegrinBondUpdater "binds" a clutch, the harmonic restoring force
     is along the bond axis from integrin to ligand.

6. **Measurement-protocol consistency**
   - The KU-2.4 biphasic / KU-2.5 catch-peak / KU-2.17 FA-growth gates
     all measure per-bond forces via ``md.bond.Harmonic``'s force
     accessor (or its Born-virial extension, mirroring H.1's
     ``_bond_virial_per_bond``). The construction r0 = 0 means a
     "zero-stretch" baseline coincides with the unbound state — the
     gate readout is therefore the *engaged* clutch force only, no
     baseline subtraction needed.

References
----------
- H.4 brief ``ffn_sim/docs/briefs/H4_fa_motor_clutch.md``.
- PHASE_0_3_DECISIONS.md §D2 / §D5 / §D6 / §D7.
- ``ffn_sim/validation/oracles/common/sanity_gate.py``:
  ``gate_unit2_1_motor_clutch``, ``gate_unit2_2_fa_growth``,
  ``gate_minifilament_topology``, ``gate_wca_cutoff``,
  ``gate_bell_evans_batch_cfl``, ``gate_bead_budget_per_cell``.
- v1 reference (acs_kb/bridge/motor_clutch.py): Chan-Odde quasi-static
  force balance + Pereverzev. v2 inverts: HOOMD particle dynamics is
  the runtime, Pereverzev is the validation oracle.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

import gsd.hoomd
import hoomd
import hoomd.md as md

from ffn_sim.integrator.baoab import make_baoab_updater
from ffn_sim.validation.pereverzev import (
    PereverzevParams,
    pereverzev_F_star,
    pereverzev_k_off,
    pereverzev_lifetime,
)


# Mathematical convention (Weeks-Chandler-Andersen): WCA r_cut/σ = 2^(1/6).
WCA_CUTOFF_FACTOR: float = 2.0 ** (1.0 / 6.0)

# Particle type names used in the H.4 frame. Externalised so the
# integrin_bonds Updater + motor module agree on the strings.
TYPE_INTEGRIN: str = "integrin"
TYPE_LIGAND: str = "ligand"
TYPE_MOTOR_BACKBONE: str = "motor_backbone"
TYPE_MOTOR_HEAD: str = "motor_head"
TYPE_ACTIN_BUNDLE: str = "actin_bundle"   # placeholder for week-3 actin coupling

BOND_TYPE_INTEGRIN: str = "integrin_ligand"
BOND_TYPE_MOTOR_HEAD: str = "motor_head_actin"
BOND_TYPE_MOTOR_BACKBONE: str = "motor_backbone"


@dataclass(frozen=True, slots=True)
class ResolvedH4:
    """Frozen H.4 parameter set (filled by ``resolve_h4``).

    All quantities SI; ``derived/*`` come from ``resolve_h4``.
    """

    # Primary scales
    kT: float
    water_viscosity: float
    seed: int
    L_box: float

    # FA topology
    n_nascent_per_cell: int
    n_mature_per_cell: int
    n_total_per_fa: int
    A_nascent: float
    A_mature: float
    h_integrin_above_substrate: float
    capture_radius_R_FA: float

    # Integrin
    k_int_bare: float
    integrin_r0: float
    alpha_vin: float
    integrin_radius: float
    ligand_radius: float

    # Vinculin / talin (toggles)
    vinculin_enabled: bool
    vinculin: dict
    talin_enabled: bool
    talin: dict

    # Pereverzev
    pereverzev: PereverzevParams
    k_on: float

    # D5 motor
    motor: dict
    n_particles_per_minifilament: int
    n_motors_per_fa: int

    # D7 LJ
    lj_enabled: bool
    lj_epsilon: float
    lj_sigma: float
    lj_r_cut: float

    # D3 + D2
    cfl_safety_factor: float
    integrin_batch_steps: int
    bond_event_rate_max: float
    cfl_safe_dt: float
    dt: float

    # FA growth oracle anchors
    fa_growth: dict

    # Derived landmarks
    catch_peak_force: float
    catch_peak_lifetime: float
    gamma_integrin: float
    gamma_ligand: float

    # Acceptance bands (passed through for tests/gates)
    acceptance: dict
    demo_mode: bool


def _stokes_drag(eta: float, R: float) -> float:
    """γ = 6π η R  [N·s/m]  (Stokes drag for a sphere of radius R)."""
    return 6.0 * math.pi * eta * R


def resolve_h4(cfg: Mapping[str, Any]) -> ResolvedH4:
    """Resolve the H.4 yaml config → frozen dataclass with derived fields.

    Mirrors ``ffn_sim.ecm.mikado.resolve_derived`` for H.1: enforces every
    Sanity-Gate-relevant invariant at construction time so downstream
    builders can assume the values are valid SI.
    """
    b = cfg["bridge"] if "bridge" in cfg else cfg

    kT = float(b["kT"])
    eta = float(b["water_viscosity"])
    seed = int(b["seed"])

    fa = b["fa"]
    n_nascent = int(fa["n_nascent_per_cell"])
    n_mature = int(fa["n_mature_per_cell"])
    n_total_per_fa = int(fa["n_total_per_fa"])
    A_nascent = float(fa["A_nascent"])
    A_mature = float(fa["A_mature"])
    h_int = float(fa["h_integrin_above_substrate"])
    R_FA = float(fa["capture_radius_R_FA"])
    L_box = float(fa["L_box"])

    if (n_nascent + n_mature) <= 0:
        raise ValueError(
            "Phase 1 H.4 requires at least one FA per cell; "
            f"got n_nascent={n_nascent}, n_mature={n_mature}."
        )
    if R_FA <= h_int:
        raise ValueError(
            f"capture_radius_R_FA={R_FA:.3e} m ≤ h_integrin_above_"
            f"substrate={h_int:.3e} m; an integrin would never reach a "
            "ligand at z=0. Adjust either parameter."
        )
    if n_total_per_fa <= 0:
        raise ValueError(
            f"n_total_per_fa={n_total_per_fa}; must be ≥ 1."
        )

    integrin = b["integrin"]
    k_int = float(integrin["k_int_bare"])
    if not (np.isfinite(k_int) and k_int > 0.0):
        raise ValueError(f"integrin.k_int_bare must be > 0; got {k_int!r}.")
    integrin_r0 = float(integrin["r0"])

    vinculin = b.get("vinculin", {})
    vinculin_enabled = bool(vinculin.get("enabled", False))
    talin = b.get("talin", {})
    talin_enabled = bool(talin.get("enabled", False))

    pereverzev = PereverzevParams.from_config(cfg)
    k_on = float(b["catch_bond"]["k_on"])

    motor = b["motor"]
    n_backbone = int(motor["n_backbone_beads"])
    n_heads_per_side = int(motor["n_heads_per_side"])
    n_part_mini = n_backbone + 2 * n_heads_per_side
    n_motors_per_fa = int(motor["n_motors_per_fa"])

    # D7 LJ
    ev = b["excluded_volume"]
    lj_enabled = bool(ev["enabled"])
    integrin_radius = float(ev["integrin_radius"])
    ligand_radius = float(ev["ligand_radius"])
    eps_kT = float(ev["epsilon_kT"])
    lj_epsilon = eps_kT * kT
    lj_sigma = integrin_radius + ligand_radius
    lj_r_cut = WCA_CUTOFF_FACTOR * lj_sigma

    # D2 batch CFL: pick dt such that batch_steps · dt · k_off_max ≲ 1e-3
    # with k_off_max = max(k_on, k_off(F_high)).
    # F_high taken at the brief's grid upper bound 60 pN (KU-2.5 sweep).
    F_high = 60.0e-12
    k_off_high = float(pereverzev_k_off(F_high, pereverzev))
    bond_event_rate_max = max(k_on, k_off_high)
    cfl_safety_factor = float(b["dynamics"]["cfl_safety_factor"])
    cfl_safe_dt = cfl_safety_factor / bond_event_rate_max
    # Integrator dt also bounded by the BAOAB CFL on the integrin-ligand
    # spring: τ = γ_int / k_int_eff (eff with vinculin allostery off ⇒
    # k_int_eff = k_int_bare). We pick the smaller of the two bounds.
    gamma_integrin = _stokes_drag(eta, integrin_radius)
    gamma_ligand = _stokes_drag(eta, ligand_radius)
    tau_spring_int = gamma_integrin / k_int
    dt_spring = cfl_safety_factor * tau_spring_int
    dt = min(cfl_safe_dt, dt_spring)
    integrin_batch_steps = int(b["dynamics"]["integrin_batch_steps"])
    # D2 batch CFL gate: events_per_batch = batch_steps · dt · k_off_max
    # must stay below 1e-3 (gate_bell_evans_batch_cfl ceiling).
    events_per_batch = integrin_batch_steps * dt * bond_event_rate_max
    if events_per_batch > 1.0e-3:
        # Shrink batch_steps to satisfy the gate at this dt rather than
        # silently violating it.
        new_batch = max(1, int(1.0e-3 / (dt * bond_event_rate_max)))
        integrin_batch_steps = new_batch

    # Pereverzev analytic landmarks.
    F_star = float(pereverzev_F_star(pereverzev))
    tau_F_star = (
        float(pereverzev_lifetime(F_star, pereverzev))
        if math.isfinite(F_star)
        else float("nan")
    )

    acceptance = b["acceptance"]
    demo_mode = bool(b.get("demo_mode", False))

    return ResolvedH4(
        kT=kT,
        water_viscosity=eta,
        seed=seed,
        L_box=L_box,
        n_nascent_per_cell=n_nascent,
        n_mature_per_cell=n_mature,
        n_total_per_fa=n_total_per_fa,
        A_nascent=A_nascent,
        A_mature=A_mature,
        h_integrin_above_substrate=h_int,
        capture_radius_R_FA=R_FA,
        k_int_bare=k_int,
        integrin_r0=integrin_r0,
        alpha_vin=float(integrin["alpha_vin"]),
        integrin_radius=integrin_radius,
        ligand_radius=ligand_radius,
        vinculin_enabled=vinculin_enabled,
        vinculin=dict(vinculin),
        talin_enabled=talin_enabled,
        talin=dict(talin),
        pereverzev=pereverzev,
        k_on=k_on,
        motor=dict(motor),
        n_particles_per_minifilament=n_part_mini,
        n_motors_per_fa=n_motors_per_fa,
        lj_enabled=lj_enabled,
        lj_epsilon=lj_epsilon,
        lj_sigma=lj_sigma,
        lj_r_cut=lj_r_cut,
        cfl_safety_factor=cfl_safety_factor,
        integrin_batch_steps=integrin_batch_steps,
        bond_event_rate_max=bond_event_rate_max,
        cfl_safe_dt=cfl_safe_dt,
        dt=dt,
        fa_growth=dict(b["fa_growth"]),
        catch_peak_force=F_star,
        catch_peak_lifetime=tau_F_star,
        gamma_integrin=gamma_integrin,
        gamma_ligand=gamma_ligand,
        acceptance=dict(acceptance),
        demo_mode=demo_mode,
    )


# ---------------------------------------------------------------------------
# FA layout helpers
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class FALayout:
    """Per-FA construction metadata persisted into the Simulation.

    Each FA owns ``n_total`` consecutive integrin particle tags
    [tag_start, tag_start + n_total) and a single ligand particle tag
    (under the FA centre at z=0). The IntegrinBondUpdater uses
    ``ligand_tag`` as the binding partner for every integrin in this FA;
    multi-ligand sampling enters at week-3 when H.1 ECM ligands replace
    the synthetic substrate.
    """

    fa_id: int
    is_mature: bool
    centre_xy: np.ndarray            # shape (2,) m
    integrin_tag_start: int
    n_total: int
    ligand_tag: int
    initial_engaged: np.ndarray      # shape (n_total,) bool  (week-1 = all False)


def _build_fa_layout(
    p: ResolvedH4, rng: np.random.Generator
) -> tuple[list[FALayout], np.ndarray, np.ndarray]:
    """Scatter FAs across the box face, allocate integrin + ligand tags.

    Returns
    -------
    layouts : list[FALayout]
        Per-FA metadata, length n_FAs_total.
    integrin_positions : np.ndarray, shape (n_integrin_total, 3)
    ligand_positions : np.ndarray, shape (n_ligand_total, 3)
    """
    n_FAs = p.n_nascent_per_cell + p.n_mature_per_cell
    # Distribute FA centres uniformly across the [-L/2 + R_FA, L/2 - R_FA]²
    # box face. R_FA margin keeps each FA's capture disk inside the box.
    margin = p.capture_radius_R_FA
    half = 0.5 * p.L_box - margin
    centres_xy = rng.uniform(-half, half, size=(n_FAs, 2))

    layouts: list[FALayout] = []
    integrin_positions: list[np.ndarray] = []
    ligand_positions: list[np.ndarray] = []

    n_int_per_fa = p.n_total_per_fa
    next_tag = 0
    # Tag block layout: [integrins of FA0, integrins of FA1, ..., ligands]
    # The integrin updater needs to know the ligand tag separately, so we
    # accumulate integrin tags first, then ligands. (The motor particles
    # come after, see ``_build_motor_particles``.)
    total_integrins = n_FAs * n_int_per_fa
    n_ligands = n_FAs  # one ligand per FA, under the centre

    for i in range(n_FAs):
        is_mature = i >= p.n_nascent_per_cell
        cxy = centres_xy[i]
        # Scatter integrins inside a disk of radius ~ √(A/π) around the
        # centre, at z = h_integrin_above_substrate. Disk radius derived
        # from FA area so mature FAs (A_mature > A_nascent) are
        # geometrically larger.
        A = p.A_mature if is_mature else p.A_nascent
        r_disk = math.sqrt(A / math.pi)
        # Uniform disk sampling: r = R·√u, θ ~ U(0, 2π).
        u = rng.uniform(0.0, 1.0, size=n_int_per_fa)
        theta = rng.uniform(0.0, 2.0 * math.pi, size=n_int_per_fa)
        r = r_disk * np.sqrt(u)
        dx = r * np.cos(theta)
        dy = r * np.sin(theta)
        pos_int = np.zeros((n_int_per_fa, 3), dtype=np.float64)
        pos_int[:, 0] = cxy[0] + dx
        pos_int[:, 1] = cxy[1] + dy
        pos_int[:, 2] = p.h_integrin_above_substrate
        integrin_positions.append(pos_int)

        tag_start = next_tag
        next_tag += n_int_per_fa

        # Ligand sits directly under the FA centre at z=0.
        layouts.append(
            FALayout(
                fa_id=i,
                is_mature=is_mature,
                centre_xy=cxy.copy(),
                integrin_tag_start=tag_start,
                n_total=n_int_per_fa,
                ligand_tag=-1,    # filled in after integrins+ligands stacked
                initial_engaged=np.zeros(n_int_per_fa, dtype=bool),
            )
        )

    # Now allocate ligand tags starting after the last integrin.
    ligand_tag_start = next_tag
    for i, lay in enumerate(layouts):
        lig_pos = np.array(
            [lay.centre_xy[0], lay.centre_xy[1], 0.0], dtype=np.float64
        )
        ligand_positions.append(lig_pos[None, :])
        lay.ligand_tag = ligand_tag_start + i

    integrin_pos = np.concatenate(integrin_positions, axis=0)
    ligand_pos = np.concatenate(ligand_positions, axis=0)
    return layouts, integrin_pos, ligand_pos


def build_h4_state(
    p: ResolvedH4, *, with_motors: bool = False
) -> tuple[gsd.hoomd.Frame, list[FALayout], dict]:
    """Build the construction-time GSD Frame for the H.4 simulation.

    Parameters
    ----------
    p : ResolvedH4
    with_motors : bool, default False
        If True, append D5 Stam-Hocky bipolar minifilaments (one rigid
        backbone + 2 × n_heads_per_side cross-bridge heads per motor).
        Disabled by default in H.4 week 1–2 (clutch-only isolation tests);
        enabled in week 3 when motor stepping enters the dynamics.

    Returns
    -------
    snap : gsd.hoomd.Frame
    layouts : list[FALayout]
        One entry per FA; tag ranges + ligand partner indices.
    meta : dict
        ``{'n_integrin', 'n_ligand', 'n_motor_particles', 'particle_types',
        'bond_types'}`` for downstream gate checks.
    """
    rng = np.random.default_rng(p.seed)
    layouts, integrin_pos, ligand_pos = _build_fa_layout(p, rng)

    particle_types = [TYPE_INTEGRIN, TYPE_LIGAND]
    typeid_int = np.zeros(integrin_pos.shape[0], dtype=np.uint32)        # → "integrin"
    typeid_lig = np.ones(ligand_pos.shape[0], dtype=np.uint32)           # → "ligand"

    if with_motors:
        from ffn_sim.bridge.motor import build_minifilaments_for_fas

        motor_pos, motor_typeids, motor_bonds, motor_bond_types, motor_meta = \
            build_minifilaments_for_fas(p, layouts, rng)
        particle_types += [TYPE_MOTOR_BACKBONE, TYPE_MOTOR_HEAD]
        # Shift motor typeids so they land at indices 2, 3 in particle_types.
        motor_typeids = motor_typeids + 2
    else:
        motor_pos = np.zeros((0, 3), dtype=np.float64)
        motor_typeids = np.zeros(0, dtype=np.uint32)
        motor_bonds = np.zeros((0, 2), dtype=np.uint32)
        motor_bond_types: list[str] = []
        motor_meta = {"n_motor_particles": 0, "n_motor_bonds": 0}

    pos_all = np.concatenate([integrin_pos, ligand_pos, motor_pos], axis=0)
    typeid_all = np.concatenate([typeid_int, typeid_lig, motor_typeids])

    n_part = pos_all.shape[0]

    # Wrap sanity: all positions inside the cubic [-L/2, L/2) box.
    L_half = 0.5 * p.L_box
    if (np.abs(pos_all[:, :2]) > L_half + 1e-9).any():
        bad = np.argwhere(np.abs(pos_all[:, :2]) > L_half + 1e-9)
        raise RuntimeError(
            "H.4 construction placed a particle outside [-L/2, L/2)² on "
            f"the box face; first offending indices: {bad[:5].tolist()}."
        )

    snap = gsd.hoomd.Frame()
    snap.particles.N = n_part
    snap.particles.types = particle_types
    snap.particles.typeid = typeid_all
    snap.particles.position = pos_all
    snap.particles.mass = np.ones(n_part, dtype=np.float64)

    # ---- Bond topology ----
    # Week-1: NO integrin↔ligand bonds at construction. They are added
    # dynamically by IntegrinBondUpdater. We still register the bond
    # type so the bond.Harmonic force can be defined at sim build time.
    bond_groups: list[np.ndarray] = []
    bond_typeids: list[np.ndarray] = []
    bond_types: list[str] = [BOND_TYPE_INTEGRIN]

    # Motor bonds, if any.
    if with_motors and motor_bonds.shape[0] > 0:
        bond_types += motor_bond_types
        bond_groups.append(motor_bonds)
        # Motor bond typeids: 1 for motor_head_actin, 2 for motor_backbone
        # (offset by 1 because integrin_ligand is typeid 0).
        # The motor builder writes typeids relative to its own list; we
        # offset by 1 here so 0→0 from the motor list collides with the
        # integrin_ligand index. ``build_minifilaments_for_fas`` returns
        # typeids in motor_bond_types ordering — we re-index against
        # the merged bond_types here.
        offset_typeids = np.zeros(motor_bonds.shape[0], dtype=np.uint32)
        for raw_id, name in enumerate(motor_bond_types):
            offset_typeids[motor_meta["motor_bond_typeids"] == raw_id] = (
                bond_types.index(name)
            )
        bond_typeids.append(offset_typeids)

    if bond_groups:
        snap.bonds.N = int(sum(g.shape[0] for g in bond_groups))
        snap.bonds.types = bond_types
        snap.bonds.typeid = np.concatenate(bond_typeids).astype(np.uint32)
        snap.bonds.group = np.concatenate(bond_groups, axis=0).astype(np.uint32)
    else:
        snap.bonds.N = 0
        snap.bonds.types = bond_types
        snap.bonds.typeid = np.zeros(0, dtype=np.uint32)
        snap.bonds.group = np.zeros((0, 2), dtype=np.uint32)

    snap.configuration.box = [p.L_box, p.L_box, p.L_box, 0.0, 0.0, 0.0]

    meta = {
        "n_integrin": int(integrin_pos.shape[0]),
        "n_ligand": int(ligand_pos.shape[0]),
        "n_motor_particles": int(motor_pos.shape[0]),
        "particle_types": list(particle_types),
        "bond_types": list(bond_types),
        **motor_meta,
    }
    return snap, layouts, meta


def build_h4_simulation(
    p: ResolvedH4,
    *,
    device: hoomd.device.Device | None = None,
    with_motors: bool = False,
    with_baoab: bool = True,
    with_integrin_updater: bool = True,
) -> tuple[hoomd.Simulation, list[FALayout], dict]:
    """Construct + wire a HOOMD Simulation for H.4.

    Parameters
    ----------
    with_motors : bool, default False
        Append D5 minifilaments. Week-1 disabled.
    with_baoab : bool, default True
        Attach D3 BAOAB-limit Updater.
    with_integrin_updater : bool, default True
        Attach the Pereverzev catch-slip dynamic-bond Updater.

    Returns
    -------
    sim : hoomd.Simulation
    layouts : list[FALayout]
    meta : dict
        ``{'baoab_action', 'baoab_updater', 'integrin_updater',
        'integrin_action', 'snap_meta'}``.
    """
    snap, layouts, snap_meta = build_h4_state(p, with_motors=with_motors)

    sim = hoomd.Simulation(device=device or hoomd.device.CPU(), seed=p.seed)
    sim.create_state_from_snapshot(snap)

    # Bond force: integrin-ligand harmonic. Week-1 has no bonds at
    # construction; the bond.Harmonic force still needs its params
    # registered so the IntegrinBondUpdater can populate bonds dynamically.
    bond = md.bond.Harmonic()
    bond.params[BOND_TYPE_INTEGRIN] = dict(k=p.k_int_bare, r0=p.integrin_r0)
    if with_motors:
        from ffn_sim.bridge.motor import (
            BACKBONE_R0_M, HEAD_REST_LENGTH_M, register_motor_bond_params,
        )

        register_motor_bond_params(bond, p)

    # D7 LJ WCA repulsive. Required pairs: integrin × ligand (always),
    # integrin × motor_head (when motors present). HOOMD 7 requires
    # params for EVERY (typeA, typeB) combination in the state, so we
    # populate zero-strength inert defaults for the unused pairs and
    # then overwrite the active ones.
    #
    # Use md.nlist.Tree (BVH) over md.nlist.Cell — H.4 single-cell box
    # is L_box=20 μm with σ=10 nm, giving L/σ ~ 2000. A Cell list
    # pre-allocates O((L/buffer)³) ≈ 8 × 10⁹ cells and segfaults on
    # allocation (same failure mode as H.1 ECM Mikado per ecm/mikado.py
    # — that file documents the original v2 Cell-list crash).
    nlist = md.nlist.Tree(buffer=0.5 * p.lj_sigma)
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    all_types = list(snap.particles.types)
    for i, t1 in enumerate(all_types):
        for t2 in all_types[i:]:
            lj.params[(t1, t2)] = dict(epsilon=0.0, sigma=p.lj_sigma)
            lj.r_cut[(t1, t2)] = 0.0
    if p.lj_enabled:
        lj.params[(TYPE_INTEGRIN, TYPE_LIGAND)] = dict(
            epsilon=p.lj_epsilon, sigma=p.lj_sigma
        )
        lj.r_cut[(TYPE_INTEGRIN, TYPE_LIGAND)] = p.lj_r_cut
        if with_motors:
            lj.params[(TYPE_INTEGRIN, TYPE_MOTOR_HEAD)] = dict(
                epsilon=p.lj_epsilon, sigma=p.lj_sigma
            )
            lj.r_cut[(TYPE_INTEGRIN, TYPE_MOTOR_HEAD)] = p.lj_r_cut
    lj.mode = "shift"

    ig = md.Integrator(dt=p.dt)
    ig.forces.append(bond)
    ig.forces.append(lj)
    sim.operations.integrator = ig

    out_meta = {"snap_meta": snap_meta, "layouts": layouts}

    if with_baoab:
        # γ for every type present in the snapshot. Integrin / ligand
        # share the same Stokes drag if their radii match; motor
        # particles take a separate γ derived from their bead radius
        # (motor module owns this).
        gamma = {
            TYPE_INTEGRIN: p.gamma_integrin,
            TYPE_LIGAND: p.gamma_ligand,
        }
        if with_motors:
            from ffn_sim.bridge.motor import motor_gammas
            gamma.update(motor_gammas(p))
        action, updater = make_baoab_updater(
            kT=p.kT, gamma=gamma, dt=p.dt, seed=p.seed
        )
        sim.operations.updaters.append(updater)
        out_meta["baoab_action"] = action
        out_meta["baoab_updater"] = updater
    else:
        out_meta["baoab_action"] = None
        out_meta["baoab_updater"] = None

    if with_integrin_updater:
        from ffn_sim.bridge.integrin_bonds import (
            IntegrinBondUpdater,
            make_integrin_updater,
        )
        action, updater = make_integrin_updater(
            sim=sim, p=p, layouts=layouts, bond=bond
        )
        sim.operations.updaters.append(updater)
        out_meta["integrin_action"] = action
        out_meta["integrin_updater"] = updater
    else:
        out_meta["integrin_action"] = None
        out_meta["integrin_updater"] = None

    return sim, layouts, out_meta
