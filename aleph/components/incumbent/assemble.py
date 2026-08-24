r"""Assemble the physiological t0 Active Cell — the ONE composed native mechanical + fluid state (P5).

This is the lead-owned integration point: it builds the resting MCF7 cell by COMPOSING the per-track
force primitives of the new engine into a single device state, so ``driver.py`` can advance them under one
outer-physical / inner-mechanical clock. It edits NO track module and NO ``ff/`` file — it only consumes the
frozen ``accumulate(...)`` primitives read-only and wires the reused ``ff/`` device kernels (bending / link
spring / reshape) that the new engine keeps.

Composed at the physiological setpoint (MCF7 R=7.5 µm, Π₀=40 Pa, the I0-B closed params):
  * cortex mechanical state — ``ac.weave.weave_cell([CORTEX])`` → the unified woven network →
    ``to_crosslinked_cortex()`` → positions + fiber topology + crosslinks (Gate-1: bit-identical to the ff
    cortex). Bending = the reused Cytosim ``cytosim_bending_kernel``; inextensibility = the reused NF2007
    ``reshape_kernel``; crosslinks = the reused ``link_spring_kernel`` (Hookean α-actinin/filamin).
  * MyosinForce (I3) — head-resolved Stam-Hocky minifilaments seeded at the cortex's own myosin sites, heads
    UNBOUND at rest (allocate_hand_state default: bound=0). The active crossbridge is off while unbound, so a
    resting minifilament adds only its passive backbone-rod + head-arm springs (~0 at the rest geometry) — a
    genuinely composed-but-resting motor. walk_dir carries the bipolar-geometry default (the I4-weave attach
    hand-off overwrites it per head on binding; irrelevant while unbound).
  * StericForce (I2b) — all-filament WCA excluded volume over a device HashGrid (SUBSUMES the MT-tip-only
    ``soft_contact_kernel``). Myosin/nucleus/membrane particles are marked dormant (active<0) so EV acts on
    actin nodes only.
  * NucleusCompartment (P4/I2) — the DEFORMABLE-MESH nuclear envelope (``ac.cell.compartments``): reused
    Helfrich bending + framework-#6 lamina areal tension + device-resident nucleoplasm volume + LINC tether
    (nucleus↔cortex, the single load path). Its live oblate mesh is the fluid domain's relative-no-flux inner
    inclusion (``MeshNucleusMaskProvider``), NOT a static sphere. Retires ``nucleus_shell_kernel`` (below).
  * MembraneCompartment (P4) — the independent plasma-membrane Helfrich sheet (``ff.membrane_surface``) at the
    MEASURED cell outer radius R_CELL_UM (the membrane is the outer boundary; PI 2026-07-17 Option B). The
    actin cortex is a shell just INSIDE it at R_CORTEX_UM, so the membrane CONTAINS the cortex (no filament-tip
    leak). Reused Helfrich bending + γ_mem area tension (K_A reservoir OFF, hard-truth #8) + ERM tether — a
    short RADIAL membrane→cortex link ≈ the submembranous gap, force-free at rest.
  * PressureCoupling + the conservative Biot substrate (I1a) — a masked FieldGrid + live
    membrane-minus-nucleus Domain at the resting mean turgor ``p_bar = p = Π₀ = 40 Pa``.  A uniform field has
    zero bulk ``grad(p)`` but still exerts the required ``p n`` boundary traction on the live membrane; the
    latter is assembled explicitly by ``MembranePressureTraction``.  The native NG-3 gate, rather than the
    initializer, decides whether the composed membrane/cortex is actually equilibrated at that preload.

ff/ RETIREMENTS applied BY NON-USE (the driver never launches these lumped kernels):
  I1a  6πηR clock / scalar turgor_kernel   → the scheduler physical clock + Biot p_bar=Π₀
  I2b  soft_contact_kernel (MT-tip only)   → StericForce (all-fiber WCA), k_EV folded into the CFL kmax
  I3   myosin_kernel / f_myo / _myosin_force→ MyosinForce (head-resolved)
  I2   nucleus_shell_kernel (radial bead ball) → NucleusCompartment (deformable mesh, P4). RETIRED SAME-COMMIT
       as this wiring (nucleus INTEGRATION.md §1 same-commit guard: no run ever composes both the bead ball and
       the mesh; the driver launches ONLY the mesh force families, never ``ff.network_warp.nucleus_shell_kernel``).

GAP magnitudes (I0-B): every value whose literature grounding is NOT closed is passed as a clearly-labeled
provisional/TEST value and NEVER tuned to an outcome. At the resting baseline the GAP'd motor magnitudes
(k_xb / f_stall / v0 / κ) do not enter the result (heads unbound), and sigma_EV=7 nm (the physically faithful
actin candidate) leaves EV inactive at the resting shell spacing — both facts are reported, not engineered.

Runtime: Warp-CUDA only (I0-A). ``build_cell`` constructs device state and therefore runs on the gbook A5000;
the dev Mac can only import / syntax-check it (FieldGrid rejects a non-CUDA device).
"""

from __future__ import annotations

import dataclasses
import warnings
from dataclasses import dataclass, field
from types import SimpleNamespace

import numpy as np
import warp as wp

# ── ac/ track primitives (consumed read-only via their §1.4 accumulate contracts) ────────────────────
from aleph.components.incumbent.compartments import (
    MembraneCompartment,
    NucleusCompartment,
    build_membrane_compartment,
    build_nucleus_compartment,
)
from aleph.components.incumbent.erm_tether import ERMBellKinetics
from aleph.components.incumbent.live_mesh_domain import LiveMeshMembraneMaskProvider
from aleph.components.incumbent.membrane_pressure import MembranePressureTraction, signed_volume
from aleph.components.incumbent.preload_contract import evaluate_preload_capacity, triangle_surface_area
from aleph.components.fluid.biot_substrate import BiotSubstrate, PressureCoupling
from aleph.components.fluid.boundary import MembraneFluxBC
from aleph.components.fluid.domain import Domain, StaticSphereMembraneProvider
from aleph.components.fluid.field_grid import FieldGrid
from aleph.components.motor.backbone_warp import (
    arm_orientation_k_theta,
    backbone_bending_k_theta,
    backbone_bending_kappa,
)
from aleph.components.motor.hand import NMIIHandParams
from aleph.components.motor.minifilament_topology import MinifilamentTopology, straddle_frame
from aleph.components.motor.minifilament_warp import MyosinForce, build_minifilament_nodes
from aleph.components.motor.resting_setpoint import (
    RestingBoundMyosinSetpoint,
    apply_resting_bound_heads,
    plan_resting_bound_heads,
    resting_tangential_load_reference,
)
from aleph.components.solid.steric_warp import StericForce
from aleph.components.weave.regions import CORTEX_REGION, RegionSpec, cortex_arp23_region, derive_cortex_arp23_split
from aleph.components.weave.walk_dir import barbed_end_node
from aleph.components.weave.woven_cell import actin_segment_topology, weave_cell

# ── engine-agnostic parameter gates (instrumentation; no new biology — ac/cell is feature-frozen) ────
from aleph.laws.turgor_pi0 import (
    PI0_CLAIM_HELA_PROXY_40PA,
    convenience_labels_record,
    resolve_pi0_pa,
)

# ── reused ff/ device kernels the new engine KEEPS (bending / link spring / reshape / step) ──────────
from aleph.laws.forces_warp import _per_triple_alpha

__all__ = ["CellConfig", "AssembledCell", "build_cell"]

# ── physiological setpoints ──────────────────────────────────────────────────────────────────────────
# R_CELL_UM = 7.5 µm is the MCF-7 SUSPENDED single-cell radius anchored to Wagner et al. 2011 (Free Radic Biol
# Med 51(3):700, PMC3147247; Coulter Counter, MEASURED dia 14.8 µm / vol 1.70 pL → R ≈ 7.4 µm). This anchor is
# used consistently across dcm/geometry.py, dcm/nondim.py, the test suite, and docs/LAYER2_ANCHORS /
# MCF7_CONVERSION_SPEC. Prior versions of THIS comment warned "do NOT assert Wagner 2011 (absent from the KB)":
# that was a KB-REGISTRY gap, not a fabrication — the paper is real and standard; it was flagged un-registered
# in SE_REGISTRATION_CANDIDATES_2026-06-30.md:131. PI 2026-07-21: register Wagner2011 as SourceEvidence +
# a suspended-radius KnowledgeClaim (reconciliation in docs/v2_audit/CORTEX_COUNT_CROSSCHECK_2026-07-21.md).
# The TAG/KB registry previously held only PARAM-cell_radius = 10 µm (a Phase-1 2D-disk placeholder, NO
# citation), now superseded by the 7.5 µm Wagner anchor. NOTE the 70,686 count itself (= 100 µm⁻² × 4π(7.5 µm)²)
# uses a GENERIC (non-MCF7) cortex areal density. Its MAGNITUDE is order-correct — it lands in the native
# 30–50 nm cortex-mesh regime (Bovellan2014); nobody counts cortical filaments directly (too dense). The real
# gap is filament LENGTH: native cortical filaments are ~hundreds of nm, so these 70,686 × 3 µm are COARSE
# representative filaments (Cytosim-style), not a short-filament census. See the CORTEX_COUNT_CROSSCHECK doc.
R_CELL_UM = 7.5           # I0-A-ratified cell OUTER radius = the MEMBRANE (PI 2026-07-17: "measurement basis
#                           was the membrane"). Membrane + fluid osmotic envelope at R_CELL_UM; the actin
#                           cortex is a shell just INSIDE it, so the membrane CONTAINS the cortex.
# The actin cortex is a ~h_cortex-thick layer just beneath the membrane. h_cortex ≈ 200 nm is KB-SOURCED and
# verified (KB-3.1 / KB-3.5, Salbreux2012_TCB · Charras2008_BJ · Clark2013 DOI 10.1016/j.bpj.2013.05.057). The
# single representative cortex shell is placed at ≈ ½·h_cortex below the membrane (the cortical-layer mid-depth;
# the finer submembranous-gap-vs-half-thickness split is not separately KB-resolved and is immaterial at the
# resting baseline — ERM is force-free, a ~1% cortex-radius change is mechanically negligible). This clears the
# cortex radial fluctuation band (±~12 nm) so the membrane CONTAINS every cortex node.
CORTEX_MEMBRANE_GAP_UM = 0.10      # ≈ ½·h_cortex (h_cortex≈200 nm, KB-3.1/3.5) — membrane→cortex-shell offset [µm]
R_CORTEX_UM = R_CELL_UM - CORTEX_MEMBRANE_GAP_UM   # 7.40 µm — the cortex shell sits ~½ a cortex-thickness inside
# ⚠ Π₀ IS GATED (PI 2026-07-25 — aleph/laws/params_turgor.yaml holds `value: null`,
# `evidence_status: "GAP — PI"`). `PI_0_PA` below is NOT a silent default: it is the EXPLICIT,
# provenance-classified 40 Pa CONVENIENCE claim (Fischer-Friedrich 2014 HeLa interphase — measured,
# WRONG CELL LINE), resolved through the gate at import so the value can never be handed out without
# its provenance. `PI_0_PROVENANCE` is emitted into every build_cell artifact.
# THREE values are live in the tree and NONE is an MCF7 measurement: 40 Pa (here, HeLa proxy),
# 133 Pa (configs/mcf7_baseline.yaml:81 + common/production_policy.py, band-implied ⇒ CIRCULAR),
# ~72 Pa (48–107, MCF7-geometry estimate, unwired). Π₀ sets the ENTIRE resting cortical tension via
# γ = ΔP·R/2, so this ~3.3× spread is a ~3.3× uncertainty on the resting baseline itself
# (PARAM_PROVENANCE_AUDIT_2026-07-24.md row 1: "the worst convenience default").
_PI_0 = resolve_pi0_pa(PI0_CLAIM_HELA_PROXY_40PA, claim="claim_hela_proxy_40pa",
                       caller="ac.cell.assemble")
PI_0_PA = _PI_0.value_pa   # resting osmotic turgor Π₀ [Pa]; 1 Pa == 1 pN/µm² so 40 in engine pressure units
#: Serializable Π₀ provenance — recorded in the build_cell ledger so no artifact carries the number
#: without "CONVENIENCE / not sourced for MCF7" attached to it.
PI_0_PROVENANCE = _PI_0.as_artifact_record()
# I0-B1 fluid closure (params_i0b1.yaml, engine units): c_v = mobility/storage_S = 50 µm²/s.
BIOT_STORAGE_S = 1.0e-4   # S = 1/M  [µm²/pN]  (M≈1e4 Pa; derived from the c_v anchor)
BIOT_MOBILITY = 5.0e-3    # k/µ  (so mobility/S = 50 µm²/s = c_v)  [µm⁴/(pN·s)]
BIOT_ALPHA = 1.0          # Biot-Willis coupling (PI-ratified 2026-07-16)
L_P = 1.6e-8              # membrane hydraulic conductivity [µm/(s·Pa)] (gamma_floor code value; draft/MCF7)

# ── I0-B2b steric (params_i0b2b.yaml) — sigma_EV GAP; use the physically faithful actin candidate ────
SIGMA_EV_UM = 0.007       # F-actin steric diameter (7 nm) — the physical candidate (PI discretization GAP)
K_EV_PROVISIONAL = 1.0e3  # WCA contact stiffness [pN/µm] — PROVISIONAL Magic-Number-Block scale (close at
#                           the native gate from the measured passive F_op; NOT tuned to a crowding outcome)

# ── I0-B3 NMII minifilament reference topology (minifilament_topology docstring; magnitudes are GAP) ──
NMII_N_BB = 14            # backbone beads (reference; count-invariant gate holds for any n_bb≥2)
NMII_N_SIDE = 10          # heads per anti-parallel side (AFINES claim_a; N_side is a PI GAP)
NMII_L_BB_UM = 0.301      # backbone contour (Billington 2013 EM; L_bb is a PI GAP)
NMII_HEAD_OFFSET_UM = 0.200   # head↔backbone arm offset (= r0_head; PI GAP, arm rests force-free at build)
# provisional/TEST mechanical magnitudes (do NOT enter the resting result — heads are UNBOUND at t0):
NMII_K_XB_TEST = 1.0e3        # crossbridge stiffness [pN/µm], mid physical band 100–1000 (MASTER knob, GAP)
NMII_K_HEAD_ARM_TEST = 1.0e2  # head↔backbone arm stiffness [pN/µm] (provisional; arm at rest ⇒ ~0 force)
NMII_K_BACKBONE_TEST = 1.0e3  # rigid-rod backbone stiffness = 10× arm (k_backbone_factor=10, numerical)
NMII_F_STALL_TEST = 0.5       # per-head stall [pN] (AFINES claim_a; GAP) — unused while unbound
NMII_V0_TEST = 0.12           # unloaded velocity [µm/s] (GAP) — unused while unbound
NMII_KAPPA_TEST = 0.5         # Hill curvature (Kovács archived; GAP) — unused while unbound
NMII_KON_TEST = 50.0          # attach rate [1/s] (config value; GAP) — unused while unbound
NMII_KOFF0 = 0.35             # Bell prefactor [1/s] (draft) — unused while unbound
NMII_F0 = 4.28e-3 / 0.6e-3    # Bell force f0 = kBT/x_beta ≈ 7.13 pN (physical)
NMII_CAPTURE_UM = 0.210       # head→actin perpendicular capture reach [µm] = r0_head (0.200) + ~10 nm slack.
# RECONCILED 2026-07-24 (PARAM_INCONSISTENCY_RECONCILE): this is the SAME physical quantity as ff/hand_kmc.py
# NMIIA_MYOSIN.capture_radius_um = 0.210 ("head_actin_capture_perp", configs/phase1_h3.yaml:386 = head_rest_length
# + 10 nm slack, KU-3.5 binding fix 2026-05-29, docs/KU35_myosin_binding_diagnosis.md). The prior 0.05 was the
# EXPLICITLY-RETIRED legacy value — phase1_h3.yaml:382-384: "Legacy 50nm was for the wrong radial-offset
# placement that left heads ~824nm from any actin." It also contradicted this module's own NMII_HEAD_OFFSET_UM =
# 0.200. Gates a point-to-segment perpendicular distance (segment_query); heads are UNBOUND at the resting
# baseline so this does not enter GATE-A. NOT the loose 0.6 diagnostic TEST reach some seed-probe scripts pass
# via resting_bound_myosin_capture_um (that is a deliberate binding-forcing knob, not this physiological reach).
# Head-on-actin STRADDLE placement (defect#3): pick an anti-parallel partner filament ~2·head_offset away so
# the ± heads land ON the two filaments (perpendicular residual → 0 at exactly 2·offset). Geometric selection
# windows (build-time; not physics magnitudes tuned to an outcome).
NMII_STRADDLE_SEP_TOL = 0.25                  # accept |partner−anchor| within ±25% of 2·head_offset
NMII_STRADDLE_ANTIPARALLEL_MAX_DOT = -0.2     # require t_A·t_B < this (the two filaments run anti-parallel)
# ── F6 angle-harmonic bending (Round-2 motor-fix): the backbone rigid-rod + head-arm orientation stiffnesses
# that make the crossbridge tension TRANSMIT collinearly (NG-1). DERIVED, not tuned (backbone_warp
# Magic-Number-Block): k_θ,bb = κ/a with κ = L_p·kBT (rigid-rod limit, L_p GAP→PI, convergence-checked);
# k_θ,arm = k_xb·r0_head² (rigid-lever limit). Unused while heads are UNBOUND at t0 (bending is internal,
# ~0 at the straight/perpendicular build rest). ─────────────────────────────────────────────────────────
NMII_BACKBONE_LP_DIAGNOSTIC_UM = 1.0  # convergence-sweep fixture only; never a production/runtime default


@dataclass
class CellConfig:
    """Assembly configuration (population + numerics). Physiological magnitudes are module constants."""

    n_filaments: int = 70686        # FULL native cortical F-actin population (100 µm⁻²·4π·7.5²)
    # ── CONFIGURABLE CORTEX GEOMETRY (mesh-fidelity study, 2026-07-24; CORTEX_MESH_FIDELITY memo §7) ──
    # These surface the previously HARD-CODED cortex geometry (architecture_spec.CORTEX) so a physiologically
    # finer-mesh cortex can be built and compared to the coarse baseline. DEFAULTS REPRODUCE TODAY EXACTLY
    # (bit-parity): 0.5/3.0/20.0 are the current architecture_spec.CORTEX values — additive, not a new default.
    # The mesh RESOLUTION is set by cortex_seg_um (crosslinks land on nodes seg_um apart → pore size ≈ seg_um);
    # the crosslink SPACING along the contour is cortex_length_um / cortex_density_per_fil. Physiology
    # (Morone 2006; Bovellan 2014; Chugh & Paluch 2018): cortical mesh ≈ crosslink spacing ≈ 50–100 nm, so a
    # ~75 nm mesh needs cortex_seg_um ≈ 0.075 µm and cortex_density_per_fil ≈ length/0.075 ≈ 40 (both DERIVED
    # from the sourced mesh, not tuned). See the memo §7 derived fine config + its node/byte cost.
    # ⚠⚠ ALL THREE DEFAULTS BELOW ARE **CONVENIENCE**, NOT SOURCED (PARAM_PROVENANCE_AUDIT_2026-07-24.md
    # rows 2/3/4). PI 2026-07-25: LABEL them, do NOT gate them, do NOT change the values. Machine-readable
    # record: aleph/laws/params_turgor.yaml `convenience_labels:` (emitted into every build_cell
    # artifact via `convenience_labels_record()`), so a run's own ledger states the deviation.
    #   • cortex_seg_um  = 0.5 µm  → CONVENIENCE (coarse). Sourced cortical mesh / crosslink spacing is
    #     50–100 nm (Morone 2006 · Bovellan 2014 · Chugh & Paluch 2018) ⇒ **5–10× TOO COARSE**. 0 KB rows
    #     support 0.5 µm. This is the deeper cause of the mesh-fidelity gap.
    #   • cortex_length_um = 3.0 µm → CONVENIENCE (coarse). Sourced ~1 µm formin / ~0.1 µm Arp2/3
    #     (Bovellan 2014, KB-3.18) ⇒ **3–30× TOO LONG**; with nucleator="formin" only, the Arp2/3
    #     short-filament population (~⅓ of cortical actin by mass) is **entirely ABSENT** — a missing
    #     ARCHITECTURE, not just a wrong number.
    #   • cortex_density_per_fil = 20.0 → CONVENIENCE (percolation). No sourced value. ⚠ two DIFFERENT
    #     justifications are on the record for the same 20.0: "DERIVED L/δ_filamin = 3.0 µm/150 nm" and
    #     "chosen so the native network stays single-spanning". It is also DOWNSTREAM of cortex_seg_um,
    #     so it cannot be closed before the mesh is sourced.
    cortex_seg_um: float = 0.5          # CONVENIENCE (coarse): 5–10× above the sourced 50–100 nm mesh
    cortex_length_um: float = 3.0       # CONVENIENCE (coarse): 3–30× above sourced formin/Arp2/3 lengths
    cortex_density_per_fil: float = 20.0  # CONVENIENCE (percolation): no sourced value; downstream of seg_um
    # ── MIXED formin + Arp2/3 CORTEX (FLAGGED prototype, 2026-07-24; CORTEX_ARP23_POPULATION memo) ──
    # The default cortex is formin-only; native cortex is ~2/3 formin (long) + ~1/3 Arp2/3-nucleated (short,
    # branched) BY MASS (Bovellan 2014, KB-3.18). DEFAULT-OFF (fraction 0.0) ⇒ today's pure-formin cortex,
    # BIT-IDENTICAL. When >0, the fixed n_filaments cortex BUDGET is split into (1−)·formin + Arp2/3 branched
    # so total areal density (~100/µm²) is PRESERVED, with the count split DERIVED from the mass fraction
    # (mass ∝ n·length; ac.weave.regions.derive_cortex_arp23_split). The Arp2/3 length/seg + mother fraction are
    # FLAGGED modeling choices (a "PI A/B/C" decision — a sourced OPTION, not a silent default); see the memo.
    cortex_arp23_fraction: float = 0.0      # Arp2/3 fraction of cortical actin BY MASS (~0.33 Bovellan; 0 = OFF)
    cortex_arp23_length_um: float = 0.15    # representative SHORT Arp2/3 filament length [µm] (0.1–0.2; PI A/B/C)
    cortex_arp23_seg_um: float = 0.05       # Arp2/3 segment rest length [µm] (≥3 nodes/filament ⇒ branchable)
    cortex_arp23_mother_fraction: float = 0.2  # fraction of the Arp2/3 pop that are branch ROOTS (rest = daughters)
    with_myosin: bool = True        # compose the head-resolved MyosinForce (resting, heads unbound)
    with_steric: bool = True        # compose the all-fiber WCA StericForce
    overlap_free_cortex: bool = True    # build the cortex overlap-free (radial ~0.2µm thickness + WCA relaxation)
                                    # so the shell starts with ZERO steric force — removes the ~64k build-
                                    # interpenetration artifact (steric-only max 2634 pN) that dominates the
                                    # resting-convergence residual, and gives a real 3D shell (R∈[7.20,7.40],
                                    # std≈0.06µm) instead of a zero-thickness sphere. Default ON = the
                                    # PHYSIOLOGICAL BASELINE (h_cortex~0.2µm, KB-3.1/3.5; density_per_fil=20 keeps
                                    # the network single-spanning under the −44% near-pair loss, 2.1× headroom).
    cortex_overlap_mode: str = "transverse"   # HOW overlap_free_cortex resolves build interpenetrations. Both
                                    # settle to the SAME param-free WCA target (shell overlap-free either way).
                                    # "transverse" (DEFAULT, byte-identical to history): shove crossing nodes
                                    # apart IN the shell tangent plane — but that knocks a node off its smooth
                                    # great-circle arc = a per-node bending KINK whose cytosim force F=α·d
                                    # (α=κ/seg³) blows up ×288 at the fine 75nm mesh (fine-cortex GATE-A plateau).
                                    # "radial_span": separate crossing fibers OUT-OF-PLANE (sphere-radial) with
                                    # the displacement spread over a ±cortex_overlap_span-node smooth cosine bump
                                    # along each fiber, so both fibers stay low-curvature arcs while separating in
                                    # 3D. Opt-in smoothness-preserving fix for the fine mesh; leaves the coarse
                                    # (0.5µm) build unaffected in intent (kink invisible there).
    cortex_overlap_span: int = 2    # radial_span half-window: bump spans ±this many arc-neighbours (0=single node)
                                    # Pass overlap_free_cortex=False to reproduce the legacy zero-thickness
                                    # γ-validation / Gate-1 RNG-parity build.
    with_nucleus: bool = True       # compose the deformable-mesh NucleusCompartment (P4/I2) + no-flux mask
    with_membrane: bool = True      # compose the plasma-membrane Helfrich sheet (P4)
    membrane_subdivisions: int = 3  # icosphere level for the plasma membrane. Default 3=642 verts (~1µm) is the
                                    # smooth RESTING sphere (backward-compat / γ-validation). PI-ratified
                                    # 2026-07-22 for DYNAMIC morphology: **6 = 40,962 verts (~131nm) is the
                                    # production baseline** (blebs ~0.5µm → ~4 nodes/neck; ERM 58/µm²), **7 =
                                    # 163,842 (~66nm) is the validation resolution** (~8 nodes/neck; ERM
                                    # 235/µm² physiological). Higher densifies the per-node ERM coupling
                                    # (n_erm=n_verts). Actual dynamic runs may instead use adaptive bleb-site
                                    # remeshing (the fidelity endgame, deferred). Memory is a non-issue (subdiv 6
                                    # full-native = 592 MB on the A5000 16 GB).
                                    # PI-ratified 2026-07-24 for the RESTING baseline: **8 = 655,362 verts (~33nm)
                                    # is the production RESTING config**. GATE A block-split localized the clean
                                    # 0.776 residual to the MEMBRANE; a subdiv 6→8 sweep drops it 10× (0.7766 →
                                    # 0.0768, 0 of 655,362 over the 0.21 gate) — textbook GRID CONVERGENCE of the
                                    # Helfrich/Young-Laplace residual (edge→0), NOT gate-loosening. At subdiv 8 the
                                    # CLEAN resting baseline PASSES GATE A (cortex 0.036 / membrane 0.077 / nucleus
                                    # 0.16, all < 0.21). A5000 builds it fine (~+614k nodes). See
                                    # docs/v2_audit/GATE_A_RESTING_CONVERGENCE_2026-07-23.md §2/§4.4.
    nucleus_subdivisions: int = 3   # icosphere level for the nuclear envelope (same 642→40k→163k scaling; PI
                                    # 2026-07-22: 6 dynamic baseline / 7 validation, as for the membrane).
    with_pressure: bool = True      # compose the Biot substrate + PressureCoupling
    #: Π₀ as a DECLARED AXIS rather than a constant (PI standing ruling 2026-08-15: every physiological
    #: value is physiological for some particular cell and is not certain, so it is an axis with a band).
    #: ``None`` takes the gated ``PI_0_PA`` and is byte-identical to the pre-axis behaviour — the default
    #: is not a second value, it is the same one. A non-null value REPLACES the turgor everywhere it acts
    #: (grid seed, CFL term, preload, membrane osmotic balance) so the arms cannot drift apart.
    #: ``0.0`` is the load-free axis point that separates the turgor LOAD from the Biot COUPLING, which
    #: ``with_pressure=False`` cannot do because it removes both at once.
    turgor_pi0_pa: float | None = None
    # No MCF7 ERM areal-density datum is registered. ``None`` keeps a diagnostic one-per-membrane-node
    # topology but NG-3 must remain open. A non-null production value requires an explicit provenance label.
    erm_density_per_um2: float | None = None
    erm_density_source: str = ""
    erm_density_mcf7_production: bool = False  # explicit provenance latch; cross-cell proxies must keep False
    erm_radial_pairing: bool = False  # pair each membrane node to its most-RADIAL cortex node so the tether
                                      # transmits the radial turgor load between the shells (resting baseline)
    # Bell on/off kinetics are an all-or-none, source-gated contract.  They intentionally have no biological
    # defaults because MCF7 k_on/k_off0/F0/capture are absent from the Contract-Graph (2026-07-21 audit).
    erm_k_on_s: float | None = None
    erm_k_off0_s: float | None = None
    erm_bell_force_pn: float | None = None
    erm_capture_radius_um: float | None = None
    erm_kinetics_source: str = ""
    erm_rebind_rest_policy: str = "formation_length"
    # The backbone persistence length is a physical magnitude, not a numerical rigid-rod convenience.  The
    # Contract-Graph currently has no ratified NMII-minifilament value, so production-form bending is opt-in
    # and source-gated; ``None`` leaves the gap explicit instead of running the diagnostic sweep value.
    nmii_backbone_lp_um: float | None = None
    nmii_backbone_lp_source: str = ""
    # NMII minifilament head-on-actin PLACEMENT (defect#3 fidelity, 2026-07-23). Each bipolar minifilament is
    # oriented to STRADDLE two anti-parallel actin filaments ~2·head_offset apart, so its ± heads land on actin
    # (~nm crossbridge extension) instead of the fixed ±head_offset perpendicular projection into a mesh void
    # (which left the heads 0.13–0.6µm off actin — only ~17% seedable at native). ON = the PHYSIOLOGICAL
    # BASELINE (heads touch the actin they will bind; the bipolar dipole anti-parallel-anchors, so it sustains
    # a balanced contractile tension). Per-minifilament fallback to the legacy midpoint placement whenever no
    # anti-parallel partner is found within tolerance. Pass False to reproduce the legacy myo geometry.
    nmii_straddle_placement: bool = True
    # RESTING BOUND-MYOSIN SETPOINT (candidate 1 — the fine-grained-faithful resting cortical-tension source;
    # see docs/v2_audit/RESTING_BASELINE_DIAGNOSIS_2026-07-23{c,d}.md and ac.motor.resting_setpoint). Additive,
    # DEFAULT-OFF for parity: with both fields None the motor stays unbound at t0 (bit-identical build) and the
    # PI-GAP is surfaced. BOTH magnitudes are PI-GAPs absent from the Contract-Graph — never defaulted here; the
    # fraction × per-head force sets the resting cortical tension γ_cortex = ΔP·R/2 − γ_mem the ERM transmits.
    # When ON, the resting path also engages RADIAL ERM pairing (so the tether transmits the radial turgor load).
    resting_bound_myosin_fraction: float | None = None
    resting_bound_myosin_force_pn: float | None = None
    resting_bound_myosin_source: str = ""
    resting_bound_myosin_capture_um: float | None = None  # bind reach [µm]; None ⇒ the motor's params.capture_radius
    dx_um: float = 0.5              # conservative field-grid resolution (numerical sizing; not literature)
    grid_pad_um: float = 3.0        # field grid half-extent beyond R_cell
    steric_grid_dim: int = 256      # hash-grid hash dimension per axis (accelerator only, grid-invariant)
    motor_segment_grid_dim: int = 128  # segment-midpoint HashGrid dimension (accelerator only)
    steric_force_cap: float | None = 1.0e3  # WCA near-core CFL guard [pN] (params_i0b2b force_cap_EV): bounds
    #   the divergent r→0 core at the handful of build-time cross-filament coincidences the physical σ_EV=7 nm
    #   node-node EV finds on the dense ff cortex (the σ_EV discretization GAP). A declared numerical guard
    #   (>> the ~0.02 pN physical operating load; never a physics magnitude tuned to an outcome), NOT the
    #   steric verdict — the load-bearing no-interpenetration verdict defers to the σ_EV/segment-segment
    #   resolution (PI GAP) + the post-I3 contractile re-run.
    device: str | None = None          # resolve the current/default CUDA device; never hard-code an ordinal
    seed: int = 0


@dataclass
class AssembledCell:
    """The composed device cell state consumed by :mod:`ac.cell.driver`."""

    cfg: CellConfig
    device: str
    # combined mechanical node array (actin nodes [0:n_actin) + myosin particles [n_actin:n_total))
    n_actin: int
    n_total: int
    pos_d: wp.array
    f_d: wp.array
    state: SimpleNamespace               # duck-typed (.pos, .node_pos, .node_volume) for the primitives
    # bending (actin)
    tri_d: wp.array
    alpha_d: wp.array
    n_tri: int
    # crosslinks (link_spring)
    xl_d: wp.array
    kxl_d: wp.array
    r0xl_d: wp.array
    n_xl: int
    # reshape (actin fibers only)
    foff_d: wp.array
    toff_d: wp.array
    soff_d: wp.array
    srest_d: wp.array
    n_fibers: int
    max_fiber_nodes: int
    # composed primitives
    steric: StericForce | None
    pressure: PressureCoupling | None
    myosin: MyosinForce | None
    nucleus: NucleusCompartment | None
    membrane: MembraneCompartment | None
    membrane_pressure: MembranePressureTraction | None  # spatial p -> live membrane normal traction (NG-3)
    nucleus_mask_provider: object | None   # MeshNucleusMaskProvider — the live no-flux inclusion (or None)
    membrane_mask_provider: object | None  # LiveMeshMembraneMaskProvider — the live outer boundary (or None)
    node_volume_d: wp.array
    solid_active_d: wp.array             # (n_total,) int32; actin skeleton >=0, non-skeleton <0
    # I4 walk_dir hand-off arrays (bound heads walk the ACTUAL actin polarity, not the nearest node)
    node_fiber_d: wp.array | None       # (n_actin,) int32 — fiber id of every actin node
    barbed_node_d: wp.array | None      # (n_fibers,) int32 — global barbed-end node of each fiber
    # fluid substrate + scheduler pieces
    grid: FieldGrid | None
    domain: Domain | None
    substrate: BiotSubstrate | None
    membrane_bc: MembraneFluxBC | None
    # numerics
    kmax: float
    dt_mu: float
    mechanical_kmax: float         # all non-pressure stiffnesses; live |Delta p| term is device-composed
    pressure_edge_um: float        # live membrane mean edge in d(pA)/dx ~ |Delta p| ell
    convergence_length_um: float       # actin segment discretization length for sqrt(eps64)*ell tolerance
    # I4 Arp2/3 70° branch-angle junctions (mixed formin+Arp2/3 cortex only; None for the formin-only default,
    # so no branch kernel launches and the build stays bit-identical). Same angle-harmonic the lamellipodium uses.
    branch_triples_d: wp.array | None = None    # (n_branch, 3) int32 — [m_after, branch_vertex, daughter_arm]
    branch_active_d: wp.array | None = None      # (n_branch,) int32 — 1 = live junction, 0 = dormant (no force)
    n_branch: int = 0
    # bookkeeping
    ledger: dict = field(default_factory=dict)

    def com(self) -> np.ndarray:
        """Centre of mass of the ACTIN nodes (host read; out-of-hot-loop diagnostic only)."""
        return self.pos_d.numpy()[: self.n_actin].mean(axis=0)


def _cortex_region(
    n_filaments: int,
    *,
    seg_um: float = 0.5,            # CONVENIENCE (coarse) — see CellConfig.cortex_seg_um + params_turgor.yaml
    length_um: float = 3.0,         # CONVENIENCE (coarse) — see CellConfig.cortex_length_um
    density_per_fil: float = 20.0,  # CONVENIENCE (percolation) — see CellConfig.cortex_density_per_fil
) -> RegionSpec:
    """A CORTEX region spec at a chosen filament count + GEOMETRY, with the cortex shell placed just INSIDE the
    plasma membrane at ``R_CORTEX_UM`` (= R_CELL_UM − CORTEX_MEMBRANE_GAP_UM). The membrane (the measured MCF7
    outer radius, PI 2026-07-17 Option B) therefore CONTAINS the cortex, and the ERM tether becomes a short
    RADIAL membrane→cortex link ≈ the submembranous gap (not a lateral tether at a coincident radius).

    ``seg_um`` / ``length_um`` / ``density_per_fil`` surface the previously hard-coded cortex geometry
    (``architecture_spec.CORTEX``) so a physiologically finer-mesh cortex can be built. They are threaded onto a
    COPY of ``CORTEX_REGION.arch`` via :func:`dataclasses.replace` (the module-level ``CORTEX`` constant is NOT
    mutated). The defaults (0.5 / 3.0 / 20.0) are the current ``CORTEX`` values, so the default call reproduces
    the historical arch bit-for-bit. Downstream (``ff.weave.weave``): the mesh discretization is
    ``beads_per_filament = round(length_um/seg_um) + 1`` and the crosslink count is
    ``n_xl = round(n_filaments·density_per_fil)``."""
    base = CORTEX_REGION.arch
    fil = dataclasses.replace(base.filament, n_filaments=n_filaments, seg_um=seg_um, length_um=length_um)
    xl = dataclasses.replace(base.crosslinker, density_per_fil=density_per_fil)
    arch = dataclasses.replace(base, filament=fil, crosslinker=xl, R_um=R_CORTEX_UM)
    return dataclasses.replace(CORTEX_REGION, arch=arch)


def _cortex_regions(cfg: "CellConfig") -> list[RegionSpec]:
    """The cortex region LIST woven into the cell (one region today, two under the mixed-cortex option).

    ``cfg.cortex_arp23_fraction == 0.0`` (default) ⇒ a single formin cortex region IDENTICAL to the historical
    ``_cortex_region(cfg.n_filaments, …)`` call — so ``weave_cell`` stays RNG- and byte-identical (Gate-1
    parity). ``> 0`` ⇒ split the fixed ``cfg.n_filaments`` cortex BUDGET into (n_formin long-formin + n_arp23
    short-branched-Arp2/3) with the count split DERIVED from the Arp2/3 MASS fraction
    (:func:`ac.weave.regions.derive_cortex_arp23_split`), preserving total areal density. The Arp2/3 region
    REUSES the Arp2/3 angle-harmonic branch path (sphere-native placement is the only new geometry)."""
    if cfg.cortex_arp23_fraction <= 0.0:
        return [_cortex_region(
            cfg.n_filaments, seg_um=cfg.cortex_seg_um, length_um=cfg.cortex_length_um,
            density_per_fil=cfg.cortex_density_per_fil)]
    n_formin, n_arp23 = derive_cortex_arp23_split(
        cfg.n_filaments, cfg.cortex_arp23_fraction, cfg.cortex_length_um, cfg.cortex_arp23_length_um)
    formin = _cortex_region(
        n_formin, seg_um=cfg.cortex_seg_um, length_um=cfg.cortex_length_um,
        density_per_fil=cfg.cortex_density_per_fil)
    arp23 = cortex_arp23_region(
        n_arp23, length_um=cfg.cortex_arp23_length_um, seg_um=cfg.cortex_arp23_seg_um,
        mother_fraction=cfg.cortex_arp23_mother_fraction, R_um=R_CORTEX_UM)
    return [formin, arp23]


def _actin_segment_topology(
    fiber_offsets: np.ndarray,
    polarity: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return flat adjacent-node segments plus the direction sign toward each fiber's barbed end.

    Delegates to the canonical I4-weave joint :func:`aleph.components.weave.woven_cell.actin_segment_topology` so the
    production cell binds the head-resolved NMII to the SAME actin↔motor topology the woven network exposes
    (single source of truth; the polarity/barbed-direction gates live on the weave side).
    """
    return actin_segment_topology(fiber_offsets, polarity)


def _actin_node_tangents(
    pos_actin: np.ndarray, fiber_offsets: np.ndarray, polarity: np.ndarray
) -> np.ndarray:
    """Per-actin-node unit tangent ``(N, 3)`` oriented toward each filament's barbed end (vectorised host).

    Central difference ``pos[next] − pos[prev]`` clamped at fiber boundaries (endpoints use the one available
    segment), signed by the per-fiber barbed flag so the tangent points barbed-ward — the polarity the
    straddle placement uses to pick an ANTI-parallel partner filament for the bipolar minifilament.
    """
    off = np.asarray(fiber_offsets, np.int64)
    n = int(pos_actin.shape[0])
    nfib = off.size - 1
    node_fiber = np.repeat(np.arange(nfib), np.diff(off))
    idx = np.arange(n)
    prev = idx - 1
    nxt = idx + 1
    starts = off[:-1]
    ends = off[1:] - 1
    prev[starts] = starts                            # first node of a fiber has no earlier neighbour
    nxt[ends] = ends                                 # last node has no later neighbour
    d = pos_actin[nxt] - pos_actin[prev]
    sign = np.where(np.asarray(polarity)[node_fiber] >= 0, 1.0, -1.0)
    d = d * sign[:, None]
    nrm = np.linalg.norm(d, axis=1, keepdims=True)
    nrm[nrm < 1e-12] = 1.0
    return d / nrm


def _select_antiparallel_partner(
    tree, pos_actin: np.ndarray, tang: np.ndarray, node_fiber: np.ndarray,
    anchor: int, head_offset_um: float,
) -> int:
    """Nearest anti-parallel actin node ~``2·head_offset`` from ``anchor`` on a DIFFERENT fiber, else ``-1``.

    Scores candidates by ``−(t_A·t_B) − |sep − 2·offset|/offset`` so the pick is strongly anti-parallel AND
    close to the ideal ``2·offset`` separation (where the ± heads land exactly on the two filament lines).
    """
    target = 2.0 * head_offset_um
    sep_lo = target * (1.0 - NMII_STRADDLE_SEP_TOL)
    sep_hi = target * (1.0 + NMII_STRADDLE_SEP_TOL)
    a = pos_actin[anchor]
    t_a = tang[anchor]
    f_a = int(node_fiber[anchor])
    cand = tree.query_ball_point(a, sep_hi)
    best, best_score = -1, -np.inf
    for j in cand:
        if int(node_fiber[j]) == f_a:
            continue
        sep = float(np.linalg.norm(pos_actin[j] - a))
        if sep < sep_lo:
            continue
        dot = float(np.dot(t_a, tang[j]))
        if dot > NMII_STRADDLE_ANTIPARALLEL_MAX_DOT:  # not anti-parallel enough
            continue
        score = -dot - abs(sep - target) / head_offset_um
        if score > best_score:
            best_score, best = score, j
    return best


def _build_myosin(
    pos_actin: np.ndarray,
    myo_i: np.ndarray,
    myo_j: np.ndarray,
    device: str,
    *,
    backbone_lp_um: float | None,
    fiber_offsets: np.ndarray | None = None,
    polarity: np.ndarray | None = None,
    straddle_placement: bool = False,
):
    """Build head-resolved minifilaments at the cortex myosin sites → (myo_positions, MyosinForce, counts).

    All particle/bond/head indices are offset to GLOBAL indices in the combined [actin | myosin] node array.
    Heads are allocated UNBOUND (allocate_hand_state default) — a resting, passive-only motor at t0.

    Placement (``straddle_placement``, defect#3 fix): when ON (needs ``fiber_offsets`` + ``polarity``), each
    bipolar minifilament is oriented to STRADDLE two anti-parallel actin filaments ~2·head_offset apart via
    :func:`aleph.components.motor.minifilament_topology.straddle_frame`, so its ± heads land ON the actin they will
    bind (~nm crossbridge extension) and the dipole anti-parallel-anchors. Per-minifilament fallback to the
    legacy geometry (centre = myo-link midpoint, axis = link direction) whenever no partner is found. OFF
    reproduces the legacy midpoint placement for every minifilament.
    """
    topo = MinifilamentTopology(n_bb=NMII_N_BB, n_heads_per_side=NMII_N_SIDE,
                                backbone_length_um=NMII_L_BB_UM, head_offset_um=NMII_HEAD_OFFSET_UM)
    n_mf = int(myo_i.shape[0])
    n_part = topo.n_particles                       # n_bb + 2·N_side
    n_actin = pos_actin.shape[0]
    n_side = topo.n_heads_per_side

    do_straddle = bool(straddle_placement) and fiber_offsets is not None and polarity is not None
    tree = tang = node_fiber = None
    n_straddled = 0
    if do_straddle:
        try:
            from scipy.spatial import cKDTree
        except ImportError:  # graceful degradation — production build must not crash on a missing accelerator
            warnings.warn(
                "scipy unavailable; NMII straddle placement disabled → legacy midpoint geometry "
                "(heads-off-actin defect#3 NOT fixed). Install scipy to enable head-on-actin placement.",
                RuntimeWarning, stacklevel=2)
            do_straddle = False
        else:
            tree = cKDTree(pos_actin)
            tang = _actin_node_tangents(pos_actin, fiber_offsets, polarity)
            node_fiber = np.repeat(np.arange(np.asarray(fiber_offsets).size - 1),
                                   np.diff(np.asarray(fiber_offsets)))

    pos_parts, bb_parts, hb_parts, hn_parts, wd_parts = [], [], [], [], []
    ba_parts, ha_parts = [], []                     # F6 backbone + head-arm angle triples (global indices)
    for m in range(n_mf):
        centre = 0.5 * (pos_actin[myo_i[m]] + pos_actin[myo_j[m]])
        axis = pos_actin[myo_j[m]] - pos_actin[myo_i[m]]
        d = build_minifilament_nodes(topo, centre.astype(np.float64), axis.astype(np.float64))
        if do_straddle:
            partner = _select_antiparallel_partner(
                tree, pos_actin, tang, node_fiber, int(myo_i[m]), NMII_HEAD_OFFSET_UM)
            if partner >= 0:
                a_idx = int(myo_i[m])
                c_s, e_x, e_y = straddle_frame(
                    pos_actin[a_idx], tang[a_idx], pos_actin[partner], tang[partner])
                wd = np.zeros((topo.n_heads, 3), np.float64)
                wd[:n_side] = -e_x                    # + arm walks −axis (outward); − arm walks +axis: contractile
                wd[n_side:] = +e_x
                d = {**d, "positions": topo.placed_positions(c_s, e_x, e_y), "walk_dir": wd}
                n_straddled += 1
        base = n_actin + m * n_part                 # global offset of this minifilament's first particle
        pos_parts.append(d["positions"])
        bb_parts.append(d["backbone_bonds"] + base)
        hb_parts.append(d["head_bonds"] + base)
        hn_parts.append(d["head_node"] + base)
        wd_parts.append(d["walk_dir"])
        ba_parts.append(d["backbone_angles"] + base)
        ha_parts.append(d["head_arm_angles"] + base)
    myo_pos = np.concatenate(pos_parts, axis=0).astype(np.float64)
    backbone_bonds = np.concatenate(bb_parts, axis=0).astype(np.int32)
    head_bonds = np.concatenate(hb_parts, axis=0).astype(np.int32)
    head_node = np.concatenate(hn_parts, axis=0).astype(np.int32)
    walk_dir = np.concatenate(wd_parts, axis=0).astype(np.float64)
    backbone_angles = np.concatenate(ba_parts, axis=0).astype(np.int32) if ba_parts else np.zeros((0, 3), np.int32)
    head_arm_angles = np.concatenate(ha_parts, axis=0).astype(np.int32)
    n_hands = int(head_node.shape[0])               # = n_mf · 2 · N_side

    params = NMIIHandParams()
    params.k_on = wp.float64(NMII_KON_TEST)
    params.k_off0 = wp.float64(NMII_KOFF0)
    params.f0 = wp.float64(NMII_F0)
    params.v0 = wp.float64(NMII_V0_TEST)
    params.f_stall = wp.float64(NMII_F_STALL_TEST)
    params.kappa = wp.float64(NMII_KAPPA_TEST)
    params.k_xb = wp.float64(NMII_K_XB_TEST)
    params.r0_head = wp.float64(NMII_HEAD_OFFSET_UM)
    params.r0_xb = wp.float64(0.0)
    params.capture_radius = wp.float64(NMII_CAPTURE_UM)

    n_seg = max(NMII_N_BB - 1, 1)
    r0_backbone = NMII_L_BB_UM / n_seg
    # Backbone bond stiffness = the params_i0b3 PER-SEGMENT value (k_backbone = k_backbone_factor·k_head_spring),
    # applied by harmonic_bond_kernel to each of the (n_bb−1) bonds (rest = L_bb/(n_bb−1)). The earlier ×(n_bb−1)
    # "rigid-rod / end-to-end" reinterpretation is WITHDRAWN (Codex cross-audit + PI: the YAML/kernel contract is
    # per-bond, not end-to-end — silently scaling ×13 over-stiffens an archived parameter). ⚠ OPEN (PI/KB): a
    # fixed per-bond k is NOT grid-invariant (end-to-end = k/(n_bb−1) varies with n_bb); if PI intends a fixed
    # end-to-end / EA contract it must be stated EXPLICITLY (k_bond = (n_bb−1)·k_end, or EA/rest_segment), not
    # applied here. The analytic reference is corrected to k_eff = k_bond/(n_bb−1) to MATCH this soft chain.
    k_backbone_seg = NMII_K_BACKBONE_TEST
    # F6 bending stiffnesses — DERIVED (backbone_warp Magic-Number-Block), not tuned to a transmission band.
    k_theta_bb = (
        backbone_bending_k_theta(backbone_bending_kappa(backbone_lp_um), r0_backbone)
        if backbone_lp_um is not None else 0.0
    )
    k_theta_arm = arm_orientation_k_theta(NMII_K_XB_TEST, NMII_HEAD_OFFSET_UM)
    bb_d = wp.array(backbone_bonds, dtype=wp.int32, device=device)
    hb_d = wp.array(head_bonds, dtype=wp.int32, device=device)
    hn_d = wp.array(head_node, dtype=wp.int32, device=device)
    ba_d = (
        wp.array(backbone_angles, dtype=wp.int32, device=device)
        if backbone_lp_um is not None and backbone_angles.shape[0] else None
    )
    ha_d = wp.array(head_arm_angles, dtype=wp.int32, device=device)
    myo = MyosinForce(bb_d, hb_d, hn_d, params,
                      k_backbone=k_backbone_seg, r0_backbone=r0_backbone,
                      k_head_spring=NMII_K_HEAD_ARM_TEST, n_hands=n_hands, walk_dir=walk_dir, device=device,
                      backbone_angles=ba_d, head_arm_angles=ha_d,
                      k_theta_backbone=k_theta_bb, k_theta_arm=k_theta_arm, segment_len_um=r0_backbone)
    counts = {"n_minifilaments": n_mf, "n_myosin_particles": myo_pos.shape[0], "n_heads": n_hands,
              "n_backbone_beads": n_mf * NMII_N_BB,
              "n_straddle_placed": n_straddled,
              "straddle_placed_frac": (n_straddled / n_mf) if n_mf else 0.0}
    return myo_pos, myo, counts


def build_cell(cfg: CellConfig | None = None) -> AssembledCell:
    """Compose the physiological t0 cell on the CUDA device (I0-A: runs on the gbook A5000)."""
    cfg = cfg or CellConfig()
    if cfg.erm_density_per_um2 is not None and not cfg.erm_density_source.strip():
        raise ValueError("a non-null ERM density requires an explicit literature/provenance label")
    if cfg.erm_density_mcf7_production and cfg.erm_density_per_um2 is None:
        raise ValueError("MCF7 production ERM provenance requires an explicit density")
    kinetic_values = (
        cfg.erm_k_on_s,
        cfg.erm_k_off0_s,
        cfg.erm_bell_force_pn,
        cfg.erm_capture_radius_um,
    )
    if any(value is not None for value in kinetic_values) and not all(value is not None for value in kinetic_values):
        raise ValueError("ERM Bell kinetics are all-or-none: k_on, k_off0, Bell force, and capture radius required")
    erm_bell_kinetics = None
    if all(value is not None for value in kinetic_values):
        erm_bell_kinetics = ERMBellKinetics(
            k_on_s=float(cfg.erm_k_on_s),
            k_off0_s=float(cfg.erm_k_off0_s),
            bell_force_pn=float(cfg.erm_bell_force_pn),
            capture_radius_um=float(cfg.erm_capture_radius_um),
            source=cfg.erm_kinetics_source,
            rebind_rest_policy=cfg.erm_rebind_rest_policy,
        )
    if cfg.erm_density_mcf7_production and erm_bell_kinetics is None:
        raise ValueError("MCF7 production ERM requires a complete source-gated Bell on/off kinetics contract")
    if cfg.nmii_backbone_lp_um is None and cfg.nmii_backbone_lp_source.strip():
        raise ValueError("NMII backbone persistence-length source requires an explicit value")
    if cfg.nmii_backbone_lp_um is not None:
        if not np.isfinite(cfg.nmii_backbone_lp_um) or cfg.nmii_backbone_lp_um <= 0.0:
            raise ValueError("NMII backbone persistence length must be finite and positive")
        if not cfg.nmii_backbone_lp_source.strip():
            raise ValueError("NMII backbone persistence length requires an explicit provenance source")
    # Resting bound-myosin setpoint: all-or-none, source-gated, NO physiological default (fraction + per-head
    # force are Contract-Graph GAPs). Both None ⇒ OFF (unbound-at-rest parity); either set ⇒ both + source
    # required. Build the frozen contract now so the build refuses a partial/defaulted setpoint before any state.
    resting_setpoint: RestingBoundMyosinSetpoint | None = None
    _resting_fields = (cfg.resting_bound_myosin_fraction, cfg.resting_bound_myosin_force_pn)
    if any(v is not None for v in _resting_fields) and not all(v is not None for v in _resting_fields):
        raise ValueError(
            "resting bound-myosin setpoint is all-or-none: both a fraction AND a per-head force are required"
        )
    if all(v is not None for v in _resting_fields):
        if not cfg.with_myosin:
            raise ValueError("resting bound-myosin setpoint requires with_myosin=True (heads to bind)")
        resting_setpoint = RestingBoundMyosinSetpoint(
            fraction=float(cfg.resting_bound_myosin_fraction),
            per_head_force_pn=float(cfg.resting_bound_myosin_force_pn),
            source=cfg.resting_bound_myosin_source,
            capture_radius_um=cfg.resting_bound_myosin_capture_um,
        )
    wp.init()
    dev = str(wp.get_device(cfg.device))
    if not wp.get_device(cfg.device).is_cuda:
        raise RuntimeError(f"ac.cell requires a CUDA GPU (I0-A). Resolved {dev!r} is not CUDA.")
    rng = np.random.default_rng(cfg.seed)

    # ── 1. cortex mechanical state (weave_cell → CrosslinkedCortex; Gate-1 delegation to ff.weave) ────
    wc = weave_cell(
        _cortex_regions(cfg),
        rng=rng,
        overlap_free=cfg.overlap_free_cortex,
        overlap_mode=cfg.cortex_overlap_mode,
        overlap_span=cfg.cortex_overlap_span,
    )
    cortex = wc.to_crosslinked_cortex()
    net = cortex.net
    pos_actin = np.ascontiguousarray(net.pos, np.float64)
    n_actin = int(net.n_nodes)

    tri = np.ascontiguousarray(net.bend_triples, np.int32)
    alpha = np.ascontiguousarray(_per_triple_alpha(net), np.float64)
    fiber_off = np.ascontiguousarray(net.fiber_offsets, np.int32)
    nodes_per_fiber = np.diff(fiber_off)
    triple_per = np.maximum(nodes_per_fiber - 2, 0)
    triple_off = np.concatenate([[0], np.cumsum(triple_per)]).astype(np.int32)
    if int(triple_off[-1]) != int(tri.shape[0]):
        raise RuntimeError("actin bending triples are not contiguous by fiber")
    seg_per = np.diff(fiber_off) - 1
    seg_off = np.concatenate([[0], np.cumsum(seg_per)]).astype(np.int32)
    srest = np.ascontiguousarray(net.seg_rest, np.float64)
    seg_mean = float(net.seg_rest.mean()) if net.seg_rest.size else 0.5

    xl = np.ascontiguousarray(np.stack([cortex.xl_i, cortex.xl_j], axis=1), np.int32) if cortex.xl_i.size \
        else np.zeros((0, 2), np.int32)
    kxl = np.ascontiguousarray(cortex.xl_k, np.float64)
    r0xl = np.ascontiguousarray(cortex.xl_rest, np.float64)

    # per-node fiber id (for steric same-filament exclusion) — actin nodes
    node_fiber_actin = np.repeat(np.arange(net.n_fibers), np.diff(fiber_off)).astype(np.int32)
    seg_node_a, seg_node_b, seg_polarity = _actin_segment_topology(fiber_off, wc.polarity)
    if seg_node_a.shape[0] != srest.shape[0]:
        raise RuntimeError("actin segment topology does not match the NF2007 segment-rest ledger")

    # ── 2. compose the global node array [actin | myosin | nucleus | membrane] ────────────────────────
    # Each non-actin block is appended at its running offset so its device force kernels address GLOBAL
    # indices; all are marked dormant for EV (active<0) with a distinct fiber_id.
    blocks = [pos_actin]

    # 2a. myosin minifilaments (optional; heads unbound at rest)
    myo = None
    myo_counts = {"n_minifilaments": 0, "n_myosin_particles": 0, "n_heads": 0, "n_backbone_beads": 0}
    n_myo_part = 0
    if cfg.with_myosin and cortex.myo_i.size:
        myo_pos, myo, myo_counts = _build_myosin(
            pos_actin,
            cortex.myo_i,
            cortex.myo_j,
            dev,
            backbone_lp_um=cfg.nmii_backbone_lp_um,
            fiber_offsets=fiber_off,
            polarity=wc.polarity,
            straddle_placement=cfg.nmii_straddle_placement,
        )
        blocks.append(myo_pos)
        n_myo_part = int(myo_pos.shape[0])

    # 2b. deformable-mesh nucleus (retires nucleus_shell_kernel; its mesh is the fluid no-flux inclusion)
    nucleus = nuc_provider = None
    n_nuc = 0
    if cfg.with_nucleus:
        nucleus, nuc_verts, nuc_provider = build_nucleus_compartment(
            pos_actin, n_actin + n_myo_part, dev, params={"subdivisions": cfg.nucleus_subdivisions})
        blocks.append(nuc_verts)
        n_nuc = int(nuc_verts.shape[0])

    # 2c. plasma-membrane Helfrich sheet (rides on the cortex via ERM; the osmotic envelope)
    membrane = None
    n_mem = 0
    if cfg.with_membrane:
        membrane, mem_verts = build_membrane_compartment(
            pos_actin,
            n_actin + n_myo_part + n_nuc,
            dev,
            R_mem=R_CELL_UM,
            params={
                "subdivisions": cfg.membrane_subdivisions,
                "erm_density_per_um2": cfg.erm_density_per_um2,
                "erm_density_source": cfg.erm_density_source,
                # The resting cortical tension can only reach the membrane through a RADIAL ERM (§2/§6): a
                # laterally-offset tether transmits none of the radial turgor load. Turn radial pairing ON with
                # the setpoint so the actively-tensioned cortex holds the membrane through the ERM (no preload).
                "erm_radial_pairing": cfg.erm_radial_pairing or (resting_setpoint is not None),
                "erm_bell_kinetics": erm_bell_kinetics,
            },
        )
        blocks.append(mem_verts)
        n_mem = int(mem_verts.shape[0])

    pos_all = np.concatenate(blocks, axis=0) if len(blocks) > 1 else pos_actin
    n_total = int(pos_all.shape[0])
    n_extra = n_myo_part + n_nuc + n_mem                      # all non-actin nodes

    # steric fiber_id + active mask over ALL nodes: non-actin blocks get distinct ids AND active<0 (skip EV)
    fiber_id = np.concatenate([
        node_fiber_actin,
        np.full(n_myo_part, -2, np.int32), np.full(n_nuc, -3, np.int32), np.full(n_mem, -4, np.int32),
    ]).astype(np.int32)
    active = np.concatenate([
        np.zeros(n_actin, np.int32), np.full(n_extra, -1, np.int32),
    ]).astype(np.int32)   # <0 dormant → EV skips myosin/nucleus/membrane

    # I4 walk_dir hand-off: fiber id of every actin node + each fiber's global barbed-end node
    barbed_actin = barbed_end_node(fiber_off, wc.polarity).astype(np.int32)     # (n_fibers,) global actin idx
    node_fiber_d = wp.array(node_fiber_actin.astype(np.int32), dtype=wp.int32, device=dev)
    barbed_node_d = wp.array(barbed_actin, dtype=wp.int32, device=dev)

    # ── 3. device upload of the combined mechanical state ────────────────────────────────────────────
    pos_d = wp.array(np.ascontiguousarray(pos_all, np.float64), dtype=wp.vec3d, device=dev)
    if myo is not None:
        myo.enable_segment_runtime(
            wp.array(seg_node_a, dtype=wp.int32, device=dev),
            wp.array(seg_node_b, dtype=wp.int32, device=dev),
            wp.array(seg_polarity, dtype=wp.int32, device=dev),
            max_segment_length_um=float(np.max(srest)),
            grid_dim=cfg.motor_segment_grid_dim,
            device=dev,
        )

    # ── 3b. resting bound-myosin SETPOINT (candidate 1) — seed a fraction of heads bound + isometrically loaded
    # onto the actin network at t0, so the cortex carries the resting active tension the ERM transmits to hold
    # the membrane turgor. OFF (default) ⇒ heads unbound at rest (parity), and the PI-GAP is surfaced. The seed
    # plan is pure-host construction (out-of-hot-loop, like ERM/LINC pairing); it is written once into the device
    # segment-hand state and is GPU-resident thereafter (I0-A). ──────────────────────────────────────────────
    resting_seed_ledger: dict[str, object] = {
        "resting_bound_myosin_status": "OFF_UNBOUND_AT_REST",
        "resting_bound_myosin_gap": (
            "PI-GAP: resting bound-myosin fraction + per-head force are absent from the Contract-Graph; "
            "set cfg.resting_bound_myosin_fraction and .resting_bound_myosin_force_pn (with a source) to "
            "activate the resting cortical-tension source. UNSET ⇒ mechanism OFF (heads unbound at t0)."
        ),
        "resting_bound_myosin_n_bound": 0,
    }
    if resting_setpoint is not None and myo is not None:
        capture_um = (
            float(cfg.resting_bound_myosin_capture_um)
            if cfg.resting_bound_myosin_capture_um is not None
            else float(NMII_CAPTURE_UM)
        )
        # GLOBAL head-node indices (build-time host read, out-of-hot-loop) + the segment topology already in host
        # NumPy: the seed planner is a pure-host construction, exactly like the ERM/LINC pairing above.
        head_node_global = myo.head_node.numpy().astype(np.int64)
        seed_plan = plan_resting_bound_heads(
            pos_all, head_node_global, seg_node_a.astype(np.int64), seg_node_b.astype(np.int64),
            seg_polarity.astype(np.int64), resting_setpoint,
            k_xb=float(NMII_K_XB_TEST), r0_xb=0.0, capture_radius_um=capture_um,
        )
        n_seeded = apply_resting_bound_heads(myo.segment_runtime.state, seed_plan, device=dev)
        realized_loads = resting_tangential_load_reference(
            pos_all, head_node_global, seed_plan, k_xb=float(NMII_K_XB_TEST), r0_xb=0.0,
        )
        max_load_err = float(np.max(np.abs(realized_loads - resting_setpoint.per_head_force_pn))) \
            if seed_plan.n_bound else 0.0
        resting_seed_ledger = {
            "resting_bound_myosin_status": "ON_SEEDED_ISOMETRIC",
            "resting_bound_myosin_source": resting_setpoint.source,
            "resting_bound_myosin_fraction_requested": float(resting_setpoint.fraction),
            "resting_bound_myosin_per_head_force_pn": float(resting_setpoint.per_head_force_pn),
            "resting_bound_myosin_n_bound": int(n_seeded),
            "resting_bound_myosin_realized_load_err_pn": max_load_err,
            "resting_bound_myosin_erm_pairing": "RADIAL",
            **{f"resting_bound_myosin_{k}": v for k, v in seed_plan.diagnostics.items()},
        }
    if nuc_provider is not None and cfg.with_pressure:    # mask tracks the LIVE envelope (refresh each remap)
        nuc_provider.bind_live(pos_d, n_actin + n_myo_part, n_nuc)
    mem_provider = None
    if membrane is not None and cfg.with_pressure:       # outer fluid boundary tracks the LIVE membrane mesh
        mem_provider = LiveMeshMembraneMaskProvider(
            pos_d, membrane.node_off, membrane.n_verts, membrane.mesh.faces,
            reference_verts=membrane.mesh.verts)
    f_d = wp.zeros(n_total, dtype=wp.vec3d, device=dev)
    tri_d = wp.array(tri, dtype=wp.int32, device=dev)
    alpha_d = wp.array(alpha, dtype=wp.float64, device=dev)
    foff_d = wp.array(fiber_off, dtype=wp.int32, device=dev)
    toff_d = wp.array(triple_off, dtype=wp.int32, device=dev)
    soff_d = wp.array(seg_off, dtype=wp.int32, device=dev)
    srest_d = wp.array(srest, dtype=wp.float64, device=dev)
    xl_d = wp.array(xl, dtype=wp.int32, device=dev)
    kxl_d = wp.array(kxl, dtype=wp.float64, device=dev)
    r0xl_d = wp.array(r0xl, dtype=wp.float64, device=dev)

    # ── Arp2/3 70° branch-angle junctions (mixed formin+Arp2/3 cortex only) ───────────────────────────
    # to_crosslinked_cortex() carries the weave's (m_after, branch_vertex, daughter_arm) branch triples; here
    # they become a device array so driver._accumulate_all launches the angle-harmonic kernel over them (the
    # SAME 70° Arp2/3 harmonic the lamellipodium's LamellipodiumBranchAngleMechanics enforces — CLAUDE.md
    # worked-example: angle-harmonic + thermal fluctuation, NOT a rigid clamp). The formin-only cortex weaves
    # ZERO branch_triples → these stay None → no branch launch → bit-identical to the pre-wiring build. Node
    # indices are actin-block indices, valid unchanged in the [actin|myosin|nucleus|membrane] global array.
    branch_triples_np = np.ascontiguousarray(cortex.branch_triples, np.int32)
    n_branch = int(branch_triples_np.shape[0])
    if n_branch:
        if int(wc.branch_active.shape[0]) != n_branch:
            raise RuntimeError("branch_active length does not match the cortex branch_triples count")
        branch_triples_d = wp.array(branch_triples_np, dtype=wp.int32, device=dev)
        branch_active_d = wp.array(
            np.ascontiguousarray(wc.branch_active.astype(np.int32)), dtype=wp.int32, device=dev)
    else:
        branch_triples_d = None
        branch_active_d = None

    # Per-node representative control volume — the ACTIN solid skeleton only (pressure -α∇p·V acts on the
    # skeleton; myosin/nucleus/membrane nodes carry their OWN mechanics, so V=0 ⇒ they feel no Biot pressure).
    # The live triangular membrane, not the analytic sphere used to seed it, is the actual discrete closed
    # domain.  Its exact polyhedral volume must therefore be the pressure/FSI quadrature volume; mixing the
    # analytic sphere volume with the triangle pressure load creates a deterministic internal net force.
    analytic_cell_vol = (4.0 / 3.0) * np.pi * R_CELL_UM**3
    cell_vol = (
        signed_volume(mem_verts, membrane.mesh.faces)
        if membrane is not None
        else analytic_cell_vol
    )
    node_vol = np.zeros(n_total, np.float64)
    node_vol[:n_actin] = cell_vol / max(n_actin, 1)
    node_volume_d = wp.array(node_vol, dtype=wp.float64, device=dev)
    solid_active_d = wp.array(active, dtype=wp.int32, device=dev)

    state = SimpleNamespace(pos=pos_d, node_pos=pos_d, node_volume=node_volume_d)

    # ── 4. StericForce (all-fiber WCA; SUBSUMES soft_contact_kernel) ─────────────────────────────────
    steric = None
    if cfg.with_steric:
        gd = (cfg.steric_grid_dim,) * 3
        steric = StericForce(fiber_id=fiber_id, sigma=SIGMA_EV_UM, k_ev=K_EV_PROVISIONAL,
                             device=dev, active=active, force_cap=cfg.steric_force_cap, grid_dim=gd)

    # ── 5. conservative Biot substrate + PressureCoupling (fluid ON from I1 at Π₀) ───────────────────
    # ONE resolved Π₀ for the whole build: the grid seed, the CFL term, the preload and the driver's
    # membrane osmotic balance all read this, so an axis point cannot be applied to some of them and not
    # others. `None` resolves to the gated constant, so the default build is unchanged.
    pi0_pa = PI_0_PA if cfg.turgor_pi0_pa is None else float(cfg.turgor_pi0_pa)
    grid = domain = substrate = membrane_bc = pressure = membrane_pressure = None
    if cfg.with_pressure:
        half = R_CELL_UM + cfg.grid_pad_um
        nx = int(round(2.0 * half / cfg.dx_um)) + 1
        origin = (-half, -half, -half)
        grid = FieldGrid((nx, nx, nx), cfg.dx_um, origin, device=dev)
        # Seed the resting mean turgor.  Uniform p has zero BULK gradient but non-zero membrane traction;
        # MembranePressureTraction below supplies that surface term and NG-3 checks the resulting preload.
        p0 = np.full(grid.shape, pi0_pa, np.float64)
        grid.set_pressure(p0)
        with wp.ScopedDevice(dev):
            grid.p_bar = wp.array(p0, dtype=wp.float64)
        # Production outer boundary = the LIVE Helfrich mesh. The sphere remains only for the explicit
        # --no-membrane diagnostic, never as a runtime substitute when the membrane compartment is present.
        fluid_membrane_provider = mem_provider or StaticSphereMembraneProvider(
            centre=(0.0, 0.0, 0.0), radius=R_CELL_UM)
        domain = Domain(grid, membrane=fluid_membrane_provider, nucleus=None)
        if nuc_provider is not None:               # inner boundary: the LIVE deformable oblate nucleus mesh
            domain.set_nucleus_boundary(nuc_provider)
        domain.classify()
        substrate = BiotSubstrate(grid, mobility=BIOT_MOBILITY, storage_S=BIOT_STORAGE_S, alpha=BIOT_ALPHA)
        domain.bind_storage(substrate.storage_S)
        membrane_bc = MembraneFluxBC(
            grid, L_p=L_P, p_ext=0.0, sigma_refl=1.0,
            pos=pos_d if membrane is not None else None,
            faces=membrane.faces_d if membrane is not None else None,
            n_faces=membrane.n_faces if membrane is not None else 0,
        )
        pressure = PressureCoupling(substrate)
        if membrane is not None:
            membrane_pressure = MembranePressureTraction(
                grid=grid, faces_d=membrane.faces_d, n_faces=membrane.n_faces, p_ext=0.0)

    # ── 6. CFL: kmax over every composed stiffness → dt_mu = 0.1/kmax (inner-solve pseudo-step) ───────
    kmax = float(net.kappa.max()) / seg_mean**3
    if kxl.size:
        kmax = max(kmax, float(kxl.max()))
    if cfg.with_steric:
        from aleph.components.solid.steric_reference import cfl_stiffness_term
        kmax = max(kmax, float(cfl_stiffness_term(SIGMA_EV_UM, steric.epsilon, cfg.steric_force_cap)))
    if myo is not None:
        kmax = max(kmax, myo.cfl_stiffness)   # per-segment rigid-rod backbone + k_xb + F6 bending eff. stiffness
    if nucleus is not None:                                   # bending κ̃/ℓ³ + areal + volume + LINC
        ell = max(nucleus.mean_edge_um, 1.0e-6)
        kmax = max(kmax, nucleus.kappa_tilde / ell**3, nucleus.k_soft, nucleus.k_ac,
                   nucleus.k_vol, nucleus.k_linc)
    pressure_edge_um = 0.0
    if membrane is not None:                                  # bending κ̃/ℓ³ + area tension + ERM
        ell = max(membrane.mean_edge_um, 1.0e-6)
        kmax = max(kmax, membrane.kappa_tilde / ell**3, membrane.gamma_mem, membrane.k_erm)
        pressure_edge_um = ell
    mechanical_kmax = kmax
    if membrane is not None and membrane_pressure is not None:
        # d(p*A)/dx ~ p*ell [pN/um], derived from the triangle pressure force Jacobian at resting Pi_0.
        kmax = max(kmax, pi0_pa * ell)
    dt_mu = 0.1 / kmax

    # ── 7. physiological preload capacity + population ledger + peak GPU bytes ──────────────────────
    preload_ledger: dict[str, object] = {}
    if membrane is not None:
        membrane_area = triangle_surface_area(mem_verts, membrane.mesh.faces)
        preload = evaluate_preload_capacity(
            pressure_pa=pi0_pa if cfg.with_pressure else 0.0,
            radius_um=R_CELL_UM,
            membrane_area_um2=membrane_area,
            membrane_tension_pn_per_um=membrane.gamma_mem if membrane.with_area_tension else 0.0,
            erm_tether_count=membrane.n_erm,
            erm_rupture_force_pn=membrane.f_rupt,
        )
        preload_ledger.update(preload.ledger_fields())
        preload_ledger.update({
            # ``membrane_surface.erm_rupture_force`` is the continuum membrane tube/tether extraction scale,
            # not a measured single-ezrin rupture force.  Preserve the historical arithmetic as a diagnostic,
            # but never let it authorize a molecular population or NG-3 production capacity.
            "preload_capacity_force_basis": "DIAGNOSTIC_MEMBRANE_TETHER_FORCE_NOT_SINGLE_ERM",
            "erm_pairing_mode": membrane.erm_pairing_mode,
            "erm_density_target_per_um2": membrane.erm_density_target_per_um2,
            "erm_density_source": membrane.erm_density_source or "ABSENT_FROM_CONTRACT_GRAPH",
            "erm_density_mcf7_production": bool(cfg.erm_density_mcf7_production),
            "erm_kinetics_mode": (
                "BELL_SLIP_ON_OFF" if membrane.erm_bell_kinetics is not None
                else "DIAGNOSTIC_HARD_THRESHOLD_NOT_PRODUCTION"
            ),
            "erm_kinetics_source": (
                membrane.erm_bell_kinetics.source if membrane.erm_bell_kinetics is not None
                else "ABSENT_FROM_CONTRACT_GRAPH"
            ),
            "erm_bell_k_on_s": (
                membrane.erm_bell_kinetics.k_on_s if membrane.erm_bell_kinetics is not None else None
            ),
            "erm_bell_k_off0_s": (
                membrane.erm_bell_kinetics.k_off0_s if membrane.erm_bell_kinetics is not None else None
            ),
            "erm_bell_force_pn": (
                membrane.erm_bell_kinetics.bell_force_pn if membrane.erm_bell_kinetics is not None else None
            ),
            "erm_capture_radius_um": (
                membrane.erm_bell_kinetics.capture_radius_um
                if membrane.erm_bell_kinetics is not None else None
            ),
            "erm_rebind_rest_policy": (
                membrane.erm_bell_kinetics.rebind_rest_policy
                if membrane.erm_bell_kinetics is not None else "UNSPECIFIED"
            ),
            "physiological_preload_status": "BLOCKED_UNSOURCED_SINGLE_ERM_CAPACITY_BASIS",
        })
    ledger = wc.population_ledger()
    ledger.update({
        "n_actin_nodes": n_actin, "n_total_nodes": n_total, "n_myosin_particles": n_myo_part,
        "n_actin_segments": int(seg_node_a.shape[0]),
        "n_nucleus_nodes": n_nuc, "n_membrane_nodes": n_mem,
        "nucleus_faces": nucleus.n_faces if nucleus is not None else 0,
        "nucleus_hinges": nucleus.n_hinges if nucleus is not None else 0,
        "nucleus_linc_tethers": nucleus.n_linc if nucleus is not None else 0,
        "membrane_faces": membrane.n_faces if membrane is not None else 0,
        "membrane_erm_tethers": membrane.n_erm if membrane is not None else 0,
        "n_bending_triples": int(tri.shape[0]), "n_crosslink_links": int(xl.shape[0]),
        "myosin_attachment_primitive": "POINT_TO_SEGMENT_BARYCENTRIC" if myo is not None else "ABSENT",
        "myosin_query_segment_count": int(seg_node_a.shape[0]) if myo is not None else 0,
        "myosin_segment_query_grid_dim": int(cfg.motor_segment_grid_dim) if myo is not None else 0,
        **{f"myosin_{k}": v for k, v in myo_counts.items()},
        "myosin_anchor_mode": "POINT_TO_SEGMENT_DYNAMIC" if myo is not None else "ABSENT",
        "nmii_backbone_lp_um": cfg.nmii_backbone_lp_um,
        "nmii_backbone_lp_source": cfg.nmii_backbone_lp_source or "ABSENT_FROM_CONTRACT_GRAPH",
        "nmii_backbone_bending_status": (
            "SOURCE_GROUNDED" if cfg.nmii_backbone_lp_um is not None
            else "BLOCKED_UNSOURCED_PHYSICAL_MAGNITUDE"
        ),
        **resting_seed_ledger,
        "field_grid_cells": int(np.prod(grid.shape)) if grid is not None else 0,
        "resting_pressure_Pa": pi0_pa if grid is not None else 0.0,
        # The axis point actually built, beside the gated claim. When these differ, PI_0_PROVENANCE below
        # describes the CLAIM and not this run — a reader must see both or the provenance misleads.
        "turgor_pi0_axis_pa": None if cfg.turgor_pi0_pa is None else float(cfg.turgor_pi0_pa),
        "turgor_pi0_is_gated_default": bool(cfg.turgor_pi0_pa is None),
        # Π₀ + convenience-default provenance travel WITH the artifact (PI 2026-07-25 claim hygiene):
        # a reader of this ledger can no longer see 40 Pa / 0.5 µm / 3.0 µm / 20 xl-per-fil without
        # also seeing that each is CONVENIENCE and what sourced value it deviates from.
        **PI_0_PROVENANCE,
        **convenience_labels_record(),
        "pressure_control_volume_um3": cell_vol,
        "analytic_sphere_volume_um3": analytic_cell_vol,
        "kmax_pN_per_um": kmax, "mechanical_kmax_pN_per_um": mechanical_kmax,
        "pressure_edge_um": pressure_edge_um, "dt_mu": dt_mu,
        **preload_ledger,
    })
    try:
        wp.synchronize_device(dev)
        d = wp.get_device(dev)
        ledger["gpu_bytes_used"] = int(d.total_memory - d.free_memory)
        ledger["gpu_bytes_total"] = int(d.total_memory)
        from warp._src.context import runtime as warp_runtime
        ledger["warp_mempool_used_current_bytes"] = int(
            warp_runtime.core.wp_cuda_device_get_mempool_used_mem_current(d.ordinal))
        ledger["warp_mempool_used_high_bytes"] = int(
            warp_runtime.core.wp_cuda_device_get_mempool_used_mem_high(d.ordinal))
        ledger["gpu_memory_accounting_status"] = "PARTIAL_MEMPOOL_ONLY"
        ledger["gpu_memory_accounting_reason"] = (
            "Warp mempool high-water excludes non-mempool CUDA allocations; exact whole-device peak remains open"
        )
    except Exception as exc:  # noqa: BLE001 - a HARD ledger may be open, but it may never disappear silently
        ledger["gpu_memory_accounting_status"] = "UNAVAILABLE"
        ledger["gpu_memory_accounting_error"] = f"{type(exc).__name__}: {exc}"

    return AssembledCell(
        cfg=cfg, device=dev, n_actin=n_actin, n_total=n_total, pos_d=pos_d, f_d=f_d, state=state,
        tri_d=tri_d, alpha_d=alpha_d, n_tri=int(tri.shape[0]),
        xl_d=xl_d, kxl_d=kxl_d, r0xl_d=r0xl_d, n_xl=int(xl.shape[0]),
        foff_d=foff_d, toff_d=toff_d, soff_d=soff_d, srest_d=srest_d, n_fibers=int(net.n_fibers),
        max_fiber_nodes=int(nodes_per_fiber.max()) if nodes_per_fiber.size else 0,
        steric=steric, pressure=pressure, myosin=myo, nucleus=nucleus, membrane=membrane,
        membrane_pressure=membrane_pressure,
        # Hand out the nucleus provider ONLY when it was bound. `bind_live` above is gated on
        # `cfg.with_pressure`, so passing it unconditionally gave every `is not None` guard in the
        # driver a provider whose `_live` is None, and the first property read RAISED
        # ("MeshNucleusMaskProvider has no live device binding"). That made
        # `with_nucleus=True, with_pressure=False` unrunnable — the exact configuration the
        # pressure A/B needs. The membrane already does this correctly at `mem_provider` above;
        # this is the nucleus catching up, and it is byte-identical whenever pressure is on.
        nucleus_mask_provider=nuc_provider if cfg.with_pressure else None,
        membrane_mask_provider=mem_provider,
        node_volume_d=node_volume_d, solid_active_d=solid_active_d,
        node_fiber_d=node_fiber_d, barbed_node_d=barbed_node_d,
        grid=grid, domain=domain, substrate=substrate, membrane_bc=membrane_bc,
        kmax=kmax, dt_mu=dt_mu, mechanical_kmax=mechanical_kmax,
        pressure_edge_um=pressure_edge_um, convergence_length_um=seg_mean,
        branch_triples_d=branch_triples_d, branch_active_d=branch_active_d, n_branch=n_branch,
        ledger=ledger)
