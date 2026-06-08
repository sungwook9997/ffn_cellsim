"""Compartment registry + recipe spine for the explicit single-cell platform.

H.7 platform unit (2026-06-08). ONE declarative catalogue of every
mechanobiology compartment ffn_cellsim knows how to assemble, plus a recipe
layer that composes a subset of them into a runnable
``configs/mcf7_baseline.yaml``-shaped manifest for the existing
:func:`ffn_sim.cell.manifest.build_baseline_cell` path.

Design intent (additive overlay, NOT a rewrite)
-----------------------------------------------
The existing build path (``manifest.resolve_baseline`` -> ``Cell.build``) is a
hand-wired if/branch assembler that works and is bit-validated. This module does
**not** replace it. It is a metadata + orchestration layer that:

1. Names every compartment once, with its default ON/OFF, hot-path priority,
   *Compartment Performance Contract*, GPU readiness, gamma-contamination
   contract, and literature citations / open PI decisions.
2. Translates a high-level **recipe** (``suspended_round``, ``adherent_passive``,
   ...) into a manifest the *unchanged* loader consumes — but only for
   compartments that are **LIVE** (already wired into ``Cell.build``).
3. **Refuses to silently wire un-ratified physics.** Enabling an ``EXPERIMENTAL``
   or ``STUB`` compartment through a recipe raises :class:`UnratifiedCompartmentError`
   (the same posture as ``resolve_baseline`` raising ``NotImplementedError`` on a
   PI-gated optional). The recipe still *declares* the intent (so a recipe is the
   one place the whole intended cell is listed) — composition just defers it.

Hard rules honoured (CLAUDE.md / AGENTS.md):
- Every compartment is **default-OFF** except the physiological baseline set
  (cytoplasm, enclosed-volume turgor, nucleus, membrane-surface) and the
  always-assembled actomyosin core (cortex, crosslinkers, myosin). Disabling a
  baseline-required compartment via a recipe raises (no silent baseline drop),
  except the single sanctioned ``allow_no_nucleus`` cortical-tension exception.
- No compartment may inject tension into the cortical-gamma estimator. Each spec
  declares ``gamma_contaminating`` + the ``denylist_bond_types`` its bonds must
  occupy; :func:`CompartmentRegistry.gamma_denylist` returns the union, which a
  test cross-checks against ``cortex.cortical_tension.ADHESION_BOND_TYPES``.
- The surface manifold is ``geometry_only`` — it must add **zero** force-bearing
  edges/area springs (mesh-as-physics is rejected). The registry encodes that as
  a hard invariant a test asserts.

This module imports only stdlib + PyYAML; it has **no** HOOMD dependency, so it is
cheap to import in tests and tooling. The ``resolve_ref`` field is an informational
dotted path to the real ``resolve_*`` function — the registry never calls it
(resolution still flows through ``manifest.py`` for LIVE compartments).
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

import yaml

_CONFIG_DIR = Path(__file__).resolve().parents[1] / "configs"
_RECIPE_DIR = _CONFIG_DIR / "recipes"


# ---------------------------------------------------------------------------
# Enumerations
# ---------------------------------------------------------------------------
class HotPathPriority(Enum):
    """Performance criticality of a compartment's per-step / per-batch work.

    P0 — always-on production per-step force/updater. MUST have a GPU-resident or
         native path planned before merge.
    P1 — recipe-on per-step force/updater likely used in production. CPU allowed
         only with explicit optimization debt + microbench.
    P2 — batch updater / binding / broad-phase search. Avoid per-step
         ``cpu_local_snapshot``; benchmark if CPU retained.
    P3 — one-time layout / diagnostics / visualization. CPU/NumPy acceptable.
    """

    P0 = "P0"
    P1 = "P1"
    P2 = "P2"
    P3 = "P3"


class GpuPath(Enum):
    """Where a compartment's hot work runs *today*."""

    NONE = "none"                 # no per-step/per-batch hot work
    CPU = "cpu"                   # cpu_local_snapshot / NumPy host path
    CUPY = "cupy"                 # gpu_local_snapshot + cupy, device-resident
    NATIVE_PLUGIN = "native"      # compiled C++/CUDA HOOMD plugin
    HOOMD_BUILTIN = "builtin"     # plain HOOMD bond/angle/pair force (already GPU)


class CompartmentStatus(Enum):
    """Wiring maturity of a compartment."""

    CORE = "core"                 # always assembled (cortex actomyosin spine)
    LIVE = "live"                 # wired into Cell.build, togglable, validated
    GEOMETRY = "geometry"         # geometry/broad-phase only, force-free, not in gamma
    EXPERIMENTAL = "experimental" # module exists, default-OFF, enabling raises (gate pending)
    STUB = "stub"                 # resolver/contract authored, mechanism behind a PI decision


# ---------------------------------------------------------------------------
# Performance + GPU contracts
# ---------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class PerformanceContract:
    """The mandatory Compartment Performance Contract (per-compartment, pre-merge).

    Counts are at the ratified ×40 mesoscale (``n_fil≈1000``) and at native scale
    (``n_fil≈38000``). Integer counts may be given as ``int`` or as a short string
    formula (e.g. ``"~2*n_xl"``) when they scale with a runtime parameter.
    """

    particle_types_added: tuple[str, ...]
    n_particles_mesoscale: int | str          # at n_fil ~ 1000
    n_particles_native: int | str             # at n_fil ~ 38000
    n_bonds: int | str
    n_angles: int | str
    per_step_force: bool
    per_batch_updater: bool
    uses_cpu_local_snapshot: bool
    uses_broad_phase: bool                     # cKDTree / global neighbour search
    expected_on_recipes: tuple[str, ...]
    hot_path_priority: HotPathPriority
    gpu_path_now: GpuPath
    native_forcecompute_candidate: bool
    bottleneck_risk: str


@dataclass(slots=True, frozen=True)
class GpuReadiness:
    """GPU/native readiness note for a compartment."""

    device_resident_now: bool
    env_flag: str | None = None               # e.g. "FFN_GPU_DEVICE_COMPARTMENTS"
    parity_test: str | None = None            # path/name of CPU<->GPU parity test
    optimization_debt: str | None = None      # explicit debt if CPU-only on a hot path
    microbench_hook: str | None = None        # profiling entry point if P0/P1


# ---------------------------------------------------------------------------
# Compartment spec
# ---------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class CompartmentSpec:
    """Declarative spec for one mechanobiology compartment.

    The registry is the single source of truth for these. ``manifest_path`` is the
    dotted key in a ``mcf7_baseline.yaml``-shaped manifest whose ``enabled`` flag
    this compartment toggles (``None`` for core/geometry/experimental compartments
    that have no manifest slot). ``resolve_ref`` is informational only.
    """

    name: str
    category: str                              # cortex|adhesion|membrane|internal|volume|multicell|geometry
    status: CompartmentStatus
    enabled_default: bool                       # True ONLY for CORE + baseline-required
    summary: str
    manifest_path: tuple[str, ...] | None
    resolve_ref: str | None
    performance_contract: PerformanceContract
    gpu_readiness: GpuReadiness
    gamma_contaminating: bool                   # do its bonds bias cortical gamma if counted?
    denylist_bond_types: tuple[str, ...] = ()   # bond types/prefixes that MUST be on the gamma denylist
    geometry_only: bool = False                 # manifold: zero force-bearing edges/area springs
    baseline_required: bool = False             # part of the physiological baseline (cannot silently drop)
    requires: tuple[str, ...] = ()              # other compartments that must be ON
    citations: tuple[str, ...] = ()
    pi_decisions: tuple[str, ...] = ()          # open PI decisions blocking activation
    sanity_gate_ref: str | None = None
    notes: str = ""

    @property
    def is_togglable(self) -> bool:
        """True if a recipe can flip this compartment via a manifest slot."""
        return self.manifest_path is not None and self.status in (
            CompartmentStatus.LIVE,
        )


class UnratifiedCompartmentError(NotImplementedError):
    """Raised when a recipe tries to wire an EXPERIMENTAL/STUB compartment.

    Mirrors ``resolve_baseline`` raising on a PI-gated optional: the physics is
    coded (or stubbed) and default-OFF, but turning it on through the production
    loader is a gate-contract change that needs PI sign-off — not a silent wire-up.
    """


class BaselineDropError(ValueError):
    """Raised when a recipe tries to disable a physiological-baseline compartment."""


# ===========================================================================
# The catalogue
# ===========================================================================
# Counts/risks for LIVE compartments are from the 2026-06-08 build-path survey;
# counts for EXPERIMENTAL/STUB compartments are design estimates to be confirmed
# by the authoring agent's microbench (flagged in bottleneck_risk).

def _live(**kw: Any) -> PerformanceContract:
    return PerformanceContract(**kw)


_SPECS: tuple[CompartmentSpec, ...] = (
    # ----------------------- CORE actomyosin spine -----------------------
    CompartmentSpec(
        name="cortex",
        category="cortex",
        status=CompartmentStatus.CORE,
        enabled_default=True,
        summary="Explicit actin cortex shell (semiflexible filament backbone + angle bending).",
        manifest_path=None,
        resolve_ref="ffn_sim.cortex.cortex.resolve_h3_derived",
        performance_contract=_live(
            particle_types_added=("actin_cortex",),
            n_particles_mesoscale="~n_fil*beads_per_fil (~1e4)",
            n_particles_native="~3.8e5",
            n_bonds="~cortex backbone bonds",
            n_angles="~cortex bending angles",
            per_step_force=False,            # bond/angle forces are HOOMD-builtin
            per_batch_updater=False,
            uses_cpu_local_snapshot=False,
            uses_broad_phase=False,
            expected_on_recipes=("all",),
            hot_path_priority=HotPathPriority.P0,
            gpu_path_now=GpuPath.HOOMD_BUILTIN,
            native_forcecompute_candidate=False,
            bottleneck_risk="HOOMD-builtin harmonic bond/angle — GPU-native already.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=True, optimization_debt=None),
        gamma_contaminating=False,           # cortex bonds ARE the gamma signal (counted)
        baseline_required=False,             # CORE (always assembled), not a baseline toggle
        citations=("KU-1.24 WLC discrete H", "AFINES_ALGORITHM_NOTES.md"),
        sanity_gate_ref="ffn_sim/tests/test_cortex.py",
        notes="The cortical-gamma estimator measures over these bonds — never denylist them.",
    ),
    CompartmentSpec(
        name="crosslinkers",
        category="cortex",
        status=CompartmentStatus.CORE,
        enabled_default=True,
        summary="Explicit alpha-actinin (slip) + filamin (catch-slip) crosslink head-pairs.",
        manifest_path=None,
        resolve_ref="ffn_sim.cortex.crosslinkers.resolve_crosslinkers",
        performance_contract=_live(
            particle_types_added=("xlink_head",),
            n_particles_mesoscale="2*n_xl",
            n_particles_native="2*n_xl(native)",
            n_bonds="n_xl intra + dynamic xlink_attach_b{i}",
            n_angles=0,
            per_step_force=False,
            per_batch_updater=True,          # XlinkBondUpdater (Bell-Evans bind/unbind)
            uses_cpu_local_snapshot=True,
            uses_broad_phase=True,           # nearest-actin candidate search
            expected_on_recipes=("all",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Batch updater on cpu_local_snapshot; P3 binder GPU-port pending.",
        ),
        gpu_readiness=GpuReadiness(
            device_resident_now=False, optimization_debt="binder host-sync; GPU port = P3 (GPU_MAIN_PORT_2026-05-31.md)",
        ),
        gamma_contaminating=False,           # xlink_intra/attach ARE counted in gamma
        baseline_required=False,             # CORE (always assembled), not a baseline toggle
        citations=("Ferrer 2008 alpha-actinin k_off0=0.066/s", "Pereverzev filamin catch-slip"),
        sanity_gate_ref="ffn_sim/tests/test_crosslinkers.py",
    ),
    CompartmentSpec(
        name="myosin",
        category="cortex",
        status=CompartmentStatus.CORE,
        enabled_default=True,
        summary="Stam-Hocky multi-head bipolar minifilament motors (grip-walk stepping, Hill F-v).",
        manifest_path=None,
        resolve_ref="ffn_sim.cortex.myosin.resolve_cortex_myosin",
        performance_contract=_live(
            particle_types_added=("cortex_myosin_backbone", "cortex_myosin_head"),
            n_particles_mesoscale="n_motors*(n_backbone+2*n_heads)",
            n_particles_native="(native motor count)",
            n_bonds="backbone + head-attach bins",
            n_angles=0,
            per_step_force=False,
            per_batch_updater=True,          # MyosinStepUpdater
            uses_cpu_local_snapshot=True,
            uses_broad_phase=True,
            expected_on_recipes=("all",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Batch updater on cpu_local_snapshot; binder GPU-port P3.",
        ),
        gpu_readiness=GpuReadiness(
            device_resident_now=False, optimization_debt="binder host-sync; GPU port = P3",
        ),
        gamma_contaminating=False,           # myosin bonds ARE the active-gamma signal
        baseline_required=False,             # CORE (always assembled), not a baseline toggle
        citations=("Chugh 2017 2-5 pN/motor", "Billington 2013 301nm contour", "Nie 2015 0.625/um2"),
        sanity_gate_ref="ffn_sim/tests/test_myosin.py",
        notes="Active-gamma = gamma(myosin ON) - gamma(myosin OFF). Never fold turgor in (Gate-B contract).",
    ),
    # ----------------------- Physiological baseline ----------------------
    CompartmentSpec(
        name="cytoplasm",
        category="volume",
        status=CompartmentStatus.LIVE,
        enabled_default=True,
        summary="Tier-1 effective-viscosity drag (MCF7 65.9 Pa.s, NOT water) via per-type BAOAB gamma.",
        manifest_path=("compartments", "cytoplasm"),
        resolve_ref="ffn_sim.cell.cytoplasm.resolve_cytoplasm",
        performance_contract=_live(
            particle_types_added=(),
            n_particles_mesoscale=0,
            n_particles_native=0,
            n_bonds=0,
            n_angles=0,
            per_step_force=False,            # modifies BAOAB gamma_map, no force compute
            per_batch_updater=False,
            uses_cpu_local_snapshot=False,
            uses_broad_phase=False,
            expected_on_recipes=("all",),
            hot_path_priority=HotPathPriority.P0,
            gpu_path_now=GpuPath.HOOMD_BUILTIN,
            native_forcecompute_candidate=False,
            bottleneck_risk="None — scales the BAOAB drag coefficient; folded into the integrator.",
        ),
        gpu_readiness=GpuReadiness(
            device_resident_now=True, env_flag="FFN_GPU_DEVICE_BAOAB",
            parity_test="ffn_sim/tests/test_baoab_device.py",
        ),
        gamma_contaminating=False,
        baseline_required=True,
        citations=("Hu 2024 MCF7 cytoplasm 65.9 Pa.s",),
        sanity_gate_ref="ffn_sim/tests/test_h10_cytoplasm.py",
        notes="Physiological-baseline HARD rule: production must NOT silently drop to water drag.",
    ),
    CompartmentSpec(
        name="enclosed_volume",
        category="volume",
        status=CompartmentStatus.LIVE,
        enabled_default=True,
        summary="Static osmotic turgor (Pi_0=133 Pa) pre-tensions the cortex shell (Young-Laplace).",
        manifest_path=("compartments", "enclosed_volume"),
        resolve_ref="ffn_sim.cortex.enclosed_volume.resolve_enclosed_volume",
        performance_contract=_live(
            particle_types_added=(),
            n_particles_mesoscale=0,
            n_particles_native=0,
            n_bonds=0,
            n_angles=0,
            per_step_force=True,             # EnclosedVolumePressure md.force.Custom
            per_batch_updater=False,
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,
            expected_on_recipes=("all",),
            hot_path_priority=HotPathPriority.P0,
            gpu_path_now=GpuPath.CUPY,
            native_forcecompute_candidate=True,
            bottleneck_risk="Per-step Custom force; GPU/native swap exists (1.07x comp, 2.89x full).",
        ),
        gpu_readiness=GpuReadiness(
            device_resident_now=True, env_flag="FFN_GPU_DEVICE_COMPARTMENTS",
            parity_test="bit-exact vs CPU (H7_NATIVE_FULLCELL_GO_2026-06-07.md)",
            microbench_hook="ffn_sim/scripts/compartment_force_profile.py",
        ),
        gamma_contaminating=False,           # passive turgor reported as a SEPARATE gamma channel
        baseline_required=True,
        citations=("Young-Laplace Pi_0=2*gamma/R", "Stewart 2011 (lower alt 40 Pa, superseded)"),
        sanity_gate_ref="ffn_sim/tests/test_enclosed_volume.py",
        notes="Gate-B MUST report active vs passive(turgor) gamma decomposition separately.",
    ),
    CompartmentSpec(
        name="nucleus",
        category="internal",
        status=CompartmentStatus.LIVE,
        enabled_default=True,
        summary="Soft chromatin + stiff-lamin two-regime radial confinement (E_nuc 4.7 kPa).",
        manifest_path=("compartments", "nucleus"),
        resolve_ref="ffn_sim.cell.nucleus.resolve_nucleus",
        performance_contract=_live(
            particle_types_added=("nucleus_bead",),
            n_particles_mesoscale="n_beads (3000)",
            n_particles_native="n_beads",
            n_bonds=0,
            n_angles=0,
            per_step_force=True,             # NucleusConfinement md.force.Custom
            per_batch_updater=False,
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,
            expected_on_recipes=("all but sanctioned no-nucleus cortical-tension build",),
            hot_path_priority=HotPathPriority.P0,
            gpu_path_now=GpuPath.CUPY,
            native_forcecompute_candidate=True,
            bottleneck_risk="Stiff lamin sets the CFL bottleneck; GPU swap exists.",
        ),
        gpu_readiness=GpuReadiness(
            device_resident_now=True, env_flag="FFN_GPU_DEVICE_COMPARTMENTS",
            parity_test="ffn_sim/cell/nucleus_confinement_gpu.py",
        ),
        gamma_contaminating=False,
        baseline_required=True,
        citations=("KU-3.B2.1 in-situ band [1e3,1e4] Pa", "ratio_lamin>=1.4"),
        sanity_gate_ref="ffn_sim/tests/test_nucleus.py",
        notes="Sanctioned no-nucleus exception ONLY for the suspended cortical-tension observable (allow_no_nucleus).",
    ),
    CompartmentSpec(
        name="membrane_surface",
        category="membrane",
        status=CompartmentStatus.LIVE,
        enabled_default=True,
        summary="Young-Laplace area-elastic plasma-membrane surface tension (gamma_mem 0.10 mN/m).",
        manifest_path=("compartments", "membrane_surface"),
        resolve_ref="ffn_sim.cell.membrane_surface.resolve_membrane_surface",
        performance_contract=_live(
            particle_types_added=(),
            n_particles_mesoscale=0,
            n_particles_native=0,
            n_bonds=0,
            n_angles=0,
            per_step_force=True,             # MembraneSurfaceTension md.force.Custom
            per_batch_updater=False,
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,
            expected_on_recipes=("all",),
            hot_path_priority=HotPathPriority.P0,
            gpu_path_now=GpuPath.CUPY,
            native_forcecompute_candidate=True,
            bottleneck_risk="Per-step Custom force; GPU swap exists.",
        ),
        gpu_readiness=GpuReadiness(
            device_resident_now=True, env_flag="FFN_GPU_DEVICE_COMPARTMENTS",
            parity_test="ffn_sim/cell/membrane_surface_gpu.py",
        ),
        gamma_contaminating=False,           # membrane tension = separate gamma_membrane channel (~0.10 mN/m)
        baseline_required=True,
        citations=("KU-3.B1 mid-band 0.10 mN/m",),
        sanity_gate_ref="ffn_sim/tests/test_membrane_surface.py",
        notes="Gate-B should report gamma_membrane as the 4th channel (currently sometimes omitted).",
    ),
    # ----------------------- Adhesion (LIVE optionals) -------------------
    CompartmentSpec(
        name="fa",
        category="adhesion",
        status=CompartmentStatus.LIVE,
        enabled_default=False,
        summary="Focal-adhesion molecular clutch: integrin-ligand Pereverzev catch + fa_actin_clutch.",
        manifest_path=("optional_subsystems", "fa"),
        resolve_ref="ffn_sim.bridge.fa.resolve_h4",
        performance_contract=_live(
            particle_types_added=("integrin", "substrate_ligand"),
            n_particles_mesoscale="n_int + n_lig",
            n_particles_native="(scaled)",
            n_bonds="integrin_ligand_b{i} (256 r0-bins) + fa_actin_clutch_b{i}",
            n_angles=0,
            per_step_force=False,
            per_batch_updater=True,          # IntegrinBondUpdater + SubstrateLigandPin
            uses_cpu_local_snapshot=True,
            uses_broad_phase=True,
            expected_on_recipes=("adherent_passive", "adherent_active_spread", "confined_migration"),
            hot_path_priority=HotPathPriority.P1,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Batch updater + ligand pin on cpu_local_snapshot; CPU debt on P1 path.",
        ),
        gpu_readiness=GpuReadiness(
            device_resident_now=False,
            optimization_debt="IntegrinBondUpdater + SubstrateLigandPin host-sync; needs P1 microbench before production scale.",
        ),
        gamma_contaminating=True,
        denylist_bond_types=("integrin_ligand", "fa_actin_clutch"),
        requires=(),
        citations=("Chan-Odde 2008 clutch", "Pereverzev catch-bond", "FA_INTEGRIN_OVERLOAD_FIX r0-bins=256"),
        pi_decisions=("Adhered build needs equilibration prelude (raw free run trips BAOAB guard).",),
        sanity_gate_ref="ffn_sim/tests/test_fa_integration.py",
        notes="Per-bond r0-bin at bind => force-free birth (overload fix). Clutch bonds are denylisted from gamma.",
    ),
    CompartmentSpec(
        name="rigid_ligand_coating",
        category="adhesion",
        status=CompartmentStatus.LIVE,
        enabled_default=False,
        summary="Rigid (glass/thin col-I coating) ligand anchor: SubstrateLigandPin, k_sub->inf.",
        manifest_path=None,                  # implied by fa ON + substrate OFF
        resolve_ref="ffn_sim.cell.cell.SubstrateLigandPin",
        performance_contract=_live(
            particle_types_added=(),
            n_particles_mesoscale=0,
            n_particles_native=0,
            n_bonds=0,
            n_angles=0,
            per_step_force=False,
            per_batch_updater=True,          # SubstrateLigandPin resets ligand z each step
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,
            expected_on_recipes=("adherent_passive", "adherent_active_spread"),
            hot_path_priority=HotPathPriority.P1,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Per-step ligand pin on cpu_local_snapshot.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="ligand pin host-sync"),
        gamma_contaminating=False,           # pins ligand positions; carries no cortical bonds
        requires=("fa",),
        citations=("glass E~64-70 GPa => rigid limit; col-I in PBS is a thin 2D coating",),
        notes="Default rigid-dish adhesion. The compliant finite-k_sub spring is the separate `substrate` compartment.",
    ),
    CompartmentSpec(
        name="substrate",
        category="adhesion",
        status=CompartmentStatus.LIVE,
        enabled_default=False,
        summary="Compliant ECM substrate spring (finite k_sub) replacing the rigid ligand pin (soft-gel variant).",
        manifest_path=("optional_subsystems", "substrate"),
        resolve_ref="ffn_sim.ecm.substrate.resolve_substrate",
        performance_contract=_live(
            particle_types_added=(),
            n_particles_mesoscale=0,
            n_particles_native=0,
            n_bonds=0,
            n_angles=0,
            per_step_force=True,             # SubstrateAnchorSpring md.force.Custom
            per_batch_updater=False,
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,
            expected_on_recipes=("(soft-gel ECM variant only)",),
            hot_path_priority=HotPathPriority.P1,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=True,
            bottleneck_risk="Per-step spring on cpu_local_snapshot.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="anchor spring host-sync"),
        gamma_contaminating=False,
        requires=("fa",),
        citations=("Chan-Odde clutch k_sub ~ 1 pN/nm per kPa",),
        pi_decisions=("k_sub REQUIRED (no magic default) — only for a defined soft gel (kPa PAA).",),
        sanity_gate_ref="ffn_sim/tests/test_substrate_spring.py",
    ),
    # ----------------------- Active spreading (LIVE optionals) -----------
    CompartmentSpec(
        name="lamellipodium",
        category="adhesion",
        status=CompartmentStatus.LIVE,
        enabled_default=False,
        summary="Arp2/3 branched leading edge: WAVE pin, barbed-end elongation, 70deg branching, capping.",
        manifest_path=("optional_subsystems", "lamellipodium"),
        resolve_ref="ffn_sim.cell.lamellipodium.resolve_h5_lamellipodium",
        performance_contract=_live(
            particle_types_added=("lamel_wave", "lamel_actin"),
            n_particles_mesoscale="n_WAVE + runtime-grown actin",
            n_particles_native="(scaled)",
            n_bonds="lamel_wave_anchor + lamel_mother + lamel_arp_branch + lamel_barbed_load",
            n_angles="arp branch angles",
            per_step_force=True,             # WaveMembranePin md.force.Custom
            per_batch_updater=True,          # 3 updaters: elongation, branching, capping
            uses_cpu_local_snapshot=True,
            uses_broad_phase=True,
            expected_on_recipes=("adherent_active_spread",),
            hot_path_priority=HotPathPriority.P1,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="3 batch updaters + runtime particle growth on cpu_local_snapshot.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="3 updaters host-sync; runtime topology growth"),
        # Lamellipodial bonds are explicit ACTIN structure (Arp2/3 branched
        # network); the cortical-gamma method-of-planes legitimately COUNTS them
        # (they are cortical actin, not an adhesion/junction load path). They are
        # therefore NOT denylisted. In the suspended Gate-A/B observable the
        # lamellipodium is OFF anyway.
        gamma_contaminating=False,
        requires=(),
        citations=("Bieling/Funk barbed-end kinetics", "Arp2/3 70deg"),
        pi_decisions=("Barbed-end rates inert at F=0 until membrane_load wires KU-5.x.",),
        sanity_gate_ref="ffn_sim/tests/test_lamellipodium.py",
        notes="Geometry dispatch: flat_plane | basal_ring (default single-cell) | polarized_patch (migration).",
    ),
    CompartmentSpec(
        name="membrane_load",
        category="membrane",
        status=CompartmentStatus.LIVE,
        enabled_default=False,
        summary="KU-5.2 Brownian-ratchet membrane brake loading lamellipodial barbed ends.",
        manifest_path=("optional_subsystems", "membrane_load"),
        resolve_ref="ffn_sim.cell.membrane.resolve_membrane",
        performance_contract=_live(
            particle_types_added=(),
            n_particles_mesoscale=0,
            n_particles_native=0,
            n_bonds=0,
            n_angles=0,
            per_step_force=False,
            per_batch_updater=True,          # MembraneLoadUpdater
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,
            expected_on_recipes=("adherent_active_spread",),
            hot_path_priority=HotPathPriority.P1,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Batch updater on barbed-end positions.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="updater host-sync"),
        gamma_contaminating=False,
        requires=("lamellipodium",),
        citations=("KU-5.2 PI-APPROVED 2026-06-06",),
        sanity_gate_ref="ffn_sim/tests/test_membrane.py",
    ),
    CompartmentSpec(
        name="erm",
        category="membrane",
        status=CompartmentStatus.LIVE,
        enabled_default=False,
        summary="Ezrin/radixin/moesin membrane-cortex tether: per-bead radial harmonic to R_cell (KU-3.18).",
        manifest_path=("optional_subsystems", "erm"),
        resolve_ref="ffn_sim.cortex.erm.resolve_erm",
        performance_contract=_live(
            particle_types_added=(),
            n_particles_mesoscale=0,
            n_particles_native=0,
            n_bonds=0,
            n_angles=0,
            per_step_force=True,             # ERMHarmonic md.force.Custom
            per_batch_updater=False,
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,
            expected_on_recipes=("(membrane-cortex tether studies)",),
            hot_path_priority=HotPathPriority.P1,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=True,
            bottleneck_risk="Per-step radial spring on cpu_local_snapshot.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="ERMHarmonic host-sync; native ForceCompute candidate"),
        gamma_contaminating=False,           # external radial field, not a cortical bond
        citations=("KU-3.18 k_ERM=0.1 N/m (brief anchor)",),
        pi_decisions=("ERM vs manifold soft-confinement U_conf overlap is an open PI decision (do not double-pin).",),
        sanity_gate_ref="ffn_sim/tests/test_erm.py",
    ),
    CompartmentSpec(
        name="turnover",
        category="cortex",
        status=CompartmentStatus.LIVE,
        enabled_default=False,
        summary="Cofilin-driven actin severing + re-annealing of cortex backbone bonds (k_sev=ln2/tau_half).",
        manifest_path=("optional_subsystems", "turnover"),
        resolve_ref="ffn_sim.cortex.turnover.resolve_turnover",
        performance_contract=_live(
            particle_types_added=(),
            n_particles_mesoscale=0,
            n_particles_native=0,
            n_bonds="modifies existing backbone bonds",
            n_angles=0,
            per_step_force=False,
            per_batch_updater=True,          # ActinTurnoverUpdater
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,
            expected_on_recipes=("(cortex remodelling studies)",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Batch topology edits on cpu_local_snapshot.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="topology-edit updater host-sync"),
        gamma_contaminating=False,           # edits cortex bonds in place (still counted as cortex)
        citations=("Chugh 2017 tau_half anchor",),
        sanity_gate_ref="ffn_sim/tests/test_turnover.py",
    ),
    # ----------------------- Geometry only -------------------------------
    CompartmentSpec(
        name="surface_manifold",
        category="geometry",
        status=CompartmentStatus.GEOMETRY,
        enabled_default=False,
        summary="Triangulated icosphere: local frames, signed normal distance, patch broad-phase, diagnostics.",
        manifest_path=None,
        resolve_ref="ffn_sim.cortex.surface_manifold.SurfaceManifold",
        performance_contract=_live(
            particle_types_added=(),
            n_particles_mesoscale=0,         # mesh vertices are NOT HOOMD particles
            n_particles_native=0,
            n_bonds=0,                       # NO mesh edges as bonds (mesh-as-physics forbidden)
            n_angles=0,
            per_step_force=False,            # geometry only; at most ONE soft normal confinement (default-OFF)
            per_batch_updater=False,
            uses_cpu_local_snapshot=False,
            uses_broad_phase=True,           # cKDTree bead->triangle, geodesic k-ring
            expected_on_recipes=("(broad-phase acceleration; opt-in)",),
            hot_path_priority=HotPathPriority.P3,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="NumPy/scipy geometry; not in the force budget.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="geometry refit is P3 (not hot)"),
        gamma_contaminating=False,
        geometry_only=True,
        citations=("H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX_2026-06-07.md",),
        pi_decisions=(
            "Soft normal confinement U_conf=1/2 k_conf(d-d0)^2 needs a Magic-Number Block (k_conf in [3e-6,3e-5] N/m) + PI sign-off.",
            "Mesh-resolution grid-invariance master gate must pass before any manifold force is enabled.",
        ),
        sanity_gate_ref="ffn_sim/tests/test_surface_manifold_fit.py",
        notes="HARD: no edge/area springs, no bending DOF, no hard projection, no gamma-as-mesh-property. Bright line: no manifold DOF in the gamma budget.",
    ),
    # ===================== MISSING compartments (default-OFF) ============
    # These are authored default-OFF as explicit-mechanism modules; enabling one
    # through a recipe raises UnratifiedCompartmentError pending its pairwise gate.
    CompartmentSpec(
        name="ventral_stress_fibers",
        category="adhesion",
        status=CompartmentStatus.EXPERIMENTAL,
        enabled_default=False,
        summary="Contractile ventral actomyosin bundles spanning FA to FA (alpha-actinin-banded, NMII-loaded).",
        manifest_path=None,
        resolve_ref="ffn_sim.cell.stress_fibers.resolve_stress_fibers",
        performance_contract=_live(
            particle_types_added=("sf_actin", "sf_myosin"),
            n_particles_mesoscale="~n_SF * beads_per_SF",
            n_particles_native="(scaled)",
            n_bonds="sf backbone + sf_alpha_actinin + sf_myosin bipolar",
            n_angles="sf bending (low — bundles are stiff)",
            per_step_force=False,
            per_batch_updater=True,
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,          # FA anchors are known endpoints
            expected_on_recipes=("adherent_active_spread (optional)",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Reuses myosin/xlink updater pattern; microbench at authoring.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="reuses cortex binder host-sync (P3 port)"),
        gamma_contaminating=True,
        denylist_bond_types=("sf_", "cortex_myosin_"),
        requires=("fa",),
        citations=(
            "Kojima/Gittes 1993 / Kojima 1994 (F-actin EA_single~4.3e-8 N)",
            "Cramer 1997 (vSF bundling N_filaments~10-30)",
            "Kumar 2006 single-SF ~10-30 nN (VALIDATION band only)",
            "Tojkander 2012 / Hotulainen-Lappalainen 2006 (sarcomeric periodicity)",
        ),
        pi_decisions=(
            "ACTIVATION BLOCKER: SF NMII currently reuses the cortex 'cortex_myosin_*' bond "
            "types (shared D5/D6 builder); those ARE the active-gamma signal and are NOT "
            "sf_-denylisted — a LIVE SF build would contaminate cortical gamma. A distinct "
            "'sf_myosin_*' prefix is required before wiring (or denylist cortex_myosin_* on SF builds).",
            "mu_SF = N_filaments*EA_single (DERIVED from F-actin EA, not back-solved from the "
            "Kumar tension band) but defaults to None: N_filaments has no single x40-mesoscale "
            "default, so the enabled build HALTS (NotImplementedError) until PI ratifies N_filaments. "
            "(2026-06-09 audit fix: retired the non-physical 1e-2 N/m placeholder.)",
            "Bundle bending rigidity (angle/buckling EI) not modelled — straight FA->FA chord only.",
        ),
        sanity_gate_ref="ffn_sim/cell/stress_fibers.py",
        notes="Anchors at FA clutches; bundle tension is a basal-plane observable, NOT cortical hoop gamma.",
    ),
    CompartmentSpec(
        name="linc",
        category="internal",
        status=CompartmentStatus.EXPERIMENTAL,
        enabled_default=False,
        summary="LINC complex: nesprin-SUN bridges coupling nuclear lamina to cyto cytoskeleton (cortex/SF/MT).",
        manifest_path=None,
        resolve_ref="ffn_sim.cell.linc.resolve_linc",
        performance_contract=_live(
            particle_types_added=("linc_anchor",),
            n_particles_mesoscale="~n_LINC",
            n_particles_native="(scaled)",
            n_bonds="linc_nesprin (nucleus_bead <-> cytoskeleton)",
            n_angles=0,
            per_step_force=False,
            per_batch_updater=False,         # static linker bonds; optional dynamic binding
            uses_cpu_local_snapshot=False,
            uses_broad_phase=False,
            expected_on_recipes=("confined_migration",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.HOOMD_BUILTIN,
            native_forcecompute_candidate=False,
            bottleneck_risk="Static harmonic linker bonds — HOOMD-builtin; low cost.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=True),
        gamma_contaminating=True,
        denylist_bond_types=("linc_",),
        requires=("nucleus",),
        sanity_gate_ref="ffn_sim/cell/linc.py",
    ),
    CompartmentSpec(
        name="intermediate_filaments",
        category="internal",
        status=CompartmentStatus.EXPERIMENTAL,
        enabled_default=False,
        summary="Keratin/vimentin IF network: extensible, strain-stiffening cage around the nucleus.",
        manifest_path=None,
        resolve_ref="ffn_sim.cell.intermediate_filaments.resolve_intermediate_filaments",
        performance_contract=_live(
            particle_types_added=("if_bead",),
            n_particles_mesoscale="~n_IF * beads_per_IF",
            n_particles_native="(scaled)",
            n_bonds="if backbone (nonlinear strain-stiffening) + if_crosslink",
            n_angles="if bending (Lp ~ 0.3-1 um, very flexible)",
            per_step_force=False,
            per_batch_updater=False,
            uses_cpu_local_snapshot=False,
            uses_broad_phase=False,
            expected_on_recipes=("confined_migration",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.HOOMD_BUILTIN,
            native_forcecompute_candidate=False,
            bottleneck_risk="Nonlinear bond may need a tabulated/custom potential (then per-step force).",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=True),
        gamma_contaminating=True,
        denylist_bond_types=("if_",),
        sanity_gate_ref="ffn_sim/cell/intermediate_filaments.py",
        notes="Very low persistence length + strain-stiffening (the load-bearing-at-large-strain network).",
    ),
    CompartmentSpec(
        name="microtubules",
        category="internal",
        status=CompartmentStatus.LIVE,   # wired 2026-06-09 (PI 소유권 허용); aster on shared forces, mt_ γ-denylisted
        enabled_default=False,
        summary="Stiff MTs radiating from the MTOC (Lp~5 mm), dynamic instability + compressive load-bearing.",
        manifest_path=("optional_subsystems", "microtubules"),
        resolve_ref="ffn_sim.cell.microtubules.resolve_microtubules",
        performance_contract=_live(
            particle_types_added=("mt_bead", "mtoc"),
            n_particles_mesoscale="~n_MT * beads_per_MT",
            n_particles_native="(scaled)",
            n_bonds="mt backbone (stiff) + mtoc anchor",
            n_angles="mt bending (very stiff — sets the CFL)",
            per_step_force=False,
            per_batch_updater=True,          # dynamic instability (catastrophe/rescue) optional
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,
            expected_on_recipes=("confined_migration",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Very high bending stiffness => tight CFL; may dominate dt. Microbench at authoring.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="dynamic-instability updater host-sync; stiff-MT CFL"),
        gamma_contaminating=True,
        denylist_bond_types=("mt_",),
        sanity_gate_ref="ffn_sim/cell/microtubules.py",
        notes="High bending modulus is a CFL risk — flag dt impact before enabling on the full cell.",
    ),
    CompartmentSpec(
        name="osmotic_regulation",
        category="volume",
        status=CompartmentStatus.LIVE,   # wired 2026-06-09 (PI 소유권 허용); manifest+post-build attach
        enabled_default=False,
        summary="Dynamic volume regulation: time-dependent osmolarity / ion-pump+aquaporin flux on enclosed volume.",
        manifest_path=("optional_subsystems", "osmotic_regulation"),
        resolve_ref="ffn_sim.cortex.osmotic_regulation.resolve_osmotic_regulation",
        performance_contract=_live(
            particle_types_added=(),
            n_particles_mesoscale=0,
            n_particles_native=0,
            n_bonds=0,
            n_angles=0,
            per_step_force=False,            # modulates the enclosed_volume setpoint, not a new force
            per_batch_updater=True,          # updates V0/Pi_0 over time
            uses_cpu_local_snapshot=False,
            uses_broad_phase=False,
            expected_on_recipes=("(volume-regulation / RVD-RVI studies)",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Cheap scalar updater modulating the existing turgor force setpoint.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="scalar setpoint updater"),
        gamma_contaminating=False,
        requires=("enclosed_volume",),
        sanity_gate_ref="ffn_sim/cortex/osmotic_regulation.py",
        notes="Extends the static turgor (enclosed_volume) with a dynamic V0(t)/Pi(t); same force compute, time-varying setpoint.",
    ),
    CompartmentSpec(
        name="membrane_reservoir",
        category="membrane",
        status=CompartmentStatus.EXPERIMENTAL,
        enabled_default=False,
        summary="Membrane excess-area reservoir + bleb machinery: cortex-membrane detachment & area buffering.",
        manifest_path=None,
        resolve_ref="ffn_sim.cell.membrane_reservoir.resolve_membrane_reservoir",
        performance_contract=_live(
            particle_types_added=("mem_node",),
            n_particles_mesoscale="(reservoir nodes; estimate at authoring)",
            n_particles_native="(scaled)",
            n_bonds="mem-cortex tether (breakable) + reservoir area buffer",
            n_angles=0,
            per_step_force=True,
            per_batch_updater=True,          # detachment/bleb nucleation kinetics
            uses_cpu_local_snapshot=True,
            uses_broad_phase=False,
            expected_on_recipes=("(blebbing / confinement studies)",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=True,
            bottleneck_risk="Couples to membrane_surface + ERM; tether-break updater host-sync.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="tether-break/bleb updater host-sync"),
        gamma_contaminating=True,
        denylist_bond_types=("mem_",),       # module GAMMA_DENYLIST_PREFIX='mem_' (mem_tether, ...)
        requires=("membrane_surface",),
        pi_decisions=("Bleb nucleation threshold + reservoir excess-area fraction need literature anchors or PI call.",),
        sanity_gate_ref="ffn_sim/cell/membrane_reservoir.py",
    ),
    CompartmentSpec(
        name="cadherin_junction",
        category="multicell",
        status=CompartmentStatus.EXPERIMENTAL,
        enabled_default=False,
        summary="Explicit E-cadherin trans-dimer catch-bond cell-cell junction (Rakshit sliding-rebinding).",
        manifest_path=None,
        resolve_ref="ffn_sim.junction.cadherin.resolve_cadherin_junction",
        performance_contract=_live(
            particle_types_added=("cadherin",),
            n_particles_mesoscale="~n_cad per cell-cell interface",
            n_particles_native="(scaled)",
            n_bonds="cadherin_trans_b{i} (dynamic catch-bond, r0-binned)",
            n_angles=0,
            per_step_force=False,
            per_batch_updater=True,          # CadherinBondUpdater (catch-slip bind/unbind)
            uses_cpu_local_snapshot=True,
            uses_broad_phase=True,           # cross-cell partner search
            expected_on_recipes=("multicell_junction",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Cross-cell partner search + catch-slip updater; adapts spheroid CadherinBondUpdater.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="cross-cell binder host-sync"),
        gamma_contaminating=True,
        denylist_bond_types=("cadherin_",),
        citations=(
            "Rakshit 2012 X-dimer sliding-rebinding (validation/cadherin_sliding_rebinding.py)",
            "spheroid/cadherin_bonds.CadherinBondUpdater (existing runtime catch-bond)",
        ),
        pi_decisions=("No slip-only shortcut, no lumped line tension (full KU-4.2 catch-bond per architectural rule).",),
        sanity_gate_ref="ffn_sim/tests/test_cadherin_catch_bond.py",
        notes="Adapts the existing single-cell-aggregate cadherin catch-bond to an explicit two-cell interface.",
    ),
    CompartmentSpec(
        name="junctional_actin",
        category="multicell",
        status=CompartmentStatus.STUB,
        enabled_default=False,
        summary="Junctional actin belt coupling the cadherin junction to each cell's cortex (alpha-catenin/vinculin).",
        manifest_path=None,
        resolve_ref="ffn_sim.junction.junctional_actin.resolve_junctional_actin",
        performance_contract=_live(
            particle_types_added=("junc_actin",),
            n_particles_mesoscale="(belt beads; estimate at authoring)",
            n_particles_native="(scaled)",
            n_bonds="cadherin <-> cortex actin coupling (alpha-catenin clutch)",
            n_angles=0,
            per_step_force=False,
            per_batch_updater=True,
            uses_cpu_local_snapshot=True,
            uses_broad_phase=True,
            expected_on_recipes=("multicell_junction",),
            hot_path_priority=HotPathPriority.P2,
            gpu_path_now=GpuPath.CPU,
            native_forcecompute_candidate=False,
            bottleneck_risk="Couples cadherin to cortex; reuses clutch pattern.",
        ),
        gpu_readiness=GpuReadiness(device_resident_now=False, optimization_debt="coupling clutch host-sync"),
        gamma_contaminating=True,
        denylist_bond_types=("junc_actin",),
        requires=("cadherin_junction",),
        pi_decisions=(
            "alpha-catenin force-dependent vinculin recruitment (catch) constants need a literature anchor or PI call.",
        ),
        sanity_gate_ref="ffn_sim/junction/junctional_actin.py",
        notes="STUB: the cadherin<->cortex mechanical coupling constants are not yet anchored — keep disabled pending PI.",
    ),
)


# ===========================================================================
# Registry
# ===========================================================================
class CompartmentRegistry:
    """Ordered catalogue of every compartment + recipe composition logic."""

    def __init__(self, specs: tuple[CompartmentSpec, ...] = _SPECS) -> None:
        self._specs: dict[str, CompartmentSpec] = {}
        for s in specs:
            if s.name in self._specs:
                raise ValueError(f"duplicate compartment name: {s.name}")
            self._specs[s.name] = s

    # ---- introspection ----
    def names(self) -> tuple[str, ...]:
        return tuple(self._specs.keys())

    def all(self) -> tuple[CompartmentSpec, ...]:
        return tuple(self._specs.values())

    def get(self, name: str) -> CompartmentSpec:
        try:
            return self._specs[name]
        except KeyError as exc:
            raise KeyError(
                f"unknown compartment '{name}'. Known: {', '.join(self.names())}"
            ) from exc

    def by_category(self, category: str) -> tuple[CompartmentSpec, ...]:
        return tuple(s for s in self._specs.values() if s.category == category)

    def by_status(self, status: CompartmentStatus) -> tuple[CompartmentSpec, ...]:
        return tuple(s for s in self._specs.values() if s.status == status)

    def baseline_required(self) -> tuple[str, ...]:
        return tuple(s.name for s in self._specs.values() if s.baseline_required)

    def core(self) -> tuple[str, ...]:
        return tuple(s.name for s in self._specs.values() if s.status is CompartmentStatus.CORE)

    def default_on(self) -> tuple[str, ...]:
        return tuple(s.name for s in self._specs.values() if s.enabled_default)

    # ---- gamma-contamination contract ----
    def gamma_denylist(self) -> frozenset[str]:
        """Union of bond types/prefixes that must be excluded from cortical gamma.

        A test cross-checks this against ``cortical_tension.ADHESION_BOND_TYPES``
        and ``ADHESION_BOND_TYPE_PREFIXES`` so that any compartment declaring it
        contaminates gamma also has its bond types denylisted in the estimator.
        """
        deny: set[str] = set()
        for s in self._specs.values():
            if s.gamma_contaminating:
                deny.update(s.denylist_bond_types)
        return frozenset(deny)

    # ---- invariants (asserted by tests) ----
    def validate_defaults(self) -> None:
        """Every non-core / non-baseline compartment must be default-OFF."""
        for s in self._specs.values():
            allowed_on = s.status is CompartmentStatus.CORE or s.baseline_required
            if s.enabled_default and not allowed_on:
                raise AssertionError(
                    f"compartment '{s.name}' is enabled_default but is neither core "
                    f"nor baseline_required (violates default-OFF rule)."
                )
            if not s.enabled_default and (s.status is CompartmentStatus.CORE or s.baseline_required):
                raise AssertionError(
                    f"compartment '{s.name}' is core/baseline_required but default-OFF."
                )

    def validate_manifold_geometry_only(self) -> None:
        """Geometry compartments must add zero force-bearing bonds/particles."""
        for s in self._specs.values():
            if s.geometry_only:
                pc = s.performance_contract
                if pc.n_bonds not in (0, "0") or pc.particle_types_added:
                    raise AssertionError(
                        f"geometry-only compartment '{s.name}' declares bonds/particles "
                        f"(mesh-as-physics is forbidden)."
                    )

    def validate_gamma_contract(self) -> None:
        """Any gamma-contaminating compartment must declare its denylisted bonds."""
        for s in self._specs.values():
            if s.gamma_contaminating and not s.denylist_bond_types:
                raise AssertionError(
                    f"compartment '{s.name}' is gamma_contaminating but declares no "
                    f"denylist_bond_types."
                )

    # ---- recipe composition ----
    def compose_manifest(
        self,
        recipe: dict,
        *,
        base_manifest: dict,
        strict: bool = True,
    ) -> tuple[dict, list[str]]:
        """Translate a recipe into a runnable manifest + a list of deferred names.

        Parameters
        ----------
        recipe : dict
            Resolved recipe (see :func:`load_recipe`): has ``enable`` (set of
            compartment names) and optional ``observables``.
        base_manifest : dict
            A ``mcf7_baseline.yaml``-shaped manifest to clone and toggle.
        strict : bool
            If True (default), a recipe that *requires* an EXPERIMENTAL/STUB
            compartment raises :class:`UnratifiedCompartmentError`. If False, such
            compartments are returned in ``deferred`` and left OFF in the manifest.

        Returns
        -------
        (manifest, deferred)
            ``manifest`` is a deep copy of ``base_manifest`` with LIVE compartment
            ``enabled`` flags set to match the recipe. ``deferred`` lists the
            EXPERIMENTAL/STUB compartments the recipe wanted but could not wire.
        """
        enable = set(recipe.get("enable", ()))
        unknown = enable - set(self._specs)
        if unknown:
            raise KeyError(f"recipe enables unknown compartment(s): {sorted(unknown)}")

        manifest = deepcopy(base_manifest)
        deferred: list[str] = []

        # Baseline-required compartments may never be silently dropped.
        for name in self.baseline_required():
            if name not in enable:
                # nucleus has the single sanctioned exception, surfaced explicitly.
                hint = (
                    " (the only sanctioned drop is the no-nucleus cortical-tension "
                    "build via allow_no_nucleus)" if name == "nucleus" else ""
                )
                raise BaselineDropError(
                    f"recipe '{recipe.get('name', '?')}' omits baseline-required "
                    f"compartment '{name}'{hint}."
                )

        for name in enable:
            spec = self._specs[name]
            # Dependency check.
            for dep in spec.requires:
                if dep not in enable:
                    raise ValueError(
                        f"compartment '{name}' requires '{dep}', which the recipe "
                        f"'{recipe.get('name', '?')}' does not enable."
                    )
            if spec.status in (CompartmentStatus.CORE,):
                continue  # always assembled; no manifest toggle
            if spec.is_togglable:
                _set_path(manifest, spec.manifest_path, "enabled", True)
                continue
            if spec.status in (CompartmentStatus.LIVE, CompartmentStatus.GEOMETRY):
                # LIVE-but-no-manifest-slot (rigid_ligand_coating) or geometry —
                # implied/handled outside the manifest toggles; nothing to set.
                continue
            # EXPERIMENTAL / STUB: cannot be wired through the loader yet.
            if strict:
                raise UnratifiedCompartmentError(
                    f"recipe '{recipe.get('name', '?')}' enables compartment "
                    f"'{name}' (status={spec.status.value}), which is not yet wired "
                    f"into Cell.build. Open PI decisions: "
                    f"{'; '.join(spec.pi_decisions) or 'pairwise gate pending'}."
                )
            deferred.append(name)

        # Also turn OFF any togglable LIVE compartment NOT in the recipe so the
        # manifest reflects exactly the recipe (default-OFF preserved).
        for spec in self._specs.values():
            if spec.is_togglable and spec.name not in enable:
                _set_path(manifest, spec.manifest_path, "enabled", False)

        return manifest, deferred


def _set_path(d: dict, path: tuple[str, ...], leaf_key: str, value: Any) -> None:
    """Set ``d[path...][leaf_key] = value``, creating intermediate dicts."""
    node = d
    for key in path:
        nxt = node.get(key)
        if not isinstance(nxt, dict):
            nxt = {}
            node[key] = nxt
        node = nxt
    node[leaf_key] = value


# ---------------------------------------------------------------------------
# Recipe loading
# ---------------------------------------------------------------------------
def load_recipe(path_or_name: str | Path, *, _seen: set[str] | None = None) -> dict:
    """Load a recipe YAML, resolving an optional ``extends`` parent chain.

    A recipe has: ``name``, ``description``, optional ``extends`` (parent recipe
    name), ``enable`` (list of compartment names), optional ``disable`` (removed
    from the inherited enable set), ``observables`` (diagnostic targets), ``scope``,
    ``notes``. The resolved recipe's ``enable`` is
    ``parent.enable - disable + enable``.
    """
    path = Path(path_or_name)
    if not path.is_absolute() and not path.exists():
        candidate = _RECIPE_DIR / path_or_name
        if candidate.exists():
            path = candidate
        else:
            path = _RECIPE_DIR / f"{path_or_name}.yaml"
    with open(path) as handle:
        recipe = yaml.safe_load(handle)

    seen = _seen or set()
    name = recipe.get("name", path.stem)
    if name in seen:
        raise ValueError(f"recipe extends cycle detected at '{name}'")
    seen.add(name)

    enable: list[str] = list(recipe.get("enable", []) or [])
    parent_name = recipe.get("extends")
    if parent_name:
        parent = load_recipe(parent_name, _seen=seen)
        merged = list(parent.get("enable", []))
        for n in enable:
            if n not in merged:
                merged.append(n)
        enable = merged
    for n in recipe.get("disable", []) or []:
        if n in enable:
            enable.remove(n)

    recipe["enable"] = enable
    recipe.setdefault("name", name)
    return recipe


def list_recipes() -> tuple[str, ...]:
    """Names of the recipe YAMLs available under ``configs/recipes/``."""
    if not _RECIPE_DIR.exists():
        return ()
    return tuple(sorted(p.stem for p in _RECIPE_DIR.glob("*.yaml")))


# Module-level default registry instance.
REGISTRY = CompartmentRegistry()
