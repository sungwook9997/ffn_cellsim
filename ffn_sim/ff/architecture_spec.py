"""Unified actin-architecture spec (task a, increment 1) — one weaving descriptor per structure.

PI 2026-06-04 ([[project-unified-actin-architecture]]): cortex / lamellipodium / filopodium /
stress-fibers / microvilli / traction are ONE category — actin structures that differ only in how the
SAME primitives (filaments + nucleator + crosslinker + motor) are WOVEN onto a manifold. This module is
that table IN CODE: an :class:`ArchitectureSpec` is one row; :func:`ff.weave.weave` builds it.

SCOPE (increment 1, PI-gated): only CORTEX and FILOPODIUM are defined here — the smallest pair that
PROVES the unification (an isotropic shell network AND a parallel bundle from ONE builder, no
structure-specific branch). The remaining rows (lamellipodium / stress-fiber / microvilli / traction)
+ their unsourced bundler kinetics (fascin, espin, …) are a PI sign-off item — see
`references/SE_REGISTRATION_CANDIDATES_2026-06-30.md` §5. Every constant here is either lit-anchored or
flagged; none is tuned to an outcome (CLAUDE.md hard rule).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ffn_sim.ff.hand_kmc import ALPHA_ACTININ, FILAMIN, NMIIA_MYOSIN, HandParams

# Arp2/3 branch geometry — lit-anchored, Magic-Number-Blocked (configs/phase1_h5.yaml, audit C5
# 2026-05-30; Fäßler et al. 2020 EMBO J 39:e104254 in-cell cryo-ET branch junction 68 ± 9°).
ARP23_BRANCH_ANGLE_RAD = 1.2217304764     # rad = 70.0° rest branch angle θ₀
ARP23_BRANCH_SIGMA_DEG = 9.0              # thermal angular spread σ_θ (Fäßler 2020)
ARP23_BRANCH_K = 0.173                    # pN·µm/rad² = 1.73e-19 J/rad² = kT/Var(θ), σ_θ=9° (equipartition)


@dataclass(frozen=True, slots=True)
class FilamentSpec:
    """How the actin filaments are nucleated + laid down."""

    nucleator: str          # "arp23" (branched) | "formin" (linear/parallel)
    length_um: float        # filament contour length L
    seg_um: float           # segment rest length ℓ₀
    n_filaments: int        # ×40-mesoscale count
    orientation: str        # "isotropic" | "parallel" | "branched_twomode"
    polarity: str = "mixed" # "mixed" | "uniform"
    # increment-2 additions (used only for orientation == "branched_twomode"; additive, default None)
    branch_angle_rad: float | None = None      # Arp2/3 mother-daughter rest angle θ₀ (70°, Fäßler 2020)
    branch_sigma_deg: float | None = None       # thermal angular width σ_θ (9°, Fäßler 2020)
    mode_axis_deg: float | None = None          # ±mode peak about the protrusion axis (35° = θ₀/2; Mueller 2017)
    branch_per_um: float | None = None          # linear Arp2/3 branch density (1.25/µm; Vinzenz 2012, PI-gated cite)
    # increment-3 additions (stress-fiber / microvillus bundles; additive, default None)
    sarcomere_um: float | None = None           # SF α-actinin Z-body / NMIIA band period (PI-gated swept; ~1 µm)
    bundle_spacing_um: float | None = None      # explicit inter-filament c2c for a bundle (SF/MV ~0.012; PI-gated)


@dataclass(frozen=True, slots=True)
class CrosslinkerSpec:
    """The crosslinker Hand + how it binds (the geometric weaving)."""

    hand: HandParams        # lit-anchored Hand preset (hand_kmc)
    bind_mode: str          # "any" (isotropic, cortex) | "parallel" (bundle) | "perp"
    density_per_fil: float  # n_crosslinks / n_filaments
    # ⚠️ fascin/espin/fimbrin bundlers are NOT yet sourced (PI-gated); increment 1 uses α-actinin /
    # filamin (already anchored) as the bundle crosslinker stand-in, FLAGGED for the filopodium row.


@dataclass(frozen=True, slots=True)
class MotorSpec:
    """The motor Hand + its arrangement (None = no motor, e.g. filopodium/microvillus core)."""

    hand: HandParams | None  # NMIIA_MYOSIN | None
    mode: str = "none"       # "none" | "bipolar" | "sarcomeric"
    density_per_fil: float = 0.0


@dataclass(frozen=True, slots=True)
class ArchitectureSpec:
    """One actin-structure weaving = manifold + filament + crosslinker + motor. ``weave(spec)`` builds it."""

    name: str
    manifold: str                       # "sphere" | "bundle" (increment 1); later "patch"|"finger"|"plane"
    R_um: float                         # manifold scale (sphere radius / bundle length)
    filament: FilamentSpec
    crosslinker: CrosslinkerSpec
    motor: MotorSpec = field(default_factory=lambda: MotorSpec(hand=None))
    notes: str = ""


# ── The table (increment 1: CORTEX + FILOPODIUM) ────────────────────────────────────────────────

# CORTEX — isotropic cross-linked shell on a sphere. Mirrors the H.3 production cortex
# (gamma_floor.build_crosslinked_cortex / CortexParams) so weave(CORTEX) reproduces the γ-floor.
CORTEX = ArchitectureSpec(
    name="cortex",
    manifold="sphere",
    R_um=7.5,                                                   # MCF7 radius (Wagner 2011; was 10µm mislabeled)
    filament=FilamentSpec(nucleator="formin", length_um=3.0, seg_um=0.5, n_filaments=21375,
                          orientation="isotropic", polarity="mixed"),   # NATIVE at R=7.5 (×40 retired)
    crosslinker=CrosslinkerSpec(hand=ALPHA_ACTININ, bind_mode="any", density_per_fil=1.0),
    motor=MotorSpec(hand=NMIIA_MYOSIN, mode="bipolar", density_per_fil=0.1),
    notes="×40 mesoscale cortical shell; the γ-floor cell. α-actinin/filamin + NMIIA all lit-anchored.",
)

# FILOPODIUM — tight PARALLEL formin bundle (no branching, no motor in the core). The architectural
# metrics (bundle count 10–30, parallel order ≈1, ~7–8 nm fascin spacing) are GEOMETRIC; the bundle
# crosslinker uses filamin (anchored) as a stand-in — fascin kinetics are PI-gated (SE §5).
FILOPODIUM = ArchitectureSpec(
    name="filopodium",
    manifold="bundle",
    R_um=3.0,                                                   # bundle length [µm]
    filament=FilamentSpec(nucleator="formin", length_um=3.0, seg_um=0.5, n_filaments=20,
                          orientation="parallel", polarity="uniform"),
    crosslinker=CrosslinkerSpec(hand=FILAMIN, bind_mode="parallel", density_per_fil=3.0),
    motor=MotorSpec(hand=None, mode="none"),
    notes="Parallel formin bundle (10–30 filaments, ~7–8 nm crosslink spacing). fascin = PI-gated.",
)

# LAMELLIPODIUM — Arp2/3 DENDRITIC ±35° two-mode array on a flat protrusion patch. Mothers seeded at
# ±35° about the protrusion axis; Arp2/3 daughters branch at θ₀=70° (flipping to the other ∓35° mode,
# so all filaments sit at ±35° two-mode AND every junction is at 70°). The branch angle is maintained
# by the angle-harmonic Arp2/3 kernel (network_warp.branch_angle_kernel, NOT a rigid constraint). Sparse
# filamin (the lamellipodial mesh is loosely bound vs the cortex shell); no contractile motor in the core.
LAMELLIPODIUM = ArchitectureSpec(
    name="lamellipodium",
    manifold="patch",
    R_um=8.0,                                                   # patch half-extent [µm]
    filament=FilamentSpec(nucleator="arp23", length_um=1.0, seg_um=0.5, n_filaments=200,
                          orientation="branched_twomode", polarity="uniform",
                          branch_angle_rad=ARP23_BRANCH_ANGLE_RAD,   # 70° (Fäßler 2020) [GROUNDED]
                          branch_sigma_deg=ARP23_BRANCH_SIGMA_DEG,    # σ_θ=9° (Fäßler 2020) [PI-gated SD]
                          mode_axis_deg=35.0,                         # ±35° = θ₀/2 (Mueller 2017) [GROUNDED]
                          branch_per_um=1.25),                        # Vinzenz 2012 1/0.80µm [PI-gated cite]
    crosslinker=CrosslinkerSpec(hand=FILAMIN, bind_mode="any", density_per_fil=0.3),  # SPARSE vs cortex 1.0
    motor=MotorSpec(hand=None, mode="none"),
    notes=("Arp2/3 dendritic ±35° two-mode patch; branch angle 70° (k_angle=0.173 pN·µm/rad²) via "
           "branch_angle_kernel. σ_θ=9° SD + Vinzenz branch density = PI-gated (KB/SE). fascin/espin "
           "bundlers PI-gated; uses anchored filamin as the sparse crosslinker stand-in."),
)

# STRESS-FIBER — ventral FA–FA actomyosin bundle with SARCOMERIC organization. formin antiparallel
# GRADED polarity (barbed/+ ends → the two FA ends, pointed/− → center; Hotulainen-Lappalainen 2006 —
# this rectification is REQUIRED, random mixed polarity does NOT contract); α-actinin Z-bodies PERIODIC
# at the sarcomere planes; NMIIA bands at the anti-registered band centers (Murrell 2015 alternation).
# ⚠️ sarcomere_um (~1µm) + N_filaments (7–30) are PI-GATED swept (not KB-grounded as point values).
STRESS_FIBER = ArchitectureSpec(
    name="stress_fiber",
    manifold="bundle",
    R_um=10.0,                                                  # fiber length L between two FAs [µm]
    filament=FilamentSpec(nucleator="formin", length_um=10.0, seg_um=0.5, n_filaments=20,
                          orientation="parallel", polarity="mixed",      # graded antiparallel
                          sarcomere_um=1.0,                              # PI-gated swept (0.5–1.4; Hotulainen 2006)
                          bundle_spacing_um=0.012),                       # actin packing ~12 nm (PI-gated)
    crosslinker=CrosslinkerSpec(hand=ALPHA_ACTININ, bind_mode="parallel", density_per_fil=2.0),  # Z-bodies, periodic
    motor=MotorSpec(hand=NMIIA_MYOSIN, mode="sarcomeric", density_per_fil=1.0),
    notes=("Ventral SF: FA–FA bundle, formin antiparallel graded polarity, α-actinin Z-bodies PERIODIC "
           "@sarcomere_um, NMIIA bands anti-registered. ACTIVE myosin target 5–6 nN (Kassianidou 2017); "
           "10–30 nN total EMERGES from FA anchoring, NOT back-solved. EA=4.37e-8 N (Kojima 1994). "
           "sarcomere_um + N_filaments PI-GATED (no KB point value)."),
)

# MICROVILLUS — apical brush-border finger: short uniform-polarity parallel bundle (no core motor).
# ⚠️ espin/fimbrin(I-plastin)/villin bundlers are UN-SOURCED (KB gap) → uses anchored FILAMIN as the
# stand-in (like the filopodium fascin case); ~12 nm packing + 20–30 count PI-GATED. The rootlet is a
# clamped-end BC, not a separate structure.
MICROVILLUS = ArchitectureSpec(
    name="microvillus",
    manifold="bundle",
    R_um=1.5,                                                   # core length 1–2 µm
    filament=FilamentSpec(nucleator="formin", length_um=1.5, seg_um=0.5, n_filaments=25,
                          orientation="parallel", polarity="uniform",     # all barbed-up (apex)
                          bundle_spacing_um=0.012),                        # ~12 nm hex packing (PI-gated)
    crosslinker=CrosslinkerSpec(hand=FILAMIN, bind_mode="parallel", density_per_fil=3.0),  # espin/fimbrin STAND-IN
    motor=MotorSpec(hand=None, mode="none"),
    notes=("Brush-border microvillus: finger bundle, formin parallel ~25 (20–30), uniform barbed-up. "
           "espin/fimbrin/villin UN-SOURCED → anchored FILAMIN stand-in. ~12 nm packing + count PI-GATED. "
           "(Lateral c2c ~12 nm; the actin helical repeat ~33 nm is a DIFFERENT axial period — see metrics.)"),
)

TABLE = {s.name: s for s in (CORTEX, FILOPODIUM, LAMELLIPODIUM, STRESS_FIBER, MICROVILLUS)}
