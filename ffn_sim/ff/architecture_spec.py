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


@dataclass(frozen=True, slots=True)
class FilamentSpec:
    """How the actin filaments are nucleated + laid down."""

    nucleator: str          # "arp23" (branched) | "formin" (linear/parallel)
    length_um: float        # filament contour length L
    seg_um: float           # segment rest length ℓ₀
    n_filaments: int        # ×40-mesoscale count
    orientation: str        # "isotropic" | "parallel" | "branched"
    polarity: str = "mixed" # "mixed" | "uniform"


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
    R_um=10.0,                                                  # MCF7, Wagner 2011 (KU-3.17)
    filament=FilamentSpec(nucleator="formin", length_um=3.0, seg_um=0.5, n_filaments=1000,
                          orientation="isotropic", polarity="mixed"),
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

TABLE = {s.name: s for s in (CORTEX, FILOPODIUM)}
