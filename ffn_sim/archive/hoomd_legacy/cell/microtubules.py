r"""Microtubule asters — stiff hollow tubes radiating from the MTOC (KU-MT).

ADDITIVE, DEFAULT-OFF compartment (EXTEND pattern; mirrors
``cortex/erm.py`` for the resolve+attach+Sanity-Gate shape and
``cell/nucleus.py``/``cell/lamellipodium.py`` for the pure
snapshot-extension builder). The cell otherwise has no explicit
microtubule cytoskeleton, so the compression-load-bearing / centrosomal
radial-array element of the composite mechanics is absent. This module
supplies it as an EXPLICIT fine-grained mechanism — no lumped proxy
force, no mesh-as-physics.

One-line purpose
----------------
Explicit ``mt_bead`` filament chains anchored at a single ``mtoc``
(centrosome) particle, each chain carrying a STIFF harmonic backbone bond
and a STIFF harmonic bending angle (high flexural rigidity ⇒ rod-like,
load-bearing in compression), with an OPTIONAL dynamic-instability batch
updater (growth / catastrophe / rescue, Mitchison-Kirschner 1984).

Explicit mechanism (what particles / bonds / angles)
----------------------------------------------------
* particles
    - ``mtoc``  : ONE centrosome particle at the seed point (the aster hub).
    - ``mt_bead``: the tubulin-segment beads; ``n_mt`` chains of
      ``beads_per_mt`` beads each, laid out radially from the MTOC along
      ``n_mt`` even directions (deterministic Fibonacci sphere).
* bonds (load path is NON-cortical — every bond type starts ``mt_``)
    - ``mt_backbone`` : harmonic stretch bond between consecutive beads of
      a chain (and MTOC→first-bead), stiffness ``k_backbone`` [N/m].
* angles
    - ``mt_bending``  : harmonic bending angle on every interior triple of a
      chain, rest angle ``t0 = π`` (straight rod), stiffness ``k_angle``
      [N·m/rad²] = [J/rad²]. This is the flexural-rigidity term — the whole
      mechanistic point (a microtubule is ~100× stiffer in bending than
      F-actin).

The backbone+angle harmonics are the SAME ``md.bond.Harmonic`` /
``md.angle.Harmonic`` machinery H.2 single-filament uses; only the
constants differ (microtubule EI ≫ actin EI). Compression load-bearing is
EMERGENT from the stiff angle term resisting buckling — it is NOT a lumped
"compression spring".

k_angle derivation (grid-invariant, NOT a magic number)
-------------------------------------------------------
The bending modulus (flexural rigidity) is ``EI = κ_B = k_B T · L_p``
(worm-like-chain identity; Gittes 1993). The discrete harmonic-angle
stiffness that reproduces a continuum WLC of rigidity ``EI`` on a chain of
segment length ``ℓ_0`` is the SAME bridge H.2 uses
(``ffn_sim/scripts/h2_single_filament.py``: ``angle_k = bending_modulus /
rest_length``)::

    k_angle = EI / ℓ_0                                  [N·m²]/[m] = [N·m/rad²]

Dimensions: [N·m²]/[m] = [N·m] = [J] per rad². ✓  Grid-invariant: a finer
discretisation (smaller ℓ_0, more interior angles) raises ``k_angle`` per
angle but the continuum bending response ``EI`` it reproduces is fixed —
``k_angle`` is DERIVED from ``EI`` and ``ℓ_0``, never chosen to pass a gate.
The backbone stretch stiffness is likewise ``k_backbone = Y_stretch / ℓ_0``
from the 1-D stretch modulus ``Y_stretch`` [N] (axial rigidity ``E·A``).

CFL — the AXIAL STRETCH term binds, not bending (FLAGGED, P-prominent)
---------------------------------------------------------------------
Two overdamped relaxation times are computed; the BINDING (tighter) one
sets the dt. The bending angle term gives a bending relaxation time (H.2
derivation, transverse deflection ``δ`` of an interior bead gives
``δθ ≈ δ/ℓ_0``, ``U = ½ k_angle (δ/ℓ_0)²`` so the effective
position-stiffness is ``k_angle/ℓ_0²``); the backbone stretch term gives an
axial relaxation time from the position-stiffness ``k_backbone = Y_stretch/ℓ_0``::

    τ_bend    = γ_b · ℓ_0² / k_angle = γ_b · ℓ_0³ / EI         [s]
    τ_stretch = γ_b / k_backbone     = γ_b · ℓ_0 / Y_stretch   [s]
    dt_cfl    = cfl_safety_factor · min(τ_bend, τ_stretch)     [s]

**The AXIAL STRETCH term is the tighter (binding) one at every production
``ℓ_0``, NOT bending.** The ratio is purely geometric::

    τ_bend / τ_stretch = ℓ_0² · Y_stretch / EI = (ℓ_0 / r_g)²

where ``r_g = √(EI/Y_stretch) = √(I/A)`` is the MT cross-section radius of
gyration ≈ 10 nm (the stretch↔bend crossover ``ℓ_0``). For ANY ``ℓ_0 > r_g``
the stretch DOF is stiffer (``k_backbone = Y/ℓ_0 ≫ EI/ℓ_0³``) so it sets the
tighter dt. Production ``ℓ_0`` is ~0.1–1 µm — far ABOVE the ~10 nm crossover
— so stretch binds by a large factor (e.g. ℓ_0 = 1.25 µm ⇒ ratio ≈ 1.6·10⁴,
stretch ~16,000× tighter). The EI ≈ 2.2·10⁻²³ N·m² ≈ 300× actin's
``EI_actin ≈ 7.3·10⁻²⁶ N·m²`` is a CROSS-COMPARTMENT statement (MT bend
vs actin bend) — it explains why MT bending is tighter than *other
compartments'* BEND CFLs, but bending is NOT the intra-module binding term.
The attach helper COMPUTES both, gates on ``min(τ_bend, τ_stretch)``, and
RAISES if the host ``dt`` exceeds it (``cfl_strict=True``); the resolver
stores both on the dataclass so a driver can read them BEFORE building.
**The dominant numerical risk is the STRETCH CFL, driven by the
``Y_stretch`` placeholder (PI_DECISIONS) — the correct PI action is to
ratify a ``Y_stretch`` KU (the dt lever), NOT an EI/bending-driven global
dt cut. Surface to PI before lowering the global dt to accommodate it.**

Dynamic instability (OPTIONAL batch updater — Mitchison-Kirschner 1984)
-----------------------------------------------------------------------
Microtubules stochastically switch between growth (rate ``v_g``) and shrink
(rate ``v_s``) phases, with catastrophe (growth→shrink, ``f_cat``) and
rescue (shrink→growth, ``f_res``) transitions. When ``dynamic_instability``
is enabled this is a ``hoomd.custom.Action`` (periodic ``CustomUpdater``,
the same batched-updater convention as ``cortex/crosslinkers.py`` /
``cortex/turnover.py``) that grows/shrinks each chain's plus-end bead count.
It is DEFAULT-OFF (a stable aster is the baseline); the off-rate parameters
have literature anchors below. Growth/shrink that adds/removes a bead is a
topology edit (snapshot rebuild), never an in-loop force injection.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_microtubules.py``.*

1. **Dimensional analysis**
   - ``EI = k_B T · L_p`` [J·m] = [N·m²]. ✓  ``k_angle = EI/ℓ_0`` [N·m]. ✓
   - ``k_backbone = Y_stretch/ℓ_0`` [N]/[m] = [N/m]. ✓
   - backbone ``U = ½ k_backbone Δr²`` [J]; angle ``U = ½ k_angle Δθ²`` [J]. ✓
   - ``τ_bend = γ_b ℓ_0³/EI`` → [N·s/m · m³ / (N·m²)] = [s]. ✓
2. **Boundary cases**
   - disabled (enabled=False / missing): resolver returns
     ``ResolvedMicrotubules(enabled=False, …zeros…)``; the builder
     EARLY-RETURNS the input snapshot UNCHANGED (0 particles, 0 bonds, 0
     angles added). The OFF-identity invariant.
   - ``beads_per_mt < 2`` → no backbone possible → ValueError.
   - ``beads_per_mt < 3`` → no interior bead for a bending angle → the
     builder adds backbone bonds but ZERO angles (documented; ≥3 needed for
     the flexural term). Resolver requires ≥3 when bending is on.
   - ``n_mt < 1``, ``L_mt ≤ 0``, ``L_p ≤ 0``, ``Y_stretch ≤ 0`` → ValueError.
   - ``r=0`` MTOC at origin: chains radiate outward; the MTOC is the chain
     ENDPOINT of each arm's first bending triple ``(MTOC, b0, b1)`` — so it
     participates in ``n_mt`` bending angles as an ENDPOINT but is never the
     VERTEX, hence no restoring torque acts at the hub. Because every arm's
     vertices are distinct particles, the arms bend INDEPENDENTLY (the MTOC
     injects no spurious inter-arm coupling). It also carries backbone bonds.
3. **Conservation / no-net-force**
   - ``md.bond.Harmonic`` / ``md.angle.Harmonic`` are Newton-3rd-law pair /
     triple potentials: each bond/angle injects ZERO net momentum. The aster
     as a whole therefore injects zero net force from its internal elements
     (it is an INTERNAL compartment, free to translate). STATIC: built on a
     straight rest configuration ⇒ every bond at ℓ_0, every angle at π ⇒ all
     internal forces are identically 0 (force-free construction).
   - The compartment adds NO external field (unlike ERM's lab-frame anchor).
4. **Numerical sanity (CFL — the AXIAL STRETCH term binds, FLAGGED)**
   - positions float64. BOTH timescales are COMPUTED and exposed:
     ``dt_cfl_bend = cfl_safety_factor·γ_b·ℓ_0³/EI`` and the (TIGHTER)
     ``dt_cfl_stretch = cfl_safety_factor·γ_b·ℓ_0/Y_stretch``; the binding
     gate uses ``min(bend, stretch)`` and the attach helper raises if ``dt``
     exceeds it.
   - At every production ``ℓ_0`` (> the ~10 nm crossover ``r_g = √(EI/Y_stretch)``)
     the STRETCH term is the tighter/binding one (``τ_bend/τ_stretch =
     (ℓ_0/r_g)² ≫ 1``); bending is NOT the intra-module bottleneck. The dt
     lever is therefore the ``Y_stretch`` placeholder (PI_DECISIONS), not EI.
5. **Sign / sense**
   - a STRETCHED backbone bond (``|Δr| > ℓ_0``) pulls the two beads TOWARD
     each other (restoring, ``F ∝ −(|Δr|−ℓ_0) r̂``). STATIC via HOOMD on a
     2-bead stretched rod: net force points inward along the bond.
   - a BENT chain (interior angle ``θ < π``) feels a restoring torque toward
     the straight rod ``θ = π`` (resists buckling — the compression
     load-bearing sense). STATIC: a kinked 3-bead chain's interior bead is
     pushed to straighten.
6. **Measurement-protocol consistency**
   - The flexural rigidity the chain reproduces is ``EI`` (and hence
     ``L_p = EI/k_BT``), recoverable by the H.2 tangent-correlation oracle
     (NOT re-run here — NOTED; same ``angle_k = EI/ℓ_0`` bridge that gate
     validated for actin).

Compartment Performance Contract
--------------------------------
* particle types added: ``mtoc`` (1) + ``mt_bead``.
* particle count: ``1 + n_mt · beads_per_mt``. At the project mesoscale a
  typical aster ``n_mt ≈ 20`` chains, ``beads_per_mt ≈ 25`` (25 × ℓ_0 ≈ a
  few-µm MT) ⇒ ``≈ 501`` particles. There is ONE MTOC regardless of cell
  scale; the count does NOT grow with ``n_fil`` (the ×40 actin scale): at
  ``n_fil = 1000`` and at native ``~38000`` the aster is the SAME ``≈ 501``
  particles (asters are centrosomal, not per-actin-filament). Scaling
  ``n_mt``/``beads_per_mt`` is a biology choice, not a discretisation knob.
* bond count: ``n_mt · beads_per_mt`` (``mt_backbone``: beads_per_mt per
  chain incl. the MTOC→first link). angle count:
  ``n_mt · (beads_per_mt − 1)`` interior ``mt_bending`` triples.
* per-step force: YES, but BUILT-IN HOOMD ``md.bond.Harmonic`` +
  ``md.angle.Harmonic`` (native C++/GPU ForceComputes) — NO Python
  ``md.force.Custom`` per step for the elastic terms, so no per-step
  ``cpu_local_snapshot`` for the mechanics. (Contrast ERM/nucleus, which
  ARE custom per-step forces.)
* per-batch updater: ONLY when ``dynamic_instability`` is enabled — a
  periodic ``CustomUpdater`` (snapshot rebuild) every ``mt_batch_steps``.
  Default OFF ⇒ no updater.
* uses cpu_local_snapshot: NO for the elastic mechanics (built-in forces).
  The optional DI updater reads a full ``get_snapshot()`` per batch (out of
  the hot loop) — flagged in its docstring.
* uses cKDTree / broad-phase: NO. The aster topology is fixed at build; DI
  edits only plus-end beads (no neighbour search).
* hot-path priority: **P2** (a stiff aster of ~500 particles with built-in
  forces is cheap relative to the ~38k actin LJ + binder hot path; the DI
  updater is out-of-loop).
* GPU path now: **builtin** (HOOMD ``md.bond.Harmonic`` / ``md.angle.Harmonic``
  run native on the selected device — CUDA or CPU — with NO port work). The
  DI updater is CPU snapshot-rebuild (default-off, out of loop).
* native ForceCompute candidate: N/A for the elastic terms (already native
  built-ins). The DI updater could become a native plugin if ever hot
  (unlikely — out of loop).
* bottleneck risk: **the axial STRETCH CFL, not the force cost.** The
  binding timescale is ``dt_cfl_stretch = cfl_safety_factor·γ_b·ℓ_0/Y_stretch``
  (tighter than ``dt_cfl_bend`` by ``(ℓ_0/r_g)² ≫ 1`` at production ℓ_0); the
  risk is that turning this compartment ON forces a global dt reduction
  (multiplicative slowdown), and the lever is the ``Y_stretch`` placeholder
  (PI_DECISIONS), NOT the MT EI. (The MT ``EI`` ≈ 300× actin only makes MT
  BENDING tighter than OTHER compartments' bend CFLs — a cross-compartment
  note, not the intra-module bottleneck.) FLAGGED prominently; PI-gated.

References
----------
- Gittes, Mickey, Nettleton, Howard (1993) *J. Cell Biol.* 120(4):923–934,
  "Flexural rigidity of microtubules and actin filaments measured from
  thermal fluctuations in shape." Microtubule EI = 2.2·10⁻²³ N·m²,
  L_p = 5.2 mm (at room T); actin EI = 7.3·10⁻²⁶ N·m², L_p ≈ 17.7 µm.
- Mitchison & Kirschner (1984) *Nature* 312:237–242, "Dynamic instability of
  microtubule growth" (the growth/catastrophe/rescue mechanism).
- Walker et al. (1988) *J. Cell Biol.* 107:1437–1448 — in-vitro DI rate
  constants (growth v_g ≈ 2 µm/min, shrink v_s ≈ 17 µm/min, f_cat ≈
  0.005 s⁻¹, f_res ≈ 0.044 s⁻¹) used as the OPTIONAL DI anchors.
- ``ffn_sim/scripts/h2_single_filament.py`` — the ``angle_k = bending_modulus
  / rest_length`` (EI/ℓ_0) bridge + ``τ_bend = γ_b ℓ_0³/EI`` CFL this module
  reuses (validated by the H.2 persistence-length gate).
- ``ffn_sim/cortex/erm.py`` — resolve+attach+Sanity-Gate shape, CFL attach
  gate convention (cfl_safety_factor, cfl_strict).
- ``ffn_sim/cell/lamellipodium.py`` — the snapshot extension w/ particles +
  bonds + angles pattern.
- HOOMD 7 ``md.bond.Harmonic`` / ``md.angle.Harmonic`` API.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

import hoomd
import hoomd.md as md


# ---------------------------------------------------------------------------
# Module-level declarations (gamma denylist + open PI decisions)
# ---------------------------------------------------------------------------
#: Every bond type this compartment creates starts with this prefix so the
#: cortical-tension (γ) estimator can DENYLIST the microtubule load path (it
#: is a NON-cortical compression-bearing element and must not contaminate the
#: cortex-tension measurement). See ``cortex/cortical_tension.py``.
GAMMA_DENYLIST_PREFIX: str = "mt_"

#: Open design decisions surfaced to PI (empty when none outstanding).
PI_DECISIONS: list[str] = [
    # Y_stretch: microtubule axial 1-D stretch modulus (E·A). A microtubule
    # is far stiffer in stretch than in bend, but the project has no ratified
    # KU anchor for the ×40-mesoscale MT axial rigidity. We default it to a
    # large value RELATIVE to the bending term (k_backbone ≫ k_angle/ℓ_0² so
    # the rod is inextensible on the bending timescale) and FLAG it: a precise
    # E·A (Y_M ≈ 1.2 GPa · ~190 nm² cross-section ⇒ E·A ≈ 2.3·10⁻⁷ N, Kis 2002
    # / Pampaloni 2006) needs a PI-ratified KU before production. Until then
    # the default is a documented stiff-rod placeholder, not a tuned number.
    "Y_stretch (MT axial 1-D stretch modulus E·A): no ratified KB anchor at "
    "the ×40 mesoscale; default 2.3e-7 N (E·A from Kis 2002 / Pampaloni 2006 "
    "E≈1.2 GPa, A≈190 nm²) is a stiff-rod placeholder — PI to ratify a KU.",
    # Dynamic-instability rate constants are in-vitro (Walker 1988); the in-vivo
    # / cancer-cell rates differ (MAPs, +TIPs, taxol). DI is DEFAULT-OFF and
    # the rates are flagged as in-vitro provenance until a cell-type KU lands.
    "Dynamic-instability rates (v_g, v_s, f_cat, f_res): in-vitro (Walker "
    "1988); in-vivo/cancer rates differ — DI stays default-OFF until a "
    "cell-type-specific KU is ratified.",
]


# ---------------------------------------------------------------------------
# Literature anchors (Magic-Number Block) — every constant cited
# ---------------------------------------------------------------------------
# Gittes 1993 JCB 120:923 — microtubule flexural rigidity & persistence length.
_EI_MT_DEFAULT: float = 2.2e-23        # N·m²   Gittes 1993 (EI microtubule)
_LP_MT_DEFAULT: float = 5.2e-3         # m      Gittes 1993 (L_p = 5.2 mm)
_KT_DEFAULT: float = 4.28e-21          # J      project k_B·T (310 K)
# Consistency: EI = kT·L_p ⇒ 4.28e-21 · 5.2e-3 = 2.23e-23 N·m² ≈ Gittes EI. ✓
# Kis 2002 PRL 89:248101 / Pampaloni 2006 PNAS 103:10248 — MT Young's modulus
# E ≈ 1.2 GPa, hollow-tube cross-section A ≈ 190 nm² ⇒ E·A ≈ 2.3e-7 N (PLACEHOLDER,
# see PI_DECISIONS).
_Y_STRETCH_DEFAULT: float = 2.3e-7     # N      E·A axial rigidity (PLACEHOLDER)
# Walker 1988 JCB 107:1437 — in-vitro dynamic-instability rates (DI default-OFF).
_V_GROW_DEFAULT: float = 2.0e-6 / 60.0     # m/s   v_g ≈ 2 µm/min
_V_SHRINK_DEFAULT: float = 17.0e-6 / 60.0  # m/s   v_s ≈ 17 µm/min
_F_CAT_DEFAULT: float = 0.005          # 1/s    catastrophe frequency
_F_RES_DEFAULT: float = 0.044          # 1/s    rescue frequency

# Literature bands (range-check guards — reject silent tuned overrides).
_LP_MT_BAND_M = (1.0e-3, 8.0e-3)       # Gittes/Howard MT L_p ≈ 1–8 mm band
#                                        (Gittes 5.2 mm default; Pampaloni 2006
#                                        length-dependent up to ~6–8 mm)


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class ResolvedMicrotubules:
    """Resolved microtubule-aster parameters (all SI).

    Stiffnesses are DERIVED from the continuum anchors via the grid-invariant
    bridges (``k_angle = EI/ℓ_0``, ``k_backbone = Y_stretch/ℓ_0``); the raw
    anchors are retained as provenance.

    A disabled compartment is ``ResolvedMicrotubules(enabled=False, …)`` with
    every numeric field zeroed — the OFF-identity baseline.
    """

    enabled: bool                      # —     compartment master switch

    # Topology
    n_mt: int                          # —     number of MT chains (aster arms)
    beads_per_mt: int                  # —     beads per chain (≥2; ≥3 for bend)
    L_mt: float                        # m     contour length per chain
    l0: float                          # m     segment (rest) length = L_mt/(beads_per_mt-1)
    mtoc_center: tuple[float, float, float]  # m  aster hub position

    # Continuum anchors (Magic-Number Block)
    EI: float                          # N·m²  flexural rigidity (Gittes 1993)
    L_p: float                         # m     persistence length (= EI/kT)
    Y_stretch: float                   # N     axial 1-D stretch modulus (E·A)
    kT: float                          # J     k_B·T

    # Derived discrete stiffnesses (grid-invariant bridges)
    k_angle: float                     # N·m/rad²  = EI/ℓ_0   (bending)
    k_backbone: float                  # N/m       = Y_stretch/ℓ_0 (stretch)
    angle_t0: float                    # rad       rest bending angle (= π)

    # Derived CFL timescales (γ_b supplied at resolve; STRETCH is the tight one
    # at production ℓ_0 — τ_bend/τ_stretch = (ℓ_0/r_g)² ≫ 1, r_g≈10 nm crossover)
    gamma_b: float                     # N·s/m  per-bead Stokes drag
    tau_bend: float                    # s      γ_b·ℓ_0³/EI
    tau_stretch: float                 # s      γ_b·ℓ_0/Y_stretch (STIFF — flagged)
    dt_cfl_bend: float                 # s      cfl_safety·τ_bend
    dt_cfl_stretch: float              # s      cfl_safety·τ_stretch (dominant risk)
    cfl_safety_factor: float           # —

    # Optional dynamic instability (Mitchison-Kirschner 1984 / Walker 1988)
    dynamic_instability: bool          # —     DI updater on? (default False)
    v_grow: float                      # m/s   growth velocity (Walker 1988)
    v_shrink: float                    # m/s   shrink velocity (Walker 1988)
    f_cat: float                       # 1/s   catastrophe frequency
    f_res: float                       # 1/s   rescue frequency

    # Bond / angle type names (γ denylist prefix-tagged)
    backbone_bond_type: str = "mt_backbone"
    bending_angle_type: str = "mt_bending"
    mtoc_type: str = "mtoc"
    bead_type: str = "mt_bead"

    @property
    def n_beads_total(self) -> int:
        """Particle count added: 1 MTOC + n_mt·beads_per_mt mt_beads."""
        if not self.enabled:
            return 0
        return 1 + self.n_mt * self.beads_per_mt

    @property
    def n_backbone_bonds(self) -> int:
        """``mt_backbone`` bond count (beads_per_mt per chain incl. MTOC link)."""
        if not self.enabled:
            return 0
        return self.n_mt * self.beads_per_mt

    @property
    def n_bending_angles(self) -> int:
        """``mt_bending`` interior-triple count (0 if beads_per_mt < 3)."""
        if not self.enabled or self.beads_per_mt < 3:
            return 0
        # Each chain: MTOC + beads_per_mt beads = (beads_per_mt+1) nodes in a
        # path ⇒ (beads_per_mt+1) - 2 = beads_per_mt - 1 interior angles.
        return self.n_mt * (self.beads_per_mt - 1)

    @property
    def dt_cfl(self) -> float:
        """The binding (tighter) CFL dt = min(bend, stretch) [s]."""
        if not self.enabled:
            return float("inf")
        return min(self.dt_cfl_bend, self.dt_cfl_stretch)


def _disabled() -> ResolvedMicrotubules:
    """OFF-identity baseline: every numeric field zeroed, ``enabled=False``."""
    return ResolvedMicrotubules(
        enabled=False,
        n_mt=0,
        beads_per_mt=0,
        L_mt=0.0,
        l0=0.0,
        mtoc_center=(0.0, 0.0, 0.0),
        EI=0.0,
        L_p=0.0,
        Y_stretch=0.0,
        kT=0.0,
        k_angle=0.0,
        k_backbone=0.0,
        angle_t0=math.pi,
        gamma_b=0.0,
        tau_bend=0.0,
        tau_stretch=0.0,
        dt_cfl_bend=0.0,
        dt_cfl_stretch=0.0,
        cfl_safety_factor=0.0,
        dynamic_instability=False,
        v_grow=0.0,
        v_shrink=0.0,
        f_cat=0.0,
        f_res=0.0,
    )


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def resolve_microtubules(
    cfg: dict,
    *,
    kT: float = _KT_DEFAULT,
    gamma_b: float,
    cfl_safety_factor: float = 0.1,
    mtoc_center: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> ResolvedMicrotubules:
    r"""Resolve the ``microtubules`` config block into SI aster parameters.

    DEFAULT-OFF: if ``cfg`` is missing the block, or it has
    ``enabled: false`` (or no ``enabled`` key ⇒ treated as False), this
    returns :func:`_disabled` — ``ResolvedMicrotubules(enabled=False, …zeros…)``
    — and NOTHING downstream builds anything.

    When enabled, the bending stiffness is DERIVED from the flexural rigidity
    ``EI`` via the grid-invariant bridge ``k_angle = EI/ℓ_0`` (the same bridge
    H.2 validated), the backbone stiffness from ``k_backbone = Y_stretch/ℓ_0``,
    and the (tight) bending CFL ``dt_cfl_bend = cfl_safety_factor·γ_b·ℓ_0³/EI``
    is computed and stored for the driver to read BEFORE building.

    Args:
        cfg: Config mapping. Accepts the root, the ``cell`` sub-dict, or the
            ``cell.microtubules`` / ``microtubules`` sub-dict. Recognised keys
            (all optional except as noted): ``enabled`` (bool, default False),
            ``n_mt`` (int), ``beads_per_mt`` (int ≥3 for bending), ``L_mt`` [m]
            (REQUIRED when enabled — host MT contour length), ``EI`` [N·m²],
            ``L_p`` [m], ``Y_stretch`` [N], ``dynamic_instability`` (bool),
            ``v_grow`` / ``v_shrink`` [m/s], ``f_cat`` / ``f_res`` [1/s].
        kT: Thermal energy [J] (default project ``4.28e-21``). Used for the
            ``EI ↔ L_p`` consistency check.
        gamma_b: Per-bead Stokes drag [N·s/m] (= 6π η R_bead). REQUIRED — the
            host supplies it; needed for the bending CFL. (No default — a
            wrong drag silently mis-sizes the CFL.)
        cfl_safety_factor: CFL safety factor (D3 BAOAB convention). Default 0.1.
        mtoc_center: Aster hub position [m] (host cell centroid). Default origin.

    Returns:
        Resolved microtubule parameters (SI). ``enabled=False`` ⇒ all zeros.

    Raises:
        ValueError: on a non-positive / out-of-band parameter when enabled
            (Sanity Gate §2).
    """
    # Unwrap nested config locations.
    if "cell" in cfg and isinstance(cfg["cell"], dict):
        cfg = cfg["cell"]
    if "microtubules" in cfg and isinstance(cfg["microtubules"], dict):
        cfg = cfg["microtubules"]

    enabled = bool(cfg.get("enabled", False))
    if not enabled:
        return _disabled()

    # γ_b is always needed for the CFL (even before the topology checks).
    _require_finite_positive("gamma_b", float(gamma_b))
    _require_finite_positive("kT", float(kT))
    if not (math.isfinite(cfl_safety_factor) and 0.0 < cfl_safety_factor <= 1.0):
        raise ValueError(
            f"cfl_safety_factor must be in (0, 1]; got {cfl_safety_factor!r}"
        )

    n_mt = int(cfg.get("n_mt", 20))
    beads_per_mt = int(cfg.get("beads_per_mt", 25))
    if n_mt < 1:
        raise ValueError(f"n_mt must be ≥ 1 when enabled; got {n_mt!r}")
    if beads_per_mt < 2:
        raise ValueError(
            f"beads_per_mt must be ≥ 2 (need a backbone bond); got "
            f"{beads_per_mt!r}"
        )
    if beads_per_mt < 3:
        # ≥3 beads are needed for an interior bending angle — the flexural term
        # is the whole mechanistic point of a microtubule. Refuse the degenerate
        # no-bending aster rather than silently dropping the rigidity.
        raise ValueError(
            f"beads_per_mt must be ≥ 3 for a bending angle (flexural rigidity "
            f"is the microtubule's defining mechanism); got {beads_per_mt!r}"
        )

    L_mt = cfg.get("L_mt", None)
    if L_mt is None:
        raise ValueError(
            "L_mt (MT contour length [m]) is REQUIRED when microtubules are "
            "enabled — set it from the host cell geometry; do NOT invent it."
        )
    L_mt = float(L_mt)
    _require_finite_positive("L_mt", L_mt)

    EI = float(cfg.get("EI", _EI_MT_DEFAULT))
    L_p = float(cfg.get("L_p", _LP_MT_DEFAULT))
    Y_stretch = float(cfg.get("Y_stretch", _Y_STRETCH_DEFAULT))
    _require_finite_positive("EI", EI)
    _require_finite_positive("L_p", L_p)
    _require_finite_positive("Y_stretch", Y_stretch)

    # KU band guard on L_p (reject silent tuned overrides).
    lo, hi = _LP_MT_BAND_M
    if not (lo <= L_p <= hi):
        raise ValueError(
            f"L_p = {L_p:.3e} m outside the Gittes/Howard microtubule "
            f"persistence-length band [{lo:.1e}, {hi:.1e}] m "
            f"({lo * 1e3:.0f}–{hi * 1e3:.0f} mm). Per CLAUDE.md "
            "no-magic-number: keep within band or surface to PI."
        )
    # EI ↔ L_p ↔ kT consistency: EI should equal kT·L_p to ~10 % (both are
    # measured independently; flag a gross mismatch — a tuned EI would break it).
    EI_from_lp = kT * L_p
    if not (0.5 * EI_from_lp <= EI <= 2.0 * EI_from_lp):
        raise ValueError(
            f"EI = {EI:.3e} N·m² inconsistent with kT·L_p = {EI_from_lp:.3e} "
            f"N·m² (factor {EI / EI_from_lp:.2f}); the WLC identity EI = kT·L_p "
            "must hold to ~2×. One of EI / L_p / kT is mis-set — surface to PI."
        )

    # Segment (rest) length: contour length over the number of springs.
    l0 = L_mt / float(beads_per_mt - 1)
    _require_finite_positive("l0", l0)

    # Grid-invariant bridges (NOT magic numbers — derived from EI / Y_stretch).
    k_angle = EI / l0                  # N·m/rad²  (== bending_modulus/rest_length)
    k_backbone = Y_stretch / l0        # N/m

    # CFL timescales (H.2 derivation). At production ℓ_0 (> the ~10 nm
    # crossover r_g = √(EI/Y_stretch)) the axial STRETCH term is the tighter /
    # binding one: τ_bend/τ_stretch = (ℓ_0/r_g)² ≫ 1. The gate uses min(...).
    tau_bend = float(gamma_b) * l0**3 / EI
    tau_stretch = float(gamma_b) * l0 / Y_stretch
    dt_cfl_bend = cfl_safety_factor * tau_bend
    dt_cfl_stretch = cfl_safety_factor * tau_stretch

    # Optional dynamic instability.
    di = bool(cfg.get("dynamic_instability", False))
    v_grow = float(cfg.get("v_grow", _V_GROW_DEFAULT))
    v_shrink = float(cfg.get("v_shrink", _V_SHRINK_DEFAULT))
    f_cat = float(cfg.get("f_cat", _F_CAT_DEFAULT))
    f_res = float(cfg.get("f_res", _F_RES_DEFAULT))
    if di:
        for nm, val in (
            ("v_grow", v_grow), ("v_shrink", v_shrink),
            ("f_cat", f_cat), ("f_res", f_res),
        ):
            _require_finite_positive(nm, val)

    return ResolvedMicrotubules(
        enabled=True,
        n_mt=n_mt,
        beads_per_mt=beads_per_mt,
        L_mt=L_mt,
        l0=l0,
        mtoc_center=tuple(float(c) for c in mtoc_center),  # type: ignore[arg-type]
        EI=EI,
        L_p=L_p,
        Y_stretch=Y_stretch,
        kT=float(kT),
        k_angle=k_angle,
        k_backbone=k_backbone,
        angle_t0=math.pi,
        gamma_b=float(gamma_b),
        tau_bend=tau_bend,
        tau_stretch=tau_stretch,
        dt_cfl_bend=dt_cfl_bend,
        dt_cfl_stretch=dt_cfl_stretch,
        cfl_safety_factor=cfl_safety_factor,
        dynamic_instability=di,
        v_grow=v_grow,
        v_shrink=v_shrink,
        f_cat=f_cat,
        f_res=f_res,
    )


# ---------------------------------------------------------------------------
# PURE topology builder (computes & returns arrays; NEVER mutates a sim/file)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class MTTopology:
    """Computed aster topology arrays for the Lead to append to a snapshot.

    All arrays are in the LOCAL aster index space: index 0 is the MTOC, then
    ``n_mt`` chains of ``beads_per_mt`` ``mt_bead`` particles in order. The
    Lead OFFSETS these by the existing particle count before merging into the
    host snapshot (single-writer ``cell.py``).
    """

    positions: np.ndarray              # (n_beads_total, 3) float64 [m]
    type_names: list[str]              # per-particle type name
    mtoc_local_index: int              # local index of the MTOC (0)
    backbone_bonds: np.ndarray         # (n_bonds, 2) int64 local indices
    bending_angles: np.ndarray         # (n_angles, 3) int64 local indices
    backbone_bond_type: str
    bending_angle_type: str


def build_mt_topology(p: ResolvedMicrotubules) -> MTTopology:
    r"""Compute the aster topology (PURE — returns arrays, mutates nothing).

    Lays out ``n_mt`` straight chains radiating from a single MTOC along even
    Fibonacci-sphere directions, each chain ``beads_per_mt`` beads at segment
    length ``ℓ_0``. Built STRAIGHT (every bond at ℓ_0, every interior angle at
    π) so the aster is force-free at construction (Sanity Gate §3).

    Returns an empty topology (no particles/bonds/angles) when the compartment
    is disabled — so callers can build unconditionally and the OFF case is a
    no-op.

    Args:
        p: Resolved microtubule parameters.

    Returns:
        :class:`MTTopology` with local-index arrays. Index 0 is the MTOC.

    Raises:
        ValueError: if ``p`` is internally inconsistent (should not happen for
            a resolver output).
    """
    if not p.enabled:
        return MTTopology(
            positions=np.empty((0, 3), dtype=np.float64),
            type_names=[],
            mtoc_local_index=-1,
            backbone_bonds=np.empty((0, 2), dtype=np.int64),
            bending_angles=np.empty((0, 3), dtype=np.int64),
            backbone_bond_type=p.backbone_bond_type,
            bending_angle_type=p.bending_angle_type,
        )

    n_mt = p.n_mt
    npb = p.beads_per_mt
    l0 = p.l0
    c = np.asarray(p.mtoc_center, dtype=np.float64).reshape(3)

    # Even radial directions via the Fibonacci sphere (deterministic, no RNG
    # clumping) — the same scheme nucleus.build_nucleus_beads uses.
    idx = np.arange(n_mt, dtype=np.float64) + 0.5
    phi = math.pi * (3.0 - math.sqrt(5.0))            # golden angle
    cos_theta = 1.0 - 2.0 * idx / float(n_mt)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    sin_theta = np.sqrt(np.maximum(0.0, 1.0 - cos_theta * cos_theta))
    az = phi * idx
    dirs = np.stack(
        [sin_theta * np.cos(az), sin_theta * np.sin(az), cos_theta], axis=1
    )

    positions = [c.copy()]                            # local index 0 = MTOC
    type_names = [p.mtoc_type]
    backbone_bonds: list[tuple[int, int]] = []
    bending_angles: list[tuple[int, int, int]] = []

    mtoc_idx = 0
    for arm in range(n_mt):
        d = dirs[arm]
        prev = mtoc_idx                               # chain starts at the MTOC
        chain_nodes = [mtoc_idx]
        for b in range(npb):
            r = c + (b + 1) * l0 * d                  # bead b at (b+1)·ℓ_0
            cur = len(positions)
            positions.append(r)
            type_names.append(p.bead_type)
            backbone_bonds.append((prev, cur))        # mt_backbone link
            chain_nodes.append(cur)
            prev = cur
        # Interior bending angles: every consecutive triple along the path
        # (MTOC, b0, b1, …). path has (npb+1) nodes ⇒ (npb-1) interior angles.
        if npb >= 3:
            for k in range(len(chain_nodes) - 2):
                bending_angles.append(
                    (chain_nodes[k], chain_nodes[k + 1], chain_nodes[k + 2])
                )

    pos_arr = np.asarray(positions, dtype=np.float64).reshape(-1, 3)
    bb_arr = (
        np.asarray(backbone_bonds, dtype=np.int64).reshape(-1, 2)
        if backbone_bonds else np.empty((0, 2), dtype=np.int64)
    )
    ang_arr = (
        np.asarray(bending_angles, dtype=np.int64).reshape(-1, 3)
        if bending_angles else np.empty((0, 3), dtype=np.int64)
    )

    return MTTopology(
        positions=pos_arr,
        type_names=type_names,
        mtoc_local_index=0,
        backbone_bonds=bb_arr,
        bending_angles=ang_arr,
        backbone_bond_type=p.backbone_bond_type,
        bending_angle_type=p.bending_angle_type,
    )


# ---------------------------------------------------------------------------
# Snapshot-extension builder (callable in isolation; no-op when disabled)
# ---------------------------------------------------------------------------
def extend_snapshot_with_microtubules(
    snap,
    p: ResolvedMicrotubules,
):
    r"""Return a snapshot extended with the MT aster (MTOC + chains).

    DEFAULT-OFF / OFF-IDENTITY: when ``p.enabled`` is False this is a
    NO-OP — it returns the INPUT ``snap`` UNCHANGED (same object, 0 particles,
    0 bonds, 0 angles added). When enabled, it builds the aster topology
    (:func:`build_mt_topology`) and returns a FRESH ``hoomd.Snapshot`` with the
    ``mtoc`` + ``mt_bead`` particles appended, the ``mt_backbone`` bonds and
    ``mt_bending`` angles registered (γ-denylist-prefixed types), and all
    existing particles/bonds/angles preserved.

    Mirrors ``cell/lamellipodium._extend_snapshot_with_new_actins`` (particles
    + bonds + angles). Callable in ISOLATION on any gsd/HOOMD snapshot — it
    does NOT need a built Cell.

    Args:
        snap: A ``hoomd.Snapshot`` or ``gsd.hoomd.Frame`` with ``.particles``,
            ``.bonds``, ``.angles``, ``.configuration``.
        p: Resolved microtubule parameters.

    Returns:
        The unchanged input ``snap`` when disabled; otherwise a new
        ``hoomd.Snapshot`` with the aster appended.
    """
    # ---- OFF-IDENTITY early return (the input is returned untouched). ----
    if not p.enabled:
        return snap

    def _grp(arr, width: int) -> np.ndarray:
        # gsd.hoomd.Frame returns None for an empty bonds/angles group; HOOMD
        # Snapshots return a (0, width) array. Normalise both to (n, width).
        if arr is None:
            return np.empty((0, width), dtype=np.int64)
        a = np.asarray(arr, dtype=np.int64)
        return a.reshape(-1, width) if a.size else np.empty((0, width), np.int64)

    def _tid(arr) -> np.ndarray:
        if arr is None:
            return np.empty((0,), dtype=np.uint32)
        return np.asarray(arr, dtype=np.uint32).reshape(-1)

    topo = build_mt_topology(p)
    n_old = int(snap.particles.N)
    n_new = int(topo.positions.shape[0])
    n_total = n_old + n_new

    def _pf(arr, default: np.ndarray) -> np.ndarray:
        # A minimally-populated gsd.hoomd.Frame returns None for unset particle
        # fields (typeid / velocity / mass / image); HOOMD fills them on load.
        # Mirror those auto-population defaults here so the in-isolation gsd
        # path (docstring promise) does not crash on np.asarray(None).
        return default if arr is None else np.asarray(arr)

    write = hoomd.Snapshot()
    write.particles.N = n_total

    # --- particle types: ensure mtoc + mt_bead registered. ---
    types = list(snap.particles.types)
    for tn in (p.mtoc_type, p.bead_type):
        if tn not in types:
            types.append(tn)
    write.particles.types = types
    type_index = {tn: i for i, tn in enumerate(types)}

    old_typeid = _pf(
        snap.particles.typeid, np.zeros(n_old, dtype=np.uint32)
    ).astype(np.uint32).reshape(-1)
    new_typeid = np.array(
        [type_index[tn] for tn in topo.type_names], dtype=np.uint32
    )
    write.particles.typeid[:] = np.concatenate([old_typeid, new_typeid])

    old_pos = _pf(
        snap.particles.position, np.zeros((n_old, 3), dtype=np.float64)
    ).astype(np.float64).reshape(-1, 3)
    write.particles.position[:] = np.concatenate(
        [old_pos, topo.positions], axis=0
    )

    # velocity / mass / image — append zeros / ones (mirror lamellipodium);
    # HOOMD auto-population defaults are typeid=0, mass=1, velocity/image=0.
    old_vel = _pf(
        snap.particles.velocity, np.zeros((n_old, 3), dtype=np.float64)
    ).astype(np.float64).reshape(-1, 3)
    write.particles.velocity[:] = np.concatenate(
        [old_vel, np.zeros((n_new, 3), dtype=np.float64)], axis=0
    )
    old_mass = _pf(
        snap.particles.mass, np.ones(n_old, dtype=np.float64)
    ).astype(np.float64).reshape(-1)
    write.particles.mass[:] = np.concatenate(
        [old_mass, np.ones(n_new, dtype=np.float64)]
    )
    old_image = _pf(
        snap.particles.image, np.zeros((n_old, 3), dtype=np.int32)
    ).astype(np.int32).reshape(-1, 3)
    write.particles.image[:] = np.concatenate(
        [old_image, np.zeros((n_new, 3), dtype=np.int32)], axis=0
    )

    write.configuration.box = list(snap.configuration.box)

    # --- bonds: existing + mt_backbone (offset local indices by n_old). ---
    # snap.bonds.types is None on a bare gsd Frame ⇒ normalise to [].
    bond_types = list(snap.bonds.types or [])
    if topo.backbone_bond_type not in bond_types:
        bond_types.append(topo.backbone_bond_type)
    bb_typeid = bond_types.index(topo.backbone_bond_type)

    old_bg = _grp(snap.bonds.group, 2)
    old_bt = _tid(snap.bonds.typeid)
    new_bg = topo.backbone_bonds + n_old              # local → global offset
    new_bt = np.full(new_bg.shape[0], bb_typeid, dtype=np.uint32)
    merged_bg = np.concatenate([old_bg, new_bg], axis=0).astype(np.uint32)
    merged_bt = np.concatenate([old_bt, new_bt]).astype(np.uint32)
    write.bonds.types = bond_types
    write.bonds.N = int(merged_bg.shape[0])
    if merged_bg.shape[0] > 0:
        write.bonds.group[:] = merged_bg
        write.bonds.typeid[:] = merged_bt

    # --- angles: existing + mt_bending. ---
    # snap.angles.types is None on a bare gsd Frame ⇒ normalise to [].
    angle_types = list(snap.angles.types or [])
    old_ag = _grp(snap.angles.group, 3)
    old_at = _tid(snap.angles.typeid)
    if topo.bending_angles.shape[0] > 0:
        if topo.bending_angle_type not in angle_types:
            angle_types.append(topo.bending_angle_type)
        ba_typeid = angle_types.index(topo.bending_angle_type)
        new_ag = topo.bending_angles + n_old
        new_at = np.full(new_ag.shape[0], ba_typeid, dtype=np.uint32)
        merged_ag = np.concatenate([old_ag, new_ag], axis=0).astype(np.uint32)
        merged_at = np.concatenate([old_at, new_at]).astype(np.uint32)
    else:
        merged_ag = old_ag.astype(np.uint32)
        merged_at = old_at.astype(np.uint32)
    if merged_ag.shape[0] > 0:
        write.angles.types = angle_types
        write.angles.N = int(merged_ag.shape[0])
        write.angles.group[:] = merged_ag
        write.angles.typeid[:] = merged_at
    elif angle_types:
        write.angles.types = angle_types
        write.angles.N = 0

    # --- pass through dihedrals / impropers if present. ---
    for grp_name in ("dihedrals", "impropers"):
        src = getattr(snap, grp_name, None)
        if src is None:
            continue
        src_n = getattr(src, "N", 0) or 0
        if int(src_n) > 0:
            dst = getattr(write, grp_name)
            dst.N = int(src.N)
            dst.types = list(src.types)
            dst.group[:] = np.asarray(src.group)
            dst.typeid[:] = np.asarray(src.typeid)

    return write


# ---------------------------------------------------------------------------
# Force attach helper (built-in Harmonic forces + the STIFF bending CFL gate)
# ---------------------------------------------------------------------------
def attach_microtubule_forces(
    sim: hoomd.Simulation,
    p: ResolvedMicrotubules,
    *,
    cfl_strict: bool = True,
) -> tuple:
    r"""Attach ``mt_backbone`` / ``mt_bending`` built-in Harmonic forces.

    DEFAULT-OFF: when ``p.enabled`` is False this is a NO-OP — it returns
    ``(None, None)`` and touches nothing (the integrator is unchanged, so a
    pre-MT run is bit-for-bit identical).

    When enabled it builds a ``md.bond.Harmonic`` for ``mt_backbone`` (rest
    length ``ℓ_0``, ``k = k_backbone``) and a ``md.angle.Harmonic`` for
    ``mt_bending`` (rest angle ``π``, ``k = k_angle``), appends them to the
    Integrator's force list, and ENFORCES the stiff bending CFL.

    **CFL — the dominant risk is the axial STRETCH term.** Both timescales are
    computed; at production ``ℓ_0`` (> the ~10 nm crossover ``r_g``) the
    binding one is ``dt_cfl_stretch = cfl_safety_factor·γ_b·ℓ_0/Y_stretch``
    (tighter than ``dt_cfl_bend`` by ``(ℓ_0/r_g)² ≫ 1``). The gate uses
    ``min(dt_cfl_bend, dt_cfl_stretch)`` and RAISES (``cfl_strict``) if
    ``dt`` exceeds it — DO NOT silently lower the global dt to accommodate this;
    surface to PI (a global dt cut is a multiplicative slowdown of every other
    compartment), and the lever is the ``Y_stretch`` placeholder
    (``PI_DECISIONS``), NOT the MT EI.

    Args:
        sim: Already-built simulation (must have an Integrator with ``dt``, and
            the ``mt_backbone`` / ``mt_bending`` types present in the state).
        p: Resolved microtubule parameters.
        cfl_strict: If True, raise ``RuntimeError`` on a bending-CFL violation;
            False only for diagnostic runs (caller's responsibility).

    Returns:
        ``(backbone_force, bending_force)`` — the attached HOOMD ForceComputes
        (``(None, None)`` when disabled). Kept for introspection in tests.

    Raises:
        RuntimeError: if no Integrator is set, or the bending CFL is violated
            (and ``cfl_strict``).
    """
    if not p.enabled:
        return (None, None)

    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "sim.operations.integrator must be set before attaching the "
            "microtubule forces."
        )

    # --- CFL gate (the flagged risk). The BINDING term is whichever of
    #     τ_bend / τ_stretch is smaller; at production ℓ_0 it is STRETCH. ---
    dt = float(ig.dt)
    dt_cfl = p.dt_cfl                                  # min(bend, stretch)
    if dt > dt_cfl and cfl_strict:
        stretch_binds = p.tau_stretch <= p.tau_bend
        binding = "AXIAL STRETCH" if stretch_binds else "BENDING"
        # r_g = √(EI/Y_stretch) is the stretch↔bend crossover ℓ_0 (≈ MT
        # radius of gyration ~10 nm); ratio = (ℓ_0/r_g)² = τ_bend/τ_stretch.
        r_g = math.sqrt(p.EI / p.Y_stretch)
        ratio = p.tau_bend / p.tau_stretch
        raise RuntimeError(
            f"Microtubule CFL violated: dt = {dt:.3e} s > "
            f"{p.cfl_safety_factor:.2f} · τ_min = {dt_cfl:.3e} s, where "
            f"τ_min = min(τ_bend = γ_b·ℓ_0³/EI = {p.tau_bend:.3e} s, "
            f"τ_stretch = γ_b·ℓ_0/Y_stretch = {p.tau_stretch:.3e} s). "
            f"The {binding} term is the BINDING one here. For any ℓ_0 > the MT "
            f"radius of gyration r_g = √(EI/Y_stretch) ≈ {r_g:.2e} m (~10 nm) "
            f"the AXIAL STRETCH term (k_backbone = Y_stretch/ℓ_0 = "
            f"{p.k_backbone:.3e} N/m) is tighter, by a factor "
            f"τ_bend/τ_stretch = (ℓ_0/r_g)² = {ratio:.0f}× at this ℓ_0. "
            f"Need dt ≤ {dt_cfl:.3e} s — {dt / dt_cfl:.1f}× smaller than the "
            "current global dt. DO NOT silently cut the global dt (it slows "
            "EVERY compartment multiplicatively) — surface to PI; the dt lever "
            "is the Y_stretch placeholder (see microtubules.PI_DECISIONS), NOT "
            "the MT EI / bending. Pass cfl_strict=False for a diagnostic run only."
        )

    backbone = md.bond.Harmonic()
    backbone.params[p.backbone_bond_type] = dict(k=p.k_backbone, r0=p.l0)
    ig.forces.append(backbone)

    bending = md.angle.Harmonic()
    bending.params[p.bending_angle_type] = dict(k=p.k_angle, t0=p.angle_t0)
    ig.forces.append(bending)

    return (backbone, bending)


# ---------------------------------------------------------------------------
# Optional dynamic-instability batch updater (Mitchison-Kirschner 1984)
# ---------------------------------------------------------------------------
class MTDynamicInstability(hoomd.custom.Action):
    r"""OPTIONAL plus-end growth/catastrophe/rescue updater (default-OFF).

    Mitchison & Kirschner (1984) dynamic instability as a batched
    ``hoomd.custom.Action`` (the same out-of-loop convention as
    ``cortex/crosslinkers.py`` / ``cortex/turnover.py``): each batch tick, each
    chain is in a growth or shrink phase; with phase-transition probabilities
    ``1 − exp(−f·Δt_batch)`` it switches phase, and grows/shrinks the plus-end
    bead count by ``round(v·Δt_batch/ℓ_0)`` beads (a topology edit via a
    snapshot rebuild — NEVER an in-loop force).

    NOT IMPLEMENTED in this minimal version: the plus-end add/remove snapshot
    rebuild is a Lead integration concern (it must coordinate with the
    single-writer ``cell.py`` topology and the BAOAB ``gamma_map``). This class
    documents the mechanism and the rate anchors; constructing it raises
    ``NotImplementedError`` so a disabled-by-default driver cannot accidentally
    run an unfinished updater. DI is gated OFF in :func:`resolve_microtubules`
    by default; enabling it requires the rebuild to be wired (PI-gated; see
    ``PI_DECISIONS``).

    Args:
        p: Resolved microtubule parameters (must have ``dynamic_instability``).
        dt: Host integrator timestep [s].
        mt_batch_steps: Steps between DI ticks (Δt_batch = mt_batch_steps·dt).
        seed: RNG seed offset.

    Raises:
        NotImplementedError: always (the plus-end rebuild is unfinished —
            honesty over completeness; DI stays default-OFF).
    """

    def __init__(
        self,
        p: ResolvedMicrotubules,
        *,
        dt: float,
        mt_batch_steps: int = 100,
        seed: int = 0,
    ) -> None:
        super().__init__()
        self.p = p
        self.dt = float(dt)
        self.mt_batch_steps = int(mt_batch_steps)
        self.dt_batch = self.dt * self.mt_batch_steps
        self.rng = np.random.default_rng(int(seed))
        raise NotImplementedError(
            "MTDynamicInstability plus-end growth/shrink snapshot-rebuild is "
            "not implemented in this minimal version (the topology edit must "
            "be coordinated with cell.py's single-writer rule + the BAOAB "
            "gamma_map). Dynamic instability is DEFAULT-OFF; wiring it is a "
            "PI-gated Lead integration step (see microtubules.PI_DECISIONS). "
            "The growth/catastrophe/rescue rate anchors (Walker 1988) are "
            "resolved and stored on ResolvedMicrotubules for that work."
        )

    def act(self, timestep: int) -> None:  # pragma: no cover - unreachable
        raise NotImplementedError
