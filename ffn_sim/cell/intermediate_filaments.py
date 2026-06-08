"""Intermediate-filament (keratin / vimentin) perinuclear cage — explicit chains.

One-line purpose
----------------
ADDITIVE, DEFAULT-OFF module supplying the cell's third cytoskeletal network:
an EXPLICIT keratin (epithelial MCF7) / vimentin intermediate-filament (IF)
perinuclear cage — flexible (low persistence length) bead chains with backbone
bonds + inter-filament crosslink bonds that carry a NON-cortical, extensible
load path around the nucleus. Every IF segment is an explicit ``if_bead``
particle joined by explicit ``if_backbone`` harmonic bonds, with explicit
``if_crosslink`` bonds bridging neighbouring filaments — NO lumped proxy force,
NO mesh-as-physics.

Why a separate network (biology)
--------------------------------
IFs are mechanically distinct from F-actin and microtubules: they are very
FLEXIBLE (persistence length ``Lp`` ~ 0.3–1 µm, ~3 orders below F-actin's
~17 µm) yet EXTENSIBLE and strain-STIFFENING — a single IF extends 2–3.5× its
rest length before rupture (Kreplak 2005; Block 2018), unlike actin which is
near-inextensible and breaks at a few-% strain. In an epithelial cell the IF
network forms a cage around the nucleus and couples to desmosomes / the LINC
complex, bearing large deformations and protecting the nucleus. Per the H.7
direction review (Guo 2013 [IF-10]): VIFs contribute LITTLE to cortical
stiffness but are critical for INTRACELLULAR mechanics — hence this is a
perinuclear/cytoplasm compartment with its OWN bond prefix (``if_``) so it is
denylisted from the cortical-tension (γ) estimator.

Explicit mechanism (what particles / bonds)
-------------------------------------------
* ``if_bead`` — IF segment beads. Each filament is a chain of ``beads_per_fil``
  beads at rest spacing ``l_seg``.
* ``if_backbone`` — permanent harmonic bond between consecutive beads on a
  filament. Small-strain LINEAR spring (this version). Stiffness DERIVED from
  the IF axial Young's modulus and the filament cross-section,
  ``k_bb = E_if · A_if / l_seg`` (grid-invariant rod-segment axial spring).
* ``if_crosslink`` — permanent harmonic cross-bridge between beads on DIFFERENT
  filaments within ``crosslink_reach`` (plectin / filaggrin-type IF
  cross-bridges). Soft relative to the backbone (cross-bridge compliance), and
  born at its as-built cross-filament spacing (per-r0 bin, ``if_crosslink_b{i}``)
  so it is FORCE-FREE at the resting geometry — a true separation-resisting
  tether, NOT a contractile r0=0 spring (see Sanity Gate item 4).

There is NO per-step custom force and NO per-batch updater in this version:
all IF bonds are PERMANENT HOOMD bonds (``md.bond.Harmonic``), evaluated by the
builtin bonded-force kernel on whatever device the integrator runs on. This is
the cheapest correct mechanistic representation; the dynamic (plectin
turnover / strain-stiffening) extensions are documented TODO/PI items below.

NONLINEAR strain-stiffening — NOT faked here (PI item)
------------------------------------------------------
The faithful IF constitutive law is strongly NONLINEAR: a soft small-strain
modulus (few-MPa equivalent) followed by pronounced strain-STIFFENING and very
large extensibility (2–3.5×) before rupture (Kreplak 2005; Block 2018;
Lichtenstern 2012). HOOMD's builtin ``md.bond.Harmonic`` is LINEAR, so this
module ships ONLY the small-strain linear harmonic backbone (with a cited
modulus) and FLAGS the nonlinear strain-stiffening + finite-extensibility +
rupture as a documented TODO (``md.bond.Table`` tabulated potential or a
custom FENE-plus-stiffening bond). We do NOT approximate the nonlinearity with
a tuned harmonic — that would be a magic-number fit. See ``PI_DECISIONS``.

Sanity Gate
-----------
1. **Dimensional analysis**
   - ``A_if = π · (d_if/2)²`` [m²]; ``k_bb = E_if · A_if / l_seg``
     [Pa·m²/m = N/m]. ✓  ``U = ½ k_bb (r − l_seg)²`` [J]. ✓
   - Crosslink stiffness ``k_xl = ratio_xl · k_bb`` [N/m] (dimensionless
     ratio). ✓
   - Persistence length ``Lp = κ / kT`` [J·m / J = m]; bending stiffness
     ``κ = Lp · kT`` [J·m] (retained as provenance; this LINEAR version uses
     no explicit angle term — flexible-chain limit, see boundary cases). ✓
2. **Boundary cases**
   - ``enabled=False`` / missing config → ``ResolvedIntermediateFilaments
     (enabled=False, …zeros…)``; builder early-returns the input snapshot
     unchanged (zero particles, zero bonds). OFF-IDENTITY.
   - Chain of length 1 bead → 0 backbone bonds (no consecutive pair). Builder
     handles ``beads_per_fil = 1`` by adding the bead with no backbone bond.
   - ``E_if`` / ``d_if`` / ``l_seg`` / ``n_filaments`` / ``beads_per_fil`` ≤ 0
     → ``resolve_intermediate_filaments`` ValueError.
   - No crosslink acceptor within ``crosslink_reach`` → that bead gets no
     crosslink bond (the cage is still backbone-connected). Builder handles.
3. **Conservation / no-net-force**
   - All IF bonds are PAIRWISE harmonic (``md.bond.Harmonic``): equal-and-
     opposite forces → the IF network injects ZERO net momentum (internal
     element). No external anchor, no lab-frame pin. Particle count after the
     builder is exactly ``N_old + n_filaments · beads_per_fil``; bond count is
     ``N_old_bonds + n_backbone + n_crosslink``.
   - IF forces act ONLY through ``if_`` bonds among ``if_bead`` particles; no
     other particle type is touched (no custom force, no mask needed).
4. **Sign / sense + force-free construction**
   - A STRETCHED backbone bond (``r > l_seg``) produces a RESTORING force that
     pulls the two beads back TOGETHER (inward along the bond) — standard
     harmonic ``F = −k_bb (r − l_seg) r̂``. Sign test asserts a stretched IF
     pair is pulled toward each other (the explicit two-bead force, computed
     analytically and via a tiny HOOMD bonded eval).
   - CROSSLINK FORCE-FREE: each ``if_crosslink`` bridges beads on DIFFERENT
     filaments at a FINITE as-built separation. A degenerate ``r0 = 0`` on such
     a bond is a CONTRACTILE spring (``F = +k_xl r`` for all ``r > 0``, energy
     minimum at coincidence) — it would inject ~1e6 kT of spurious pre-tension
     and self-contract the cage at t=0, violating the physiological-baseline
     rule. Instead each crosslink is binned (``if_crosslink_b{i}``, N_XL_BINS
     bins over ``(0, crosslink_reach]``) to a rest length ≈ its born separation,
     so the cage construction strain energy is ~O(kT)·n_xl (thermal), NOT
     ~1e6 kT — a true separation-resisting tether. Sign/energy test asserts
     mean ``|r − r0| ≈ 0`` at build (force-free) and a stretched crosslink
     restores inward.
5. **CFL note (stiff)**
   - The stiffest IF spring is the backbone ``k_bb``. Overdamped relax time
     ``τ = γ_if / k_bb`` must satisfy ``dt ≤ cfl_safety_factor · τ`` (same D3
     BAOAB convention as ``cortex/erm.py`` / ``cell/nucleus.py``). The attach
     helper raises if violated (using the stiffest spring, conservative). At
     the cited small-strain modulus IF is SOFT (few MPa over a ~10 nm
     filament → k_bb ~ 1e-3 N/m, ~100× softer than the 0.1 N/m ERM pin), so
     the IF CFL is not the binding step constraint — but the gate is still
     enforced.
6. **Shared bond force (no second Harmonic)**
   - HOOMD 7.0.1 ``md.bond.Harmonic`` demands params for EVERY system bond
     type. Two coexisting Harmonics each crash on the other's types. So for the
     integrated cell the IF bond params are REGISTERED onto the cell's single
     shared ``md.bond.Harmonic`` via ``register_if_bond_params`` (mirroring
     ``stress_fibers`` / ``linc``), NOT a standalone force.
     ``attach_if_bonds_to_simulation`` builds a standalone force and is
     RESTRICTED to IF-only systems — it refuses to attach if the state already
     carries any non-``if_`` bond type.

Compartment Performance Contract
--------------------------------
* Particle types added: ``if_bead`` (1 new type).
* Particle count: ``n_filaments · beads_per_fil``. At the mesoscopic n_fil=1000
  cell, a perinuclear IF cage is a SUBSET of the network (the cage, not all
  ~38k native filaments): default ``n_filaments=200``, ``beads_per_fil=10`` →
  ~2000 ``if_bead``. At native ~38000-filament scale the same relative cage
  fraction is ~7600 filaments → set by the host; the builder is O(N) in the
  beads it is asked to add. The 200×10 default is a host-overridable cage size,
  not a hard count.
* Bond/angle count: backbone ``n_filaments · (beads_per_fil − 1)`` (≈ 1800 at
  default), crosslink ``≤ n_filaments · beads_per_fil`` (cKDTree, one nearest
  cross-filament acceptor per bead), split across ``N_XL_BINS`` per-r0 bin
  subtypes (``if_crosslink_b{i}``) so each is force-free at its born separation.
  NO angles in this linear version (flexible chain limit; bending/strain-
  stiffening is the TODO angle/table term).
* Per-step force: NO custom per-step force — all IF load is via PERMANENT
  ``md.bond.Harmonic`` bonds (builtin bonded kernel).
* Per-batch updater: NO (permanent bonds; no dynamic turnover this version).
* Uses ``cpu_local_snapshot``: NO (snapshot-extension builder uses the global
  ``get_snapshot``/passed Snapshot; no per-step CPU sync).
* Uses cKDTree / broad-phase: YES — once, at BUILD time only, to find
  cross-filament crosslink acceptors within ``crosslink_reach`` (not in the hot
  loop).
* Hot-path priority: P2 (per task). IF bonds add a fixed, modest count of
  builtin harmonic bonds — they ride the existing bonded kernel; no new
  per-step Python.
* GPU path now: builtin (``md.bond.Harmonic`` runs natively on the integrator
  device — CPU or CUDA — no port needed). Build-time cKDTree is CPU/one-shot.
* Native ForceCompute candidate: N/A for the linear backbone (already a builtin
  native bonded force). The TODO nonlinear bond, if not expressible as
  ``md.bond.Table``, would be a native ForceCompute candidate.
* Bottleneck risk: LOW. Fixed bond count, builtin kernel, no per-step Python,
  no per-step snapshot sync. Build-time cKDTree is one-shot O(N log N).

References
----------
- Mücke, N. et al. (2004) "Assessing the flexibility of intermediate filaments
  by atomic force microscopy." J. Mol. Biol. 335(5):1241–1250 — vimentin
  persistence length Lp ≈ 0.4–1 µm (AFM contour analysis).
- Lichtenstern, T. et al. (2012) "Complex formation and kinetics of
  filament assembly exhibited by the simple epithelial keratins K8 and K18."
  J. Struct. Biol. 177(1):54–62 — keratin K8/K18 assembly / flexibility
  (epithelial keratin Lp lower than vimentin, ~0.3–0.5 µm range).
- Kreplak, L. et al. (2005) "Exploring the mechanical behavior of single
  intermediate filaments." J. Mol. Biol. 354(3):569–577 — single IF extends
  ~2–3.5× rest length before rupture; strain-stiffening (the nonlinearity this
  module FLAGS, does not fake).
- Block, J. et al. (2018) "Viscoelastic properties of vimentin originate from
  nonequilibrium conformational changes." Sci. Adv. 4(6):eaat1161,
  DOI 10.1126/sciadv.aat1161 — single-vimentin tensile memory / nonequilibrium
  α-helix unfolding (the strain-stiffening nonlinearity this module FLAGS).
- Guo, M. et al. (2013) "The role of vimentin intermediate filaments in
  cortical and cytoplasmic mechanics." Biophys. J. 105(7):1562–1568 — VIFs
  contribute little to CORTICAL stiffness, dominate INTRACELLULAR mechanics
  (motivates the perinuclear/cytoplasm placement + ``if_`` γ-denylist).
- ``ffn_sim/cell/nucleus.py`` — additive default-off compartment template
  (resolve + dataclass + Sanity Gate + CFL attach gate).
- ``ffn_sim/cell/lamellipodium.py`` — snapshot-extension (particle+bond append)
  builder pattern.
- HOOMD 7 ``md.bond.Harmonic`` / ``md.bond.Table`` API.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.cortex.crosslinkers import (
    xlink_attach_bin_names,
    xlink_attach_bin_rest_lengths,
)


# ---------------------------------------------------------------------------
# Module-level contract constants
# ---------------------------------------------------------------------------
# Every dynamic / load-bearing bond this module creates carries the non-cortical
# IF load path; the cortical-tension (γ) estimator denylists this prefix so IF
# bonds never contaminate the cortical-tension budget.
GAMMA_DENYLIST_PREFIX: str = "if_"

# Open design decisions surfaced to PI (honesty over completeness).
PI_DECISIONS: list[str] = [
    "NONLINEAR strain-stiffening + finite extensibility (2-3.5x, Kreplak 2005 / "
    "Block 2018): this version ships ONLY the small-strain LINEAR harmonic "
    "backbone (cited modulus). The faithful nonlinear strain-stiffening law "
    "needs an md.bond.Table tabulated potential or a custom FENE+stiffening "
    "bond and is NOT faked with a tuned harmonic. TODO/PI: build the tabulated "
    "IF constitutive curve from Kreplak 2005 / Block 2018 force-extension data; "
    "until then build_intermediate_filament_bond_potential(nonlinear=True) "
    "raises NotImplementedError.",
    "E_if (IF axial small-strain Young's modulus): the literature spread is "
    "wide (single-filament few-MPa to network-dependent). Default 6.0 MPa is "
    "the canonical single-IF small-strain order (Kreplak 2005 / Guo 2013 "
    "regime); range-checked [1e6, 1e7] Pa. Tighten per cell type (keratin "
    "vs vimentin) with PI sign-off.",
    "d_if (filament diameter, sets cross-section A_if): 10 nm full assembled "
    "IF diameter (canonical TEM, Mucke 2004 / Herrmann reviews). For the x40 "
    "mesoscopic bundle scale a host may rescale; left at the physical 10 nm "
    "single-IF value (the modulus is per single-IF cross-section). PI to "
    "ratify any mesoscale bundle re-derivation.",
    "ratio_xl (IF crosslink stiffness / backbone): plectin/filaggrin "
    "cross-bridge compliance has no single clean number at this scale; default "
    "0.1 (cross-bridges softer than the backbone, same convention as the "
    "actin xlink << backbone separation). Flagged; PI to anchor.",
]


# ---------------------------------------------------------------------------
# Magic-Number Block (literature anchors / range-check bands)
# ---------------------------------------------------------------------------
# IF axial small-strain Young's modulus band (single-IF, Kreplak 2005 / Guo
# 2013 regime; few MPa). Outside this band the resolver raises.
_E_IF_BAND_PA: tuple[float, float] = (1.0e6, 1.0e7)
# IF persistence length band [m]: keratin ~0.3-0.5 µm (Lichtenstern 2012),
# vimentin ~0.4-1 µm (Mucke 2004). Used to range-check the provided Lp.
_LP_IF_BAND_M: tuple[float, float] = (0.2e-6, 1.2e-6)
# Full assembled IF diameter [m]: canonical ~10 nm (TEM; Mucke 2004 /
# Herrmann reviews). Range-checked 8-12 nm.
_D_IF_BAND_M: tuple[float, float] = (8.0e-9, 12.0e-9)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class ResolvedIntermediateFilaments:
    """Intermediate-filament cage parameters (all SI).

    The backbone stiffness is DERIVED from the IF axial Young's modulus and
    the filament cross-section (``k_bb = E_if · A_if / l_seg``); the crosslink
    stiffness from the dimensionless ``ratio_xl``. Raw anchors are retained as
    provenance.

    Attributes:
        enabled: Master switch. When False the builder is a no-op.
        n_filaments: Number of IF chains in the perinuclear cage (count).
        beads_per_fil: Beads per IF chain (count; ≥ 1).
        l_seg: Backbone bond rest length / bead spacing [m].
        d_if: Assembled filament diameter [m] (sets cross-section).
        A_if: Filament cross-sectional area ``π (d_if/2)²`` [m²] (derived).
        E_if: IF axial small-strain Young's modulus [Pa] (provenance).
        Lp: IF persistence length [m] (provenance; flexible-chain anchor).
        kappa_bend: Bending stiffness ``Lp · kT`` [J·m] (derived PROVENANCE
            ONLY — consumed by NO force in this linear version, which carries no
            angle term; it is the input the documented WLC bending-angle TODO
            would wire as ``k_angle = kappa_bend / l_seg`` into an
            ``md.angle.Harmonic`` (θ0 = π)).
        k_bb: Backbone harmonic stiffness ``E_if A_if / l_seg`` [N/m] (derived).
        ratio_xl: Crosslink stiffness as a fraction of the backbone [—].
        k_xl: Crosslink harmonic stiffness ``ratio_xl · k_bb`` [N/m] (derived).
        crosslink_reach: Max cross-filament distance to place a crosslink [m].
        R_cage_inner: Inner radius of the perinuclear cage shell [m] (seed geo).
        R_cage_outer: Outer radius of the perinuclear cage shell [m] (seed geo).
        max_stretch_ratio: Documented finite-extensibility limit (2-3.5×;
            Kreplak 2005). NOT enforced by the linear bond — provenance for the
            nonlinear TODO.
    """

    enabled: bool
    n_filaments: int
    beads_per_fil: int
    l_seg: float            # m
    d_if: float             # m
    A_if: float             # m²  derived = π (d_if/2)²
    E_if: float             # Pa  provenance
    Lp: float               # m   provenance
    kappa_bend: float       # J·m derived = Lp · kT
    k_bb: float             # N/m derived = E_if A_if / l_seg
    ratio_xl: float         # —
    k_xl: float             # N/m derived = ratio_xl · k_bb
    crosslink_reach: float  # m
    R_cage_inner: float     # m
    R_cage_outer: float     # m
    max_stretch_ratio: float = 3.0   # Kreplak 2005 (2-3.5×); provenance only

    @property
    def n_beads_total(self) -> int:
        """Total ``if_bead`` particles the cage adds."""
        return int(self.n_filaments) * int(self.beads_per_fil)

    @property
    def n_backbone_bonds(self) -> int:
        """Total ``if_backbone`` bonds (chain has beads_per_fil − 1 each)."""
        return int(self.n_filaments) * max(0, int(self.beads_per_fil) - 1)


def _disabled_resolved() -> ResolvedIntermediateFilaments:
    """OFF-identity resolved object: enabled=False with zeros."""
    return ResolvedIntermediateFilaments(
        enabled=False,
        n_filaments=0,
        beads_per_fil=0,
        l_seg=0.0,
        d_if=0.0,
        A_if=0.0,
        E_if=0.0,
        Lp=0.0,
        kappa_bend=0.0,
        k_bb=0.0,
        ratio_xl=0.0,
        k_xl=0.0,
        crosslink_reach=0.0,
        R_cage_inner=0.0,
        R_cage_outer=0.0,
        max_stretch_ratio=0.0,
    )


def resolve_intermediate_filaments(
    cfg: dict,
    *,
    kT: float,
    R_cell: float,
    R_nuc: float | None = None,
) -> ResolvedIntermediateFilaments:
    """Resolve the ``intermediate_filaments`` config block into SI params.

    DEFAULT-OFF: when ``cfg`` is missing the block, or the block has
    ``enabled: false`` (or ``enabled`` absent / falsey), this returns
    ``ResolvedIntermediateFilaments(enabled=False, …zeros…)`` and NOTHING is
    built downstream.

    The backbone stiffness is DERIVED (not invented):
    ``k_bb = E_if · A_if / l_seg`` with ``A_if = π (d_if/2)²`` — the axial
    spring constant of a linear-elastic rod segment of modulus ``E_if``,
    cross-section ``A_if``, length ``l_seg``. Grid-invariant in the sense that
    it is set by the material modulus + geometry, not chosen to pass a gate
    (the sign / dimensional tests do not use ``E_if`` at all).

    Args:
        cfg: Config mapping. May be the YAML root, a ``cell`` sub-dict, or the
            ``intermediate_filaments`` sub-dict directly.
        kT: Thermal energy [J] (sets ``kappa_bend = Lp · kT`` provenance).
        R_cell: Cell radius [m] (REQUIRED; range-checks the cage geometry).
        R_nuc: Nuclear radius [m] (optional). If given, the default cage shell
            is seeded just outside the nucleus ``[R_nuc, R_nuc + thickness]``;
            else a default perinuclear band inside ``R_cell`` is used.

    Returns:
        Resolved IF parameters (SI). ``enabled=False`` when disabled.

    Raises:
        ValueError: on a non-positive / out-of-band parameter when ENABLED.
    """
    # Drill into the relevant sub-dict if a parent mapping was passed.
    if isinstance(cfg, dict) and "cell" in cfg and "intermediate_filaments" in cfg.get("cell", {}):
        cfg = cfg["cell"]
    if isinstance(cfg, dict) and "intermediate_filaments" in cfg:
        cfg = cfg["intermediate_filaments"]
    if cfg is None:
        cfg = {}

    # DEFAULT-OFF gate.
    if not bool(cfg.get("enabled", False)):
        return _disabled_resolved()

    _require_finite_positive("kT", kT)
    _require_finite_positive("R_cell", R_cell)

    # --- Geometry / counts ---
    n_filaments = int(cfg.get("n_filaments", 200))
    beads_per_fil = int(cfg.get("beads_per_fil", 10))
    if n_filaments < 1:
        raise ValueError(f"n_filaments must be ≥ 1 when enabled; got {n_filaments!r}")
    if beads_per_fil < 1:
        raise ValueError(f"beads_per_fil must be ≥ 1; got {beads_per_fil!r}")

    # --- Material anchors (Magic-Number Block) ---
    E_if = float(cfg.get("E_if", 6.0e6))      # Pa  Kreplak 2005 / Guo 2013 (few MPa)
    Lp = float(cfg.get("Lp", 0.5e-6))         # m   Mucke 2004 / Lichtenstern 2012
    d_if = float(cfg.get("d_if", 10.0e-9))    # m   canonical ~10 nm assembled IF
    # Backbone rest length / bead spacing. Default = one persistence length so
    # the chain is flexible. NOTE: the Kuhn length is b = 2·Lp, so l_seg = Lp is
    # HALF a Kuhn segment, not one. This LINEAR version carries NO explicit
    # bending angle term (flexible-chain / FJC limit), so the cage is
    # INTENTIONALLY floppier than a stiffness-matched WLC (an angle-free chain
    # with sub-Kuhn bonds under-reports ⟨R²⟩ by up to ~2×). The bending
    # angle / strain-stiffening term is the documented TODO (see the module
    # docstring "NONLINEAR strain-stiffening" + PI_DECISIONS). Set
    # ``l_seg = 2·Lp`` if WLC large-scale stats must be matched by an angle-free
    # FJC. Grid-derivable, not tuned.
    l_seg = float(cfg.get("l_seg", Lp))
    ratio_xl = float(cfg.get("ratio_xl", 0.1))
    max_stretch_ratio = float(cfg.get("max_stretch_ratio", 3.0))

    _require_finite_positive("E_if", E_if)
    _require_finite_positive("Lp", Lp)
    _require_finite_positive("d_if", d_if)
    _require_finite_positive("l_seg", l_seg)
    _require_finite_positive("ratio_xl", ratio_xl)
    _require_finite_positive("max_stretch_ratio", max_stretch_ratio)

    # Band guards (no silent tuned override).
    elo, ehi = _E_IF_BAND_PA
    if not (elo <= E_if <= ehi):
        raise ValueError(
            f"E_if = {E_if:.3e} Pa outside the IF small-strain modulus band "
            f"[{elo:.1e}, {ehi:.1e}] Pa (few-MPa single-IF; Kreplak 2005 / "
            "Guo 2013). Keep within band or surface to PI."
        )
    llo, lhi = _LP_IF_BAND_M
    if not (llo <= Lp <= lhi):
        raise ValueError(
            f"Lp = {Lp:.3e} m outside the IF persistence-length band "
            f"[{llo:.1e}, {lhi:.1e}] m (keratin ~0.3-0.5 µm Lichtenstern 2012; "
            "vimentin ~0.4-1 µm Mucke 2004). Keep within band or surface to PI."
        )
    dlo, dhi = _D_IF_BAND_M
    if not (dlo <= d_if <= dhi):
        raise ValueError(
            f"d_if = {d_if:.3e} m outside the assembled-IF diameter band "
            f"[{dlo:.1e}, {dhi:.1e}] m (canonical ~10 nm; Mucke 2004). Keep "
            "within band or surface to PI."
        )

    # --- Derived SI quantities ---
    A_if = math.pi * (0.5 * d_if) ** 2               # m²  cross-section
    k_bb = E_if * A_if / l_seg                        # N/m rod-segment axial spring
    k_xl = ratio_xl * k_bb                            # N/m crosslink (softer)
    kappa_bend = Lp * kT                              # J·m bending stiffness (PROVENANCE ONLY; no force consumes it — future angle-term input)

    # --- Cage shell geometry ---
    default_thickness = float(cfg.get("cage_thickness", 0.3 * R_cell))
    if R_nuc is not None and R_nuc > 0.0:
        R_cage_inner = float(cfg.get("R_cage_inner", R_nuc))
        R_cage_outer = float(cfg.get("R_cage_outer", R_nuc + default_thickness))
    else:
        # Perinuclear band sitting in the inner ~[0.3, 0.6] R_cell by default.
        R_cage_inner = float(cfg.get("R_cage_inner", 0.30 * R_cell))
        R_cage_outer = float(cfg.get("R_cage_outer", 0.60 * R_cell))
    if not (0.0 < R_cage_inner < R_cage_outer <= R_cell):
        raise ValueError(
            f"cage shell must satisfy 0 < R_cage_inner ({R_cage_inner:.3e}) "
            f"< R_cage_outer ({R_cage_outer:.3e}) ≤ R_cell ({R_cell:.3e}) [m]."
        )

    # Crosslink reach defaults to one segment (nearest-neighbour filaments).
    crosslink_reach = float(cfg.get("crosslink_reach", l_seg))
    _require_finite_positive("crosslink_reach", crosslink_reach)

    return ResolvedIntermediateFilaments(
        enabled=True,
        n_filaments=n_filaments,
        beads_per_fil=beads_per_fil,
        l_seg=l_seg,
        d_if=d_if,
        A_if=A_if,
        E_if=E_if,
        Lp=Lp,
        kappa_bend=kappa_bend,
        k_bb=k_bb,
        ratio_xl=ratio_xl,
        k_xl=k_xl,
        crosslink_reach=crosslink_reach,
        R_cage_inner=R_cage_inner,
        R_cage_outer=R_cage_outer,
        max_stretch_ratio=max_stretch_ratio,
    )


# ---------------------------------------------------------------------------
# Bond-type names (all carry the GAMMA_DENYLIST_PREFIX)
# ---------------------------------------------------------------------------
BACKBONE_BOND_TYPE: str = "if_backbone"
CROSSLINK_BOND_TYPE: str = "if_crosslink"
IF_BEAD_TYPE: str = "if_bead"

# Crosslink rest-length binning. Each cross-bridge is force-free at its as-built
# cross-filament separation: instead of a single degenerate r0=0 (which would be
# a CONTRACTILE spring across already-separated beads, injecting spurious
# construction pre-stress — see Sanity Gate item 4 below), the crosslink bond is
# split into N_XL_BINS per-r0 subtypes over (0, crosslink_reach], and each
# crosslink is assigned the bin whose centre ≈ its born separation. This mirrors
# the platform convention in cortex/crosslinkers.py:xlink_attach_bin_rest_lengths
# (and cortex/connected_mesh.py / cortex/cortex.py), so the IF cage is force-free
# (energy ~ O(kT)·n_xl, thermal) at t=0 rather than ~1e6 kT pre-tensioned.
N_XL_BINS: int = 16


def if_crosslink_bin_names(n_bins: int = N_XL_BINS) -> list[str]:
    """HOOMD bond-type names for the per-r0 IF crosslink bins (γ-denylisted)."""
    return [f"if_crosslink_b{i}" for i in range(int(n_bins))]


assert BACKBONE_BOND_TYPE.startswith(GAMMA_DENYLIST_PREFIX)
assert CROSSLINK_BOND_TYPE.startswith(GAMMA_DENYLIST_PREFIX)
assert all(n.startswith(GAMMA_DENYLIST_PREFIX) for n in if_crosslink_bin_names())


# ---------------------------------------------------------------------------
# PURE layout helper (computes & returns; NEVER mutates a sim/file/snapshot)
# ---------------------------------------------------------------------------
def build_if_cage_layout(
    p: ResolvedIntermediateFilaments,
    centroid: tuple[float, float, float] | np.ndarray,
    *,
    gamma_if: float,
    seed: int = 0,
) -> dict[str, object]:
    """Compute the perinuclear IF cage layout (positions + bond pairs).

    PURE helper: it computes ``if_bead`` positions, the per-filament backbone
    bond pairs, and the cross-filament crosslink bond pairs, and RETURNS them
    as local-index arrays. It MUST NOT mutate any simulation / snapshot / file.
    The builder / Lead offsets these local indices by the existing particle
    count when appending.

    Each filament is a short chain of ``beads_per_fil`` beads laid along a
    random great-circle-tangent direction, seeded in the perinuclear shell
    ``[R_cage_inner, R_cage_outer]`` about ``centroid`` (deterministic RNG).
    Crosslinks are placed by a one-shot cKDTree query: for each bead, the
    nearest bead on a DIFFERENT filament within ``crosslink_reach`` (one
    crosslink per bead, de-duplicated).

    Args:
        p: Resolved IF parameters (must be ``enabled``).
        centroid: Cell / nucleus centroid ``(x, y, z)`` [m].
        gamma_if: Per-bead Stokes drag [N·s/m] (= 6π η R_bead; > 0).
        seed: Deterministic RNG seed offset.

    Returns:
        Dict with:
          * ``"positions"``: ndarray ``(n_beads_total, 3)`` float64 [m].
          * ``"type_name"``: str (``"if_bead"``).
          * ``"gamma"``: ndarray ``(n_beads_total,)`` float64 [N·s/m].
          * ``"backbone_bonds"``: ndarray ``(n_bb, 2)`` int64 LOCAL indices.
          * ``"crosslink_bonds"``: ndarray ``(n_xl, 2)`` int64 LOCAL indices.
          * ``"crosslink_dists"``: ndarray ``(n_xl,)`` float64 — each crosslink's
            as-built cross-filament separation [m] (its force-free rest length).
          * ``"crosslink_bins"``: ndarray ``(n_xl,)`` int64 — per-crosslink
            per-r0 bin index (assigned so the bond is force-free at its born
            separation; see ``N_XL_BINS``).
          * ``"crosslink_bin_rest_lengths"``: ndarray ``(N_XL_BINS,)`` float64 —
            the bin-centre rest lengths [m].
          * ``"crosslink_bin_types"``: list[str] — the per-bin bond-type names.
          * ``"backbone_bond_type"`` / ``"crosslink_bond_type"``: str.

    Raises:
        RuntimeError: if ``p`` is disabled (caller must gate on ``p.enabled``).
        ValueError: on non-positive ``gamma_if``.
    """
    if not p.enabled:
        raise RuntimeError(
            "build_if_cage_layout called on a disabled ResolvedIntermediate"
            "Filaments; gate on p.enabled before calling."
        )
    _require_finite_positive("gamma_if", gamma_if)

    c = np.asarray(centroid, dtype=np.float64).reshape(3)
    if not np.isfinite(c).all():
        raise ValueError(f"centroid must be finite; got {centroid!r}")

    rng = np.random.default_rng(int(seed) + 9173)
    nf = int(p.n_filaments)
    bpf = int(p.beads_per_fil)

    positions = np.empty((nf * bpf, 3), dtype=np.float64)
    fil_id = np.empty(nf * bpf, dtype=np.int64)
    backbone_bonds: list[tuple[int, int]] = []

    # Even angular seeding of filament centres on the perinuclear shell
    # (Fibonacci directions), each chain laid along a tangent direction.
    idx = np.arange(nf, dtype=np.float64) + 0.5
    golden = math.pi * (3.0 - math.sqrt(5.0))
    cos_t = np.clip(1.0 - 2.0 * idx / float(nf), -1.0, 1.0)
    sin_t = np.sqrt(np.maximum(0.0, 1.0 - cos_t * cos_t))
    az = golden * idx
    seed_dirs = np.stack(
        [sin_t * np.cos(az), sin_t * np.sin(az), cos_t], axis=1
    )
    # Radii spread across the shell (each filament centre at a random shell r).
    shell_r = rng.uniform(p.R_cage_inner, p.R_cage_outer, size=nf)

    bead = 0
    for f in range(nf):
        n_hat = seed_dirs[f]
        centre = c + shell_r[f] * n_hat
        # A tangent direction (perpendicular to the radial normal): build an
        # arbitrary perpendicular and rotate randomly about n_hat.
        ref = np.array([0.0, 0.0, 1.0]) if abs(n_hat[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        t1 = np.cross(n_hat, ref)
        t1 /= np.linalg.norm(t1)
        t2 = np.cross(n_hat, t1)
        ang = rng.uniform(0.0, 2.0 * math.pi)
        tangent = math.cos(ang) * t1 + math.sin(ang) * t2
        # Lay the chain centred on `centre` along `tangent` at spacing l_seg.
        offsets = (np.arange(bpf, dtype=np.float64) - 0.5 * (bpf - 1)) * p.l_seg
        first = bead
        for b in range(bpf):
            positions[bead] = centre + offsets[b] * tangent
            fil_id[bead] = f
            if b > 0:
                backbone_bonds.append((bead - 1, bead))
            bead += 1
        _ = first  # (chain start index; kept for clarity)

    # --- Crosslinks: one-shot cKDTree, nearest cross-filament bead in reach ---
    # Each crosslink is recorded with its as-built cross-filament separation
    # ``d_in`` so the builder can place it FORCE-FREE at construction (per-r0
    # bin), not as a contractile r0=0 spring (Sanity Gate item 4).
    crosslink_bonds: list[tuple[int, int]] = []
    crosslink_dists: list[float] = []
    if p.crosslink_reach > 0.0 and positions.shape[0] > 1:
        try:
            from scipy.spatial import cKDTree  # build-time only, not hot loop

            tree = cKDTree(positions)
            # Query a few neighbours so we can skip same-filament beads.
            kq = min(positions.shape[0], 6)
            dists, nbrs = tree.query(positions, k=kq)
            seen: set[tuple[int, int]] = set()
            for i in range(positions.shape[0]):
                for d_in, j in zip(np.atleast_1d(dists[i]), np.atleast_1d(nbrs[i])):
                    j = int(j)
                    if j == i:
                        continue
                    if fil_id[j] == fil_id[i]:
                        continue  # same filament → backbone, not a crosslink
                    if d_in > p.crosslink_reach:
                        continue
                    key = (i, j) if i < j else (j, i)
                    if key in seen:
                        continue
                    seen.add(key)
                    crosslink_bonds.append(key)
                    crosslink_dists.append(float(d_in))
                    break  # one crosslink per bead (nearest cross-filament)
        except Exception:  # scipy unavailable / degenerate → cage is still
            crosslink_bonds = []  # backbone-connected; crosslinks optional.
            crosslink_dists = []

    xl_arr = (
        np.array(crosslink_bonds, dtype=np.int64).reshape(-1, 2)
        if crosslink_bonds else np.empty((0, 2), dtype=np.int64)
    )
    xl_dist = (
        np.array(crosslink_dists, dtype=np.float64)
        if crosslink_dists else np.empty((0,), dtype=np.float64)
    )
    # Assign each crosslink the per-r0 bin whose centre is nearest its born
    # separation → the bond is (near-)force-free at construction. Bin edges span
    # (0, crosslink_reach]; clamp any tiny float overshoot into the last bin.
    bin_centres = xlink_attach_bin_rest_lengths(N_XL_BINS, float(p.crosslink_reach))
    if xl_dist.shape[0] > 0:
        bin_width = float(p.crosslink_reach) / float(N_XL_BINS)
        xl_bin = np.clip(
            (xl_dist / bin_width).astype(np.int64), 0, N_XL_BINS - 1
        )
    else:
        xl_bin = np.empty((0,), dtype=np.int64)

    return {
        "positions": positions,
        "type_name": IF_BEAD_TYPE,
        "gamma": np.full(positions.shape[0], float(gamma_if), dtype=np.float64),
        "backbone_bonds": (
            np.array(backbone_bonds, dtype=np.int64).reshape(-1, 2)
            if backbone_bonds else np.empty((0, 2), dtype=np.int64)
        ),
        "crosslink_bonds": xl_arr,
        "crosslink_dists": xl_dist,         # born cross-filament separation [m]
        "crosslink_bins": xl_bin,           # per-r0 bin index per crosslink
        "crosslink_bin_rest_lengths": bin_centres,  # bin-centre r0 [m]
        "backbone_bond_type": BACKBONE_BOND_TYPE,
        "crosslink_bond_type": CROSSLINK_BOND_TYPE,
        "crosslink_bin_types": if_crosslink_bin_names(N_XL_BINS),
    }


# ---------------------------------------------------------------------------
# Snapshot-extension builder (callable in isolation; NO-OP when disabled)
# ---------------------------------------------------------------------------
def extend_snapshot_with_if_cage(
    read_snap,
    p: ResolvedIntermediateFilaments,
    centroid: tuple[float, float, float] | np.ndarray,
    *,
    gamma_if: float,
    seed: int = 0,
):
    """Return a snapshot extended with the explicit IF cage (or the input).

    DEFAULT-OFF / OFF-IDENTITY: when ``p.enabled`` is False this returns
    ``read_snap`` UNCHANGED (same object, zero particles, zero bonds added).

    When enabled it builds a FRESH ``hoomd.Snapshot`` that appends:
      * ``n_beads_total`` ``if_bead`` particles (perinuclear cage),
      * ``if_backbone`` harmonic bonds (chain backbone),
      * per-r0-bin ``if_crosslink_b{i}`` harmonic bonds (cross-filament
        bridges, each force-free at its as-built separation),
    leaving every pre-existing particle / bond / angle untouched (mirrors the
    ``cell/lamellipodium.py`` ``_extend_snapshot_with_new_actins`` pattern).

    Args:
        read_snap: A ``hoomd.Snapshot`` (or ``sim.state.get_snapshot()`` result)
            to extend. NOT mutated; a new snapshot is returned when enabled.
        p: Resolved IF parameters.
        centroid: Cell / nucleus centroid ``(x, y, z)`` [m].
        gamma_if: Per-bead Stokes drag [N·s/m] (> 0). NOTE: the snapshot itself
            does not store γ (the host BAOAB ``gamma_map`` does); ``gamma_if``
            is validated and passed through to ``build_if_cage_layout`` so the
            caller can read it from the returned layout if needed.
        seed: Deterministic RNG seed offset.

    Returns:
        The extended ``hoomd.Snapshot`` (enabled) or ``read_snap`` (disabled).
    """
    if not p.enabled:
        return read_snap  # OFF-IDENTITY: input unchanged, zero particles/bonds.

    layout = build_if_cage_layout(p, centroid, gamma_if=gamma_if, seed=seed)
    new_pos = np.asarray(layout["positions"], dtype=np.float64).reshape(-1, 3)
    bb = np.asarray(layout["backbone_bonds"], dtype=np.int64).reshape(-1, 2)
    xl = np.asarray(layout["crosslink_bonds"], dtype=np.int64).reshape(-1, 2)
    xl_bin = np.asarray(layout["crosslink_bins"], dtype=np.int64).reshape(-1)

    n_old = int(read_snap.particles.N)
    n_new = new_pos.shape[0]
    n_total = n_old + n_new

    write_snap = hoomd.Snapshot()
    write_snap.particles.N = n_total

    # Particle types: ensure if_bead registered.
    types = list(read_snap.particles.types)
    if IF_BEAD_TYPE not in types:
        types.append(IF_BEAD_TYPE)
    if_typeid = types.index(IF_BEAD_TYPE)
    write_snap.particles.types = types

    # A build-time gsd.hoomd.Frame returns None for unset per-particle fields
    # (typeid/velocity/mass/image); HOOMD fills them on load. Mirror those
    # auto-population defaults so the in-build path (snapshot extended BEFORE
    # create_state_from_snapshot) does not crash on np.asarray(None) (0-d). The
    # standalone HOOMD-snapshot path already has these populated. (Mirrors the
    # microtubules / cadherin extenders' None-guards.)
    def _pf(arr, default: np.ndarray) -> np.ndarray:
        return default if arr is None else np.asarray(arr)

    old_typeid = _pf(
        read_snap.particles.typeid, np.zeros(n_old, dtype=np.uint32)
    ).reshape(-1)
    new_typeid = np.full(n_new, if_typeid, dtype=np.uint32)
    write_snap.particles.typeid[:] = np.concatenate([old_typeid, new_typeid])

    old_pos = _pf(
        read_snap.particles.position, np.zeros((n_old, 3), dtype=np.float64)
    ).reshape(-1, 3)
    write_snap.particles.position[:] = np.concatenate([old_pos, new_pos], axis=0)

    old_vel = _pf(
        read_snap.particles.velocity, np.zeros((n_old, 3), dtype=np.float64)
    ).reshape(-1, 3)
    write_snap.particles.velocity[:] = np.concatenate(
        [old_vel, np.zeros((n_new, 3), dtype=np.float64)], axis=0
    )
    old_mass = _pf(
        read_snap.particles.mass, np.ones(n_old, dtype=np.float64)
    ).reshape(-1)
    write_snap.particles.mass[:] = np.concatenate(
        [old_mass, np.ones(n_new, dtype=np.float64)]
    )
    old_image = _pf(
        read_snap.particles.image, np.zeros((n_old, 3), dtype=np.int32)
    ).reshape(-1, 3)
    write_snap.particles.image[:] = np.concatenate(
        [old_image, np.zeros((n_new, 3), dtype=np.int32)], axis=0
    )
    write_snap.configuration.box = list(read_snap.configuration.box)

    # Bonds: existing + if_backbone + per-r0-bin if_crosslink (offset
    # local→global). Crosslinks use one bond subtype per r0 bin so each is
    # FORCE-FREE at its born cross-filament separation (Sanity Gate item 4), not
    # a single contractile r0=0 type.
    bond_types = list(read_snap.bonds.types)
    if BACKBONE_BOND_TYPE not in bond_types:
        bond_types.append(BACKBONE_BOND_TYPE)
    xl_bin_types = if_crosslink_bin_names(N_XL_BINS)
    for name in xl_bin_types:
        if name not in bond_types:
            bond_types.append(name)
    bb_typeid = bond_types.index(BACKBONE_BOND_TYPE)
    xl_bin_typeids = np.array(
        [bond_types.index(name) for name in xl_bin_types], dtype=np.uint32
    )

    old_bg = np.asarray(read_snap.bonds.group, dtype=np.int64).reshape(-1, 2)
    old_bt = np.asarray(read_snap.bonds.typeid, dtype=np.uint32).reshape(-1)
    bb_global = bb + n_old if bb.shape[0] else bb
    xl_global = xl + n_old if xl.shape[0] else xl
    # Map each crosslink's bin index → its per-bin bond typeid.
    xl_typeids = (
        xl_bin_typeids[xl_bin] if xl_global.shape[0] else
        np.empty((0,), dtype=np.uint32)
    )
    bg_pieces = [old_bg, bb_global, xl_global]
    bt_pieces = [
        old_bt,
        np.full(bb_global.shape[0], bb_typeid, dtype=np.uint32),
        xl_typeids,
    ]
    merged_bg = np.concatenate(bg_pieces, axis=0).astype(np.uint32)
    merged_bt = np.concatenate(bt_pieces).astype(np.uint32)
    write_snap.bonds.N = int(merged_bg.shape[0])
    write_snap.bonds.types = bond_types
    if merged_bg.shape[0] > 0:
        write_snap.bonds.group[:] = merged_bg
        write_snap.bonds.typeid[:] = merged_bt

    # Pass-through angles / dihedrals / impropers untouched. A build-time gsd
    # Frame returns None for N / types of an empty group → guard both (the
    # standalone HOOMD-snapshot path has them as 0 / []).
    for grp_name in ("angles", "dihedrals", "impropers"):
        src = getattr(read_snap, grp_name)
        dst = getattr(write_snap, grp_name)
        src_n = int(src.N) if getattr(src, "N", None) else 0
        src_types = list(src.types) if getattr(src, "types", None) else []
        if src_n > 0:
            dst.N = src_n
            dst.types = src_types
            dst.group[:] = np.asarray(src.group)
            dst.typeid[:] = np.asarray(src.typeid)
        elif src_types:
            dst.N = 0
            dst.types = src_types

    return write_snap


# ---------------------------------------------------------------------------
# Force attach: register the IF harmonic bonds on an existing simulation
# ---------------------------------------------------------------------------
def _raise_if_nonlinear(nonlinear: bool) -> None:
    """Guard: the faithful nonlinear IF law is never silently faked."""
    if nonlinear:
        raise NotImplementedError(
            "Nonlinear IF strain-stiffening / finite-extensibility (Kreplak "
            "2005 2-3.5×; Block 2018) is NOT implemented. The faithful law "
            "needs an md.bond.Table tabulated potential or a custom FENE+"
            "stiffening bond built from force-extension data — it is NOT "
            "approximated with a tuned harmonic. See PI_DECISIONS."
        )


def register_if_bond_params(
    bond: "md.bond.Harmonic",
    p: ResolvedIntermediateFilaments,
    *,
    nonlinear: bool = False,
) -> "md.bond.Harmonic":
    """Register the IF bond params ON the cell's SHARED ``md.bond.Harmonic``.

    Platform shared-force pattern (mirrors
    ``stress_fibers.register_stress_fiber_bond_params`` /
    ``linc.configure_linc_bond_potential``): the Lead owns the single
    ``md.bond.Harmonic`` instance for the cell, and each compartment writes its
    own bond-type params onto it — it does NOT construct a SECOND competing
    ``md.bond.Harmonic`` (HOOMD 7.0.1 demands params for EVERY system bond type
    on every Harmonic ForceCompute, so two coexisting Harmonics each crash on
    the other's types — Sanity Gate item 6).

    Writes:
      * ``if_backbone``      → ``dict(k=k_bb, r0=l_seg)`` (force-free at the
        as-built bead spacing).
      * ``if_crosslink_b{i}`` (i in ``range(N_XL_BINS)``) → ``dict(k=k_xl,
        r0=bin_centre_i)`` — a harmonic cross-bridge at its as-built rest length
        (per-r0 bin), force-free at the resting cross-filament spacing (NOT a
        contractile r0=0 spring; Sanity Gate item 4).

    No-op (returns ``bond`` unchanged) when ``p`` is disabled.

    Args:
        bond: The cell's shared ``md.bond.Harmonic`` to write onto.
        p: Resolved IF parameters.
        nonlinear: If True, raises NotImplementedError (see PI_DECISIONS).

    Returns:
        The same ``bond`` object (params populated when enabled).

    Raises:
        NotImplementedError: when ``nonlinear=True`` (documented PI/TODO item).
    """
    _raise_if_nonlinear(nonlinear)
    if not p.enabled:
        return bond
    # if_backbone: rest length = bead spacing l_seg, stiffness k_bb (force-free
    # at the as-built spacing).
    bond.params[BACKBONE_BOND_TYPE] = dict(k=p.k_bb, r0=p.l_seg)
    # if_crosslink_b{i}: per-r0 cross-bridge. Each crosslink rides the bin whose
    # centre ≈ its born cross-filament separation, so the bond is force-free at
    # the resting spacing (a true separation-resisting tether, NOT r0=0).
    bin_r0 = xlink_attach_bin_rest_lengths(N_XL_BINS, float(p.crosslink_reach))
    for i, name in enumerate(if_crosslink_bin_names(N_XL_BINS)):
        bond.params[name] = dict(k=p.k_xl, r0=float(bin_r0[i]))
    return bond


def build_intermediate_filament_bonds(
    p: ResolvedIntermediateFilaments,
    *,
    nonlinear: bool = False,
) -> "md.bond.Harmonic | None":
    """Build a STANDALONE ``md.bond.Harmonic`` for an IF-ONLY system.

    For the integrated cell, prefer :func:`register_if_bond_params`, which writes
    onto the cell's single shared ``md.bond.Harmonic`` (a second standalone
    Harmonic would crash any system that carries other bond types — HOOMD 7.0.1
    demands params for every bond type on every Harmonic; Sanity Gate item 6).
    This standalone constructor is for diagnostic / unit IF-only snapshots whose
    only bond types are ``if_backbone`` + ``if_crosslink_b{i}``.

    Small-strain LINEAR backbone + per-r0-bin crosslink springs (this version).
    Returns a configured ``md.bond.Harmonic`` with params set for ``if_backbone``
    (``k_bb``, rest length ``l_seg``) and each ``if_crosslink_b{i}`` (``k_xl``,
    rest length = the bin-centre cross-filament spacing — force-free at its
    as-built separation, NOT a contractile r0=0 tether).

    Args:
        p: Resolved IF parameters.
        nonlinear: If True, request the faithful strain-stiffening +
            finite-extensibility constitutive law — NOT implemented (would be a
            magic-number fit with a plain harmonic). Raises NotImplementedError
            so the nonlinearity is never silently faked (see ``PI_DECISIONS``).

    Returns:
        A configured ``md.bond.Harmonic`` (enabled) or ``None`` (disabled).

    Raises:
        NotImplementedError: when ``nonlinear=True`` (documented PI/TODO item).
    """
    _raise_if_nonlinear(nonlinear)
    if not p.enabled:
        return None
    harmonic = md.bond.Harmonic()
    return register_if_bond_params(harmonic, p, nonlinear=False)


def attach_if_bonds_to_simulation(
    sim: hoomd.Simulation,
    p: ResolvedIntermediateFilaments,
    *,
    gamma_if: float | None = None,
    cfl_safety_factor: float = 0.1,
    cfl_strict: bool = True,
) -> "md.bond.Harmonic | None":
    """CFL-gate + append a STANDALONE IF harmonic-bond force (no-op off).

    DEFAULT-OFF: returns ``None`` and touches nothing when ``p.enabled`` is
    False.

    ⚠️ For the INTEGRATED cell, register the IF params onto the cell's single
    shared ``md.bond.Harmonic`` via :func:`register_if_bond_params` instead —
    this helper appends a SECOND standalone Harmonic, which is valid ONLY for an
    IF-ONLY system (no other bond types). HOOMD 7.0.1 demands params for every
    system bond type on every Harmonic, so a standalone IF force coexisting with
    e.g. a cortex bond crashes (Sanity Gate item 6). It refuses to attach if the
    sim already carries non-IF bond types.

    CFL gate (mirrors ``cortex/erm.py`` / ``cell/nucleus.py``): the stiffest IF
    spring is the backbone ``k_bb``; ``τ = γ_if / k_bb`` and ``dt`` must satisfy
    ``dt ≤ cfl_safety_factor · τ``. Raises ``RuntimeError`` if violated (when
    ``gamma_if`` is supplied and ``cfl_strict``).

    Args:
        sim: Already-built simulation (must have an Integrator with ``dt``).
        p: Resolved IF parameters.
        gamma_if: Per-bead Stokes drag [N·s/m] for the CFL gate (optional; if
            None the gate is skipped, caller's responsibility).
        cfl_safety_factor: Same convention as the bond/angle/ERM CFL gates.
        cfl_strict: If True, raise on CFL violation.

    Returns:
        The attached ``md.bond.Harmonic`` (enabled) or ``None`` (disabled).

    Raises:
        RuntimeError: if no Integrator is set, or the CFL gate is violated.
    """
    if not p.enabled:
        return None
    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "sim.operations.integrator must be set before attaching IF bonds."
        )
    # Coexistence guard (Sanity Gate item 6): a standalone IF Harmonic is only
    # safe in an IF-ONLY system. If the state already carries non-IF bond types,
    # a second Harmonic would crash (HOOMD demands params for every type on every
    # Harmonic) — route through register_if_bond_params on the shared force.
    try:
        sys_bond_types = list(sim.state.bond_types)
    except Exception:
        sys_bond_types = []
    foreign = [
        t for t in sys_bond_types
        if not str(t).startswith(GAMMA_DENYLIST_PREFIX)
    ]
    if foreign:
        raise RuntimeError(
            "attach_if_bonds_to_simulation builds a STANDALONE md.bond.Harmonic, "
            f"but the system already carries non-IF bond types {foreign!r}. Two "
            "coexisting Harmonics each demand params for the other's bond types "
            "and crash (HOOMD 7.0.1). For the integrated cell, register IF params "
            "onto the cell's shared md.bond.Harmonic via register_if_bond_params "
            "instead of attaching a second force."
        )
    if gamma_if is not None and p.k_bb > 0.0:
        _require_finite_positive("gamma_if", gamma_if)
        dt = float(ig.dt)
        tau_if = gamma_if / p.k_bb
        dt_cfl_if = cfl_safety_factor * tau_if
        if dt > dt_cfl_if and cfl_strict:
            raise RuntimeError(
                f"IF backbone CFL violated: dt = {dt:.3e} s > "
                f"{cfl_safety_factor:.2f} · τ_if = {dt_cfl_if:.3e} s "
                f"(τ_if = γ_if / k_bb = {tau_if:.3e} s, k_bb = {p.k_bb:.3e} "
                "N/m). Reduce dt or soften k_bb (PI sign-off vs the cited "
                "E_if). Pass cfl_strict=False to skip for diagnostic runs."
            )
    harmonic = build_intermediate_filament_bonds(p, nonlinear=False)
    ig.forces.append(harmonic)
    return harmonic
