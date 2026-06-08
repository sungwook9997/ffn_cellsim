"""H.8 plasma-membrane RESERVOIR + bleb-nucleation machinery (KU-3.B1 / KU-3.5).

ADDITIVE, DEFAULT-OFF module. The plasma membrane carries EXCESS area (folds,
microvilli, caveolae) that buffers in-plane tension, and is held against the
actomyosin cortex by a finite **membrane–cortex adhesion** (the MCA / ezrin-
radixin-moesin link). When the cytoplasmic hydrostatic pressure locally exceeds
the adhesion, the membrane DETACHES from the cortex and a **bleb** nucleates
(Charras & Paluch 2008; Tinevez et al. 2009). This module supplies the missing
fine-grained machinery for BOTH facts:

1. **Membrane–cortex tether bonds** (``mem_tether``): explicit breakable harmonic
   bonds, one per anchored membrane bead, joining a membrane-surface bead to its
   nearest cortex bead. Each tether carries a real adhesion energy / rupture
   load. A bond ruptures (Bell-Evans) when the local detachment load exceeds the
   adhesion — that rupture IS the bleb-nucleation event (an emergent, explicit
   topology change, not a lumped switch).
2. **Area reservoir** (a bookkeeping reserve of excess membrane area): the
   construction shell holds a fraction ``f_excess`` of slack area folded away;
   as the membrane is stretched, reservoir area is RELEASED to grow the elastic
   reference area ``A0``, buffering tension before the area-elastic term of
   :mod:`ffn_sim.cell.membrane_surface` engages. This is the explicit area-store
   that lets tension stay near-constant over the buffering regime
   (Raucher & Sheetz 1999; Figard & Sokac 2014).

This file is DISTINCT from :mod:`ffn_sim.cell.membrane_surface` (the H.8 in-plane
γ_membrane + K_A surface module). **It REQUIRES membrane_surface**: the membrane
beads it tethers and the reference area ``A0`` it grows are the membrane_surface
shell + reference area. It is conceptually coupled to :mod:`ffn_sim.cortex.erm`
(the ERM radial tether is the lab-frame analog of the per-bead membrane–cortex
link; the adhesion energy here is the energetic counterpart of that pinning).

Explicit mechanism — what particles / bonds (per CLAUDE.md fine-grained rule)
----------------------------------------------------------------------------
* Particles added: **NONE** at construction in the default build — the tethers
  reuse EXISTING membrane-surface beads and EXISTING cortex beads (the membrane
  rides on the same shell tags as the cortex in the H.8 surface model). The
  bleb-membrane patch that detaches is the same membrane bead set, now unbonded
  from the cortex. (A future variant that gives the membrane its OWN bead layer
  offset radially outward from the cortex is a documented TODO; it would add
  ``n_mem`` ``mem_bead`` particles. NOT done here to stay strictly additive and
  avoid an un-anchored geometry.)
* Bonds added: ``mem_tether`` — one explicit breakable harmonic bond per anchored
  membrane bead → nearest cortex bead. Every bond type starts with the
  :data:`GAMMA_DENYLIST_PREFIX` ``'mem_'`` so the cortical-tension estimator
  (:func:`ffn_sim.cortex.cortical_tension`) DENYLISTS it — the tether is a
  non-cortical, radial detachment load path and must NOT contaminate γ_soft.
* Updater: a per-batch **Bell-Evans rupture** Action (:class:`MembraneTetherUpdater`)
  removes ``mem_tether`` bonds whose detachment load exceeds the adhesion. Rupture
  = bleb nucleation. (Re-attachment / bleb retraction by ERM re-binding is a
  documented TODO — the first version is rupture-only, matching the bleb-growth
  experiments where nucleation is the studied event.)

The bleb-nucleation THRESHOLD and the reservoir EXCESS-AREA FRACTION are uncertain
for MCF7 (see PI_DECISIONS): both are set to ``None`` by default and the code
paths that need them raise :class:`NotImplementedError` (the rupture updater and
the reservoir release are kept DISABLED until a PI-anchored value is supplied).
The breakable-tether TOPOLOGY itself is provided with the cited membrane–cortex
adhesion energy where a literature value exists (the MCA band, Dai & Sheetz 1999;
Diz-Muñoz et al. 2010), so the static, force-free tether mesh is buildable now.

Adhesion energy → rupture force bridge (no invented number)
-----------------------------------------------------------
The membrane–cortex adhesion is an energy per unit area ``W_MCA`` [J/m²]
(KU-3.B1.3 band 1e-6 .. 1e-4 J/m²; Dai & Sheetz 1999; Diz-Muñoz 2010). Each
membrane bead carries an area share ``A_bead = 4π R_cell² / n_anchored``; its
tether holds an adhesion energy ``E_tether = W_MCA · A_bead`` [J]. The tether is
a harmonic bond of stiffness ``k_tether`` and rest length ``r0`` (the
construction membrane↔cortex separation). The bond ruptures when the elastic
energy stored by detachment exceeds the adhesion energy, i.e. at a critical
extension ``Δ_c`` with ``½ k_tether Δ_c² = E_tether`` → the critical rupture
FORCE is::

    F_c = k_tether · Δ_c = √(2 k_tether · W_MCA · A_bead)            [N]      (R)

Dimensions: √([N/m]·[J/m²]·[m²]) = √([N/m]·[J]) = √([N/m]·[N·m]) = √(N²) = N. ✓
``F_c`` is grid-aware: ``A_bead ∝ 1/n_anchored`` so the per-bead rupture force
falls as the membrane is discretised finer, while the TOTAL adhesion energy
``Σ E_tether = W_MCA · 4π R_cell²`` is N-independent (intensive material adhesion)
— exactly the per-bead area-share bridge ``cortex/enclosed_volume.py`` and
``cell/membrane_surface.py`` use. ``F_c`` is NOT tuned to pass a gate; it is the
adhesion energy expressed as a rupture load.

The Tinevez (2009) **critical-tension** ``σ_crit`` (the cortical tension below
which a bleb cannot expand) is the COMPLEMENTARY, cell-scale threshold; its MCF7
value is uncertain (PI decision) and is carried as ``None`` — the rupture updater
that would use it stays disabled. The static tether topology only needs ``W_MCA``.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_membrane_reservoir.py``.*

1. **Dimensional analysis**
   - ``W_MCA`` [J/m²]; ``A_bead = 4π R²/n`` [m²]; ``E_tether = W_MCA·A_bead`` [J]. ✓
   - ``k_tether`` [N/m]; rupture force ``F_c = √(2 k_tether W_MCA A_bead)``
     = √([N/m]·[J]) = [N]. ✓  Critical extension ``Δ_c = F_c/k_tether`` [m]. ✓
   - Reservoir: ``f_excess`` dimensionless; released area
     ``A_release = f_excess · A0`` [m²] grows the elastic reference area. ✓
   - At construction (membrane bead at its cortex anchor, ``Δr = r0``): tether
     force = 0 and ``U = 0`` — force-free tether mesh (Sanity Gate boundary).
2. **Boundary cases**
   - ``enabled=False`` (or missing): resolver returns ``ResolvedMembraneReservoir``
     with ``.enabled=False`` and zeros; the topology builder EARLY-RETURNS the
     input snapshot unchanged (no particles, no bonds added). STATIC test.
   - ``n_anchored = 0`` (no membrane beads in tether range, or no cortex acceptor
     within ``max_tether_dist``): builder adds zero ``mem_tether`` bonds; trivially
     correct. STATIC test.
   - ``W_MCA``, ``k_tether``, ``R_cell``, ``max_tether_dist`` ≤ 0 or non-finite:
     :func:`resolve_membrane_reservoir` raises ``ValueError``. STATIC test.
   - ``sigma_crit_bleb`` is ``None`` (MCF7 unknown — PI decision): the rupture
     updater's ``__init__`` raises ``NotImplementedError`` (bleb path disabled),
     never invents a threshold. STATIC test.
   - ``f_excess`` is ``None`` (MCF7 unknown — PI decision): :meth:`released_area`
     raises ``NotImplementedError`` (reservoir release disabled). STATIC test.
3. **Conservation / no-net-force**
   - ``mem_tether`` is a symmetric ``md.bond.Harmonic`` pair force (Newton's 3rd
     law) → it injects ZERO net momentum on the membrane+cortex system. STATIC
     reasoning; a tether between bead pairs is force-symmetric by construction.
   - The tether load path is RADIAL detachment (membrane vs cortex), NOT in-plane
     cortical tension; named ``mem_*`` so the cortical-tension denylist drops it
     (no γ contamination). STATIC test asserts the prefix + denylist match.
   - Particle count is INVARIANT (default build reuses existing beads); only the
     ``mem_tether`` bond family is appended. STATIC test asserts ΔN_particles = 0.
4. **Numerical sanity (CFL)**
   - Positions / energies float64. The tether harmonic has stiffness ``k_tether``;
     its relax time ``τ = γ_b / k_tether`` gates ``dt ≤ cfl_safety_factor·τ`` — the
     attach helper raises like ``erm.py`` / ``crosslinkers.py`` if violated.
   - The rupture updater is a per-BATCH Action; its batch CFL
     ``batch_steps·dt·k_off_max ≤ 1e-3`` mirrors the D2 crosslinker contract — but
     it is only constructible once ``sigma_crit_bleb`` is PI-anchored (disabled).
5. **Sign / sense**
   - A STRETCHED tether (membrane bead pulled OUTWARD from cortex, ``Δr > r0``)
     pulls the membrane bead back INWARD toward its cortex anchor (restoring) —
     the adhesion resists detachment. STATIC test on the harmonic force sign.
   - Larger adhesion ``W_MCA`` → larger rupture force ``F_c`` (harder to bleb):
     ``F_c`` monotonically INCREASES with ``W_MCA``. STATIC test.
6. **Measurement-protocol consistency**
   - The MCA adhesion ``W_MCA`` is the same energy density the membrane tether-
     force oracle of ``membrane_surface.py`` uses (γ_MCA); the rupture-force bridge
     (eq. R) is the per-bead expression of that adhesion. NO claim that a bleb gate
     PASSES — this supplies the explicit breakable topology + the adhesion anchor.

Compartment Performance Contract
--------------------------------
* Particle types added: NONE (default build reuses membrane + cortex beads). A
  future own-layer variant would add ``mem_bead`` (TODO, not built here).
* Particle count added @ n_fil=1000 mesoscale: 0. @ native ~38000 scale: 0.
* Bond/angle count: ``n_anchored`` ``mem_tether`` bonds (≤ one per membrane bead;
  ~ n_cortex_beads, so O(7e3) @ n_fil=1000, O(2.6e5) @ native). 0 angles.
* Per-step force: NO. The tether is a standard ``md.bond.Harmonic`` evaluated by
  HOOMD's native bonded ForceCompute (GPU-resident, builtin) — no Python per-step
  custom force.
* Per-batch updater: YES (rupture = bleb nucleation) — but DISABLED until
  ``sigma_crit_bleb`` is PI-anchored; when enabled it runs every ``batch_steps``
  like the D2 crosslinker / H.4 integrin updaters.
* Uses cpu_local_snapshot: NO (builder is a pure snapshot transform; the rupture
  updater, when enabled, uses ``sim.state.get_snapshot()`` like the D2 updater).
* Uses cKDTree / broad-phase: YES, ONCE at construction — nearest-cortex-bead
  acceptor query for the tether mesh (O(n_mem · log n_cortex)); none in the hot
  loop (the tether pairs are fixed once built).
* Hot-path priority: **P2** (per the task spec). The static tether mesh is a
  builtin bonded force (cheap); the rupture updater is batched and rare.
* GPU path now: **builtin** for the tether force (HOOMD native ``bond.Harmonic``);
  the (disabled) rupture updater would be CPU like the other batched binders until
  the GPU-main binder port lands.
* Native ForceCompute candidate: N/A for the tether (already builtin). The rupture
  updater is a candidate for the shared CUDA dynamic-bond plugin (future).
* Bottleneck risk: LOW. Construction KDTree is one-shot; the per-step cost is a
  native bonded force over ~n_cortex bonds, negligible vs the BAOAB integrator.

References
----------
- Charras, G. & Paluch, E. (2008) "Blebs lead the way: how to migrate without
  lamellipodia." Nat. Rev. Mol. Cell Biol. 9:730–736. doi:10.1038/nrm2453
  (bleb nucleation when hydrostatic pressure exceeds membrane–cortex adhesion;
  membrane reservoir of folds/microvilli).
- Tinevez, J.-Y. et al. (2009) "Role of cortical tension in bleb growth."
  Proc. Natl. Acad. Sci. USA 106:18581–18586. doi:10.1073/pnas.0903353106
  (a CRITICAL cortical tension below which blebs cannot expand — the σ_crit
  threshold; value cell-line specific, MCF7 uncertain → PI decision).
- Dai, J. & Sheetz, M.P. (1999) "Membrane tether formation from blebbing cells."
  Biophys. J. 77:3363–3370. doi:10.1016/S0006-3495(99)77168-7 (membrane–cortex
  adhesion energy from tether pulling; W_MCA in the 1e-6..1e-4 J/m² band).
- Diz-Muñoz, A. et al. (2010) "Control of directed cell migration in vivo by
  membrane-to-cortex attachment." PLoS Biol. 8:e1000544.
  doi:10.1371/journal.pbio.1000544 (membrane-to-cortex attachment / MCA energy).
- Raucher, D. & Sheetz, M.P. (1999) "Characteristics of a membrane reservoir
  buffering membrane tension." Biophys. J. 77:1992–2002.
  doi:10.1016/S0006-3495(99)77040-2 (the membrane area reservoir that buffers
  tension).
- Figard, L. & Sokac, A.M. (2014) "A membrane reservoir at the cell surface:
  unfolding the plasma membrane to fuel cell shape change." Bioarchitecture
  4:39–46. doi:10.4161/bioa.29069 (excess-area reservoir, a few % to tens of %;
  MCF7 value uncertain → PI decision).
- KU-3.B1.3 membrane–cortex adhesion energy band 1e-6 .. 1e-4 J/m² (see
  ``cell/membrane_surface.py`` GAMMA_MCA_BAND; same anchor).
- Structural analogs (in-tree templates): ``cortex/erm.py`` (radial membrane–
  cortex pinning + Sanity Gate + CFL attach gate), ``cortex/crosslinkers.py``
  (breakable Bell-Evans dynamic-bond family + per-batch Action + bond binning),
  ``cell/membrane_surface.py`` (the REQUIRED surface module: shell beads + A0),
  ``cell/nucleus.py`` (per-bead area-share bridge), ``cell/cell.py``
  SubstrateLigandPin (custom Action) + the FA snapshot bond-extension pattern.
- HOOMD 7 ``md.bond.Harmonic`` / ``hoomd.custom.Action`` APIs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd


# ---------------------------------------------------------------------------
# GAMMA denylist + open PI decisions (module-level constants per task spec)
# ---------------------------------------------------------------------------
#: Every bond type this module creates starts with this prefix so the cortical-
#: tension estimator (``ffn_sim.cortex.cortical_tension``) can DENYLIST it — the
#: membrane–cortex tether is a RADIAL detachment load path, NOT in-plane cortical
#: tension, and must never contaminate the γ_soft method-of-planes sum.
GAMMA_DENYLIST_PREFIX: str = "mem_"

#: Open, PI-gated decisions (uncertain MCF7 constants kept ``None`` + disabled).
PI_DECISIONS: list[str] = [
    "sigma_crit_bleb (Tinevez 2009 critical cortical tension for bleb growth) is "
    "cell-line specific; the MCF7 value is uncertain. Set to None; the Bell-Evans "
    "bleb-rupture updater (MembraneTetherUpdater) raises NotImplementedError until "
    "PI anchors a value. Do NOT invent a threshold to make blebs nucleate.",
    "f_excess (membrane reservoir excess-area fraction; Raucher-Sheetz 1999 / "
    "Figard 2014 report a few % to tens of %) is uncertain for MCF7. Set to None; "
    "released_area() raises NotImplementedError so the reservoir-release buffering "
    "stays disabled until PI anchors the MCF7 excess-area fraction.",
]

#: Membrane–cortex adhesion energy band, KU-3.B1.3 (Dai & Sheetz 1999; Diz-Muñoz
#: 2010). Same anchor as ``membrane_surface.py`` GAMMA_MCA_BAND. [J/m²]
W_MCA_BAND: tuple[float, float] = (1.0e-6, 1.0e-4)
#: Default membrane–cortex adhesion energy: mid-band 1e-5 J/m² (cited, in-band).
DEFAULT_W_MCA: float = 1.0e-5  # J/m²  KU-3.B1.3 (Dai & Sheetz 1999)

#: Canonical breakable tether bond-type family name (starts with the denylist
#: prefix so it is excluded from γ_soft).
MEM_TETHER_BOND: str = "mem_tether"


def _require_finite_positive(name: str, x: float) -> None:
    """Raise ValueError unless ``x`` is finite and strictly positive."""
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedMembraneReservoir:
    """Membrane reservoir + bleb-tether parameters (all SI).

    Attributes:
        enabled: Master switch. When False this module is a no-op: the topology
            builder early-returns the input snapshot unchanged.
        W_MCA: Membrane–cortex adhesion energy density [J/m²]. KU-3.B1.3 band
            1e-6..1e-4 (Dai & Sheetz 1999; Diz-Muñoz 2010). Sets the tether
            adhesion energy / rupture force (eq. R).
        k_tether: Membrane–cortex tether harmonic stiffness [N/m]. The same
            order as the ERM radial pinning ``k_ERM`` (KU-3.18, 0.1 N/m); the
            tether is the per-bead bond form of that membrane–cortex link.
        R_cell: Nominal shell radius [m] (sets the per-bead area share
            ``A_bead = 4π R_cell²/n_anchored``).
        max_tether_dist: Acceptor search radius [m] — a membrane bead tethers to
            a cortex bead within this distance. Default = the tether construction
            reach; no membrane bead farther than this is anchored.
        sigma_crit_bleb: Tinevez 2009 critical cortical tension for bleb growth
            [N/m]. ``None`` for MCF7 (PI decision) → the rupture updater is
            disabled (raises NotImplementedError). NEVER invented.
        f_excess: Reservoir excess-area fraction [dimensionless]. ``None`` for
            MCF7 (PI decision) → reservoir release is disabled (raises
            NotImplementedError). NEVER invented.
        batch_steps: BAOAB steps between rupture-updater ticks (D2 batch
            convention). Only used when the rupture path is PI-enabled.
        dt: Host sim timestep [s] (for the batch CFL when the updater is enabled).
        seed: Deterministic RNG seed for the (PI-gated) rupture updater.
    """

    enabled: bool
    # Adhesion / tether mechanics.
    W_MCA: float                 # J/m²  membrane–cortex adhesion (KU-3.B1.3)
    k_tether: float              # N/m   tether harmonic stiffness
    R_cell: float                # m     shell radius (per-bead area share)
    max_tether_dist: float       # m     acceptor search radius
    # PI-gated uncertain MCF7 constants — None keeps their code path disabled.
    sigma_crit_bleb: float | None = None   # N/m  Tinevez 2009 σ_crit (MCF7 ?)
    f_excess: float | None = None          # —    reservoir excess area (MCF7 ?)
    # Batch (rupture updater; only used when PI-enabled).
    batch_steps: int = 100
    dt: float = 0.0
    seed: int = 42
    # Derived diagnostics (per-bead, populated for n_anchored if known).
    A_bead_at_R: float = 0.0     # m²   4π R_cell² (n_anchored=1 reference share)
    extras: dict[str, Any] = field(default_factory=dict)

    # ---- derived helpers (pure; safe to call on enabled config) ----
    def area_share(self, n_anchored: int) -> float:
        """Per-bead membrane area share ``A_bead = 4π R_cell²/n_anchored`` [m²]."""
        if int(n_anchored) < 1:
            raise ValueError(f"n_anchored must be ≥ 1; got {n_anchored!r}")
        return 4.0 * math.pi * self.R_cell * self.R_cell / float(n_anchored)

    def tether_adhesion_energy(self, n_anchored: int) -> float:
        """Per-tether adhesion energy ``E = W_MCA · A_bead`` [J]."""
        return self.W_MCA * self.area_share(n_anchored)

    def rupture_force(self, n_anchored: int) -> float:
        """Critical tether rupture force ``F_c = √(2 k_tether W_MCA A_bead)`` [N].

        Eq. (R): the adhesion energy expressed as a harmonic rupture load. Grid-
        aware (``A_bead ∝ 1/n_anchored``); the TOTAL adhesion
        ``Σ E = W_MCA·4π R_cell²`` is N-independent.
        """
        return math.sqrt(
            2.0 * self.k_tether * self.tether_adhesion_energy(n_anchored)
        )

    def rupture_extension(self, n_anchored: int) -> float:
        """Critical tether extension ``Δ_c = F_c / k_tether`` [m]."""
        return self.rupture_force(n_anchored) / self.k_tether

    def released_area(self, A0: float) -> float:
        """Reservoir area released to grow the reference area ``A0`` [m²].

        ``A_release = f_excess · A0``. DISABLED until ``f_excess`` is PI-anchored
        for MCF7 (see PI_DECISIONS) — raises :class:`NotImplementedError` so the
        reservoir buffering is never run from an invented excess-area fraction.
        """
        if self.f_excess is None:
            raise NotImplementedError(
                "membrane reservoir release is disabled: f_excess (excess-area "
                "fraction) is uncertain for MCF7 and intentionally None (see "
                "PI_DECISIONS). Surface to PI to anchor the MCF7 excess-area "
                "fraction (Raucher-Sheetz 1999 / Figard 2014 report a few %–"
                "tens of %) before enabling reservoir buffering."
            )
        _require_finite_positive("A0", A0)
        return self.f_excess * float(A0)


def resolve_membrane_reservoir(
    cfg: dict,
    *,
    R_cell: float,
    dt: float = 0.0,
) -> ResolvedMembraneReservoir:
    """Resolve the ``membrane_reservoir`` config block.

    ``cfg`` may be the YAML root, a ``cell:`` sub-dict, or the
    ``membrane_reservoir`` sub-dict. ``R_cell`` (shell radius) comes from the
    host cortex / membrane-surface resolved params. ``dt`` is the host sim
    timestep (only needed when the PI-gated rupture updater is later enabled).

    DEFAULT-OFF: when ``cfg`` has ``enabled: false`` (or the key is missing),
    returns ``ResolvedMembraneReservoir(enabled=False, ...)`` with zeros — the
    topology builder then early-returns its input unchanged.

    The uncertain MCF7 constants ``sigma_crit_bleb`` and ``f_excess`` default to
    ``None`` (PI decisions); supplying them in ``cfg`` is allowed (and enables the
    corresponding path) but no default is invented.

    Args:
        cfg: Config mapping (root, ``cell``, or ``membrane_reservoir`` sub-dict).
        R_cell: Shell radius [m] (host geometry; REQUIRED, > 0).
        dt: Host sim timestep [s] (for the rupture-updater batch CFL; optional).

    Returns:
        Resolved membrane-reservoir parameters (SI).

    Raises:
        ValueError: on a non-positive / non-finite enabled parameter, or an
            out-of-band ``W_MCA`` (Sanity Gate §2).
    """
    if "cell" in cfg:
        cfg = cfg["cell"]
    if "membrane_reservoir" in cfg:
        cfg = cfg["membrane_reservoir"]

    enabled = bool(cfg.get("enabled", False))

    if not enabled:
        # DEFAULT-OFF identity: zeros, both PI-gated constants None.
        return ResolvedMembraneReservoir(
            enabled=False,
            W_MCA=0.0,
            k_tether=0.0,
            R_cell=0.0,
            max_tether_dist=0.0,
            sigma_crit_bleb=None,
            f_excess=None,
            batch_steps=int(cfg.get("batch_steps", 100)),
            dt=0.0,
            seed=int(cfg.get("seed", 42)),
            A_bead_at_R=0.0,
        )

    # ---- enabled: read + validate parameters ----
    W_MCA = float(cfg.get("W_MCA", DEFAULT_W_MCA))
    # k_tether order = ERM radial pinning (KU-3.18, 0.1 N/m) — the per-bead bond
    # form of the membrane–cortex link; a caller may override within reason.
    k_tether = float(cfg.get("k_tether", 0.1))
    # Acceptor reach: default = the mesoscale ~ membrane-cortex gap scale. The
    # cortex is ~200 nm thick (KU-3.17); the membrane sits within a tether length
    # of it, so 200 nm is the physical acceptor reach (NOT a tuned bind_scale).
    max_tether_dist = float(cfg.get("max_tether_dist", 200.0e-9))

    # Optional PI-anchored uncertain constants (default None — never invented).
    sigma_crit_bleb = cfg.get("sigma_crit_bleb", None)
    sigma_crit_bleb = (
        None if sigma_crit_bleb is None else float(sigma_crit_bleb)
    )
    f_excess = cfg.get("f_excess", None)
    f_excess = None if f_excess is None else float(f_excess)

    batch_steps = int(cfg.get("batch_steps", 100))
    seed = int(cfg.get("seed", 42))

    # §2 boundary checks (enabled path).
    _require_finite_positive("R_cell", R_cell)
    _require_finite_positive("W_MCA", W_MCA)
    _require_finite_positive("k_tether", k_tether)
    _require_finite_positive("max_tether_dist", max_tether_dist)
    if batch_steps < 1:
        raise ValueError(f"batch_steps must be ≥ 1; got {batch_steps!r}")
    if not (math.isfinite(dt) and dt >= 0.0):
        raise ValueError(f"dt must be finite and ≥ 0; got {dt!r}")

    # KU-3.B1.3 adhesion-energy band guard (no magic number / tuned override).
    lo, hi = W_MCA_BAND
    if not (lo <= W_MCA <= hi):
        raise ValueError(
            f"W_MCA = {W_MCA:.3e} J/m² outside the KU-3.B1.3 membrane–cortex "
            f"adhesion band [{lo:.1e}, {hi:.1e}] J/m² (Dai & Sheetz 1999; "
            "Diz-Muñoz 2010). Per CLAUDE.md no-magic-number: keep within band "
            "or surface to PI."
        )
    # Optional anchored constants must be physical when supplied.
    if sigma_crit_bleb is not None:
        _require_finite_positive("sigma_crit_bleb", sigma_crit_bleb)
    if f_excess is not None:
        if not (math.isfinite(f_excess) and f_excess >= 0.0):
            raise ValueError(
                f"f_excess must be finite and ≥ 0; got {f_excess!r}"
            )

    return ResolvedMembraneReservoir(
        enabled=True,
        W_MCA=W_MCA,
        k_tether=k_tether,
        R_cell=float(R_cell),
        max_tether_dist=max_tether_dist,
        sigma_crit_bleb=sigma_crit_bleb,
        f_excess=f_excess,
        batch_steps=batch_steps,
        dt=float(dt),
        seed=seed,
        A_bead_at_R=4.0 * math.pi * float(R_cell) * float(R_cell),
    )


# ---------------------------------------------------------------------------
# Tether topology builder (snapshot extension — pure transform)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class MembraneTetherLayout:
    """Bookkeeping for the seeded membrane–cortex tether mesh.

    Attributes:
        tether_pairs: ``(n_tether, 2)`` int64 — ``(membrane_bead_tag,
            cortex_bead_tag)`` per tether (force-free at construction).
        tether_r0: ``(n_tether,)`` float64 — per-tether rest length [m] (the
            construction membrane↔cortex separation).
        n_tether: Number of seeded tethers.
        rupture_force: Per-tether critical rupture force ``F_c`` [N] (eq. R; the
            same for every tether given a fixed ``n_anchored`` area share).
    """

    tether_pairs: np.ndarray
    tether_r0: np.ndarray
    n_tether: int
    rupture_force: float


def build_membrane_tethers(
    snap: Any,
    p: ResolvedMembraneReservoir,
    *,
    membrane_tag_range: tuple[int, int],
    cortex_tag_range: tuple[int, int],
) -> tuple[Any, MembraneTetherLayout]:
    """Append explicit breakable ``mem_tether`` bonds to a snapshot/frame.

    For each membrane bead (tags in ``membrane_tag_range``) the nearest cortex
    bead (tags in ``cortex_tag_range``) within ``p.max_tether_dist`` is found via
    a one-shot KDTree and joined by a ``mem_tether`` harmonic bond at the
    construction separation (force-free). Rupture of such a bond IS bleb
    nucleation (handled by :class:`MembraneTetherUpdater`, PI-gated).

    DEFAULT-OFF / no-op contract: if ``p.enabled`` is False, OR no membrane bead
    has a cortex acceptor in range, the snapshot is returned UNCHANGED (no
    particles, no bonds added) with an empty layout. This builder NEVER mutates
    a live simulation/state — it is a pure snapshot→snapshot transform (the Lead
    sets the resulting snapshot, single-writer convention).

    Args:
        snap: A ``gsd.hoomd.Frame`` or HOOMD ``Snapshot`` with ``particles``
            (position, tag/typeid) and ``bonds`` (types/typeid/group/N). Treated
            read-mostly; a NEW frame is returned (input not mutated when bonds are
            added; the SAME object is returned unchanged on the no-op path).
        p: Resolved membrane-reservoir parameters.
        membrane_tag_range: ``[start, end)`` tag range of membrane beads to tether.
        cortex_tag_range: ``[start, end)`` tag range of cortex acceptor beads.

    Returns:
        ``(out_snap, layout)``. On the no-op path ``out_snap is snap`` and
        ``layout.n_tether == 0``.

    Raises:
        ValueError: on an inverted tag range.
    """
    m_lo, m_hi = int(membrane_tag_range[0]), int(membrane_tag_range[1])
    c_lo, c_hi = int(cortex_tag_range[0]), int(cortex_tag_range[1])
    if m_hi < m_lo or c_hi < c_lo:
        raise ValueError(
            "tag ranges must satisfy end ≥ start; got membrane "
            f"({m_lo}, {m_hi}), cortex ({c_lo}, {c_hi})."
        )

    empty_layout = MembraneTetherLayout(
        tether_pairs=np.empty((0, 2), dtype=np.int64),
        tether_r0=np.empty((0,), dtype=np.float64),
        n_tether=0,
        rupture_force=0.0,
    )

    # ---- DEFAULT-OFF identity: return the input unchanged ----
    if not p.enabled:
        return snap, empty_layout

    pos = np.asarray(snap.particles.position, dtype=np.float64).reshape(-1, 3)
    n_part = pos.shape[0]
    tags = np.arange(n_part, dtype=np.int64)  # gsd/HOOMD: row index == tag
    mem_mask = (tags >= m_lo) & (tags < m_hi)
    cor_mask = (tags >= c_lo) & (tags < c_hi)
    mem_idx = np.flatnonzero(mem_mask)
    cor_idx = np.flatnonzero(cor_mask)

    if mem_idx.size == 0 or cor_idx.size == 0:
        # No beads in range → no tethers; snapshot unchanged.
        return snap, empty_layout

    # One-shot nearest-cortex-acceptor query (construction only; no hot loop).
    from scipy.spatial import cKDTree

    tree = cKDTree(pos[cor_idx])
    d, j = tree.query(pos[mem_idx], k=1, distance_upper_bound=p.max_tether_dist)
    in_range = np.isfinite(d) & (d <= p.max_tether_dist)
    if not in_range.any():
        return snap, empty_layout

    mem_anchored = mem_idx[in_range]
    cor_anchored = cor_idx[j[in_range]]
    r0 = d[in_range].astype(np.float64)
    n_tether = int(mem_anchored.size)
    n_anchored = n_tether
    rupture_force = p.rupture_force(n_anchored)

    tether_pairs = np.stack([mem_anchored, cor_anchored], axis=1).astype(np.int64)

    # ---- append mem_tether bond family to a NEW frame ----
    old_types = list(snap.bonds.types) if snap.bonds.N or snap.bonds.types else []
    old_N = int(snap.bonds.N)
    old_group = (
        np.asarray(snap.bonds.group, dtype=np.int64).reshape(old_N, 2)
        if old_N > 0
        else np.empty((0, 2), dtype=np.int64)
    )
    old_typeid = (
        np.asarray(snap.bonds.typeid, dtype=np.uint32)
        if old_N > 0
        else np.empty((0,), dtype=np.uint32)
    )

    new_types = list(old_types)
    if MEM_TETHER_BOND not in new_types:
        new_types.append(MEM_TETHER_BOND)
    tether_typeid = np.uint32(new_types.index(MEM_TETHER_BOND))

    new_group = np.concatenate(
        [old_group, tether_pairs], axis=0
    ).astype(np.uint32)
    new_typeid = np.concatenate(
        [old_typeid, np.full(n_tether, tether_typeid, dtype=np.uint32)]
    )

    out = _clone_frame_with_bonds(snap, new_types, new_typeid, new_group)

    layout = MembraneTetherLayout(
        tether_pairs=tether_pairs,
        tether_r0=r0,
        n_tether=n_tether,
        rupture_force=rupture_force,
    )
    return out, layout


def _clone_frame_with_bonds(
    snap: Any,
    bond_types: list[str],
    bond_typeid: np.ndarray,
    bond_group: np.ndarray,
) -> Any:
    """Return a new ``gsd.hoomd.Frame`` copy of ``snap`` with new bonds.

    Particles, angles and box are passed through unchanged; only the bond block
    is replaced. Mirrors the FA snapshot bond-extension pattern in
    ``cell/cell.py`` (a fresh frame, single-writer). Used only on the enabled,
    has-tethers path (the no-op path returns the input object directly).
    """
    import gsd.hoomd

    out = gsd.hoomd.Frame()
    out.particles.N = int(np.asarray(snap.particles.position).reshape(-1, 3).shape[0])
    out.particles.types = list(snap.particles.types)
    out.particles.typeid = np.asarray(snap.particles.typeid, dtype=np.uint32)
    out.particles.position = np.asarray(
        snap.particles.position, dtype=np.float64
    ).reshape(-1, 3)
    if getattr(snap.particles, "mass", None) is not None:
        try:
            out.particles.mass = np.asarray(snap.particles.mass, dtype=np.float64)
        except (TypeError, ValueError):
            pass

    out.bonds.N = int(bond_group.shape[0])
    out.bonds.types = list(bond_types)
    out.bonds.typeid = bond_typeid.astype(np.uint32)
    out.bonds.group = bond_group.astype(np.uint32)

    n_ang = int(getattr(snap.angles, "N", 0) or 0)
    if n_ang > 0:
        out.angles.N = n_ang
        out.angles.types = list(snap.angles.types)
        out.angles.typeid = np.asarray(snap.angles.typeid, dtype=np.uint32)
        out.angles.group = np.asarray(snap.angles.group, dtype=np.uint32)

    out.configuration.box = list(snap.configuration.box)
    return out


# ---------------------------------------------------------------------------
# Attach the mem_tether harmonic force to a built simulation (builtin bond)
# ---------------------------------------------------------------------------
def attach_membrane_tether_force(
    sim: hoomd.Simulation,
    p: ResolvedMembraneReservoir,
    layout: MembraneTetherLayout,
    *,
    gamma_b: float | None = None,
    cfl_safety_factor: float = 0.1,
    cfl_strict: bool = True,
) -> Any | None:
    """Attach the ``mem_tether`` HOOMD builtin harmonic bond to an Integrator.

    The tether is a standard ``md.bond.Harmonic`` (native, GPU-resident) — NOT a
    Python per-step custom force. Per-bond rest length ``r0`` is the construction
    separation; stiffness ``k = p.k_tether``. The bond CFL ``dt ≤ α·γ_b/k_tether``
    is gated like ``erm.py`` / ``crosslinkers.py``.

    No-op when ``p.enabled`` is False OR ``layout.n_tether == 0``: returns None
    without touching the integrator.

    Args:
        sim: Built simulation with an Integrator (and the ``mem_tether`` bond type
            already present in its state, from :func:`build_membrane_tethers`).
        p: Resolved membrane-reservoir parameters.
        layout: The tether layout from :func:`build_membrane_tethers`.
        gamma_b: Per-bead Stokes drag [N·s/m] for the CFL gate (None → skip gate).
        cfl_safety_factor: Same convention as the bond/ERM CFL gates. Default 0.1.
        cfl_strict: Raise on CFL violation if True.

    Returns:
        The attached ``md.bond.Harmonic`` force, or None on the no-op path.

    Raises:
        RuntimeError: if no Integrator is set, or the CFL gate is violated.
    """
    import hoomd.md as md

    if not p.enabled or layout.n_tether == 0:
        return None

    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "sim.operations.integrator must be set before attaching the "
            "mem_tether bond force."
        )
    if gamma_b is not None:
        dt = float(ig.dt)
        tau = gamma_b / p.k_tether
        dt_cfl = cfl_safety_factor * tau
        if dt > dt_cfl and cfl_strict:
            raise RuntimeError(
                f"mem_tether CFL violated: dt = {dt:.3e} s > "
                f"{cfl_safety_factor:.2f}·τ = {dt_cfl:.3e} s "
                f"(τ = γ_b/k_tether = {tau:.3e} s). Reduce dt or soften "
                "k_tether (PI sign-off). Pass cfl_strict=False to skip."
            )

    bond = md.bond.Harmonic()
    # Construction-separation rest length; force-free at t=0. All seeded tethers
    # share the canonical MEM_TETHER_BOND type — use the median r0 as the single
    # type rest length (the per-bond construction spread is sub-nm; the tether is
    # a soft adhesion link, residual mismatch ≪ kT, same convention as erm.py).
    r0 = float(np.median(layout.tether_r0)) if layout.n_tether else 0.0
    bond.params[MEM_TETHER_BOND] = dict(k=p.k_tether, r0=r0)
    ig.forces.append(bond)
    return bond


# ---------------------------------------------------------------------------
# Bleb-rupture updater (Bell-Evans) — PI-GATED, disabled until σ_crit anchored
# ---------------------------------------------------------------------------
class MembraneTetherUpdater(hoomd.custom.Action):
    """Per-batch Bell-Evans rupture of ``mem_tether`` bonds (= bleb nucleation).

    A ``mem_tether`` ruptures when its detachment load exceeds the membrane–cortex
    adhesion; rupture is the EXPLICIT bleb-nucleation event (topology change, not a
    lumped switch). The off-rate is Bell-Evans against the rupture force
    ``F_c`` (eq. R) referenced to the Tinevez (2009) critical tension.

    DISABLED in this version: the MCF7 critical tension ``sigma_crit_bleb`` is
    uncertain (PI decision), so ``__init__`` raises :class:`NotImplementedError`
    when ``p.sigma_crit_bleb is None`` — the bleb path NEVER runs from an invented
    threshold. Supply a PI-anchored ``sigma_crit_bleb`` to enable it.

    Args:
        p: Resolved membrane-reservoir parameters.
        layout: The seeded tether layout.
        kT: Thermal energy [J] for the Bell-Evans exponent.
        seed_offset: RNG offset added to ``simulation.seed`` (default 5, to avoid
            collision with the BAOAB/xlink/integrin updaters at 0/1/2).

    Raises:
        NotImplementedError: when ``p.sigma_crit_bleb`` is None (MCF7 unknown).
    """

    def __init__(
        self,
        p: ResolvedMembraneReservoir,
        layout: MembraneTetherLayout,
        *,
        kT: float,
        seed_offset: int = 5,
    ) -> None:
        super().__init__()
        if p.sigma_crit_bleb is None:
            raise NotImplementedError(
                "MembraneTetherUpdater (bleb-rupture) is disabled: "
                "sigma_crit_bleb (Tinevez 2009 critical cortical tension for "
                "bleb growth) is uncertain for MCF7 and intentionally None (see "
                "PI_DECISIONS). Surface to PI to anchor the MCF7 critical tension "
                "before enabling the bleb-nucleation path — do not invent it."
            )
        self.p = p
        self.layout = layout
        self.kT = float(kT)
        self._rng = np.random.default_rng(p.seed + seed_offset)
        self._sim_ref: hoomd.Simulation | None = None

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        # Unreachable while disabled (__init__ raises). The rupture loop mirrors
        # the D2 crosslinker break path: read bonds, compute mem_tether load,
        # sample Bell-Evans p_break against F_c, drop ruptured tethers. Left as a
        # documented TODO contingent on the PI-anchored σ_crit.
        raise NotImplementedError(
            "MembraneTetherUpdater.act requires a PI-anchored sigma_crit_bleb."
        )
