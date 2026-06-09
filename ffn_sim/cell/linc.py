"""H.x LINC complex — explicit nesprin-SUN bridges across the nuclear envelope.

ADDITIVE, DEFAULT-OFF compartment module. The LINC (Linker of Nucleoskeleton
and Cytoskeleton) complex mechanically COUPLES the nuclear lamina to the
cytoplasmic cytoskeleton: a nesprin (KASH-domain) outer-membrane protein binds
SUN (inner-membrane) proteins through the perinuclear space, and its
cytoplasmic spectrin-repeat tail engages cortical actin / perinuclear actin
caps / stress fibers / microtubules. Force generated at the cortex (or a
substrate-anchored stress fiber) is transmitted through LINC into the lamina;
removing LINC mechanically decouples the nucleus (Lombardi 2011 JBC).

This module supplies the EXPLICIT mechanism per CLAUDE.md: each LINC bridge is
a real HOOMD **harmonic bond** (``linc_nesprin``) between one ``nucleus_bead``
particle tag and one cytoskeletal particle tag (cortex actin / stress-fiber /
MT bead). There is NO lumped coupling force and NO mesh-as-physics — a LINC
bridge is one bond carrying one molecular spring's load, exactly as a real
nesprin-2 giant spectrin-repeat tail is a molecular spring (Autore 2013 PLoS
ONE; Arsenovic 2016 Biophys J showed nesprin-2G bears myosin-dependent tension
by FRET — TSMod sensor, calibrated range ~1-5 pN — but reports relative FRET,
not an absolute pN value; the resting-tension magnitude ~8 pN comes from the
separate mini-nesprin-2G/CB construct of Déjardin 2020 JCB). The topology is
built geometrically: each nucleus bead bonds to its
nearest eligible cytoskeletal bead within a capture radius, producing
``n_bridges`` radial nesprin springs distributed over the envelope.

Requires the nucleus compartment (``cell/nucleus.py``) to be ON — the bonds
attach to ``nucleus_bead`` tags. With no nucleus present the resolver still
works but the builder finds zero acceptor pairs and is a no-op.

Mechanism (fine-grained per CLAUDE.md)
--------------------------------------
For each formed LINC bridge between nucleus-bead position ``r_n`` and
cytoskeleton-bead position ``r_c``::

    U   = ½ k_linc (|r_c − r_n| − r0)²
    F_n = +k_linc (|Δr| − r0) · Δr̂          (Δr̂ = (r_c − r_n)/|r_c − r_n|)
    F_c = −F_n                                (Newton 3rd law — built-in Harmonic)

implemented with HOOMD's **builtin** ``md.bond.Harmonic`` (NOT a custom force):
the integrator and the native constrained-BAOAB plugin both already handle
``md.bond.Harmonic`` on the GPU, so LINC adds zero per-step Python and zero
``cpu_local_snapshot`` cost — it is pure topology + a builtin potential.

The per-bond stiffness ``k_linc`` of a single nesprin spectrin-repeat spring is
genuinely UNCERTAIN at this coarse scale (see PI_DECISIONS): nesprin-2 giant
has ~56 spectrin repeats whose force-extension is highly nonlinear (repeats
unfold sequentially), so there is no single well-defined Hookean stiffness, and
the published single-molecule numbers (Autore 2013 mechanics; AFM repeat
unfolding ~25–35 pN plateaus) do not pin a linear ``N/m``. Therefore
``k_linc`` defaults to ``None`` and the enabled topology-build path RAISES
``NotImplementedError`` until the PI supplies a literature-anchored value. The
resting tension ``f_rest`` (~8 pN, Déjardin 2020 JCB) IS known and is exposed
as a provenance field / diagnostic, but a tension is not a stiffness and is not
used to invent one. (Note: Arsenovic 2016 showed nesprin-2G bears
myosin-dependent tension by FRET but reported only relative ratiometric FRET
with a TSMod sensor of calibrated range ~1-5 pN — it did NOT report an absolute
pN value; the ~8 pN resting magnitude is Déjardin 2020's CB/mini-nesprin-2G
result, not Arsenovic's.)

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_linc.py``.*

1. **Dimensional analysis**
   - ``k_linc`` [N/m]; ``(|Δr| − r0)`` [m] → ``F = k_linc·Δ`` [N]. ✓
   - ``U = ½ k_linc Δ²`` [J]. ✓
   - ``f_rest`` [N] is a force, NOT a stiffness — it is provenance only and is
     never divided by a length to fabricate ``k_linc`` (would need an extension
     scale that is itself unknown).
   - At ``|Δr| = r0`` (bond at rest length): ``F = 0``, ``U = 0``.

2. **Boundary cases**
   - No nucleus beads present, or no cytoskeleton bead within ``capture_radius``
     of any nucleus bead: zero bridges formed → builder returns the snapshot
     unchanged (no-op). RUNTIME pass.
   - ``enabled = False``: resolver returns ``.enabled = False`` with zeroed
     fields; the builder early-returns the input snapshot unchanged.
   - ``k_linc = None`` while enabled: the topology builder raises
     ``NotImplementedError`` (honest — no invented number). The resolver itself
     does NOT raise (resolving a config is allowed; only BUILDING the load-
     bearing bonds requires the stiffness).
   - ``k_linc ≤ 0`` (if a value is supplied), ``r0 ≤ 0``, ``capture_radius ≤
     0``, ``n_bridges_max < 0``: ``resolve_linc`` ValueError.

3. **Conservation invariants (no net force)**
   - Each LINC bond is a builtin ``md.bond.Harmonic`` pair → equal-and-opposite
     forces on the nucleus bead and the cytoskeleton bead (Newton 3rd law).
     The bond family therefore injects EXACTLY zero net momentum into the
     system (unlike the lab-frame-anchored ERM external field). The force flows
     between two explicit internal particles — a genuine load path, not a field.
   - Bridges only touch ``nucleus_bead`` ↔ cytoskeleton pairs in the supplied
     tag sets; no other particles gain bonds.

4. **Numerical sanity (CFL)**
   - All geometry in float64; pairing distances in metres.
   - The bond is stiff: relax time ``τ = γ_bead / k_linc`` (γ from the host
     cytoskeleton/nucleus drag) must satisfy ``dt ≤ cfl_safety_factor · τ``.
     :func:`linc_cfl_dt_max` returns the gate; the attach helper raises if a
     supplied ``dt`` violates it (same convention as erm.py / nucleus.py). The
     gate cannot be evaluated until ``k_linc`` is known (None → skipped with a
     warning that the build path is disabled anyway).

5. **Sign / sense**
   - A STRETCHED LINC bond (``|Δr| > r0``) pulls the nucleus bead TOWARD the
     cytoskeleton bead and vice-versa (restoring / attractive) — i.e. the
     nucleus is pulled OUTWARD toward a cortex anchor that sits outside it, and
     the cortex anchor is pulled INWARD toward the nucleus. STATIC test on
     :func:`linc_bond_force_along_axis` asserts the sign.
   - A COMPRESSED bond (``|Δr| < r0``) pushes them apart. Harmonic (not a
     cable) by construction — documented; a tension-only cable variant is a TODO
     (nesprin tails buckle under compression, carrying ~no compressive load).

6. **Measurement-protocol consistency**
   - The resting per-nesprin tension ``f_rest`` (~8 pN, Déjardin 2020 JCB) is
     the validation oracle: at equilibrium the mean LINC bond tension should sit
     near ``f_rest`` once ``k_linc`` and the pre-strain are set by the host. We
     expose ``f_rest`` and the per-bond tension formula but do NOT claim any
     gate PASSes here (that is a Lead integration / production concern). NOTE:
     Arsenovic 2016 measured myosin-dependent FRET on nesprin-2G (TSMod sensor,
     calibrated range ~1-5 pN) but reported relative FRET, not an absolute pN
     value — the ~8 pN magnitude is Déjardin 2020's CB/mini-nesprin-2G result.

Compartment Performance Contract
--------------------------------
* **Particle types added**: NONE. LINC reuses existing ``nucleus_bead`` and
  cytoskeleton particles; it only adds bonds. (No new particle cloud — this is
  the key difference from nucleus.py.)
* **Particle count** at n_fil=1000 / native ~38000: +0 particles at both
  scales. LINC is topology-only.
* **Bond/angle count**: +``n_bridges`` bonds (one ``linc_nesprin`` bond per
  formed bridge). No angles. The GEOMETRIC ceiling on ``n_bridges`` is the
  nucleus bead count ``n_nuc`` (one LINC per envelope bead via the k=1
  nearest-acceptor pairing); ``n_bridges_max`` is a SAFETY/no-cap sentinel that
  by default (1e6) never binds — it is present only so a host CAN cap explicitly,
  so in practice ``n_bridges = min(n_nuc, n_bridges_max) = n_nuc``. Biology:
  thousands of LINC per nucleus (Lombardi 2011); at the ×40 mesoscale this maps
  to O(n_nuc) (tens–hundreds) — set by the nucleus bead resolution, NOT by the
  1e6 sentinel. Zero per-bond runtime cost beyond the builtin harmonic
  evaluation.
* **Per-step force**: YES, but via the BUILTIN ``md.bond.Harmonic`` only — no
  Python ``md.force.Custom`` is attached. No per-step Python.
* **Per-batch updater**: NO. LINC bonds are STATIC topology in this minimal
  version (formed once at build). A force-dependent nesprin unfolding / LINC
  turnover updater (Bell-Evans-style) is a documented TODO, not in the hot loop
  now.
* **uses cpu_local_snapshot**: NO (topology built once on a CPU gsd snapshot at
  construction; no per-step host sync).
* **uses cKDTree/broad-phase**: YES, ONCE at build time only (the geometric
  nucleus↔cytoskeleton nearest-acceptor pairing uses ``scipy.spatial.cKDTree``;
  O(n_nuc · log n_cyto)). NOT in the hot loop.
* **hot-path priority**: P2 (per the task spec).
* **GPU path now**: builtin (``md.bond.Harmonic`` runs GPU-resident; the native
  constrained-BAOAB plugin already evaluates harmonic bonds).
* **native ForceCompute candidate**: NO — already a HOOMD-builtin bond; nothing
  to port.
* **bottleneck risk**: LOW. Build-time cKDTree pairing is one-off; runtime is a
  handful of extra builtin harmonic bonds. The only stiffness-driven cost is the
  CFL on ``k_linc`` (could force a smaller ``dt`` if LINC is stiffer than the
  cortex/nucleus springs — gated explicitly).

References
----------
- Lombardi M.L. et al. (2011) *J. Biol. Chem.* 286:26743 — LINC (nesprin-SUN)
  transmits cytoskeletal force to the nucleus; disruption mechanically
  decouples the nucleus.
- Autore F. et al. (2013) *PLoS ONE* 8:e63633 — nesprin-2 spectrin-repeat
  structure / mechanics (spectrin-repeat molecular spring; sequential
  repeat-unfolding nonlinearity → no single Hookean stiffness).
- Arsenovic P.T. et al. (2016) *Biophys. J.* 110:34 (DOI 10.1016/j.bpj.2015.11.014)
  — TSMod FRET tension sensor showing nesprin-2G bears myosin-dependent tension
  in living cells; reports RELATIVE ratiometric FRET (sensor calibrated range
  ~1-5 pN), NOT an absolute resting pN value.
- Déjardin T. et al. (2020) *J. Cell Biol.* 219(10):e201908036
  (DOI 10.1083/jcb.201908036) — mini-nesprin-2G / CB TSMod tension biosensor;
  nesprin-2G (CB construct) is held under ~8 pN of resting tension generated by
  the actomyosin + microtubule cytoskeletons, balanced by cell–cell adhesion.
  This is the source of the ``f_rest`` ~8 pN magnitude (not Arsenovic 2016).
- Crisp M. et al. (2006) *J. Cell Biol.* 172:41 — coining of LINC; SUN-KASH
  bridge across the nuclear envelope (LINC density: thousands per nucleus).
- ``ffn_sim/cell/nucleus.py`` — nucleus_bead cloud this module bonds to.
- ``ffn_sim/cortex/crosslinkers.py`` — dynamic-bond / cKDTree pairing analog
  (here STATIC topology, builtin harmonic, no per-step updater).
- HOOMD 7 ``md.bond.Harmonic`` API.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np

import hoomd

try:
    from scipy.spatial import cKDTree
except Exception:  # pragma: no cover - scipy is a hard dep of the env
    cKDTree = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Module-level contracts (read by the registry + the cortical-tension denylist)
# ---------------------------------------------------------------------------
#: Every bond type this module creates starts with this prefix so the cortical-
#: tension (γ) estimator can denylist the LINC nesprin load path (it is a
#: nucleus↔cytoskeleton bridge, NOT a cortical-tension bond).
GAMMA_DENYLIST_PREFIX: str = "linc_"

#: Canonical LINC bond type name (nesprin-SUN spectrin-repeat molecular spring).
LINC_BOND_NESPRIN: str = "linc_nesprin"


def linc_nesprin_exact_type_names(n: int) -> list[str]:
    """Per-bond EXACT-r0 ``linc_nesprin`` bond-type names (γ-denylisted).

    Mirrors the FA molecular-clutch convention (``fa_actin_clutch`` /
    ``fa_actin_clutch_b{i}`` in :func:`cell.cell.build_cortex_full_simulation`):
    one bond type PER formed LINC bridge, each carrying that bridge's EXACT
    as-built nucleus↔cytoskeleton separation as its rest length, so every bond
    is born FORCE-FREE at the resting geometry — even though the per-nesprin
    stiffness ``k_linc`` (route-A folded-rod secant ~1e-2 N/m) is ~100× stiffer
    than the soft IF crosslink, where a coarse per-r0 bin (cf.
    ``intermediate_filaments.if_crosslink_b{i}``) would leave hundreds of kT of
    spurious construction pre-stress. The first name is the canonical
    ``linc_nesprin``; the rest are ``linc_nesprin_b1 ..`` (every name keeps the
    ``linc_`` γ-denylist prefix). Used by :func:`compute_linc_layout` /
    :func:`extend_snapshot_with_linc_layout` for the integrated cell (Option A:
    nucleus → perinuclear IF cage), where the construction separations vary
    (CV ~24%) so a single ``r0`` cannot be force-free.
    """
    return [
        LINC_BOND_NESPRIN if i == 0 else f"{LINC_BOND_NESPRIN}_b{i}"
        for i in range(int(n))
    ]

#: Open PI decisions (CLAUDE.md no-magic-number protocol). Empty when none.
PI_DECISIONS: list[str] = [
    "k_linc (single-nesprin spectrin-repeat spring stiffness, N/m): UNKNOWN. "
    "Nesprin-2 giant is a ~56-repeat spectrin spring whose force-extension is "
    "nonlinear (sequential repeat unfolding ~25-35 pN plateaus; Autore 2013), "
    "so there is no single Hookean N/m in the literature. Resting tension is "
    "known (~8 pN, Déjardin 2020 JCB mini-nesprin-2G/CB sensor; Arsenovic 2016 "
    "showed myosin-dependent FRET but reported only relative ratiometric FRET, "
    "not an absolute pN value) but a tension is not a stiffness and is not "
    "divided by an unknown extension scale to fabricate one. Default None; the "
    "enabled topology builder raises NotImplementedError until the PI supplies "
    "a literature-anchored k_linc.",
]


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class ResolvedLINC:
    """Resolved LINC-complex parameters (all SI).

    Attributes:
        enabled: Master on/off. When False the builder is a strict no-op and
            every other field is zeroed / None.
        k_linc: Per-nesprin harmonic spring stiffness [N/m]. ``None`` when
            UNKNOWN (default) — the enabled topology build raises
            ``NotImplementedError``. See module ``PI_DECISIONS``.
        r0: Bond rest length [m] (nesprin-SUN bridge span across the perinuclear
            space + the cytoplasmic tail reach at rest). Default derived from
            the perinuclear gap ~50 nm (Crisp 2006) — a GEOMETRIC length, not a
            tuned force constant.
        capture_radius: Geometric pairing radius [m]: a nucleus bead bonds to
            the nearest cytoskeleton bead within this radius. No pair beyond it.
        n_bridges_max: SAFETY/no-cap sentinel on the number of LINC bonds
            formed, NOT the binding ceiling. The true geometric ceiling is the
            nucleus bead count ``n_nuc`` (one LINC per envelope bead, k=1
            nearest-acceptor pairing); ``n_bridges_max`` defaults to 1e6, which
            never binds — it exists only so a host CAN cap explicitly. Biology
            has thousands of LINC per nucleus (Lombardi 2011); the mesoscale
            count is set by ``n_nuc``, not by this sentinel.
        f_rest: Resting per-nesprin tension [N] (~8 pN, Déjardin 2020 JCB
            mini-nesprin-2G/CB sensor). PROVENANCE / diagnostic only — NOT used
            to derive ``k_linc``. (Arsenovic 2016 reported relative FRET, not an
            absolute pN value, so the magnitude is Déjardin 2020's, not his.)
        bond_type_name: Bond type name (always starts with ``linc_``).
    """

    enabled: bool
    k_linc: float | None
    r0: float
    capture_radius: float
    n_bridges_max: int
    f_rest: float
    bond_type_name: str = LINC_BOND_NESPRIN


# ---------------------------------------------------------------------------
# Literature anchors (Magic-Number Block)
# ---------------------------------------------------------------------------
#: Perinuclear-space gap across which the SUN-KASH bridge spans (Crisp 2006
#: J Cell Biol — the nuclear-envelope luminal space is ~30-50 nm). Used as the
#: GEOMETRIC default rest length (a length, derivable, not a tuned force knob).
_PERINUCLEAR_GAP_M: float = 50.0e-9          # m   (50 nm, Crisp 2006)

#: Resting per-nesprin tension on the mini-nesprin-2G / CB TSMod tension sensor
#: (Déjardin 2020 JCB, DOI 10.1083/jcb.201908036 — nesprin-2G/CB is held under
#: ~8 pN of resting tension generated by the actomyosin + microtubule
#: cytoskeletons, balanced by cell-cell adhesion). Provenance only.
#: NOTE: Arsenovic 2016 (DOI 10.1016/j.bpj.2015.11.014) measured myosin-dependent
#: FRET on nesprin-2G (TSMod, calibrated range ~1-5 pN) but reported only
#: RELATIVE ratiometric FRET, not an absolute pN value — the ~8 pN magnitude is
#: Déjardin 2020's, not Arsenovic's. ⚠️ This is an ORACLE-TARGET value: the
#: 2 pN → 8 pN re-anchor + re-attribution must be PI-signed-off (no-magic-number /
#: no-gate-loosening) and the KB SourceEvidence row refreshed before the LINC
#: gate is treated as live.
_NESPRIN_REST_TENSION_N: float = 8.0e-12     # N   (8 pN, Déjardin 2020 JCB)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def resolve_linc(
    cfg: dict,
    *,
    R_cell: float | None = None,
    R_nuc: float | None = None,
) -> ResolvedLINC:
    """Resolve a ``linc`` config block into SI parameters.

    Accepts the YAML root, a ``cell`` sub-dict, or a ``linc`` sub-dict. When
    ``enabled`` is False or absent, returns a zeroed ``ResolvedLINC`` with
    ``.enabled = False`` (strict DEFAULT-OFF — the builder then no-ops).

    The rest length ``r0`` defaults to the perinuclear-gap geometric length
    (~50 nm, Crisp 2006). The default ``capture_radius`` (3× that gap = 150 nm)
    is the ANATOMICALLY correct nesprin span across the perinuclear space, but
    it is NOT the bead-to-bead pairing scale in this mesoscale geometry: at the
    platform's bead representation the nucleus is a filled cloud out to
    R_nuc ≈ 0.25·R_cell and the cortex/SF acceptors sit near r ≈ R_cell, so the
    nearest nucleus↔cortex bead separation is R_cell − R_nuc ≈ 5.6 µm (≈37× the
    150 nm default) — ZERO LINC bonds form at construction with the default
    (the identical zero-bonds-at-real-geometry lesson the FA path already
    learned, cell.py:1040-1050). The literal load path therefore needs EITHER
    (a) a perinuclear/cap acceptor population seeded near the envelope, OR (b) a
    caller-supplied ``capture_radius`` on the R_cell − R_nuc scale (just as FA
    exposes ``fa_clutch_capture_radius``); when LINC is wired into
    ``build_baseline_cell`` the production ``capture_radius`` must come from
    cell geometry (or perinuclear acceptors), per the physiological-baseline
    rule, NOT this 3×-gap default. Keep both numbers and label which is which:
    150 nm = real nesprin span; R_cell − R_nuc = bead-pairing scale.
    ``k_linc`` is UNKNOWN by default (``None``) — resolving is allowed, but
    BUILDING the load-bearing bonds with ``k_linc is None`` raises
    ``NotImplementedError`` (no invented number).

    Args:
        cfg: Config mapping (root, ``cell``, or ``linc`` sub-dict).
        R_cell: Host cell radius [m] (optional; only used to sanity-bound the
            capture radius if provided).
        R_nuc: Host nuclear radius [m] (optional; reserved for a future
            radius-aware default rest length — unused in this minimal version).

    Returns:
        Resolved LINC parameters (SI). ``.enabled = False`` ⇒ all-zero / None.

    Raises:
        ValueError: on a non-positive / non-finite supplied parameter, or a
            capture radius larger than the host cell (a clearly wrong pairing
            scale) when ``R_cell`` is provided.
    """
    if isinstance(cfg, dict) and "cell" in cfg and "linc" not in cfg:
        cfg = cfg["cell"]
    if isinstance(cfg, dict) and "linc" in cfg:
        cfg = cfg["linc"]
    if not isinstance(cfg, dict):
        cfg = {}

    enabled = bool(cfg.get("enabled", False))

    if not enabled:
        # Strict DEFAULT-OFF: zeroed fields, None stiffness, .enabled False.
        return ResolvedLINC(
            enabled=False,
            k_linc=None,
            r0=0.0,
            capture_radius=0.0,
            n_bridges_max=0,
            f_rest=0.0,
            bond_type_name=LINC_BOND_NESPRIN,
        )

    # k_linc: None unless the PI supplies a literature-anchored value.
    k_raw = cfg.get("k_linc", None)
    k_linc: float | None = None if k_raw is None else float(k_raw)

    r0 = float(cfg.get("r0", _PERINUCLEAR_GAP_M))
    # Capture radius default: 3× the perinuclear gap = 150 nm. This is the
    # anatomically correct real nesprin span across the perinuclear space, NOT
    # the bead-to-bead pairing scale in this mesoscale geometry: at the platform
    # bead representation the nucleus↔cortex bead separation is R_cell - R_nuc
    # ≈ 5.6 µm (≈37× this default), so ZERO LINC bonds form at construction with
    # the default — see the resolve_linc docstring + the R_cell-vs-R_nuc guard
    # below + the FA lesson at cell.py:1040-1050. A caller wiring LINC at the
    # real op-point MUST pass a capture_radius on the R_cell - R_nuc scale (or
    # seed perinuclear acceptors). GEOMETRIC, grid-invariant (multiple of a length).
    capture_radius = float(cfg.get("capture_radius", 3.0 * _PERINUCLEAR_GAP_M))
    # Sentinel: effectively no cap. The true ceiling on n_bridges is n_nuc (one
    # LINC per envelope bead, k=1 nearest-acceptor pairing); this 1e6 default
    # never binds — it is here only so a host CAN cap the bridge COUNT explicitly.
    n_bridges_max = int(cfg.get("n_bridges_max", 1_000_000))
    f_rest = float(cfg.get("f_rest", _NESPRIN_REST_TENSION_N))

    # §2 boundary checks (only on values that ARE supplied / used).
    _require_finite_positive("r0", r0)
    _require_finite_positive("capture_radius", capture_radius)
    _require_finite_positive("f_rest", f_rest)
    if n_bridges_max < 0:
        raise ValueError(f"n_bridges_max must be ≥ 0; got {n_bridges_max!r}")
    if k_linc is not None:
        _require_finite_positive("k_linc", k_linc)
    if R_cell is not None:
        _require_finite_positive("R_cell", R_cell)
        if capture_radius > R_cell:
            raise ValueError(
                f"capture_radius = {capture_radius:.3e} m exceeds the cell "
                f"radius R_cell = {R_cell:.3e} m — LINC bridges the nuclear "
                "envelope (a short span), this pairing scale is wrong."
            )
        # Symmetric lower-bound geometry check: at this platform's bead
        # representation the nucleus is a filled cloud out to R_nuc and the
        # cortex/SF acceptors sit near R_cell, so the nearest nucleus↔cortex
        # bead separation is ~ (R_cell - R_nuc). If the capture radius is below
        # that gap AND the only acceptors are the cortex shell, ZERO bonds form
        # at construction (the latent geometry trap; mirrors the FA lesson at
        # cell.py:1040-1050). Warn loudly rather than silently no-op. We warn
        # (not raise) because a caller MAY seed perinuclear/cap acceptors near
        # the envelope that close the gap — the resolver cannot see the acceptor
        # classes (those come from the build-time tag sets), so a hard raise
        # would be wrong for that legitimate path.
        if R_nuc is not None:
            _require_finite_positive("R_nuc", R_nuc)
            nuc_cortex_gap = R_cell - R_nuc
            if nuc_cortex_gap > 0.0 and capture_radius < nuc_cortex_gap:
                warnings.warn(
                    f"LINC capture_radius = {capture_radius:.3e} m is below the "
                    f"nucleus↔cortex bead gap R_cell - R_nuc = "
                    f"{nuc_cortex_gap:.3e} m (ratio "
                    f"{nuc_cortex_gap / capture_radius:.1f}×). At this mesoscale "
                    "bead geometry ZERO LINC bonds will form at construction if "
                    "the only acceptor class is the cortex shell. Pass a "
                    "capture_radius on the R_cell - R_nuc scale, or seed "
                    "perinuclear/cap acceptors near the envelope (cf. "
                    "fa_clutch_capture_radius; FA lesson cell.py:1040-1050).",
                    stacklevel=2,
                )

    return ResolvedLINC(
        enabled=True,
        k_linc=k_linc,
        r0=r0,
        capture_radius=capture_radius,
        n_bridges_max=n_bridges_max,
        f_rest=f_rest,
        bond_type_name=LINC_BOND_NESPRIN,
    )


# ---------------------------------------------------------------------------
# Scalar force law (pure — testable without a sim)
# ---------------------------------------------------------------------------
def linc_bond_force_along_axis(k_linc: float, r0: float, sep: float) -> float:
    """Signed harmonic force projected on the nucleus→cytoskeleton axis [N].

    Convention: returns the force component felt by the NUCLEUS bead along the
    unit vector pointing FROM the nucleus bead TO the cytoskeleton bead
    (``Δr̂``). For a stretched bond (``sep > r0``) this is POSITIVE — the
    nucleus is pulled toward the cytoskeleton anchor (restoring / attractive).
    For a compressed bond (``sep < r0``) it is NEGATIVE — pushed apart.

    This mirrors the per-bond force a builtin ``md.bond.Harmonic`` applies
    (``F = −k (sep − r0)`` on each end along the bond; here expressed as the
    component on the nucleus bead toward the cytoskeleton bead, hence the sign
    flips to ``+k (sep − r0)``). Exposed as a scalar for the sign-sense
    Sanity-Gate check.

    Args:
        k_linc: Per-nesprin stiffness [N/m] (> 0).
        r0: Bond rest length [m] (> 0).
        sep: Current nucleus↔cytoskeleton separation [m] (≥ 0).

    Returns:
        Force component on the nucleus bead along ``Δr̂`` [N].

    Raises:
        ValueError: on a non-positive stiffness / rest length / negative sep.
    """
    _require_finite_positive("k_linc", k_linc)
    _require_finite_positive("r0", r0)
    if not (math.isfinite(sep) and sep >= 0.0):
        raise ValueError(f"sep must be finite and ≥ 0; got {sep!r}")
    return k_linc * (sep - r0)


def linc_bond_potential(k_linc: float, r0: float, sep: float) -> float:
    """Non-negative harmonic LINC bond potential ``½ k (sep − r0)²`` [J]."""
    _require_finite_positive("k_linc", k_linc)
    _require_finite_positive("r0", r0)
    if not (math.isfinite(sep) and sep >= 0.0):
        raise ValueError(f"sep must be finite and ≥ 0; got {sep!r}")
    d = sep - r0
    return 0.5 * k_linc * d * d


def linc_cfl_dt_max(
    k_linc: float, gamma_bead: float, cfl_safety_factor: float = 0.1
) -> float:
    """Max stable ``dt`` for the stiff LINC bond [s] (overdamped CFL gate).

    ``dt ≤ cfl_safety_factor · τ`` with ``τ = γ_bead / k_linc`` — the same
    convention erm.py / nucleus.py use for their stiff radial springs.

    Args:
        k_linc: Per-nesprin stiffness [N/m] (> 0).
        gamma_bead: Per-bead Stokes drag of the bonded beads [N·s/m] (> 0).
        cfl_safety_factor: Safety factor (default 0.1, D3 BAOAB convention).

    Returns:
        Maximum stable integration step ``dt`` [s].
    """
    _require_finite_positive("k_linc", k_linc)
    _require_finite_positive("gamma_bead", gamma_bead)
    _require_finite_positive("cfl_safety_factor", cfl_safety_factor)
    return cfl_safety_factor * (gamma_bead / k_linc)


# ---------------------------------------------------------------------------
# Geometric topology pairing (pure — returns bond pairs, mutates nothing)
# ---------------------------------------------------------------------------
def pair_linc_bridges(
    nucleus_positions: np.ndarray,
    cytoskeleton_positions: np.ndarray,
    *,
    capture_radius: float,
    n_bridges_max: int,
    unique_acceptor: bool = False,
) -> np.ndarray:
    """Geometrically pair nucleus beads to nearest cytoskeleton beads.

    PURE helper: computes which nucleus bead bonds to which cytoskeleton bead
    and RETURNS the index pairs ``(nuc_local_idx, cyto_local_idx)``. It mutates
    no simulation/snapshot/file. Each nucleus bead bonds to its single nearest
    cytoskeleton bead within ``capture_radius`` (one LINC per envelope bead);
    nucleus beads with no acceptor in range form no bond. The result is capped
    at ``n_bridges_max`` (keeping the SHORTEST bridges — the tightest, most
    load-bearing couplings — first).

    Args:
        nucleus_positions: ``(n_nuc, 3)`` nucleus-bead positions [m].
        cytoskeleton_positions: ``(n_cyto, 3)`` cytoskeleton-bead positions [m].
        capture_radius: Max pairing distance [m] (> 0).
        n_bridges_max: Cap on the number of bridges (≥ 0).
        unique_acceptor: If True, each cytoskeleton acceptor bead is used by AT
            MOST ONE bridge (greedy shortest-first matching) — "one nesprin per
            IF anchor point". This BOUNDS the per-acceptor bond degree to +1,
            which is REQUIRED on the integrated cell: without it many nucleus
            surface beads share a single nearest ``if_bead`` and concentrate
            LINC bonds on it, overflowing the HOOMD nlist per-particle exclusion
            cap (the same class as the single-hub MTOC degree limit). Default
            False preserves the one-LINC-per-nucleus-bead semantics for the
            standalone smoke / unit tests.

    Returns:
        ``(n_bridges, 2)`` int array of LOCAL index pairs
        ``[nuc_local_idx, cyto_local_idx]`` (empty ``(0, 2)`` if none).

    Raises:
        ValueError: on a non-positive capture radius or negative cap.
        RuntimeError: if scipy is unavailable.
    """
    _require_finite_positive("capture_radius", capture_radius)
    if n_bridges_max < 0:
        raise ValueError(f"n_bridges_max must be ≥ 0; got {n_bridges_max!r}")

    nuc = np.asarray(nucleus_positions, dtype=np.float64).reshape(-1, 3)
    cyto = np.asarray(cytoskeleton_positions, dtype=np.float64).reshape(-1, 3)
    if nuc.shape[0] == 0 or cyto.shape[0] == 0 or n_bridges_max == 0:
        return np.empty((0, 2), dtype=np.int64)
    if cKDTree is None:  # pragma: no cover
        raise RuntimeError("scipy.spatial.cKDTree is required for LINC pairing.")

    tree = cKDTree(cyto)
    # Nearest cytoskeleton acceptor for each nucleus bead.
    dist, idx = tree.query(nuc, k=1)
    dist = np.atleast_1d(dist)
    idx = np.atleast_1d(idx)

    in_range = dist <= capture_radius
    nuc_idx = np.nonzero(in_range)[0]
    if nuc_idx.size == 0:
        return np.empty((0, 2), dtype=np.int64)

    cyto_idx = idx[nuc_idx]
    d_sel = dist[nuc_idx]
    # Keep the shortest (tightest) bridges first, then cap.
    order = np.argsort(d_sel, kind="stable")
    nuc_ord = nuc_idx[order]
    cyto_ord = cyto_idx[order]
    if not unique_acceptor:
        nuc_keep = nuc_ord[:n_bridges_max]
        cyto_keep = cyto_ord[:n_bridges_max]
        return np.stack([nuc_keep, cyto_keep], axis=1).astype(np.int64)

    # Greedy 1:1 acceptor matching (shortest-first): each acceptor used once so
    # the per-acceptor bond degree gains at most +1 (nlist exclusion-cap safe).
    seen: set[int] = set()
    keep_n: list[int] = []
    keep_c: list[int] = []
    for ni, ci in zip(nuc_ord.tolist(), cyto_ord.tolist()):
        if ci in seen:
            continue
        seen.add(ci)
        keep_n.append(ni)
        keep_c.append(ci)
        if len(keep_n) >= n_bridges_max:
            break
    if not keep_n:
        return np.empty((0, 2), dtype=np.int64)
    return np.stack(
        [np.array(keep_n, dtype=np.int64), np.array(keep_c, dtype=np.int64)],
        axis=1,
    )


# ---------------------------------------------------------------------------
# Per-bond exact-r0 LINC layout (Option A: nucleus → perinuclear IF cage)
# ---------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class LINCLayout:
    """Resolved per-bond LINC bridge topology (EXACT-r0, force-free at build).

    Returned by :func:`compute_linc_layout` for the integrated cell. Each formed
    LINC bridge gets its OWN bond type carrying its EXACT as-built separation as
    the rest length, so the cage couples force-free at the resting geometry (the
    ~8 pN resting tension is generated by actomyosin AFTER equilibration, not
    injected as construction pre-stress — physiological-baseline rule).

    Attributes:
        pairs: ``(n_bridges, 2)`` uint32 GLOBAL particle tags
            ``[nucleus_tag, cytoskeleton_tag]`` (tag == row index on a
            tag-ordered build snapshot).
        r0: ``(n_bridges,)`` float64 EXACT as-built separation per bridge [m]
            (the bond's force-free rest length).
        type_names: per-bridge bond-type names (one type per bridge; every name
            keeps the ``linc_`` γ-denylist prefix).
    """

    pairs: np.ndarray
    r0: np.ndarray
    type_names: tuple[str, ...]

    @property
    def n_bridges(self) -> int:
        return int(self.pairs.shape[0])


def compute_linc_layout(
    snap: "object",
    p: ResolvedLINC,
    *,
    nucleus_tags: np.ndarray | list[int],
    cytoskeleton_tags: np.ndarray | list[int],
) -> "LINCLayout | None":
    """Geometric per-bond EXACT-r0 LINC layout (pure; mutates nothing).

    Pairs each nucleus bead to its nearest cytoskeletal acceptor (Option A: the
    perinuclear ``if_bead`` cage) within ``p.capture_radius`` and records the
    EXACT as-built separation as each bond's force-free rest length. Returns
    ``None`` (a no-op signal) when LINC is disabled, no tags are supplied, or
    the geometry forms zero bridges. Does NOT require ``k_linc`` (pairing is
    pure geometry; the stiffness is consumed only at potential-configure time).

    The snapshot is consumed TAG-ORDERED (global tag == row index, cell.py
    extender convention); ``.particles.tag`` is never read (a build-time
    ``gsd.hoomd.Frame`` does not expose it).

    Args:
        snap: A ``gsd.hoomd.Frame`` / ``hoomd.Snapshot`` with
            ``.particles.position`` / ``.N``.
        p: Resolved LINC parameters.
        nucleus_tags: ``nucleus_bead`` tags (== row indices).
        cytoskeleton_tags: acceptor (``if_bead`` / cortex / SF / MT) tags.

    Returns:
        A :class:`LINCLayout` or ``None`` (no bridges / disabled).
    """
    if not p.enabled:
        return None
    nuc_tags = np.asarray(nucleus_tags, dtype=np.int64).reshape(-1)
    cyto_tags = np.asarray(cytoskeleton_tags, dtype=np.int64).reshape(-1)
    if nuc_tags.size == 0 or cyto_tags.size == 0:
        return None

    pos = np.asarray(snap.particles.position, dtype=np.float64).reshape(-1, 3)
    n_part = int(snap.particles.N)
    nuc_rows = nuc_tags[(nuc_tags >= 0) & (nuc_tags < n_part)]
    cyto_rows = cyto_tags[(cyto_tags >= 0) & (cyto_tags < n_part)]
    if nuc_rows.size == 0 or cyto_rows.size == 0:
        return None

    pairs_local = pair_linc_bridges(
        pos[nuc_rows],
        pos[cyto_rows],
        capture_radius=p.capture_radius,
        n_bridges_max=p.n_bridges_max,
        # One nesprin per IF anchor: bound the per-acceptor bond degree to +1 so
        # the integrated cell does not overflow the nlist exclusion cap (many
        # nucleus surface beads otherwise share one nearest if_bead).
        unique_acceptor=True,
    )
    if pairs_local.shape[0] == 0:
        return None

    # Local-index pairs → global tags (rows ARE the global tags, tag-ordered).
    nuc_bond_tags = nuc_rows[pairs_local[:, 0]]
    cyto_bond_tags = cyto_rows[pairs_local[:, 1]]
    pairs = np.stack([nuc_bond_tags, cyto_bond_tags], axis=1).astype(np.uint32)
    sep = np.linalg.norm(
        pos[nuc_bond_tags] - pos[cyto_bond_tags], axis=1
    ).astype(np.float64)
    names = tuple(linc_nesprin_exact_type_names(int(pairs.shape[0])))
    return LINCLayout(pairs=pairs, r0=sep, type_names=names)


# ---------------------------------------------------------------------------
# Fresh-snapshot bond-append helper (shared; handles gsd.Frame + hoomd.Snapshot)
# ---------------------------------------------------------------------------
def _rebuild_snapshot_appending_linc_bonds(
    snap: "object",
    new_pairs: np.ndarray,
    per_bond_type_names: list[str],
) -> "hoomd.Snapshot":
    """Return a FRESH ``hoomd.Snapshot`` = ``snap`` + appended LINC bonds.

    LINC adds NO particles (bonds-only), so every particle / angle / dihedral /
    improper is copied UNCHANGED; only the supplied ``linc_`` bonds are appended.
    Handles both a build-time ``gsd.hoomd.Frame`` (None-valued unset fields) and
    a ``hoomd.Snapshot``, mirroring ``microtubules.extend_snapshot_with_micro
    tubules`` — an in-place mutation of ``snap.bonds`` is NOT used because a
    ``hoomd.Snapshot`` BondDataSnapshot has no whole-array setter.

    Args:
        new_pairs: ``(n_new, 2)`` GLOBAL bond tags.
        per_bond_type_names: ``len == n_new`` bond-type name per appended bond
            (all ``linc_``-prefixed). Distinct names are registered once.
    """
    def _grp(arr, width: int) -> np.ndarray:
        if arr is None:
            return np.empty((0, width), dtype=np.int64)
        a = np.asarray(arr, dtype=np.int64)
        return a.reshape(-1, width) if a.size else np.empty((0, width), np.int64)

    def _tid(arr) -> np.ndarray:
        if arr is None:
            return np.empty((0,), dtype=np.uint32)
        return np.asarray(arr, dtype=np.uint32).reshape(-1)

    def _pf(arr, default: np.ndarray) -> np.ndarray:
        return default if arr is None else np.asarray(arr)

    n = int(snap.particles.N)
    write = hoomd.Snapshot()
    write.particles.N = n
    write.particles.types = list(snap.particles.types)
    write.particles.typeid[:] = _pf(
        snap.particles.typeid, np.zeros(n, dtype=np.uint32)
    ).astype(np.uint32).reshape(-1)
    write.particles.position[:] = _pf(
        snap.particles.position, np.zeros((n, 3), dtype=np.float64)
    ).astype(np.float64).reshape(-1, 3)
    write.particles.velocity[:] = _pf(
        snap.particles.velocity, np.zeros((n, 3), dtype=np.float64)
    ).astype(np.float64).reshape(-1, 3)
    write.particles.mass[:] = _pf(
        snap.particles.mass, np.ones(n, dtype=np.float64)
    ).astype(np.float64).reshape(-1)
    write.particles.image[:] = _pf(
        snap.particles.image, np.zeros((n, 3), dtype=np.int32)
    ).astype(np.int32).reshape(-1, 3)
    # A minimal build-time gsd.hoomd.Frame may leave the box unset (None); a
    # hoomd.Snapshot always carries one. Default to a unit box only when absent
    # (the unit tests never create a state from the no-box frame).
    _box = snap.configuration.box
    write.configuration.box = (
        list(_box) if _box is not None else [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    )

    # Bonds: existing + appended LINC bonds. Register each distinct LINC type
    # once, then resolve the per-bond typeid.
    bond_types = list(snap.bonds.types or [])
    for tn in per_bond_type_names:
        if tn not in bond_types:
            bond_types.append(tn)
    old_bg = _grp(snap.bonds.group, 2)
    old_bt = _tid(snap.bonds.typeid)
    new_bg = np.asarray(new_pairs, dtype=np.int64).reshape(-1, 2)
    new_bt = np.array(
        [bond_types.index(tn) for tn in per_bond_type_names], dtype=np.uint32
    )
    merged_bg = np.concatenate([old_bg, new_bg], axis=0).astype(np.uint32)
    merged_bt = np.concatenate([old_bt, new_bt]).astype(np.uint32)
    write.bonds.types = bond_types
    write.bonds.N = int(merged_bg.shape[0])
    if merged_bg.shape[0] > 0:
        write.bonds.group[:] = merged_bg
        write.bonds.typeid[:] = merged_bt

    # Pass through angles / dihedrals / impropers untouched (None-guarded).
    for grp_name in ("angles", "dihedrals", "impropers"):
        src = getattr(snap, grp_name, None)
        if src is None:
            continue
        src_n = int(getattr(src, "N", 0) or 0)
        src_types = list(getattr(src, "types", None) or [])
        dst = getattr(write, grp_name)
        if src_n > 0:
            dst.N = src_n
            dst.types = src_types
            dst.group[:] = np.asarray(src.group)
            dst.typeid[:] = np.asarray(src.typeid)
        elif src_types:
            dst.N = 0
            dst.types = src_types

    return write


# ---------------------------------------------------------------------------
# Snapshot-extension builders (add linc_nesprin bonds; no-op when disabled)
# ---------------------------------------------------------------------------
def extend_snapshot_with_linc(
    snap: "object",
    p: ResolvedLINC,
    *,
    nucleus_tags: np.ndarray | list[int],
    cytoskeleton_tags: np.ndarray | list[int],
) -> "object":
    """Append SINGLE-r0 ``linc_nesprin`` harmonic bonds (nucleus↔cytoskeleton).

    The simple uniform-``r0`` path: every formed bridge gets the one canonical
    ``linc_nesprin`` bond type at the resolved ``p.r0`` rest length. It is the
    right path when the acceptors are SEEDED at ``r0`` (e.g. the perinuclear-cap
    smoke topology) so every bond is born force-free. For the integrated cell's
    Option A (nucleus → the pre-existing IF cage), where the as-built
    separations vary (CV ~24%) and a single ``r0`` cannot be force-free, use
    :func:`compute_linc_layout` + :func:`extend_snapshot_with_linc_layout`
    (per-bond EXACT-r0) instead.

    Returns a FRESH ``hoomd.Snapshot`` (LINC adds bonds only; all particles /
    angles are copied unchanged). Works in isolation on a ``gsd.hoomd.Frame`` or
    a ``hoomd.Snapshot``, consumed TAG-ORDERED (global tag == row index;
    ``.particles.tag`` is never read). STRICT NO-OP — returns the input snapshot
    UNCHANGED (same object) — when:

    * ``p.enabled`` is False (DEFAULT-OFF), or
    * no nucleus / cytoskeleton tags are supplied, or
    * the geometric pairing forms zero bridges.

    Args:
        snap: A ``gsd.hoomd.Frame`` / ``hoomd.Snapshot`` (``.particles.position``
            / ``.N`` / ``.bonds.*``).
        p: Resolved LINC parameters. ``.enabled = False`` ⇒ no-op.
        nucleus_tags: ``nucleus_bead`` tags (== row indices).
        cytoskeleton_tags: acceptor bead tags (== row indices).

    Returns:
        A fresh extended ``hoomd.Snapshot`` (bonds added) or the input ``snap``
        unchanged (no-op).

    Raises:
        NotImplementedError: if ``p.enabled`` and ``p.k_linc is None`` — the
            load-bearing stiffness is UNKNOWN (see module ``PI_DECISIONS``).
    """
    # DEFAULT-OFF strict no-op.
    if not p.enabled:
        return snap
    # Honest disabled-build: stiffness genuinely unknown (raise BEFORE pairing,
    # matching the historical contract — k_linc None always raises on an enabled
    # build attempt).
    if p.k_linc is None:
        raise NotImplementedError(
            "LINC topology build requires k_linc (per-nesprin spectrin-repeat "
            "spring stiffness, N/m), which is UNKNOWN — see "
            "ffn_sim.cell.linc.PI_DECISIONS. Nesprin-2 giant has no single "
            "Hookean stiffness in the literature (nonlinear repeat unfolding); "
            "the ~8 pN resting tension (Déjardin 2020 JCB; Arsenovic 2016 "
            "reported relative FRET, not an absolute pN value) is a force, not a "
            "stiffness. Surface to PI for a literature-anchored k_linc before "
            "enabling LINC. (Resolving the config is allowed; building the "
            "bonds is not, until k_linc is set.)"
        )

    layout = compute_linc_layout(
        snap, p, nucleus_tags=nucleus_tags, cytoskeleton_tags=cytoskeleton_tags
    )
    if layout is None:
        return snap  # no tags / no acceptor in range → no-op

    # SINGLE-r0: every bridge rides the one canonical linc_nesprin type.
    names = [p.bond_type_name] * layout.n_bridges
    return _rebuild_snapshot_appending_linc_bonds(snap, layout.pairs, names)


def extend_snapshot_with_linc_layout(
    snap: "object",
    p: ResolvedLINC,
    layout: "LINCLayout | None",
) -> "object":
    """Append per-bond EXACT-r0 ``linc_nesprin`` bonds from a precomputed layout.

    The integrated-cell Option-A path: each LINC bridge gets its own bond type
    (``linc_nesprin`` / ``linc_nesprin_b{i}``) so it is born force-free at its
    EXACT as-built separation (the :func:`compute_linc_layout` ``r0`` per
    bridge), configured by :func:`configure_linc_bond_potential`. Returns a
    FRESH ``hoomd.Snapshot``; strict no-op (returns ``snap``) when LINC is
    disabled or ``layout`` is ``None``.

    Raises:
        NotImplementedError: if ``p.enabled`` and ``p.k_linc is None`` (the
            stiffness is UNKNOWN — see ``PI_DECISIONS``). Raised here too so the
            build halts BEFORE appending bonds it cannot give a potential.
    """
    if not p.enabled or layout is None:
        return snap
    if p.k_linc is None:
        raise NotImplementedError(
            "LINC topology build requires k_linc (UNKNOWN) — see "
            "ffn_sim.cell.linc.PI_DECISIONS. Surface to PI before enabling LINC."
        )
    return _rebuild_snapshot_appending_linc_bonds(
        snap, layout.pairs, list(layout.type_names)
    )


def configure_linc_bond_potential(
    harmonic: "object",
    p: ResolvedLINC,
    *,
    layout: "LINCLayout | None" = None,
) -> "object":
    """Set the ``linc_nesprin`` parameters on the cell's shared ``md.bond.Harmonic``.

    The Lead owns the single shared ``md.bond.Harmonic``; this helper writes the
    LINC bond-type params onto it (it does NOT construct a competing potential).
    No-op when LINC is disabled.

    * ``layout=None`` (SINGLE-r0): registers the one canonical ``linc_nesprin``
      type at ``k=k_linc`` / ``r0=p.r0``.
    * ``layout`` given (per-bond EXACT-r0): registers EACH ``layout.type_names``
      at ``k=k_linc`` / ``r0=layout.r0[i]`` — every bond force-free at its
      as-built separation (matches :func:`extend_snapshot_with_linc_layout`).

    Args:
        harmonic: The cell's shared ``hoomd.md.bond.Harmonic``.
        p: Resolved LINC parameters.
        layout: Optional per-bond EXACT-r0 layout from
            :func:`compute_linc_layout`.

    Returns:
        The same ``harmonic`` (for chaining / introspection).

    Raises:
        NotImplementedError: if enabled and ``k_linc is None`` (unknown
            stiffness — see ``PI_DECISIONS``).
    """
    if not p.enabled:
        return harmonic
    if p.k_linc is None:
        raise NotImplementedError(
            "Cannot configure the linc_nesprin bond potential: k_linc is "
            "UNKNOWN (see ffn_sim.cell.linc.PI_DECISIONS). Surface to PI."
        )
    k = float(p.k_linc)
    if layout is None:
        harmonic.params[p.bond_type_name] = dict(k=k, r0=float(p.r0))
        return harmonic
    for name, r0 in zip(layout.type_names, np.asarray(layout.r0, dtype=float)):
        harmonic.params[name] = dict(k=k, r0=float(r0))
    return harmonic
