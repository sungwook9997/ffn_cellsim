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
ONE; Arsenovic 2016 Biophys J measured ~2 pN resting tension per nesprin by
FRET). The topology is built geometrically: each nucleus bead bonds to its
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
resting tension ``f_rest`` (~2 pN, Arsenovic 2016) IS known and is exposed as a
provenance field / diagnostic, but a tension is not a stiffness and is not used
to invent one.

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
   - The resting per-nesprin tension ``f_rest`` (~2 pN, Arsenovic 2016) is the
     validation oracle: at equilibrium the mean LINC bond tension should sit
     near ``f_rest`` once ``k_linc`` and the pre-strain are set by the host. We
     expose ``f_rest`` and the per-bond tension formula but do NOT claim any
     gate PASSes here (that is a Lead integration / production concern).

Compartment Performance Contract
--------------------------------
* **Particle types added**: NONE. LINC reuses existing ``nucleus_bead`` and
  cytoskeleton particles; it only adds bonds. (No new particle cloud — this is
  the key difference from nucleus.py.)
* **Particle count** at n_fil=1000 / native ~38000: +0 particles at both
  scales. LINC is topology-only.
* **Bond/angle count**: +``n_bridges`` bonds (one ``linc_nesprin`` bond per
  formed bridge), ``n_bridges = min(n_nucleus_beads, n_bridges_max)`` capped.
  No angles. Biology: thousands of LINC per nucleus (Lombardi 2011); at the
  ×40 mesoscale ``n_bridges`` is O(n_nucleus_beads), i.e. tens–hundreds, with
  ``n_bridges_max`` capping it. Zero per-bond runtime cost beyond the builtin
  harmonic evaluation.
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
- Arsenovic P.T. et al. (2016) *Biophys. J.* 110:34 — FRET tension sensor:
  nesprin-2G bears ~2 pN of resting tension per molecule in living cells.
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

#: Open PI decisions (CLAUDE.md no-magic-number protocol). Empty when none.
PI_DECISIONS: list[str] = [
    "k_linc (single-nesprin spectrin-repeat spring stiffness, N/m): UNKNOWN. "
    "Nesprin-2 giant is a ~56-repeat spectrin spring whose force-extension is "
    "nonlinear (sequential repeat unfolding ~25-35 pN plateaus; Autore 2013), "
    "so there is no single Hookean N/m in the literature. Resting tension is "
    "known (~2 pN, Arsenovic 2016) but a tension is not a stiffness and is not "
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
        n_bridges_max: Hard cap on the number of LINC bonds formed (one per
            nucleus bead at most; biology has thousands per nucleus —
            Lombardi 2011 — capped here at the ×40 mesoscale).
        f_rest: Resting per-nesprin tension [N] (~2 pN, Arsenovic 2016).
            PROVENANCE / diagnostic only — NOT used to derive ``k_linc``.
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

#: Resting per-nesprin tension measured by FRET in living cells
#: (Arsenovic 2016 Biophys J — ~2 pN on nesprin-2G). Provenance only.
_NESPRIN_REST_TENSION_N: float = 2.0e-12     # N   (2 pN, Arsenovic 2016)


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
    (~50 nm, Crisp 2006); the capture radius defaults to a small multiple of
    that gap (the nucleus↔cortex separation at the envelope is short, since
    LINC physically spans the membrane). ``k_linc`` is UNKNOWN by default
    (``None``) — resolving is allowed, but BUILDING the load-bearing bonds with
    ``k_linc is None`` raises ``NotImplementedError`` (no invented number).

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
    # Capture radius default: 3× the perinuclear gap (the nucleus-bead ↔
    # cortex/SF-bead separation at the envelope is short; LINC spans the
    # membrane). A GEOMETRIC choice, grid-invariant (multiple of a length).
    capture_radius = float(cfg.get("capture_radius", 3.0 * _PERINUCLEAR_GAP_M))
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
    nuc_idx = nuc_idx[order][:n_bridges_max]
    cyto_idx = cyto_idx[order][:n_bridges_max]
    return np.stack([nuc_idx, cyto_idx], axis=1).astype(np.int64)


# ---------------------------------------------------------------------------
# Snapshot-extension builder (adds linc_nesprin bonds; no-op when disabled)
# ---------------------------------------------------------------------------
def extend_snapshot_with_linc(
    snap: "object",
    p: ResolvedLINC,
    *,
    nucleus_tags: np.ndarray | list[int],
    cytoskeleton_tags: np.ndarray | list[int],
) -> "object":
    """Add ``linc_nesprin`` harmonic bonds between nucleus and cytoskeleton.

    Mirrors the ``cell/cell.py`` bond-type-extension pattern: carries existing
    bonds across, registers the ``linc_nesprin`` bond type (if absent), pairs
    nucleus↔cytoskeleton beads geometrically, and appends one bond per formed
    bridge. The function is callable IN ISOLATION on a ``gsd.hoomd.Frame`` (or
    any snapshot exposing ``.particles.position/.tag`` and ``.bonds.*``) and is
    a STRICT NO-OP — returning the input snapshot unchanged — when:

    * ``p.enabled`` is False (DEFAULT-OFF), or
    * no nucleus / cytoskeleton tags are supplied, or
    * the geometric pairing forms zero bridges.

    It mutates the passed snapshot's ``bonds`` arrays in place (matching the
    single-writer gsd-Frame convention in cell.py) and returns it.

    Args:
        snap: A ``gsd.hoomd.Frame`` (or HOOMD snapshot) with positions, tags,
            and ``bonds`` arrays.
        p: Resolved LINC parameters. ``.enabled = False`` ⇒ no-op.
        nucleus_tags: Particle tags of the ``nucleus_bead`` cloud.
        cytoskeleton_tags: Particle tags of the cytoskeletal acceptor beads
            (cortex actin / stress fiber / MT) LINC may bond to.

    Returns:
        The (possibly extended) snapshot.

    Raises:
        NotImplementedError: if ``p.enabled`` and ``p.k_linc is None`` — the
            load-bearing stiffness is UNKNOWN (see module ``PI_DECISIONS``).
            Honesty over completeness: no invented number.
    """
    # DEFAULT-OFF strict no-op.
    if not p.enabled:
        return snap

    # Honest disabled-build: stiffness genuinely unknown.
    if p.k_linc is None:
        raise NotImplementedError(
            "LINC topology build requires k_linc (per-nesprin spectrin-repeat "
            "spring stiffness, N/m), which is UNKNOWN — see "
            "ffn_sim.cell.linc.PI_DECISIONS. Nesprin-2 giant has no single "
            "Hookean stiffness in the literature (nonlinear repeat unfolding); "
            "the ~2 pN resting tension (Arsenovic 2016) is a force, not a "
            "stiffness. Surface to PI for a literature-anchored k_linc before "
            "enabling LINC. (Resolving the config is allowed; building the "
            "bonds is not, until k_linc is set.)"
        )

    nuc_tags = np.asarray(nucleus_tags, dtype=np.int64).reshape(-1)
    cyto_tags = np.asarray(cytoskeleton_tags, dtype=np.int64).reshape(-1)
    if nuc_tags.size == 0 or cyto_tags.size == 0:
        return snap  # nothing to bond → no-op

    pos = np.asarray(snap.particles.position, dtype=np.float64).reshape(-1, 3)
    tags = np.asarray(snap.particles.tag, dtype=np.int64).reshape(-1)
    # Map global tag → row index.
    tag_to_row = {int(t): i for i, t in enumerate(tags)}

    nuc_rows = np.array(
        [tag_to_row[int(t)] for t in nuc_tags if int(t) in tag_to_row],
        dtype=np.int64,
    )
    cyto_rows = np.array(
        [tag_to_row[int(t)] for t in cyto_tags if int(t) in tag_to_row],
        dtype=np.int64,
    )
    if nuc_rows.size == 0 or cyto_rows.size == 0:
        return snap

    pairs_local = pair_linc_bridges(
        pos[nuc_rows],
        pos[cyto_rows],
        capture_radius=p.capture_radius,
        n_bridges_max=p.n_bridges_max,
    )
    if pairs_local.shape[0] == 0:
        return snap  # no acceptor in range → no-op

    # Local-index pairs → global tags (HOOMD bond groups are tags).
    nuc_bond_tags = tags[nuc_rows[pairs_local[:, 0]]]
    cyto_bond_tags = tags[cyto_rows[pairs_local[:, 1]]]
    new_pairs = np.stack([nuc_bond_tags, cyto_bond_tags], axis=1).astype(
        np.uint32
    )
    n_new = int(new_pairs.shape[0])

    # Register the bond type (carry existing bond types/groups/typeids across).
    old_types = list(snap.bonds.types) if snap.bonds.types is not None else []
    old_N = int(snap.bonds.N)
    old_group = (
        np.asarray(snap.bonds.group, dtype=np.int64).reshape(old_N, 2)
        if old_N > 0
        else np.empty((0, 2), dtype=np.int64)
    ).astype(np.uint32)
    old_typeid = (
        np.asarray(snap.bonds.typeid, dtype=np.uint32)
        if old_N > 0
        else np.empty((0,), dtype=np.uint32)
    )

    new_types = list(old_types)
    if p.bond_type_name not in new_types:
        new_types.append(p.bond_type_name)
    linc_typeid = new_types.index(p.bond_type_name)

    group_all = np.concatenate([old_group, new_pairs], axis=0)
    typeid_all = np.concatenate(
        [old_typeid, np.full(n_new, linc_typeid, dtype=np.uint32)]
    )

    snap.bonds.N = int(group_all.shape[0])
    snap.bonds.types = new_types
    snap.bonds.group = group_all.astype(np.uint32)
    snap.bonds.typeid = typeid_all.astype(np.uint32)
    return snap


def configure_linc_bond_potential(
    harmonic: "object",
    p: ResolvedLINC,
) -> "object":
    """Set the ``linc_nesprin`` parameters on a builtin ``md.bond.Harmonic``.

    The Lead owns the single ``md.bond.Harmonic`` instance for the cell; this
    helper writes the LINC bond-type params into it (``k`` = ``k_linc``,
    ``r0`` = ``p.r0``) without constructing a competing potential. It is a
    no-op when LINC is disabled.

    Args:
        harmonic: A ``hoomd.md.bond.Harmonic`` whose ``.params`` mapping accepts
            the LINC bond type.
        p: Resolved LINC parameters.

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
    harmonic.params[p.bond_type_name] = dict(k=float(p.k_linc), r0=float(p.r0))
    return harmonic
