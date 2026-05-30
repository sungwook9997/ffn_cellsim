"""H.9 nucleus mechanics — stiff deformable nuclear element (KU-3.B2).

ADDITIVE, DEFAULT-OFF module (EXTEND memo Part B Template 2 + Template 1;
``docs/briefs/H9_nucleus.md``). The cell otherwise has NO explicit nucleus,
so the composite KU-3.1 rounding / confined-migration mechanics are missing
their mechanically co-dominant, rate-limiting element. This module supplies
a stiff nuclear element: a ``nucleus_bead`` particle cloud near the cell
centroid, confined by a two-regime radial harmonic (chromatin interior +
lamin-shell strain-stiffening). NO integrator change; when the Lead does NOT
attach :class:`NucleusConfinement`, pre-H.9 runs are bit-for-bit identical.

The brief (``docs/briefs/H9_nucleus.md``) is a DRAFT skeleton, NOT a
contract; this is the correct MINIMAL first version (radial confinement +
strain-stiffening knee + the E_nuc→k_eff dimensional bridge + a PURE
bead-builder helper). Deeper physics (a lamin SHELL of explicit bonded beads,
chromatin internal crosslinks, nuclear-envelope rupture under pore transit,
the MMP pore-enlargement rescue branch) is documented TODO below.

Mechanism (fine-grained per CLAUDE.md — minimal but real)
---------------------------------------------------------
The nucleus is represented by ``n_beads`` ``nucleus_bead`` particles seeded
in a roughly-spherical cloud of radius ``R_nuc`` about the cell centroid
``r_c``. A radial harmonic confinement holds each bead near the nuclear
radius ``R_nuc`` measured from the (live) bead centroid::

    d   = r − R_nuc        (signed radial displacement from the shell)
    F_i = −k_eff(d) · d · r̂_i          r̂_i = (r_i − r_c)/|r_i − r_c|

This mirrors the ``cortex/erm.py`` radial-spring SHAPE (a harmonic toward a
target radius) but, like ``cortex/enclosed_volume.py``, anchors the normals
at the BEAD CENTROID (not a fixed origin) so the confinement injects ≈ 0 net
momentum (Sanity Gate §3) — the nucleus is an INTERNAL element, free to
translate with the cell, not pinned to the lab frame.

Two-regime stiffness (chromatin + lamin shell)
----------------------------------------------
KU-3.B2 sets the nucleus as a soft chromatin interior wrapped in a
strain-stiffening lamin-A/C shell. We model this as a BILINEAR radial
spring in the signed displacement ``d`` (continuous force at the knee):

* small strain ``|d| ≤ d_knee``: stiffness ``k_chrom`` (chromatin).
* large strain ``|d| > d_knee``: stiffness ``k_chrom + k_lamin`` (the lamin
  shell engages and the response stiffens).

The force magnitude is built so it is CONTINUOUS at the knee
``d = ±d_knee`` (standard bilinear / strain-stiffening spring — value
matches across the knee by construction; the SLOPE steps up). Defining the
knee force ``F_knee = k_chrom · d_knee``::

    |d| ≤ d_knee :  F_mag(d) = −k_chrom · d
    |d| > d_knee :  F_mag(d) = −sign(d)·[ F_knee + (k_chrom+k_lamin)·(|d|−d_knee) ]

At ``|d| = d_knee`` both branches give ``F_mag = −sign(d)·F_knee`` → C0
continuous (Sanity Gate §strain-stiffening). The bilinear potential is
C1 (continuous force) but NOT C2 (slope jumps) — physically correct for an
engagement knee. This is the minimal mechanistic two-regime law; a smooth
(``tanh``-blended) knee is a documented TODO.

E_nuc → k_eff dimensional bridge (no invented number)
-----------------------------------------------------
The chromatin stiffness ``k_chrom`` is DERIVED from the continuum Young's
modulus ``E_nuc`` via the bead surface-area share — the SAME bridge
``cortex/enclosed_volume.py`` uses for ``F_i = ΔP·A_i`` (there ΔP = bulk
stress; here σ = E·ε). A radial bead displacement ``d`` is a strain
``ε = d / R_nuc``; the elastic stress is ``σ = E_nuc · ε``; the restoring
force on a bead carrying its share ``A_bead = 4π R_nuc² / n_beads`` of the
nuclear surface is::

    F = σ · A_bead = E_nuc · (d / R_nuc) · (4π R_nuc² / n_beads)
      = (4π E_nuc R_nuc / n_beads) · d

so the per-bead chromatin stiffness is the **dimensional bridge**::

    k_chrom = 4π · E_nuc · R_nuc / n_beads          [Pa·m = N/m]      (B)

Dimensions: [Pa]·[m] = [N/m]. ✓  Grid-aware: the ``1/n_beads`` share makes
the TOTAL radial confinement stiffness ``n_beads · k_chrom = 4π E_nuc R_nuc``
independent of the discretization (intensive material response), exactly as
enclosed_volume's per-bead ``A_i = S/N`` share keeps the net pressure load
N-independent. It is NOT chosen to pass a gate: the strain-stiffening /
continuity / sign tests do not use ``E_nuc`` at all (they check the SHAPE of
the law for any positive stiffness).

The lamin-shell increment ``k_lamin`` is set by the nucleus:cytoplasm
stiffness RATIO anchor (KU-3.B2.1): the high-strain stiffness is
``k_chrom + k_lamin``; we anchor the shell so the high-strain modulus is
``ratio_lamin ×`` the chromatin modulus, i.e.
``k_lamin = (ratio_lamin − 1) · k_chrom``. ``ratio_lamin`` is taken from the
KU-3.B2.1 in-situ nucleus:cytoplasm band 1.4–5× — NOT the 10× isolated value
(explicitly NOT hard-coded per the brief). Default ``ratio_lamin = 3.0``
(mid in-situ decade). This is a stated literature ratio, range-checked
against 1.4–5×; a value outside the band raises (guards a tuned override).

Map to confined migration (KU-1.V.3.3 pore limit) — coupling note
-----------------------------------------------------------------
KU-3.B2.3 / KU-1.V.3.3: confined-migration ARREST occurs when the nucleus
is compressed to ~10% cross-section / a critical pore ~7 µm² (tumor). This
module exposes :meth:`NucleusConfinement.critical_pore_area` /
:func:`nucleus_cross_section_area` so the Lead's confined-migration driver
can detect the arrest condition (nuclear cross-section ≤ critical pore). The
ARREST physics itself (an ECM pore as a geometric constraint + the lamin
shell resisting transit) is a Lead integration concern wiring this nucleus
to the ECM pore — documented TODO; here we only NOTE the coupling and supply
the geometric predicate. We do NOT claim the KU-1.V.3.3 / KU-3.1 gate PASSes.

Magic-Number Block — every physical constant
--------------------------------------------
* ``E_nuc`` (Young's modulus) — KU-3.B2.1: nuclear modulus 1–10 kPa
  (in-situ ~5, isolated ~8). Default 5.0e3 Pa (in-situ anchor). Resolver
  range-checks 1e3–1e4 Pa; outside raises.
* ``ratio_lamin`` (nucleus:cytoplasm high-strain stiffness ratio) —
  KU-3.B2.1: 1.4–5× in-situ (NOT 10× isolated — brief explicit). Default
  3.0. Resolver range-checks [1.4, 5.0].
* ``knee_strain`` (lamin engagement strain) — KU-3.B2: the lamin shell
  strain-stiffens above a threshold; the chromatin small-strain regime
  spans ``|ε| ≤ knee_strain``. No tight single literature value exists for
  the engagement strain at this ×40 coarse scale, so it is a REQUIRED-with-
  documented-default parameter (default 0.10 = 10% radial strain, the same
  order as the ~10% cross-section arrest scale KU-3.B2.3); flagged below.
* ``R_nuc`` (nuclear radius) — REQUIRED (no default): set by the host cell
  geometry (a typical mammalian nucleus is ~5 µm diameter, R_nuc ~ 2.5 µm,
  but the Lead supplies the value consistent with the cell's R_cell so it is
  not invented here). ``critical_pore_area_um2`` KU-1.V.3.3 default 7.0 µm²
  (tumor; T-cell 4, neutrophil 2 — the Lead selects per cell type).
* ``gamma_nuc`` (per-bead Stokes drag) — REQUIRED for the CFL gate (supplied
  by the host, = 6π η R_bead like the cortex ``gamma_b``); no default.

Lamin-A scaling (KU-3.B2.2) — exposed as oracle helpers, not runtime
--------------------------------------------------------------------
KU-3.B2.2: nuclear viscosity ~ [lamin-A]^(3±1), elasticity ~ [lamin-A]^(~0.5)
(lamin-A → viscosity, lamin-B → elasticity). These are exposed as TEST-ONLY
scaling oracles :func:`lamin_a_viscosity_scaling` /
:func:`lamin_a_elasticity_scaling` (a relative-modulus multiplier the Lead
can fold into ``E_nuc`` when modelling lamin-A up/down-regulation, e.g.
cancer lamin-A/C loss → softer nucleus). They are NOT called by the runtime
force compute (no per-step lamin dynamics in this minimal version — TODO).

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_nucleus.py``.*

1. **Dimensional analysis**
   - ``k_chrom = 4π E_nuc R_nuc / n_beads`` [Pa·m = N/m]; ``d`` [m] →
     ``F = k·d`` [N]. ✓  ``k_lamin = (ratio_lamin−1) k_chrom`` [N/m]. ✓
   - ``U = ½ k_chrom d²`` (interior) + lamin shell increment [J]. ✓
   - At ``r = R_nuc`` (``d = 0``): ``F = 0``, ``U = 0`` — force-free at the
     nuclear shell (the zero-force configuration).
2. **Boundary cases**
   - ``r = R_nuc`` → ``d = 0`` → ``F = 0`` (zero-force config).
   - bead AT the centroid (``r = 0``): ``r̂`` undefined → that bead's force
     set to 0 (it has no radial direction). Nucleus beads seed at
     ``r ≈ R_nuc ≫ 0`` so never hit in practice.
   - ``E_nuc``/``R_nuc``/``ratio_lamin``/``knee_strain``/``n_beads`` ≤ 0 or
     out of band: caught by ``resolve_nucleus`` / builder ValueError.
   - empty nucleus (n_beads = 0): the custom force iterates zero nucleus
     particles → trivially zero (the builder forbids n_beads < 1).
3. **Conservation invariants**
   - **Net force ≈ 0**: normals about the bead CENTROID → for a symmetric
     cloud ``Σ r̂_i ≈ 0`` so ``Σ F_i ≈ 0`` to the discretisation residual
     (the nucleus is INTERNAL; must not inject net momentum — unlike the
     lab-frame-anchored ERM field). STATIC test asserts ``|Σ F| ≪ Σ|F|``.
   - confinement acts ONLY on tags in ``nucleus_tag_range``; all other
     particles untouched (mask). STATIC.
4. **Numerical sanity (CFL)**
   - positions read float64; forces float64.
   - per-bead radial stiffness is at most the high-strain
     ``k_hi = k_chrom + k_lamin``; relax time ``τ = γ_nuc / k_hi`` must
     satisfy ``dt ≤ cfl_safety_factor · τ``; the attach helper raises (like
     erm.py / enclosed_volume.py) if violated, using ``k_hi`` (the stiffest
     branch) so the gate is conservative.
5. **Sign / sense**
   - ``r > R_nuc`` (bead outside the shell): ``d > 0`` → ``F ∝ −r̂`` pulls
     INWARD (confines toward the shell). STATIC.
   - ``r < R_nuc`` (bead inside): ``d < 0`` → ``F ∝ +r̂`` pushes OUTWARD
     (toward the shell). STATIC.
   - strain-stiffening MONOTONE and CONTINUOUS at the knee. STATIC.
6. **Measurement protocol**
   - E_nuc→k_eff bridge stated (eq. B) and unit-checked.
   - couples to KU-1.V.3.3 pore limit via
     :func:`nucleus_cross_section_area` / ``critical_pore_area`` (NOTED, not
     a claimed PASS).
   - BAOAB smoke: nucleus beads + confinement-on runs a few hundred steps
     with no NaN/Inf and the cloud stays bounded near R_nuc.

References
----------
- Brief: ``ffn_sim/docs/briefs/H9_nucleus.md`` (DRAFT skeleton).
- KU-3.B2.1 nuclear modulus E_nuc 1–10 kPa; nucleus:cytoplasm 1.4–5× in-situ
  / 3–10× isolated (technique-dependent; Notion KU v2 KU-3.B2).
- KU-3.B2.2 lamin-A scaling: viscosity ~ [lamin-A]^(3±1), elasticity ~
  [lamin-A]^(~0.5) (lamin-A→viscosity, lamin-B→elasticity).
- KU-3.B2.3 / KU-1.V.3.3 confined-migration arrest: ~10% cross-section /
  critical pore ~7 µm² (tumor) / 4 (T-cell) / 2 (neutrophil).
- ``ffn_sim/cortex/erm.py`` — radial-spring SHAPE analog (md.force.Custom
  radial field on a tag range + §1–6 Sanity Gate + CFL attach gate).
- ``ffn_sim/cortex/enclosed_volume.py`` — centroid-relative normals (net
  force ≈ 0 for an internal element) + per-bead area-share bridge.
- HOOMD 7 ``md.force.Custom`` API.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import hoomd
import hoomd.md as md


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedNucleus:
    """Nucleus confinement parameters (all SI).

    Stiffnesses are DERIVED from the continuum anchors via the E_nuc→k_eff
    bridge (eq. B in the module docstring); the raw anchors are retained as
    provenance.
    """

    # Derived per-bead radial stiffnesses [N/m]
    k_chrom: float               # chromatin interior small-strain stiffness
    k_lamin: float               # lamin-shell strain-stiffening increment
    R_nuc: float                 # m   nuclear radius (target shell radius)
    n_beads: int                 # nucleus bead count (sets the area share)
    d_knee: float                # m   knee displacement = knee_strain · R_nuc

    # Provenance anchors (Magic-Number Block)
    E_nuc: float = 5.0e3         # Pa   KU-3.B2.1 in-situ nuclear modulus
    ratio_lamin: float = 3.0     # —    KU-3.B2.1 nucleus:cyto in-situ (1.4–5×)
    knee_strain: float = 0.10    # —    KU-3.B2 lamin engagement strain
    critical_pore_area_um2: float = 7.0   # µm² KU-1.V.3.3 tumor pore limit

    @property
    def k_hi(self) -> float:
        """High-strain (lamin-engaged) per-bead stiffness [N/m] (CFL gate)."""
        return self.k_chrom + self.k_lamin

    @property
    def F_knee(self) -> float:
        """Force magnitude at the knee, ``k_chrom · d_knee`` [N]."""
        return self.k_chrom * self.d_knee


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


# KU-3.B2 / KU-1.V.3.3 literature bands (range-check guards).
_E_NUC_BAND_PA = (1.0e3, 1.0e4)        # KU-3.B2.1 1–10 kPa
_RATIO_LAMIN_BAND = (1.4, 5.0)         # KU-3.B2.1 in-situ (NOT 10× isolated)


def resolve_nucleus(
    cfg: dict,
    *,
    R_nuc: float,
    n_beads: int,
) -> ResolvedNucleus:
    """Resolve the ``nucleus`` config block into SI per-bead stiffnesses.

    ``cfg`` is the YAML root, the ``cell:`` sub-dict, or the ``cell.nucleus``
    / ``nucleus`` sub-dict. ``R_nuc`` (nuclear radius) and ``n_beads`` come
    from the host cell geometry (REQUIRED — not invented here).

    The chromatin stiffness is DERIVED from ``E_nuc`` via the area-share
    bridge ``k_chrom = 4π E_nuc R_nuc / n_beads`` (eq. B); the lamin
    increment from the nucleus:cytoplasm ratio
    ``k_lamin = (ratio_lamin − 1) · k_chrom``. Anchors are range-checked
    against the KU-3.B2.1 bands (a value outside the band raises — guards a
    silent tuned override).

    Args:
        cfg: Config mapping (root, ``cell``, or ``nucleus`` sub-dict).
        R_nuc: Nuclear radius [m] (host cell geometry; REQUIRED).
        n_beads: Number of nucleus beads (sets the per-bead area share).

    Returns:
        Resolved nucleus parameters (SI).

    Raises:
        ValueError: on a non-positive / out-of-band anchor (Sanity Gate §2).
    """
    if "cell" in cfg:
        cfg = cfg["cell"]
    if "nucleus" in cfg:
        cfg = cfg["nucleus"]

    E_nuc = float(cfg.get("E_nuc", 5.0e3))
    ratio_lamin = float(cfg.get("ratio_lamin", 3.0))
    knee_strain = float(cfg.get("knee_strain", 0.10))
    critical_pore = float(cfg.get("critical_pore_area_um2", 7.0))

    # §2 boundary checks (REQUIRED host geometry + anchors).
    _require_finite_positive("R_nuc", R_nuc)
    if int(n_beads) < 1:
        raise ValueError(f"n_beads must be ≥ 1; got {n_beads!r}")
    _require_finite_positive("E_nuc", E_nuc)
    _require_finite_positive("ratio_lamin", ratio_lamin)
    _require_finite_positive("critical_pore_area_um2", critical_pore)
    if not (math.isfinite(knee_strain) and 0.0 < knee_strain < 1.0):
        raise ValueError(
            f"knee_strain must be finite and in (0, 1); got {knee_strain}"
        )

    # KU-3.B2.1 band guards (no magic number / no tuned override).
    lo, hi = _E_NUC_BAND_PA
    if not (lo <= E_nuc <= hi):
        raise ValueError(
            f"E_nuc = {E_nuc:.3e} Pa outside the KU-3.B2.1 nuclear-modulus "
            f"band [{lo:.1e}, {hi:.1e}] Pa (1–10 kPa). Per CLAUDE.md "
            "no-magic-number: keep within band or surface to PI."
        )
    rlo, rhi = _RATIO_LAMIN_BAND
    if not (rlo <= ratio_lamin <= rhi):
        raise ValueError(
            f"ratio_lamin = {ratio_lamin:.3f} outside the KU-3.B2.1 in-situ "
            f"nucleus:cytoplasm band [{rlo}, {rhi}]× (do NOT hard-code the "
            "10× isolated value — brief explicit). Surface to PI to change."
        )

    # E_nuc → k_eff bridge (eq. B): per-bead chromatin radial stiffness.
    k_chrom = 4.0 * math.pi * E_nuc * R_nuc / float(n_beads)
    # Lamin-shell increment from the nucleus:cytoplasm stiffness ratio.
    k_lamin = (ratio_lamin - 1.0) * k_chrom

    return ResolvedNucleus(
        k_chrom=k_chrom,
        k_lamin=k_lamin,
        R_nuc=float(R_nuc),
        n_beads=int(n_beads),
        d_knee=knee_strain * float(R_nuc),
        E_nuc=E_nuc,
        ratio_lamin=ratio_lamin,
        knee_strain=knee_strain,
        critical_pore_area_um2=critical_pore,
    )


# ---------------------------------------------------------------------------
# Two-regime radial force law (pure scalar — testable without a sim)
# ---------------------------------------------------------------------------
def confinement_force_magnitude(p: ResolvedNucleus, d: float) -> float:
    """Signed radial force magnitude for displacement ``d = r − R_nuc`` [N].

    BILINEAR strain-stiffening spring (continuous at the knee):

    * ``|d| ≤ d_knee``: ``F = −k_chrom · d``  (chromatin interior).
    * ``|d| > d_knee``: the lamin shell engages — the SLOPE steepens to
      ``k_chrom + k_lamin`` while the force VALUE stays continuous::

          F = −sign(d) · [ F_knee + (k_chrom + k_lamin)·(|d| − d_knee) ]

    Returns the radial force component along ``+r̂`` (negative = inward /
    restoring for ``d > 0``). This is the SAME law the per-bead custom force
    applies; exposed as a scalar so the strain-stiffening / continuity / sign
    Sanity-Gate checks run analytically (no sim).

    Args:
        p: Resolved nucleus parameters.
        d: Signed radial displacement from the shell, ``r − R_nuc`` [m].

    Returns:
        Radial force component along ``+r̂`` [N].
    """
    ad = abs(d)
    if ad <= p.d_knee:
        return -p.k_chrom * d
    sign = 1.0 if d > 0.0 else -1.0
    return -sign * (p.F_knee + (p.k_chrom + p.k_lamin) * (ad - p.d_knee))


def confinement_potential(p: ResolvedNucleus, d: float) -> float:
    """Non-negative bilinear confinement potential ``U(d)`` [J].

    The integral of ``−F_mag(d)`` (a restoring spring): a quadratic well of
    curvature ``k_chrom`` up to the knee, continued with curvature
    ``k_chrom + k_lamin`` beyond (C1: value AND slope of U continuous at the
    knee since ``F`` is continuous). ``U(0) = 0``; ``U ≥ 0`` everywhere.

    Args:
        p: Resolved nucleus parameters.
        d: Signed radial displacement [m].

    Returns:
        Confinement potential energy [J].
    """
    ad = abs(d)
    if ad <= p.d_knee:
        return 0.5 * p.k_chrom * d * d
    # Interior well up to the knee, plus the stiffened shell beyond.
    u_knee = 0.5 * p.k_chrom * p.d_knee * p.d_knee
    e = ad - p.d_knee
    return u_knee + p.F_knee * e + 0.5 * (p.k_chrom + p.k_lamin) * e * e


# ---------------------------------------------------------------------------
# KU-1.V.3.3 / KU-3.B2.3 confined-migration coupling (geometric predicate)
# ---------------------------------------------------------------------------
def nucleus_cross_section_area(positions: np.ndarray, normal_axis: int = 1) -> float:
    """Projected nuclear cross-section area ⟂ ``normal_axis`` [m²].

    Sphere-equivalent cross-section from the in-plane gyration radius of the
    nucleus-bead cloud projected onto the plane ⟂ to ``normal_axis`` (the
    confinement-pore transit axis). Used by the Lead's confined-migration
    driver to detect the KU-3.B2.3 / KU-1.V.3.3 ARREST condition (nuclear
    cross-section ≤ critical pore). NOTE-only coupling — the arrest physics
    (ECM pore as a geometric constraint) is a Lead integration concern.

    Args:
        positions: Nucleus-bead positions, shape ``(n, 3)`` [m].
        normal_axis: Pore-transit axis (0=x, 1=y, 2=z). Default 1 (``+ŷ``,
            the leading-edge advance axis the lamellipodium already uses).

    Returns:
        Projected cross-section area ``π R_g_plane²`` [m²]; 0 for < 1 bead.
    """
    pos = np.asarray(positions, dtype=np.float64).reshape(-1, 3)
    if pos.shape[0] == 0:
        return 0.0
    in_plane = [a for a in range(3) if a != int(normal_axis)]
    centroid = pos.mean(axis=0)
    d = pos[:, in_plane] - centroid[in_plane]
    # In-plane radius of gyration → sphere-equivalent projected radius.
    r_plane = float(np.sqrt((d * d).sum(axis=1).mean()))
    return math.pi * r_plane * r_plane


# ---------------------------------------------------------------------------
# Lamin-A scaling oracles (KU-3.B2.2; TEST/host helpers, NOT runtime)
# ---------------------------------------------------------------------------
def lamin_a_viscosity_scaling(lamin_a_rel: float, exponent: float = 3.0) -> float:
    """Relative nuclear-viscosity multiplier ~ [lamin-A]^(3±1) (KU-3.B2.2).

    Lamin-A → viscosity. Oracle helper for the host to fold lamin-A
    up/down-regulation into a viscosity scale; NOT called by the runtime.

    Args:
        lamin_a_rel: Lamin-A level relative to baseline (1.0 = baseline).
        exponent: Scaling exponent, KU-3.B2.2 ``3 ± 1`` (default 3.0).

    Returns:
        Relative viscosity multiplier ``lamin_a_rel ** exponent``.
    """
    _require_finite_positive("lamin_a_rel", lamin_a_rel)
    return float(lamin_a_rel) ** float(exponent)


def lamin_a_elasticity_scaling(lamin_a_rel: float, exponent: float = 0.5) -> float:
    """Relative nuclear-elasticity multiplier ~ [lamin-A]^(~0.5) (KU-3.B2.2).

    Oracle helper (the host may fold this into ``E_nuc`` for cancer lamin-A/C
    loss → softer nucleus). NOT called by the runtime.

    Args:
        lamin_a_rel: Lamin-A level relative to baseline (1.0 = baseline).
        exponent: Scaling exponent, KU-3.B2.2 ``~0.5`` (default 0.5).

    Returns:
        Relative elasticity multiplier ``lamin_a_rel ** exponent``.
    """
    _require_finite_positive("lamin_a_rel", lamin_a_rel)
    return float(lamin_a_rel) ** float(exponent)


# ---------------------------------------------------------------------------
# PURE bead-builder helper (computes & returns; NEVER mutates a sim/file)
# ---------------------------------------------------------------------------
def build_nucleus_beads(
    centroid: tuple[float, float, float] | np.ndarray,
    R_nuc: float,
    n_beads: int,
    *,
    gamma_nuc: float,
    type_name: str = "nucleus_bead",
    seed: int = 0,
    fill: bool = True,
) -> dict[str, object]:
    """Compute ``nucleus_bead`` snapshot arrays for the Lead to append.

    PURE helper (EXTEND Template 2): it computes positions / type-name /
    gamma arrays and RETURNS them. It MUST NOT mutate any simulation, state,
    snapshot, or file — the Lead appends these to the snapshot and the
    ``gamma_map`` later (single-writer ``cell.py``).

    Beads are seeded on (or filling) a sphere of radius ``R_nuc`` about
    ``centroid`` using a deterministic Fibonacci-sphere distribution (even
    angular coverage, no RNG clumping), so the confinement starts near its
    force-free configuration (``r ≈ R_nuc`` → ``d ≈ 0``). With ``fill=True``
    the radii are spread across a few shells filling the interior (a solid
    nuclear cloud); with ``fill=False`` all beads sit on the ``R_nuc`` shell.

    Args:
        centroid: Cell/nucleus centroid ``(x, y, z)`` [m].
        R_nuc: Nuclear radius [m] (> 0).
        n_beads: Number of nucleus beads (≥ 1).
        gamma_nuc: Per-bead Stokes drag [N·s/m] (= 6π η R_bead; > 0). Goes
            into the returned ``gamma`` array (one entry per bead) for the
            Lead to merge into the BAOAB ``gamma_map``.
        type_name: Particle type name (default ``"nucleus_bead"``).
        seed: Deterministic offset for the Fibonacci spiral (default 0).
        fill: If True, fill the interior across shells (solid cloud); if
            False, place all beads on the ``R_nuc`` shell.

    Returns:
        Dict with:
          * ``"positions"``: ndarray ``(n_beads, 3)`` float64 [m].
          * ``"type_name"``: str.
          * ``"gamma"``: ndarray ``(n_beads,)`` float64 [N·s/m].
          * ``"radii"``: ndarray ``(n_beads,)`` float64 [m] (seed radius
            from the centroid; diagnostic).

    Raises:
        ValueError: on non-positive ``R_nuc`` / ``gamma_nuc`` or ``n_beads<1``.
    """
    _require_finite_positive("R_nuc", R_nuc)
    _require_finite_positive("gamma_nuc", gamma_nuc)
    n = int(n_beads)
    if n < 1:
        raise ValueError(f"n_beads must be ≥ 1; got {n_beads!r}")

    c = np.asarray(centroid, dtype=np.float64).reshape(3)
    if not np.isfinite(c).all():
        raise ValueError(f"centroid must be finite; got {centroid!r}")

    # Deterministic Fibonacci-sphere directions (even angular coverage).
    idx = np.arange(n, dtype=np.float64) + 0.5
    phi = math.pi * (3.0 - math.sqrt(5.0))            # golden angle
    cos_theta = 1.0 - 2.0 * idx / float(n)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    sin_theta = np.sqrt(np.maximum(0.0, 1.0 - cos_theta * cos_theta))
    az = phi * idx + float(seed)
    dirs = np.stack(
        [sin_theta * np.cos(az), sin_theta * np.sin(az), cos_theta], axis=1
    )

    if fill and n > 1:
        # Spread radii so the cloud fills the interior with ~uniform volume
        # density (r ∝ k^{1/3}), the outermost bead landing on R_nuc.
        k = np.arange(1, n + 1, dtype=np.float64)
        radii = R_nuc * (k / float(n)) ** (1.0 / 3.0)
        # Decouple radius from the Fibonacci DIRECTION index. The spiral orders
        # directions monotonically in z (cos_theta = 1 - 2k/n), so leaving radii
        # in index order clusters small-radius beads at +z and large-radius beads
        # at -z; the radial confinement-force magnitudes then correlate with
        # direction and do NOT cancel -> a spurious NET force / COM drift (the
        # §3 internal-element / momentum-conservation violation, even though
        # set_forces references the instantaneous centroid). A deterministic
        # permutation decorrelates radius from direction -> isotropic cloud ->
        # net force ~ 0.
        radii = radii[np.random.default_rng(int(seed) + 1).permutation(n)]
    else:
        radii = np.full(n, float(R_nuc), dtype=np.float64)

    positions = c[None, :] + radii[:, None] * dirs

    return {
        "positions": positions.astype(np.float64),
        "type_name": str(type_name),
        "gamma": np.full(n, float(gamma_nuc), dtype=np.float64),
        "radii": radii.astype(np.float64),
    }


# ---------------------------------------------------------------------------
# Custom force compute
# ---------------------------------------------------------------------------
class NucleusConfinement(md.force.Custom):
    """Two-regime radial confinement on the nucleus-bead cloud (KU-3.B2).

    Each step:

    1. Read nucleus-bead positions (tags in ``nucleus_tag_range``).
    2. Centroid ``r_c`` of the nucleus beads (live).
    3. ``d_i = |r_i − r_c| − R_nuc``; radial force from the BILINEAR law
       :func:`confinement_force_magnitude` (continuous at the knee).
    4. ``F_i = F_mag(d_i) · r̂_i`` with ``r̂_i = (r_i − r_c)/|r_i − r_c|``.

    A non-negative confinement potential :func:`confinement_potential` is
    reported per-bead so HOOMD's ``.energy`` is meaningful.

    Mirrors ``cortex/erm.py`` (radial-spring SHAPE) but anchors normals at
    the BEAD CENTROID (like ``cortex/enclosed_volume.py``) so the
    confinement injects ≈ 0 net momentum: the nucleus is an INTERNAL element
    free to translate with the cell (Sanity Gate §3). Computed float64 from
    ``cpu_local_snapshot`` → ``cpu_local_force_arrays`` (row-indexed,
    single-rank, same guarantee erm.py documents).

    Args:
        p: Resolved nucleus parameters.
        nucleus_tag_range: Tag range ``[start, end)`` of the nucleus beads.
            Particles in this range receive the confinement; others
            untouched.
        aniso: HOOMD anisotropic flag (default False — central force only).
    """

    def __init__(
        self,
        p: ResolvedNucleus,
        nucleus_tag_range: tuple[int, int],
        *,
        aniso: bool = False,
    ) -> None:
        super().__init__(aniso=aniso)
        self.p = p
        self.tag_start = int(nucleus_tag_range[0])
        self.tag_end = int(nucleus_tag_range[1])
        if self.tag_end < self.tag_start:
            raise ValueError(
                f"nucleus_tag_range must satisfy end ≥ start; "
                f"got ({self.tag_start}, {self.tag_end})"
            )
        # Last-computed diagnostics (populated each set_forces).
        self.last_n: int = 0
        self.last_max_strain: float = float("nan")
        self.last_cross_section_area: float = float("nan")

    def critical_pore_area(self) -> float:
        """KU-1.V.3.3 critical confinement pore area [m²].

        ``critical_pore_area_um2`` (default 7 µm² tumor) converted to m².
        The Lead's confined-migration driver compares
        :func:`nucleus_cross_section_area` against this to detect arrest.
        NOTE-only coupling (no arrest physics enforced here).
        """
        return self.p.critical_pore_area_um2 * 1.0e-12

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()

        mask = (tag >= self.tag_start) & (tag < self.tag_end)
        F_vec = np.zeros_like(pos)
        U_per = np.zeros(pos.shape[0], dtype=np.float64)

        n_nuc = int(mask.sum())
        self.last_n = n_nuc
        if n_nuc > 0:
            nuc_pos = pos[mask]
            centroid = nuc_pos.mean(axis=0)
            dx = nuc_pos - centroid
            r = np.linalg.norm(dx, axis=1)
            r_safe = np.where(r > 0.0, r, 1.0)
            r_hat = dx / r_safe[:, None]
            r_hat[r <= 0.0] = 0.0          # singular centroid bead → no dir

            d = r - self.p.R_nuc           # signed radial displacement
            # Vectorised bilinear law (continuous at the knee).
            ad = np.abs(d)
            sign = np.sign(d)
            interior = ad <= self.p.d_knee
            F_mag = np.where(
                interior,
                -self.p.k_chrom * d,
                -sign * (
                    self.p.F_knee
                    + (self.p.k_chrom + self.p.k_lamin) * (ad - self.p.d_knee)
                ),
            )
            F_shell = F_mag[:, None] * r_hat
            F_vec[mask] = F_shell

            # Non-negative bilinear potential per bead.
            u_knee = 0.5 * self.p.k_chrom * self.p.d_knee * self.p.d_knee
            e = ad - self.p.d_knee
            U = np.where(
                interior,
                0.5 * self.p.k_chrom * d * d,
                u_knee
                + self.p.F_knee * e
                + 0.5 * (self.p.k_chrom + self.p.k_lamin) * e * e,
            )
            U_per[mask] = U

            self.last_max_strain = (
                float(ad.max() / self.p.R_nuc) if n_nuc else float("nan")
            )
            self.last_cross_section_area = nucleus_cross_section_area(nuc_pos)

        with self.cpu_local_force_arrays as arrays:
            arrays.force[:] = F_vec
            arrays.potential_energy[:] = U_per


# ---------------------------------------------------------------------------
# Public helper: attach nucleus confinement to an existing simulation
# ---------------------------------------------------------------------------
def attach_nucleus_confinement(
    sim: hoomd.Simulation,
    p_nuc: ResolvedNucleus,
    *,
    nucleus_tags: tuple[int, int],
    gamma_nuc: float,
    cfl_safety_factor: float = 0.1,
    cfl_strict: bool = True,
) -> NucleusConfinement:
    """Append a :class:`NucleusConfinement` custom force to an Integrator.

    CFL gate (mirrors ``cortex/erm.py`` / ``cortex/enclosed_volume.py``):
    the stiffest branch is ``k_hi = k_chrom + k_lamin``; the relax time is
    ``τ = γ_nuc / k_hi`` and ``dt`` must satisfy
    ``dt ≤ cfl_safety_factor · τ``. Raises ``RuntimeError`` if violated
    (using ``k_hi`` so the gate is conservative).

    Args:
        sim: Already-built simulation (must have an Integrator with ``dt``).
        p_nuc: Resolved nucleus parameters.
        nucleus_tags: Tag range ``[start, end)`` of the nucleus beads.
        gamma_nuc: Per-bead Stokes drag [N·s/m] (= 6π η R_bead; > 0). REQUIRED
            for the CFL gate (no default — the host supplies it).
        cfl_safety_factor: Same convention as the bond/angle/ERM CFL gates
            (D3 BAOAB). Default 0.1.
        cfl_strict: If True, raise on CFL violation; False for diagnostic
            runs (caller's responsibility).

    Returns:
        The attached :class:`NucleusConfinement` (kept for introspection).

    Raises:
        RuntimeError: if no Integrator is set, or the CFL gate is violated
            (and ``cfl_strict``).
        ValueError: on non-positive ``gamma_nuc``.
    """
    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "sim.operations.integrator must be set before attaching "
            "the nucleus confinement force."
        )
    _require_finite_positive("gamma_nuc", gamma_nuc)

    nuc = NucleusConfinement(p_nuc, nucleus_tags)

    k_hi = p_nuc.k_hi
    if k_hi > 0.0:
        dt = float(ig.dt)
        tau_nuc = gamma_nuc / k_hi
        dt_cfl_nuc = cfl_safety_factor * tau_nuc
        if dt > dt_cfl_nuc and cfl_strict:
            raise RuntimeError(
                f"Nucleus-confinement CFL violated: dt = {dt:.3e} s > "
                f"{cfl_safety_factor:.2f} · τ_nuc = {dt_cfl_nuc:.3e} s "
                f"(τ_nuc = γ_nuc / k_hi = {tau_nuc:.3e} s, k_hi = "
                f"k_chrom + k_lamin = {k_hi:.3e} N/m). Either soften the "
                "nucleus (lower E_nuc/ratio_lamin — requires PI sign-off vs "
                "the KU-3.B2.1 anchors) OR reduce dt (need dt ≤ "
                f"{dt_cfl_nuc:.3e} s). Pass cfl_strict=False to skip for "
                "diagnostic runs."
            )

    ig.forces.append(nuc)
    return nuc
