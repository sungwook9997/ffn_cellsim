"""H.8 plasma-membrane SURFACE mechanics — in-plane tension + area elasticity.

ADDITIVE, DEFAULT-OFF module (H.8 brief ``docs/briefs/H8_membrane_surface.md``;
EXTEND memo ``docs/CELL_MECHANICS_EXTEND_VS_REBUILD.md`` §B.4). The H.3 cortex
supplies a cortical tension γ_cortex, but the *plasma membrane* itself carries an
in-plane surface tension that the cortex shell does not represent. Without it, the
whole-cell apparent tension probed by micropipette aspiration / KU-3.5 is attributed
entirely to the cortex, which is physically wrong: the measured tension is
**composite**, ``γ_total = γ_membrane + γ_cortex`` (Sens & Plastino 2015; Fischer-
Friedrich 2014). This module supplies the missing γ_membrane term so KU-3.5 can be
re-derived as a composite, *without* touching the cortex modules or the integrator.

This file is the H.8 SURFACE module. It is DISTINCT from ``cell/membrane.py``
(the H.5 leading-edge load reservoir) — different concept, different concern.

Mechanism (fine-grained per CLAUDE.md, mirrors ``cortex/enclosed_volume.py``)
----------------------------------------------------------------------------
The same cortex shell beads (tags ``[0, n_cortex_actin)``) carry the plasma
membrane. Two in-plane contributions, both expressed as an inward Young-Laplace
pressure on the closed shell — the standard thin-shell relation ``ΔP = 2γ/R`` for a
sphere of radius ``R`` and surface tension ``γ`` (Laplace 1806; Evans & Skalak
1980):

1. **Bare in-plane tension** ``γ_mem`` (constant) → inward pressure ``2 γ_mem / R``.
2. **Area-expansion elasticity** with modulus ``K_A`` → an elastic tension that
   grows with areal strain, ``γ_area(A) = K_A · (A − A0) / A0`` (the in-plane area-
   stretch law; Evans & Rawicz 1990; Rawicz 2000). A stretched membrane (A > A0)
   develops a positive restoring tension that pulls the shell inward to shrink the
   area; a compressed membrane (A < A0) develops a negative (outward) tension.

These combine into a single composite membrane tension that is a function of the
instantaneous shell area::

    γ_tot(A) = γ_mem + K_A · (A − A0) / A0                                  [N/m]
    ΔP       = 2 · γ_tot(A) / R_mean                                       [Pa]
    F_i      = −ΔP · (S / N) · n̂_i                                         [N]

with ``S = 4π R_mean²`` the shell area, ``A_i = S/N`` the per-bead area share,
``n̂_i = (r_i − r_c)/|r_i − r_c|`` the OUTWARD radial unit vector about the shell
centroid ``r_c``, and the leading minus sign making a positive tension pull the
beads INWARD (a surface tension over a curved closed surface compresses it). This is
exactly the per-bead shell-force shape of ``cortex/enclosed_volume.py``; the *only*
difference is the pressure law (constant + area-elastic surface tension here, vs the
bulk/osmotic ΔP = −K_vol(V−V0)/V0 there). The two modules are complementary:
``enclosed_volume`` resists VOLUME change (cytoplasm incompressibility), this resists
AREA change / supplies the surface tension. Both attach to the same shell tags and
sum into HOOMD's net force.

Composite-tension reading (the H.8 purpose)
-------------------------------------------
Because both this membrane tension and the cortex's effective γ_cortex act as an
inward Laplace pressure on the SAME shell, the equilibrium Young-Laplace balance the
shell settles into reflects the SUM ``γ_total = γ_mem + γ_cortex``. The KU-3.5 tension
gate measures γ from exactly this Laplace balance (the cortex driver extracts an
effective γ from the shell's pressure–curvature relation / aspiration response — see
``outputs/h3/REPORT.md`` KU-3.5 protocol). Attaching this module therefore shifts the
measured γ by + γ_mem, making the reading composite. THIS MODULE DOES NOT CLAIM the
KU-3.5 gate passes — it supplies the missing additive term; re-validation is a Lead
integration step.

Helfrich bending κ_m — documented TODO (NOT implemented here)
-------------------------------------------------------------
The Helfrich curvature-elastic energy ``U_bend = ½ κ_m ∮ (2H − c0)² dA`` requires a
per-bead estimate of the local mean curvature ``H`` on an unstructured particle shell
(no fixed mesh connectivity at this ×40 coarse scale). A defensible discrete mean-
curvature operator (cotangent Laplacian / Meyer 2003, or a local quadric fit) is a
non-trivial later refinement and is **deliberately deferred** to keep this first
version correct-and-minimal (brief: "mean-curvature estimation on a particle shell is
non-trivial"). What IS provided is the bending modulus as a literature-anchored
constant and the **tether-force oracle** ``f_t = 2π √(2 κ_m (T_m + γ_MCA))`` (Hochmuth
1996; Derényi 2002) as a TEST-ONLY closed form, used to sanity-check that the κ_m / T
constants land the membrane tether force in the KU-3.B1 band 5–40 pN. κ_m enters NO
runtime force in this version. See :func:`tether_force`.

Magic-Number Block — physical constants
----------------------------------------
Every constant below is a KU-3.B1 literature anchor (see ``References`` and the
annotated module-level constants). NONE is tuned to make a gate pass:

* ``γ_mem`` (membrane in-plane tension) is a REQUIRED parameter of the resolver — it
  has NO default, because the whole point of H.8 is that the *decomposition*
  ``T = T_m + γ_MCA`` of the apparent tension is exactly what is being studied; baking
  in a default would prejudge it. The caller supplies ``T_m`` (membrane tension) from
  the KU-3.B1 band. (``DEFAULT_T_APPARENT`` / ``DEFAULT_GAMMA_MCA`` are provided as
  *annotated anchors* a caller may compose, not as a silent default.)
* ``K_A`` defaults to the KU-3.B1 area-expansion modulus 0.24 N/m (Rawicz 2000) — an
  intensive material modulus [N/m], grid-invariant (does not depend on N, ℓ₀, dt; the
  per-bead 1/N area share carries the count dependence so the TOTAL area force is
  N-independent — verified in test §CFL/grid).
* ``κ_m`` defaults to 1×10⁻¹⁹ J (≈ 24 k_BT; KU-3.B1, Rawicz 2000) — used only by the
  tether oracle in this version.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_membrane_surface.py``.*

1. **Dimensional analysis**
   - ``γ_mem``, ``γ_area``, ``K_A`` [N/m]; ``(A−A0)/A0`` dimensionless →
     ``γ_tot`` [N/m]. ``ΔP = 2 γ_tot / R`` [N/m / m = N/m² = Pa]. ``A_i = S/N``
     [m²]; ``F_i = ΔP·A_i`` [Pa·m² = N]. ✓
   - Tether ``f_t = 2π √(2 κ_m (T_m + γ_MCA))`` = √(J · (N/m)) = √(N·m · N/m) =
     √(N²) = N. ✓
   - A non-negative diagnostic potential is reported (constant-tension surface work
     ``γ_mem·(A−A0)`` + area-elastic ``½ K_A (A−A0)²/A0`` [J]) so ``.energy`` is
     meaningful; the per-bead force is the radial gradient of this surface work.

2. **Boundary cases**
   - ``A = A0`` AND ``γ_mem = 0``: ``γ_tot = 0`` → ΔP = 0 → F = 0 (the zero-force
     configuration). With ``γ_mem > 0`` the bare tension still pulls inward at
     A = A0 (a real membrane under tension on a curved surface is never force-free) —
     the area-elastic PART is what vanishes at A = A0. STATIC test asserts the
     area-only contribution is zero at A = A0.
   - Empty shell (n = 0): force loop over zero beads → trivially zero.
   - Bead at centroid (``r_i = r_c``): n̂ undefined → that bead's force set to 0
     (cortex beads sit at r ≈ R_cell ≫ 0; never hit in practice).
   - ``K_A ≤ 0``, ``R_cell ≤ 0``, ``A0 ≤ 0``, non-finite γ_mem: caught at resolve.

3. **Conservation invariants**
   - **Net force ≈ 0 (no spurious COM drift)**: n̂_i are centroid-relative, so for a
     symmetric shell ``Σ F_i = −ΔP·A·Σ n̂_i ≈ 0`` to the discretisation residual.
     This is an INTERNAL surface-tension potential and must not inject net momentum
     (same invariant, and same Fibonacci-lattice residual bound, as
     ``enclosed_volume.py`` §3). STATIC test asserts |Σ F_i| ≪ Σ|F_i|.
   - The force acts ONLY on tags in ``shell_tag_range``; all other particles
     untouched (mask by tag).

4. **Numerical sanity / CFL**
   - Positions read in float64; estimator + forces in float64.
   - Effective per-bead radial stiffness ``k_eff = (8π γ_mem + 16π K_A)/N²`` (derived
     in :meth:`MembraneSurfaceTension.effective_bead_stiffness`): both terms ∝ 1/N²,
     so for N ~ 7000 with γ_mem ~ 1e-4 N/m, K_A ~ 0.24 N/m this is ~1e-7 N/m → the
     relax time τ = γ_b/k_eff ≫ cortex dt_CFL. The membrane surface term is FAR
     softer than bonds/angles/crosslinkers and never tightens CFL — but the attach
     helper enforces the ``dt ≤ α·τ`` gate anyway (parity with erm.py / enclosed_
     volume.py). STATIC test confirms τ ≫ dt_cfl.

5. **Sign / sense**
   - ``γ_tot > 0`` (membrane under net tension): ``F_i ∝ −n̂_i`` → INWARD; a surface
     tension over a closed curved surface COMPRESSES it (Laplace). STATIC.
   - Stretched membrane ``A > A0`` (with γ_mem = 0): ``γ_area > 0`` → inward
     restoring force shrinks the area (area-expansion RESISTS stretch). STATIC.
   - Compressed ``A < A0`` (γ_mem = 0): ``γ_area < 0`` → outward, re-expands. STATIC.

6. **Measurement protocol**
   - KU-3.5 measures γ from the shell's Young-Laplace pressure–curvature balance
     (aspiration / pressure response). This module adds an inward Laplace pressure
     ``2 γ_mem / R`` of EXACTLY that form, so the γ it contributes is read on the same
     footing as γ_cortex → the measured tension becomes ``γ_mem + γ_cortex`` (test
     §composite). The tether oracle ``f_t`` (KU-3.B1.4) is checked to land in
     5–40 pN for the cited κ_m, T (test §tether). NO claim that KU-3.5 PASSES.

References
----------
- Brief: ``ffn_sim/docs/briefs/H8_membrane_surface.md`` (DRAFT skeleton; KU-3.B1).
- EXTEND memo: ``ffn_sim/docs/CELL_MECHANICS_EXTEND_VS_REBUILD.md`` §B.4 (attach
  contract, Template 1; composite γ_total = γ_membrane + γ_cortex).
- Structural analogs (verified in-tree template): ``ffn_sim/cortex/enclosed_volume.py``
  (EnclosedVolumePressure shell-force shape + centroid-relative normals + CFL gate)
  and ``ffn_sim/cortex/erm.py`` (md.force.Custom radial field + Sanity Gate + attach
  gate).
- KU-3.B1 anchors (apparent tension T 0.03–0.3 mN/m; K_A ≈ 0.24 N/m, lysis 3–10
  mN/m, max strain 2–5%; κ_m ≈ 1e-19 J / 10–30 k_BT; γ_MCA ≈ 1e-5 J/m²; tether
  f_t 5–40 pN): Rawicz et al. 2000 (Biophys. J. 79:328, bending & area modulus);
  Evans & Rawicz 1990 (PRL 64:2094, area-stretch & bending); Sens & Plastino 2015
  (J. Phys. Condens. Matter 27:273103, composite cortex+membrane tension); Fischer-
  Friedrich et al. 2014 (Sci. Rep. 4:6213, cortex tension); Hochmuth 1996 / Derényi
  et al. 2002 (PRL 88:238101, membrane tether force f_t = 2π√(2κ(T+γ_MCA))).
- HOOMD 7 ``md.force.Custom`` API.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import hoomd
import hoomd.md as md


# ---------------------------------------------------------------------------
# KU-3.B1 literature anchors (annotated module-level constants — NO silent
# defaults for γ_mem; see Magic-Number Block in the module docstring).
# ---------------------------------------------------------------------------
# Apparent (whole-cell) membrane tension band and its decomposition
# T = T_m + γ_MCA. KU-3.B1.1 (Sens & Plastino 2015; Fischer-Friedrich 2014).
DEFAULT_T_APPARENT: float = 3.0e-5      # N/m   apparent tension, band 0.03–0.3 mN/m
T_APPARENT_BAND: tuple[float, float] = (3.0e-5, 3.0e-4)  # N/m  (0.03–0.3 mN/m)

# Membrane–cortex adhesion energy density. KU-3.B1.3 (band 1e-6 .. 1e-4 J/m²).
DEFAULT_GAMMA_MCA: float = 1.0e-5       # J/m²  MCA adhesion
GAMMA_MCA_BAND: tuple[float, float] = (1.0e-6, 1.0e-4)   # J/m²

# Area-expansion (stretch) modulus. KU-3.B1.2 (Rawicz 2000); lysis 3–10 mN/m,
# max strain 2–5% → K_A·strain_max ≈ lysis tension.
DEFAULT_K_A: float = 0.24               # N/m   area-expansion modulus
LYSIS_TENSION_BAND: tuple[float, float] = (3.0e-3, 1.0e-2)  # N/m (3–10 mN/m)
MAX_AREAL_STRAIN_BAND: tuple[float, float] = (0.02, 0.05)   # 2–5 %

# Helfrich bending modulus (TODO stub — used only by the tether oracle here).
# KU-3.B1.2 (Rawicz 2000); 1e-19 J ≈ 24 k_BT at 300 K.
DEFAULT_KAPPA_M: float = 1.0e-19        # J     bending modulus, 10–30 k_BT

# Membrane tether-force acceptance band. KU-3.B1.4 (Hochmuth 1996; Derényi 2002).
TETHER_FORCE_BAND: tuple[float, float] = (5.0e-12, 4.0e-11)  # N  (5–40 pN)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedMembraneSurface:
    """Plasma-membrane surface-mechanics parameters (all SI).

    Attributes:
        gamma_mem: Bare in-plane membrane tension γ_mem [N/m]. REQUIRED (no
            default — the apparent-tension decomposition T = T_m + γ_MCA is the
            object of study; see Magic-Number Block). Typically the membrane part
            T_m of the KU-3.B1 apparent tension.
        K_A: Area-expansion (stretch) modulus [N/m]. KU-3.B1.2 (Rawicz 2000).
        A0: Reference (construction) shell area [m²]; γ_area = K_A·(A−A0)/A0.
        R_cell: Nominal shell radius [m] (sets A0 = 4π R_cell² when not given,
            and the CFL effective-stiffness scale).
        kappa_m: Helfrich bending modulus [J] (TODO stub; tether oracle only).
        gamma_mca: Membrane–cortex adhesion energy density [J/m²] (tether oracle).
    """

    gamma_mem: float                 # N/m   bare in-plane membrane tension (REQUIRED)
    K_A: float                       # N/m   area-expansion modulus, KU-3.B1.2
    A0: float                        # m²    reference (construction) shell area
    R_cell: float                    # m     nominal shell radius

    # κ_m / γ_MCA carried for the tether oracle + later Helfrich refinement;
    # NOT used by the runtime force in this version.
    kappa_m: float = DEFAULT_KAPPA_M       # J     bending modulus (TODO stub)
    gamma_mca: float = DEFAULT_GAMMA_MCA   # J/m²  membrane–cortex adhesion

    # Derived diagnostic.
    tether_force: float = 0.0        # N   f_t = 2π√(2 κ_m (γ_mem + γ_MCA))


def resolve_membrane_surface(
    cfg: dict,
    *,
    R_cell: float,
    A0: float | None = None,
) -> ResolvedMembraneSurface:
    """Resolve the ``membrane_surface`` config block.

    ``cfg`` may be the YAML root, a ``cell:`` sub-dict, or the
    ``membrane_surface`` sub-dict. ``R_cell`` comes from the host cortex resolved
    params. The reference area ``A0`` defaults to the sphere area ``4π R_cell²``
    (construction area) unless an explicit ``A0`` is supplied (e.g. measured from
    the construction snapshot).

    ``gamma_mem`` is REQUIRED in ``cfg`` (key ``gamma_mem`` or ``T_m``) — there is
    no default, per the Magic-Number Block (the apparent-tension decomposition is
    the object of study, so a silent default would prejudge it). ``K_A``, ``kappa_m``,
    ``gamma_mca`` default to their KU-3.B1 anchors.

    Raises:
        ValueError: if ``gamma_mem``/``T_m`` is absent, or any parameter is
            non-finite/non-positive (Sanity Gate §2 boundary checks).
    """
    if "cell" in cfg:
        cfg = cfg["cell"]
    if "membrane_surface" in cfg:
        cfg = cfg["membrane_surface"]

    # γ_mem is REQUIRED (accept either key name). No default.
    if "gamma_mem" in cfg:
        gamma_mem = float(cfg["gamma_mem"])
    elif "T_m" in cfg:
        gamma_mem = float(cfg["T_m"])
    else:
        raise ValueError(
            "membrane_surface requires an explicit 'gamma_mem' (or 'T_m') "
            "membrane in-plane tension [N/m]; no default is provided because "
            "the apparent-tension decomposition T = T_m + γ_MCA (KU-3.B1) is "
            "the object of study. Supply T_m from the KU-3.B1 band "
            f"{T_APPARENT_BAND} N/m."
        )

    K_A = float(cfg.get("K_A", DEFAULT_K_A))
    kappa_m = float(cfg.get("kappa_m", DEFAULT_KAPPA_M))
    gamma_mca = float(cfg.get("gamma_mca", DEFAULT_GAMMA_MCA))

    # §2 boundary checks.
    _require_finite_positive("R_cell", R_cell)
    _require_finite_positive("K_A", K_A)
    _require_finite_positive("kappa_m", kappa_m)
    _require_finite_positive("gamma_mca", gamma_mca)
    if not math.isfinite(gamma_mem):
        raise ValueError(f"gamma_mem must be finite; got {gamma_mem!r}")
    if gamma_mem < 0.0:
        # A negative bare tension is unphysical for a resting membrane.
        raise ValueError(f"gamma_mem must be ≥ 0; got {gamma_mem!r}")

    if A0 is None:
        A0 = 4.0 * math.pi * R_cell**2
    _require_finite_positive("A0", A0)

    p = ResolvedMembraneSurface(
        gamma_mem=gamma_mem,
        K_A=K_A,
        A0=A0,
        R_cell=R_cell,
        kappa_m=kappa_m,
        gamma_mca=gamma_mca,
    )
    # Derived diagnostic: membrane tether force (KU-3.B1.4 oracle).
    p.tether_force = tether_force(kappa_m, gamma_mem + gamma_mca)
    return p


# ---------------------------------------------------------------------------
# Closed-form oracles (TEST-ONLY helpers; NOT used by the runtime force)
# ---------------------------------------------------------------------------
def membrane_tension(p: ResolvedMembraneSurface, A: float) -> float:
    """Composite membrane surface tension at shell area ``A`` [N/m].

    ``γ_tot(A) = γ_mem + K_A·(A − A0)/A0`` — the SAME law the runtime per-bead
    force uses, exposed as a scalar so tests can check it analytically. This is the
    bare in-plane tension plus the in-plane area-stretch tension (Evans & Rawicz
    1990; Rawicz 2000).
    """
    return p.gamma_mem + p.K_A * (A - p.A0) / p.A0


def laplace_pressure(gamma: float, R: float) -> float:
    """Young-Laplace inward pressure for a sphere: ΔP = 2γ/R [Pa].

    ACCEPTANCE ORACLE (KU-3.5 measurement form). Used in tests to confirm the
    runtime per-bead force realises this pressure. NEVER called by the runtime.
    """
    return 2.0 * gamma / R


def tether_force(kappa_m: float, tension: float) -> float:
    """Membrane tether (nanotube) pulling force [N] — KU-3.B1.4 oracle.

    ``f_t = 2π √(2 κ_m (T_m + γ_MCA))`` (Hochmuth 1996; Derényi et al. 2002), the
    equilibrium force to hold a membrane tether against bending κ_m and the
    effective tension ``tension = T_m + γ_MCA``. TEST-ONLY closed form: κ_m enters
    NO runtime force in this version (Helfrich bending is a documented TODO stub).

    Args:
        kappa_m: Bending modulus [J].
        tension: Effective tension ``T_m + γ_MCA`` [N/m].

    Returns:
        Tether force [N]; lands in the KU-3.B1 band 5–40 pN for the cited anchors.
    """
    return 2.0 * math.pi * math.sqrt(2.0 * kappa_m * tension)


# ---------------------------------------------------------------------------
# Custom force compute
# ---------------------------------------------------------------------------
class MembraneSurfaceTension(md.force.Custom):
    """Per-bead plasma-membrane surface-tension force on the cortex shell (H.8).

    Each step:

    1. Read shell-bead positions (tags in ``shell_tag_range``).
    2. Centroid ``r_c`` of the shell beads.
    3. ``R_mean = ⟨|r_i − r_c|⟩``; shell area ``S = 4π R_mean²``.
    4. Composite tension ``γ_tot = γ_mem + K_A·(S − A0)/A0``; inward Laplace
       pressure ``ΔP = 2 γ_tot / R_mean``.
    5. ``F_i = −ΔP · (S / N) · n̂_i`` with ``n̂_i = (r_i − r_c)/|r_i − r_c|`` the
       outward radial unit vector (the leading minus makes positive tension pull
       INWARD).

    A non-negative diagnostic surface-work potential ``γ_mem·(S − A0) +
    ½ K_A (S − A0)²/A0`` (clamped at 0 for the small negative constant-tension part
    when S < A0) is reported onto the shell beads so ``.energy`` is meaningful.

    Args:
        p: Resolved membrane-surface config.
        shell_tag_range: Tag range ``[start, end)`` of the cortex-actin shell beads.
            Particles whose tag is in this range receive the membrane force; all
            others are untouched.
        aniso: Passed to ``md.force.Custom`` (isotropic; default False).

    Notes:
        Implemented as ``md.force.Custom`` (mirrors ``cortex/erm.py`` and
        ``cortex/enclosed_volume.py``) so the contribution participates in HOOMD's
        ``net_force`` accumulator, which the L-M BAOAB Updater reads each step.
        Positions read via ``cpu_local_snapshot`` and forces written via
        ``cpu_local_force_arrays`` — both ROW-indexed; HOOMD guarantees row ``i`` is
        the same particle in both on a single rank.

        Centroid-relative normals keep the net surface-tension force ≈ 0 (Sanity
        Gate §3): this is an INTERNAL surface potential and must not inject net
        momentum, unlike the ERM external field.

        κ_m (Helfrich bending) is NOT used here — documented TODO stub (see module
        docstring); only γ_mem + K_A enter the runtime force.
    """

    def __init__(
        self,
        p: ResolvedMembraneSurface,
        shell_tag_range: tuple[int, int],
        *,
        aniso: bool = False,
    ) -> None:
        super().__init__(aniso=aniso)
        self.p = p
        self.tag_start = int(shell_tag_range[0])
        self.tag_end = int(shell_tag_range[1])
        if self.tag_end < self.tag_start:
            raise ValueError(
                f"shell_tag_range must satisfy end ≥ start; "
                f"got ({self.tag_start}, {self.tag_end})"
            )
        # Last-computed diagnostics (populated each set_forces).
        self.last_area: float = float("nan")
        self.last_tension: float = float("nan")
        self.last_pressure: float = float("nan")
        self.last_R_mean: float = float("nan")

    # -- estimator (static so tests can call it on a bare position array) --
    @staticmethod
    def estimate_area(
        positions: np.ndarray,
    ) -> tuple[float, float, np.ndarray, np.ndarray]:
        """Sphere-equivalent shell area + radius from a shell-bead array.

        Mirrors ``EnclosedVolumePressure.estimate_volume`` (same centroid-relative,
        translation-invariant estimator) but returns the AREA (this module resists
        area change, not volume).

        Args:
            positions: Shell-bead positions [m], shape (n, 3).

        Returns:
            A tuple ``(S, R_mean, centroid, radii)`` where ``S = 4π R_mean²`` is the
            shell area [m²], ``R_mean = ⟨|r_i − r_c|⟩`` [m], ``centroid`` is the bead
            centroid [m] (shape (3,)), and ``radii`` are per-bead distances from the
            centroid [m] (shape (n,)).
        """
        pos = np.asarray(positions, dtype=np.float64).reshape(-1, 3)
        n = pos.shape[0]
        if n == 0:
            return 0.0, 0.0, np.zeros(3, dtype=np.float64), np.empty(0)
        centroid = pos.mean(axis=0)
        radii = np.linalg.norm(pos - centroid, axis=1)
        R_mean = float(radii.mean())
        S = 4.0 * math.pi * R_mean**2
        return S, R_mean, centroid, radii

    def effective_bead_stiffness(self, n_shell: int) -> float:
        """Per-bead radial stiffness of the membrane surface term [N/m] (CFL gate).

        Moving ONE shell bead radially by δr changes ``R_mean`` by δr/N, hence the
        shell area ``S = 4π R_mean²`` by ``dS = 8π R_mean·(δr/N)`` and the inward
        pressure force on that bead. Linearising the per-bead force
        ``F_i = −2 γ_tot(S)/R_mean · (S/N)`` about the reference sphere
        (R_mean = R_cell, S = A0) gives the dominant radial stiffnesses:

            constant-tension part:  k_γ = 8π γ_mem / N²
            area-elastic part:      k_A = 16π K_A / N²
            k_eff = (8π γ_mem + 16π K_A) / N²

        Both ∝ 1/N² → tiny for N ~ 7000 (≪ bond/angle/crosslinker stiffness), so the
        membrane surface term never tightens CFL (verified in test §CFL). The attach
        helper enforces ``dt ≤ α·γ_b/k_eff`` anyway for parity with erm.py.

        Args:
            n_shell: Number of shell beads N.

        Returns:
            Effective per-bead radial stiffness [N/m].
        """
        N = float(n_shell)
        return (8.0 * math.pi * self.p.gamma_mem
                + 16.0 * math.pi * self.p.K_A) / (N * N)

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()

        mask = (tag >= self.tag_start) & (tag < self.tag_end)
        F_vec = np.zeros_like(pos)
        U_per = np.zeros(pos.shape[0], dtype=np.float64)

        n_shell = int(mask.sum())
        if n_shell > 0:
            shell_pos = pos[mask]
            S, R_mean, centroid, radii = self.estimate_area(shell_pos)
            if R_mean > 0.0:
                gamma_tot = self.p.gamma_mem + self.p.K_A * (S - self.p.A0) / self.p.A0
                dP = 2.0 * gamma_tot / R_mean          # inward Laplace pressure [Pa]
                # Per-bead outward normal about the centroid.
                dx = shell_pos - centroid
                r_safe = np.where(radii > 0.0, radii, 1.0)
                n_hat = dx / r_safe[:, None]
                n_hat[radii <= 0.0] = 0.0
                A_i = S / n_shell
                # Minus sign: positive tension → INWARD (−n̂) compressive force.
                F_shell = (-(dP * A_i)) * n_hat
                F_vec[mask] = F_shell

                # Diagnostic non-negative surface-work potential, distributed evenly:
                #   constant-tension work γ_mem·(S − A0) (clamped ≥ 0)
                # + area-elastic        ½ K_A (S − A0)²/A0  (always ≥ 0).
                U_const = max(self.p.gamma_mem * (S - self.p.A0), 0.0)
                U_area = 0.5 * self.p.K_A * (S - self.p.A0) ** 2 / self.p.A0
                U_total = U_const + U_area
                U_per[mask] = U_total / n_shell

                self.last_area = S
                self.last_tension = gamma_tot
                self.last_pressure = dP
                self.last_R_mean = R_mean

        with self.cpu_local_force_arrays as arrays:
            arrays.force[:] = F_vec
            arrays.potential_energy[:] = U_per


# ---------------------------------------------------------------------------
# Public helper: attach membrane surface tension to an existing simulation
# ---------------------------------------------------------------------------
def attach_membrane_surface(
    sim: hoomd.Simulation,
    p_mem: ResolvedMembraneSurface,
    *,
    shell_tag_range: tuple[int, int],
    gamma_b: float | None = None,
    n_shell: int | None = None,
    cfl_safety_factor: float = 0.1,
    cfl_strict: bool = True,
) -> MembraneSurfaceTension:
    """Append a MembraneSurfaceTension custom force to an existing Integrator.

    Mirrors ``attach_enclosed_volume_to_simulation`` / ``attach_erm_to_simulation``:
    builds the force over ``shell_tag_range`` and applies the CFL gate
    ``dt ≤ cfl_safety_factor · γ_b / k_eff`` before appending. DEFAULT-OFF contract:
    if the Lead does NOT call this helper, nothing is attached and existing runs are
    bit-for-bit identical (this module registers nothing globally on import).

    Args:
        sim: Already-built cortex simulation (must have an Integrator).
        p_mem: Resolved membrane-surface config.
        shell_tag_range: Tag range ``[start, end)`` of the cortex-actin shell beads.
        gamma_b: Per-bead Stokes drag [N·s/m]. Required to gate the membrane CFL
            ``dt ≤ cfl_safety_factor · γ_b / k_eff``. If None, the CFL gate is
            skipped (caller responsible). The membrane surface k_eff is extremely
            soft (≪ bond/angle/xl stiffness), so this gate essentially never fires —
            enforced for parity with erm.py.
        n_shell: Number of shell beads (= ``shell_tag_range[1] − shell_tag_range[0]``
            if omitted). Used for the CFL per-bead stiffness.
        cfl_safety_factor: Same convention as the bond/angle/ERM/EV CFL gates
            (D3 BAOAB). Default 0.1.
        cfl_strict: If True, raise on CFL violation; set False for diagnostic runs.

    Returns:
        The attached :class:`MembraneSurfaceTension` force compute (kept for
        introspection in tests).

    Raises:
        RuntimeError: if ``sim.operations.integrator`` is unset, or the CFL gate is
            violated and ``cfl_strict`` is True.
    """
    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "sim.operations.integrator must be set before attaching "
            "the membrane surface tension force."
        )
    if n_shell is None:
        n_shell = int(shell_tag_range[1]) - int(shell_tag_range[0])

    mem = MembraneSurfaceTension(p_mem, shell_tag_range)

    if gamma_b is not None and n_shell > 0:
        dt = float(ig.dt)
        k_eff = mem.effective_bead_stiffness(n_shell)
        if k_eff > 0.0:
            tau_mem = gamma_b / k_eff
            dt_cfl_mem = cfl_safety_factor * tau_mem
            if dt > dt_cfl_mem and cfl_strict:
                raise RuntimeError(
                    f"Membrane-surface CFL violated: dt = {dt:.3e} s > "
                    f"{cfl_safety_factor:.2f} · τ_mem = {dt_cfl_mem:.3e} s "
                    f"(τ_mem = γ_b / k_eff = {tau_mem:.3e} s, k_eff = "
                    f"{k_eff:.3e} N/m). Reduce γ_mem / K_A (softer membrane — "
                    "requires PI sign-off vs the KU-3.B1 anchors) OR reduce dt "
                    f"(need dt ≤ {dt_cfl_mem:.3e} s). Pass cfl_strict=False to "
                    "skip for diagnostic runs."
                )

    ig.forces.append(mem)
    return mem
