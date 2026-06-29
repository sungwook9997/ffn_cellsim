"""H.3 enclosed-volume / intracellular-pressure term (KU-3.1).

ADDITIVE, DEFAULT-OFF module (mechanism audit 2026-05-30 §5 #2;
``docs/briefs/H4_FA_INTEGRATION_DESIGN.md`` §4). The cortex is otherwise
an OPEN shell of tangent filaments + radial ERM pins, so it can be
deflated to a point at no energetic cost — KU-3.1 mitotic rounding is
therefore non-mechanistic. Real cells maintain an intracellular
hydrostatic pressure ΔP (osmotic + actomyosin), and rounding is driven by
the Young-Laplace balance ΔP = 2γ/R between this pressure and the cortical
tension γ. This module supplies the missing pressure term so a closed
cortex shell resists volume change and ΔP = 2γ/R can EMERGE.

Mechanism (fine-grained per CLAUDE.md)
--------------------------------------
The cortex shell beads enclose a volume V. Define a reference volume
``V0`` (the construction volume). The pressure is a linear bulk/osmotic
response::

    ΔP = −K_vol · (V − V0) / V0

with ``K_vol`` a cytoplasm bulk-response coefficient (Pa).

**Baseline osmotic turgor (PI-ratified 2026-06-04).** A real rounded cell is an
osmotically-PRESSURISED drop: it carries a resting intracellular hydrostatic
pressure ``Π₀`` (Stewart 2011: ~40 Pa interphase → ~400 Pa metaphase) that
pre-tensions the cortex via Young-Laplace ``γ ≈ Π₀·R/2`` BEFORE any myosin.
The full pressure is therefore::

    ΔP = Π₀ − K_vol · (V − V0) / V0

so at the construction volume (V = V0) ``ΔP = Π₀ > 0`` — the cortex bonds carry
a hoop tension that balances the outward turgor (the physically-correct resting
state). Myosin then MODULATES on top (Bohec 2026). ``turgor_dP0 = 0`` (default)
recovers the prior force-free-at-construction behaviour (additive / opt-in).
Distribute the pressure as a normal force over the shell beads::

    F_i = ΔP · A_i · n̂_i

with ``A_i = S(V) / N`` the per-bead share of the shell surface area and
``n̂_i`` the OUTWARD radial unit vector from the cell centroid. Positive
ΔP (V < V0, cell compressed) → outward (restoring) force; negative ΔP
(V > V0, cell inflated) → inward (restoring) force.

The closed-form Young-Laplace ΔP = 2γ/R is the ACCEPTANCE ORACLE used in
``ffn_sim/tests/test_enclosed_volume.py`` only — NEVER the runtime. The
runtime is the explicit per-bead force above.

Volume estimator
----------------
At this ×40 coarse scale the cortex is a roughly-spherical bead cloud, so
the simplest defensible estimator is the **sphere-equivalent volume from
the mean bead radius about the centroid**::

    R_mean = ⟨|r_i − r_c|⟩      (mean over shell beads)
    V      = (4/3) π R_mean³
    S      = 4 π R_mean²         (shell area, for the per-bead share A_i)

where ``r_c`` is the centroid of the shell beads (the COM-free estimator
makes V translation-invariant, so the pressure force introduces NO
spurious net force / COM drift — see Sanity Gate §3). This is documented
in the FA brief §4 as acceptable for the spherical cortex; a convex-hull
estimator is a later refinement (not needed for the spherical mitotic
shell and far costlier per step). ``R_mean`` (mean) rather than an RMS or
max radius keeps the estimator a smooth, well-conditioned function of the
bead positions and makes the Young-Laplace recovery (test §Young-Laplace)
exact for a perfect sphere.

Magic-Number Block — K_vol (no empirical magic number)
------------------------------------------------------
``K_vol`` is an OSMOTIC-MODULUS anchor, NOT a tuned spring. It is fixed so
the construction-pressure scale matches the KU-3.1 interphase osmotic
anchor ΔP ≈ 40 Pa (Stewart 2011) at a 1 % radius compression:

* Derivable: for a spherical shell a fractional radius change δR/R gives a
  fractional volume change ``δV/V0 = 3 · δR/R`` (V ∝ R³, linearised). The
  bulk law ΔP = −K_vol · δV/V0 then gives, at a compression δR/R = 0.01
  (δV/V0 = 0.03):

      K_vol = ΔP_anchor / (δV/V0) = ΔP_anchor / (3 · 0.01)
            = 40 Pa / 0.03 ≈ 1.333 × 10³ Pa.

  So ``K_vol = ΔP_ref / (3 · strain_ref)`` with ΔP_ref = 40 Pa (Stewart
  2011 interphase) and strain_ref = 0.01 (1 % radius compression — the
  linear-regime reference point, not a fit). Both are stated literature /
  geometric anchors; the resolver recomputes K_vol from them.
* Grid-invariant: K_vol is an intensive material modulus [Pa]. It does
  NOT depend on n_filaments, beads_per_filament, ℓ₀, or dt — the per-bead
  force F_i = ΔP · (S/N) · n̂ carries the 1/N share so the TOTAL outward
  pressure force Σ|F_i| = ΔP · S is N-independent (the physical
  pressure × area). Doubling the bead count halves A_i and doubles the
  count → identical net pressure load. Verified in test §6.
* Not chosen to pass a gate: the Young-Laplace test (§Young-Laplace) does
  NOT use K_vol at all — it recovers the equilibrium ΔP from the
  force-balance and checks it equals 2γ/R for an externally imposed γ.
  K_vol only sets WHERE that equilibrium sits, not the law it obeys.

Cross-check (consistency, not a fit): at the cortex default R_cell =
10 μm and KU-3.5 cortical tension γ = 1 mN/m (Salbreux 2012), the
Young-Laplace pressure is ΔP = 2γ/R = 2·1e-3/1e-5 = 200 Pa — squarely
between the 40 Pa interphase and 400 Pa metaphase anchors (Stewart 2011;
Fischer-Friedrich 2014). The K_vol anchor and the γ anchor are mutually
consistent.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_enclosed_volume.py``.*

1. **Dimensional analysis**
   - ``K_vol`` [Pa = N/m²]; ``(V−V0)/V0`` dimensionless → ΔP [Pa]. ✓
   - ``A_i = S/N`` [m²]; ``F_i = ΔP·A_i`` [Pa·m² = N]. ✓
   - ``U`` is the bulk elastic energy ``½ (K_vol/V0) (V−V0)²`` [Pa·m³ = J]
     (reported as a diagnostic potential; the per-bead force is the
     gradient of this pressure work to leading order). ✓
   - At construction (V = V0): ΔP = 0 → F = 0, U = 0. Force-free at
     construction.

2. **Boundary cases**
   - ``K_vol ≤ 0``: caught by ``resolve_enclosed_volume`` ValueError.
   - ``V0 ≤ 0`` / ``R_cell ≤ 0``: caught at resolve time.
   - Empty shell (n beads = 0): the custom force iterates zero shell
     particles; ΔP from V=0 is clamped (R_mean=0 → V=0 → guarded), force
     trivially zero. (In practice the cortex always has beads.)
   - A bead AT the centroid (r_i = r_c): n̂ undefined; we set that bead's
     pressure force to 0 (it contributes 0 area share direction). Cortex
     beads sit at r ≈ R_cell ≫ 0 so this is never hit in practice.
   - ``strain_ref`` ≤ 0 or ≥ 1, or ``dP_ref`` ≤ 0: caught at resolve.

3. **Conservation invariants**
   - **Net force ≈ 0 (no spurious COM drift)**: the outward radial unit
     vectors n̂_i are computed about the bead CENTROID r_c. For a
     symmetric (centroid-centered) shell Σ n̂_i ≈ 0, so Σ F_i = ΔP·A·Σn̂_i
     ≈ 0 to the discretisation residual. This is the key difference from
     the ERM field (which is anchored at a FIXED origin and intentionally
     breaks translational symmetry): the pressure term is INTERNAL
     (gradient of an enclosed-volume potential) and must not inject net
     momentum. STATIC test asserts |Σ F_i| ≪ Σ|F_i|.
   - **Sign of work**: ΔP·δV ≤ 0 about V0 (restoring) — the potential is a
     non-negative quadratic in (V−V0). STATIC.
   - The pressure force acts ONLY on particles whose tag is in
     ``shell_tag_range`` (cortex actin). Other particles untouched.

4. **Numerical sanity**
   - All positions read in float64; estimator + forces computed in
     float64.
   - CFL: the pressure term has an effective per-bead radial stiffness
     ``k_eff = ΔP-slope w.r.t. bead radius``. For the sphere-equivalent
     estimator, moving ONE bead radially by δr changes R_mean by δr/N,
     hence V by (S/N)·(δr/N)·... → the per-bead restoring stiffness is
     ``k_eff = 3 K_vol S² / (N² V0)`` (derived in
     :meth:`EnclosedVolumePressure.effective_bead_stiffness`). The relax
     time ``τ_vol = γ_b / k_eff`` must satisfy ``dt ≤ α · τ_vol``; the
     attach helper raises (like erm.py) if violated. With the default
     K_vol ≈ 1.3e3 Pa, R=10 μm, N=7000, this k_eff is ~1e-12 N/m → τ_vol
     ~ms ≫ cortex dt_CFL=13 ns, so the pressure term NEVER tightens CFL
     (it is far softer than even the crosslinkers). STATIC test confirms.

5. **Sign / sense**
   - V < V0 (cell compressed): (V−V0) < 0 → ΔP = −K_vol·(neg)/V0 > 0 →
     F_i = ΔP·A·n̂ points OUTWARD (+n̂): restoring (re-inflates). STATIC.
   - V > V0 (cell inflated): ΔP < 0 → F_i points INWARD (−n̂): restoring
     (re-deflates). STATIC.
   - V = V0: ΔP = 0, F = 0. STATIC.

6. **Measurement protocol**
   - **Young-Laplace recovery (the KU-3.1 mechanism oracle)**: place a
     perfect sphere of radius R; impose an inward cortical-tension-like
     line load equivalent to a surface tension γ (the test applies the
     analytic Laplace inward pressure 2γ/R as an external balancing
     pressure). The enclosed-volume term's outward ΔP balances it at the
     volume where ``K_vol·(V0−V)/V0 = 2γ/R``; the recovered ΔP equals
     2γ/R within tolerance. ANALYTICAL test (no long sim needed).
   - **BAOAB smoke**: cortex + enclosed-volume-on runs a few hundred
     steps with no NaN/Inf and the enclosed volume stays bounded near V0.

References
----------
- Mechanism audit: ``ffn_sim/docs/v2_audit/MECHANISM_AUDIT_2026-05-30.md``
  §5 #2 (intracellular pressure / enclosed-volume — the actual KU-3.1
  mechanism) + §6 row 13 (mitotic rounding ΔP=2γ/R).
- Design: ``ffn_sim/docs/briefs/H4_FA_INTEGRATION_DESIGN.md`` §4.
- KU-3.1 osmotic / pressure anchor: Stewart et al. 2011 (Nature,
  "Hydrostatic pressure and the actomyosin cortex drive mitotic cell
  rounding"), ΔP ≈ 40 Pa interphase → ~400 Pa metaphase;
  Fischer-Friedrich et al. 2014.
- KU-3.5 cortical tension γ ≈ 1 mN/m (Salbreux 2012) — the consistency
  cross-check above and the test oracle γ.
- ``ffn_sim/cortex/erm.py`` — the structural analog (md.force.Custom
  radial field on cortex beads + §1-6 Sanity Gate + CFL attach gate).
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
class ResolvedEnclosedVolume:
    """Enclosed-volume / intracellular-pressure parameters (all SI)."""

    K_vol: float                 # Pa   cytoplasm bulk-response coefficient
    V0: float                    # m³   reference (construction) enclosed volume
    R_cell: float                # m    nominal shell radius (V0 = 4/3 π R³)

    # Provenance of K_vol (Magic-Number Block anchors; recomputed → K_vol)
    dP_ref: float = 40.0         # Pa   KU-3.1 interphase anchor (Stewart 2011)
    strain_ref: float = 0.01     # —    1 % radius compression (linear ref pt)

    # Baseline osmotic TURGOR Π₀ [Pa] — the RESTING intracellular hydrostatic
    # pressure (Stewart 2011: ~40 Pa interphase → ~400 Pa metaphase). With
    # turgor_dP0 > 0 the shell is PRE-TENSIONED at construction (ΔP = Π₀ at
    # V=V0), so the cortex carries a Young-Laplace hoop tension γ ≈ Π₀·R/2
    # WITHOUT myosin — the physically-correct rounded-cell state (a real cell is
    # an osmotically-pressurised drop; myosin MODULATES on top, Bohec 2026).
    # Default 0.0 reproduces the prior force-free-at-construction behaviour
    # (additive/opt-in per CLAUDE.md). Anchored to the SAME Stewart 2011 datum
    # the K_vol block already cites — NOT a tuned value (PI-ratified 2026-06-04).
    turgor_dP0: float = 0.0      # Pa   baseline osmotic turgor (resting ΔP)

    # Derived diagnostic
    dP_at_strain_ref: float = 0.0   # ΔP at strain_ref (≈ dP_ref by construction)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def resolve_enclosed_volume(
    cfg: dict,
    *,
    R_cell: float,
    V0: float | None = None,
) -> ResolvedEnclosedVolume:
    """Resolve the ``enclosed_volume`` config block.

    ``cfg`` is the YAML root, the ``cortex:`` sub-dict, or the
    ``cortex.enclosed_volume`` sub-dict. ``R_cell`` comes from the host
    cortex resolved params; the reference volume ``V0`` defaults to the
    sphere-equivalent ``4/3 π R_cell³`` (the construction volume) unless an
    explicit ``V0`` is supplied (e.g. measured from the construction
    snapshot centroid).

    K_vol is recomputed from first principles (Magic-Number Block):
    ``K_vol = dP_ref / (3 · strain_ref)`` — the osmotic-modulus anchor that
    makes ΔP at a ``strain_ref`` radius compression equal the KU-3.1
    interphase anchor ``dP_ref`` (Stewart 2011). If the config supplies an
    explicit ``K_vol`` it is taken verbatim and checked for consistency
    against the anchors (a >1 % mismatch raises — guards against a tuned
    value silently overriding the derivation).
    """
    if "cortex" in cfg:
        cfg = cfg["cortex"]
    if "enclosed_volume" in cfg:
        cfg = cfg["enclosed_volume"]

    dP_ref = float(cfg.get("dP_ref", 40.0))
    strain_ref = float(cfg.get("strain_ref", 0.01))
    # Baseline osmotic turgor (opt-in; 0 = prior force-free-at-construction).
    turgor_dP0 = float(cfg.get("turgor_dP0", 0.0))

    # §2 boundary checks on anchors.
    _require_finite_positive("R_cell", R_cell)
    _require_finite_positive("dP_ref", dP_ref)
    if not (math.isfinite(turgor_dP0) and turgor_dP0 >= 0.0):
        raise ValueError(f"turgor_dP0 must be finite and >= 0; got {turgor_dP0!r}")
    if not (math.isfinite(strain_ref) and 0.0 < strain_ref < 1.0):
        raise ValueError(
            f"strain_ref must be finite and in (0, 1); got {strain_ref}"
        )

    # First-principles K_vol from the Magic-Number Block.
    #   δV/V0 = 3 · δR/R = 3 · strain_ref  (V ∝ R³, linearised)
    #   ΔP = K_vol · δV/V0  →  K_vol = dP_ref / (3 · strain_ref)
    K_vol_derived = dP_ref / (3.0 * strain_ref)

    K_vol_cfg = cfg.get("K_vol", None)
    if K_vol_cfg is not None:
        K_vol = float(K_vol_cfg)
        _require_finite_positive("K_vol", K_vol)
        # Consistency guard: an explicit K_vol must match the anchors.
        if not math.isclose(K_vol, K_vol_derived, rel_tol=1.0e-2):
            raise ValueError(
                f"Explicit K_vol = {K_vol:.4e} Pa disagrees with the "
                f"first-principles anchor dP_ref/(3·strain_ref) = "
                f"{K_vol_derived:.4e} Pa (dP_ref={dP_ref} Pa, "
                f"strain_ref={strain_ref}). Per CLAUDE.md no-magic-number: "
                "either remove the explicit K_vol (let it derive) or update "
                "the anchors and surface to PI."
            )
    else:
        K_vol = K_vol_derived

    if V0 is None:
        V0 = (4.0 / 3.0) * math.pi * R_cell**3
    _require_finite_positive("V0", V0)
    _require_finite_positive("K_vol", K_vol)

    p = ResolvedEnclosedVolume(
        K_vol=K_vol,
        V0=V0,
        R_cell=R_cell,
        dP_ref=dP_ref,
        strain_ref=strain_ref,
        turgor_dP0=turgor_dP0,
    )
    # Diagnostic: ΔP at strain_ref must equal dP_ref by construction.
    p.dP_at_strain_ref = K_vol * (3.0 * strain_ref)
    return p


# ---------------------------------------------------------------------------
# Closed-form oracles (TEST-ONLY helpers; NOT used by the runtime force)
# ---------------------------------------------------------------------------
def young_laplace_pressure(gamma: float, R: float) -> float:
    """Young-Laplace pressure for a sphere: ΔP = 2γ/R [Pa].

    ACCEPTANCE ORACLE (KU-3.1 mechanism). Used in tests to check the
    runtime per-bead force recovers the correct equilibrium ΔP. NEVER
    called by the runtime force compute.
    """
    return 2.0 * gamma / R


def equilibrium_pressure_from_volume(p: ResolvedEnclosedVolume, V: float) -> float:
    """Runtime bulk law ΔP = −K_vol·(V−V0)/V0 [Pa] for a given V.

    This is the SAME law the per-bead force uses (exposed as a scalar so
    tests can check the Young-Laplace balance analytically without a sim).
    """
    return -p.K_vol * (V - p.V0) / p.V0


def sphere_volume_for_pressure(p: ResolvedEnclosedVolume, dP: float) -> float:
    """Volume V at which the bulk law yields pressure ``dP``.

    Inverse of :func:`equilibrium_pressure_from_volume`:
    ``V = V0 · (1 − dP/K_vol)``. Test helper.
    """
    return p.V0 * (1.0 - dP / p.K_vol)


# ---------------------------------------------------------------------------
# Custom force compute
# ---------------------------------------------------------------------------
class EnclosedVolumePressure(md.force.Custom):
    """Per-bead intracellular-pressure force on the cortex shell (KU-3.1).

    Each step:

    1. Read shell-bead positions (tags in ``shell_tag_range``).
    2. Centroid ``r_c`` of the shell beads.
    3. ``R_mean = ⟨|r_i − r_c|⟩``; enclosed volume ``V = 4/3 π R_mean³``;
       shell area ``S = 4 π R_mean²``.
    4. ``ΔP = −K_vol (V − V0)/V0``.
    5. ``F_i = ΔP · (S / N) · n̂_i`` with ``n̂_i`` the outward radial unit
       vector ``(r_i − r_c)/|r_i − r_c|``.

    A non-negative diagnostic potential ``½ (K_vol/V0)(V−V0)²`` is reported
    (summed onto the shell beads) so HOOMD's ``.energy`` is meaningful.

    Parameters
    ----------
    p : ResolvedEnclosedVolume
        Resolved enclosed-volume config.
    shell_tag_range : tuple[int, int]
        Tag range [start, end) of the cortex-actin shell beads. Particles
        whose tag is in this range receive the pressure force; all others
        are untouched.

    Notes
    -----
    Implemented as ``md.force.Custom`` (mirrors ``cortex/erm.py``) so the
    contribution participates in HOOMD's ``net_force`` accumulator, which
    the L-M BAOAB Updater reads each step. Computed in float64 from
    positions via ``cpu_local_snapshot`` / ``cpu_local_force_arrays`` —
    row-indexed, same single-rank guarantee erm.py documents.

    The centroid-relative outward normals make the net pressure force
    ≈ 0 (Sanity Gate §3): unlike the ERM external field, this is an
    INTERNAL enclosed-volume potential and must not inject net momentum.
    """

    def __init__(
        self,
        p: ResolvedEnclosedVolume,
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
        # Last-computed diagnostics (populated each set_forces; handy for
        # tests / monitors). Defaults make the pre-run state explicit.
        self.last_volume: float = float("nan")
        self.last_pressure: float = float("nan")
        self.last_R_mean: float = float("nan")

    # -- estimator (static so tests can call it on a bare position array) --
    @staticmethod
    def estimate_volume(positions: np.ndarray) -> tuple[float, float, np.ndarray, np.ndarray]:
        """Sphere-equivalent enclosed volume from a shell-bead array.

        Parameters
        ----------
        positions : ndarray, shape (n, 3)
            Shell-bead positions [m].

        Returns
        -------
        V : float
            Enclosed volume ``4/3 π R_mean³`` [m³].
        S : float
            Shell area ``4 π R_mean²`` [m²].
        centroid : ndarray, shape (3,)
            Bead centroid [m].
        radii : ndarray, shape (n,)
            Per-bead distance from the centroid [m].
        """
        pos = np.asarray(positions, dtype=np.float64).reshape(-1, 3)
        n = pos.shape[0]
        if n == 0:
            return 0.0, 0.0, np.zeros(3, dtype=np.float64), np.empty(0)
        centroid = pos.mean(axis=0)
        radii = np.linalg.norm(pos - centroid, axis=1)
        R_mean = float(radii.mean())
        V = (4.0 / 3.0) * math.pi * R_mean**3
        S = 4.0 * math.pi * R_mean**2
        return V, S, centroid, radii

    def effective_bead_stiffness(self, n_shell: int) -> float:
        """Per-bead radial stiffness of the pressure term [N/m] (CFL gate).

        Moving ONE shell bead radially outward by δr increases R_mean by
        δr/N, hence V by (dV/dR_mean)·(δr/N) = S·(δr/N). The restoring
        change in that bead's outward force per unit δr is

            k_eff = |dF_bead/dr| = A_i · |d(ΔP)/dr|
                  = (S/N) · (K_vol/V0) · (dV/dr)
                  = (S/N) · (K_vol/V0) · S·(1/N)
                  = K_vol · S² / (N² · V0).

        (S = 4πR², V0 = 4/3 πR³ → S²/V0 = 12πR, so
         k_eff = 12 π R K_vol / N² — tiny for N ~ 7000.)
        """
        S = 4.0 * math.pi * self.p.R_cell**2
        return self.p.K_vol * S * S / (float(n_shell) ** 2 * self.p.V0)

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
            V, S, centroid, radii = self.estimate_volume(shell_pos)
            # Baseline osmotic turgor Π₀ (resting intracellular pressure) PLUS
            # the elastic bulk response. At V=V0 the elastic term is 0 so
            # ΔP = Π₀ > 0 → outward pre-tension (the cortex carries a
            # Young-Laplace hoop tension γ ≈ Π₀·R/2 without myosin). Π₀=0
            # recovers the prior force-free-at-construction behaviour.
            dP = self.p.turgor_dP0 - self.p.K_vol * (V - self.p.V0) / self.p.V0
            # Per-bead outward normal about the centroid.
            dx = shell_pos - centroid
            r_safe = np.where(radii > 0.0, radii, 1.0)
            n_hat = dx / r_safe[:, None]
            # Zero the (singular) centroid-coincident beads' direction.
            n_hat[radii <= 0.0] = 0.0
            A_i = S / n_shell
            F_shell = (dP * A_i) * n_hat
            F_vec[mask] = F_shell
            # Diagnostic non-negative bulk potential, distributed evenly.
            U_total = 0.5 * (self.p.K_vol / self.p.V0) * (V - self.p.V0) ** 2
            U_per[mask] = U_total / n_shell

            self.last_volume = V
            self.last_pressure = dP
            self.last_R_mean = float(radii.mean())

        with self.cpu_local_force_arrays as arrays:
            arrays.force[:] = F_vec
            arrays.potential_energy[:] = U_per


# ---------------------------------------------------------------------------
# Public helper: attach enclosed-volume pressure to an existing simulation
# ---------------------------------------------------------------------------
def attach_enclosed_volume_to_simulation(
    sim: hoomd.Simulation,
    p_ev: ResolvedEnclosedVolume,
    *,
    shell_tag_range: tuple[int, int],
    n_shell: int | None = None,
    gamma_b: float | None = None,
    cfl_safety_factor: float = 0.1,
    cfl_strict: bool = True,
) -> EnclosedVolumePressure:
    """Append an EnclosedVolumePressure custom force to an existing Integrator.

    Parameters
    ----------
    sim : hoomd.Simulation
        Already-built cortex simulation (must have an Integrator).
    p_ev : ResolvedEnclosedVolume
    shell_tag_range : tuple[int, int]
        Tag range of the cortex-actin shell beads.
    n_shell : int, optional
        Number of shell beads (= ``shell_tag_range[1] − shell_tag_range[0]``
        if omitted). Used for the CFL gate's per-bead stiffness.
    gamma_b : float, optional
        Per-bead Stokes drag [N·s/m]. Required to gate the pressure CFL
        ``dt ≤ cfl_safety_factor · γ_b / k_eff``. If None, the CFL gate is
        skipped (caller responsible). The enclosed-volume k_eff is
        extremely soft (≪ bond/angle/xl stiffness), so this gate essentially
        never fires — but it is enforced for parity with erm.py.
    cfl_safety_factor : float, default 0.1
        Same convention as the bond/angle/ERM CFL gates (D3 BAOAB).
    cfl_strict : bool, default True
        If True, raise on CFL violation; set False for diagnostic runs.

    Returns
    -------
    ev_force : EnclosedVolumePressure
        The attached force compute (kept for introspection in tests).
    """
    ig = sim.operations.integrator
    if ig is None:
        raise RuntimeError(
            "sim.operations.integrator must be set before attaching "
            "the enclosed-volume pressure force."
        )
    if n_shell is None:
        n_shell = int(shell_tag_range[1]) - int(shell_tag_range[0])

    ev = EnclosedVolumePressure(p_ev, shell_tag_range)

    if gamma_b is not None and n_shell > 0:
        dt = float(ig.dt)
        k_eff = ev.effective_bead_stiffness(n_shell)
        if k_eff > 0.0:
            tau_vol = gamma_b / k_eff
            dt_cfl_vol = cfl_safety_factor * tau_vol
            if dt > dt_cfl_vol and cfl_strict:
                raise RuntimeError(
                    f"Enclosed-volume CFL violated: dt = {dt:.3e} s > "
                    f"{cfl_safety_factor:.2f} · τ_vol = {dt_cfl_vol:.3e} s "
                    f"(τ_vol = γ_b / k_eff = {tau_vol:.3e} s, k_eff = "
                    f"{k_eff:.3e} N/m). Reduce K_vol (softer pressure — "
                    "requires PI sign-off vs the KU-3.1 osmotic anchor) OR "
                    f"reduce dt (need dt ≤ {dt_cfl_vol:.3e} s). Pass "
                    "cfl_strict=False to skip for diagnostic runs."
                )

    ig.forces.append(ev)
    return ev
