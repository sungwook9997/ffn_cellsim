"""A1: ligand identity → ACTIVE edge-traction magnitude (the faithful Bare/Pre/Lam4 driver).

L2.6 mapped the three PI ligand conditions (Bare / Pre / Lam4) onto the COARSE *passive*
substrate-adhesion depth (``substrate.resolve_substrate(adhesion_ratio=…)``) and found the
emergent A/A₀(R₀) curves barely separate — a cohesive MCF7 forms a 3D cap regardless of the
passive well depth (L2.6 REPORT). The L2.6 finding *and* the validated traction axis (A/A₀
rises monotonically with active edge-traction) point to the real driver: the ligand modulates
the cell's **active traction** through its integrin-clutch kinetics, not the passive adhesion.

This module is that faithful upgrade. It maps each PI condition → an edge-cell traction
magnitude (the ``f_traction`` of ``spheroid.spreading.edge_outward_forces`` /
``proliferation.run_growth_pooled``) derived from the per-species integrin catch-slip kinetics
in ``bridge.ligand_species`` (the literature SoT), so the three conditions yield three distinct
emergent a/b/c curves *from the mechanism*, never fit to the PI poster.

Mechanism (engaged-clutch → traction; PI-exp map §"b = bound-integrin count × clutch force")
--------------------------------------------------------------------------------------------
A cell-substrate edge contact transmits traction through an ensemble of integrin clutches:

    f_traction(condition) = T_ref · density(condition) · clutch_strength(ligand)

* **clutch_strength(ligand)** — the ANCHORED ligand-IDENTITY factor (col-I vs laminin-111),
  derived entirely from measured kinetics. Per engaged clutch the steady transmitted force
  scales as the bound occupancy × the bond's characteristic force:

      strength(ligand) ∝ φ(ligand) · F_s(ligand)
      φ      = k_on / (k_on + k_off0)            steady bound fraction (k_on = KU-2.4 0.3 s⁻¹)
      F_s    = k_BT / x_β                         Bell-Evans slip characteristic force

  reported RELATIVE to the col-I reference (=1.0). With the registry anchors (col-I k_off0
  1.3 s⁻¹ / x_β 0.23 nm; laminin-111 k_off0 1.85 s⁻¹ / x_β 0.28 nm, α7β1-invasin slip proxy):
  col-I φ·F_s = 0.1875·18.6 pN, laminin = 0.1395·15.3 pN ⇒ **laminin ≈ 0.61× col-I** — i.e.
  laminin is the WEAKER clutch (lower occupancy AND lower per-bond force), exactly the measured
  direction (breast-epithelial traction lower on LN-111, P=0.016–0.028; Taubenberger 2007 col-I
  vs the laminin proxy). The ordering is NOT guessed — it falls out of the registry kinetics.

* **density(condition)** — the ligand-AVAILABILITY factor (Bare vs Pre, both col-I). Bare =
  untreated pV4D4 (low adsorbed density → diffuse β1, low n_engaged); Pre = surface-adsorbed
  col-I (higher density → peripheral β1); Lam4 = + soluble laminin-111 (uniform β1).
  ⚠️ The pV4D4 col-I adsorption DENSITY is a documented literature GAP (PI-exp map: "no
  surviving primary quantitative claim … route to the collaborator's lab" — Im Sung Gap, KAIST,
  who made the surface). So the Bare<Pre ratio is a FLAGGED modeling axis (qualitative anchor:
  Bare low-availability < Pre higher-density; LIGAND_PRESENTATION_MECHANISM.md), a sweepable
  argument, NOT a measured constant — pending the collaborator datum. The ligand-IDENTITY axis
  above is the anchored result of A1; the density axis is the flagged knob.

* **T_ref** — the absolute traction scale of the reference (col-I, reference density) clutch.
  There is NO direct MCF7 single-cell traction in the literature (the 15–25 nN micropillar
  value was REFUTED in verification), so the absolute scale is unanchored; it is set inside the
  B1 stable band (edge traction ≤ ~3 nN/cell or the overdamped CBM detaches the cell) and the
  SCIENCE is in the relative ordering + the emergent curve shapes, not the absolute magnitude.

Single-cell ↔ collective note (honest, expected): this maps the SINGLE-CELL traction, where
laminin is the weaker clutch. The PI poster's *collective* Lam4 enhancement (small-size penalty
removed) is a documented single-cell↔collective split (PI-exp map). Whether the emergent
*collective* A/A₀ reproduces the poster ordering is the A2 overlay question — A1 must NOT
engineer it; it reports what the measured single-cell clutch kinetics produce.

Sanity Gate
-----------
- Dimensional: k_off0 [s⁻¹]; x_β [m]; F_s,f_traction [N]; φ,density,strength [–].
- Boundary: density→0 ⇒ f_traction→0 (no ligand, no active traction); col-I strength ≡ 1.0.
- Sign/sense: lower k_off0 ⇒ higher φ ⇒ stronger clutch ⇒ more traction; laminin (higher
  k_off0, smaller F_s) ⇒ weaker than col-I. Higher density ⇒ more engaged ⇒ more traction.
- Range: every resolved f_traction ≤ the B1 stable ceiling (asserted; surfaced if exceeded).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from ffn_sim.bridge.ligand_species import (
    DEFAULT_LIGAND_FOR_CONDITION,
    K_BT,
    LIGAND_REGISTRY,
)

if TYPE_CHECKING:
    from ffn_sim.spheroid.params import ResolvedL2

__all__ = [
    "ResolvedLigandTraction",
    "KU_2_4_K_ON",
    "T_REF_DEFAULT",
    "DENSITY_FACTORS_DEFAULT",
    "STABLE_TRACTION_CEILING",
    "LP_EDGE",
    "LP_UNIFORM",
    "BETA1_DISTRIBUTION",
    "clutch_strength",
    "resolve_ligand_traction",
    "V0_PROTRUSION_BIELING",
    "V0_WHOLECELL_MCF7",
    "ResolvedActiveTraction",
    "resolve_active_traction",
]

# KU-2.4 generic integrin on-rate (motor-clutch base set); PI-exp validation map §Layer-1.
KU_2_4_K_ON: float = 0.3  # s⁻¹

# Reference edge-traction for the col-I, reference-density clutch (N). UNANCHORED absolute
# scale (no MCF7 single-cell traction in lit; the 15–25 nN micropillar value was REFUTED).
# Set inside the B1 stable band so the strongest condition stays ≤ STABLE_TRACTION_CEILING;
# the RELATIVE ordering, not this magnitude, is the A1 result. A documented numerical policy.
T_REF_DEFAULT: float = 2.5e-9  # N

# B1 stable-band ceiling: edge traction above this detaches an edge cell in the overdamped
# large-dt CBM (REPORT §B1). Resolved tractions are asserted ≤ this.
STABLE_TRACTION_CEILING: float = 3.0e-9  # N

# Ligand-availability (density) factor per condition. FLAGGED modeling axis (Bare<Pre is the
# qualitative anchor; the pV4D4 col-I adsorption density is a literature gap pending the
# collaborator's surface-chemistry data — see module docstring). Pre = the reference (1.0).
DENSITY_FACTORS_DEFAULT: dict[str, float] = {
    "Bare": 0.6,   # low adsorbed density (diffuse β1)            — FLAGGED (pending anchor)
    "Pre": 1.0,    # higher surface-adsorbed col-I (reference)    — FLAGGED (pending anchor)
    "Lam4": 1.0,   # + soluble laminin-111, uniform β1 (= Pre availability; isolates identity)
}

# --- A4: traction DISTRIBUTION axis (β1 uniformity → traction screening length Lp) ----------
# A1 set the traction MAGNITUDE from clutch kinetics. A4 adds the orthogonal DISTRIBUTION axis
# from the measured β1 IF pattern (LIGAND_PRESENTATION_MECHANISM.md): Bare "diffuse" / Pre
# "peripheral" / Lam4 "uniform". In the active-wetting model the traction is screened to within
# Lp of the cluster edge (`spreading.edge_outward_forces`): a SMALL Lp = edge/peripheral
# engagement (the rim ∝ 1/R drives the law), a LARGE Lp = whole-footprint UNIFORM engagement.
# So β1 distribution maps onto Lp: peripheral/diffuse → edge Lp; uniform → Lp ≫ spheroid.
LP_EDGE: float = 11.0e-6        # m  edge screening length (Pérez-González 2019; A1 default)
# "uniform-β1 sentinel": Lp ≫ any spheroid radius ⇒ exp(−depth/Lp) ≈ 1 ⇒ every basal cell
# transmits traction (uniform engagement), the measured Lam4 phenotype. NOT a fitted length —
# it is the Lp→∞ (uniform) limit expressed as a value far above the cluster scale (~0.1–0.4 mm).
LP_UNIFORM: float = 1.0e-3      # m

# β1 IF distribution per condition (LIGAND_PRESENTATION_MECHANISM.md). "uniform" (Lam4) is the
# A4 collective axis; "diffuse"/"peripheral" (col-I) are edge-localized for traction transmission
# (the Bare diffuse-vs-Pre peripheral sub-distinction is carried by the density factor, A1).
BETA1_DISTRIBUTION: dict[str, str] = {"Bare": "diffuse", "Pre": "peripheral", "Lam4": "uniform"}
_LP_FOR_DISTRIBUTION: dict[str, float] = {
    "diffuse": LP_EDGE, "peripheral": LP_EDGE, "uniform": LP_UNIFORM,
}


@dataclass(frozen=True)
class ResolvedLigandTraction:
    """Resolved ligand-condition → active edge-traction (SI), with its anchored provenance."""

    condition: str          # "Bare" | "Pre" | "Lam4"
    ligand: str             # registry key (collagen_I | laminin_111)
    integrin: str           # receptor (alpha2beta1 | alpha6beta1)
    engaged_fraction: float # φ = k_on/(k_on+k_off0)                          [–]
    F_s: float              # Bell-Evans slip characteristic force k_BT/x_β   [N]
    clutch_strength: float  # φ·F_s relative to the col-I reference (=1.0)    [–]
    density_factor: float   # ligand availability (FLAGGED modeling axis)     [–]
    f_traction: float       # edge-cell traction = T_ref·density·strength     [N]
    beta1_distribution: str # measured β1 IF pattern: diffuse|peripheral|uniform
    Lp: float               # traction screening length (m): edge vs uniform (A4 axis)
    uniform_beta1: bool     # True for Lam4 ("uniform β1" → whole-footprint traction, A4)
    proxy: bool             # True if the ligand kinetics are a literature proxy (laminin)


def _phi_Fs(ligand: str, k_on: float) -> tuple[float, float]:
    """(engaged fraction φ, slip force F_s) for a slip ligand from the registry kinetics."""
    lig = LIGAND_REGISTRY[ligand]
    if lig.bond_class != "slip":
        # A1 maps the PI conditions (col-I, laminin) which are slip; catch-slip (FN) would
        # need the Pereverzev occupancy, not the Bell φ — out of scope for these conditions.
        raise ValueError(
            f"resolve_ligand_traction: {ligand!r} is {lig.bond_class!r}; A1 covers slip "
            "ligands (col-I, laminin). FN catch-slip is not a PI condition."
        )
    assert lig.k_off0 is not None and lig.x_beta is not None
    phi = k_on / (k_on + lig.k_off0)
    F_s = K_BT / lig.x_beta
    return phi, F_s


def clutch_strength(ligand: str, *, k_on: float = KU_2_4_K_ON, reference: str = "collagen_I") -> float:
    """Relative clutch strength φ·F_s of ``ligand`` vs the ``reference`` ligand (=1.0).

    The ANCHORED ligand-identity factor: occupancy × per-bond force, both from measured
    kinetics. < 1 means a weaker clutch (less transmitted traction) than the reference.
    """
    phi, F_s = _phi_Fs(ligand, k_on)
    phi_r, F_s_r = _phi_Fs(reference, k_on)
    return (phi * F_s) / (phi_r * F_s_r)


def resolve_ligand_traction(
    condition: str,
    *,
    t_ref: float = T_REF_DEFAULT,
    k_on: float = KU_2_4_K_ON,
    density_factors: dict[str, float] | None = None,
    uniform_beta1: bool = False,
) -> ResolvedLigandTraction:
    """Resolve a PI condition (Bare/Pre/Lam4) → active edge-traction from clutch kinetics.

    Args:
        condition: one of ``"Bare"``, ``"Pre"``, ``"Lam4"`` (DEFAULT_LIGAND_FOR_CONDITION).
        t_ref: reference edge-traction (N) of the col-I reference-density clutch (absolute
            scale; unanchored — set in the B1 stable band, see module docstring).
        k_on: generic integrin on-rate (KU-2.4 default 0.3 s⁻¹).
        density_factors: per-condition availability factors (FLAGGED modeling axis; defaults
            to ``DENSITY_FACTORS_DEFAULT``). Pass a custom dict to sweep the density axis (A3).
        uniform_beta1: A4 axis. ``False`` (default) → the traction screening length ``Lp`` is
            taken from the condition's measured β1 IF pattern (``BETA1_DISTRIBUTION``: Lam4
            "uniform" → ``LP_UNIFORM``, col-I edge → ``LP_EDGE``). ``True`` forces the uniform
            (whole-footprint) limit regardless of condition — used to A/B the uniform-β1
            mechanism against the A1 edge-localized baseline.

    Returns:
        ``ResolvedLigandTraction`` with the resolved ``f_traction`` and its provenance.

    Raises:
        ValueError: on an unknown condition, a non-finite/negative t_ref, or a resolved
            traction above the B1 stable ceiling (surfaced, never silently clamped).
    """
    if condition not in DEFAULT_LIGAND_FOR_CONDITION:
        raise ValueError(
            f"unknown condition {condition!r}; expected one of "
            f"{sorted(DEFAULT_LIGAND_FOR_CONDITION)}."
        )
    if not (t_ref > 0.0):
        raise ValueError("t_ref must be > 0.")
    dens = DENSITY_FACTORS_DEFAULT if density_factors is None else density_factors
    if condition not in dens:
        raise ValueError(f"density_factors missing condition {condition!r}.")

    ligand = DEFAULT_LIGAND_FOR_CONDITION[condition]
    lig = LIGAND_REGISTRY[ligand]
    phi, F_s = _phi_Fs(ligand, k_on)
    strength = clutch_strength(ligand, k_on=k_on)
    density = float(dens[condition])
    f_traction = float(t_ref * density * strength)

    distribution = BETA1_DISTRIBUTION[condition]
    is_uniform = bool(uniform_beta1 or distribution == "uniform")
    Lp = LP_UNIFORM if is_uniform else _LP_FOR_DISTRIBUTION[distribution]

    if f_traction > STABLE_TRACTION_CEILING:
        raise ValueError(
            f"resolve_ligand_traction({condition!r}) → {f_traction*1e9:.2f} nN exceeds the "
            f"B1 stable ceiling {STABLE_TRACTION_CEILING*1e9:.1f} nN (edge cells would detach). "
            "Lower t_ref or the density factor (REPORT §B1)."
        )

    return ResolvedLigandTraction(
        condition=condition,
        ligand=ligand,
        integrin=lig.integrin,
        engaged_fraction=float(phi),
        F_s=float(F_s),
        clutch_strength=float(strength),
        density_factor=density,
        f_traction=f_traction,
        beta1_distribution=distribution,
        Lp=float(Lp),
        uniform_beta1=is_uniform,
        proxy=bool(lig.proxy),
    )


# ============================================================================================
# Lamellipodium → CBM active-traction SCALE-BRIDGE (2026-06-04, PI direction)
# ============================================================================================
# The A1 `T_REF_DEFAULT` (2.5 nN) is an UNANCHORED heuristic placed inside the B1 stable band
# so the RELATIVE ligand ordering is the science (the magnitude is admittedly a knob). This
# bridge replaces that heuristic magnitude with one ANCHORED to the single-cell H.5
# lamellipodium — the same move that anchored cohesion to the Iturri-2020 de-adhesion force.
#
# Overdamped map (the FROZEN integrator's own relation): a per-cell self-propulsion VELOCITY v0
# is sustained by a force f_active = v0 · γ_cell (since the overdamped step is r += (F/γ)·dt, a
# steady drift v0 needs F = v0·γ). γ_cell is the clutch-ensemble migration drag (KU-2.18
# Bangasser 2013; resolve_layer2), NOT water-Stokes — so this is the crawling-relevant drag.
#
# Two anchor sources (the ~6× between them IS the magnitude question — overlay-only, never fit):
#  (A) PROTRUSION anchor — the platform's OWN runtime mechanism (literature-first, in-tree):
#      v0_prot = k_elong0 · δ_elong = 11.6 s⁻¹ · 2.7 nm = 31.3 nm/s = 1.88 µm/min  (Bieling 2016
#      slip force-velocity, phase1_h5.yaml; the unloaded barbed-end translocation the KU-5.2
#      oracle inverts). ⚠️ NOT k_elong0·rest_length (5.8 µm/s) — that is a network overgrowth
#      rate, not a translocation, and is unphysical as a cell speed.
#  (B) WHOLE-CELL anchor — an MCF7-specific single-cell speed used as a VALIDATION OVERLAY (not
#      the platform's mechanism): v0_wc ≈ 0.32 µm/min (≈19 µm/h, breast-epithelial single-cell
#      migration). FLAGGED as overlay/proxy: single-cell ≠ collective, and it is an MCF7 number
#      we compare against, not a mechanism we run.
V0_PROTRUSION_BIELING: float = 11.6 * 2.7e-9      # m/s  k_elong0·δ_elong (Bieling 2016, in-tree)
V0_WHOLECELL_MCF7: float = 0.32e-6 / 60.0          # m/s  ≈0.32 µm/min (overlay proxy, FLAGGED)


@dataclass(frozen=True)
class ResolvedActiveTraction:
    """Lamellipodium-anchored per-cell active traction f_active = v0·γ_cell (SI), both anchors."""

    gamma_cell: float            # N·s/m  clutch-ensemble migration drag (KU-2.18)
    v0_protrusion: float         # m/s    Bieling barbed-end translocation (platform mechanism)
    v0_wholecell: float          # m/s    MCF7 single-cell speed (overlay proxy)
    f_protrusion: float          # N      v0_protrusion·γ_cell (the literature-first anchor)
    f_wholecell: float           # N      v0_wholecell·γ_cell (the overlay-anchored value)
    ceiling: float               # N      B1 stable-band ejection ceiling (STABLE_TRACTION_CEILING)
    protrusion_exceeds_ceiling: bool  # f_protrusion > ceiling → overdamped CBM ejects (D3/B1 gate)


def resolve_active_traction(resolved: "ResolvedL2") -> ResolvedActiveTraction:
    """Resolve the lamellipodium-anchored per-cell active traction (no sim — pure arithmetic).

    Returns BOTH anchors (protrusion = the platform's Bieling mechanism; whole-cell = the MCF7
    overlay proxy) so the magnitude question is explicit, and flags whether the protrusion
    anchor exceeds the B1 ejection ceiling (the regime where the overdamped large-dt CBM cannot
    hold an edge cell — the REPORT B1/D3 sub-stepping question, surfaced not silently capped).

    Args:
        resolved: resolved Layer-2 params (provides ``gamma_cell``, the clutch-ensemble drag).
    """
    g = float(resolved.gamma_cell)
    f_prot = V0_PROTRUSION_BIELING * g
    f_wc = V0_WHOLECELL_MCF7 * g
    return ResolvedActiveTraction(
        gamma_cell=g,
        v0_protrusion=V0_PROTRUSION_BIELING,
        v0_wholecell=V0_WHOLECELL_MCF7,
        f_protrusion=float(f_prot),
        f_wholecell=float(f_wc),
        ceiling=STABLE_TRACTION_CEILING,
        protrusion_exceeds_ceiling=bool(f_prot > STABLE_TRACTION_CEILING),
    )
