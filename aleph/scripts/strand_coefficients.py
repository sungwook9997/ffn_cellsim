#!/usr/bin/env python
r"""Derive the axial and bending coefficients for the strand populations that carry NO law today.

**Why this exists.** Measured 2026-08-24 on the tau record: eleven populations stand and **two carry
force**. `microtubule` sits in the arena with 73,500 segments and 73,000 bending triples and no law
bound to any of them; the same for `filopodium`, `stress_fiber`, `sf_arc`, `lamellipodium`,
`intermediate_filament` and `chromatin`. The driver's `--bind POP=k_axial,kappa` has existed the whole
time and the tau runs recorded `extra_bindings: []` — **nobody passed the flag.** So this is not new
machinery; it is the seven pairs of numbers the existing flag refuses to default.

**PI ruling 2026-08-24, option (다) — MIXED.** Derive where the literature is firm; declare an axis
where it is not. The two are DIFFERENT GRADES and this module refuses to blur them: a derived value
carries the sourced constants it came from, and a gap raises rather than returning a plausible number.

**The convention, and it is already fixed by the kernels.**

* ``kappa_pn_um2`` is the flexural rigidity **κ = k_B·T·L_p** [pN·µm²]. It is a MATERIAL property and
  grid-invariant; ``accumulate_bending`` forms the per-triple stiffness itself as ``α = κ / seg³``.
  So the number passed does not change with the discretisation, which is what the charter's
  "derivable and grid-invariant" requires. ⚠ The step enters CUBED, so it is required rather than
  inferred — a 2× error in the step is an 8× error in the force.
* ``k_axial_pn_per_um`` is **EA / L_seg** [pN/µm], and that one IS grid-dependent by construction: a
  chain of stiffer, shorter springs has the same continuum stiffness. Reported with its step.

⚠ **SIX OF THE SEVEN STRAND POPULATIONS ARE DISCRETISED AT THE CORTEX STEP — AND CHROMATIN IS NOT.**
``populations.build_remaining_populations`` takes one ``seg_um`` and hands it to six builders, its own
docstring calling it *"the cortex constant the provenance audit"* names. ``build_chromatin`` takes
``subunit_um`` instead and stands its chain at **0.6 µm**, the source's Mbp domain — twelve times
coarser. An earlier draft of this paragraph said "every", and the difference is not cosmetic: a step
enters bending CUBED, so twelve times is one thousand seven hundred times.

⚠ And the shared step is itself an unratified modelling decision. 50 nm is a fine step for actin
(L_p 17 µm) and resolves a microtubule (L_p ≈ 4,673 µm, derived below) roughly a hundred thousand
times finer than its own persistence length — resolution the physics cannot use, bought at 73,500
segments. REPORTED here rather than changed, because changing it changes the cell.

⚠ **WHAT THIS IS NOT.** It is not a claim that these coefficients are MCF7's. It is not a claim that
these laws belong on these populations — whether a chromatin subunit chain should carry a worm-like
bending term at all is a modelling question with its own answer. And no magnitude produced by a run
that uses these becomes quotable: `STATE.md` (c) stands.

Sanity Gate:
    * dimensional — κ [pN·µm²] = [pN·µm]·[µm]; EA [pN] = [pN/µm²]·[µm²]; k_axial [pN/µm] = [pN]/[µm].
      Each is asserted against its factors' units in the table below rather than assumed.
    * boundary cases — a population whose material constant is absent RAISES with the name of the
      missing constant; it never falls back to actin, which would silently give a microtubule the
      stiffness of an actin filament.
    * conservation — the two independent routes to a bending rigidity (persistence length, and
      Young's modulus × second moment) are BOTH computed where both inputs exist, and their
      disagreement is reported. They are not averaged.
    * CFL/precision — none; this is closed-form arithmetic on host floats.
    * sign sense — every coefficient must be strictly positive; a non-positive one raises.
    * measurement protocol — pure arithmetic over declared constants. No device, no run, no fit.

engine units: length µm, force pN, energy pN·µm. Runtime: host arithmetic; CPU-importable.

Usage:
    python aleph/scripts/strand_coefficients.py
    python aleph/scripts/strand_coefficients.py --bind-flags
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

from aleph.laws import units as U

#: The discretisation every strand population is built at. ⚠ NOT a choice made here — read from the
#: same axis the builders read, so this file cannot drift from the cell it describes.
from aleph.world.build import NATIVE_AXES

SEG_UM: float = float(NATIVE_AXES["cortex_seg_um"]["value"])   # type: ignore[index]

#: Actin filament axial stiffness [pN]. Kojima 1994, via ``laws.kim_network.EA_ACTIN_PN``.
EA_ACTIN_PN: float = 4.4e4

#: Intermediate-filament elastic constants, from ``laws.intermediate_filaments``.
E_IF_PN_UM2: float = 6.0e6      # Kreplak 2005 / Guo 2013; band [1, 10] MPa
D_IF_UM: float = 0.010          # Mücke 2004, canonical TEM
LP_IF_BAND_UM: tuple[float, float] = (0.3, 1.0)

#: Chromatin intersubunit spring [pN/µm] — Stephens et al. 2017 MBoC 28:1984 Table 1 p7, 1.6 nN/µm.
#: ⚠ Already a SPRING CONSTANT per link, so it is NOT divided by a discretisation step: the source's
#: element is the Mbp domain itself. `build/chromatin` stands 552 of them, one chain, 551 segments.
K_CHROMATIN_PN_PER_UM: float = 1.6e3

#: Chromatin subunit diameter [µm], same table. ⚠ This population is discretised at 0.6 µm while the
#: other six are at the cortex's 0.05 µm — a twelve-fold difference that matters wherever a step
#: enters a coefficient, which for bending is CUBED.
SUBUNIT_UM_CHROMATIN: float = 0.6

#: Microtubule outer radius [µm] — ``laws.microtubule.MT_RADIUS_UM``, 25 nm diameter.
R_MT_OUTER_UM: float = 0.0125

#: Microtubule LUMEN radius [µm]. ⚠ **DECLARED 2026-08-24 by PI ruling, and it is not a free choice.**
#: The repository declared only the outer radius, so EA had no annular area and this module refused it.
#: 8.5 nm is the standard 13-protofilament lumen (wall ≈ 4 nm), and what makes it a DERIVATION rather
#: than a guess is that it is over-determined by a constant already sourced here: with it,
#: ``E = KAPPA_MT / I`` comes out at 1.33 GPa against a literature range of 1-2 GPa. The bending
#: rigidity and the geometry check each other; a wall thickness chosen to make something pass would
#: put E outside that range and be visible. `_demo` asserts the check rather than the value.
R_MT_LUMEN_UM: float = 0.0085

#: Cortex axial stiffness [pN/µm] — PI ruling 2026-08-24. ⚠ **The production runs used 88,000, which
#: is 10x SOFTER than the same actin this module derives for every other actin population.** It was
#: honestly recorded as a TEST POINT, so it is not a hidden defect; but the charter requires starting
#: at the physiological operating point, and leaving it would give ONE MATERIAL TWO STIFFNESSES —
#: stress fibres at 880,000 and the cortex at 88,000 — which makes the label rather than the material
#: set the physics. The cost of changing it is lower than it sounds: every gamma and tau taken on the
#: old value is already non-quotable under `STATE.md` (c) 3 and (c) 17.
K_CORTEX_PN_PER_UM: float = EA_ACTIN_PN / 0.05


@dataclass(frozen=True, slots=True)
class Coefficient:
    """One population's pair, with the grade and the inputs that produced it."""

    population: str
    kappa_pn_um2: float | None
    k_axial_pn_per_um: float | None
    grade: str                    # DERIVED | SOURCED | PI_GAP
    basis: str

    def flag(self) -> str | None:
        """The driver argument for this population, or ``None`` when the axial factor is missing.

        ⚠ Three shapes, and getting them wrong is a run that refuses at launch or one that binds
        nothing:
          * ``cortex`` does NOT go through ``--bind``. ``build_all`` stands it, not
            ``build_remaining_populations``, and `world_phase4_native` refuses a ``--bind`` naming a
            population PHASE 1 did not build. It takes ``--cortex-k-axial`` / ``--cortex-kappa``.
          * a population with no bending triple emits ``POP=k_axial`` alone — supplying a kappa there
            is refused, because nothing would launch over it.
          * everything else emits ``POP=k_axial,kappa``.
        """
        if self.k_axial_pn_per_um is None:
            return None
        if self.population == "cortex":
            if self.kappa_pn_um2 is None:
                return None
            return (f"--cortex-k-axial {self.k_axial_pn_per_um:.6g} "
                    f"--cortex-kappa {self.kappa_pn_um2:.6g}")
        if self.kappa_pn_um2 is None:
            return f"--bind {self.population}={self.k_axial_pn_per_um:.6g}"
        return f"--bind {self.population}={self.k_axial_pn_per_um:.6g},{self.kappa_pn_um2:.6g}"


def kappa_from_persistence_length(lp_um: float) -> float:
    """κ = k_B·T·L_p [pN·µm²]. The definition, not a fit."""
    if not (lp_um > 0.0):
        raise ValueError(f"persistence length must be positive; got {lp_um!r}")
    return float(U.KBT) * float(lp_um)


def kappa_from_rod(young_pn_um2: float, diameter_um: float) -> float:
    """κ = E·I for a solid circular rod, ``I = π d⁴ / 64`` [pN·µm²].

    ⚠ SOLID. A microtubule is a hollow tube and this form does not apply to it; the caller must not
    reach for this as a generic fallback, which is why it takes a diameter rather than a population.
    """
    second_moment = math.pi * diameter_um**4 / 64.0
    return float(young_pn_um2) * second_moment


def ea_axial(young_pn_um2: float, diameter_um: float) -> float:
    """EA for a solid circular rod [pN]."""
    return float(young_pn_um2) * math.pi * (diameter_um / 2.0) ** 2


def derive(seg_um: float = SEG_UM) -> list[Coefficient]:
    """The seven pairs, each with its grade. Raises on a non-positive step."""
    if not (seg_um > 0.0):
        raise ValueError(f"seg_um must be positive; it divides EA and cubes into the bending: {seg_um!r}")

    kappa_actin = kappa_from_persistence_length(U.LP_ACTIN_UM)
    k_actin = EA_ACTIN_PN / seg_um

    # ── IF: over-determined, and the two routes AGREE. Reported as a cross-check, never averaged.
    kappa_if_rod = kappa_from_rod(E_IF_PN_UM2, D_IF_UM)
    lp_lo, lp_hi = (kappa_from_persistence_length(x) for x in LP_IF_BAND_UM)
    if_agrees = lp_lo <= kappa_if_rod <= lp_hi
    ea_if = ea_axial(E_IF_PN_UM2, D_IF_UM)

    # ── MT: κ is SOURCED outright. The axial route needs the lumen radius, which this repo does not
    #    declare, so it refuses. E is BACK-DERIVED from the sourced κ only to report the cross-check.
    kappa_mt = float(U.KAPPA_MT)
    i_mt = math.pi * (R_MT_OUTER_UM**4 - R_MT_LUMEN_UM**4) / 4.0
    e_mt = kappa_mt / i_mt
    a_mt = math.pi * (R_MT_OUTER_UM**2 - R_MT_LUMEN_UM**2)
    ea_mt = e_mt * a_mt

    out = [
        Coefficient("microtubule", kappa_mt, ea_mt / seg_um, "SOURCED + DERIVED",
                    f"κ = KAPPA_MT = {kappa_mt} pN·µm² (Nédélec & Foethke 2007, NJP 9:427 p9) — "
                    f"SOURCED, implying L_p = κ/k_BT = {kappa_mt / U.KBT:,.0f} µm. "
                    f"k_axial DERIVED from the sourced κ plus the tube geometry: I = π(r_o⁴-r_i⁴)/4 = "
                    f"{i_mt:.4g} µm⁴, so E = κ/I = {e_mt / 1e9:.2f} GPa — ⚠ CROSS-CHECK against the "
                    f"literature 1-2 GPa, which is what makes the lumen radius a derivation and not a "
                    f"guess. EA = E·π(r_o²-r_i²) = {ea_mt:.4g} pN. "
                    f"⚠ AND IT BECOMES THE STIFFEST ELEMENT IN THE CELL: the explicit overdamped bound "
                    f"1/(μ·k) falls to {1.0 / (1e-6 * ea_mt / seg_um):.4g} s at μ=1e-6, so dt=0.05 s "
                    f"keeps only {1.0 / (1e-6 * ea_mt / seg_um) / 0.05:.1f}x margin against 227x today. "
                    f"⚠ A SEPARATE QUESTION FOLLOWS: at L_p = {kappa_mt / U.KBT:,.0f} µm a 50 nm step "
                    f"resolves a microtubule {kappa_mt / U.KBT / seg_um:,.0f} times finer than its own "
                    f"persistence length — resolution the physics cannot use, bought at 73,500 segments."),
        Coefficient("intermediate_filament", kappa_if_rod, ea_if / seg_um, "DERIVED",
                    f"κ = E·πd⁴/64 = {kappa_if_rod:.4g} pN·µm² from E={E_IF_PN_UM2:.1e} pN/µm², "
                    f"d={D_IF_UM} µm. ⚠ CROSS-CHECK: the independent persistence-length band "
                    f"L_p ∈ {LP_IF_BAND_UM} µm gives κ ∈ [{lp_lo:.4g}, {lp_hi:.4g}] and the rod value "
                    f"{'FALLS INSIDE IT' if if_agrees else 'DOES NOT — do not use either'}. "
                    f"EA = {ea_if:.4g} pN."),
        Coefficient("stress_fiber", kappa_actin, k_actin, "DERIVED",
                    f"actin: κ = k_BT·L_p = {kappa_actin:.4g}, EA = {EA_ACTIN_PN:.1e} pN (Kojima 1994). "
                    f"⚠ SINGLE-FILAMENT values on a structure that is a BUNDLE — the builder stands one "
                    f"chain per fibre, so the bundling factor is absent here and is a modelling gap."),
        Coefficient("sf_arc", kappa_actin, k_actin, "DERIVED",
                    "same actin constants and the same bundle caveat as stress_fiber."),
        Coefficient("lamellipodium", kappa_actin, k_actin, "DERIVED",
                    "single actin filaments; the dendritic branch angle is a separate law, unbound."),
        Coefficient("filopodium", kappa_actin, k_actin, "DERIVED",
                    "⚠ actin SINGLE-filament values. A filopodium is a fascin-crosslinked bundle of "
                    "~20 and its effective rigidity is super-additive in the crosslink density; cards "
                    "F1-F4 hold exactly those constants and none is bound. Using these makes a "
                    "filopodium as floppy as one filament, which is a stated understatement."),
        Coefficient("chromatin", None, K_CHROMATIN_PN_PER_UM, "SOURCED",
                    f"⚠ I REFUSED THIS FIRST AND THE REFUSAL WAS WRONG. k_axial = kp = "
                    f"{K_CHROMATIN_PN_PER_UM:,.0f} pN/µm is SOURCED — Stephens et al. 2017 MBoC "
                    f"28:1984 Table 1 p7, 1.6 nN/µm, the same table build/lamina reads. It is passed "
                    f"AS THE SPRING, not divided by a step: the source's segment IS the Mbp domain "
                    f"(σp = {SUBUNIT_UM_CHROMATIN} µm), so this population is NOT on the cortex "
                    f"discretisation the other six share. "
                    f"⚠ κ is REFUSED for a structural reason, not a missing number: build/chromatin "
                    f"claims NO ANGLE3 range — the source models chromatin as *extensible springs of "
                    f"ZERO bending modulus* and the builder says so in its own words. There is no "
                    f"bending triple to launch over, and adding a worm-like term would be adding "
                    f"physics the source does not have."),
    ]
    out.append(Coefficient(
        "cortex", kappa_actin, K_CORTEX_PN_PER_UM, "DERIVED",
        f"⚠ ALREADY BOUND, and at 88,000 pN/µm — 10x softer than this same actin. PI ruling "
        f"2026-08-24 raises it to EA/L_seg = {K_CORTEX_PN_PER_UM:,.0f}. Listed here so the cell has "
        f"ONE actin stiffness rather than one per population. ⚠ Every gamma and tau measured on the "
        f"old value describes a cortex 10x softer than actin; all of them are already non-quotable "
        f"under STATE.md (c) 3 and (c) 17, which is why the change is affordable."))

    for c in out:
        for value in (c.kappa_pn_um2, c.k_axial_pn_per_um):
            if value is not None and not (value > 0.0 and math.isfinite(value)):
                raise ValueError(f"{c.population}: non-positive or non-finite coefficient {value!r}")
    return out


def _demo() -> None:
    """Self-check without pytest.

    ⚠ Pins the CROSS-CHECK and the REFUSALS, not the values — a value pinned here would have to be
    edited whenever a sourced constant is corrected, which is the wrong direction of authority.
    """
    rows = derive()
    by = {c.population: c for c in rows}

    # The two independent IF routes must agree, or the derivation is not a derivation.
    kappa_rod = kappa_from_rod(E_IF_PN_UM2, D_IF_UM)
    lo, hi = (kappa_from_persistence_length(x) for x in LP_IF_BAND_UM)
    assert lo <= kappa_rod <= hi, (lo, kappa_rod, hi)

    # Actin's declared kappa must reproduce laws.units' own constant, or the two have drifted.
    assert abs(kappa_from_persistence_length(U.LP_ACTIN_UM) - U.KAPPA_ACTIN) < 1e-12

    # Refusals stay refusals.
    # ⚠ Chromatin's axial is SOURCED and its bending is structurally absent (no ANGLE3 range).
    # An earlier draft refused both and was wrong on both counts.
    assert by["chromatin"].kappa_pn_um2 is None, "chromatin has no bending triple to launch over"
    assert by["chromatin"].k_axial_pn_per_um == K_CHROMATIN_PN_PER_UM
    assert by["chromatin"].grade == "SOURCED"
    # MT axial is DERIVED now, and what licenses it is the cross-check, so pin THAT.
    i_mt = math.pi * (R_MT_OUTER_UM**4 - R_MT_LUMEN_UM**4) / 4.0
    e_gpa = (U.KAPPA_MT / i_mt) / 1e9
    assert 1.0 <= e_gpa <= 2.0, (
        f"the lumen radius implies E = {e_gpa:.2f} GPa, outside the literature 1-2 GPa. The geometry "
        f"and the sourced kappa no longer check each other, so neither may be used.")
    assert by["microtubule"].k_axial_pn_per_um is not None
    assert by["microtubule"].kappa_pn_um2 == U.KAPPA_MT

    # A microtubule must not end up with actin's rigidity.
    assert by["microtubule"].kappa_pn_um2 > 100.0 * by["stress_fiber"].kappa_pn_um2

    n_flags = sum(1 for c in rows if c.flag() is not None)
    print(f"[strand_coefficients] self-check OK — {n_flags} of {len(rows)} populations bindable, "
          f"seg {SEG_UM} µm")


def main(argv: list[str] | None = None) -> int:
    """Print the table, or the flags for the driver."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bind-flags", action="store_true", help="print only the --bind arguments")
    ap.add_argument("--seg-um", type=float, default=SEG_UM)
    args = ap.parse_args(argv)

    rows = derive(args.seg_um)
    if args.bind_flags:
        print(" ".join(c.flag() for c in rows if c.flag()))
        return 0

    print(f"discretisation: seg = {args.seg_um} µm, the SAME step for every strand population "
          f"(populations.build_remaining_populations passes one value to all seven)\n")
    for c in rows:
        kap = "REFUSED" if c.kappa_pn_um2 is None else f"{c.kappa_pn_um2:.6g}"
        axi = "REFUSED" if c.k_axial_pn_per_um is None else f"{c.k_axial_pn_per_um:.6g}"
        print(f"{c.population:24s} κ {kap:>12s} pN·µm²   k_axial {axi:>12s} pN/µm   [{c.grade}]")
    print()
    for c in rows:
        print(f"── {c.population}\n   {c.basis}\n")
    bindable = [c for c in rows if c.flag()]
    print(f"{len(bindable)} of {len(rows)} populations have both coefficients. "
          f"{len(rows) - len(bindable)} REFUSE rather than defaulting.")
    print("⚠ No magnitude from a run using these is quotable. These are physiological starting "
          "points, not MCF7 measurements, and STATE.md (c) stands.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
