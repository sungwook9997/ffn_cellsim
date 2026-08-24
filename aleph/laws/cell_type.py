"""FF cell-type presets — a named, KB-grounded bundle of the existing migration knobs (NOT new physics).

The epithelial MCF7 and a motile mesenchymal/EMT cell run the SAME fine-grained engine (directed front
polymerization, Stam-Hocky minifilament contraction, catch-slip clutches, collagen re-grip); they differ only in
how those molecular processes are SPATIALLY ORGANISED. A ``CellTypeProfile`` records that organisation as a set
of already-existing knobs, so ``--cell-type`` is pure composition — no lumped migration law.

Grounding (each knob → a verified KB claim):
- **KB-3.11** (migration modes): mesenchymal = strong FA/traction, MMP-dependent, stiff ECM, SINGLE leading edge;
  MCF7 = "mesenchymal+collective hybrid" but poorly-migratory as a single cell.
- **SE248 Betorz2023**: motile cells set a front (high Rac1/Cdc42 → F-actin protrusion) and a rear (LESS actin,
  MORE ACTIVE MYOSIN motors); the cell migrates toward the front, retracting the rear. → ``polarize`` (front
  nascent adhesion / rear de-adhesion) + ``myo_rear_bias`` (rear-weighted minifilaments).
- **KB-4.12** (EMT): E-cadherin 50-100%↓, N-cadherin/vimentin↑ → collective→single-cell switch (grounds the
  single-cell mesenchymal arm; no cell-cell junction here).
- **KB-2.2** (FA maturation) + **KB-2.12** (per-clutch load): mesenchymal keeps the fast nascent-adhesion
  turnover cycle (fewer, front-biased, higher-load FAs) rather than maturing to stable long-lived FBs → ``n_fa``
  capped to a discrete count (KB-2.2 ~5-10 mature + nascent; per-clutch lands in the KB-2.12 5-20 pN band).
- **KB-1.23** (FA captures fibres, R_FA≈1.5µm): the ``ecm_regrip`` walking-adhesion cycle on collagen.

⚠️ PI-flagged (qualitative, not a single measured MCF7-vs-mesenchymal number): ``front_frac`` 0.6,
``myo_rear_bias`` 0.7, ``n_fa`` 50 are grounded in DIRECTION (SE248/KB-3.11/KB-2.2) but the exact values are
modeling choices — they must NOT be swept to hit a target speed. ``contractility_mult`` stays 1.0 (KB-3.14
licenses 2-10× EMT myosin upregulation but that is a separate PI-gated lever, off by default).
"""
from __future__ import annotations

from dataclasses import dataclass, replace


@dataclass(frozen=True, slots=True)
class CellTypeProfile:
    """A named bundle of migration knobs. ``mcf7_epithelial`` == the current default behaviour (back-compat)."""

    name: str
    polarize: bool           # front-Rac nascent adhesion + rear-Rho de-adhesion (--polarize)
    ecm_regrip: bool         # walking adhesion: front clutches re-grip the nearest current collagen (KB-1.23)
    front_frac: float        # leading-edge cap fraction: 0.5 broad (epithelial) vs 0.6 single dominant front
    myo_rear_bias: float     # fraction of minifilaments relocated to the REAR cap (SE248 rear-myosin); 0 = uniform
    n_fa: int                # discrete FA sites (0 = one clutch/node, epithelial); ~50 fewer/dynamic (mesenchymal)
    contractility_mult: float  # EMT myosin upregulation multiplier (PI-gated; 1.0 default)
    # --- MECHANICAL setpoints (None → engine MCF7 default). Unlike the knobs above these change the cell's
    #     PHYSICS (drag, geometry, nucleus), so a preset that sets them is a real cell-type mechanical build, not
    #     MCF7 mechanics migrating in a mesenchymal pattern. Derived/sourced values only — never swept to a target.
    R_um: float | None = None          # suspended cell radius [µm]  (None → CortexParams 7.5, MCF7 Wagner2011)
    eta_Pa_s: float | None = None      # cytoplasm viscosity [Pa·s]  (None → units.ETA_CYTOPLASM 65.9, MCF7 Dessard2024)
    E_nuc_Pa: float | None = None      # nucleus Young's modulus [Pa] (None → resolve_nucleus 399, MCF7 Fischer2020)


_PROFILES = {
    # DEFAULT — byte-for-byte the current run when no other flags are given (the rear-bias hook is a no-op at 0.0)
    "mcf7_epithelial": CellTypeProfile("mcf7_epithelial", polarize=False, ecm_regrip=False,
                                       front_frac=0.5, myo_rear_bias=0.0, n_fa=0, contractility_mult=1.0),
    # MOTILE — front-Rac protrusion / walking adhesions / single leading edge (SE248, KB-3.11, KB-4.12).
    # Native-validated stable config (2026-07-11 autonomous decisions, plan autonomous-log):
    #  • myo_rear_bias 0.0: the 70%-rear-concentrated minifilament relocation is NUMERICALLY UNSTABLE at native
    #    (OverflowError in the from-resting relaxation) AND added no migration benefit — dropped from the preset
    #    (the --myo-rear-bias flag + build hook remain for PI experimentation). Contraction = the existing uniform
    #    myosin; the migration front-back ASYMMETRY is set by polarize (adhesion turnover), the validated mechanism.
    #  • n_fa 100: the only STABLE+running native config. n_fa=0 (one/node, 10956 clutches) + ecm_regrip + com_drag
    #    OVERFLOWS (numerical, the many-clutch re-grip on compliant collagen); n_fa=50 COLLAPSES traction. n_fa=100
    #    runs stably (v≈0.22 nm/s). HONEST LIMITATION (2026-07-11, not a tuning target): native migration of a
    #    strongly-adherent cell on physiological collagen is sub-physiological (v~0.2-0.57) and numerically fragile
    #    across these mechanisms — the model reproduces MCF7's poor motility; full-speed invasion is an open
    #    adhesion-drag/solver item, NOT reachable by parameter choice. MMP needs more front traction than this config
    #    sustains, so proteolysis is under-engaged here (flagged). The migration is directionally correct (biphasic
    #    in adhesion/matrix density, compliant helps) but modest — reported honestly rather than tuned to a target.
    "mesenchymal": CellTypeProfile("mesenchymal", polarize=True, ecm_regrip=True,
                                   front_frac=0.6, myo_rear_bias=0.0, n_fa=100, contractility_mult=1.0),
}
_PROFILES["emt"] = replace(_PROFILES["mesenchymal"], name="emt")
# MDA-MB-231 — the mesenchymal ORGANISATION (identical migration knobs) PLUS 3 of the cell's ~7 mechanical axes
# set to 231-specific values (R, cytoplasm η, nucleus E_nuc). This is a PARTIAL 231 mechanical build: turgor ΔP
# (40 Pa MCF7 proxy), the ΔP·R/2-emergent cortical tension γ, cortical/actin density, contractility (NMII) and
# adhesion/cadherin identity all remain MCF7/generic and PI-gated (unresolved) — see
# docs/v2_audit/MDA_MB_231_BASELINE_CONTRACT_2026-07-13.md §3C. NOT yet a complete 231 cell. Mechanics wired
# (2026-07-13, primary-source verified in wf_35f1b9e4-d78):
#   • eta 10.7 Pa·s — KB-PIV-7 (Dessard/Manneville/Berret 2024, 10.1039/d4na00003j) Table-1 primary value
#     (L=3±1µm; the all-wires pooled figure is 12.0±5.7). PI pick 2026-07-13: paper's quotable primary. ~6× < MCF7 65.9.
#   • R 7.5 µm — direct suspended single-cell radius (Cognart/Viovy/Villard 2020, 10.1038/s41598-020-63316-w, dia ~15µm;
#     PI pick 2026-07-13 over the earlier 8.0). CORRECTS the vol-derived 10.5µm (RodriguezCruz 3D volume was
#     protrusion-inflated) → 231 suspended size ≈ MCF7 7.5µm (231 is a SMALL mesenchymal cell). band 7.5-10.
#   • E_nuc 157.7 Pa — direct 231 in-situ over-nucleus modulus (Fischer, Hayn & Mierke 2020, 10.3389/fcell.2020.00393
#     — SAME paper as MCF7's 399; 231 nucleus ~2.5× softer, lamin-A/C↓; corroborated by Poertner2025 live-MCF7-nucleus
#     150-170 Pa). resolve_nucleus lower band extended 200→150 Pa (PI-ratified 2026-07-13) to admit it. ratio_lamin
#     stays 3.0 (no quantitative 231 lamin-A/C fold-change exists → not invented). Low-impact for basal-substrate crawl.
# ⚠️ Inputs for a DIAGNOSTIC run; R/eta/E_nuc all sourced+ratified. Never swept to a target.
_PROFILES["mda_mb_231"] = replace(_PROFILES["mesenchymal"], name="mda_mb_231",
                                  R_um=7.5, eta_Pa_s=10.7, E_nuc_Pa=157.7)
_ALIASES = {"mcf7": "mcf7_epithelial", "epithelial": "mcf7_epithelial", "mda": "mda_mb_231",
            "mda231": "mda_mb_231", "231": "mda_mb_231", "mesench": "mesenchymal", "motile": "mesenchymal"}


def resolve_cell_type(name: str) -> CellTypeProfile:
    """Return the :class:`CellTypeProfile` for ``name`` (registry key or alias); raise on unknown."""
    key = _ALIASES.get(name, name)
    if key not in _PROFILES:
        raise ValueError(f"unknown --cell-type {name!r}; choose from {sorted(_PROFILES)} (+ aliases {sorted(_ALIASES)})")
    return _PROFILES[key]


__all__ = ["CellTypeProfile", "resolve_cell_type"]
