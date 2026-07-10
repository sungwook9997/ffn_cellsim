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
    #  • n_fa 100: physiological FA count (KB-2.12 per-cell 10-100 nN / single FA 1-10 nN → ~10-100 FAs). n_fa=50
    #    COLLAPSED traction to 0 on sparse physiological collagen (too few grips = the biphasic soft-arm); n_fa=100
    #    sits in the physiological + working-traction regime. NOT tuned-to-speed — a KB-2.12 count that avoids the
    #    traction-collapse failure mode. The dynamics come from polarize's front-form/rear-release clutch turnover.
    "mesenchymal": CellTypeProfile("mesenchymal", polarize=True, ecm_regrip=True,
                                   front_frac=0.6, myo_rear_bias=0.0, n_fa=100, contractility_mult=1.0),
}
_PROFILES["emt"] = replace(_PROFILES["mesenchymal"], name="emt")
_ALIASES = {"mcf7": "mcf7_epithelial", "epithelial": "mcf7_epithelial", "mda": "mesenchymal",
            "mesench": "mesenchymal", "motile": "mesenchymal"}


def resolve_cell_type(name: str) -> CellTypeProfile:
    """Return the :class:`CellTypeProfile` for ``name`` (registry key or alias); raise on unknown."""
    key = _ALIASES.get(name, name)
    if key not in _PROFILES:
        raise ValueError(f"unknown --cell-type {name!r}; choose from {sorted(_PROFILES)} (+ aliases {sorted(_ALIASES)})")
    return _PROFILES[key]


__all__ = ["CellTypeProfile", "resolve_cell_type"]
