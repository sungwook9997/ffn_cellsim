"""ECM ligand identity → integrin clutch kinetics (col-I / laminin-111 / FN).

STANDALONE, IMPORT-ISOLATED data module — imported by **no** existing runtime file,
so attaching it leaves the live build and the existing test suite bit-for-bit
unaffected. It is the literature-anchored source-of-truth that the (Lead-owned)
``bridge/fa.py`` ligand-identity extension consumes; the wiring is a Lead sequencing
step documented in ``docs/LIGAND_IDENTITY_FA_SPEC.md``. NOT wired here.

The PI experiment (MCF7 spheroid on iCVD pV4D4, col-I coated; PI-exp validation map)
presents three ligand conditions; the platform clutch (``bridge/integrin_bonds.py``)
currently models a single FN-like α5β1 catch-slip bond (KU-2.18). This module gives
the clutch a **ligand identity** as a per-species parameter set.

Design — ligand identity is a ``PereverzevParams`` swap, not a new mechanism
----------------------------------------------------------------------------
The runtime off-rate is ``pereverzev_k_off(F, PereverzevParams(k_s, F_s, k_c, F_c))``
with ``k_off(F) = k_s·e^{F/F_s} + k_c·e^{−F/F_c}`` (slip term + catch term).

* **col-I → α2β1** and **laminin-111 → α6β1** are **slip** bonds. A Bell-Evans slip
  ``k_off(F) = k_off0·e^{F·x_β/kT}`` maps onto the Pereverzev form with
  ``k_s = k_off0`` and ``F_s = k_BT/x_β``. The catch pathway is *disabled* by setting
  ``k_c = CATCH_DISABLED_FRACTION · k_s`` (negligible to 1e-9 across the tensile
  regime) rather than 0 — ``PereverzevParams`` requires every field finite and > 0.
  ``pereverzev_F_star`` then correctly returns **NaN** (a slip bond has no catch
  peak), which the KU-2.5 catch-peak gate must learn to skip for slip ligands.
* **fibronectin → α5β1** keeps the full catch-slip = the platform KU-2.18 default.

All constants are published acceptance values (NOT fitting targets); every laminin
value is ``proxy``-flagged (direct α6β1–LN-111 force kinetics are genuinely absent
from the literature — PI-exp map dossier ``wf_1494c786-e79``).

Magic-Number Block
------------------
No tuned constants. Each kinetic value is literature-sourced (``source`` field).
``F_s = k_BT/x_β`` is the Bell-Evans definition (no free parameter); ``k_BT`` uses the
simulation body temperature ``T = 310.15 K``. AFM SMFS x_β were measured at RT
(~298 K, a ~4 % systematic on F_s, inside the bands). ``CATCH_DISABLED_FRACTION``
(1e-9) is a numerical sentinel meaning "no catch pathway", not a physics tuning knob.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ffn_sim.validation.pereverzev import (
    DEFAULT_PARAMS as KU_2_18_FN_PARAMS,  # FN α5β1 full catch-slip (Bangasser 2013)
    PereverzevParams,
    pereverzev_F_star,
)

# --- Physical constants (SI) --------------------------------------------------
K_B: float = 1.380649e-23          # J/K   Boltzmann (CODATA)
T_BODY: float = 310.15             # K     simulation body temperature (37 °C)
K_BT: float = K_B * T_BODY         # J     ≈ 4.2816e-21 J

# Numerical sentinel: a slip bond has NO catch pathway. PereverzevParams forbids
# k_c = 0 (requires every field finite and > 0), so the catch term is made
# negligible: k_c = CATCH_DISABLED_FRACTION · k_s gives catch/slip ≤ 1e-9 for all
# F ≥ 0, and pereverzev_F_star() → NaN (no peak). NOT a physics tuning constant.
CATCH_DISABLED_FRACTION: float = 1.0e-9


@dataclass(frozen=True, slots=True)
class LigandSpecies:
    """One ECM ligand presentation and its integrin-clutch kinetics.

    Attributes
    ----------
    name : str
        Registry key (e.g. ``"collagen_I"``).
    integrin : str
        Receptor (e.g. ``"alpha2beta1"``).
    bond_class : str
        ``"slip"`` (Bell-Evans, catch disabled) or ``"catch-slip"`` (Pereverzev).
    k_off0 : float | None
        Bell-Evans slip prefactor ``k_off(F=0)`` [s⁻¹]; ``None`` for catch-slip.
    x_beta : float | None
        Bell-Evans slip distance [m]; ``None`` for catch-slip.
    k_off0_range : tuple[float, float] | None
        Literature range for ``k_off0`` [s⁻¹] where reported.
    catch_params : PereverzevParams | None
        Full catch-slip parameter set; ``None`` for slip bonds.
    k_on : float | None
        On-rate [s⁻¹]; ``None`` = use the generic KU-2.4 default (off-rate-only
        ligand identity until a ligand-specific K_D→k_on calibration lands).
    source : str
        Primary literature anchor.
    confidence : str
        Adversarial-verification confidence (``"high (3-0)"`` / ``"med (proxy)"``).
    proxy : bool
        True if the kinetics are a documented proxy (not a direct measurement).
    notes : str
        Provenance / caveats.
    """

    name: str
    integrin: str
    bond_class: str
    k_off0: float | None
    x_beta: float | None
    k_off0_range: tuple[float, float] | None
    catch_params: PereverzevParams | None
    k_on: float | None
    source: str
    confidence: str
    proxy: bool
    notes: str

    @property
    def F_s(self) -> float | None:
        """Bell-Evans slip characteristic force ``F_s = k_BT / x_β`` [N].

        ``None`` for catch-slip bonds (use ``catch_params.F_s`` instead).
        """
        if self.x_beta is None:
            return None
        return K_BT / self.x_beta


# --- The ligand registry (literature-anchored; see docs/LIGAND_IDENTITY_FA_SPEC.md) ---
#
# F_s = k_BT/x_β at T = 310.15 K (k_BT = 4.2816e-21 J):
#   col-I single-cell  x_β=0.23 nm -> F_s = 18.6 pN
#   col-I locked-open  x_β=0.70 nm -> F_s =  6.1 pN
#   laminin α6β1       x_β=0.28 nm -> F_s = 15.3 pN  (≈ reported f_b 14.1±1.3 pN ✓)

LIGAND_REGISTRY: dict[str, LigandSpecies] = {
    # collagen-I → α2β1, slip. DEFAULT col-I anchor = the single-cell (cell-context)
    # measurement (Taubenberger 2007 MBoC 18:1634): k_off0 1.3 s⁻¹, x_β 2.3 Å,
    # lifetime 0.8±0.7 s, single-bond rupture 47±13 pN @500 pN/s.
    "collagen_I": LigandSpecies(
        name="collagen_I",
        integrin="alpha2beta1",
        bond_class="slip",
        k_off0=1.3,
        x_beta=0.23e-9,
        k_off0_range=(1.3 - 1.3, 1.3 + 1.3),  # 1.3±1.3 s⁻¹ as reported
        catch_params=None,
        k_on=None,
        source="Taubenberger 2007 MBoC 18:1634 (e06-09-0777), single-cell SMFS",
        confidence="high (3-0)",
        proxy=False,
        notes=(
            "Cell-context α2β1–col-I slip. Rupture 38→90 pN over 180–8800 pN/s "
            "(state loading rate). DEFAULT col-I anchor; see locked-open variant."
        ),
    ),
    # collagen-I locked-open variant (isolated I-domain SMFS): k_off0 0.44 s⁻¹,
    # x_β 0.7 nm (GPP control 11 s⁻¹ / 0.37 nm). Lower-affinity context.
    "collagen_I_locked_open": LigandSpecies(
        name="collagen_I_locked_open",
        integrin="alpha2beta1",
        bond_class="slip",
        k_off0=0.44,
        x_beta=0.70e-9,
        k_off0_range=None,
        catch_params=None,
        k_on=None,
        source="PMC3588017 (α2 I-domain AFM SMFS, locked-open high-affinity)",
        confidence="high (3-0)",
        proxy=False,
        notes="Isolated-domain SMFS variant of col-I; not the default (cell-context preferred).",
    ),
    # laminin-111 → α6β1, slip PROXY. Direct α6β1–LN-111 force kinetics are
    # GENUINELY ABSENT from the literature (PI-exp dossier wf_1494c786-e79). Shape
    # proxy = α7β1–invasin (PMC3882471): x_β 0.28 nm, f_b 12–15 pN (14.1±1.3),
    # k_off0 1.4–2.3 s⁻¹, single-barrier slip. Scaled WEAKER than the FN clutch
    # (rel-traction anchor: breast-epithelial traction lower on LN-111, P=0.016–0.028).
    "laminin_111": LigandSpecies(
        name="laminin_111",
        integrin="alpha6beta1",
        bond_class="slip",
        k_off0=1.85,  # representative midpoint of the 1.4–2.3 s⁻¹ proxy range
        x_beta=0.28e-9,
        k_off0_range=(1.4, 2.3),
        catch_params=None,
        k_on=None,
        source="PROXY: α7β1–invasin PMC3882471 (NOT laminin); rel-scaling PMC3391238",
        confidence="med (proxy)",
        proxy=True,
        notes=(
            "α6β1–LN-111 direct force kinetics absent in lit; α7β1-invasin slip shape, "
            "scaled weaker than FN. Affinity context: K_D ≈ 1–20 nM (LN-111 = α6β1's "
            "3rd-rank laminin). Catch-vs-slip unclassified; slip assumed (default)."
        ),
    ),
    # fibronectin → α5β1, full catch-slip = the platform KU-2.18 default (Bangasser
    # 2013; consistent with Kong 2009 catch 10–30 pN, peak 2–10 s @20–25 pN).
    "fibronectin": LigandSpecies(
        name="fibronectin",
        integrin="alpha5beta1",
        bond_class="catch-slip",
        k_off0=None,
        x_beta=None,
        k_off0_range=None,
        catch_params=KU_2_18_FN_PARAMS,
        k_on=None,
        source="KU-2.18 Bangasser 2013 BiophysJ 105:581; Kong 2009 JCB 185:1275",
        confidence="high (3-0)",
        proxy=False,
        notes="Platform baseline/benchmark; the existing clutch — bit-for-bit default.",
    ),
}

# PI-experiment condition → default ligand (PI-exp validation map §Layer-1).
DEFAULT_LIGAND_FOR_CONDITION: dict[str, str] = {
    "Bare": "collagen_I",   # untreated pV4D4, col-I coated
    "Pre": "collagen_I",    # surface-adsorbed col-I
    "Lam4": "laminin_111",  # + soluble laminin-111 (4 µg/mL) → α6β1
}


def pereverzev_params_for(species: str) -> PereverzevParams:
    """Return the ``PereverzevParams`` for a ligand, ready for ``IntegrinBondUpdater``.

    Catch-slip ligands (FN) return their full parameter set. Slip ligands (col-I,
    laminin) return a Bell-Evans slip mapped onto the Pereverzev form with the catch
    pathway disabled (``k_s = k_off0``, ``F_s = k_BT/x_β``,
    ``k_c = CATCH_DISABLED_FRACTION·k_s``, ``F_c = F_s``) — so ``pereverzev_k_off``
    runs them unchanged and ``pereverzev_F_star`` → NaN (no catch peak).

    Parameters
    ----------
    species : str
        A key in :data:`LIGAND_REGISTRY`.

    Returns
    -------
    PereverzevParams
    """
    lig = LIGAND_REGISTRY[species]
    if lig.bond_class == "catch-slip":
        assert lig.catch_params is not None
        return lig.catch_params
    # slip → Bell-Evans mapped to Pereverzev with catch disabled
    assert lig.k_off0 is not None and lig.x_beta is not None
    F_s = K_BT / lig.x_beta
    k_s = lig.k_off0
    return PereverzevParams(
        k_s=k_s,
        F_s=F_s,
        k_c=CATCH_DISABLED_FRACTION * k_s,
        F_c=F_s,
    )


def bell_evans_k_off(F_newton: float, species: str) -> float:
    """Bell-Evans slip off-rate ``k_off0·exp(F·x_β/k_BT)`` [s⁻¹] for a slip ligand.

    The faithful closed form for slip bonds, used by the Layer-1 acceptance gate
    (emergent k_off vs this oracle). Raises for a catch-slip species.

    Parameters
    ----------
    F_newton : float
        Tensile bond force [N], must be ≥ 0.
    species : str
        A slip-class key in :data:`LIGAND_REGISTRY`.
    """
    lig = LIGAND_REGISTRY[species]
    if lig.bond_class != "slip":
        raise ValueError(
            f"bell_evans_k_off: {species!r} is {lig.bond_class!r}, not slip; "
            "use pereverzev_k_off for catch-slip ligands."
        )
    if F_newton < 0.0:
        raise ValueError("bell_evans_k_off: F must be ≥ 0 (tensile loading).")
    assert lig.k_off0 is not None and lig.x_beta is not None
    return lig.k_off0 * math.exp(F_newton * lig.x_beta / K_BT)


# Catch peak (NaN for slip ligands — they have none) for convenience/tests.
def catch_peak_force_for(species: str) -> float:
    """Analytic catch-peak force F* [N] for a ligand, or NaN if slip (no peak)."""
    return float(pereverzev_F_star(pereverzev_params_for(species)))
