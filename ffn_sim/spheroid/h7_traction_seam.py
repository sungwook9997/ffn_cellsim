"""H.7 fine-grained single-cell traction → Layer-2 CBM active-traction SCALE-BRIDGE (consumer).

WHY THIS EXISTS
---------------
The Layer-2 center-based CBM line is CLOSED (2026-06-05): it reproduced the PI spreading law
A/A0 = a + b/R + c/R² in FORM (r² = 0.998, zero-calibration) but its MAGNITUDE is bounded ~5–9×
under the PI medians. The ``b/R`` term is the **single-cell active traction** that drives
spreading; Layer-2 has only ever set it from HEURISTIC anchors:
  * ``motility_bridge.resolve_active_traction``  — whole-cell migration speed × γ ⇒ ~1.6 nN,
                                                   protrusion velocity × γ ⇒ ~9.4 nN (a 6× fork);
  * ``ligand_traction.py``                       — ``T_ref = 2.5 nN`` clutch-scale GUESS, with the
                                                   in-code caveat *"no direct MCF7 single-cell
                                                   traction in the literature"*.

H.7 has now PRODUCED that missing quantity mechanistically. A graded-polarity SARCOMERIC ventral
stress fiber, built from §9-corrected continuous-stroke myosin on FA anchors, rectifies the motor
into a same-seed coherent contractile traction of **+131 ± 8 pN per single SF** (16σ; +2.67 pN per
engaged head ≈ the native cross-bridge F_stall). This module is the CONSUMER side of the
H.7→Layer-2 seam: it turns that fine-grained, particle/bond traction into the coarse CBM per-cell
``f_active`` parameter — replacing the heuristic with a mechanistic anchor, **never fitting to the
PI A/A0 (overlay-only HARD rule)**.

THE BRIDGE (dimensional, no new mechanism, no integrator change)
---------------------------------------------------------------
H.7 measures a per-SF traction in pN; Layer-2's ``f_active`` is a per-CELL traction in nN. The
conversion is the cell's stress-fiber count::

    f_cell = per_SF_traction · N_SF                         [pN · (–) = pN → nN]

with ``N_SF`` ANCHORED, not tuned: each ventral SF spans two focal adhesions (barbed/+ ends at the
FAs), so ``N_SF ≈ FA_count / 2`` (Hotulainen & Lappalainen 2006 JCB; FA_count = 43 from KU-2.4,
``configs/phase1_h4.yaml``) ⇒ N_SF ≈ 21.5, consistent with the ~20 ventral SF/cell reported by
Hotulainen-Lappalainen. This yields the TEST-DENSITY aggregate (the H.7 fiber is a deliberately
under-dense test unit: 16 minifilaments / 3 sarcomeres / 6 µm, ~1–2 orders under physiological
per-SF minifilament density — H7_SF_ARRAY_TRACTION_SCALEUP §3).

The PHYSIOLOGICAL single-cell traction is anchored independently to the nearest literature TFM
datum — Gil-Redondo 2023 (Microsc Res Tech, DOI 10.1002/jemt.24368): MCF-7 total contractility
**102 ± 59 nN** over a 1822 ± 886 µm² contact area (~63 Pa mean stress). H.7's scaling law
(F_aggregate ∝ engaged-head count × 2.67 pN) reaches this at physiological minifilament density;
the ratio physiological ÷ test-density = the documented density gap (~×30–40 ≈ the sanctioned ×40
mesoscale coarse-graining), the SAME overlap/density physics as the cortical-γ density datum, NOT a
sign or mechanism error (the sign + rectification are 16σ-confirmed).

So the seam supplies TWO mechanistic anchors for the decisive fork, both *measured/derived*, never
fitted:
  1. ``f_cell_test``      — raw H.7 mechanism at the validated TEST-fiber density (~2.8 nN);
  2. ``f_cell_physio``    — the same mechanism scaled to physiological density (Gil-Redondo, 102 nN).

SEAM RECORD (when H.7 lands the array aggregate)
------------------------------------------------
H.7's ``H7_SF_ARRAY_TRACTION_SCALEUP §5`` designs a seam record
``outputs/h7/production/h7_traction_seam.json`` written by the H.7 SF-array GPU run. When present
(dropped into ``outputs/layer2/seam_inputs/`` or readable in the H.7 production dir), this consumer
reads the MEASURED array aggregate from it instead of the provisional per-SF × N_SF estimate. Until
then it is PROVISIONAL from the decisive single-SF value (status flagged on the result).

Sanity Gate
-----------
- Dimensional: per_SF [N] · N_SF [–] = f_cell [N]; physiological anchor [N]; ratio [–]. ✓
- Boundary: per_SF = 0 ⇒ f_cell = 0 (no traction, reduces to G1 cohesion-locked). N_SF ≥ 1.
- Sign: per_SF > 0 is CONTRACTILE (the sarcomeric rectification); mixed-polarity SMOKE was −
  (slackening) and is NOT consumed — only the decisive sarcomeric (+) value bridges.
- Provenance: per_SF from the vendored H.7 artifact (or seam record); N_SF from the FA-pairing
  anchor; physiological from Gil-Redondo 2023. NO value fitted to the PI A/A0.
- Measurement: ``density_gap_factor = f_cell_physio / f_cell_test`` is reported (the honest
  fine-grained-prediction residual), never absorbed into a tuned magnitude.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "ResolvedH7Traction",
    "resolve_h7_traction",
    "FA_COUNT_PER_CELL",
    "N_SF_PER_CELL",
    "GIL_REDONDO_TOTAL_CONTRACTILITY_N",
    "GIL_REDONDO_CONTACT_AREA_M2",
]

# ---- Literature anchors (HARD rule: derivable / cited, never tuned to a gate) ----

# Focal adhesions per cell — KU-2.4 (configs/phase1_h4.yaml; brief 20–50 nascent + 5–10 mature →
# 43 = 35 + 8 in the H.7 SF-array doc §2). Used only to DERIVE the SF count.
FA_COUNT_PER_CELL: int = 43

# Ventral stress fibers per cell. Each ventral SF spans TWO focal adhesions (barbed/+ ends pinned
# at the FAs), so N_SF ≈ FA_count / 2 (Hotulainen & Lappalainen 2006, JCB 173:383). = 21.5,
# consistent with the ~20 ventral SF/cell that paper reports. DERIVED, not a free knob.
N_SF_PER_CELL: float = FA_COUNT_PER_CELL / 2.0

# Physiological single-cell traction anchor — nearest literature TFM datum to the PI platform.
# Gil-Redondo 2023 (Microsc Res Tech, DOI 10.1002/jemt.24368), Table 1, MCF-7 control n=37:
# total contractility 102 ± 59 nN over a 1822 ± 886 µm² contact area (≈ 63 Pa mean stress).
# OVERLAY anchor (the PI platform is MCF-7 on pV4D4/col-I integrin-β1; Gil-Redondo is MCF-7 on
# fibronectin/PDMS) — used to set the physiological SCALE, never fitted to the PI A/A0.
GIL_REDONDO_TOTAL_CONTRACTILITY_N: float = 102.0e-9
GIL_REDONDO_CONTACT_AREA_M2: float = 1822.0e-12

# Decisive single sarcomeric-SF traction (provisional fallback when no seam record / array
# aggregate is present). Source: H.7 h7_sf_2c_sarcomeric_gpu.json (loop23), vendored read-only at
# outputs/layer2/seam_inputs/. Read from that file at runtime; this is the documented value.
_DECISIVE_PER_SF_pN: float = 131.08

_SEAM_INPUTS = Path(__file__).resolve().parents[1] / "outputs" / "layer2" / "seam_inputs"
_H7_PROD = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "production"


@dataclass(frozen=True)
class ResolvedH7Traction:
    """H.7-anchored CBM single-cell active traction (SI). See module docstring for the bridge."""

    per_sf_traction: float        # N    decisive single sarcomeric-SF coherent contractile traction
    per_sf_sem: float             # N    same-seed coherent differential SEM (±8 pN floor)
    n_sf: float                   # –    ventral SF / cell (FA-pairing anchor)
    f_cell_test: float            # N    per_sf · n_sf  (raw mechanism, TEST-fiber density)
    f_cell_physio: float          # N    physiological-density anchor (Gil-Redondo total contractility)
    contact_area: float           # m²   MCF-7 spread footprint (Gil-Redondo)
    stress_test_Pa: float         # Pa   f_cell_test / contact_area
    stress_physio_Pa: float       # Pa   f_cell_physio / contact_area (≈ 63 Pa, Gil-Redondo)
    density_gap_factor: float     # –    f_cell_physio / f_cell_test (honest density residual)
    source: str                   # "seam_record" | "provisional_single_sf"
    provisional: bool             # True until the H.7 array aggregate / seam record lands
    provenance: str


def _find_seam_record() -> Path | None:
    """Locate the H.7-supplied seam record, if it has landed (vendored copy preferred)."""
    for d in (_SEAM_INPUTS, _H7_PROD):
        p = d / "h7_traction_seam.json"
        if p.exists():
            return p
    return None


def _read_decisive_per_sf() -> tuple[float, float]:
    """Read the decisive single-SF coherent traction (pN) + SEM from the vendored H.7 artifact.

    Falls back to the documented constant if the vendored file is absent (so the bridge never
    silently breaks), but the vendored artifact is the source of truth.
    """
    p = _SEAM_INPUTS / "h7_sf_2c_sarcomeric_gpu.json"
    if p.exists():
        d = json.loads(p.read_text())
        return float(d["coherent_differential_pN"]), float(d.get("coherent_differential_sem_pN", 0.0))
    return _DECISIVE_PER_SF_pN, 8.2


def resolve_h7_traction(*, n_sf: float = N_SF_PER_CELL) -> ResolvedH7Traction:
    """Resolve the H.7-anchored CBM single-cell active traction (pure arithmetic bridge, no sim).

    Reads the H.7-supplied seam record ``h7_traction_seam.json`` if present (measured array
    aggregate); otherwise computes the per-cell aggregate PROVISIONALLY from the decisive
    single-SF traction × N_SF. The physiological anchor (Gil-Redondo) is always supplied.

    Args:
        n_sf: ventral stress fibers per cell (default = FA-pairing anchor ~21.5). Override only to
            test another anchored count.

    Returns:
        ``ResolvedH7Traction`` with the test-density and physiological per-cell tractions, the
        contact-area stress, the density-gap factor, and the provenance/provisional status.
    """
    if n_sf <= 0.0:
        raise ValueError("n_sf must be > 0.")
    area = GIL_REDONDO_CONTACT_AREA_M2
    f_physio = GIL_REDONDO_TOTAL_CONTRACTILITY_N

    seam = _find_seam_record()
    if seam is not None:
        d = json.loads(seam.read_text())
        per_sf = float(d["per_SF_traction_N"])
        per_sf_sem = float(d.get("per_SF_traction_sem_N", 0.0))
        # When H.7 supplies the measured array aggregate, prefer it for f_cell_test.
        f_test = float(d.get("single_cell_total_traction_N", per_sf * n_sf))
        n_sf_used = float(d.get("n_SF", n_sf))
        source, provisional = "seam_record", False
        prov = str(d.get("provenance", "H.7 seam record h7_traction_seam.json"))
    else:
        per_sf_pN, per_sf_sem_pN = _read_decisive_per_sf()
        per_sf = per_sf_pN * 1e-12
        per_sf_sem = per_sf_sem_pN * 1e-12
        f_test = per_sf * n_sf
        n_sf_used = n_sf
        source, provisional = "provisional_single_sf", True
        prov = (
            "PROVISIONAL: H.7 decisive single sarcomeric-SF +131 pN "
            "(h7_sf_2c_sarcomeric_gpu.json, loop23) × N_SF (FA-pairing anchor); "
            "awaiting H.7 SF-array GPU aggregate / h7_traction_seam.json. "
            "Physiological scale = Gil-Redondo 2023 MCF-7 102 nN / 63 Pa."
        )

    return ResolvedH7Traction(
        per_sf_traction=per_sf,
        per_sf_sem=per_sf_sem,
        n_sf=n_sf_used,
        f_cell_test=f_test,
        f_cell_physio=f_physio,
        contact_area=area,
        stress_test_Pa=f_test / area,
        stress_physio_Pa=f_physio / area,
        density_gap_factor=f_physio / f_test if f_test > 0 else float("inf"),
        source=source,
        provisional=provisional,
        provenance=prov,
    )
