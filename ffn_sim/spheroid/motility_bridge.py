"""Single-cell H.5 lamellipodium → Layer-2 CBM active-motility SCALE-BRIDGE.

The Layer-2 CBM reproduces the PI spreading-law FORM (A/A0 = a + b/R + c/R², r²≈0.998 across the
full PI R₀ range) but UNDER-SPREADS in magnitude ~5–9× vs the PI MCF7 data. The cohesion side is
anchored TWO independent ways (Iturri-2020 de-adhesion AND the D2 surface-tension γ=0.57 mN/m in
the KU-3.5 band), so it is NOT "over-strong by an unanchored amount." The remaining free scale is
the **active per-cell traction/self-propulsion**, whose absolute magnitude was an unanchored
heuristic (the B1 ``STABLE_TRACTION_CEILING`` ≤3 nN; ``spreading.py``/``ligand_traction.py`` flag
this explicitly). This module ANCHORS it the way cohesion was anchored — from a measured
single-cell quantity, never by fitting to the PI A/A0 (overlay-only HARD rule).

The bridge (no new mechanism, no integrator change)
---------------------------------------------------
A CBM cell is overdamped: the FROZEN BAOAB step is ``r(t+dt) = r(t) + (F/γ)·dt`` ⇒ a cell driven
by a constant active force moves at terminal velocity ``v = f_active / γ_cell``. Inverting, an
anchored cell speed ``v0`` sets the active force the existing ``ActiveMotility`` / edge-traction
knobs consume::

    f_active = v0 · γ_cell        [m/s · N·s/m = N]

``γ_cell`` is the per-cell MIGRATION drag = clutch-ensemble friction ``n_eng·κ_clutch/k_off``
(KU-2.18 Bangasser 2013; ``params.resolve_layer2`` ⇒ 0.30 N·s/m), the SAME single-cell anchor the
cohesion bridge used — NOT water-Stokes.

Two candidate single-cell speed anchors (their RATIO is the question)
--------------------------------------------------------------------
1. **Protrusion velocity (fine-grained, in-tree, the H.5 mechanism's own constant):** the Bieling
   2016 unloaded barbed-end elongation velocity ``v0 = k_elong⁰·δ_elong`` — read here from the
   SINGLE-CELL config ``configs/phase1_h5.yaml`` (the source of truth the H.5 runtime uses), so
   the bridge tracks the platform's own protrusion mechanism. With k_elong⁰=11.6 s⁻¹,
   δ_elong=2.7 nm ⇒ v0 ≈ 31.3 nm/s = 1.88 µm/min ⇒ **f_active ≈ 9.4 nN**.
   ⚠️ NOT ``k_elong⁰·rest_length`` (=5.8 µm/s) — that is a network OVERGROWTH rate, not a
   barbed-end translocation, and is nonsense as a cell speed.
2. **Whole-cell migration speed (MCF7-specific, OVERLAY-validation):** the measured MCF7
   single-cell speed 19.0±6.6 µm/h ≈ 0.32 µm/min (DeepBIT 2026, N=2600; in the focal-adhesion
   epithelial band 0.2–1 µm/min) ⇒ **f_active ≈ 1.6 nN**, in the B1 stable band.

The two differ by ~6× — the **protrusion machinery is ~6× faster than whole-cell translocation**
(the expected motility/clutch loss). That ~6× ≈ the observed ~5–9× magnitude gap: i.e. the
under-spread is consistent with the model effectively running at the slow whole-cell speed while
the experiment's collective spread carries closer to the protrusion scale. Which speed is THE
anchor (vs the validation overlay), and whether to lift the 3 nN ejection ceiling to admit the
protrusion value, is a PI call (anchor choice = science; ceiling lift brushes the integrator
freeze). This resolver computes both and asserts nothing it shouldn't — it surfaces the fork.

Sanity Gate
-----------
- Dimensional: v0 [m/s]; γ_cell [N·s/m]; f_active = v0·γ_cell [N]. ✓
- Boundary: v0=0 ⇒ f_active=0 (no spreading, reduces to G1). Both speeds > 0.
- Sign/sense: larger v0 ⇒ larger f_active ⇒ more spread. Monotone.
- Provenance: k_elong⁰, δ_elong read from the H.5 SoT config (no duplicated constant); the MCF7
  speed is a cited literature overlay anchor. NO value is fitted to the PI A/A0.
- Measurement: ``active_ratio = f_protrusion/f_wholecell`` is the active-side scale ambiguity,
  reported against the observed magnitude gap (overlay).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from ffn_sim.spheroid.ligand_traction import STABLE_TRACTION_CEILING
from ffn_sim.spheroid.params import ResolvedL2

__all__ = ["ResolvedMotility", "resolve_active_traction", "V0_MCF7_WHOLECELL"]

# Single-cell H.5 config = the source of truth for the Bieling protrusion constants (so the
# bridge reads the SAME k_elong⁰, δ_elong the H.5 lamellipodium runtime uses — a true bridge).
_H5_CFG = Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml"

# MCF7 whole-cell single-cell migration speed — OVERLAY VALIDATION anchor (NOT fitted to PI A/A0).
# 19.0 ± 6.6 µm/h (DeepBIT 2026, N=2600 MCF7 cells; bioRxiv 2025.04.01.646705), inside the
# focal-adhesion epithelial band 0.2–1 µm/min (PMC7071959). A WHOLE-CELL translocation speed —
# distinct from the sub-cellular protrusion velocity (anchor 1).
V0_MCF7_WHOLECELL: float = 19.0e-6 / 3600.0  # m/s  (= 0.317 µm/min)


@dataclass(frozen=True)
class ResolvedMotility:
    """Anchored CBM active-traction candidates (SI). See module docstring for the bridge."""

    gamma_cell: float                  # N·s/m  per-cell migration drag (clutch-ensemble)
    v0_protrusion: float               # m/s    Bieling unloaded barbed-end protrusion velocity
    f_protrusion: float                # N      = v0_protrusion · gamma_cell (fine-grained anchor)
    v0_wholecell: float                # m/s    MCF7 whole-cell migration (overlay validation)
    f_wholecell: float                 # N      = v0_wholecell · gamma_cell
    ceiling: float                     # N      B1 STABLE_TRACTION_CEILING
    protrusion_exceeds_ceiling: bool   # f_protrusion > ceiling (→ needs sub-stepping/dt, PI/D3)
    active_ratio: float                # f_protrusion / f_wholecell ≈ the active-side scale gap


def _load_bieling_v0() -> float:
    """Unloaded barbed-end protrusion velocity v0 = k_elong⁰·δ_elong from the H.5 SoT config."""
    cfg = yaml.safe_load(_H5_CFG.read_text())
    # phase1_h5.yaml nests the elongation constants under the lamellipodium/dynamics block; find
    # k_elong_0 [1/s] and delta_elong [m] wherever they live (single source of truth).
    def _find(key: str):
        stack = [cfg]
        while stack:
            node = stack.pop()
            if isinstance(node, dict):
                if key in node:
                    return float(node[key])
                stack.extend(node.values())
        raise KeyError(f"{key} not found in {_H5_CFG}")
    k_elong_0 = _find("k_elong_0")     # 1/s   (Bieling 2016)
    delta_elong = _find("delta_elong")  # m     (Bieling 2016 G-actin attachment length)
    return k_elong_0 * delta_elong      # m/s


def resolve_active_traction(
    resolved: ResolvedL2,
    *,
    v0_wholecell: float = V0_MCF7_WHOLECELL,
) -> ResolvedMotility:
    """Resolve the two anchored CBM active-traction candidates (no sim; pure arithmetic bridge).

    Args:
        resolved: resolved Layer-2 params (provides ``gamma_cell``, the clutch-ensemble drag).
        v0_wholecell: whole-cell migration speed for the overlay-validation arm (default MCF7
            DeepBIT 2026). Override only to test another cell line's measured speed.

    Returns:
        ``ResolvedMotility`` with the protrusion-anchored and whole-cell-anchored f_active, the
        B1 ceiling, the exceed flag, and the protrusion/whole-cell ratio (the active-side gap).
    """
    g = float(resolved.gamma_cell)
    if g <= 0.0:
        raise ValueError("gamma_cell must be > 0 to map a speed to a force.")
    v0_prot = _load_bieling_v0()
    if v0_prot <= 0.0 or v0_wholecell <= 0.0:
        raise ValueError("anchor speeds must be strictly positive.")
    f_prot = v0_prot * g
    f_whole = v0_wholecell * g
    return ResolvedMotility(
        gamma_cell=g,
        v0_protrusion=v0_prot,
        f_protrusion=f_prot,
        v0_wholecell=float(v0_wholecell),
        f_wholecell=f_whole,
        ceiling=float(STABLE_TRACTION_CEILING),
        protrusion_exceeds_ceiling=bool(f_prot > STABLE_TRACTION_CEILING),
        active_ratio=f_prot / f_whole,
    )
