"""Thermal worm-like-chain (Marko-Siggia) force-extension for FF fibrillar segments — the single source of law.

PI-gated, DEFAULT-OFF, PROVISIONAL. The athermal linear-spring collagen reproduces the reference-concentration
modulus but the G'(c) exponent is n≈1 (density-linear) vs the literature n≈2 — KB-diagnosed as a THERMAL
semiflexible effect (KB-1.30 G0~kB·T·L_p²/ξ⁵, KB-1.7 ξ~c^(−1/2) ⇒ c^2.5; CSCALING_REGIME_FINDING.md). This module
is the shared, host-side + kernel-mirrored WLC law used by the relax kernel, the virial (`ecm_material_stress`),
and the energy route (`_elastic_energy`), so all three read ONE tension.

The extensible-WLC (Marko-Siggia entropic branch below the enthalpic crossover x_max, a finite EA backbone wall
above it — C¹-continuous, so the 1/(1−x)² divergence is removed with finite max stiffness EA/Lc):

    x = L / Lc                                              (fractional extension; L end-to-end, Lc contour)
    F(x) = (kB·T / L_p) · [ 1/(4(1−x)²) − 1/4 + x ]         for x ≤ x_max   [pN]
    F(x) = F(x_max) + (EA/Lc)·(L − x_max·Lc)                for x > x_max   (enthalpic wall)

The crux is the CONTOUR LENGTH per segment: the thermal slack that carries the c-dependence is the EMERGENT
geometric mesh ξ ~ c^(−1/2) (NOT the c-independent 0.5 µm discretization seg_rest — with L_p≫seg a per-segment
slack from seg_rest is nearly rigid and gives c^1, the classic failure). Kratky-Porod stiff-limit:
    seg_Lc = seg_rest · (1 + ξ/(6·L_p)),   ξ = (V / L_total)^(1/2) [3D],   x₀ = 1/(1+ξ/(6L_p))  (baseline op-point).

Constants are all derived/literature-anchored (kB·T=U.KBT 310K; L_p, EA from the spec; ξ emergent; 1/6 the
Kratky-Porod coefficient; x_max from f'_WLC=EA) — NOT tuned to any exponent.

References: Marko & Siggia 1995 (Macromolecules 28:8759); MacKintosh/Storm 2005 (K~σ^{3/2}); KB-1.30/1.7/1.3.
"""

from __future__ import annotations

import numpy as np


def wlc_x_crossover(EA: float, Lp: float, kBT: float) -> float:
    """Entropic→enthalpic crossover x_max where the WLC differential stiffness f'_WLC(x) equals the backbone EA.

    f'_WLC(x) = (kB·T/L_p)·[1/(2(1−x)³) + 1] = EA  ⇒  (1−x_max)³ = kB·T / (2(EA·L_p − kB·T)).
    Above x_max the force is the linear EA wall, so the max stiffness is finite (= EA/Lc) — no divergence.
    Physics-fixed by (EA, L_p, kB·T); not a tuned number."""
    denom = 2.0 * (EA * Lp - kBT)
    if denom <= 0.0:
        return 0.99                                        # degenerate (EA≲kBT/Lp) → soft cap; not reached for collagen
    x = 1.0 - (kBT / denom) ** (1.0 / 3.0)
    return float(min(max(x, 0.5), 0.999))                  # x_hard=0.999 numerical guard (non-binding for collagen)


def wlc_tension_np(L, Lc, Lp: float, EA: float, kBT: float, x_max: float) -> np.ndarray:
    """Extensible-WLC segment tension F(L) [pN] (host oracle; the kernel mirrors this exactly)."""
    L = np.asarray(L, float); Lc = np.asarray(Lc, float)
    x = L / np.maximum(Lc, 1e-12)
    xe = np.minimum(x, x_max)                              # entropic branch capped at the crossover
    F_ent = (kBT / Lp) * (1.0 / (4.0 * (1.0 - xe) ** 2) - 0.25 + xe)
    F_max = (kBT / Lp) * (1.0 / (4.0 * (1.0 - x_max) ** 2) - 0.25 + x_max)
    F_wall = F_max + (EA / np.maximum(Lc, 1e-12)) * (L - x_max * Lc)   # enthalpic backbone above x_max
    return np.where(x <= x_max, F_ent, F_wall)


def wlc_energy_np(L, Lc, Lp: float, EA: float, kBT: float, x_max: float) -> np.ndarray:
    """Extensible-WLC strain energy U(L) [pN·µm] = ∫F dR (host; for the energy-route G=2ΔU/(Vγ²))."""
    L = np.asarray(L, float); Lc = np.asarray(Lc, float)
    x = L / np.maximum(Lc, 1e-12)
    xe = np.minimum(x, x_max)
    # ∫F dR = Lc·∫F dx ; ∫[1/(4(1−x)²) − 1/4 + x]dx = 1/(4(1−x)) − x/4 + x²/2
    U_ent = Lc * (kBT / Lp) * (1.0 / (4.0 * (1.0 - xe)) - xe / 4.0 + xe ** 2 / 2.0)
    F_max = (kBT / Lp) * (1.0 / (4.0 * (1.0 - x_max) ** 2) - 0.25 + x_max)
    R_max = x_max * Lc
    U_at_max = Lc * (kBT / Lp) * (1.0 / (4.0 * (1.0 - x_max)) - x_max / 4.0 + x_max ** 2 / 2.0)
    dR = L - R_max
    U_wall = U_at_max + F_max * dR + 0.5 * (EA / np.maximum(Lc, 1e-12)) * dR ** 2
    return np.where(x <= x_max, U_ent, U_wall)


def geometric_mesh_um(seg_rest: np.ndarray, box_lo, box_hi, dim: int = 3) -> float:
    """Emergent geometric mesh ξ = (V/L_total)^(1/2) [3D] or (A_xy/L_total) [2D] — the confinement/slack scale.

    L_total = Σ seg_rest ∝ n_fibers ∝ c, so ξ ~ c^(−1/2) EMERGES (matches spec.conc_exponent=0.5). Measured from
    the BUILT network (least gameable), NOT an input."""
    L_total = float(np.sum(seg_rest))
    if L_total <= 0.0:
        return 0.0
    lo = np.asarray(box_lo, float); hi = np.asarray(box_hi, float)
    if dim == 2:
        A = float((hi[0] - lo[0]) * (hi[1] - lo[1]))
        return A / L_total
    V = float(np.prod(hi - lo))
    return (V / L_total) ** 0.5


def segment_contour_lengths(seg_rest: np.ndarray, Lp: float, xi: float) -> np.ndarray:
    """Per-segment contour length seg_Lc = seg_rest·(1 + ξ/(6·L_p)) (Kratky-Porod stiff-limit thermal slack).

    Distributes the confinement span's uniform contour-excess over its discretization segments. Baseline operating
    point x₀ = seg_rest/seg_Lc = 1/(1+ξ/(6L_p)) — the physiological resting extension (small, balanced pre-tension),
    NOT a floppy or pre-stressed null."""
    return np.asarray(seg_rest, float) * (1.0 + xi / (6.0 * Lp))


__all__ = ["wlc_x_crossover", "wlc_tension_np", "wlc_energy_np", "geometric_mesh_um", "segment_contour_lengths"]
