r"""Non-muscle myosin IIA bipolar minifilament — explicit multi-head + LINEAR force-velocity (FF, Warp).

The live FF myosin (``network_warp.myosin_kernel``) is a lumped constant-force dipole whose magnitude f_myo is
SWEPT as a controlled variable. This makes it the fine-grained Stam-Hocky minifilament with a DERIVED stall
force and the correct force-velocity law.

⚠ Force-velocity is LINEAR, not Hill (PI-ratified 2026-07-07 off the lit-ingest): a non-muscle-myosin-IIA Hill
hyperbola does NOT exist in the literature — a/F₀≈0.25 is muscle-only (Hill 1938). The sourced non-muscle law is
the linear/affine ``v(F) = v₀·(1 − F/F_s)`` (Freedman 2017 AFINES; Cytosim; Tam 2021). So the CLAUDE.md
"linear-stall = Wrong / Hill = Right" row is MUSCLE-based and does not bind non-muscle myosin — the linear law is
the correct, sourced choice here. (CLAUDE.md to be annotated.)

Stam-Hocky bipolar minifilament (KB-3.18): ~N_heads myosin split into two anti-parallel half-filaments; each
head has per-head stall ~F_head, duty ratio r (fraction bound; ~0.1 non-muscle, Kovács 2003). The per-side
ENSEMBLE stall = N_side · F_head (all engaged at isometric load) → per-minifilament 50–100 pN (KB-3.18). The
contractile force between the two actin anchors it bridges is DERIVED from this, not swept.

No magic numbers: N_heads, F_head from KB-3.18; duty from Kovács 2003; v₀ from KB-PIV-4 (v_u=120 nm/s). The
γ-floor experiment now measures the EMERGENT γ from these explicit minifilaments (the known ~530× deficit is
reported, not tuned away — see the honest-deviation clause in the build plan DoD).

Gates (analytic, non-Hill): (1) ensemble stall F_s = N_side·F_head lands in the KB-3.18 per-minifilament band
50–100 pN; (2) LINEAR FV: v(0)=v₀, v(F_s)=0, exactly affine in between (a Hill fit would show curvature a/F_s>0 —
here a/F_s→∞, i.e. straight); (3) force-free at v=v₀ (unloaded minifilament exerts ~0 contractile force).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

# KB-3.18 minifilament composition + Kovács 2003 duty + KB-PIV-4 unloaded velocity
N_HEADS_MINIFIL = 60           # MOTOR HEADS per bipolar minifilament = KB-3.18 ~30 myosin-II dimers × 2 heads
F_HEAD_PN = 2.0                # per-head stall force [pN] (KB-3.18 ~2 pN)
DUTY_NM2A = 0.1                # non-muscle IIA duty ratio (Kovács 2003; band 0.1–0.3)
V0_UM_S = 0.120               # unloaded shortening velocity [µm/s] (KB-PIV-4 v_u = 120 nm/s)


@dataclass
class MyosinMinifilament:
    n_heads: int = N_HEADS_MINIFIL
    f_head_pn: float = F_HEAD_PN
    duty: float = DUTY_NM2A
    v0_um_s: float = V0_UM_S

    @property
    def n_side(self) -> int:
        """Heads per half-filament (the two anti-parallel sides of the bipolar minifilament)."""
        return self.n_heads // 2

    @property
    def f_stall_pn(self) -> float:
        """Ensemble isometric stall force [pN] = N_side · F_head (all heads engaged under load)."""
        return self.n_side * self.f_head_pn

    @property
    def f_resting_pn(self) -> float:
        """Resting contractile force at the duty ratio (fraction bound, unloaded) [pN] = N_side·duty·F_head."""
        return self.n_side * self.duty * self.f_head_pn


def resolve_myosin() -> MyosinMinifilament:
    return MyosinMinifilament()


def linear_force_velocity(F_pn, m: MyosinMinifilament | None = None):
    """Non-muscle LINEAR force-velocity ``v(F) = v₀·(1 − F/F_s)`` [µm/s]; clamped ≥0 past stall (no lengthening)."""
    m = m or resolve_myosin()
    v = m.v0_um_s * (1.0 - np.asarray(F_pn, float) / m.f_stall_pn)
    return np.maximum(v, 0.0)


@wp.kernel
def minifilament_kernel(pos: wp.array(dtype=wp.vec3d), links: wp.array(dtype=wp.int32, ndim=2),
                        v_slide: wp.array(dtype=wp.float64), f_stall: wp.float64, v0: wp.float64,
                        force: wp.array(dtype=wp.vec3d)):
    """Bipolar-minifilament contractile force between actin anchors ``links[t]``, from the INVERTED linear
    force-velocity: at the current inter-anchor sliding velocity ``v_slide`` the minifilament exerts
    ``F = F_stall·(1 − v_slide/v₀)`` (isometric v=0 ⇒ full stall F_stall; free-running v=v₀ ⇒ 0). Pulls i,j
    together. In the quasi-static implicit solve the equilibrium is v→0 ⇒ the stall prestress (the γ source)."""
    t = wp.tid()
    i = links[t, 0]; j = links[t, 1]
    d = pos[j] - pos[i]
    L = wp.length(d)
    if L < wp.float64(1e-12):
        return
    Fmag = f_stall * (wp.float64(1.0) - v_slide[t] / v0)
    if Fmag < wp.float64(0.0):
        Fmag = wp.float64(0.0)                             # past stall: no active lengthening force
    f = (Fmag / L) * d
    wp.atomic_add(force, i, f)                             # i toward j (shorten)
    wp.atomic_add(force, j, -f)


__all__ = ["MyosinMinifilament", "resolve_myosin", "linear_force_velocity", "minifilament_kernel",
           "N_HEADS_MINIFIL", "F_HEAD_PN", "DUTY_NM2A", "V0_UM_S"]
