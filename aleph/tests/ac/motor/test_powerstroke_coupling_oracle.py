"""Self-test of the power-stroke -> force coupling oracle (pure NumPy — no Warp/CUDA).

THE gate that catches the I3 decoupling defect on the dev Mac: it asserts the walked abscissa (the myosin
power stroke) actually ENTERS the crossbridge force, so a stepping head makes a directed contractile force and
force-velocity self-limits at F_stall. The earlier suite could not catch this — it never exercised the
mechanics<->kinetics coupling (the force kernels were only run on the gbook). This oracle mirrors the fixed
device kernels bit-for-formula, so the coupling is now a LOCAL gate.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.motor.hill_fv_analytic import hill_force
from aleph.components.motor.minifilament_topology import working_stroke_strain
from aleph.components.motor.powerstroke_analytic import (
    crossbridge_force,
    stall_abscissa,
    step_engaged_head,
)

# illustrative (GAP/provisional) params — the SHAPE/coupling gates are magnitude-independent (params_i0b3.yaml)
K_XB = 200.0     # pN/um (GAP: 100-1000 band; MASTER knob)
F_STALL = 2.0    # pN    (GAP: 0.5 vs 2.0)
V0 = 0.2         # um/s  (GAP: 0.12 vs 0.2)
R0_XB = 0.0      # um    (~0: a bound head sits on its actin site)
WALK = np.array([1.0, 0.0, 0.0])   # unit barbed-end direction


# ── property (a): a walking head makes a DIRECTED contractile force of the correct sign ─────────────
def test_walking_head_directed_contractile_force() -> None:
    """s > 0 (a head that has stepped) drags actin along -walk_dir (contractile) with |F| = k_xb (s - r0_xb)."""
    s = 0.008  # um walked (8 nm), below the stall strain
    xb = crossbridge_force(head_pos=(0, 0, 0), anchor_pos=(0, 0, 0), abscissa=s,
                           walk_dir=WALK, k_xb=K_XB, r0_xb=R0_XB)
    # actin is dragged the way OPPOSITE the walk -> contraction (negative projection on walk_dir)
    assert xb["f_actin_along_walk"] < 0.0
    # magnitude equals the crossbridge tension k_xb (s - r0_xb), directed
    assert np.linalg.norm(xb["f_on_actin"]) == pytest.approx(K_XB * (s - R0_XB), rel=1e-12)
    # the force is along the axis (no spurious transverse component here)
    assert xb["f_on_actin"][1] == pytest.approx(0.0, abs=1e-12)
    assert xb["f_on_actin"][2] == pytest.approx(0.0, abs=1e-12)


def test_newton_third_law_on_crossbridge() -> None:
    """The crossbridge force is equal-and-opposite on head and actin (Newton's 3rd law) to machine precision."""
    xb = crossbridge_force(head_pos=(0.001, 0.002, 0.0), anchor_pos=(0, 0, 0), abscissa=0.006,
                           walk_dir=WALK, k_xb=K_XB)
    assert xb["newton3rd_residual"] == pytest.approx(0.0, abs=1e-13)
    assert np.allclose(xb["f_on_head"], -xb["f_on_actin"])


# ── property (c): the abscissa/power-stroke actually ENTERS the force (the decoupling arbiter) ──────
def test_zero_abscissa_gives_zero_active_force() -> None:
    """s = 0 with the head on its anchor => zero crossbridge force: the passive baseline (no power stroke)."""
    xb = crossbridge_force(head_pos=(0, 0, 0), anchor_pos=(0, 0, 0), abscissa=0.0, walk_dir=WALK, k_xb=K_XB)
    assert np.linalg.norm(xb["f_on_actin"]) == pytest.approx(0.0, abs=1e-12)
    assert xb["load"] == pytest.approx(0.0, abs=1e-12)


def test_abscissa_enters_the_force_decoupling_arbiter() -> None:
    """d(load)/d(abscissa) = k_xb and the force GROWS with the walk — a decoupled force would not move.

    This is the exact defect gate: the original kernels advanced the abscissa but never read it in the force,
    so d(force)/d(abscissa) was 0. Here it must be k_xb (nonzero).
    """
    s1, s2 = 0.003, 0.009
    xb1 = crossbridge_force((0, 0, 0), (0, 0, 0), s1, WALK, K_XB, R0_XB)
    xb2 = crossbridge_force((0, 0, 0), (0, 0, 0), s2, WALK, K_XB, R0_XB)
    # the load rises with the walked displacement, slope exactly k_xb
    assert (xb2["load"] - xb1["load"]) / (s2 - s1) == pytest.approx(K_XB, rel=1e-10)
    # the contractile force magnitude strictly increases with abscissa (the active component is real)
    assert np.linalg.norm(xb2["f_on_actin"]) > np.linalg.norm(xb1["f_on_actin"]) > 0.0


def test_zero_walk_dir_is_passive() -> None:
    """walk_dir = 0 (unset polarity) => the power stroke cannot advance the attachment => passive head."""
    xb = crossbridge_force((0, 0, 0), (0, 0, 0), abscissa=0.05, walk_dir=(0, 0, 0), k_xb=K_XB)
    assert np.linalg.norm(xb["f_on_actin"]) == pytest.approx(0.0, abs=1e-12)
    assert xb["load"] == pytest.approx(0.0, abs=1e-12)


# ── property (b): force-velocity emerges; force -> F_stall as v -> 0 ────────────────────────────────
@pytest.mark.parametrize("kappa", [0.25, 0.5, 1e8])
def test_force_velocity_settles_at_stall(kappa: float) -> None:
    """Isometric stepping self-limits: load -> F_stall, v -> 0, strain -> F_stall/k_xb + r0_xb, |F_actin| -> F_stall."""
    traj = step_engaged_head(V0, F_STALL, kappa, K_XB, r0_xb=R0_XB, walk_dir=WALK,
                             tau=1.0e-3, n_steps=60_000)
    assert traj["velocity"][-1] == pytest.approx(0.0, abs=1e-4)          # walk has stalled
    assert traj["load"][-1] == pytest.approx(F_STALL, rel=2e-3)         # tension settled at the stall force
    assert traj["abscissa"][-1] == pytest.approx(stall_abscissa(F_STALL, K_XB, R0_XB), rel=2e-3)
    # the settled contractile force magnitude is the per-head stall force (directed, contractile)
    assert abs(traj["f_actin_along_walk"][-1]) == pytest.approx(F_STALL, rel=2e-3)
    assert traj["f_actin_along_walk"][-1] < 0.0


def test_force_free_at_start_then_builds() -> None:
    """At s = 0 the head is force-free (unloaded, v = v0); the contractile force BUILDS as it walks."""
    traj = step_engaged_head(V0, F_STALL, kappa=0.5, k_xb=K_XB, tau=1.0e-3, n_steps=60_000)
    assert traj["load"][0] == pytest.approx(0.0, abs=1e-12)
    assert traj["velocity"][0] == pytest.approx(V0, rel=1e-12)          # force-free => full unloaded speed
    assert traj["f_actin_along_walk"][0] == pytest.approx(0.0, abs=1e-12)
    assert abs(traj["f_actin_along_walk"][-1]) > abs(traj["f_actin_along_walk"][0])


def test_load_rises_monotonically_and_never_overshoots_stall() -> None:
    """The walked tension rises monotonically toward F_stall and never exceeds it (clamped Hill self-limiting)."""
    traj = step_engaged_head(V0, F_STALL, kappa=0.5, k_xb=K_XB, tau=1.0e-3, n_steps=60_000)
    assert np.all(np.diff(traj["load"]) >= -1e-12)                       # non-decreasing
    assert np.all(traj["load"] <= F_STALL + 1e-9)                        # never overshoots the stall force
    assert np.all(traj["velocity"] >= 0.0)                               # never actively lengthens


def test_stall_abscissa_matches_working_stroke_strain() -> None:
    """With r0_xb = 0 the settled abscissa equals minifilament_topology.working_stroke_strain(F_stall, k_xb)."""
    assert stall_abscissa(F_STALL, K_XB, r0_xb=0.0) == pytest.approx(working_stroke_strain(F_STALL, K_XB), rel=1e-12)


def test_load_at_stall_is_hill_consistent() -> None:
    """The settled load is exactly where the Hill inverse map gives v = 0 (F(0) = F_stall) — internal consistency."""
    # F(v=0) = F_stall for any kappa; the coupled loop must settle there
    assert float(hill_force(0.0, V0, F_STALL, kappa=0.5)) == pytest.approx(F_STALL, rel=1e-12)
