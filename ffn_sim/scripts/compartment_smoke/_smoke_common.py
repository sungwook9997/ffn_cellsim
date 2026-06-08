"""Shared plumbing for compartment ENABLED-PATH smoke harnesses (2026-06-09).

PLUMBING ONLY. Each harness in this package assembles ONE default-OFF
compartment (+ a tiny inert cortex stub) into a minimal HOOMD simulation, steps
it with the project's Leimkuhler-Matthews BAOAB integrator, asserts the enabled
path does not crash, and measures a COARSE observable. These are the FIRST
actual executions of each compartment's enabled build+step path, so they catch
crash-on-enable bugs a static read / sim.run(0) coexistence test cannot.

NOT production physics
----------------------
* NO band pass/fail, NO biological claim. The observable is a plumbing sanity
  number (bonds stayed finite, rod stayed rod-like, updater stepped), never a
  validation verdict.
* None-gated constants use the ``COMPARTMENT_ACTIVATION_DETAILED_PLANS_2026-06-09``
  PI-CANDIDATE values, used here ONLY to make the plumbing run and labelled
  ``SMOKE-ONLY`` at every call site. They are NOT ratified; the production
  enabled path still raises until PI ratifies them.
* Every figure carries a SMOKE watermark.

Drag / kT are at the physiological cytoplasm setpoint (eta=65.9 Pa.s, Hu 2024
MCF7; kT=4.28e-21 J at 310 K) per the physiological-baseline rule — the smoke
runs FROM the real operating point, not water.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Mapping

import numpy as np

import gsd.hoomd
import hoomd
import hoomd.md as md

from ffn_sim.integrator.baoab import make_baoab_updater

# --------------------------------------------------------------------------
# Physiological constants (cited; NOT smoke-invented)
# --------------------------------------------------------------------------
ETA_CYTOPLASM: float = 65.9        # Pa.s   MCF7 cytoplasm (Hu 2024; registry)
KT: float = 4.28e-21               # J      project k_B.T at 310 K

# --------------------------------------------------------------------------
# Output locations (this harness package owns a dedicated subtree so it never
# collides with h7 production artifacts)
# --------------------------------------------------------------------------
_FFN = Path(__file__).resolve().parents[2]          # ffn_sim/
_OUT = _FFN / "outputs" / "compartment_smoke"
FIG_DIR = _OUT / "figs"
JSON_DIR = _OUT / "json"


def stokes_drag(r_bead: float, eta: float = ETA_CYTOPLASM) -> float:
    """Per-bead Stokes drag gamma_b = 6*pi*eta*R [N.s/m]."""
    return 6.0 * math.pi * eta * r_bead


# --------------------------------------------------------------------------
# Tiny inert cortex stub (no bonds/angles → never collides with a compartment's
# own md.bond.Harmonic type coverage; just exercises multi-type coexistence +
# the BAOAB gamma_map coverage path)
# --------------------------------------------------------------------------
def cortex_stub_frame(
    *,
    n_beads: int = 12,
    r_shell: float = 1.0e-6,
    box_l: float = 5.0e-5,
    bead_type: str = "cortex_actin",
) -> gsd.hoomd.Frame:
    """A minimal gsd Frame of ``n_beads`` inert beads on a Fibonacci shell.

    No bonds, no angles — purely a host the compartment extends. Dense tags
    [0, n_beads). Box is large enough that the compartment fits inside.
    """
    idx = np.arange(n_beads, dtype=np.float64) + 0.5
    phi = math.pi * (3.0 - math.sqrt(5.0))
    cos_t = np.clip(1.0 - 2.0 * idx / float(n_beads), -1.0, 1.0)
    sin_t = np.sqrt(np.maximum(0.0, 1.0 - cos_t * cos_t))
    az = phi * idx
    pos = np.stack(
        [r_shell * sin_t * np.cos(az),
         r_shell * sin_t * np.sin(az),
         r_shell * cos_t], axis=1
    )
    fr = gsd.hoomd.Frame()
    fr.particles.N = int(n_beads)
    fr.particles.types = [bead_type]
    fr.particles.typeid = [0] * int(n_beads)
    fr.particles.position = pos.tolist()
    fr.particles.mass = [1.0] * int(n_beads)
    fr.particles.velocity = [[0.0, 0.0, 0.0]] * int(n_beads)
    fr.particles.image = [[0, 0, 0]] * int(n_beads)
    fr.configuration.box = [box_l, box_l, box_l, 0.0, 0.0, 0.0]
    fr.bonds.N = 0
    fr.bonds.types = []
    fr.angles.N = 0
    fr.angles.types = []
    return fr


# --------------------------------------------------------------------------
# Minimal stepping simulation (project BAOAB, faithful enabled path)
# --------------------------------------------------------------------------
def make_sim(snap, *, seed: int = 1) -> hoomd.Simulation:
    """Create a CPU Simulation with state from ``snap``."""
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=seed)
    sim.create_state_from_snapshot(snap)
    return sim


def set_integrator(sim: hoomd.Simulation, dt: float) -> md.Integrator:
    """Attach an empty-methods Integrator (BAOAB Action does the position step)."""
    ig = md.Integrator(dt=dt, methods=[])
    sim.operations.integrator = ig
    return ig


def auto_gamma_map(
    sim: hoomd.Simulation,
    *,
    radius_by_type: Mapping[str, float] | None = None,
    default_radius: float = 1.0e-7,
    eta: float = ETA_CYTOPLASM,
) -> dict[str, float]:
    """Build a gamma_map covering EVERY particle type in the state.

    Per-type Stokes drag at cytoplasm viscosity, using ``radius_by_type`` where
    given and ``default_radius`` (cortex-bead scale, 100 nm) otherwise. Auto-
    coverage guarantees the BAOAB attach never trips its missing-type guard from
    the harness side, so any failure isolates to the compartment build/step.
    """
    radius_by_type = dict(radius_by_type or {})
    gamma: dict[str, float] = {}
    for t in sim.state.particle_types:
        r = radius_by_type.get(t, default_radius)
        gamma[t] = stokes_drag(r, eta)
    return gamma


def attach_baoab(
    sim: hoomd.Simulation,
    *,
    dt: float,
    gamma_map: Mapping[str, float],
    seed: int = 2,
):
    """Wrap + append the L-M BAOAB per-step CustomUpdater. Returns the Action."""
    action, updater = make_baoab_updater(
        kT=KT, gamma=dict(gamma_map), dt=dt, seed=seed
    )
    sim.operations.updaters.append(updater)
    return action


def positions_by_tag(sim: hoomd.Simulation) -> np.ndarray:
    """(N, 3) positions ordered by particle tag (stable identity)."""
    with sim.state.cpu_local_snapshot as snap:
        tags = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position, dtype=np.float64).copy()
    out = np.empty_like(pos)
    out[tags] = pos
    return out


def all_finite(sim: hoomd.Simulation) -> bool:
    """True iff every particle position is finite (no NaN/Inf blow-up)."""
    return bool(np.all(np.isfinite(positions_by_tag(sim))))


# --------------------------------------------------------------------------
# Figure + JSON helpers
# --------------------------------------------------------------------------
def _watermark(fig) -> None:
    fig.text(
        0.5, 0.5, "SMOKE",
        fontsize=80, color="0.85", alpha=0.5,
        ha="center", va="center", rotation=30, zorder=0,
    )
    fig.text(
        0.5, 0.015,
        "PLUMBING SMOKE — not a physics claim; SMOKE-ONLY candidate constants; "
        "no band pass/fail",
        fontsize=7, color="0.4", ha="center", va="bottom",
    )


def save_fig(fig, name: str) -> Path:
    """Save a figure (with SMOKE watermark) under outputs/compartment_smoke/figs."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    _watermark(fig)
    path = FIG_DIR / (name if name.endswith(".png") else f"{name}.png")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    return path


def save_json(name: str, data: dict) -> Path:
    """Save a result dict under outputs/compartment_smoke/json."""
    JSON_DIR.mkdir(parents=True, exist_ok=True)
    path = JSON_DIR / (name if name.endswith(".json") else f"{name}.json")
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, default=str)
    return path
