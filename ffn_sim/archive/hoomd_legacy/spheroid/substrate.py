"""L2.6 substrate confinement: a z=0 adhesive wall → quasi-2D wetting (the experiment geometry).

The PI experiment is a spheroid spreading ON a substrate (pV4D4 / col-I coated dish); the
spread area A is a quasi-2D footprint. In free 3D the CBM grows a ball and the projected area
conflates 3D growth with 2D spread (and the convex hull spans any 3D roughness). Confining the
cells to the substrate makes the spheroid spread as a quasi-2D cap/monolayer, so the projected
area IS the spread footprint and the proliferating-rim fraction scales as the 2D edge ratio
~1/R — the clean geometric origin of the A/A0 = a + b/R + c/R² law.

Mechanism (a Morse adhesive wall, the cell-substrate analog of the cell-cell cohesion)
--------------------------------------------------------------------------------------
A planar wall at z = 0 with a ``md.external.wall.Morse`` potential gives each cell:
- a repulsive floor (cells cannot sink through the substrate), and
- an adhesion well of depth ``D_sub`` at standoff ``r0_sub = R_cell`` (a cell rests one radius
  above the dish, adhered) — the cell-substrate adhesion that drives wetting/spreading.

The substrate-adhesion vs cell-cohesion balance sets the wetting (3D cap → flat monolayer),
the platform analog of the experiment's contact angle. ``D_sub`` is the natural place to encode
the Bare / Pre / Lam4 ligand conditions (stronger substrate adhesion → more spreading); for now
it is anchored to the measured cohesion scale and flagged as the ligand-condition axis.

Pool-ready: the wall acts per particle type; ``void`` (the L2.4b parked pool) gets r_cut 0.

Sanity Gate
-----------
- Dimensional: z, r0_sub [m]; D_sub [J]; alpha_sub [1/m]; r_cut [m].
- Boundary: D_sub=0 → no substrate (free 3D, reduces to the L2.x bulk model). A cell far above
  the wall (z ≫ r0_sub + r_cut) feels nothing.
- Sign/sense: larger D_sub → flatter spheroid (more wetting, larger footprint A); larger
  cohesion → rounder (less spreading). The wall adds NO net in-plane force (xy isotropic).
- Conservation: a conservative external potential; BAOAB integrates as before; count invariant.
"""

from __future__ import annotations

from dataclasses import dataclass

import hoomd
from hoomd import md

from ffn_sim.archive.hoomd_legacy.spheroid.params import ResolvedL2

__all__ = ["ResolvedSubstrate", "resolve_substrate", "add_substrate_wall"]

# Numerical-policy: truncate the Morse wall tail this many decay lengths beyond r0_sub.
_WALL_CUTOFF_N_RANGES = 5.0


@dataclass(frozen=True)
class ResolvedSubstrate:
    """Resolved z=0 substrate-adhesion wall parameters (SI)."""

    D_sub: float        # J     substrate-adhesion well depth (cell-substrate; ligand axis)
    alpha_sub: float    # 1/m   inverse adhesive range (= cell-cell morse_alpha by default)
    r0_sub: float       # m     cell standoff = R_cell (cell rests one radius above the dish)
    r_cut: float        # m     wall potential cutoff


def resolve_substrate(
    resolved: ResolvedL2,
    *,
    adhesion_ratio: float = 1.0,
) -> ResolvedSubstrate:
    """Derive the substrate wall from the cell scales.

    Args:
        resolved: resolved CBM params (R_cell, contact_zone, D_e cohesion scale).
        adhesion_ratio: substrate adhesion as a multiple of the cell-cell cohesion energy D_e
            (the wetting knob / ligand-condition axis). 1.0 = balanced wetting; >1 spreads more.
    """
    alpha_sub = 1.0 / resolved.contact_zone_width
    return ResolvedSubstrate(
        D_sub=float(adhesion_ratio * resolved.D_e),
        alpha_sub=float(alpha_sub),
        r0_sub=float(resolved.R_cell),
        r_cut=float(resolved.R_cell + _WALL_CUTOFF_N_RANGES / alpha_sub),
    )


def add_substrate_wall(
    sim: hoomd.Simulation, sub: ResolvedSubstrate, *, cell_type: str = "cell"
) -> "md.external.wall.Morse":
    """Append a z=0 adhesive Morse wall to the simulation's integrator. Returns the wall force.

    The wall acts on ``cell_type`` particles (adhesion + floor); other types (e.g. ``void``)
    get r_cut 0 (no interaction). Must be called after ``sim.operations.integrator`` is set.
    """
    integrator = sim.operations.integrator
    if integrator is None:
        raise RuntimeError("add_substrate_wall: set sim.operations.integrator first.")
    plane = hoomd.wall.Plane(origin=(0.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0))
    wall = md.external.wall.Morse(walls=[plane])
    type_names = list(sim.state.particle_types)
    for t in type_names:
        if t == cell_type:
            wall.params[t] = dict(
                D0=sub.D_sub, alpha=sub.alpha_sub, r0=sub.r0_sub, r_cut=sub.r_cut
            )
        else:
            wall.params[t] = dict(D0=0.0, alpha=sub.alpha_sub, r0=sub.r0_sub, r_cut=0.0)
    integrator.forces.append(wall)
    return wall
