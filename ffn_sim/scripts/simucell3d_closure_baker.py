#!/usr/bin/env python3
"""ffn_cellsim -> SimuCell3D closure-baker (pipeline P0).

Bakes ffn-derived single-cell mechanical *closures* into a SimuCell3D parameter
XML + an MCF7-sized seed mesh, so the SimuCell3D deformable-cell **chassis** runs
an MCF7 spheroid driven by ffn closures. This is the offline foundry->engine
coupling: ffn (Python/HOOMD, the mechanistic foundry) computes the reduced laws;
SimuCell3D (C++ DCM, the morphology engine) consumes them. The Python binding of
SimuCell3D is one-shot (no live force injection), so coupling is by baked values,
which is exactly correct under the cortex/tissue timescale separation.

P0 scope = the mature, FORK-FREE closures that map onto native SimuCell3D params:
  gamma  cortical tension  -> face_type surface_tension   (NATIVE, 1:1)
  omega  cohesion          -> face_type adherence_strength (native contact param)
  growth/division          -> cell growth_rate + division_volume (native)
The C++ plugin forks come later: substrate adhesion (P1, the PI spreading expt),
lamellipodium traction (P2), active anisotropic cortical stress (P3, Gate-B).

Closure provenance (ffn module audit, 2026-06-09):
  gamma = Hosseini 2020 AdvSci MCF7 rounded cortical tension 0.27 mN/m
          (ffn cortex/cortical_tension.measure_cortical_tension; literature anchor
           used now, ffn-derived active channel pending Gate-B)
  omega = E-cadherin catch-bond ensemble, de-adhesion peak ~6.5 nN/contact
          (Iturri 2020; ffn validation/cadherin_sliding_rebinding + spheroid/cadherin_bonds)
  R_cell = 7.5 um  (MCF7, Wagner 2011)
  K      = 2.5e3 Pa  cytoplasm/volume incompressibility (SimuCell3D cell default; MCF7-plausible)
  cycle  = 30 h MCF7 doubling (ffn spheroid/proliferation) -> growth ACCELERATED for tractability

All SI units. SimuCell3D unit tags: length [m], surface_tension [N/m]=[M/T^2],
adherence/repulsion [M/(T^2 L^2)], bulk_modulus [Pa]=[M/(L T^2)], growth [m^3/s].

Usage:
    python simucell3d_closure_baker.py [--sc3d-root /path] [--seed-mesh sphere.vtk]
Writes  <root>/data/input_meshes/mcf7_sphere.vtk  and  <root>/mcf7_p0.xml .
Run it:  cd <root>/build && ./simucell3d ../mcf7_p0.xml
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

# --------------------------------------------------------------------------- #
# ffn-derived MCF7 closures (single-cell -> chassis), with provenance above.
# --------------------------------------------------------------------------- #
R_CELL = 7.5e-6                                   # m, MCF7 radius (Wagner 2011)
V_CELL = 4.0 / 3.0 * math.pi * R_CELL**3          # m^3, ~1.767e-15
GAMMA_CORTICAL = 0.27e-3                           # N/m, Hosseini 2020 MCF7 (0.27 mN/m)
COHESION_PEAK_N = 6.5e-9                           # N/contact, cadherin de-adhesion (Iturri 2020)
BULK_MODULUS = 2.5e3                               # Pa, cell incompressibility
CYCLE_TIME_S = 30 * 3600.0                         # s, MCF7 doubling (physiological reference)

# --- closure -> SimuCell3D param mapping ----------------------------------- #
# surface tension is a 1:1 native map.
SURFACE_TENSION = GAMMA_CORTICAL                   # [N/m]
# adherence_strength [M/(T^2 L^2)] is a force-density, NOT the raw nN; the contact
# model multiplies it by geometric (area/length) factors. The principled
# ffn(6.5 nN) -> adherence calibration needs the contact-model force law (P0.1
# refinement). For P0 we set cohesion at the repulsion scale so daughters cohere
# into a spheroid (repulsion default 1e8); flagged for calibration.
ADHERENCE_STRENGTH = 1.0e8                         # [M/(T^2 L^2)]  (P0 cohesive; calibrate vs 6.5 nN)
REPULSION_STRENGTH = 1.0e8                         # [M/(T^2 L^2)]  excluded volume (SimuCell3D default)
DIVISION_VOLUME = 2.0 * V_CELL                     # m^3, divide at doubled volume -> two MCF7 cells
MIN_VOL = 0.1 * V_CELL                             # m^3, delete shrunken cells

# Temporal params: ACCELERATED for tractability (cannot simulate 30 h at dt=1e-7 s;
# mechanics equilibrate in ~ms while growth is the slow variable -> the standard DCM
# acceleration. Physiological anchoring lives in gamma/omega/size, not wall-clock growth.)
GROWTH_RATE = 1.0e-11                              # m^3/s (accelerated; ratio to relaxation is what matters)
GROWTH_RATE_STD = 1.0e-12
DAMPING = 3.0e-9                                   # [M/T] overdamped drag (accelerated; ffn eta=65.9 Pa.s maps here, P0.1)
TIME_STEP = 1.0e-7                                 # s
SIM_DURATION = 8.0e-4                              # s  (~4-5 division generations -> small spheroid)
SAMPLING_PERIOD = 1.0e-5                           # s  (~80 morphology frames)
MIN_EDGE_LENGTH = 7.5e-7                           # m  (~10% of cell radius)
CONTACT_CUTOFF = 5.0e-7                            # m


def scale_seed_mesh(seed_path: Path, out_path: Path, target_radius: float) -> dict:
    """Recenter a VTK-legacy sphere mesh to the origin and scale it to ``target_radius``.

    Only the POINTS block is rewritten; connectivity is preserved verbatim.
    Returns geometry diagnostics.
    """
    lines = seed_path.read_text().split("\n")
    i = 0
    while not lines[i].startswith("POINTS"):
        i += 1
    n_pts = int(lines[i].split()[1])
    header_end = i  # POINTS line index
    # gather coordinate tokens until we have 3*n
    vals: list[float] = []
    j = i + 1
    while len(vals) < 3 * n_pts:
        toks = lines[j].split()
        vals.extend(float(t) for t in toks)
        j += 1
    body_start = j  # first line after the point coordinates

    pts = [(vals[3 * k], vals[3 * k + 1], vals[3 * k + 2]) for k in range(n_pts)]
    cx = sum(p[0] for p in pts) / n_pts
    cy = sum(p[1] for p in pts) / n_pts
    cz = sum(p[2] for p in pts) / n_pts
    r_mean = sum(math.dist(p, (cx, cy, cz)) for p in pts) / n_pts
    factor = target_radius / r_mean

    scaled = [
        ((p[0] - cx) * factor, (p[1] - cy) * factor, (p[2] - cz) * factor) for p in pts
    ]
    # rebuild: header (through POINTS line) + scaled coords (3 per line) + rest
    out: list[str] = lines[: header_end + 1]
    for p in scaled:
        out.append(f"{p[0]:.6e} {p[1]:.6e} {p[2]:.6e}")
    out.extend(lines[body_start:])
    out_path.write_text("\n".join(out))
    return {
        "n_pts": n_pts,
        "seed_radius_um": r_mean * 1e6,
        "scale_factor": factor,
        "target_radius_um": target_radius * 1e6,
    }


def _face_type(fid: int, name: str) -> str:
    return f"""            <face_type>
                <global_face_id>{fid}</global_face_id>
                <face_type_name>{name}</face_type_name>
                <adherence_strength>{ADHERENCE_STRENGTH:.6e}</adherence_strength>
                <repulsion_strength>{REPULSION_STRENGTH:.6e}</repulsion_strength>
                <surface_tension>{SURFACE_TENSION:.6e}</surface_tension>
                <bending_modulus>0</bending_modulus>
            </face_type>"""


def build_xml(mesh_rel: str, out_folder: str) -> str:
    faces = "\n".join(
        _face_type(k, name) for k, name in enumerate(("apical", "lateral", "basal"))
    )
    return f"""<?xml version="1.0"?>
<!-- MCF7 P0 spheroid, baked from ffn closures by simucell3d_closure_baker.py.
     Requires compile flags: POLARIZATION_MODE_INDEX 1, CONTACT_MODEL_INDEX 0, DYNAMIC_MODEL_INDEX 1 (overdamped). -->
<numerical_parameters>
    <input_mesh_file_path>{mesh_rel}</input_mesh_file_path>
    <output_mesh_folder_path>{out_folder}</output_mesh_folder_path>
    <perform_initial_triangulation>1</perform_initial_triangulation>
    <enable_edge_swap_operation>0</enable_edge_swap_operation>
    <damping_coefficient>{DAMPING:.6e}</damping_coefficient>
    <simulation_duration>{SIM_DURATION:.6e}</simulation_duration>
    <sampling_period>{SAMPLING_PERIOD:.6e}</sampling_period>
    <time_step>{TIME_STEP:.6e}</time_step>
    <min_edge_length>{MIN_EDGE_LENGTH:.6e}</min_edge_length>
    <contact_cutoff_adhesion>{CONTACT_CUTOFF:.6e}</contact_cutoff_adhesion>
    <contact_cutoff_repulsion>{CONTACT_CUTOFF:.6e}</contact_cutoff_repulsion>
</numerical_parameters>

<cell_types>
    <cell_type>
        <cell_type_name>MCF7</cell_type_name>
        <global_cell_id>0</global_cell_id>
        <cell_mass_density>1.0e3</cell_mass_density>
        <cell_bulk_modulus>{BULK_MODULUS:.6e}</cell_bulk_modulus>
        <max_inner_pressure>INF</max_inner_pressure>
        <avg_growth_rate>{GROWTH_RATE:.6e}</avg_growth_rate>
        <std_growth_rate>{GROWTH_RATE_STD:.6e}</std_growth_rate>
        <target_isoperimetric_ratio>250</target_isoperimetric_ratio>
        <area_elasticity_modulus>0</area_elasticity_modulus>
        <angle_regularization_factor>0</angle_regularization_factor>
        <avg_division_volume>{DIVISION_VOLUME:.6e}</avg_division_volume>
        <std_division_volume>{0.1 * DIVISION_VOLUME:.6e}</std_division_volume>
        <surface_coupling_max_curvature>5e6</surface_coupling_max_curvature>
        <min_vol>{MIN_VOL:.6e}</min_vol>
        <face_types>
{faces}
        </face_types>
    </cell_type>
</cell_types>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description="Bake ffn closures -> SimuCell3D MCF7 P0 config")
    ap.add_argument("--sc3d-root", default="/Users/sw1/simucell3d_ref", type=Path)
    ap.add_argument("--seed-mesh", default="sphere.vtk")
    args = ap.parse_args()

    root: Path = args.sc3d_root
    seed = root / "data" / "input_meshes" / args.seed_mesh
    mesh_out = root / "data" / "input_meshes" / "mcf7_sphere.vtk"
    xml_out = root / "mcf7_p0.xml"

    geo = scale_seed_mesh(seed, mesh_out, R_CELL)
    xml = build_xml("/data/input_meshes/mcf7_sphere.vtk", "../simulation_results/mcf7_p0")
    xml_out.write_text(xml)

    print("=== ffn -> SimuCell3D closure-baker (P0) ===")
    print(f"  seed mesh        : {seed.name}  ({geo['n_pts']} pts, r={geo['seed_radius_um']:.2f} um)")
    print(f"  scaled -> MCF7    : {mesh_out}  (x{geo['scale_factor']:.3f} -> r={geo['target_radius_um']:.2f} um)")
    print(f"  XML              : {xml_out}")
    print("  --- baked closures ---")
    print(f"  gamma (surf.tens) : {SURFACE_TENSION:.3e} N/m   (Hosseini MCF7 0.27 mN/m)")
    print(f"  omega (adherence) : {ADHERENCE_STRENGTH:.3e}     (P0 cohesive; calibrate vs 6.5 nN)")
    print(f"  K (bulk modulus)  : {BULK_MODULUS:.3e} Pa")
    print(f"  V0 / V_div        : {V_CELL:.3e} / {DIVISION_VOLUME:.3e} m^3")
    print(f"  growth (accel.)   : {GROWTH_RATE:.3e} m^3/s   (phys. cycle {CYCLE_TIME_S:.1e}s acc.)")
    print(f"  sim duration      : {SIM_DURATION:.1e} s @ dt {TIME_STEP:.1e}  (~{int(SIM_DURATION/TIME_STEP)} steps)")
    print(f"\nRun:  cd {root}/build && ./simucell3d ../mcf7_p0.xml")


if __name__ == "__main__":
    main()
