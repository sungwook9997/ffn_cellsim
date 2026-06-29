"""B1 — basal contractile mesh GEOMETRY layout (KU-3.5, PI directive (a)).

Sanity-gate checks for ``ffn_sim/cell/basal_mesh.generate_basal_mesh_layout``:
planar (in the basal band), within the footprint disk, force-free chains, long
cables aligned to the long axis, short infill isotropic — plus boundary cases.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.archive.hoomd_legacy.cell.basal_mesh import (
    basal_apparatus_report,
    basal_connectivity_report,
    basal_mesh_build_report,
    basal_mesh_reach,
    build_basal_filament_network,
    connect_basal_mesh,
    generate_basal_mesh_layout,
)
from ffn_sim.archive.hoomd_legacy.cell.basal_surface import build_flat_basal_surface
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers

ELL0 = 0.5e-6
FOOT_R = 5.0e-6
Z_BASAL = -7.0e-6
BAND = 200.0e-9

_H3_CFG = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"


def _p_xl():
    return resolve_crosslinkers(yaml.safe_load(open(_H3_CFG)), dt=1e-9)


def _layout(n=400, seed=7, **kw):
    return generate_basal_mesh_layout(
        n_filaments=n, ell0=ELL0, footprint_radius=FOOT_R, z_basal=Z_BASAL,
        band_thickness=BAND, seed=seed, **kw,
    )


# --- §3 conservation / counts ---
def test_counts_consistent():
    lay = _layout()
    nb = lay.n_beads_per_filament
    assert int(lay.positions_flat.shape[0]) == int(nb.sum())
    assert int(lay.bond_groups.shape[0]) == int((nb - 1).sum())
    assert int(lay.angle_groups.shape[0]) == int(np.clip(nb - 2, 0, None).sum())
    # filament_starts are exclusive cumulative sums.
    expect = np.zeros_like(nb)
    expect[1:] = np.cumsum(nb[:-1])
    assert np.array_equal(lay.filament_starts, expect)


def test_formin_fraction_recovered():
    lay = _layout(n=2000, formin_fraction=0.15, seed=3)
    frac = float(lay.is_formin.mean())
    assert 0.12 <= frac <= 0.18  # ~0.15 ± sampling


# --- §5 planar + within footprint + force-free ---
def test_planar_in_basal_band():
    lay = _layout()
    z = lay.positions_flat[:, 2]
    assert np.all(z <= Z_BASAL + 1e-12)
    assert np.all(z >= Z_BASAL - BAND - 1e-12)


def test_within_footprint_disk():
    lay = _layout()
    radial = np.linalg.norm(lay.centers_of_mass[:, :2], axis=1)
    assert np.all(radial <= FOOT_R + 1e-12)


def test_chains_force_free():
    lay = _layout()
    bg = lay.bond_groups
    seg = lay.positions_flat[bg[:, 0]] - lay.positions_flat[bg[:, 1]]
    ln = np.linalg.norm(seg, axis=1)
    assert float(np.max(np.abs(ln - ELL0) / ELL0)) < 1e-9


# --- §6 measurement: cable alignment + infill isotropy ---
def test_cables_aligned_infill_isotropic():
    lay = _layout(n=1500, long_axis=np.array([1.0, 0.0, 0.0]), seed=11)
    rep = basal_mesh_build_report(
        lay, ell0=ELL0, footprint_radius=FOOT_R, z_basal=Z_BASAL,
        band_thickness=BAND, long_axis=np.array([1.0, 0.0, 0.0]),
    )
    c = rep["controls"]
    # perfectly aligned cables (jitter 0) → |cos| ≈ 1.
    assert c["cables_aligned"]["ok"]
    assert c["cables_aligned"]["mean_abs_cos"] > 0.99
    # isotropic infill → mean |cos| ≈ 2/π ≈ 0.637.
    assert c["infill_isotropic"]["ok"]
    assert abs(c["infill_isotropic"]["mean_abs_cos"] - 2.0 / math.pi) < 0.08


def test_cables_follow_a_rotated_long_axis():
    axis = np.array([0.0, 1.0, 0.0])  # +y
    lay = _layout(n=1500, long_axis=axis, seed=5)
    tang = lay.tangents[lay.is_formin]
    absdot = np.abs(tang[:, :2] @ axis[:2])
    assert float(np.mean(absdot)) > 0.99


def test_report_pass_verdict():
    lay = _layout(n=1000, seed=2)
    rep = basal_mesh_build_report(
        lay, ell0=ELL0, footprint_radius=FOOT_R, z_basal=Z_BASAL,
        band_thickness=BAND,
    )
    assert rep["verdict"] == "PASS", rep["controls"]


# --- §2 boundary cases ---
@pytest.mark.parametrize("kw", [
    dict(n_filaments=0),
    dict(n_filaments=10, ell0=0.0),
    dict(n_filaments=10, footprint_radius=-1.0),
    dict(n_filaments=10, formin_fraction=1.5),
    dict(n_filaments=10, band_thickness=0.0),
])
def test_boundary_raises(kw):
    base = dict(
        n_filaments=10, ell0=ELL0, footprint_radius=FOOT_R, z_basal=Z_BASAL,
        band_thickness=BAND, seed=1,
    )
    base.update(kw)
    with pytest.raises(ValueError):
        generate_basal_mesh_layout(**base)


# --- B2 connectivity ---
def test_basal_mesh_reach_disk_formula():
    r = basal_mesh_reach(FOOT_R, 600)
    assert math.isclose(r, FOOT_R * math.sqrt(math.pi / 600), rel_tol=1e-12)


def test_connect_basal_mesh_percolates():
    lay = _layout(n=600, seed=73)
    seed = connect_basal_mesh(
        lay, _p_xl(), footprint_radius=FOOT_R, z_struct=3.3, bundle_mult=2,
        rng=np.random.default_rng(1),
    )
    rep = basal_connectivity_report(seed)
    c = rep["controls"]
    assert c["giant_fraction"]["ok"], c["giant_fraction"]
    assert c["coordination_z"]["ok"], c["coordination_z"]
    assert c["L_over_lc"]["ok"], c["L_over_lc"]
    assert rep["verdict"] == "PASS"
    # the seed actually produced crosslinks bridging DIFFERENT filaments.
    assert seed.n_xl > 0
    bf = np.asarray(seed.bridge_filaments)
    assert np.all(bf[:, 0] != bf[:, 1])


def test_connect_basal_mesh_z_in_band():
    lay = _layout(n=600, seed=73)
    seed = connect_basal_mesh(
        lay, _p_xl(), footprint_radius=FOOT_R, z_struct=3.3,
        rng=np.random.default_rng(2),
    )
    z = seed.z_struct_realised
    assert 3.0 <= z <= 3.5, f"realised z={z} out of band"
    assert seed.giant_fraction >= 0.9


# --- B3 integration: F-layer ON the S-layer, cables anchored FA→FA ---
def _apparatus(n_cables=10, n_infill=400, seed=2):
    surf = build_flat_basal_surface(
        footprint_radius=FOOT_R, z_basal=Z_BASAL, n_rings=10, n_fa=60,
        rng=np.random.default_rng(1),
    )
    return build_basal_filament_network(
        surf, ell0=ELL0, n_cables=n_cables, n_infill=n_infill,
        rng=np.random.default_rng(seed),
    )


def test_apparatus_cables_anchored_force_free_planar():
    app = _apparatus()
    rep = basal_apparatus_report(app, ell0=ELL0)
    c = rep["controls"]
    assert c["cables_anchored"]                      # both ends of every cable
    assert c["n_anchor_bonds"] == 2 * app.n_cables
    assert c["anchors_force_free"]["ok"]             # cable end sits AT the FA
    assert c["in_slab"]["ok"]                        # within the basal actin slab
    assert c["chains_force_free"]["ok"]              # per-filament exact spacing
    assert rep["verdict"] == "PASS"


def test_apparatus_cables_are_the_formin_set():
    app = _apparatus(n_cables=8)
    lay = app.layout
    # the first n_cables filaments are the cables (is_formin True).
    assert lay.is_formin[:8].all()
    assert not lay.is_formin[8:].any()


def test_apparatus_anchor_ends_at_fa_positions():
    app = _apparatus()
    pos = app.layout.positions_flat
    ab = app.anchor_bonds
    sep = np.linalg.norm(pos[ab[:, 0]] - app.fa_positions[ab[:, 1]], axis=1)
    assert float(np.max(sep)) < 1e-15   # end bead AT the FA


def test_apparatus_combined_mesh_percolates():
    app = _apparatus(n_cables=10, n_infill=500)
    seed = connect_basal_mesh(
        app.layout, _p_xl(), footprint_radius=FOOT_R, z_struct=3.3,
        rng=np.random.default_rng(4),
    )
    rep = basal_connectivity_report(seed)
    assert rep["controls"]["giant_fraction"]["ok"], rep["controls"]
    assert rep["controls"]["coordination_z"]["ok"], rep["controls"]


def test_apparatus_needs_enough_fa():
    surf = build_flat_basal_surface(
        footprint_radius=FOOT_R, z_basal=Z_BASAL, n_rings=10, n_fa=10,
        rng=np.random.default_rng(1),
    )
    with pytest.raises(ValueError):
        build_basal_filament_network(surf, ell0=ELL0, n_cables=20, n_infill=50)


# --- B4 sf_myosin_ NMII placement on the apparatus ---
def _p_sf_myosin(n_motors=40):
    import yaml
    from copy import deepcopy
    from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin
    cfg = deepcopy(yaml.safe_load(open(_H3_CFG)))
    cfg["cortex"]["myosin"]["prefix"] = "sf_myosin_"
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors
    return resolve_cortex_myosin(cfg, dt=1e-9)


def test_sf_myosin_placement_gate_pass():
    from ffn_sim.archive.hoomd_legacy.cell.basal_mesh import (
        place_sf_myosin_on_apparatus, sf_myosin_placement_report,
    )
    app = _apparatus(n_cables=12, n_infill=500)
    p_sf = _p_sf_myosin(40)
    myo = place_sf_myosin_on_apparatus(app, p_sf, rng=np.random.default_rng(5))
    rep = sf_myosin_placement_report(app, p_sf, myo)
    c = rep["controls"]
    assert c["assembled"]
    assert c["force_free"]["ok"]
    assert c["heads_near_actin"]["ok"]                 # ≥80% heads binding-eligible
    assert c["heads_near_actin"]["eligible_fraction"] >= 0.8
    assert c["sf_denylisted"]                           # all sf_myosin_ types under sf_
    assert c["distinct_from_cortex_myosin"]             # cortical active-γ untouched
    assert rep["verdict"] == "PASS"


def test_sf_myosin_requires_sf_prefix():
    import yaml
    from copy import deepcopy
    from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin
    from ffn_sim.archive.hoomd_legacy.cell.basal_mesh import place_sf_myosin_on_apparatus
    app = _apparatus()
    cfg = deepcopy(yaml.safe_load(open(_H3_CFG)))  # default prefix cortex_myosin_
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = 20
    p_cortical = resolve_cortex_myosin(cfg, dt=1e-9)
    with pytest.raises(ValueError):
        place_sf_myosin_on_apparatus(app, p_cortical)   # cortex_ prefix → refused


def test_sf_myosin_types_all_sf_prefixed():
    from ffn_sim.archive.hoomd_legacy.cortex.myosin import (
        myosin_attach_bin_names, myosin_backbone_bond_name,
        myosin_head_backbone_bond_name, myosin_particle_type_names,
    )
    p_sf = _p_sf_myosin(10)
    types = [
        myosin_backbone_bond_name(p_sf.prefix),
        myosin_head_backbone_bond_name(p_sf.prefix),
        *myosin_particle_type_names(p_sf.prefix),
        *myosin_attach_bin_names(p_sf.n_bins, p_sf.prefix),
    ]
    assert all(t.startswith("sf_myosin_") for t in types)
    assert not any(t.startswith("cortex_myosin_") for t in types)


# --- bending EI anchor (the missing SF/basal flexural term) ---
def test_basal_bending_resolve_scales_with_N():
    from ffn_sim.archive.hoomd_legacy.cell.basal_mesh import resolve_basal_bending, EI_SINGLE_ACTIN
    b10 = resolve_basal_bending(N_filaments=10, ell0=ELL0)
    b20 = resolve_basal_bending(N_filaments=20, ell0=ELL0)
    # loose bundle: EI_cable = N·EI_single.
    assert math.isclose(b10["EI_cable_N_m2"], 10 * EI_SINGLE_ACTIN, rel_tol=1e-12)
    assert math.isclose(b20["EI_cable_N_m2"], 20 * EI_SINGLE_ACTIN, rel_tol=1e-12)
    assert math.isclose(b20["EI_cable_N_m2"], 2 * b10["EI_cable_N_m2"], rel_tol=1e-12)
    # tight bundle scales N².
    bt = resolve_basal_bending(N_filaments=10, ell0=ELL0, bundle_coupling="tight")
    assert math.isclose(bt["EI_cable_N_m2"], 100 * EI_SINGLE_ACTIN, rel_tol=1e-12)


def test_basal_bending_report_force_free_and_bend_before_stretch():
    app = _apparatus(n_cables=12, n_infill=500)
    rep = basal_bending_report(app, N_filaments=20, ell0=ELL0)
    c = rep["controls"]
    assert c["angle_groups_present"]
    assert c["n_angle_triplets"] > 0
    assert c["force_free"]["ok"]                         # straight chains → θ≈π
    bvs = c["bend_vs_stretch"]
    # the AFINES condition: stretch ≫ bending → filaments bend before they stretch.
    assert bvs["bend_before_stretch"]
    assert bvs["ratio_ka_over_kbend"] > 1e3              # large headroom (~1.5e5)
    assert rep["verdict"] == "PASS"


def test_bend_vs_stretch_ratio_is_N_independent():
    # both ka (=N·EA/ℓ0) and kbend (=N·EI/ℓ0³) scale with N → ratio cancels N.
    app = _apparatus()
    r7 = basal_bending_report(app, N_filaments=7, ell0=ELL0)["controls"]["bend_vs_stretch"]["ratio_ka_over_kbend"]
    r30 = basal_bending_report(app, N_filaments=30, ell0=ELL0)["controls"]["bend_vs_stretch"]["ratio_ka_over_kbend"]
    assert math.isclose(r7, r30, rel_tol=1e-9)


from ffn_sim.archive.hoomd_legacy.cell.basal_mesh import basal_bending_report  # noqa: E402


def test_jitter_spreads_cables():
    lay0 = _layout(n=2000, long_axis_jitter_deg=0.0, seed=9)
    layj = _layout(n=2000, long_axis_jitter_deg=20.0, seed=9)
    a0 = np.abs(lay0.tangents[lay0.is_formin][:, 0])
    aj = np.abs(layj.tangents[layj.is_formin][:, 0])
    # jitter lowers the mean alignment to the +x axis.
    assert float(np.mean(aj)) < float(np.mean(a0))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
