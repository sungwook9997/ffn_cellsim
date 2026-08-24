r"""CPU tests for the CONFIGURABLE cortex geometry (mesh-fidelity study, 2026-07-24).

These assert the spec/arch construction + derived counts of :func:`ac.cell.assemble._cortex_region` and the
new :class:`ac.cell.assemble.CellConfig` geometry knobs — WITHOUT a CUDA build (``build_cell`` needs the A5000;
guarded out here). The two facts under test:

  1. **Default = bit-parity.** ``CellConfig`` defaults (0.5 / 3.0 / 20.0) reproduce the historical hard-coded
     ``architecture_spec.CORTEX`` geometry EXACTLY — 7 nodes/filament (494,802 nodes), 20 crosslinks/filament
     (1,413,720 crosslinks), and the module-level ``CORTEX`` constant is never mutated.
  2. **Finer mesh scales as expected.** ``cortex_seg_um`` drives ``beads_per_filament = round(L/seg)+1`` and the
     node count; ``cortex_density_per_fil`` drives ``n_xl = round(n_fil·density)``. A ~75–100 nm physiological
     mesh (memo §7) is ``seg_um ≈ 0.075–0.1 µm``.

The downstream discretization formulas mirrored here are the SoT ones in ``ff.weave.weave``:
``beads_per_filament = int(round(fs.length_um / fs.seg_um)) + 1`` and ``n_xl = int(round(n_fil·density))``.
"""

from __future__ import annotations

import dataclasses

import pytest

from aleph.components.incumbent.assemble import CellConfig, R_CORTEX_UM, _cortex_region
from aleph.laws.architecture_spec import CORTEX

N_FIL = 70686  # FULL native cortical F-actin population (the production baseline)


def _beads_per_filament(arch) -> int:
    """Mirror ff.weave.weave: beads_per_filament = round(length_um / seg_um) + 1."""
    return int(round(arch.filament.length_um / arch.filament.seg_um)) + 1


def _n_nodes(arch) -> int:
    return _beads_per_filament(arch) * arch.filament.n_filaments


def _n_xl(arch) -> int:
    """Mirror ff.weave.weave: n_xl = round(n_filaments · density_per_fil)."""
    return int(round(arch.filament.n_filaments * arch.crosslinker.density_per_fil))


# ── 1. defaults reproduce today EXACTLY (bit-parity) ───────────────────────────────────────────────


def test_cellconfig_defaults_match_hardcoded_cortex():
    """The new CellConfig knobs default to the historical architecture_spec.CORTEX values."""
    cfg = CellConfig()
    assert cfg.cortex_seg_um == CORTEX.filament.seg_um == 0.5
    assert cfg.cortex_length_um == CORTEX.filament.length_um == 3.0
    assert cfg.cortex_density_per_fil == CORTEX.crosslinker.density_per_fil == 20.0


def test_default_cortex_region_is_bit_identical_to_legacy():
    """Default _cortex_region == the OLD logic (replace only n_filaments + R_um) — additive, no drift."""
    # Reconstruct exactly what the pre-change _cortex_region produced.
    old_fil = dataclasses.replace(CORTEX.filament, n_filaments=N_FIL)
    old_arch = dataclasses.replace(CORTEX, filament=old_fil, R_um=R_CORTEX_UM)

    got = _cortex_region(N_FIL).arch
    assert got == old_arch  # frozen dataclasses compare by value → proves bit-parity


def test_default_counts_are_todays_baseline():
    """Default geometry ⇒ 7 nodes/fil, 494,802 actin nodes, 1,413,720 crosslinks."""
    arch = _cortex_region(N_FIL).arch
    assert _beads_per_filament(arch) == 7
    assert _n_nodes(arch) == 494_802
    assert _n_xl(arch) == 1_413_720
    assert arch.R_um == R_CORTEX_UM  # cortex shell sits just inside the membrane


def test_module_level_cortex_constant_not_mutated():
    """Threading a finer geometry must NOT mutate the shared CORTEX constant."""
    _ = _cortex_region(N_FIL, seg_um=0.05, length_um=1.0, density_per_fil=60.0)
    assert CORTEX.filament.seg_um == 0.5
    assert CORTEX.filament.length_um == 3.0
    assert CORTEX.crosslinker.density_per_fil == 20.0


# ── 2. finer mesh + density flow through as expected ───────────────────────────────────────────────


def test_finer_seg_scales_beads_and_nodes():
    """cortex_seg_um=0.1 µm ⇒ 31 nodes/fil ⇒ 2,191,266 actin nodes (memo §7 upper-physiological mesh)."""
    arch = _cortex_region(N_FIL, seg_um=0.1).arch
    assert arch.filament.seg_um == 0.1
    assert _beads_per_filament(arch) == 31
    assert _n_nodes(arch) == 2_191_266  # 70,686 × 31
    # density untouched ⇒ crosslink count unchanged from baseline
    assert _n_xl(arch) == 1_413_720


def test_density_flows_to_crosslink_count():
    """cortex_density_per_fil is the crosslink knob: n_xl = round(n_fil·density)."""
    arch = _cortex_region(N_FIL, density_per_fil=40.0).arch
    assert arch.crosslinker.density_per_fil == 40.0
    assert _n_xl(arch) == 2_827_440  # 70,686 × 40


def test_derived_75nm_fine_config():
    """The memo §7 DERIVED fine config (seg=0.075 µm, density=40) reproduces its stated cost.

    seg_um = mesh_target (75 nm) ⇒ beads = round(3.0/0.075)+1 = 41; density = length/mesh_target = 3.0/0.075
    = 40 holds the ~75 nm crosslink spacing. Cost: 2,898,126 nodes (~5.86× the 494,802 baseline),
    2,827,440 crosslinks (~2.0×).
    """
    arch = _cortex_region(N_FIL, seg_um=0.075, density_per_fil=40.0).arch
    assert _beads_per_filament(arch) == 41
    assert _n_nodes(arch) == 2_898_126
    assert _n_xl(arch) == 2_827_440
    assert _n_nodes(arch) / 494_802 == pytest.approx(5.857, rel=1e-3)


def test_crosslink_contour_spacing_matches_mesh_target():
    """Crosslink spacing along the contour = length/density; hitting ~75 nm ⇒ density ≈ 40 at L=3 µm."""
    length_um, mesh_target_um = 3.0, 0.075
    density = length_um / mesh_target_um
    assert density == pytest.approx(40.0)
    arch = _cortex_region(N_FIL, seg_um=mesh_target_um, length_um=length_um, density_per_fil=density).arch
    spacing_um = arch.filament.length_um / arch.crosslinker.density_per_fil
    assert spacing_um == pytest.approx(mesh_target_um)  # 75 nm crosslink spacing ≈ mesh (Chugh & Paluch 2018)
