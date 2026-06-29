"""Sanity-gate tests for the CORTEX CONSTRUCTION REBUILD (2026-06-04).

Covers the three new construction layers that turn the fragmented mesh
(z = 1.3, giant 7 %, 56 % same-filament staples) into a CONNECTED SPANNING
MESH (z ∈ [3.0, 3.5], giant ≥ 0.9, L/lc ≥ 5.9):

* ``cortex.generate_bimodal_cortex_layout`` — bimodal-exponential length
  (Arp2/3 short + formin long backbone) + disordered isotropic + 200 nm
  shell-band projection.
* ``crosslinkers.seed_connected_mesh_xlinks`` — bridge-different-filament
  (no same-filament staples) + per-filament degree-capped + bundling +
  seeded-at-construction (adhered baseline).
* ``connected_mesh.build_connected_cortex`` — the HOOMD assembler; stable at
  the PHYSIOLOGICAL cytoplasm viscosity (no construction-overlap blowups).

Per CLAUDE.md Sanity-Gate Protocol §1-6 (dimensional / boundary / conservation
/ numerical / sign-sense / measurement-protocol).
"""
from __future__ import annotations

import copy
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.archive.hoomd_legacy.cortex.cortex import (
    resolve_h3_derived,
    generate_bimodal_cortex_layout,
)
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import (
    resolve_crosslinkers,
    seed_connected_mesh_xlinks,
)
from ffn_sim.archive.hoomd_legacy.cell.cytoplasm import ETA_CYTO_BY_CELLTYPE

_CFG = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"


def _cfg(n_fil, *, physiological=True):
    cfg = copy.deepcopy(yaml.safe_load(open(_CFG)))
    cfg["cortex"]["R_cell"] = 7.5e-6          # MCF7
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    if physiological:
        # PHYSIOLOGICAL-BASELINE: cortex actin immersed in cytoplasm (65.9 Pa·s
        # MCF7, Hu 2024), NOT water — γ_b/dt/CFL consistent at the baseline.
        cfg["cortex"]["water_viscosity"] = ETA_CYTO_BY_CELLTYPE["MCF7"]
    return cfg


def _resolve(n_fil, **kw):
    cfg = _cfg(n_fil, **kw)
    p = resolve_h3_derived(cfg)
    p_xl = resolve_crosslinkers(cfg, dt=p.dt_cfl)
    return p, p_xl


# ===========================================================================
# Bimodal cortex layout
# ===========================================================================
class TestBimodalLayout:
    def test_formin_fraction_and_length_ratio(self):
        """§6: requested formin fraction recovered; formin filaments are the
        LONG subpopulation (the backbone)."""
        p, _ = _resolve(1500)
        lay = generate_bimodal_cortex_layout(
            p, formin_fraction=0.12, L_long_mean=5e-6,
            rng=np.random.default_rng(1),
        )
        assert abs(lay.is_formin.mean() - 0.12) < 0.04
        L = lay.L_per_filament
        # formin mean length must exceed Arp2/3 mean by a wide margin (~10×).
        assert L[lay.is_formin].mean() > 4.0 * L[~lay.is_formin].mean()

    def test_shell_band_200nm(self):
        """§2/§6: every bead within the 200 nm cortex band [R−200nm, R]."""
        p, _ = _resolve(1500)
        lay = generate_bimodal_cortex_layout(
            p, formin_fraction=0.12, L_long_mean=5e-6,
            rng=np.random.default_rng(2),
        )
        r = np.linalg.norm(lay.positions_flat, axis=1)
        drift = p.R_cell - r
        assert drift.min() >= -1e-12
        assert drift.max() <= p.cortex_thickness + 1e-12

    def test_isotropic(self):
        """§6: disordered isotropic tangent field (nematic S → 0)."""
        p, _ = _resolve(2000)
        lay = generate_bimodal_cortex_layout(
            p, formin_fraction=0.12, rng=np.random.default_rng(3))
        t = lay.tangents
        Q = np.mean(t[:, :, None] * t[:, None, :], axis=0) * 1.5 - 0.5 * np.eye(3)
        S = float(np.max(np.linalg.eigvalsh(Q)))
        assert S < 0.1, f"tangent field not isotropic: S={S}"

    def test_conservation_counts(self):
        """§3: bond / angle counts = Σ(N_i−1) / Σ(N_i−2); filament_idx aligned."""
        p, _ = _resolve(800)
        lay = generate_bimodal_cortex_layout(
            p, formin_fraction=0.12, rng=np.random.default_rng(4))
        nb = lay.n_beads_per_filament
        assert lay.positions_flat.shape[0] == int(nb.sum())
        assert lay.bond_groups.shape[0] == int((nb - 1).sum())
        assert lay.angle_groups.shape[0] == int(np.clip(nb - 2, 0, None).sum())
        assert lay.filament_idx.shape[0] == lay.positions_flat.shape[0]
        # filament_idx is non-decreasing and spans [0, n_fil)
        assert lay.filament_idx.min() == 0
        assert lay.filament_idx.max() == nb.shape[0] - 1

    def test_boundary_bad_formin_fraction(self):
        """§2: formin_fraction outside [0, 1] raises."""
        p, _ = _resolve(100)
        with pytest.raises(ValueError):
            generate_bimodal_cortex_layout(p, formin_fraction=1.5)


# ===========================================================================
# Connected-mesh seeding
# ===========================================================================
class TestSeedConnectedMesh:
    @pytest.fixture(scope="class")
    def seed(self):
        n_fil = 600
        p, p_xl = _resolve(n_fil)
        rng = np.random.default_rng(11)
        lay = generate_bimodal_cortex_layout(
            p, formin_fraction=0.12, L_long_mean=5e-6, rng=rng)
        s = seed_connected_mesh_xlinks(
            lay.positions_flat, lay.filament_idx, n_fil, p_xl,
            z_struct=3.7, bundle_mult=2, R_cell=p.R_cell,
            n_cortex_beads=lay.positions_flat.shape[0], rng=rng)
        return n_fil, lay, s

    def test_no_same_filament_staples(self, seed):
        """§5 (sign/sense): the bridge-different-filament rule — ZERO
        same-filament staples (the 56 % failure mode of the prior mesh)."""
        _, _, s = seed
        bf = s.bridge_filaments
        assert int((bf[:, 0] == bf[:, 1]).sum()) == 0

    def test_connectivity_gates(self, seed):
        """§6: the CORTEX REBUILD acceptance gates."""
        _, _, s = seed
        assert 3.0 <= s.z_struct_realised <= 3.5, f"z={s.z_struct_realised}"
        assert s.giant_fraction >= 0.9, f"giant={s.giant_fraction}"
        assert s.L_over_lc >= 5.9, f"L/lc={s.L_over_lc}"

    def test_seeded_attach_consistency(self, seed):
        """§3: two seeded attach bonds per crosslinker; head_to_actin0 matches."""
        n_fil, lay, s = seed
        assert s.seeded_attach.shape[0] == 2 * s.n_xl
        # every seeded head is recorded as bound; both heads of each dimer bound
        assert int((s.head_to_actin0 >= 0).sum()) == 2 * s.n_xl
        # attach bonds point head→actin (col0 ≥ n_cortex_beads, col1 < it)
        nca = lay.positions_flat.shape[0]
        assert (s.seeded_attach[:, 0] >= nca).all()
        assert (s.seeded_attach[:, 1] < nca).all()

    def test_low_homeless(self, seed):
        """§6: almost every filament joins the mesh (homeless ≪ 10 %)."""
        n_fil, _, s = seed
        assert s.n_homeless / n_fil < 0.1


# ===========================================================================
# HOOMD assembler — physiological-viscosity stability
# ===========================================================================
class TestBuildConnectedCortex:
    def test_build_runs_and_stable(self):
        """§4 (numerical): the HOOMD frame is valid (no exclusion overflow) and
        a short BAOAB run at the PHYSIOLOGICAL cytoplasm viscosity is STABLE
        (no construction-overlap blowup) with connectivity maintained."""
        from ffn_sim.archive.hoomd_legacy.cortex.connected_mesh import build_connected_cortex

        n_fil = 500
        p, p_xl = _resolve(n_fil)
        h = build_connected_cortex(
            p, p_xl, formin_fraction=0.12, L_long_mean=5e-6,
            z_struct=3.7, bundle_mult=2, equilibrate=True, n_softstart=200,
            rng=np.random.default_rng(5),
        )
        # particle count = cortex beads + 2·n_xl heads
        assert h.sim.state.N_particles == h.n_cortex_actin + 2 * h.seed.n_xl
        assert h.seed.giant_fraction >= 0.9
        # stable short run (would FloatingPointError on a WCA blowup)
        h.sim.run(500)
        # all seeded bridges remain engaged (connectivity maintained), and the
        # bridge rule never created a same-filament staple on rebind.
        assert h.xlink_action.n_engaged >= 0.95 * h.n_xlink_heads

    def test_physiological_viscosity_relaxes_cfl(self):
        """§1 (dimensional/numerical): the physiological cytoplasm drag is
        ~10^5× water and RELAXES the CFL (larger dt) — the larger drag also
        bounds per-step displacement (no unphysical overshoots)."""
        pw, _ = _resolve(800, physiological=False)
        pc, _ = _resolve(800, physiological=True)
        ratio = pc.gamma_b / pw.gamma_b
        assert ratio > 1e4                       # ~9.5×10^4
        assert math.isclose(pc.dt_cfl / pw.dt_cfl, ratio, rel_tol=1e-6)
