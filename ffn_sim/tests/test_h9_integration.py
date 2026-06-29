"""H.9 nucleus + KU-3.1 enclosed-volume INTEGRATION smoke tests.

Covers the additive, default-off wiring of the H.9 nucleus
(:class:`ffn_sim.archive.hoomd_legacy.cell.nucleus.NucleusConfinement`) and the KU-3.1
enclosed-volume pressure (:class:`ffn_sim.archive.hoomd_legacy.cortex.enclosed_volume.EnclosedVolumePressure`)
into :func:`ffn_sim.archive.hoomd_legacy.cell.cell.build_cortex_full_simulation`.

Asserted here (the SA-1 integration contract):

- **DEFAULT-OFF / bit-for-bit identical**: with ``p_nucleus=None`` AND
  ``p_enclosed_volume=None`` the builder's snapshot (particle count, types,
  positions, bonds, angles) is byte-identical to the bare-cortex builder, and
  no nucleus/enclosed-volume force is attached (``TestDefaultOffBitIdentical``).
- **tag-APPEND invariant**: the ``nucleus_bead`` cloud is appended LAST, at tags
  ``[N_cortex_pre, N_cortex_pre + n_beads)``; nothing already present is
  re-tagged (``TestTagAppend``).
- **KU-3.1 SMOKE**: a SMALL cell (n_fil≈200) with BOTH nucleus ON and
  enclosed-volume ON builds, runs a few hundred BAOAB steps crash-free (finite
  positions), and the runtime enclosed-volume pressure is non-trivial with the
  correct Young-Laplace sign/scale (~2γ/R) under a controlled compression
  (``TestKU31Smoke``).
- **gamma_map coverage**: the appended ``nucleus_bead`` type is registered in the
  BAOAB gamma_map (else BAOAB HALTS) — verified implicitly by the crash-free run
  and explicitly via the BAOAB action's gamma dict (``TestGammaMapCoverage``).

These are SMOKE / integration checks. The nucleus LAW's shape (strain-stiffening,
sign, net-force, CFL) is unit-tested in ``test_nucleus.py``; the enclosed-volume
Young-Laplace recovery oracle is in ``test_enclosed_volume.py``. This file does NOT
claim the composite KU-3.1 / KU-3.5 gate PASSES — it asserts the wiring is correct
and runs, and surfaces the multi-timescale CFL finding (see module note below).

Multi-timescale CFL note (PI-gate flag, NOT a contract change here)
------------------------------------------------------------------
At the small-cortex ``dt_cfl`` (≈1.3e-8 s) the nucleus stiffest branch
``k_hi = k_chrom + k_lamin`` (∝ 1/n_beads) is stiffer than the cortex per-bead
timescale, so the nucleus CFL gate (``dt ≤ 0.1·γ_nuc/k_hi``) only passes for a
SUFFICIENTLY DISCRETIZED nucleus (here n_beads=300 → ~1.9× CFL margin). The
production NATIVE path needs the B1 global-dt reconciliation
(``reconcile_dt=True``) extended to the nucleus term — flagged to the Lead /PI.
The smoke deliberately uses a CFL-safe discretization so the default strict CFL
gate (parity with erm/enclosed_volume) is honoured, NOT loosened.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import math

import numpy as np
import pytest
import yaml

import hoomd

from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.enclosed_volume import (
    resolve_enclosed_volume,
    young_laplace_pressure,
)
from ffn_sim.archive.hoomd_legacy.cell.nucleus import resolve_nucleus
from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)

# A CFL-safe nucleus discretization at the small-cortex dt_cfl (see module
# note): k_hi ∝ 1/n_beads, so 300 beads keeps 0.1·γ_nuc/k_hi above dt_cfl with
# ~1.9× margin. R_nuc = 0.25·R_cell (a typical nucleus:cell radius ratio).
_N_NUC_BEADS = 300
_R_NUC_FRACTION = 0.25


def _load_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _demo_cfg(n_filaments: int = 200) -> dict:
    cfg = deepcopy(_load_cfg())
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    return cfg


@pytest.fixture(scope="module")
def p_cortex():
    return resolve_h3_derived(_demo_cfg())


@pytest.fixture(scope="module")
def p_ev(p_cortex):
    return resolve_enclosed_volume({}, R_cell=p_cortex.R_cell)


@pytest.fixture(scope="module")
def p_nuc(p_cortex):
    return resolve_nucleus(
        {}, R_nuc=_R_NUC_FRACTION * p_cortex.R_cell, n_beads=_N_NUC_BEADS
    )


def _snapshot_arrays(sim):
    """Return a hashable summary of the construction snapshot for equality."""
    snap = sim.state.get_snapshot()
    return {
        "N": int(snap.particles.N),
        "types": tuple(snap.particles.types),
        "typeid": np.asarray(snap.particles.typeid).copy(),
        "position": np.asarray(snap.particles.position, dtype=np.float64).copy(),
        "bond_N": int(snap.bonds.N),
        "bond_types": tuple(snap.bonds.types),
        "bond_group": np.asarray(snap.bonds.group).copy(),
        "angle_N": int(snap.angles.N),
        "angle_group": np.asarray(snap.angles.group).copy(),
    }


# ---------------------------------------------------------------------------
# DEFAULT-OFF — both compartments None ⇒ bit-for-bit identical builder
# ---------------------------------------------------------------------------
class TestDefaultOffBitIdentical:
    def test_nucleus_and_ev_off_snapshot_identical(self, p_cortex):
        """p_nucleus=None AND p_enclosed_volume=None ⇒ construction snapshot
        is byte-identical to the bare-cortex full-builder (additive invariant).
        """
        h_bare = build_cortex_full_simulation(p_cortex, with_baoab=True)
        h_off = build_cortex_full_simulation(
            p_cortex, p_nucleus=None, p_enclosed_volume=None, with_baoab=True
        )
        a = _snapshot_arrays(h_bare["sim"])
        b = _snapshot_arrays(h_off["sim"])
        assert a["N"] == b["N"]
        assert a["types"] == b["types"]
        np.testing.assert_array_equal(a["typeid"], b["typeid"])
        np.testing.assert_array_equal(a["position"], b["position"])
        assert a["bond_N"] == b["bond_N"]
        assert a["bond_types"] == b["bond_types"]
        np.testing.assert_array_equal(a["bond_group"], b["bond_group"])
        assert a["angle_N"] == b["angle_N"]
        np.testing.assert_array_equal(a["angle_group"], b["angle_group"])

    def test_off_handles_are_none(self, p_cortex):
        """Off-path handles for nucleus / enclosed-volume are None / 0."""
        h = build_cortex_full_simulation(p_cortex, with_baoab=True)
        assert h["nucleus_force"] is None
        assert h["nucleus_integration"] is None
        assert h["n_nucleus_beads"] == 0
        assert h["enclosed_volume_force"] is None

    def test_no_nucleus_bead_type_when_off(self, p_cortex):
        """The nucleus_bead particle type is absent when p_nucleus is None."""
        h = build_cortex_full_simulation(p_cortex, with_baoab=True)
        snap = h["sim"].state.get_snapshot()
        assert "nucleus_bead" not in list(snap.particles.types)


# ---------------------------------------------------------------------------
# tag-APPEND invariant — nucleus block is appended LAST, nothing re-tagged
# ---------------------------------------------------------------------------
class TestTagAppend:
    def test_nucleus_appended_after_cortex(self, p_cortex, p_nuc):
        """nucleus_bead tags occupy [N_cortex, N_cortex + n_beads); the cortex
        tag range [0, N_cortex) is byte-identical to the nucleus-off build."""
        N_cortex = p_cortex.n_filaments * p_cortex.beads_per_filament
        h_off = build_cortex_full_simulation(p_cortex, with_baoab=True)
        h_on = build_cortex_full_simulation(
            p_cortex, p_nucleus=p_nuc, with_baoab=True
        )
        ni = h_on["nucleus_integration"]
        assert ni is not None
        assert ni.nucleus_tag_start == N_cortex
        assert ni.n_beads == _N_NUC_BEADS
        assert h_on["sim"].state.N_particles == N_cortex + _N_NUC_BEADS

        # The pre-existing cortex particle block must be UNCHANGED (tag-APPEND:
        # appending must not perturb earlier tags' positions/types).
        off = _snapshot_arrays(h_off["sim"])
        on = _snapshot_arrays(h_on["sim"])
        np.testing.assert_array_equal(
            on["position"][:N_cortex], off["position"][:N_cortex]
        )
        np.testing.assert_array_equal(
            on["typeid"][:N_cortex], off["typeid"][:N_cortex]
        )
        # Bonds/angles reference only cortex tags and must be unchanged.
        assert on["bond_N"] == off["bond_N"]
        np.testing.assert_array_equal(on["bond_group"], off["bond_group"])
        assert on["angle_N"] == off["angle_N"]
        np.testing.assert_array_equal(on["angle_group"], off["angle_group"])

    def test_nucleus_beads_seeded_near_centroid(self, p_cortex, p_nuc):
        """Appended nucleus beads sit near the cortex centroid (r ≤ R_nuc),
        inside the cortex shell (R_cell), as the confinement expects."""
        N_cortex = p_cortex.n_filaments * p_cortex.beads_per_filament
        h = build_cortex_full_simulation(
            p_cortex, p_nucleus=p_nuc, with_baoab=True
        )
        snap = h["sim"].state.get_snapshot()
        pos = np.asarray(snap.particles.position, dtype=np.float64)
        nuc = pos[N_cortex:]
        c = nuc.mean(axis=0)
        r = np.linalg.norm(nuc - c, axis=1)
        # Seeded as a fill cloud of nominal radius R_nuc. The outermost bead
        # sits at R_nuc from the SEED centroid; the recomputed centroid of the
        # finite Fibonacci-fill cloud is marginally offset, so r.max() about
        # the recomputed centroid can exceed R_nuc by a small discretization
        # residual (~1%). The physically meaningful bound is that the whole
        # cloud stays well inside the cortex shell (R_cell).
        assert r.max() <= p_nuc.R_nuc * 1.05
        assert r.max() < 0.5 * p_cortex.R_cell


# ---------------------------------------------------------------------------
# gamma_map coverage — appended type registered (else BAOAB HALTS)
# ---------------------------------------------------------------------------
class TestGammaMapCoverage:
    def test_nucleus_in_baoab_gamma_map(self, p_cortex, p_nuc):
        """The BAOAB action's gamma dict carries nucleus_bead (γ > 0 finite);
        absence would HALT BAOAB at attach. We read it off the action."""
        h = build_cortex_full_simulation(
            p_cortex, p_nucleus=p_nuc, with_baoab=True
        )
        action = h["baoab_action"]
        gamma = None
        for attr in ("gamma", "_gamma", "gamma_map", "_gamma_map"):
            if hasattr(action, attr):
                cand = getattr(action, attr)
                if isinstance(cand, dict):
                    gamma = cand
                    break
        assert gamma is not None, (
            "could not locate the BAOAB gamma map on the action; "
            f"action attrs: {dir(action)}"
        )
        assert "nucleus_bead" in gamma
        g = float(gamma["nucleus_bead"])
        assert math.isfinite(g) and g > 0.0
        # The host default = cortex gamma_b.
        assert math.isclose(g, p_cortex.gamma_b, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# KU-3.1 SMOKE — nucleus ON + enclosed-volume ON, run crash-free, YL pressure
# ---------------------------------------------------------------------------
class TestKU31Smoke:
    def test_build_runs_crashfree(self, p_cortex, p_ev, p_nuc):
        """SMALL cell + nucleus ON + enclosed-volume ON runs a few hundred
        BAOAB steps with no NaN/Inf (the core KU-3.1 integration smoke)."""
        N_cortex = p_cortex.n_filaments * p_cortex.beads_per_filament
        h = build_cortex_full_simulation(
            p_cortex,
            p_enclosed_volume=p_ev,
            p_nucleus=p_nuc,
            with_baoab=True,
        )
        sim = h["sim"]
        assert sim.state.N_particles == N_cortex + _N_NUC_BEADS
        assert h["enclosed_volume_force"] is not None
        assert h["nucleus_force"] is not None

        sim.run(300)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position).copy()
        assert np.isfinite(pos).all(), (
            "nucleus+enclosed-volume BAOAB produced NaN/Inf"
        )

        # Nucleus stays a bounded cloud near R_nuc; enclosed volume near V0.
        ev = h["enclosed_volume_force"]
        nf = h["nucleus_force"]
        assert nf.last_n == _N_NUC_BEADS
        assert np.isfinite(ev.last_volume) and ev.last_volume > 0.0
        # Volume stays within a few % of V0 (the shell is not deflating).
        assert 0.8 < ev.last_volume / p_ev.V0 < 1.25

    def test_enclosed_volume_pressure_young_laplace_sign_scale(
        self, p_cortex, p_ev, p_nuc
    ):
        """Under a controlled 1% radius COMPRESSION the runtime enclosed-volume
        pressure is POSITIVE (outward / restoring) and of the Young-Laplace
        scale 2γ/R (the KU-3.1 mechanism), confirming the per-bead force is
        wired and non-trivial — NOT a claimed KU-3.1 gate PASS.
        """
        N_cortex = p_cortex.n_filaments * p_cortex.beads_per_filament
        h = build_cortex_full_simulation(
            p_cortex,
            p_enclosed_volume=p_ev,
            p_nucleus=p_nuc,
            with_baoab=True,
        )
        sim = h["sim"]
        ev = h["enclosed_volume_force"]

        # Compress the cortex shell by eps about its centroid (uniform inward
        # scale) so V < V0 → ΔP must be > 0 (outward restoring).
        eps = 0.01
        snap = sim.state.get_snapshot()
        pos = np.asarray(snap.particles.position, dtype=np.float64)
        shell = pos[:N_cortex]
        c = shell.mean(axis=0)
        pos[:N_cortex] = c + (1.0 - eps) * (shell - c)
        sim.state.set_snapshot(snap)
        sim.run(0)  # recompute forces at the compressed config

        dP = ev.last_pressure
        R = ev.last_R_mean
        assert np.isfinite(dP) and np.isfinite(R) and R > 0.0
        # Sign: compression (V<V0) → outward restoring pressure > 0.
        assert dP > 0.0, f"compression should give outward ΔP>0; got {dP}"

        # Scale: ΔP for a δR/R = eps compression is ~K_vol·3·eps (the K_vol
        # anchor makes this ≈ the KU-3.1 interphase ΔP_ref = 40 Pa). The
        # discrete tangent-filament shell is not a perfect sphere, so accept a
        # factor ~2 band about the ideal. The implied surface tension
        # γ = ΔP·R/2 must land in the cortical-tension decade (KU-3.5
        # γ ≈ 0.1–2 mN/m), i.e. ΔP is of the Young-Laplace 2γ/R scale.
        dP_ideal = p_ev.K_vol * 3.0 * eps          # ≈ 40 Pa
        assert 0.3 * dP_ideal < dP < 3.0 * dP_ideal, (
            f"runtime ΔP={dP:.2f} Pa outside the [0.3,3]×{dP_ideal:.1f} Pa "
            "Young-Laplace band for a 1% compression."
        )
        gamma_equiv = dP * R / 2.0                  # ΔP = 2γ/R  → γ
        assert 1.0e-4 < gamma_equiv < 2.0e-3, (
            f"implied surface tension {gamma_equiv*1e3:.3f} mN/m outside the "
            "KU-3.5 cortical-tension decade 0.1–2 mN/m."
        )
        # Cross-check the oracle: 2·γ_equiv/R reproduces ΔP by construction.
        assert math.isclose(
            young_laplace_pressure(gamma_equiv, R), dP, rel_tol=1e-9
        )

    def test_nucleus_confinement_keeps_cloud_bounded(
        self, p_cortex, p_nuc
    ):
        """With the nucleus alone ON, the confinement keeps the bead cloud
        bounded near R_nuc over 300 steps (no escape / blow-up)."""
        N_cortex = p_cortex.n_filaments * p_cortex.beads_per_filament
        h = build_cortex_full_simulation(
            p_cortex, p_nucleus=p_nuc, with_baoab=True
        )
        sim = h["sim"]
        sim.run(300)
        snap = sim.state.get_snapshot()
        pos = np.asarray(snap.particles.position, dtype=np.float64)
        nuc = pos[N_cortex:]
        c = nuc.mean(axis=0)
        r = np.linalg.norm(nuc - c, axis=1)
        assert np.isfinite(pos).all()
        assert (r < 5.0 * p_nuc.R_nuc).all(), (
            f"nucleus cloud unbounded: max r {r.max():.3e} > 5·R_nuc"
        )
