"""STATIC + analytical + Brownian-smoke sanity-gate tests for H.10 cytoplasm
Tier-1 (KU-3.B3.1) — ``ffn_sim/cell/cytoplasm.py``.

Covers the module Sanity Gate §1–6:

- §1 Dimensional / ratio law        (TestDimensional)
- §2 Boundary / DEFAULT-OFF         (TestDefaultOff — bit-for-bit identity)
- §3 FDT + viscosity-contrast GATE  (TestViscosityContrastGate — analytical D-ratio
                                      + a standalone HOOMD Brownian smoke: D = kT/γ_b
                                      per type, MCF7 ≈ 5× MDA)
- §4 Numerical / CFL re-check       (TestCFLReCheck — relaxes at the ~5e4× jump)
- §5 Sign-sense                     (TestSignSense — higher η ⇒ slower diffusion)
- resolver / provenance             (TestResolve — cell-type table, validation)

STANDALONE by design (SA-3 shared-tree discipline): NOT built via cell.py. The
Brownian smoke is a few FREE particles under the real BAOAB integrator (mirrors
``tests/test_baoab.py`` plumbing), NOT a full cell. The override is verified to
be the SAME dict BAOAB already consumes, so D = kT/γ_b stays exact per type.

Validation = reproduce the ~5× MCF7/MDA viscosity RATIO (=> D-ratio), NOT the
absolute Pa·s (which is microrheology-probe-length-dependent, η ∝ L²).
"""

from __future__ import annotations

import math
import tempfile

import numpy as np
import pytest

import gsd.hoomd
import hoomd
import hoomd.md as md

from ffn_sim.cell.cytoplasm import (
    DEFAULT_IMMERSED_TYPES,
    ETA_CYTO_BY_CELLTYPE,
    ETA_WATER,
    CFLReCheck,
    CytoplasmTier1,
    apply_cytoplasm_drag,
    cfl_relaxation_at_drag_jump,
    cytoplasm_gamma_override,
    resolve_cytoplasm,
    stokes_drag,
)
from ffn_sim.integrator.baoab import make_baoab_updater

# Reference cell-type viscosities (Hu 2024, PMC10929591; MRS L ≈ 3 µm).
_ETA_MCF7 = ETA_CYTO_BY_CELLTYPE["MCF7"]        # 65.9 Pa·s
_ETA_MDA = ETA_CYTO_BY_CELLTYPE["MDA-MB-231"]   # 12.0 Pa·s
_RATIO_MCF7_MDA = _ETA_MCF7 / _ETA_MDA          # ≈ 5.49×

# A representative water-drag gamma_map (the values cell.py builds at η=water).
# γ_b = 6π η_water R_bead with R_bead = 30 nm (cortex) ≈ 3.91e-10 N·s/m.
_GAMMA_B_WATER = 6.0 * math.pi * ETA_WATER * 30.0e-9


def _water_gamma_map() -> dict[str, float]:
    """A gamma_map keyed by the immersed types (+ one water-exposed type)."""
    gm = {t: _GAMMA_B_WATER for t in DEFAULT_IMMERSED_TYPES}
    # A type NOT in the immersed set — must be left untouched by the override.
    gm["membrane_bead"] = _GAMMA_B_WATER * 1.7
    return gm


# ---------------------------------------------------------------------------
# §1 Dimensional / ratio law
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_stokes_drag_units_value(self):
        # γ = 6π η R; check the numeric matches the closed form.
        g = stokes_drag(_ETA_MCF7, 30.0e-9)
        assert math.isclose(g, 6.0 * math.pi * _ETA_MCF7 * 30.0e-9, rel_tol=0)
        assert g > 0.0 and math.isfinite(g)

    def test_override_is_water_gamma_times_ratio(self):
        # The override scales each immersed γ_b by exactly η_eff/η_water, i.e.
        # γ_eff = 6π η_eff R — dimensionally [N·s/m], and equal to the direct
        # Stokes computation at the cortex radius.
        cyto = resolve_cytoplasm(cell_type="MCF7")
        gm = {"actin_cortex": _GAMMA_B_WATER}
        out = cytoplasm_gamma_override(gm, cyto)
        direct = stokes_drag(_ETA_MCF7, 30.0e-9)
        assert math.isclose(out["actin_cortex"], direct, rel_tol=1e-12)

    def test_diffusion_ratio_is_inverse_viscosity_ratio(self):
        # D_eff/D_water = η_water/η_eff (since D = kT/γ, γ ∝ η).
        cyto = resolve_cytoplasm(cell_type="MCF7")
        assert math.isclose(
            cyto.diffusion_ratio_vs_water(), ETA_WATER / _ETA_MCF7, rel_tol=1e-12
        )


# ---------------------------------------------------------------------------
# §2 Boundary / DEFAULT-OFF — bit-for-bit identity
# ---------------------------------------------------------------------------
class TestDefaultOff:
    def test_none_returns_unchanged_copy(self):
        gm = _water_gamma_map()
        out = apply_cytoplasm_drag(gm, None)
        assert out == gm           # exact equality
        assert out is not gm       # new object (no aliasing of the caller's dict)

    def test_water_eta_is_exact_identity(self):
        # eta_eff == eta_water => every entry is the *input float*, bit-for-bit.
        cyto = resolve_cytoplasm(eta_eff=ETA_WATER)
        assert cyto.is_noop
        gm = _water_gamma_map()
        out = apply_cytoplasm_drag(gm, cyto)
        for k in gm:
            assert out[k] == gm[k]            # bit-for-bit (==, not isclose)
            assert out[k] is gm[k] or out[k] == gm[k]

    def test_cell_type_water_is_noop(self):
        cyto = resolve_cytoplasm(cell_type="water")
        assert cyto.is_noop
        gm = _water_gamma_map()
        assert apply_cytoplasm_drag(gm, cyto) == gm

    def test_default_constructed_is_noop(self):
        assert CytoplasmTier1().is_noop
        assert resolve_cytoplasm().is_noop

    def test_non_immersed_type_untouched_when_enabled(self):
        # membrane_bead is water-exposed (not immersed) => keeps its water drag
        # even when Tier-1 is enabled.
        cyto = resolve_cytoplasm(cell_type="MCF7")
        gm = _water_gamma_map()
        out = cytoplasm_gamma_override(gm, cyto)
        assert out["membrane_bead"] == gm["membrane_bead"]   # unchanged
        # ...while the immersed actin_cortex IS scaled up.
        assert out["actin_cortex"] > gm["actin_cortex"]


# ---------------------------------------------------------------------------
# §3 FDT + viscosity-contrast GATE
# ---------------------------------------------------------------------------
class TestViscosityContrastGate:
    """KU-3.B3.1: per-type η_eff ⇒ MCF7-tagged diffusion ≈ 5× slower than MDA,
    via D = kT/γ_b. Analytical (the override math) + an empirical Brownian smoke
    (D measured from MSD of free beads under the real BAOAB integrator)."""

    def test_analytical_mcf7_about_5x_mda(self):
        # Two cell-type maps from the SAME water map. The immersed-type drag
        # ratio must equal the viscosity ratio ≈ 5.5×.
        gm = {"actin_cortex": _GAMMA_B_WATER}
        g_mcf7 = cytoplasm_gamma_override(gm, resolve_cytoplasm(cell_type="MCF7"))
        g_mda = cytoplasm_gamma_override(gm, resolve_cytoplasm(cell_type="MDA-MB-231"))
        gamma_ratio = g_mcf7["actin_cortex"] / g_mda["actin_cortex"]
        # γ ∝ η ⇒ γ-ratio == η-ratio; D = kT/γ ⇒ D-ratio == 1/γ-ratio.
        assert math.isclose(gamma_ratio, _RATIO_MCF7_MDA, rel_tol=1e-12)
        assert 5.0 <= gamma_ratio <= 6.0          # "≈ 5×" band (Hu 2024: 5.5×)
        # MCF7 diffuses ~5.5× SLOWER than MDA (D-ratio).
        D_ratio_mcf7_over_mda = g_mda["actin_cortex"] / g_mcf7["actin_cortex"]
        assert math.isclose(D_ratio_mcf7_over_mda, 1.0 / _RATIO_MCF7_MDA, rel_tol=1e-12)
        assert D_ratio_mcf7_over_mda < 1.0        # slower

    def test_brownian_smoke_D_equals_kT_over_gamma_per_type(self):
        """Standalone HOOMD Brownian smoke (NOT a cell): free beads of two
        types in one box, drag set per type via the cytoplasm override, no
        conservative force ⇒ pure diffusion. Assert measured D = kT/γ_b per
        type (FDT preserved per-type), hence the MCF7/MDA D-ratio.

        Uses HOOMD-unit scales (kT=1, γ_water=1) for a clean, fast MSD — the
        FDT identity D=kT/γ is scale-free; the cell-type η RATIO (5.5×) is the
        physical content and is applied as the drag scale.
        """
        # --- Build a 2-type free-particle state (no bonds) --------------------
        n_per = 400
        n = 2 * n_per
        box = 1.0e6                  # huge box: never wrap during the smoke
        types = ("mcf7_bead", "mda_bead")
        snap = gsd.hoomd.Frame()
        snap.particles.N = n
        snap.particles.types = list(types)
        tid = np.zeros(n, dtype=np.uint32)
        tid[n_per:] = 1             # second half is mda_bead
        snap.particles.typeid = tid
        snap.particles.position = np.zeros((n, 3), dtype=float)
        snap.particles.mass = np.ones(n)
        snap.configuration.box = [box, box, box, 0.0, 0.0, 0.0]
        fp = tempfile.NamedTemporaryFile(suffix=".gsd", delete=False).name
        with gsd.hoomd.open(fp, mode="w") as f:
            f.append(snap)

        kT = 1.0
        dt = 1.0e-3
        gamma_water = 1.0           # HOOMD-unit "water" drag

        # Per-type water gamma_map (both types start at water drag = 1.0).
        gm_water = {t: gamma_water for t in types}
        # Engage the override at the MCF7/MDA viscosity contrast. Treat the two
        # bead types as the two cell lines' immersed beads.
        cyto_mcf7 = resolve_cytoplasm(
            eta_eff=_ETA_MCF7, eta_water=ETA_WATER, immersed_types=("mcf7_bead",)
        )
        cyto_mda = resolve_cytoplasm(
            eta_eff=_ETA_MDA, eta_water=ETA_WATER, immersed_types=("mda_bead",)
        )
        gm = cytoplasm_gamma_override(gm_water, cyto_mcf7)
        gm = cytoplasm_gamma_override(gm, cyto_mda)
        gamma_mcf7 = gamma_water * (_ETA_MCF7 / ETA_WATER)
        gamma_mda = gamma_water * (_ETA_MDA / ETA_WATER)
        assert math.isclose(gm["mcf7_bead"], gamma_mcf7, rel_tol=1e-12)
        assert math.isclose(gm["mda_bead"], gamma_mda, rel_tol=1e-12)
        # γ ∝ η so the drag (hence 1/D) ratio is the η ratio.
        assert math.isclose(gm["mcf7_bead"] / gm["mda_bead"], _RATIO_MCF7_MDA,
                            rel_tol=1e-12)

        # --- Run BAOAB on free beads (forces=[], methods=[]) ------------------
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=2024)
        sim.create_state_from_gsd(filename=fp)
        ig = md.Integrator(dt=dt)   # no forces => pure diffusion
        sim.operations.integrator = ig
        action, updater = make_baoab_updater(kT=kT, gamma=gm, dt=dt, seed=7)
        sim.operations.updaters.append(updater)

        sim.run(0)                  # seed net_force (=0) before first step
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag).copy()
            tid_by_tag = np.empty(n, dtype=np.int64)
            tid_by_tag[tag] = np.asarray(s.particles.typeid)
            pos0 = np.empty((n, 3))
            pos0[tag] = np.asarray(s.particles.position)

        n_steps = 4000
        sim.run(n_steps)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag)
            posT = np.empty((n, 3))
            posT[tag] = np.asarray(s.particles.position)

        # Unwrapped displacement (box never wrapped — huge box).
        disp2 = np.sum((posT - pos0) ** 2, axis=1)           # (n,)
        t = n_steps * dt
        # D = <Δr²>/(2 d t), d=3.
        msd_mcf7 = disp2[tid_by_tag == 0].mean()
        msd_mda = disp2[tid_by_tag == 1].mean()
        D_mcf7 = msd_mcf7 / (2.0 * 3.0 * t)
        D_mda = msd_mda / (2.0 * 3.0 * t)

        D_exp_mcf7 = kT / gamma_mcf7
        D_exp_mda = kT / gamma_mda

        # FDT: measured D ≈ kT/γ_b per type (Stokes-Einstein). 400 beads ⇒
        # statistical scatter on MSD ~ √(2/N) ≈ 7%; use a 15% gate.
        assert math.isclose(D_mcf7, D_exp_mcf7, rel_tol=0.15), (
            f"MCF7 D={D_mcf7:.4g} vs kT/γ={D_exp_mcf7:.4g}")
        assert math.isclose(D_mda, D_exp_mda, rel_tol=0.15), (
            f"MDA D={D_mda:.4g} vs kT/γ={D_exp_mda:.4g}")

        # The discriminator: MCF7 diffuses ≈ 5.5× SLOWER than MDA.
        meas_ratio = D_mda / D_mcf7
        assert math.isclose(meas_ratio, _RATIO_MCF7_MDA, rel_tol=0.20), (
            f"measured D_mda/D_mcf7={meas_ratio:.3g} vs η-ratio "
            f"{_RATIO_MCF7_MDA:.3g}")
        assert 4.5 <= meas_ratio <= 6.5     # "≈ 5×" contrast, empirically


# ---------------------------------------------------------------------------
# §4 Numerical / CFL re-check at the drag jump
# ---------------------------------------------------------------------------
class TestCFLReCheck:
    def test_cfl_relaxes_at_worst_case_jump(self):
        # Worst-case: MCF7 (65.9 Pa·s) over water (6.913e-4) ≈ 9.5×10⁴× drag.
        res = cfl_relaxation_at_drag_jump(
            gamma_water=_GAMMA_B_WATER,
            k_eff=1.0e-3,                 # a representative stiff spring [N/m]
            cfl_safety=0.1,
            eta_water=ETA_WATER,
            eta_eff=_ETA_MCF7,
            bead_mass=1.0e-18,
        )
        assert isinstance(res, CFLReCheck)
        # τ = γ_b/k_eff scales with drag => the dt ceiling GROWS => relaxes.
        assert res.relaxes is True
        assert res.dt_ceiling_cyto > res.dt_ceiling_water
        assert math.isclose(
            res.dt_ceiling_cyto / res.dt_ceiling_water,
            _ETA_MCF7 / ETA_WATER, rel_tol=1e-9,
        )
        assert res.ratio > 1.0e4          # the documented ~5×10⁴× class jump
        # τ_v = m/γ SHRINKS with drag (deeper overdamped limit — informational).
        assert res.tau_v_cyto < res.tau_v_water

    def test_cfl_noop_at_water(self):
        res = cfl_relaxation_at_drag_jump(
            gamma_water=_GAMMA_B_WATER, k_eff=1.0e-3, cfl_safety=0.1,
            eta_water=ETA_WATER, eta_eff=ETA_WATER,
        )
        assert res.relaxes is True            # ≥ (equal) still counts as safe
        assert math.isclose(res.ratio, 1.0, rel_tol=0)
        assert math.isclose(res.dt_ceiling_cyto, res.dt_ceiling_water, rel_tol=0)

    def test_cfl_rejects_bad_inputs(self):
        with pytest.raises(ValueError):
            cfl_relaxation_at_drag_jump(
                gamma_water=-1.0, k_eff=1.0, cfl_safety=0.1, eta_eff=10.0
            )
        with pytest.raises(ValueError):
            cfl_relaxation_at_drag_jump(
                gamma_water=1.0, k_eff=0.0, cfl_safety=0.1, eta_eff=10.0
            )


# ---------------------------------------------------------------------------
# §5 Sign-sense
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_higher_eta_lower_diffusion(self):
        # MCF7 (more viscous) must diffuse slower than MDA: D ∝ 1/η.
        d_mcf7 = resolve_cytoplasm(cell_type="MCF7").diffusion_ratio_vs_water()
        d_mda = resolve_cytoplasm(cell_type="MDA-MB-231").diffusion_ratio_vs_water()
        assert d_mcf7 < d_mda                 # MCF7 slower (lower D)
        # Both below water (enabled cytoplasm => slower than water).
        assert d_mcf7 < 1.0 and d_mda < 1.0

    def test_enabling_increases_drag(self):
        cyto = resolve_cytoplasm(cell_type="MCF7")
        assert cyto.viscosity_ratio > 1.0     # drag goes UP (slower), right sign


# ---------------------------------------------------------------------------
# resolver / provenance
# ---------------------------------------------------------------------------
class TestResolve:
    def test_eta_table_matches_literature(self):
        # Magic-Number Block provenance (Hu 2024 / MRS).
        assert ETA_CYTO_BY_CELLTYPE["MCF10A"] == 41.6
        assert ETA_CYTO_BY_CELLTYPE["MCF7"] == 65.9
        assert ETA_CYTO_BY_CELLTYPE["MDA-MB-231"] == 12.0
        assert ETA_CYTO_BY_CELLTYPE["water"] == ETA_WATER
        # The discriminator the unit validates: MCF7 ≈ 5× MDA.
        assert 5.0 <= _RATIO_MCF7_MDA <= 6.0

    def test_eta_eff_overrides_cell_type(self):
        cyto = resolve_cytoplasm(eta_eff=99.0, cell_type="MDA-MB-231")
        assert cyto.eta_eff == 99.0           # explicit wins

    def test_unknown_cell_type_raises(self):
        with pytest.raises(KeyError):
            resolve_cytoplasm(cell_type="HeLa")

    def test_bad_eta_raises(self):
        with pytest.raises(ValueError):
            resolve_cytoplasm(eta_eff=-1.0)
        with pytest.raises(ValueError):
            resolve_cytoplasm(eta_eff=0.0)
        with pytest.raises(ValueError):
            CytoplasmTier1(eta_eff=math.inf)

    def test_custom_immersed_types(self):
        cyto = resolve_cytoplasm(cell_type="MCF7", immersed_types=("actin_cortex",))
        gm = _water_gamma_map()
        out = cytoplasm_gamma_override(gm, cyto)
        assert out["actin_cortex"] > gm["actin_cortex"]      # scaled
        # xlink_head is in DEFAULT set but NOT in this custom set => untouched.
        assert out["xlink_head"] == gm["xlink_head"]

    def test_override_is_pure_no_mutation(self):
        gm = _water_gamma_map()
        gm_snapshot = dict(gm)
        cyto = resolve_cytoplasm(cell_type="MCF7")
        _ = cytoplasm_gamma_override(gm, cyto)
        assert gm == gm_snapshot              # input dict not mutated
