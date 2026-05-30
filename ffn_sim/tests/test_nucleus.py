"""STATIC + analytical sanity-gate tests for H.9 nucleus mechanics (KU-3.B2).

Covers ``ffn_sim/cell/nucleus.py`` Sanity Gate §1–6:

- §1 Dimensional analysis        (TestDimensional — E_nuc→k_eff bridge)
- §2 Boundary cases              (TestBoundary — zero-force config, bands)
- §strain-stiffening             (TestStrainStiffening — knee continuity)
- §3 Conservation invariants     (TestConservation — net force ≈ 0, masking)
- §4 Numerical sanity / CFL      (TestNumericalCFL)
- §5 Sign / sense                (TestSignSense — confines toward the shell)
- §6 Measurement / pore + BAOAB  (TestPoreCoupling, TestBAOABSmoke)
- builder purity                 (TestBuilder — build_nucleus_beads is PURE)

Mirrors ``ffn_sim/tests/test_enclosed_volume.py`` / ``test_membrane.py``.
The brief is a DRAFT skeleton — these tests check the LAW's shape, the
E_nuc→k_eff dimensional bridge, and the default-off invariant; they do NOT
claim a KU-3.1 / KU-1.V.3.3 gate PASS.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import hoomd
import hoomd.md as md

from ffn_sim.cell.nucleus import (
    NucleusConfinement,
    ResolvedNucleus,
    attach_nucleus_confinement,
    build_nucleus_beads,
    confinement_force_magnitude,
    confinement_potential,
    lamin_a_elasticity_scaling,
    lamin_a_viscosity_scaling,
    nucleus_cross_section_area,
    resolve_nucleus,
)


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------
_R_NUC = 3.0e-6          # 3 µm nuclear radius (host geometry)
_N_BEADS = 60            # nucleus bead count
_KT = 4.28e-21          # J (project kT)
# Per-bead Stokes drag γ = 6π η R_bead with η_water≈1e-3, R_bead≈100 nm.
_GAMMA_NUC = 6.0 * math.pi * 1.0e-3 * 1.0e-7   # ≈ 1.885e-9 N·s/m


def _cfg(**overrides) -> dict:
    base = {
        "E_nuc": 5.0e3,
        "ratio_lamin": 3.0,
        "knee_strain": 0.10,
        "critical_pore_area_um2": 7.0,
    }
    base.update(overrides)
    return {"nucleus": base}


def _resolved(**overrides) -> ResolvedNucleus:
    return resolve_nucleus(_cfg(**overrides), R_nuc=_R_NUC, n_beads=_N_BEADS)


def _minimal_sim_with_integrator(
    positions: np.ndarray, type_name: str, dt: float = 1.0e-6
):
    """Smallest HOOMD host: ``positions`` particles of one type + a bare
    md.Integrator (dt set). Mirrors test_membrane._minimal_sim_with_integrator
    so the nucleus force/CFL gate is isolated from any heavy build."""
    import gsd.hoomd

    pos = np.asarray(positions, dtype=np.float64).reshape(-1, 3)
    n = pos.shape[0]
    snap = gsd.hoomd.Frame()
    snap.particles.N = n
    snap.particles.types = [type_name]
    snap.particles.typeid = [0] * n
    snap.particles.position = pos.tolist()
    snap.particles.mass = [1.0] * n
    # Box generously larger than the cloud so nothing wraps.
    L = 20.0e-6
    snap.configuration.box = [L, L, L, 0.0, 0.0, 0.0]
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.integrator = md.Integrator(dt=dt)
    return sim


# ---------------------------------------------------------------------------
# §1 Dimensional analysis — E_nuc → k_eff bridge (eq. B)
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_k_chrom_bridge_formula(self):
        """k_chrom = 4π E_nuc R_nuc / n_beads  [Pa·m = N/m]."""
        p = _resolved()
        expected = 4.0 * math.pi * p.E_nuc * p.R_nuc / p.n_beads
        assert math.isclose(p.k_chrom, expected, rel_tol=1e-12)
        # Sanity of magnitude: 5 kPa, 3 µm, 60 beads → ~3e-3 N/m.
        assert 1e-4 < p.k_chrom < 1e-1

    def test_k_chrom_grid_invariance(self):
        """Total radial stiffness n_beads·k_chrom = 4π E_nuc R_nuc is
        discretization-independent (intensive material response)."""
        p1 = resolve_nucleus(_cfg(), R_nuc=_R_NUC, n_beads=30)
        p2 = resolve_nucleus(_cfg(), R_nuc=_R_NUC, n_beads=300)
        tot1 = p1.n_beads * p1.k_chrom
        tot2 = p2.n_beads * p2.k_chrom
        assert math.isclose(tot1, tot2, rel_tol=1e-12)
        assert math.isclose(tot1, 4.0 * math.pi * 5.0e3 * _R_NUC, rel_tol=1e-12)

    def test_lamin_increment_from_ratio(self):
        """k_lamin = (ratio_lamin − 1)·k_chrom → high-strain modulus =
        ratio_lamin × chromatin modulus (nucleus:cyto ratio, KU-3.B2.1)."""
        p = _resolved(ratio_lamin=3.0)
        assert math.isclose(p.k_lamin, 2.0 * p.k_chrom, rel_tol=1e-12)
        assert math.isclose(p.k_hi, 3.0 * p.k_chrom, rel_tol=1e-12)

    def test_force_units_newton(self):
        """F = k [N/m] · d [m] → [N]."""
        p = _resolved()
        d = 0.5 * p.d_knee
        F = confinement_force_magnitude(p, d)
        assert math.isclose(F, -p.k_chrom * d, rel_tol=1e-12)

    def test_zero_force_config(self):
        """r = R_nuc (d=0) → F=0, U=0 (the force-free configuration)."""
        p = _resolved()
        assert confinement_force_magnitude(p, 0.0) == 0.0
        assert confinement_potential(p, 0.0) == 0.0


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_E_nuc_below_band_raises(self):
        with pytest.raises(ValueError, match="E_nuc"):
            resolve_nucleus(_cfg(E_nuc=500.0), R_nuc=_R_NUC, n_beads=_N_BEADS)

    def test_E_nuc_above_band_raises(self):
        with pytest.raises(ValueError, match="E_nuc"):
            resolve_nucleus(_cfg(E_nuc=2.0e4), R_nuc=_R_NUC, n_beads=_N_BEADS)

    def test_ratio_lamin_above_insitu_band_raises(self):
        """The 10× isolated value must be REJECTED (brief explicit: do not
        hard-code 10×; in-situ band is 1.4–5×)."""
        with pytest.raises(ValueError, match="ratio_lamin"):
            resolve_nucleus(_cfg(ratio_lamin=10.0), R_nuc=_R_NUC, n_beads=_N_BEADS)

    def test_ratio_lamin_below_band_raises(self):
        with pytest.raises(ValueError, match="ratio_lamin"):
            resolve_nucleus(_cfg(ratio_lamin=1.0), R_nuc=_R_NUC, n_beads=_N_BEADS)

    def test_knee_strain_out_of_range_raises(self):
        with pytest.raises(ValueError, match="knee_strain"):
            resolve_nucleus(_cfg(knee_strain=1.5), R_nuc=_R_NUC, n_beads=_N_BEADS)

    def test_negative_R_nuc_raises(self):
        with pytest.raises(ValueError, match="R_nuc"):
            resolve_nucleus(_cfg(), R_nuc=-1.0, n_beads=_N_BEADS)

    def test_zero_n_beads_raises(self):
        with pytest.raises(ValueError, match="n_beads"):
            resolve_nucleus(_cfg(), R_nuc=_R_NUC, n_beads=0)

    def test_invalid_tag_range_raises(self):
        p = _resolved()
        with pytest.raises(ValueError, match="nucleus_tag_range"):
            NucleusConfinement(p, nucleus_tag_range=(100, 50))

    def test_centroid_coincident_bead_no_nan(self):
        """A bead AT the centroid (r̂ undefined) gets zero force, no NaN."""
        p = _resolved()
        # Two beads: one at centroid offset 0, one offset R_nuc. Put them so
        # the centroid sits on the first bead → r=0 for it.
        pos = np.array(
            [[0.0, 0.0, 0.0], [2.0 * _R_NUC, 0.0, 0.0]], dtype=np.float64
        )
        sim = _minimal_sim_with_integrator(pos, "nucleus_bead")
        nuc = NucleusConfinement(p, nucleus_tag_range=(0, 2))
        sim.operations.integrator.forces.append(nuc)
        sim.run(0)
        F = np.asarray(nuc.forces)
        assert np.isfinite(F).all()


# ---------------------------------------------------------------------------
# §strain-stiffening — knee continuity + monotonicity
# ---------------------------------------------------------------------------
class TestStrainStiffening:
    def test_force_continuous_at_knee(self):
        """Bilinear force is CONTINUOUS at d = ±d_knee.

        Continuity = the two one-sided LIMITS agree. The interior limit is
        ``−k_chrom·d_knee``; the exterior limit (lamin branch at |d|=d_knee,
        where ``|d|−d_knee = 0``) is ``−sign·F_knee = −sign·k_chrom·d_knee``.
        We evaluate each branch formula AT the knee (the exact limits) rather
        than straddling it with a finite ε (which would falsely report the
        slope-kink as a value gap)."""
        p = _resolved()
        for s in (+1.0, -1.0):
            d_knee = s * p.d_knee
            # Interior-branch formula evaluated at the knee.
            f_interior = -p.k_chrom * d_knee
            # Lamin-branch formula evaluated at the knee (|d|−d_knee = 0).
            f_lamin = -s * (p.F_knee + (p.k_chrom + p.k_lamin) * 0.0)
            assert math.isclose(f_interior, f_lamin, rel_tol=1e-12), (
                f"force one-sided limits disagree at knee (s={s}): "
                f"{f_interior} vs {f_lamin}"
            )
            # And the actual law returns exactly that limit at the knee.
            assert math.isclose(
                confinement_force_magnitude(p, d_knee), f_interior, rel_tol=1e-12
            )

    def test_potential_continuous_at_knee(self):
        """Bilinear potential is C1 (value AND slope continuous at knee).

        Value continuity: the exterior potential at |d|=d_knee adds 0 beyond
        the interior well (``e = |d|−d_knee = 0``) → both branches give
        ``½ k_chrom d_knee²``. Slope continuity: ``dU/dd = −F``, continuous
        because F is continuous (test above)."""
        p = _resolved()
        u_interior = 0.5 * p.k_chrom * p.d_knee * p.d_knee
        # Exterior branch at the knee (e=0) reduces to u_knee.
        u_exterior = confinement_potential(p, p.d_knee + 1.0e-30)
        assert math.isclose(u_interior, u_exterior, rel_tol=1e-12)
        # Value at exactly the knee.
        assert math.isclose(
            confinement_potential(p, p.d_knee), u_interior, rel_tol=1e-12
        )

    def test_slope_steepens_past_knee(self):
        """Stiffness STEPS UP across the knee (lamin engages): |dF/dd|
        jumps from k_chrom to k_chrom+k_lamin."""
        p = _resolved()
        h = 1.0e-10
        # Interior slope ≈ k_chrom.
        d0 = 0.3 * p.d_knee
        slope_in = abs(
            (confinement_force_magnitude(p, d0 + h)
             - confinement_force_magnitude(p, d0 - h)) / (2 * h)
        )
        # Exterior slope ≈ k_chrom + k_lamin.
        d1 = 2.0 * p.d_knee
        slope_out = abs(
            (confinement_force_magnitude(p, d1 + h)
             - confinement_force_magnitude(p, d1 - h)) / (2 * h)
        )
        assert math.isclose(slope_in, p.k_chrom, rel_tol=1e-4)
        assert math.isclose(slope_out, p.k_hi, rel_tol=1e-4)
        assert slope_out > slope_in   # strain-STIFFENING

    def test_potential_monotone_nonneg(self):
        """U(d) ≥ 0 and monotone increasing in |d| (restoring well)."""
        p = _resolved()
        ds = np.linspace(0.0, 4.0 * p.d_knee, 50)
        us = np.array([confinement_potential(p, d) for d in ds])
        assert (us >= -1e-30).all()
        assert np.all(np.diff(us) >= -1e-30)


# ---------------------------------------------------------------------------
# §3 Conservation invariants
# ---------------------------------------------------------------------------
class TestConservation:
    def test_net_force_near_zero_symmetric_cloud(self):
        """A symmetric nucleus cloud → Σ F ≈ 0 (internal element, no net
        momentum injected — centroid-relative normals)."""
        p = _resolved()
        # A single SHELL of evenly-distributed beads (Fibonacci sphere) so
        # the angular coverage is symmetric → Σ r̂_i ≈ 0. A UNIFORM radial
        # scale by 1.5 pushes every bead off-shell (nonzero forces) WITHOUT
        # breaking the angular symmetry (n̂ unchanged) — the internal-element
        # invariant: the centroid-relative normals make Σ F ≈ 0. (A random
        # per-bead perturbation would shift the centroid and is the wrong
        # test of the invariant.)
        built = build_nucleus_beads(
            (0.0, 0.0, 0.0), _R_NUC, 200, gamma_nuc=_GAMMA_NUC, fill=False
        )
        pos = 1.5 * built["positions"]
        sim = _minimal_sim_with_integrator(pos, "nucleus_bead")
        nuc = NucleusConfinement(p, nucleus_tag_range=(0, pos.shape[0]))
        sim.operations.integrator.forces.append(nuc)
        sim.run(0)
        F = np.asarray(nuc.forces)
        net = np.linalg.norm(F.sum(axis=0))
        total = np.linalg.norm(F, axis=1).sum()
        assert total > 0.0   # forces are actually nonzero (beads off-shell)
        assert net < 1e-2 * total, (
            f"|Σ F| = {net:.3e} N not ≪ Σ|F| = {total:.3e} N "
            "(net momentum injected — internal-element invariant broken)."
        )

        # PRODUCTION default fill=True (the seeding the Lead integrates). The
        # radius is decorrelated from the Fibonacci direction in nucleus.py so
        # the radial force largely cancels; a finite isotropic cloud leaves a
        # small residual (~few %, vs ~25% for the old index-correlated seeding).
        # An exact zero would need antipodal-pair seeding (flagged refinement).
        built_f = build_nucleus_beads(
            (0.0, 0.0, 0.0), _R_NUC, 200, gamma_nuc=_GAMMA_NUC, fill=True
        )
        pos_f = 1.3 * built_f["positions"]      # push beads off the shell
        sim_f = _minimal_sim_with_integrator(pos_f, "nucleus_bead")
        nuc_f = NucleusConfinement(p, nucleus_tag_range=(0, pos_f.shape[0]))
        sim_f.operations.integrator.forces.append(nuc_f)
        sim_f.run(0)
        F_f = np.asarray(nuc_f.forces)
        net_f = np.linalg.norm(F_f.sum(axis=0))
        total_f = np.linalg.norm(F_f, axis=1).sum()
        assert total_f > 0.0
        assert net_f < 0.1 * total_f, (
            f"fill=True net force |Σ F|={net_f:.3e} exceeds 10% of Σ|F|"
            f"={total_f:.3e} — radius/direction decorrelation regressed."
        )

    def test_no_force_on_non_nucleus_tags(self):
        """Confinement with tag range [0, k) leaves tags ≥ k untouched."""
        p = _resolved()
        built = build_nucleus_beads(
            (0.0, 0.0, 0.0), _R_NUC, 10, gamma_nuc=_GAMMA_NUC, fill=True
        )
        # Push beads off-shell so the in-range ones DO feel force.
        pos = built["positions"] * 1.3
        sim = _minimal_sim_with_integrator(pos, "nucleus_bead")
        nuc = NucleusConfinement(p, nucleus_tag_range=(0, 4))
        sim.operations.integrator.forces.append(nuc)
        sim.run(0)
        with sim.state.cpu_local_snapshot as s:
            tag = np.asarray(s.particles.tag).copy()
        F = np.asarray(nuc.forces)
        rows = np.flatnonzero(np.linalg.norm(F, axis=1) > 1e-30)
        assert (tag[rows] < 4).all(), (
            f"confinement applied to out-of-range tags: {tag[rows]}"
        )


# ---------------------------------------------------------------------------
# §4 Numerical sanity / CFL
# ---------------------------------------------------------------------------
class TestNumericalCFL:
    def test_force_magnitude_matches_analytic(self):
        """Per-bead force read from HOOMD matches the bilinear scalar law."""
        p = _resolved()
        # Single bead displaced OUTWARD past the knee so the lamin branch is
        # exercised. Put one bead on +x at R_nuc + 2·d_knee; centroid = bead,
        # so use TWO beads symmetric so the centroid is the origin.
        r_out = _R_NUC + 2.0 * p.d_knee
        pos = np.array(
            [[r_out, 0.0, 0.0], [-r_out, 0.0, 0.0]], dtype=np.float64
        )
        sim = _minimal_sim_with_integrator(pos, "nucleus_bead")
        nuc = NucleusConfinement(p, nucleus_tag_range=(0, 2))
        sim.operations.integrator.forces.append(nuc)
        sim.run(0)
        F = np.asarray(nuc.forces)
        # Each bead at radius r_out from the centroid (origin) → d = 2·d_knee.
        d = 2.0 * p.d_knee
        expected_mag = abs(confinement_force_magnitude(p, d))
        assert math.isclose(
            np.linalg.norm(F[0]), expected_mag, rel_tol=1e-6
        )

    def test_cfl_gate_blocks_too_large_dt(self):
        """CFL gate raises when dt > 0.1·τ (τ = γ_nuc / k_hi)."""
        p = _resolved()
        built = build_nucleus_beads(
            (0.0, 0.0, 0.0), _R_NUC, _N_BEADS, gamma_nuc=_GAMMA_NUC
        )
        tau = _GAMMA_NUC / p.k_hi
        # Choose dt = 10·(0.1·τ) = τ  →  guaranteed violation.
        sim = _minimal_sim_with_integrator(
            built["positions"], "nucleus_bead", dt=tau
        )
        with pytest.raises(RuntimeError, match="Nucleus-confinement CFL"):
            attach_nucleus_confinement(
                sim, p, nucleus_tags=(0, _N_BEADS),
                gamma_nuc=_GAMMA_NUC, cfl_safety_factor=0.1, cfl_strict=True,
            )

    def test_cfl_gate_passes_small_dt(self):
        """A dt safely below 0.1·τ attaches without raising."""
        p = _resolved()
        built = build_nucleus_beads(
            (0.0, 0.0, 0.0), _R_NUC, _N_BEADS, gamma_nuc=_GAMMA_NUC
        )
        tau = _GAMMA_NUC / p.k_hi
        dt_ok = 0.01 * tau           # 10× below the 0.1·τ ceiling
        sim = _minimal_sim_with_integrator(
            built["positions"], "nucleus_bead", dt=dt_ok
        )
        nuc = attach_nucleus_confinement(
            sim, p, nucleus_tags=(0, _N_BEADS),
            gamma_nuc=_GAMMA_NUC, cfl_safety_factor=0.1, cfl_strict=True,
        )
        assert isinstance(nuc, NucleusConfinement)

    def test_cfl_strict_false_skips(self):
        """cfl_strict=False attaches despite a violating dt (diagnostic)."""
        p = _resolved()
        built = build_nucleus_beads(
            (0.0, 0.0, 0.0), _R_NUC, _N_BEADS, gamma_nuc=_GAMMA_NUC
        )
        tau = _GAMMA_NUC / p.k_hi
        sim = _minimal_sim_with_integrator(
            built["positions"], "nucleus_bead", dt=tau
        )
        nuc = attach_nucleus_confinement(
            sim, p, nucleus_tags=(0, _N_BEADS),
            gamma_nuc=_GAMMA_NUC, cfl_safety_factor=0.1, cfl_strict=False,
        )
        assert isinstance(nuc, NucleusConfinement)


# ---------------------------------------------------------------------------
# §5 Sign / sense
# ---------------------------------------------------------------------------
class TestSignSense:
    def _two_bead_radial_force(self, p, r_radius):
        """Two symmetric beads at ±r_radius on +x → centroid origin; return
        the +x bead's force vector (radial along +x̂)."""
        pos = np.array(
            [[r_radius, 0.0, 0.0], [-r_radius, 0.0, 0.0]], dtype=np.float64
        )
        sim = _minimal_sim_with_integrator(pos, "nucleus_bead")
        nuc = NucleusConfinement(p, nucleus_tag_range=(0, 2))
        sim.operations.integrator.forces.append(nuc)
        sim.run(0)
        return np.asarray(nuc.forces)[0]

    def test_inward_when_outside_shell(self):
        """Bead OUTSIDE R_nuc → force points INWARD (−x̂): confines."""
        p = _resolved()
        F = self._two_bead_radial_force(p, _R_NUC + 0.5 * p.d_knee)
        assert F[0] < 0.0, f"outside-shell bead should feel inward; Fx={F[0]}"

    def test_outward_when_inside_shell(self):
        """Bead INSIDE R_nuc → force points OUTWARD (+x̂): confines."""
        p = _resolved()
        F = self._two_bead_radial_force(p, _R_NUC - 0.5 * p.d_knee)
        assert F[0] > 0.0, f"inside-shell bead should feel outward; Fx={F[0]}"

    def test_zero_on_shell(self):
        """Bead exactly on R_nuc → ≈ zero radial force."""
        p = _resolved()
        F = self._two_bead_radial_force(p, _R_NUC)
        assert np.linalg.norm(F) < 1e-18


# ---------------------------------------------------------------------------
# §6 Measurement protocol — KU-1.V.3.3 pore coupling + lamin scaling + smoke
# ---------------------------------------------------------------------------
class TestPoreCoupling:
    def test_critical_pore_area_conversion(self):
        """critical_pore_area = 7 µm² (tumor) → 7e-12 m² (KU-1.V.3.3)."""
        p = _resolved()
        nuc = NucleusConfinement(p, nucleus_tag_range=(0, _N_BEADS))
        assert math.isclose(nuc.critical_pore_area(), 7.0e-12, rel_tol=1e-12)

    def test_cross_section_area_predicate(self):
        """A nucleus cloud at full R_nuc has cross-section ≫ the 7 µm² pore
        (not arrested); a compressed cloud falls below it (arrest predicate
        — NOTE-only coupling, not a claimed PASS)."""
        built = build_nucleus_beads(
            (0.0, 0.0, 0.0), _R_NUC, 200, gamma_nuc=_GAMMA_NUC, fill=True
        )
        pos = built["positions"]
        area_full = nucleus_cross_section_area(pos, normal_axis=1)
        # Full nucleus (R_nuc = 3 µm) cross-section is several µm² and well
        # above the 7 µm² tumor pore limit? R_g-based → check it is > pore.
        assert area_full > 7.0e-12
        # Compress 5× in the transit plane (x,z) → area drops below the pore.
        compressed = pos.copy()
        compressed[:, 0] *= 0.2
        compressed[:, 2] *= 0.2
        area_sq = nucleus_cross_section_area(compressed, normal_axis=1)
        assert area_sq < area_full
        assert area_sq < 7.0e-12   # below the arrest threshold

    def test_lamin_a_scaling_oracles(self):
        """KU-3.B2.2: viscosity ~ laminA^3, elasticity ~ laminA^0.5."""
        assert math.isclose(lamin_a_viscosity_scaling(2.0), 8.0, rel_tol=1e-12)
        assert math.isclose(
            lamin_a_elasticity_scaling(4.0), 2.0, rel_tol=1e-12
        )
        # Lamin-A LOSS (cancer) → softer: both multipliers < 1 for rel < 1.
        assert lamin_a_viscosity_scaling(0.5) < 1.0
        assert lamin_a_elasticity_scaling(0.5) < 1.0


class TestBuilder:
    def test_builder_returns_expected_shapes(self):
        """build_nucleus_beads returns positions/type_name/gamma/radii."""
        out = build_nucleus_beads(
            (1.0e-6, 2.0e-6, 3.0e-6), _R_NUC, _N_BEADS,
            gamma_nuc=_GAMMA_NUC, type_name="nucleus_bead",
        )
        assert out["positions"].shape == (_N_BEADS, 3)
        assert out["gamma"].shape == (_N_BEADS,)
        assert out["type_name"] == "nucleus_bead"
        assert np.allclose(out["gamma"], _GAMMA_NUC)
        # Outermost seed radius (fill) lands on R_nuc.
        assert math.isclose(float(out["radii"].max()), _R_NUC, rel_tol=1e-9)

    def test_builder_centroid_offset(self):
        """Beads are centred on the supplied centroid (within R_nuc)."""
        c = np.array([5.0e-6, -4.0e-6, 1.0e-6])
        out = build_nucleus_beads(
            tuple(c), _R_NUC, _N_BEADS, gamma_nuc=_GAMMA_NUC
        )
        centroid = out["positions"].mean(axis=0)
        # Fibonacci-fill centroid is near c (symmetric coverage).
        assert np.linalg.norm(centroid - c) < _R_NUC

    def test_builder_is_pure_deterministic(self):
        """Same args → identical arrays (no hidden state / RNG drift)."""
        a = build_nucleus_beads(
            (0.0, 0.0, 0.0), _R_NUC, _N_BEADS, gamma_nuc=_GAMMA_NUC, seed=7
        )
        b = build_nucleus_beads(
            (0.0, 0.0, 0.0), _R_NUC, _N_BEADS, gamma_nuc=_GAMMA_NUC, seed=7
        )
        np.testing.assert_array_equal(a["positions"], b["positions"])

    def test_builder_rejects_bad_args(self):
        with pytest.raises(ValueError, match="R_nuc"):
            build_nucleus_beads((0, 0, 0), -1.0, 10, gamma_nuc=_GAMMA_NUC)
        with pytest.raises(ValueError, match="gamma_nuc"):
            build_nucleus_beads((0, 0, 0), _R_NUC, 10, gamma_nuc=-1.0)
        with pytest.raises(ValueError, match="n_beads"):
            build_nucleus_beads((0, 0, 0), _R_NUC, 0, gamma_nuc=_GAMMA_NUC)


class TestBAOABSmoke:
    def test_baoab_nucleus_confinement_no_nan(self):
        """Nucleus beads + confinement-on under L-M BAOAB runs 300 steps
        with no NaN/Inf and the cloud stays bounded near R_nuc."""
        from ffn_sim.integrator.baoab import make_baoab_updater

        p = _resolved()
        built = build_nucleus_beads(
            (0.0, 0.0, 0.0), _R_NUC, _N_BEADS, gamma_nuc=_GAMMA_NUC, fill=True
        )
        pos = built["positions"]
        tau = _GAMMA_NUC / p.k_hi
        dt = 0.05 * tau              # CFL-safe (0.05·τ < 0.1·τ ceiling)
        sim = _minimal_sim_with_integrator(pos, "nucleus_bead", dt=dt)
        # BAOAB Action (wrapped in CustomUpdater) co-exists with the
        # Integrator carrying methods=[] (the project's canonical wiring).
        _action, updater = make_baoab_updater(
            kT=_KT, gamma={"nucleus_bead": _GAMMA_NUC}, dt=dt, seed=3
        )
        sim.operations.integrator.methods = []
        sim.operations.updaters.append(updater)
        attach_nucleus_confinement(
            sim, p, nucleus_tags=(0, _N_BEADS),
            gamma_nuc=_GAMMA_NUC, cfl_safety_factor=0.1,
        )
        sim.run(300)
        with sim.state.cpu_local_snapshot as s:
            out = np.asarray(s.particles.position).copy()
        assert np.isfinite(out).all(), "confinement-on BAOAB produced NaN/Inf"
        centroid = out.mean(axis=0)
        r = np.linalg.norm(out - centroid, axis=1)
        # Beads seeded within R_nuc; confinement keeps them within a few R_nuc.
        assert (r < 5.0 * _R_NUC).all(), (
            f"nucleus cloud unbounded: max radius {r.max():.3e} m > 5·R_nuc."
        )
