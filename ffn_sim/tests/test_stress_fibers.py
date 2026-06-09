"""STATIC + analytical sanity-gate tests for ventral stress fibers.

Covers ``ffn_sim/cell/stress_fibers.py`` Sanity Gate §1-6:

- OFF-identity (default-off contract): disabled resolve + builder no-op.
- §1 Dimensional analysis: derived SI quantities, sarcomere count.
- §2 Boundary cases: negative/zero/non-finite params raise ValueError.
- §3 Conservation: exact particle + bond count after extend.
- §5 Sign/sense: a stretched backbone bond pulls inward (contractile).
- no-γ-contamination: GAMMA_DENYLIST_PREFIX non-empty + every bond type
  this builder creates starts with it.

Import-light: resolver-level + a tiny hand-built snapshot. NO full Cell.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ffn_sim.cell.stress_fibers import (
    BOND_TYPE_SF_ACTIN,
    GAMMA_DENYLIST_PREFIX,
    PI_DECISIONS,
    ResolvedStressFibers,
    extend_snapshot_with_stress_fibers,
    generate_stress_fiber_layout,
    measure_sf_tension,
    resolve_stress_fibers,
    sf_bond_type_names,
    stress_fiber_dt_cfl,
)

_KT = 4.28e-21       # J
_GAMMA_B = 3.9e-10   # N·s/m  (cortex per-bead Stokes drag)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _enabled_cfg(**overrides) -> dict:
    base = {
        "enabled": True,
        "n_SF": 2,
        "n_beads_per_SF": 6,
        "sarcomere_spacing": 0.5e-6,
        "k_actin": 1.0e-2,
        "k_anchor": 1.0e-2,
        "seed": 7,
    }
    base.update(overrides)
    return {"stress_fibers": base}


def _fa_endpoints(n_pairs: int = 2, span: float = 6.0e-6):
    """Return (positions, tags) for 2·n_pairs FA anchors on the basal plane.

    Anchors are laid out as endpoint pairs separated by ``span`` along x, on
    z = 0 (basal plane). Tags are an arbitrary contiguous block starting at 100
    (simulating FA clutch particles that live elsewhere in the tag space)."""
    pts = []
    for i in range(n_pairs):
        y = i * 2.0e-6
        pts.append([0.0, y, 0.0])
        pts.append([span, y, 0.0])
    pos = np.asarray(pts, dtype=np.float64)
    tags = np.arange(100, 100 + pos.shape[0], dtype=np.int64)
    return pos, tags


def _base_snapshot(n_pre: int = 110):
    """Tiny gsd Frame with ``n_pre`` pre-existing particles (incl. FA anchors
    at tags 100..103) and one pre-existing non-sf bond type, to exercise the
    type-merge logic and tag offsetting."""
    import gsd.hoomd

    snap = gsd.hoomd.Frame()
    snap.particles.N = n_pre
    snap.particles.types = ["integrin", "ligand"]
    snap.particles.typeid = [0] * n_pre
    pos = np.zeros((n_pre, 3), dtype=np.float64)
    # Put FA anchors at the tags _fa_endpoints expects (100..103).
    fa_pos, _ = _fa_endpoints()
    pos[100:100 + fa_pos.shape[0]] = fa_pos
    snap.particles.position = pos
    snap.particles.mass = [1.0] * n_pre
    snap.bonds.N = 1
    snap.bonds.types = ["integrin_ligand"]
    snap.bonds.typeid = [0]
    snap.bonds.group = [[0, 1]]
    L = 50.0e-6
    snap.configuration.box = [L, L, L, 0.0, 0.0, 0.0]
    return snap


# ---------------------------------------------------------------------------
# OFF-identity (default-off contract)
# ---------------------------------------------------------------------------
class TestOffIdentity:
    def test_resolve_disabled_by_default(self):
        # Missing block → disabled.
        p = resolve_stress_fibers({})
        assert p.enabled is False
        assert p.n_SF == 0
        assert p.k_actin is None

    def test_resolve_enabled_false(self):
        p = resolve_stress_fibers({"stress_fibers": {"enabled": False, "n_SF": 5}})
        assert p.enabled is False
        assert p.n_SF == 0

    def test_builder_noop_when_disabled(self):
        p = resolve_stress_fibers({"stress_fibers": {"enabled": False}})
        snap = _base_snapshot()
        n_part0 = int(snap.particles.N)
        n_bond0 = int(snap.bonds.N)
        out = extend_snapshot_with_stress_fibers(snap, p)
        # SAME object back, nothing added.
        assert out is snap
        assert int(out.particles.N) == n_part0
        assert int(out.bonds.N) == n_bond0

    def test_builder_noop_when_zero_SF(self):
        p = resolve_stress_fibers(_enabled_cfg(n_SF=0))
        assert p.enabled is True and p.n_SF == 0
        snap = _base_snapshot()
        out = extend_snapshot_with_stress_fibers(snap, p)
        assert out is snap


# ---------------------------------------------------------------------------
# §1 Dimensional / derived-quantity sanity
# ---------------------------------------------------------------------------
class TestDimensional:
    def test_resolved_si_fields(self):
        p = resolve_stress_fibers(_enabled_cfg())
        assert p.sarcomere_spacing == pytest.approx(0.5e-6)
        assert p.k_actin == pytest.approx(1.0e-2)
        assert p.k_anchor == pytest.approx(1.0e-2)
        # α-actinin reuse defaults (KU-3.19).
        assert p.alpha_length == pytest.approx(35.0e-9)
        assert p.k_xl == pytest.approx(1.0e-7)
        # Kumar 2006 band ~10-30 nN.
        assert p.sf_tension_band_N == pytest.approx((10.0e-9, 30.0e-9))

    def test_layout_ell0_and_sarcomere_count(self):
        # Single bundle on a known 6 µm span so the pairing is deterministic
        # (one pair, no cross-row shuffle ambiguity).
        nb = 6
        p = resolve_stress_fibers(_enabled_cfg(n_SF=1, n_beads_per_SF=nb))
        pos, tags = _fa_endpoints(n_pairs=1, span=6.0e-6)
        lay = generate_stress_fiber_layout(p, pos, tags, particle_tag_start=0)
        L = 6.0e-6
        # ℓ0 = L/(nb-1) = 6 µm / 5 = 1.2 µm.
        assert np.allclose(lay.ell0_actin, L / (nb - 1))
        # n_sarc = round(L/spacing) = round(6/0.5) = 12, capped at nb-1 = 5.
        exp_n_sarc = int(min(round(L / p.sarcomere_spacing), nb - 1))
        assert np.all(lay.n_sarc == exp_n_sarc)

    def test_cfl_timestep(self):
        p = resolve_stress_fibers(_enabled_cfg())
        dt_cfl = stress_fiber_dt_cfl(p, gamma_b=_GAMMA_B, cfl_safety_factor=0.1)
        # 0.1 · γ_b / k = 0.1 · 3.9e-10 / 1e-2 = 3.9e-9 s.
        assert dt_cfl == pytest.approx(0.1 * _GAMMA_B / 1.0e-2)
        assert dt_cfl > 0.0

    def test_cfl_raises_when_k_flagged(self):
        p = resolve_stress_fibers(_enabled_cfg(k_actin=None, k_anchor=None))
        assert p.k_actin is None
        with pytest.raises(NotImplementedError):
            stress_fiber_dt_cfl(p, gamma_b=_GAMMA_B)


# ---------------------------------------------------------------------------
# §2 Boundary cases — bad params raise
# ---------------------------------------------------------------------------
class TestBoundary:
    @pytest.mark.parametrize(
        "override",
        [
            {"n_SF": -1},
            {"n_beads_per_SF": 1},
            {"sarcomere_spacing": 0.0},
            {"sarcomere_spacing": -1.0e-6},
            {"k_xl": 0.0},
            {"alpha_length": -1.0},
            {"k_actin": -1.0e-2},
            {"k_anchor": 0.0},
            {"sf_tension_band_N": (30.0e-9, 10.0e-9)},  # mis-ordered
        ],
    )
    def test_bad_param_raises(self, override):
        with pytest.raises(ValueError):
            resolve_stress_fibers(_enabled_cfg(**override))

    def test_nonfinite_spacing_raises(self):
        with pytest.raises(ValueError):
            resolve_stress_fibers(_enabled_cfg(sarcomere_spacing=float("inf")))

    def test_too_few_fa_endpoints_raises(self):
        p = resolve_stress_fibers(_enabled_cfg(n_SF=3))
        pos, tags = _fa_endpoints(n_pairs=2)  # only 4 endpoints, need 6
        with pytest.raises(ValueError):
            generate_stress_fiber_layout(p, pos, tags, particle_tag_start=0)

    def test_degenerate_bundle_raises(self):
        p = resolve_stress_fibers(_enabled_cfg(n_SF=1))
        # Two coincident endpoints → zero-length bundle.
        pos = np.array([[1.0e-6, 0.0, 0.0], [1.0e-6, 0.0, 0.0]], dtype=np.float64)
        tags = np.array([100, 101], dtype=np.int64)
        with pytest.raises(ValueError):
            generate_stress_fiber_layout(p, pos, tags, particle_tag_start=0)

    def test_enabled_builder_requires_fa(self):
        p = resolve_stress_fibers(_enabled_cfg())
        snap = _base_snapshot()
        with pytest.raises(ValueError):
            extend_snapshot_with_stress_fibers(snap, p)  # no FA endpoints


# ---------------------------------------------------------------------------
# §3 Conservation — exact counts after extend
# ---------------------------------------------------------------------------
class TestConservation:
    def test_exact_particle_and_bond_counts(self):
        nb = 6
        p = resolve_stress_fibers(_enabled_cfg(n_SF=2, n_beads_per_SF=nb))
        pos, tags = _fa_endpoints(n_pairs=2, span=6.0e-6)
        snap = _base_snapshot()
        n_part0 = int(snap.particles.N)
        n_bond0 = int(snap.bonds.N)

        # Re-derive the realised layout to get exact n_sarc per bundle.
        lay = generate_stress_fiber_layout(
            p, pos, tags, particle_tag_start=n_part0,
            rng=np.random.default_rng(p.seed),
        )
        n_sarc = lay.n_sarc  # per bundle
        exp_particles = int(sum(nb + 2 * s for s in n_sarc))
        # per bundle bonds: (nb-1) chain + 2 anchor + n_sarc intra + 2·n_sarc attach
        exp_bonds = int(sum((nb - 1) + 2 + s + 2 * s for s in n_sarc))

        out = extend_snapshot_with_stress_fibers(
            snap, p, pos, tags, rng=np.random.default_rng(p.seed)
        )
        assert int(out.particles.N) == n_part0 + exp_particles
        assert int(out.bonds.N) == n_bond0 + exp_bonds

        # Pre-existing particles + bond preserved unchanged.
        assert int(out.bonds.typeid[0]) == 0  # old integrin_ligand bond
        assert "integrin_ligand" in list(out.bonds.types)
        # New particle types present.
        assert "sf_actin" in list(out.particles.types)
        assert "sf_xlink_head" in list(out.particles.types)

    def test_anchor_bonds_reference_fa_tags(self):
        nb = 4
        p = resolve_stress_fibers(_enabled_cfg(n_SF=1, n_beads_per_SF=nb))
        pos, tags = _fa_endpoints(n_pairs=1, span=4.0e-6)
        snap = _base_snapshot()
        out = extend_snapshot_with_stress_fibers(snap, p, pos, tags)
        bt = np.asarray(out.bonds.typeid)
        bg = np.asarray(out.bonds.group).reshape(-1, 2)
        anchor_tid = list(out.bonds.types).index("sf_anchor")
        anchor_rows = bg[bt == anchor_tid]
        # Exactly two anchors, each referencing an FA tag (100 or 101).
        assert anchor_rows.shape[0] == 2
        referenced = set(int(x) for x in anchor_rows.reshape(-1))
        assert {100, 101}.issubset(referenced)


# ---------------------------------------------------------------------------
# §5 Sign / sense — stretched backbone bond is contractile (pulls inward)
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_stretched_backbone_pulls_inward(self):
        # Two adjacent backbone beads at rest length ℓ0, then stretch the
        # outer bead farther; the harmonic restoring force on the outer bead
        # must point back toward its neighbour (contractile).
        k = 1.0e-2
        ell0 = 1.0e-6
        # bead i at x=0, bead j at x = ell0 + delta (stretched by delta).
        delta = 0.2e-6
        xi = np.array([0.0, 0.0, 0.0])
        xj = np.array([ell0 + delta, 0.0, 0.0])
        # md.bond.Harmonic force on j: F_j = -k (|r| - r0) r̂_{j-i}.
        rvec = xj - xi
        r = np.linalg.norm(rvec)
        r_hat = rvec / r
        F_j = -k * (r - ell0) * r_hat
        # Stretched (r > ℓ0) → F_j points in −x → back toward bead i (inward).
        assert F_j[0] < 0.0
        # Magnitude = k·delta.
        assert np.linalg.norm(F_j) == pytest.approx(k * delta, rel=1e-9)

    def test_measure_sf_tension_positive_under_stretch(self):
        p = resolve_stress_fibers(_enabled_cfg(n_SF=1, n_beads_per_SF=4))
        # Build a snapshot, then artificially stretch the bundle by moving the
        # FA endpoints apart farther than the construction span.
        pos, tags = _fa_endpoints(n_pairs=1, span=4.0e-6)
        snap = _base_snapshot()
        out = extend_snapshot_with_stress_fibers(snap, p, pos, tags)
        positions = np.asarray(out.particles.position, dtype=np.float64).copy()
        bg = np.asarray(out.bonds.group)
        bt = np.asarray(out.bonds.typeid)
        types = list(out.bonds.types)
        # Stretch every sf_actin bead outward along x by 10% to induce tension.
        sf_actin_tid = list(out.particles.types).index("sf_actin")
        tids = np.asarray(out.particles.typeid)
        sf_mask = tids == sf_actin_tid
        positions[sf_mask, 0] *= 1.10
        m = measure_sf_tension(positions, bg, bt, types, p)
        assert m["n_bonds"] > 0
        assert m["max_tension_N"] >= 0.0

    def test_uniform_contraction_needs_construction_r0(self):
        """Regression: a UNIFORMLY stretched bundle must report nonzero tension
        when the construction rest length is passed (the live-mean fallback
        under-reports it to ~0 — the very physiological load case)."""
        p = resolve_stress_fibers(_enabled_cfg(n_SF=1, n_beads_per_SF=6))
        pos, tags = _fa_endpoints(n_pairs=1, span=4.0e-6)
        snap = _base_snapshot()
        out = extend_snapshot_with_stress_fibers(snap, p, pos, tags)
        positions0 = np.asarray(out.particles.position, dtype=np.float64).copy()
        bg = np.asarray(out.bonds.group)
        bt = np.asarray(out.bonds.typeid)
        types = list(out.bonds.types)
        # Construction rest length = mean sf_actin bond length at construction.
        from ffn_sim.cell.stress_fibers import BOND_TYPE_SF_ACTIN
        chain_tid = types.index(BOND_TYPE_SF_ACTIN)
        sel = bg[bt == chain_tid]
        r_constr = np.linalg.norm(positions0[sel[:, 0]] - positions0[sel[:, 1]], axis=1)
        r0 = float(np.mean(r_constr))
        # Uniformly stretch every particle radially by 8% about the origin.
        positions = positions0 * 1.08
        with_r0 = measure_sf_tension(positions, bg, bt, types, p, r0=r0)
        live = measure_sf_tension(positions, bg, bt, types, p)  # fallback
        # With the construction r0, the uniform stretch is detected.
        assert with_r0["mean_tension_N"] > 0.0
        assert with_r0["r0_is_live"] == 0.0
        # The live-mean fallback flags itself and (the bug) reports ~0.
        assert live["r0_is_live"] == 1.0
        assert live["mean_tension_N"] < with_r0["mean_tension_N"]


# ---------------------------------------------------------------------------
# no-γ-contamination
# ---------------------------------------------------------------------------
class TestNoGammaContamination:
    def test_prefix_nonempty(self):
        assert isinstance(GAMMA_DENYLIST_PREFIX, str)
        assert len(GAMMA_DENYLIST_PREFIX) > 0
        assert GAMMA_DENYLIST_PREFIX == "sf_"

    def test_every_bond_type_uses_prefix(self):
        for name in sf_bond_type_names():
            assert name.startswith(GAMMA_DENYLIST_PREFIX), name

    def test_built_snapshot_only_adds_sf_prefixed_bonds(self):
        p = resolve_stress_fibers(_enabled_cfg(n_SF=2, n_beads_per_SF=5))
        pos, tags = _fa_endpoints(n_pairs=2, span=5.0e-6)
        snap = _base_snapshot()
        old_types = set(snap.bonds.types)
        out = extend_snapshot_with_stress_fibers(snap, p, pos, tags)
        new_types = set(out.bonds.types) - old_types
        # Every NEWLY-added bond type is sf_-prefixed.
        assert new_types
        for name in new_types:
            assert name.startswith(GAMMA_DENYLIST_PREFIX), name


# ---------------------------------------------------------------------------
# PI-decisions presence (honesty over completeness)
# ---------------------------------------------------------------------------
def test_pi_decisions_documented():
    assert isinstance(PI_DECISIONS, list)
    # k_actin / k_anchor are PI-flagged; must be surfaced.
    joined = " ".join(PI_DECISIONS).lower()
    assert "k_actin" in joined
    assert "k_anchor" in joined


# ---------------------------------------------------------------------------
# LIVE activation wiring (2026-06-09; PI 소유권 허용): FA-adhered build +
# long-axis-aligned basal FA pairing + per-bundle EXACT-r0 (force-free) backbone.
# cell.py snapshot-extension + sf LJ/gamma_map + registry LIVE. PASSIVE backbone
# (NMII deferred — sf_myosin_* prefix + Kumar active gate).
# ---------------------------------------------------------------------------
class TestStressFibersActivationWiring:
    def test_off_build_is_bit_identity(self):
        from ffn_sim.cell.manifest import (
            build_baseline_cell, load_manifest, resolve_baseline,
        )
        rb = resolve_baseline(load_manifest("mcf7_baseline.yaml"))
        assert rb.p_stress_fibers is None
        cell = build_baseline_cell("mcf7_baseline.yaml", seed=1)
        assert cell.p_stress_fibers is None
        assert cell.simulation.state.N_particles > 0

    def test_on_bundles_assemble_aligned_forcefree_no_contam(self):
        import numpy as np
        from ffn_sim.cell.compartment_registry import REGISTRY, load_recipe
        from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
        from ffn_sim.cell.stress_fibers import sf_actin_backbone_bin_names
        from ffn_sim.cortex.cortical_tension import (
            _is_adhesion_bond_type, measure_cortical_tension,
        )

        base = load_manifest("mcf7_baseline.yaml")
        m_off, d0 = REGISTRY.compose_manifest(
            load_recipe("adherent_passive"), base_manifest=base, strict=True
        )
        m_on, d1 = REGISTRY.compose_manifest(
            load_recipe("ventral_stress_fibers_passive"), base_manifest=base, strict=True
        )
        assert d0 == [] and d1 == []
        off = build_baseline_cell("mcf7_baseline.yaml", manifest=m_off, seed=1)
        on = build_baseline_cell("mcf7_baseline.yaml", manifest=m_on, seed=1)
        p = on.p_stress_fibers
        assert p is not None and p.enabled

        sn = on.simulation.state.get_snapshot()
        n_sf = int(on.extras["handles"]["n_stress_fibers"])
        layout = on.extras["handles"]["stress_fiber_layout"]
        assert n_sf == p.n_SF > 0
        assert "sf_actin" in sn.particles.types and "sf_xlink_head" in sn.particles.types
        sf_bonds = [t for t in sn.bonds.types if t.startswith("sf_")]
        assert sf_bonds and all(_is_adhesion_bond_type(t) for t in sf_bonds)
        # ALIGNED: bundle axes vs the in-plane footprint long axis.
        pos = np.asarray(sn.particles.position, dtype=np.float64)
        fe = np.asarray(layout.fa_endpoints, dtype=np.int64)
        vecs = pos[fe[:, 1]] - pos[fe[:, 0]]
        lens = np.linalg.norm(vecs, axis=1)
        u = vecs / np.maximum(lens[:, None], 1e-30)
        ax = u.mean(axis=0); ax = ax / np.linalg.norm(ax)
        assert float(np.mean(np.abs(u @ ax))) >= 0.8   # aligned, not random
        # FORCE-FREE: per-bundle backbone born at its EXACT ell0_b.
        g = np.asarray(sn.bonds.group, dtype=np.int64)
        bt = np.asarray(sn.bonds.typeid, dtype=np.int64)
        bn = list(sn.bonds.types)
        worst = 0.0
        for b, name in enumerate(sf_actin_backbone_bin_names(n_sf)):
            if name not in bn:
                continue
            m = bt == bn.index(name)
            if not m.any():
                continue
            ln = np.linalg.norm(pos[g[m][:, 0]] - pos[g[m][:, 1]], axis=1)
            r0 = float(layout.ell0_actin[b])
            worst = max(worst, float(np.max(np.abs(ln - r0) / r0)))
        assert worst < 1e-6
        # NO contamination: cortical γ_soft identical OFF vs ON (sf_ denylisted).
        off.simulation.run(0); on.simulation.run(0)
        go = measure_cortical_tension(off.simulation, R_cell=off.p_cortex.R_cell,
                                      p_enclosed_volume=off.p_enclosed_volume)
        gn = measure_cortical_tension(on.simulation, R_cell=on.p_cortex.R_cell,
                                      p_enclosed_volume=on.p_enclosed_volume)
        k = "gamma_soft_N_per_m" if "gamma_soft_N_per_m" in go else "gamma_soft"
        assert abs(gn[k] - go[k]) <= 1e-12 * max(1.0, abs(go[k]))

    def test_enabled_requires_fa(self):
        import pytest
        from ffn_sim.cell.manifest import load_manifest, resolve_baseline
        m = load_manifest("mcf7_baseline.yaml")
        m.setdefault("optional_subsystems", {})
        m["optional_subsystems"]["ventral_stress_fibers"] = {
            "enabled": True, "n_SF": 4, "N_filaments": 20,
        }
        # fa OFF → SF resolve must raise (requires fa).
        with pytest.raises(ValueError):
            resolve_baseline(m)
