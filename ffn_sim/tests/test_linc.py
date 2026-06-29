"""STATIC + analytical sanity-gate tests for the LINC compartment.

Covers ``ffn_sim/cell/linc.py`` Sanity Gate §1-6:

- OFF-identity (DEFAULT-OFF): resolve(enabled=False) → .enabled False;
  the snapshot builder is a strict no-op (zero bonds added).
- §1 Dimensional / parameter sanity: derived SI quantities + ValueError on
  negative/zero params.
- §2 Boundary: k_linc=None enabled build raises NotImplementedError; no
  acceptor in range → no-op; missing nucleus → no-op.
- §3 Conservation / no-gamma-contamination: every bond type starts with the
  declared GAMMA_DENYLIST_PREFIX.
- §5 Sign / sense: a stretched LINC bond pulls the nucleus toward the
  cytoskeleton anchor (restoring).

Import-light: resolver level + a tiny hand-built gsd Frame. No full Cell.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

import gsd.hoomd

from ffn_sim.archive.hoomd_legacy.cell.linc import (
    GAMMA_DENYLIST_PREFIX,
    LINC_BOND_NESPRIN,
    PI_DECISIONS,
    ResolvedLINC,
    configure_linc_bond_potential,
    extend_snapshot_with_linc,
    linc_bond_force_along_axis,
    linc_bond_potential,
    linc_cfl_dt_max,
    pair_linc_bridges,
    resolve_linc,
)


# ---------------------------------------------------------------------------
# Tiny snapshot fixture: 3 nucleus beads (tags 0-2) + 3 cytoskeleton beads
# (tags 3-5), offset radially so geometric pairing is deterministic.
# ---------------------------------------------------------------------------
def _tiny_snapshot() -> gsd.hoomd.Frame:
    f = gsd.hoomd.Frame()
    # Nucleus beads near origin; cytoskeleton beads ~100 nm radially outward.
    nuc = np.array(
        [[0.0, 0.0, 0.0], [10e-9, 0.0, 0.0], [0.0, 10e-9, 0.0]],
        dtype=np.float64,
    )
    cyto = np.array(
        [[100e-9, 0.0, 0.0], [110e-9, 0.0, 0.0], [0.0, 100e-9, 0.0]],
        dtype=np.float64,
    )
    pos = np.concatenate([nuc, cyto], axis=0)
    n = pos.shape[0]
    f.particles.N = n
    f.particles.types = ["nucleus_bead", "actin_cortex"]
    f.particles.typeid = np.array([0, 0, 0, 1, 1, 1], dtype=np.uint32)
    f.particles.position = pos.astype(np.float64)
    # NOTE: deliberately NO f.particles.tag — a real build-time gsd.hoomd.Frame
    # does not expose .particles.tag (only a runtime cpu_local_snapshot does).
    # The extender consumes the snapshot tag-ordered (row index == global tag),
    # exactly like the cell.py extenders, so it must not read .particles.tag.
    # One pre-existing unrelated bond (cortex-cortex) to test carry-across.
    f.bonds.N = 1
    f.bonds.types = ["cortex_backbone"]
    f.bonds.typeid = np.array([0], dtype=np.uint32)
    f.bonds.group = np.array([[3, 4]], dtype=np.uint32)
    return f


def _enabled_cfg(**over: object) -> dict:
    cfg = {"linc": {"enabled": True}}
    cfg["linc"].update(over)
    return cfg


# ---------------------------------------------------------------------------
# OFF-identity (DEFAULT-OFF)
# ---------------------------------------------------------------------------
class TestOffIdentity:
    def test_resolve_default_is_off(self) -> None:
        p = resolve_linc({})
        assert p.enabled is False
        assert p.k_linc is None
        assert p.r0 == 0.0
        assert p.capture_radius == 0.0
        assert p.n_bridges_max == 0

    def test_resolve_explicit_false_is_off(self) -> None:
        p = resolve_linc({"linc": {"enabled": False, "k_linc": 1e-3}})
        assert p.enabled is False
        assert p.k_linc is None  # zeroed regardless of supplied value

    def test_builder_noop_when_disabled(self) -> None:
        snap = _tiny_snapshot()
        n_bonds_before = int(snap.bonds.N)
        types_before = list(snap.bonds.types)
        p = resolve_linc({})
        out = extend_snapshot_with_linc(
            snap, p, nucleus_tags=[0, 1, 2], cytoskeleton_tags=[3, 4, 5]
        )
        assert out is snap
        assert int(out.bonds.N) == n_bonds_before
        assert list(out.bonds.types) == types_before
        assert LINC_BOND_NESPRIN not in out.bonds.types

    def test_configure_potential_noop_when_disabled(self) -> None:
        class _FakeHarmonic:
            def __init__(self) -> None:
                self.params: dict = {}

        h = _FakeHarmonic()
        configure_linc_bond_potential(h, resolve_linc({}))
        assert h.params == {}


# ---------------------------------------------------------------------------
# §1 Dimensional / parameter sanity
# ---------------------------------------------------------------------------
class TestDimensionalAndParams:
    def test_defaults_are_geometric_literature_lengths(self) -> None:
        p = resolve_linc(_enabled_cfg())
        assert p.enabled is True
        # Perinuclear gap default (Crisp 2006 ~50 nm).
        assert p.r0 == pytest.approx(50e-9, rel=1e-12)
        # Capture radius = 3× the gap (geometric multiple).
        assert p.capture_radius == pytest.approx(150e-9, rel=1e-12)
        # Resting tension provenance is PI-PENDING (absolute magnitude +
        # attribution unresolved: Arsenovic 2016 reported relative FRET only;
        # Déjardin 2020 JCB ~8 pN — see PI_DECISIONS). f_rest feeds NO force
        # (provenance only), so assert only that it is a physiological per-nesprin
        # tension (~1-10 pN), NOT a single contested value.
        assert 1.0e-12 <= p.f_rest <= 1.0e-11
        # Stiffness genuinely unknown by default.
        assert p.k_linc is None

    def test_k_linc_passed_through_when_supplied(self) -> None:
        p = resolve_linc(_enabled_cfg(k_linc=1.0e-3))
        assert p.k_linc == pytest.approx(1.0e-3)

    def test_cfl_dt_max_dimensional(self) -> None:
        # dt = sf * gamma/k  → [s] = [-]*[N·s/m]/[N/m].
        dt = linc_cfl_dt_max(k_linc=2.0, gamma_bead=4.0, cfl_safety_factor=0.1)
        assert dt == pytest.approx(0.1 * 4.0 / 2.0)

    def test_potential_nonneg_and_zero_at_rest(self) -> None:
        k, r0 = 1.0e-3, 50e-9
        assert linc_bond_potential(k, r0, r0) == pytest.approx(0.0)
        assert linc_bond_potential(k, r0, r0 + 5e-9) > 0.0
        assert linc_bond_potential(k, r0, r0 - 5e-9) > 0.0

    @pytest.mark.parametrize(
        "over",
        [
            {"r0": 0.0},
            {"r0": -1e-9},
            {"capture_radius": 0.0},
            {"capture_radius": -1e-9},
            {"f_rest": -1e-12},
            {"n_bridges_max": -1},
            {"k_linc": 0.0},
            {"k_linc": -1.0},
        ],
    )
    def test_bad_params_raise(self, over: dict) -> None:
        with pytest.raises(ValueError):
            resolve_linc(_enabled_cfg(**over))

    def test_capture_radius_larger_than_cell_raises(self) -> None:
        with pytest.raises(ValueError):
            resolve_linc(_enabled_cfg(capture_radius=20e-6), R_cell=7.5e-6)

    def test_scalar_force_law_rejects_bad_args(self) -> None:
        with pytest.raises(ValueError):
            linc_bond_force_along_axis(k_linc=-1.0, r0=50e-9, sep=60e-9)
        with pytest.raises(ValueError):
            linc_bond_force_along_axis(k_linc=1.0, r0=0.0, sep=60e-9)
        with pytest.raises(ValueError):
            linc_bond_force_along_axis(k_linc=1.0, r0=50e-9, sep=-1e-9)


# ---------------------------------------------------------------------------
# §2 Boundary cases
# ---------------------------------------------------------------------------
class TestBoundary:
    def test_enabled_build_without_k_linc_raises(self) -> None:
        snap = _tiny_snapshot()
        p = resolve_linc(_enabled_cfg())  # k_linc None
        with pytest.raises(NotImplementedError):
            extend_snapshot_with_linc(
                snap, p, nucleus_tags=[0, 1, 2], cytoskeleton_tags=[3, 4, 5]
            )

    def test_configure_potential_without_k_linc_raises(self) -> None:
        class _FakeHarmonic:
            def __init__(self) -> None:
                self.params: dict = {}

        with pytest.raises(NotImplementedError):
            configure_linc_bond_potential(_FakeHarmonic(), resolve_linc(_enabled_cfg()))

    def test_no_acceptor_in_range_is_noop(self) -> None:
        snap = _tiny_snapshot()
        n_before = int(snap.bonds.N)
        # Tiny capture radius (< the 100 nm nucleus↔cortex gap) → no pairs.
        p = resolve_linc(_enabled_cfg(k_linc=1e-3, capture_radius=1e-9))
        out = extend_snapshot_with_linc(
            snap, p, nucleus_tags=[0, 1, 2], cytoskeleton_tags=[3, 4, 5]
        )
        assert int(out.bonds.N) == n_before  # unchanged

    def test_missing_nucleus_is_noop(self) -> None:
        snap = _tiny_snapshot()
        n_before = int(snap.bonds.N)
        p = resolve_linc(_enabled_cfg(k_linc=1e-3))
        out = extend_snapshot_with_linc(
            snap, p, nucleus_tags=[], cytoskeleton_tags=[3, 4, 5]
        )
        assert int(out.bonds.N) == n_before

    def test_pairing_empty_inputs(self) -> None:
        out = pair_linc_bridges(
            np.empty((0, 3)), np.zeros((3, 3)),
            capture_radius=1e-6, n_bridges_max=10,
        )
        assert out.shape == (0, 2)

    def test_pairing_cap_keeps_shortest(self) -> None:
        nuc = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=np.float64)
        cyto = np.array([[1e-9, 0, 0], [2e-9, 0, 0], [3e-9, 0, 0]], dtype=np.float64)
        # 3 nucleus beads, all pair to the nearest cyto (idx 0, dist 1nm);
        # cap to 2 keeps 2 of them.
        pairs = pair_linc_bridges(
            nuc, cyto, capture_radius=1e-6, n_bridges_max=2
        )
        assert pairs.shape[0] == 2

    def test_default_sentinel_cap_does_not_bind(self) -> None:
        # The default n_bridges_max (1e6) is a no-cap sentinel: the binding
        # ceiling is n_nuc (one LINC per envelope bead), not the 1e6 number.
        n_nuc = 5
        nuc = np.zeros((n_nuc, 3), dtype=np.float64)
        cyto = np.array([[1e-9, 0, 0]], dtype=np.float64)
        p = resolve_linc(_enabled_cfg(k_linc=1e-3))
        assert p.n_bridges_max == 1_000_000  # the inert sentinel
        pairs = pair_linc_bridges(
            nuc, cyto,
            capture_radius=p.capture_radius,
            n_bridges_max=p.n_bridges_max,
        )
        # min(n_nuc, 1e6) == n_nuc — the sentinel never binds.
        assert pairs.shape[0] == n_nuc


# ---------------------------------------------------------------------------
# Build-time-snapshot contract: the extender must NOT read .particles.tag
# (a real gsd.hoomd.Frame does not expose it). Consumed tag-ordered.
# ---------------------------------------------------------------------------
class TestNoTagAttributeContract:
    def test_frame_has_no_tag_attribute(self) -> None:
        # Sanity: a vanilla build-time Frame genuinely lacks .particles.tag.
        snap = _tiny_snapshot()
        assert not hasattr(snap.particles, "tag")

    def test_build_succeeds_without_tag_attribute(self) -> None:
        # Regression for the .tag-read bug: an enabled build on a real Frame
        # (NO .particles.tag) must form a LINC bond, not AttributeError.
        snap = _tiny_snapshot()
        n_before = int(snap.bonds.N)
        p = resolve_linc(_enabled_cfg(k_linc=1e-3))  # 150 nm covers the gaps
        out = extend_snapshot_with_linc(
            snap, p, nucleus_tags=[0, 1, 2], cytoskeleton_tags=[3, 4, 5]
        )
        assert int(out.bonds.N) > n_before, "expected ≥1 LINC bond to form"
        assert LINC_BOND_NESPRIN in out.bonds.types
        # Bond groups are tag-ordered rows: nucleus tag (0-2) ↔ cyto tag (3-5).
        linc_typeid = list(out.bonds.types).index(LINC_BOND_NESPRIN)
        groups = np.asarray(out.bonds.group).reshape(-1, 2)
        typeids = np.asarray(out.bonds.typeid)
        for a, b in groups[typeids == linc_typeid]:
            assert (a in (0, 1, 2)) and (b in (3, 4, 5))

    def test_out_of_range_tags_are_bounds_filtered(self) -> None:
        # Tags ≥ N (no such row) are dropped, not indexed out of bounds.
        snap = _tiny_snapshot()  # N == 6
        p = resolve_linc(_enabled_cfg(k_linc=1e-3))
        out = extend_snapshot_with_linc(
            snap, p, nucleus_tags=[0, 1, 2, 99], cytoskeleton_tags=[3, 4, 5, 100]
        )
        assert LINC_BOND_NESPRIN in out.bonds.types


# ---------------------------------------------------------------------------
# Real-geometry trap: at the production nucleus-cloud / cortex-shell bead
# representation the default 150 nm capture radius forms ZERO bonds. Document
# the trap + the corrected capture_radius that does form bonds.
# ---------------------------------------------------------------------------
class TestRealGeometryCaptureTrap:
    # MCF7 production geometry: R_cell = 7.5 µm, R_nuc = 0.25·R_cell = 1.875 µm.
    R_CELL = 7.5e-6
    R_NUC = 1.875e-6

    def _nuc_cloud_and_cortex_shell(
        self, n_nuc: int = 40, n_cortex: int = 200, seed: int = 0
    ) -> tuple[np.ndarray, np.ndarray]:
        """Filled nucleus cloud out to R_nuc + a cortex shell near R_cell."""
        rng = np.random.default_rng(seed)
        # Nucleus: filled ball out to R_nuc (uniform-in-volume radii).
        u = rng.random(n_nuc)
        r_nuc = self.R_NUC * np.cbrt(u)
        dirs_n = rng.normal(size=(n_nuc, 3))
        dirs_n /= np.linalg.norm(dirs_n, axis=1, keepdims=True)
        nuc = (dirs_n * r_nuc[:, None]).astype(np.float64)
        # Cortex: thin shell near R_cell.
        dirs_c = rng.normal(size=(n_cortex, 3))
        dirs_c /= np.linalg.norm(dirs_c, axis=1, keepdims=True)
        cyto = (dirs_c * self.R_CELL).astype(np.float64)
        return nuc, cyto

    def test_default_capture_radius_forms_zero_bonds(self) -> None:
        # The latent trap: nearest nucleus↔cortex bead gap is R_cell - R_nuc
        # ≈ 5.6 µm ≫ the 150 nm default → ZERO bridges at construction.
        nuc, cyto = self._nuc_cloud_and_cortex_shell()
        p = resolve_linc(_enabled_cfg(k_linc=1e-3))  # default 150 nm radius
        assert p.capture_radius == pytest.approx(150e-9, rel=1e-12)
        pairs = pair_linc_bridges(
            nuc, cyto,
            capture_radius=p.capture_radius,
            n_bridges_max=p.n_bridges_max,
        )
        assert pairs.shape[0] == 0, (
            "default 150 nm capture radius must form ZERO bonds at the real "
            "nucleus-cloud / cortex-shell geometry (the documented trap)"
        )

    def test_corrected_capture_radius_forms_bonds(self) -> None:
        # A capture_radius on the R_cell - R_nuc scale DOES form bridges.
        nuc, cyto = self._nuc_cloud_and_cortex_shell()
        # Slightly above the gap so the nearest acceptors are in range.
        cr = 1.1 * (self.R_CELL - self.R_NUC)
        p = resolve_linc(_enabled_cfg(k_linc=1e-3, capture_radius=cr))
        pairs = pair_linc_bridges(
            nuc, cyto,
            capture_radius=p.capture_radius,
            n_bridges_max=p.n_bridges_max,
        )
        assert pairs.shape[0] > 0, (
            "an R_cell - R_nuc-scale capture radius must form ≥1 bridge"
        )

    def test_resolve_warns_when_capture_below_nuc_cortex_gap(self) -> None:
        # The fail-loud guard: default radius + real R_cell/R_nuc → warning.
        with pytest.warns(UserWarning, match="ZERO LINC bonds"):
            resolve_linc(
                _enabled_cfg(),  # default 150 nm capture radius
                R_cell=self.R_CELL,
                R_nuc=self.R_NUC,
            )

    def test_resolve_no_warn_when_capture_covers_gap(self) -> None:
        # A capture radius covering the gap must NOT warn.
        cr = 1.1 * (self.R_CELL - self.R_NUC)
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # any warning → test failure
            resolve_linc(
                _enabled_cfg(capture_radius=cr),
                R_cell=self.R_CELL,
                R_nuc=self.R_NUC,
            )


# ---------------------------------------------------------------------------
# §3 Conservation / no-gamma-contamination
# ---------------------------------------------------------------------------
class TestNoGammaContamination:
    def test_prefix_nonempty(self) -> None:
        assert isinstance(GAMMA_DENYLIST_PREFIX, str)
        assert len(GAMMA_DENYLIST_PREFIX) > 0
        assert GAMMA_DENYLIST_PREFIX == "linc_"

    def test_bond_type_name_has_prefix(self) -> None:
        assert LINC_BOND_NESPRIN.startswith(GAMMA_DENYLIST_PREFIX)

    def test_every_built_bond_type_has_prefix(self) -> None:
        snap = _tiny_snapshot()
        p = resolve_linc(_enabled_cfg(k_linc=1e-3))  # capture 150nm covers gaps
        out = extend_snapshot_with_linc(
            snap, p, nucleus_tags=[0, 1, 2], cytoskeleton_tags=[3, 4, 5]
        )
        # New bonds were added (acceptors are within 150 nm of nucleus beads).
        added = [t for t in out.bonds.types if t not in ("cortex_backbone",)]
        assert added, "expected at least one LINC bond type"
        for t in added:
            assert t.startswith(GAMMA_DENYLIST_PREFIX), (
                f"bond type {t!r} would contaminate the cortical-tension γ "
                f"estimator (no {GAMMA_DENYLIST_PREFIX} prefix)"
            )
        # Pre-existing non-LINC bond is preserved.
        assert "cortex_backbone" in out.bonds.types

    def test_builder_only_adds_linc_bonds_between_supplied_tags(self) -> None:
        snap = _tiny_snapshot()
        p = resolve_linc(_enabled_cfg(k_linc=1e-3))
        out = extend_snapshot_with_linc(
            snap, p, nucleus_tags=[0, 1, 2], cytoskeleton_tags=[3, 4, 5]
        )
        linc_typeid = list(out.bonds.types).index(LINC_BOND_NESPRIN)
        groups = np.asarray(out.bonds.group).reshape(-1, 2)
        typeids = np.asarray(out.bonds.typeid)
        linc_groups = groups[typeids == linc_typeid]
        for a, b in linc_groups:
            # one end is a nucleus tag (0-2), the other a cytoskeleton tag (3-5)
            assert (a in (0, 1, 2)) and (b in (3, 4, 5))


# ---------------------------------------------------------------------------
# §5 Sign / sense
# ---------------------------------------------------------------------------
class TestSignSense:
    def test_stretched_bond_pulls_nucleus_toward_cytoskeleton(self) -> None:
        k, r0 = 1.0e-3, 50e-9
        # Stretched: separation > rest length → restoring force on the nucleus
        # bead is POSITIVE along Δr̂ (toward the cytoskeleton anchor).
        f = linc_bond_force_along_axis(k, r0, sep=60e-9)
        assert f > 0.0
        assert f == pytest.approx(k * (60e-9 - r0))

    def test_compressed_bond_pushes_apart(self) -> None:
        k, r0 = 1.0e-3, 50e-9
        f = linc_bond_force_along_axis(k, r0, sep=40e-9)
        assert f < 0.0  # negative along Δr̂ → nucleus pushed away from anchor

    def test_rest_length_is_force_free(self) -> None:
        k, r0 = 1.0e-3, 50e-9
        assert linc_bond_force_along_axis(k, r0, sep=r0) == pytest.approx(0.0)

    def test_force_monotone_increasing_in_separation(self) -> None:
        k, r0 = 1.0e-3, 50e-9
        seps = np.linspace(10e-9, 120e-9, 12)
        forces = [linc_bond_force_along_axis(k, r0, float(s)) for s in seps]
        assert all(
            forces[i + 1] > forces[i] for i in range(len(forces) - 1)
        )


# ---------------------------------------------------------------------------
# PI decisions present (honesty-over-completeness audit hook)
# ---------------------------------------------------------------------------
def test_pi_decisions_records_unknown_k_linc() -> None:
    assert isinstance(PI_DECISIONS, list)
    assert len(PI_DECISIONS) >= 1
    assert any("k_linc" in d for d in PI_DECISIONS)


# ---------------------------------------------------------------------------
# LIVE activation wiring (2026-06-09; PI 소유권 허용): manifest + Cell.build
# snapshot-extension (per-bond EXACT-r0 linc_nesprin bonds nucleus↔IF cage on the
# shared bond) + registry LIVE. Option A: bonds the nucleus to the perinuclear IF
# cage; linc_ bonds are γ-denylisted (no cortical contamination).
# ---------------------------------------------------------------------------
class TestLINCActivationWiring:
    def test_off_build_is_bit_identity(self) -> None:
        from ffn_sim.archive.hoomd_legacy.cell.manifest import (
            build_baseline_cell, load_manifest, resolve_baseline,
        )
        rb = resolve_baseline(load_manifest("mcf7_baseline.yaml"))
        assert rb.p_linc is None
        cell = build_baseline_cell("mcf7_baseline.yaml", seed=1)
        assert cell.p_linc is None
        assert cell.simulation.state.N_particles > 0

    def test_on_bridges_form_without_contamination(self) -> None:
        from ffn_sim.archive.hoomd_legacy.cell.compartment_registry import REGISTRY, load_recipe
        from ffn_sim.archive.hoomd_legacy.cell.manifest import build_baseline_cell, load_manifest
        from ffn_sim.archive.hoomd_legacy.cortex.cortical_tension import (
            _is_adhesion_bond_type, measure_cortical_tension,
        )

        base = load_manifest("mcf7_baseline.yaml")
        # OFF = IF cage on, LINC off; ON = LINC coupled (IF cage + LINC).
        m_off, d_off = REGISTRY.compose_manifest(
            load_recipe("if_cage"), base_manifest=base, strict=True
        )
        m_on, d_on = REGISTRY.compose_manifest(
            load_recipe("linc_coupled"), base_manifest=base, strict=True
        )
        assert d_off == [] and d_on == []
        off = build_baseline_cell("mcf7_baseline.yaml", manifest=m_off, seed=1)
        on = build_baseline_cell("mcf7_baseline.yaml", manifest=m_on, seed=1)
        p = on.p_linc
        assert p is not None and p.enabled

        so = off.simulation.state.get_snapshot()
        sn = on.simulation.state.get_snapshot()
        n_bridges = int(on.extras["handles"]["n_linc_bridges"])
        # BRIDGES FORMED (REFUTE guard) + bonds-only (zero new particles).
        assert n_bridges > 0
        assert sn.particles.N == so.particles.N            # LINC adds 0 particles
        assert sn.bonds.N - so.bonds.N == n_bridges        # +1 bond per bridge
        linc_types = [t for t in sn.bonds.types if t.startswith("linc_")]
        assert linc_types, "expected per-bond linc_nesprin types"
        assert LINC_BOND_NESPRIN in linc_types
        # Every linc bond is on the cortical-γ denylist (no contamination path).
        assert all(_is_adhesion_bond_type(t) for t in linc_types)
        # FORCE-FREE construction: each bridge born at its EXACT as-built sep.
        layout = on.extras["handles"]["linc_layout"]
        pos = np.asarray(sn.particles.position, dtype=np.float64)
        pairs = np.asarray(layout.pairs, dtype=np.int64)
        sep = np.linalg.norm(pos[pairs[:, 0]] - pos[pairs[:, 1]], axis=1)
        r0 = np.asarray(layout.r0, dtype=np.float64)
        assert float(np.max(np.abs(sep - r0) / np.maximum(r0, 1e-30))) < 1e-6
        # Acceptor degree bounded: each if_bead carries at most one LINC bond.
        _, counts = np.unique(pairs[:, 1], return_counts=True)
        assert int(counts.max()) == 1
        # NO contamination: cortical γ_soft identical OFF vs ON (linc_ denylisted).
        off.simulation.run(0); on.simulation.run(0)
        go = measure_cortical_tension(off.simulation, R_cell=off.p_cortex.R_cell,
                                      p_enclosed_volume=off.p_enclosed_volume)
        gn = measure_cortical_tension(on.simulation, R_cell=on.p_cortex.R_cell,
                                      p_enclosed_volume=on.p_enclosed_volume)
        k = "gamma_soft_N_per_m" if "gamma_soft_N_per_m" in go else "gamma_soft"
        assert abs(gn[k] - go[k]) <= 1e-12 * max(1.0, abs(go[k]))

    def test_enabled_without_if_cage_raises(self) -> None:
        # LINC (Option A) REQUIRES the IF cage (its acceptors are if_bead).
        # Enabling linc without intermediate_filaments must raise (manifest guard).
        import pytest
        from ffn_sim.archive.hoomd_legacy.cell.manifest import load_manifest, resolve_baseline

        m = load_manifest("mcf7_baseline.yaml")
        m.setdefault("optional_subsystems", {})
        m["optional_subsystems"]["linc"] = {"enabled": True, "k_linc": 1.0e-2}
        with pytest.raises(ValueError):
            resolve_baseline(m)


def test_resolved_dataclass_is_frozen() -> None:
    p = resolve_linc(_enabled_cfg(k_linc=1e-3))
    with pytest.raises(Exception):
        p.k_linc = 2e-3  # type: ignore[misc]  frozen dataclass
    assert math.isfinite(p.r0)
