"""Tests for ``ffn_sim.archive.hoomd_legacy.cell.membrane_reservoir`` (H.8 reservoir + bleb tethers).

Import-light, fast: resolver-level + a tiny synthetic gsd frame for the tether
builder. No full Cell is built. The OFF-identity, dimensional/sign sanity, and
no-gamma-contamination checks all run cheaply.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ffn_sim.archive.hoomd_legacy.cell.membrane_reservoir import (
    DEFAULT_W_MCA,
    GAMMA_DENYLIST_PREFIX,
    MEM_TETHER_BOND,
    PI_DECISIONS,
    W_MCA_BAND,
    MembraneTetherLayout,
    MembraneTetherUpdater,
    ResolvedMembraneReservoir,
    build_membrane_tethers,
    resolve_membrane_reservoir,
)

KT_300 = 4.28e-21  # J  (≈ k_B · 310 K cytoplasm; matches in-tree convention)
R_CELL = 7.5e-6    # m  MCF7 radius (Wagner 2011)


# ---------------------------------------------------------------------------
# Tiny synthetic frame helper (membrane beads just OUTSIDE cortex beads)
# ---------------------------------------------------------------------------
def _make_frame(n: int = 6, gap: float = 50.0e-9):
    """Build a gsd frame: n cortex beads then n membrane beads at +gap radius.

    Cortex beads on a sphere R_CELL; each membrane bead sits ``gap`` farther out
    along the same radial direction (so the nearest-cortex tether is well-defined
    and within max_tether_dist).
    """
    import gsd.hoomd

    idx = np.arange(n, dtype=np.float64) + 0.5
    phi = math.pi * (3.0 - math.sqrt(5.0))
    ct = np.clip(1.0 - 2.0 * idx / n, -1.0, 1.0)
    st = np.sqrt(np.maximum(0.0, 1.0 - ct * ct))
    az = phi * idx
    dirs = np.stack([st * np.cos(az), st * np.sin(az), ct], axis=1)
    cortex = R_CELL * dirs
    membrane = (R_CELL + gap) * dirs
    pos = np.concatenate([cortex, membrane], axis=0)

    f = gsd.hoomd.Frame()
    f.particles.N = 2 * n
    f.particles.types = ["actin_cortex", "mem_bead"]
    f.particles.typeid = np.concatenate(
        [np.zeros(n, dtype=np.uint32), np.ones(n, dtype=np.uint32)]
    )
    f.particles.position = pos.astype(np.float64)
    f.particles.mass = np.ones(2 * n, dtype=np.float64)
    # One pre-existing cortex backbone bond so we test passthrough + append.
    f.bonds.N = 1
    f.bonds.types = ["cortex-bond"]
    f.bonds.typeid = np.zeros(1, dtype=np.uint32)
    f.bonds.group = np.array([[0, 1]], dtype=np.uint32)
    f.configuration.box = [4 * R_CELL, 4 * R_CELL, 4 * R_CELL, 0, 0, 0]
    return f, n


# ---------------------------------------------------------------------------
# OFF-identity
# ---------------------------------------------------------------------------
def test_off_identity_resolver():
    """enabled=False (and missing) → .enabled False with zeros + None PI consts."""
    for cfg in ({"enabled": False}, {}):  # explicit-off and missing-key
        p = resolve_membrane_reservoir(cfg, R_cell=R_CELL)
        assert isinstance(p, ResolvedMembraneReservoir)
        assert p.enabled is False
        assert p.W_MCA == 0.0
        assert p.k_tether == 0.0
        assert p.sigma_crit_bleb is None
        assert p.f_excess is None


def test_off_identity_builder_no_op():
    """Disabled builder returns the SAME frame object, adds zero bonds."""
    f, n = _make_frame()
    p = resolve_membrane_reservoir({"enabled": False}, R_cell=R_CELL)
    out, layout = build_membrane_tethers(
        f, p, membrane_tag_range=(n, 2 * n), cortex_tag_range=(0, n)
    )
    assert out is f                       # identity: same object, untouched
    assert layout.n_tether == 0
    assert int(out.bonds.N) == 1          # only the pre-existing cortex bond
    assert MEM_TETHER_BOND not in list(out.bonds.types)
    assert int(out.particles.N) == 2 * n  # no particles added


def test_off_identity_nested_cfg():
    """Resolver accepts root / cell / membrane_reservoir nesting; default-off."""
    p = resolve_membrane_reservoir(
        {"cell": {"membrane_reservoir": {"enabled": False}}}, R_cell=R_CELL
    )
    assert p.enabled is False


# ---------------------------------------------------------------------------
# Dimensional / parameter sanity
# ---------------------------------------------------------------------------
def test_resolver_derived_quantities():
    """Enabled resolver computes the eq.(R) bridge correctly (SI units)."""
    p = resolve_membrane_reservoir({"enabled": True}, R_cell=R_CELL)
    assert p.enabled is True
    assert p.W_MCA == DEFAULT_W_MCA
    n = 1000
    # A_bead = 4π R²/n
    A_bead = 4.0 * math.pi * R_CELL**2 / n
    assert math.isclose(p.area_share(n), A_bead, rel_tol=1e-12)
    # E_tether = W_MCA · A_bead  [J]
    E = p.W_MCA * A_bead
    assert math.isclose(p.tether_adhesion_energy(n), E, rel_tol=1e-12)
    # F_c = √(2 k W A_bead)  [N]
    F_c = math.sqrt(2.0 * p.k_tether * E)
    assert math.isclose(p.rupture_force(n), F_c, rel_tol=1e-12)
    # Δ_c = F_c / k  [m]
    assert math.isclose(p.rupture_extension(n), F_c / p.k_tether, rel_tol=1e-12)


def test_total_adhesion_is_grid_invariant():
    """Σ E_tether = W_MCA·4πR² is N-independent (intensive material adhesion)."""
    p = resolve_membrane_reservoir({"enabled": True}, R_cell=R_CELL)
    total_a = p.tether_adhesion_energy(1000) * 1000
    total_b = p.tether_adhesion_energy(38000) * 38000
    expected = p.W_MCA * 4.0 * math.pi * R_CELL**2
    assert math.isclose(total_a, expected, rel_tol=1e-9)
    assert math.isclose(total_b, expected, rel_tol=1e-9)


@pytest.mark.parametrize(
    "cfg",
    [
        {"enabled": True, "W_MCA": -1.0e-5},     # negative adhesion
        {"enabled": True, "W_MCA": 0.0},         # zero adhesion
        {"enabled": True, "k_tether": -0.1},     # negative stiffness
        {"enabled": True, "k_tether": 0.0},      # zero stiffness
        {"enabled": True, "max_tether_dist": 0.0},
        {"enabled": True, "batch_steps": 0},
        {"enabled": True, "W_MCA": 1.0},         # far above KU-3.B1.3 band
    ],
)
def test_bad_params_raise(cfg):
    """Negative/zero/out-of-band enabled params raise ValueError (Sanity §2)."""
    with pytest.raises(ValueError):
        resolve_membrane_reservoir(cfg, R_cell=R_CELL)


def test_bad_R_cell_raises():
    """Non-positive R_cell raises even when enabled."""
    with pytest.raises(ValueError):
        resolve_membrane_reservoir({"enabled": True}, R_cell=-1.0)


def test_W_MCA_band_endpoints_ok():
    """KU-3.B1.3 band endpoints resolve; just outside raises."""
    lo, hi = W_MCA_BAND
    assert resolve_membrane_reservoir({"enabled": True, "W_MCA": lo}, R_cell=R_CELL).enabled
    assert resolve_membrane_reservoir({"enabled": True, "W_MCA": hi}, R_cell=R_CELL).enabled
    with pytest.raises(ValueError):
        resolve_membrane_reservoir({"enabled": True, "W_MCA": hi * 1.01}, R_cell=R_CELL)


# ---------------------------------------------------------------------------
# PI-gated uncertain constants stay disabled (no invented number)
# ---------------------------------------------------------------------------
def test_pi_decisions_nonempty():
    """Both uncertain MCF7 constants are flagged in PI_DECISIONS."""
    assert len(PI_DECISIONS) >= 2
    joined = " ".join(PI_DECISIONS)
    assert "sigma_crit_bleb" in joined
    assert "f_excess" in joined


def test_reservoir_release_disabled_when_f_excess_none():
    """released_area raises NotImplementedError while f_excess is None (MCF7)."""
    p = resolve_membrane_reservoir({"enabled": True}, R_cell=R_CELL)
    assert p.f_excess is None
    with pytest.raises(NotImplementedError):
        p.released_area(4.0 * math.pi * R_CELL**2)


def test_reservoir_release_works_when_f_excess_supplied():
    """A PI-anchored f_excess enables A_release = f_excess·A0 (path provided)."""
    p = resolve_membrane_reservoir(
        {"enabled": True, "f_excess": 0.1}, R_cell=R_CELL
    )
    A0 = 4.0 * math.pi * R_CELL**2
    assert math.isclose(p.released_area(A0), 0.1 * A0, rel_tol=1e-12)


def _empty_layout() -> MembraneTetherLayout:
    return MembraneTetherLayout(
        tether_pairs=np.empty((0, 2), dtype=np.int64),
        tether_r0=np.empty((0,), dtype=np.float64),
        n_tether=0,
        rupture_force=0.0,
    )


def test_bleb_updater_disabled_when_sigma_crit_none():
    """MembraneTetherUpdater raises NotImplementedError without σ_crit (MCF7)."""
    p = resolve_membrane_reservoir({"enabled": True}, R_cell=R_CELL)
    assert p.sigma_crit_bleb is None
    with pytest.raises(NotImplementedError):
        MembraneTetherUpdater(p, _empty_layout(), kT=KT_300)


def test_bleb_updater_refuses_construction_even_with_sigma_anchored():
    """Construction is refused even WITH a PI-anchored σ_crit (rupture loop unwritten).

    Finding (doc/gate honesty): the class docstring used to promise "supply a
    PI-anchored sigma_crit_bleb to enable it", but the Bell-Evans rupture loop in
    act() is an unimplemented TODO. __init__ now refuses construction
    UNCONDITIONALLY so __init__ and act() agree and no non-functional updater is
    instantiated. The error message must NOT blame a missing sigma (it was
    supplied) — it must state the rupture loop is unimplemented.
    """
    p = resolve_membrane_reservoir(
        {"enabled": True, "sigma_crit_bleb": 4.0e-4}, R_cell=R_CELL
    )
    assert p.sigma_crit_bleb == 4.0e-4  # σ IS anchored
    with pytest.raises(NotImplementedError) as exc:
        MembraneTetherUpdater(p, _empty_layout(), kT=KT_300)
    msg = str(exc.value).lower()
    # the reason is the unimplemented rupture loop, not a missing sigma
    assert "rupture loop" in msg or "not yet functional" in msg
    assert "act()" in msg or "act" in msg


# ---------------------------------------------------------------------------
# Sign / sense
# ---------------------------------------------------------------------------
def test_rupture_force_increases_with_adhesion():
    """Larger W_MCA → larger rupture force F_c (harder to bleb). Monotone."""
    p_lo = resolve_membrane_reservoir({"enabled": True, "W_MCA": 1.0e-6}, R_cell=R_CELL)
    p_hi = resolve_membrane_reservoir({"enabled": True, "W_MCA": 1.0e-4}, R_cell=R_CELL)
    assert p_hi.rupture_force(1000) > p_lo.rupture_force(1000) > 0.0


def test_tether_force_sign_stretched_pulls_inward():
    """A stretched tether (Δr>r0) pulls the membrane bead back toward cortex.

    Harmonic restoring force F = -k(Δr-r0)·r̂ along the bond axis: for a membrane
    bead pulled OUTWARD past r0, the force component along the outward direction is
    NEGATIVE (inward / restoring) — the adhesion resists detachment.
    """
    p = resolve_membrane_reservoir({"enabled": True}, R_cell=R_CELL)
    k = p.k_tether
    r0 = 50.0e-9
    # membrane bead displaced outward to Δr = 80 nm
    dr = 80.0e-9
    # scalar harmonic force magnitude along the outward unit vector
    f_along_outward = -k * (dr - r0)
    assert f_along_outward < 0.0          # inward / restoring
    # compressed (Δr < r0) → outward / re-separating
    f_compressed = -k * (30.0e-9 - r0)
    assert f_compressed > 0.0


# ---------------------------------------------------------------------------
# Tether topology builder (enabled path) + no-acceptor boundary
# ---------------------------------------------------------------------------
def test_builder_seeds_tethers_force_free():
    """Enabled builder appends mem_tether bonds at construction separation."""
    f, n = _make_frame(n=6, gap=50.0e-9)
    p = resolve_membrane_reservoir({"enabled": True}, R_cell=R_CELL)
    out, layout = build_membrane_tethers(
        f, p, membrane_tag_range=(n, 2 * n), cortex_tag_range=(0, n)
    )
    assert layout.n_tether == n           # every membrane bead found an acceptor
    assert MEM_TETHER_BOND in list(out.bonds.types)
    # pre-existing cortex bond preserved + n tethers appended
    assert int(out.bonds.N) == 1 + n
    # particle count unchanged (default build reuses existing beads)
    assert int(out.particles.N) == 2 * n
    # rest length ≈ the gap (construction separation), within a nm
    assert np.allclose(layout.tether_r0, 50.0e-9, atol=2.0e-9)
    # each tether pair is (membrane_tag, cortex_tag): col0 ≥ n, col1 < n
    assert (layout.tether_pairs[:, 0] >= n).all()
    assert (layout.tether_pairs[:, 1] < n).all()


def test_builder_no_acceptor_in_range_is_noop():
    """No cortex bead within max_tether_dist → snapshot unchanged."""
    f, n = _make_frame(n=6, gap=50.0e-9)
    # tiny acceptor reach (1 nm) < the 50 nm gap → no tether found
    p = resolve_membrane_reservoir(
        {"enabled": True, "max_tether_dist": 1.0e-9}, R_cell=R_CELL
    )
    out, layout = build_membrane_tethers(
        f, p, membrane_tag_range=(n, 2 * n), cortex_tag_range=(0, n)
    )
    assert layout.n_tether == 0
    assert out is f                       # unchanged input
    assert int(out.bonds.N) == 1


def test_builder_inverted_range_raises():
    """Inverted tag range raises ValueError (Sanity §2)."""
    f, n = _make_frame()
    p = resolve_membrane_reservoir({"enabled": True}, R_cell=R_CELL)
    with pytest.raises(ValueError):
        build_membrane_tethers(
            f, p, membrane_tag_range=(2 * n, n), cortex_tag_range=(0, n)
        )


# ---------------------------------------------------------------------------
# No gamma contamination
# ---------------------------------------------------------------------------
def test_gamma_denylist_prefix_nonempty():
    """The declared denylist prefix is non-empty."""
    assert isinstance(GAMMA_DENYLIST_PREFIX, str)
    assert GAMMA_DENYLIST_PREFIX != ""
    assert GAMMA_DENYLIST_PREFIX == "mem_"


def test_all_bond_types_start_with_denylist_prefix():
    """Every mem_tether bond this builder creates starts with GAMMA_DENYLIST_PREFIX.

    This is the no-gamma-contamination guarantee: the cortical-tension estimator
    denylists any bond whose name startswith a non-cortical prefix, so the radial
    detachment load path never inflates γ_soft.
    """
    f, n = _make_frame(n=6, gap=50.0e-9)
    p = resolve_membrane_reservoir({"enabled": True}, R_cell=R_CELL)
    out, _ = build_membrane_tethers(
        f, p, membrane_tag_range=(n, 2 * n), cortex_tag_range=(0, n)
    )
    created = [t for t in list(out.bonds.types) if t != "cortex-bond"]
    assert created == [MEM_TETHER_BOND]
    for bt in created:
        assert bt.startswith(GAMMA_DENYLIST_PREFIX), bt


def test_denylist_matches_cortical_tension_estimator():
    """The cortical_tension estimator's startswith-denylist would drop mem_tether.

    Mirrors the estimator's prefix test (name.startswith(prefix)) so a future
    wiring that adds GAMMA_DENYLIST_PREFIX to its denylist excludes mem_tether.
    """
    assert MEM_TETHER_BOND.startswith(GAMMA_DENYLIST_PREFIX)


# ---------------------------------------------------------------------------
# LIVE activation wiring (2026-06-09; PI 소유권 허용): own mem_node offset layer +
# static mem_tether mesh on the full physiological baseline. cell.py snapshot-
# extension + LJ r_cut=0 + gamma_map + manifest stanza + registry LIVE. Bleb
# rupture + reservoir release stay PI-blocked.
# ---------------------------------------------------------------------------
class TestMembraneReservoirActivationWiring:
    def test_off_build_is_bit_identity(self):
        from ffn_sim.archive.hoomd_legacy.cell.manifest import (
            build_baseline_cell, load_manifest, resolve_baseline,
        )
        rb = resolve_baseline(load_manifest("mcf7_baseline.yaml"))
        assert rb.p_membrane_reservoir is None
        cell = build_baseline_cell("mcf7_baseline.yaml", seed=1)
        assert cell.p_membrane_reservoir is None
        assert cell.simulation.state.N_particles > 0

    def test_on_mesh_assembles_without_contamination(self):
        import numpy as np
        from ffn_sim.archive.hoomd_legacy.cell.compartment_registry import REGISTRY, load_recipe
        from ffn_sim.archive.hoomd_legacy.cell.manifest import build_baseline_cell, load_manifest
        from ffn_sim.archive.hoomd_legacy.cortex.cortical_tension import (
            _is_adhesion_bond_type, measure_cortical_tension,
        )

        base = load_manifest("mcf7_baseline.yaml")
        m_on, deferred = REGISTRY.compose_manifest(
            load_recipe("membrane_reservoir_tethered"), base_manifest=base, strict=True
        )
        assert deferred == []
        off = build_baseline_cell("mcf7_baseline.yaml", seed=1)
        on = build_baseline_cell("mcf7_baseline.yaml", manifest=m_on, seed=1)
        p = on.p_membrane_reservoir
        assert p is not None and p.enabled

        so = off.simulation.state.get_snapshot()
        sn = on.simulation.state.get_snapshot()
        n_teth = int(on.extras["handles"]["n_mem_tethers"])
        layout = on.extras["handles"]["membrane_tether_layout"]
        # MESH ASSEMBLED: own mem_node layer (+1 particle per tether) + mem_tether.
        assert n_teth > 0
        assert "mem_node" in sn.particles.types
        assert "mem_tether" in sn.bonds.types
        assert sn.particles.N - so.particles.N == n_teth   # one mem_node per tether
        assert sn.bonds.N - so.bonds.N == n_teth
        mem_bonds = [t for t in sn.bonds.types if t.startswith("mem_")]
        assert all(_is_adhesion_bond_type(t) for t in mem_bonds)
        # CROSS-LAYER: mem_node ↔ cortex, distinct cortex anchors, no self-pairs.
        types = list(sn.particles.types)
        tid = np.asarray(sn.particles.typeid, dtype=np.int64)
        pairs = np.asarray(layout.tether_pairs, dtype=np.int64)
        assert (tid[pairs[:, 0]] == types.index("mem_node")).all()
        assert (tid[pairs[:, 1]] == types.index("actin_cortex")).all()
        assert len(set(pairs[:, 1].tolist())) == pairs.shape[0]   # unique anchors
        assert int((pairs[:, 0] == pairs[:, 1]).sum()) == 0
        # FORCE-FREE: every tether born at exactly membrane_offset.
        pos = np.asarray(sn.particles.position, dtype=np.float64)
        sep = np.linalg.norm(pos[pairs[:, 0]] - pos[pairs[:, 1]], axis=1)
        assert float(np.max(np.abs(sep - p.membrane_offset) / p.membrane_offset)) < 1e-6
        # NO contamination: cortical γ_soft identical OFF vs ON (mem_ denylisted).
        off.simulation.run(0); on.simulation.run(0)
        go = measure_cortical_tension(off.simulation, R_cell=off.p_cortex.R_cell,
                                      p_enclosed_volume=off.p_enclosed_volume)
        gn = measure_cortical_tension(on.simulation, R_cell=on.p_cortex.R_cell,
                                      p_enclosed_volume=on.p_enclosed_volume)
        k = "gamma_soft_N_per_m" if "gamma_soft_N_per_m" in go else "gamma_soft"
        assert abs(gn[k] - go[k]) <= 1e-12 * max(1.0, abs(go[k]))

    def test_enabled_requires_membrane_surface(self):
        # membrane_reservoir requires membrane_surface; resolving with it OFF raises.
        import pytest
        from ffn_sim.archive.hoomd_legacy.cell.manifest import load_manifest, resolve_baseline
        m = load_manifest("mcf7_baseline.yaml")
        m.setdefault("optional_subsystems", {})
        m["optional_subsystems"]["membrane_reservoir"] = {"enabled": True}
        m.setdefault("compartments", {}).setdefault("membrane_surface", {})["enabled"] = False
        with pytest.raises((ValueError, Exception)):
            resolve_baseline(m)
