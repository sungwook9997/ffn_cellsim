"""H.4 -> Cell focal-adhesion integration (alpha restart S0/S1/S2) tests.

Covers the FIRST vertical slice of FA wiring into
``build_cortex_full_simulation`` / ``Cell.build``:

* S0 substrate anchor -- a fixed/immobile ligand layer at z=0
  (``ligand`` particles held by :class:`ffn_sim.cell.cell.SubstrateLigandPin`).
* S1 integrin-ligand catch bond -- the ``integrin_ligand`` bond TYPE is
  registered and an :class:`ffn_sim.bridge.integrin_bonds.IntegrinBondUpdater`
  is attached (Pereverzev catch-slip, reused as-is; the PI-gated
  Pereverzev->Kong migration is NOT done here).
* S2 ``fa_actin_clutch`` bond -- the NEW load-path topology connecting an
  integrin to its nearest cortex-actin bead, force-free at construction
  (per-bond exact-r0).

Contract (CLAUDE.md, H4_FA_INTEGRATION_DESIGN.md):

* ADDITIVE + DEFAULT-OFF -- with ``p_fa=None`` the builder is bit-for-bit
  the pre-FA builder. :func:`test_fa_off_is_noop` asserts the particle /
  bond counts + type lists + position checksums are identical to a bare
  cortex+myosin build.
* The clutch load path is an explicit HOOMD bond, never a scalar.

Two honest, documented findings pinned by these tests
-----------------------------------------------------
1. **S5 tag-space unification -> FA coexists with myosin/xlink/lamellipodium.**
   The FA block (integrins then substrate ligands) is APPENDED LAST, in
   natural append order after cortex / myosin / xlink / lamellipodium, so
   NOTHING already present is re-tagged and every other subsystem's
   absolute-tag Updater bookkeeping stays valid. The reused
   ``IntegrinBondUpdater`` was generalized to resolve per-integrin state by
   an explicit global-tag->local-index map (built from
   ``integrin_tag_start``) and to read positions by global tag (a HOOMD
   snapshot is tag-ordered), so integrins no longer need to occupy
   ``[0, n_int)``. The previous ``NotImplementedError`` for FA +
   (myosin|xlink|lamellipodium) is removed.
   :func:`test_fa_with_myosin_builds` (+ xlink / lamellipodium smoke tests)
   pin the coexistence.

2. **Physical capture radius barely clutches.** At the real cell geometry
   the substrate plane (z~0) and the cortex shell (r~R_cell ~10 um) are
   spatially disjoint, so the physical ``capture_radius_R_FA`` (1.5 um)
   reaches almost no cortex beads -- the literal load path is essentially
   absent until the cell is equilibrated / positioned onto the substrate
   (S5/B2). :func:`test_fa_on_builds` therefore passes an explicit
   ``fa_clutch_capture_radius`` spanning the gap to exercise the real
   ``fa_actin_clutch`` topology (the governed ``resolve_h4`` default is
   unchanged), and
   :func:`test_fa_physical_capture_radius_barely_clutches` pins the
   honest near-zero-clutch finding.

Stability scope
---------------
:func:`test_fa_short_run_no_nan` runs a SHORT BAOAB warm-up (a few hundred
steps at the cortex CFL ``dt``) -- NOT production. It checks construction +
warm-up: no NaN/Inf, substrate ligands stay at z=0. Full production
stability needs the PI-gated equilibration prelude (design B2) + the global
dt reconciliation (design B1); those are out of scope for this slice.
"""

from __future__ import annotations

import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import resolve_cortex_myosin
from ffn_sim.bridge.fa import resolve_h4
from ffn_sim.cell.cell import (
    Cell,
    CellBuildOptions,
    build_cortex_full_simulation,
)

CONFIG_DIR = Path(__file__).resolve().parents[1] / "configs"

# Multi-minute / production-scale checks are opt-in (mirrors the existing
# H3_*_PRODUCTION gate convention).
RUN_PRODUCTION = os.environ.get("H4_FA_PRODUCTION", "0") == "1"


# ---------------------------------------------------------------------------
# Config fixtures -- CI-friendly demo scale (mirrors test_cell_full.py)
# ---------------------------------------------------------------------------
def _cortex_cfg() -> dict:
    cfg = deepcopy(yaml.safe_load(open(CONFIG_DIR / "phase1_h3.yaml")))
    cfg["cortex"]["n_filaments"] = 40
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = 5
    return cfg


@pytest.fixture(scope="module")
def cortex_cfg() -> dict:
    return _cortex_cfg()


@pytest.fixture(scope="module")
def p_cortex(cortex_cfg):
    return resolve_h3_derived(cortex_cfg)


@pytest.fixture(scope="module")
def p_myosin(cortex_cfg, p_cortex):
    return resolve_cortex_myosin(cortex_cfg, dt=p_cortex.dt_cfl)


@pytest.fixture(scope="module")
def p_fa(p_cortex):
    """Small FA demo resolved on configs/phase1_h4.yaml.

    Box matched to the cortex box so integrins/ligands land inside the same
    periodic cell as the cortex shell.
    """
    cfg = deepcopy(yaml.safe_load(open(CONFIG_DIR / "phase1_h4.yaml")))
    cfg["bridge"]["fa"]["n_nascent_per_cell"] = 2
    cfg["bridge"]["fa"]["n_mature_per_cell"] = 1
    cfg["bridge"]["fa"]["n_total_per_fa"] = 4
    cfg["bridge"]["fa"]["L_box"] = float(p_cortex.L_box)
    return resolve_h4(cfg)


# ---------------------------------------------------------------------------
# S0/S1/S2 -- ADDITIVE + DEFAULT-OFF
# ---------------------------------------------------------------------------
def _fingerprint(snap) -> dict:
    pos = np.asarray(snap.particles.position)
    return {
        "N": int(snap.particles.N),
        "ptypes": list(snap.particles.types),
        "bonds_N": int(snap.bonds.N),
        "btypes": list(snap.bonds.types),
        "angles_N": int(snap.angles.N),
        "pos_sum": float(pos.sum()),
        "pos_sq": float((pos * pos).sum()),
    }


def test_fa_off_is_noop(p_cortex, p_myosin):
    """``p_fa=None`` => build is bit-for-bit the pre-FA cortex+myosin build.

    Asserts identical particle / bond counts + TYPE lists + position
    checksums against a bare cortex+myosin build, and that all FA handles
    are None / 0 (no FA Updaters attached, no FA particle / bond types leak
    into the state).
    """
    h_base = build_cortex_full_simulation(p_cortex, p_myosin=p_myosin)
    h_off = build_cortex_full_simulation(p_cortex, p_myosin=p_myosin, p_fa=None)

    fp_base = _fingerprint(h_base["sim"].state.get_snapshot())
    fp_off = _fingerprint(h_off["sim"].state.get_snapshot())
    assert fp_off == fp_base, (
        "FA-off build diverged from the bare cortex+myosin build:\n"
        f"  base={fp_base}\n  off ={fp_off}"
    )

    # No FA particle / bond types registered.
    assert "integrin" not in fp_off["ptypes"]
    assert "ligand" not in fp_off["ptypes"]
    assert "integrin_ligand" not in fp_off["btypes"]
    assert "fa_actin_clutch" not in fp_off["btypes"]

    # FA handles inert.
    for key in (
        "fa_integration", "fa_layout", "integrin_action", "integrin_updater",
        "ligand_pin_action", "ligand_pin_updater",
    ):
        assert h_off[key] is None, f"{key} should be None when p_fa=None"
    assert h_off["n_fa_integrins"] == 0
    assert h_off["n_substrate_ligands"] == 0
    assert h_off["n_fa_clutch_bonds"] == 0

    # Same number of attached updaters (no extra FA updaters).
    assert (
        len(h_off["sim"].operations.updaters)
        == len(h_base["sim"].operations.updaters)
    )


def test_cell_build_fa_off_is_noop(p_cortex, p_myosin):
    """Cell.build with with_fa=False leaves the FA slot dead (regression).

    The ``fa_integrin`` tag slot stays zero-width and ``cell.fa`` is None,
    exactly as the pre-FA Cell reported.
    """
    cell = Cell.build(
        p_cortex, p_myosin=p_myosin,
        options=CellBuildOptions(with_myosin=True),
    )
    assert cell.fa is None
    tr = cell.tag_ranges()
    lo, hi = tr["fa_integrin"]
    assert lo == hi, "fa_integrin slot must be zero-width when FA off"
    assert "substrate_ligand" not in tr
    assert sum(cell.bead_count_summary().values()) == cell.n_particles


# ---------------------------------------------------------------------------
# S5 tag-space unification -- FA + (myosin/xlink/lamellipodium) now COEXIST
# ---------------------------------------------------------------------------
def test_fa_with_myosin_builds(p_cortex, p_myosin, p_fa):
    """S5: FA + myosin coexist (was NotImplementedError before tag-space
    unification).

    Replaces the old ``test_fa_with_myosin_raises``. The FA block is appended
    LAST so the myosin Updater's absolute cortex-actin tag bookkeeping is
    untouched, and the generalized IntegrinBondUpdater resolves integrin state
    by an explicit tag->local map. Asserts:

    * all three bond families present (integrin_ligand + fa_actin_clutch +
      cortex_myosin_*),
    * clutch bonds connect integrins to CORTEX-ACTIN tags (NOT myosin tags),
    * a short BAOAB warm-up runs with no NaN/Inf and substrate ligands stay
      at z=0.

    Full production stability still needs the PI-gated B2 equilibration +
    B1 dt reconciliation (out of scope); this is construction + warm-up only.
    """
    big_R = 2.0 * p_cortex.R_cell  # span the substrate<->cortex gap
    h = build_cortex_full_simulation(
        p_cortex, p_myosin=p_myosin, p_fa=p_fa,
        fa_clutch_capture_radius=big_R, with_baoab=True,
    )
    sim = h["sim"]
    snap = sim.state.get_snapshot()
    fi = h["fa_integration"]
    assert fi is not None
    assert h["n_myosin_particles"] > 0

    # All three bond families coexist.
    btypes = list(snap.bonds.types)
    assert "integrin_ligand" in btypes
    assert "fa_actin_clutch" in btypes
    assert any(b.startswith("cortex_myosin") for b in btypes), (
        f"cortex_myosin bond types missing from {btypes}"
    )

    # Particle types coexist.
    ptypes = list(snap.particles.types)
    assert "integrin" in ptypes and "ligand" in ptypes
    assert "cortex_myosin_backbone" in ptypes
    assert "cortex_myosin_head" in ptypes

    # Tag-range bookkeeping: cortex actin at the FRONT, FA appended LAST.
    n_cortex = int(h["n_cortex_actin"])
    n_int = int(h["n_fa_integrins"])
    assert fi.integrin_tag_start >= n_cortex + h["n_myosin_particles"], (
        "integrins must be appended AFTER cortex + myosin (S5)"
    )

    # S2: clutch bonds connect integrins to CORTEX-ACTIN beads, NOT myosin.
    # cortex actin keeps its original tags [0, n_cortex_actin); myosin sits in
    # [n_cortex, n_cortex + n_myosin); integrins are appended after both.
    cp = np.asarray(fi.clutch_pairs)
    assert cp.shape[0] > 0, "expected clutch bonds with the wide capture radius"
    assert np.all(cp[:, 0] >= fi.integrin_tag_start), "col 0 must be an integrin"
    assert np.all(cp[:, 0] < fi.integrin_tag_start + n_int)
    assert np.all(cp[:, 1] < n_cortex), (
        "clutch col 1 must be a cortex-actin tag (< n_cortex_actin), NOT a "
        "myosin tag"
    )

    # A clutch bond joins an integrin to a REAL cortex-actin bead, by tag +
    # particle type, AFTER unification (mechanistic-integrity check).
    typeid = np.asarray(snap.particles.typeid)
    int_tid = ptypes.index("integrin")
    actin_tid = ptypes.index("actin_cortex")
    a0, b0 = int(cp[0, 0]), int(cp[0, 1])
    assert typeid[a0] == int_tid, "clutch col 0 particle must be type integrin"
    assert typeid[b0] == actin_tid, (
        "clutch col 1 particle must be type actin_cortex (a real cortex bead)"
    )

    # Short BAOAB warm-up: no NaN/Inf, substrate ligands immobile at z=0.
    lig = np.arange(
        fi.ligand_tag_start, fi.ligand_tag_start + fi.n_substrate_ligands
    )
    pre = np.asarray(sim.state.get_snapshot().particles.position).copy()
    sim.run(0)
    sim.run(200)
    post = np.asarray(sim.state.get_snapshot().particles.position)
    assert np.all(np.isfinite(post)), "NaN/Inf after FA+myosin warm-up"
    assert np.allclose(post[lig], pre[lig]), "substrate ligands drifted (S0 pin)"
    # FLOOR-FIX (a): ligands are seeded directly UNDER a south-cap cortex bead
    # (contact-footprint seeding), so each ligand's xy coincides with a
    # cortex-actin bead (was: all ligands on the z=0 equatorial plane). The S0
    # immobility contract is asserted just above; verify the contact geometry.
    # Check CONSTRUCTION seeding (pre): the immobile ligand sits under its cortex
    # bead's construction xy. Post-warm-up the cortex bead has drifted thermally
    # ~nm while the pinned ligand stayed put, so compare against pre, not post.
    _n_cortex = int(h["n_cortex_actin"])
    _lig_xy = pre[lig][:, :2]
    _cortex_xy = pre[:_n_cortex, :2]
    _dmin = np.min(
        np.linalg.norm(_lig_xy[:, None, :] - _cortex_xy[None, :, :], axis=2), axis=1
    )
    assert np.all(_dmin < 1e-9), "each ligand must sit under a cortex-actin bead"


def _resolve_xlinks_or_skip(cortex_cfg, p_cortex):
    """Resolve the cortex crosslinker block from the h3 config (small demo)."""
    from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
    cfg = deepcopy(cortex_cfg)
    if "dynamic_crosslinkers" not in cfg.get("cortex", {}):
        pytest.skip("no dynamic_crosslinkers block in phase1_h3.yaml")
    cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = 8  # small demo
    return resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)


def test_fa_with_xlinks_builds(cortex_cfg, p_cortex, p_fa):
    """S5 smoke: FA + crosslinkers coexist + build (was NotImplementedError).

    Asserts integrin/ligand + integrin_ligand + xlink bond types all present
    and the FA block is appended after the cortex + xlink-head block.
    """
    p_xlinks = _resolve_xlinks_or_skip(cortex_cfg, p_cortex)
    big_R = 2.0 * p_cortex.R_cell
    h = build_cortex_full_simulation(
        p_cortex, p_xlinks=p_xlinks, p_fa=p_fa,
        fa_clutch_capture_radius=big_R,
    )
    snap = h["sim"].state.get_snapshot()
    btypes = list(snap.bonds.types)
    assert "integrin_ligand" in btypes
    assert "fa_actin_clutch" in btypes
    assert any(b.startswith("xlink") for b in btypes), btypes
    assert h["n_xlink_heads"] > 0
    fi = h["fa_integration"]
    n_cortex = int(h["n_cortex_actin"])
    assert fi.integrin_tag_start >= n_cortex + h["n_xlink_heads"], (
        "integrins appended after cortex + xlink heads (S5)"
    )
    # clutch col 1 stays a cortex-actin tag (not an xlink-head tag).
    cp = np.asarray(fi.clutch_pairs)
    if cp.shape[0] > 0:
        assert np.all(cp[:, 1] < n_cortex)


def _resolve_lamel_or_skip(cortex_cfg, p_cortex):
    """Resolve a small lamellipodium demo from configs/phase1_h5.yaml.

    The lamellipodium config lives in its own ``phase1_h5.yaml`` (it is not a
    sub-block of the cortex config), loaded the same way test_lamellipodium.py
    does. Override n_WAVE to a small demo count for a fast smoke build.
    """
    from ffn_sim.cell.lamellipodium import resolve_h5_lamellipodium
    h5_path = CONFIG_DIR / "phase1_h5.yaml"
    if not h5_path.exists():
        pytest.skip("configs/phase1_h5.yaml not found")
    cfg = deepcopy(yaml.safe_load(h5_path.read_text()))
    block = cfg.get("lamellipodium", cfg)
    block["n_WAVE"] = 3  # small demo
    return resolve_h5_lamellipodium(
        cfg, L_box=p_cortex.L_box, dt=p_cortex.dt_cfl, kT=p_cortex.kT,
    )


def test_fa_with_lamellipodium_builds(cortex_cfg, p_cortex, p_fa):
    """S5 smoke: FA + lamellipodium coexist + build (was NotImplementedError).

    Asserts integrin/ligand + integrin_ligand + lamellipodium types present
    and the FA block is appended after the wave/lamellipodium block.
    """
    p_lamel = _resolve_lamel_or_skip(cortex_cfg, p_cortex)
    big_R = 2.0 * p_cortex.R_cell
    h = build_cortex_full_simulation(
        p_cortex, p_lamellipodium=p_lamel, p_fa=p_fa,
        fa_clutch_capture_radius=big_R,
    )
    snap = h["sim"].state.get_snapshot()
    ptypes = list(snap.particles.types)
    btypes = list(snap.bonds.types)
    assert "integrin" in ptypes and "ligand" in ptypes
    assert "wave_particle" in ptypes
    assert "integrin_ligand" in btypes
    assert "fa_actin_clutch" in btypes
    assert h["n_wave_particles"] > 0
    fi = h["fa_integration"]
    n_cortex = int(h["n_cortex_actin"])
    # integrins appended after cortex + wave + lamellipodium actin block.
    assert fi.integrin_tag_start >= (
        n_cortex + h["n_wave_particles"] + h["n_lamellipodium_actin"]
    ), "integrins appended after cortex + lamellipodium (S5)"
    cp = np.asarray(fi.clutch_pairs)
    if cp.shape[0] > 0:
        assert np.all(cp[:, 1] < n_cortex)


# ---------------------------------------------------------------------------
# S0/S1/S2 -- ON path builds (cortex + FA, the slice scope)
# ---------------------------------------------------------------------------
def test_fa_on_builds(p_cortex, p_fa):
    """FA on (cortex + FA): substrate ligands at z=0, both bond types
    registered, integrins APPENDED last (S5) at tags
    [n_cortex_actin, +n_int), and a fa_actin_clutch bond actually connects an
    integrin to a cortex-actin tag.

    The clutch is exercised with an explicit capture radius spanning the
    construction substrate<->cortex gap (documented in the module
    docstring).
    """
    big_R = 2.0 * p_cortex.R_cell  # spans cell diameter => every integrin
    h = build_cortex_full_simulation(
        p_cortex, p_fa=p_fa, fa_clutch_capture_radius=big_R,
    )
    sim = h["sim"]
    snap = sim.state.get_snapshot()
    fi = h["fa_integration"]
    assert fi is not None

    # particle/bond types registered
    assert "integrin" in list(snap.particles.types)
    assert "ligand" in list(snap.particles.types)
    assert "integrin_ligand" in list(snap.bonds.types)
    assert "fa_actin_clutch" in list(snap.bonds.types)

    # counts
    n_int_expected = (
        (p_fa.n_nascent_per_cell + p_fa.n_mature_per_cell) * p_fa.n_total_per_fa
    )
    n_lig_expected = p_fa.n_nascent_per_cell + p_fa.n_mature_per_cell
    assert h["n_fa_integrins"] == n_int_expected
    assert h["n_substrate_ligands"] == n_lig_expected

    # S5: integrins APPENDED last, just after the cortex-actin block (no other
    # subsystems here), at global tags [n_cortex_actin, n_cortex_actin+n_int).
    n_cortex = int(h["n_cortex_actin"])
    assert fi.integrin_tag_start == n_cortex
    pos = np.asarray(snap.particles.position)
    typeid = np.asarray(snap.particles.typeid)
    int_tid = list(snap.particles.types).index("integrin")
    assert np.all(
        typeid[n_cortex:n_cortex + n_int_expected] == int_tid
    ), "S5: integrins occupy the appended block [n_cortex_actin, +n_int)"

    # FLOOR-FIX (a): ligands are seeded directly UNDER a south-cap cortex bead
    # (contact-footprint), so each ligand's xy coincides with a cortex-actin bead
    # (was: all ligands on the z=0 plane). Substrate immobility is the S0 pin's job.
    lig = np.arange(
        fi.ligand_tag_start, fi.ligand_tag_start + fi.n_substrate_ligands
    )
    _lig_xy = pos[lig][:, :2]
    _cortex_xy = pos[:n_cortex, :2]
    _dmin = np.min(
        np.linalg.norm(_lig_xy[:, None, :] - _cortex_xy[None, :, :], axis=2), axis=1
    )
    assert np.all(_dmin < 1e-9), "each ligand must sit under a cortex-actin bead"

    # S1: integrin_ligand catch updater attached, none bound at construction
    assert h["integrin_action"] is not None
    assert h["integrin_updater"] is not None
    assert h["integrin_action"].n_engaged == 0

    # S0: ligand pin attached
    assert h["ligand_pin_action"] is not None
    assert h["ligand_pin_updater"] is not None

    # S2: clutch bonds connect integrin -> cortex-actin (the literal load path).
    # S5 append-last: cortex actin keeps tags [0, n_cortex_actin); integrins
    # are appended at [n_cortex_actin, +n_int).
    assert fi.n_clutch_bonds == n_int_expected, (
        "with a cell-diameter capture radius every integrin should clutch to "
        "its nearest cortex-actin bead"
    )
    cp = np.asarray(fi.clutch_pairs)
    # S5 append-last: integrin tags >= n_cortex, cortex-actin tags < n_cortex.
    assert np.all(cp[:, 0] >= n_cortex), "clutch col 0 must be an integrin tag"
    assert np.all(cp[:, 0] < n_cortex + n_int_expected)
    assert np.all(cp[:, 1] < n_cortex), (
        "clutch col 1 must be an (unshifted) cortex-actin tag"
    )

    # S2 construction force-free: each clutch bond's r0 == its initial
    # integrin<->actin separation (per-bond exact-r0 binning).
    bond = None
    for f in sim.operations.integrator.forces:
        if "fa_actin_clutch" in getattr(f, "params", {}):
            bond = f
            break
    assert bond is not None
    r0_canon = float(bond.params["fa_actin_clutch"]["r0"])
    assert np.isclose(r0_canon, fi.clutch_bin_r0[0])
    # k reuses the integrin bare clutch stiffness.
    assert np.isclose(float(bond.params["fa_actin_clutch"]["k"]), p_fa.k_int_bare)


def test_fa_physical_capture_radius_fully_clutches(p_cortex, p_fa):
    """FLOOR-FIX (a) [KU-3.5 v4, 2026-05-31] — *** TEST-CONTRACT CHANGE, PI-FLAGGED ***.

    This test previously asserted the BUG (``..._barely_clutches``): that the
    physical ``capture_radius_R_FA`` (1.5 um) "barely clutches" because the cortex
    shell (origin-centred, south pole at z ~ -R_cell) and the z=0 substrate are
    spatially disjoint, leaving the literal load path essentially absent (the
    KU-3.5 v4 gamma floor). The contact-footprint seeding
    (``cell.py._extend_snapshot_with_fa``) now seeds each FA directly under a
    south-cap cortex bead, so the PHYSICAL radius fully clutches with NO
    cell-diameter override. Inverted to regression-guard the fix.
    """
    h = build_cortex_full_simulation(p_cortex, p_fa=p_fa)  # physical radius
    fi = h["fa_integration"]
    n_int = h["n_fa_integrins"]
    # The clutch fan-in cap (_MAX_CLUTCH_FANIN in cell.py) trades a few clutch
    # bonds for HOOMD nlist-exclusion safety (the "Too many bonds to process
    # exclusions" overflow on long production runs), so a small fraction of
    # integrins in dense FA disks stay unclutched. The contact-footprint fix is
    # proven by clutching the MAJORITY (>=70%) vs the old disjoint-geometry ~17%
    # (9/52). (PI-flagged: the clutch count is a HOOMD-capacity tradeoff.)
    assert fi.n_clutch_bonds >= 0.7 * n_int, (
        f"physical capture_radius_R_FA={p_fa.capture_radius_R_FA:.2e} m should "
        f"clutch most of the {n_int} integrins via contact-footprint seeding; "
        f"got {fi.n_clutch_bonds} (old disjoint geometry clutched ~17%)."
    )
    assert n_int > 0
    assert h["n_substrate_ligands"] > 0
    assert h["integrin_action"] is not None


def test_cell_build_with_fa_requires_p_fa(p_cortex):
    """Cell.build(with_fa=True) without p_fa raises (mirrors with_myosin)."""
    with pytest.raises(ValueError, match="requires p_fa"):
        Cell.build(p_cortex, options=CellBuildOptions(with_fa=True))


def test_cell_build_fa_on_tag_ranges(p_cortex, p_fa):
    """Cell.build(with_fa=True), cortex + FA: populates cell.fa + integrin /
    substrate_ligand tag ranges; tag ranges stay contiguous & cover N.

    FA-on layout (S5 tag-space unification) APPENDS the FA block LAST, so
    ``fa_integrin`` starts AFTER the cortex block and ``substrate_ligand``
    follows it (cortex_actin keeps tag 0).
    """
    big_R = 2.0 * p_cortex.R_cell
    cell = Cell.build(
        p_cortex, p_fa=p_fa,
        fa_clutch_capture_radius=big_R,
        options=CellBuildOptions(with_fa=True),
    )
    assert cell.fa is not None
    tr = cell.tag_ranges()
    assert "substrate_ligand" in tr
    # Cortex actin keeps tag 0 (S5 append-last); FA block comes after it.
    assert tr["cortex_actin"][0] == 0
    lo, hi = tr["fa_integrin"]
    assert lo >= cell.n_cortex_actin, "FA-on layout appends integrins last (S5)"
    assert hi - lo == cell.fa.n_integrins > 0
    lo2, hi2 = tr["substrate_ligand"]
    assert hi2 - lo2 == cell.fa.n_substrate_ligands > 0

    # contiguity + full coverage
    prev = 0
    for name, (s, e) in sorted(tr.items(), key=lambda kv: kv[1][0]):
        assert s == prev, f"tag ranges not contiguous at {name}: {s} != {prev}"
        prev = max(prev, e)
    assert prev == cell.n_particles

    # bead-count summary still sums to N (now includes integrin + ligand).
    assert sum(cell.bead_count_summary().values()) == cell.n_particles


# ---------------------------------------------------------------------------
# Short BAOAB warm-up -- construction + a few hundred steps (NOT production)
# ---------------------------------------------------------------------------
def test_fa_short_run_no_nan(p_cortex, p_fa):
    """Short BAOAB warm-up (cortex + FA): no NaN/Inf; substrate ligands stay
    immobile.

    NOTE: SHORT warm-up at the cortex CFL dt -- full production stability
    needs the PI-gated equilibration prelude (B2) + dt reconciliation (B1);
    this test validates construction + warm-up only.
    """
    big_R = 2.0 * p_cortex.R_cell  # form real clutch bonds for the warm-up
    h = build_cortex_full_simulation(
        p_cortex, p_fa=p_fa, fa_clutch_capture_radius=big_R, with_baoab=True,
    )
    sim = h["sim"]
    fi = h["fa_integration"]
    lig = np.arange(
        fi.ligand_tag_start, fi.ligand_tag_start + fi.n_substrate_ligands
    )
    pre = np.asarray(sim.state.get_snapshot().particles.position).copy()

    sim.run(0)            # seed net_force for the L-M Action
    sim.run(300)          # short warm-up (NOT production)

    post = np.asarray(sim.state.get_snapshot().particles.position)
    assert np.all(np.isfinite(post)), "NaN/Inf after FA warm-up"

    # S0 immobility: substrate ligands held at their construction position.
    assert np.allclose(post[lig], pre[lig]), (
        "substrate ligands drifted during warm-up (S0 pin failed)"
    )
    # FLOOR-FIX (a): ligands sit UNDER south-cap cortex beads (contact-footprint),
    # not on the old z=0 plane — verify the contact geometry (xy under a bead).
    # Check CONSTRUCTION seeding (pre): the immobile ligand sits under its cortex
    # bead's construction xy. Post-warm-up the cortex bead has drifted thermally
    # ~nm while the pinned ligand stayed put, so compare against pre, not post.
    _n_cortex = int(h["n_cortex_actin"])
    _lig_xy = pre[lig][:, :2]
    _cortex_xy = pre[:_n_cortex, :2]
    _dmin = np.min(
        np.linalg.norm(_lig_xy[:, None, :] - _cortex_xy[None, :, :], axis=2), axis=1
    )
    assert np.all(_dmin < 1e-9), "each ligand must sit under a cortex-actin bead"


@pytest.mark.skipif(
    not RUN_PRODUCTION,
    reason="opt-in: set H4_FA_PRODUCTION=1 (multi-minute longer warm-up)",
)
def test_fa_longer_warmup_opt_in(p_cortex, p_fa):
    """Opt-in longer warm-up (cortex + FA) -- H4_FA_PRODUCTION=1.

    Still NOT full production (no equilibration prelude); just a longer
    construction-stability soak too slow for default CI.
    """
    big_R = 2.0 * p_cortex.R_cell
    h = build_cortex_full_simulation(
        p_cortex, p_fa=p_fa,
        fa_clutch_capture_radius=big_R, with_baoab=True,
    )
    sim = h["sim"]
    sim.run(0)
    sim.run(5000)
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position).copy()
    assert np.all(np.isfinite(pos))
