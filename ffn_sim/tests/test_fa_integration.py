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
1. **Integrin-tag prepend -> FA is cortex-only in this slice.** The reused
   ``IntegrinBondUpdater`` requires integrin global tags == dense
   ``[0, n_int)`` (it indexes its per-integrin state by
   ``integrin_tag_start`` AND uses that as a snapshot-position index). So
   the FA extension PREPENDS integrins at ``[0, n_int)`` and shifts every
   other particle up by ``n_int``. That collides with the absolute-tag
   bookkeeping the myosin / xlink / lamellipodium Updaters captured at
   layout time (silently-wrong head/xlink binding), so the builder raises
   ``NotImplementedError`` for those combinations -- unifying the tag
   bookkeeping is the S5 integration step.
   :func:`test_fa_with_myosin_raises` pins that guard.

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
# S5 boundary -- FA + (myosin/xlink/lamellipodium) refused loudly
# ---------------------------------------------------------------------------
def test_fa_with_myosin_raises(p_cortex, p_myosin, p_fa):
    """FA + myosin must raise NotImplementedError in this slice.

    The integrin-tag prepend collides with the myosin Updater's absolute
    cortex-actin tag bookkeeping (silently-wrong binding); unifying that is
    S5. The builder refuses the combination rather than produce wrong
    forces.
    """
    with pytest.raises(NotImplementedError, match="not supported"):
        build_cortex_full_simulation(p_cortex, p_myosin=p_myosin, p_fa=p_fa)


# ---------------------------------------------------------------------------
# S0/S1/S2 -- ON path builds (cortex + FA, the slice scope)
# ---------------------------------------------------------------------------
def test_fa_on_builds(p_cortex, p_fa):
    """FA on (cortex + FA): substrate ligands at z=0, both bond types
    registered, integrins prepended at tags [0, n_int), and a
    fa_actin_clutch bond actually connects an integrin to a cortex-actin tag.

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

    # integrins prepended at global tags [0, n_int)
    assert fi.integrin_tag_start == 0
    pos = np.asarray(snap.particles.position)
    typeid = np.asarray(snap.particles.typeid)
    int_tid = list(snap.particles.types).index("integrin")
    assert np.all(typeid[:n_int_expected] == int_tid), (
        "the reused IntegrinBondUpdater requires integrins at tags [0, n_int)"
    )

    # substrate ligands at z=0
    lig = np.arange(
        fi.ligand_tag_start, fi.ligand_tag_start + fi.n_substrate_ligands
    )
    assert np.allclose(pos[lig, 2], 0.0), "substrate ligands must sit at z=0"

    # S1: integrin_ligand catch updater attached, none bound at construction
    assert h["integrin_action"] is not None
    assert h["integrin_updater"] is not None
    assert h["integrin_action"].n_engaged == 0

    # S0: ligand pin attached
    assert h["ligand_pin_action"] is not None
    assert h["ligand_pin_updater"] is not None

    # S2: clutch bonds connect integrin -> cortex-actin (the literal load path).
    # cortex actin lives at SHIFTED tags [n_int, n_int + n_cortex_actin).
    assert fi.n_clutch_bonds == n_int_expected, (
        "with a cell-diameter capture radius every integrin should clutch to "
        "its nearest cortex-actin bead"
    )
    cp = np.asarray(fi.clutch_pairs)
    n_cortex = int(h["n_cortex_actin"])
    assert np.all(cp[:, 0] < n_int_expected), "clutch col 0 must be an integrin tag"
    assert np.all(cp[:, 1] >= n_int_expected), (
        "clutch col 1 must be a (shifted) cortex-actin tag"
    )
    assert np.all(cp[:, 1] < n_int_expected + n_cortex)

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


def test_fa_physical_capture_radius_barely_clutches(p_cortex, p_fa):
    """Documented honest finding: the physical ``capture_radius_R_FA``
    (1.5 um) barely reaches the cortex shell, so the vast majority of
    integrins fail to form a load-path clutch at construction -- the
    literal load path is essentially absent until the cell is equilibrated
    onto the substrate (S5/B2). Asserts (a) integrins + ligands + the
    integrin_ligand catch updater still build, and (b) fewer than half the
    integrins clutch at the physical radius (vs ALL at the cell-diameter
    override in ``test_fa_on_builds``).
    """
    h = build_cortex_full_simulation(p_cortex, p_fa=p_fa)  # default radius
    fi = h["fa_integration"]
    n_int = h["n_fa_integrins"]
    assert fi.n_clutch_bonds < 0.5 * n_int, (
        f"physical capture_radius_R_FA={p_fa.capture_radius_R_FA:.2e} m "
        f"clutched {fi.n_clutch_bonds}/{n_int} integrins; the disjoint "
        "substrate/cortex geometry should leave most integrins unclutched "
        "until equilibration (S5/B2)."
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

    FA-on layout is integrin-first (the prepend contract), so
    ``fa_integrin`` starts at 0 and ``cortex_actin`` is shifted up.
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
    lo, hi = tr["fa_integrin"]
    assert lo == 0, "FA-on layout prepends integrins at tag 0"
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

    # S0 immobility: substrate ligands held at construction z=0.
    assert np.allclose(post[lig], pre[lig]), (
        "substrate ligands drifted during warm-up (S0 pin failed)"
    )
    assert np.allclose(post[lig, 2], 0.0), "substrate ligands left z=0"


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
