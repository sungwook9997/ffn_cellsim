"""Compartment registry + recipe spine contract tests (H.7 platform, 2026-06-08).

Pin the invariants of ffn_sim/cell/compartment_registry.py:
  * every compartment is default-OFF except the actomyosin core + physiological
    baseline (no silent baseline drop, no surprise-ON compartment);
  * geometry-only compartments add zero force-bearing bonds/particles
    (mesh-as-physics rejected);
  * every gamma-contaminating compartment declares its denylisted bonds, and the
    LIVE ones are actually covered by the cortical-tension estimator denylist;
  * the six recipes load (extends chains resolve) and compose onto the MCF7
    baseline manifest WITHOUT wiring un-ratified physics;
  * composition refuses to drop a baseline compartment and refuses to force an
    EXPERIMENTAL/STUB compartment into the manifest.

Import-light: the registry has no HOOMD dependency. The one test that round-trips
a composed manifest through resolve_baseline imports HOOMD-backed resolvers and is
skipped if HOOMD is unavailable.
"""

from __future__ import annotations

import pytest

from ffn_sim.cell.compartment_registry import (
    REGISTRY,
    BaselineDropError,
    CompartmentRegistry,
    CompartmentStatus,
    GpuPath,
    HotPathPriority,
    UnratifiedCompartmentError,
    list_recipes,
    load_recipe,
)
from ffn_sim.cell.manifest import load_manifest

RECIPE_NAMES = (
    "suspended_round",
    "adherent_passive",
    "adherent_active_spread",
    "multicell_junction",
    "confined_migration",
    "full_physiological",
)


# --------------------------------------------------------------------------
# Catalogue invariants
# --------------------------------------------------------------------------
def test_registry_lists_all_expected_compartments():
    names = set(REGISTRY.names())
    expected = {
        # core
        "cortex", "crosslinkers", "myosin",
        # baseline
        "cytoplasm", "enclosed_volume", "nucleus", "membrane_surface",
        # live adhesion/spreading
        "fa", "rigid_ligand_coating", "substrate", "lamellipodium",
        "membrane_load", "erm", "turnover",
        # geometry
        "surface_manifold",
        # missing (experimental/stub)
        "ventral_stress_fibers", "linc", "intermediate_filaments",
        "microtubules", "osmotic_regulation", "membrane_reservoir",
        "cadherin_junction", "junctional_actin",
    }
    assert expected <= names, f"missing compartments: {expected - names}"


def test_defaults_are_off_except_core_and_baseline():
    REGISTRY.validate_defaults()  # raises on violation
    on = set(REGISTRY.default_on())
    assert on == {
        "cortex", "crosslinkers", "myosin",          # core
        "cytoplasm", "enclosed_volume", "nucleus", "membrane_surface",  # baseline
    }


def test_baseline_required_set():
    assert set(REGISTRY.baseline_required()) == {
        "cytoplasm", "enclosed_volume", "nucleus", "membrane_surface",
    }


def test_geometry_only_adds_no_physics():
    REGISTRY.validate_manifold_geometry_only()  # raises on violation
    sm = REGISTRY.get("surface_manifold")
    assert sm.geometry_only is True
    assert sm.performance_contract.per_step_force is False
    assert sm.performance_contract.n_bonds in (0, "0")
    assert sm.performance_contract.particle_types_added == ()


def test_gamma_contract_declared():
    REGISTRY.validate_gamma_contract()  # gamma_contaminating => denylist non-empty


def test_live_contaminating_compartments_are_actually_denylisted():
    """LIVE/CORE compartments that contaminate gamma must already be excluded by
    the cortical-tension estimator's denylist (no silent contamination today)."""
    ct = pytest.importorskip("ffn_sim.cortex.cortical_tension")

    def covered(bt: str) -> bool:
        # Use the estimator's ACTUAL non-cortical filter — now registry-driven
        # (ADHESION_BOND_TYPES + ADHESION_BOND_TYPE_PREFIXES +
        # NONCORTICAL_COMPARTMENT_PREFIXES), so a LIVE compartment's declared
        # prefix (e.g. 'mt_') is honored without a hardcoded estimator edit.
        return ct._is_adhesion_bond_type(bt) or ct._is_adhesion_bond_type(bt + "x")

    for spec in REGISTRY.all():
        if not spec.gamma_contaminating:
            continue
        if spec.status in (CompartmentStatus.CORE, CompartmentStatus.LIVE):
            for bt in spec.denylist_bond_types:
                assert covered(bt), (
                    f"LIVE compartment '{spec.name}' bond '{bt}' is not in the "
                    f"cortical-tension denylist — it would contaminate gamma."
                )


def test_experimental_contaminating_compartments_declare_future_denylist():
    """Experimental/stub compartments that will contaminate gamma must declare the
    bond prefixes that MUST be added to the estimator denylist when they go LIVE."""
    for spec in REGISTRY.all():
        if spec.gamma_contaminating and spec.status in (
            CompartmentStatus.EXPERIMENTAL, CompartmentStatus.STUB,
        ):
            assert spec.denylist_bond_types, spec.name


def test_p0_compartments_have_gpu_path():
    """P0 = always-on per-step work => must have a GPU-resident/native path
    planned (HARD rule: GPU before merge for P0)."""
    for spec in REGISTRY.all():
        pc = spec.performance_contract
        if pc.hot_path_priority is HotPathPriority.P0:
            ok = (
                spec.gpu_readiness.device_resident_now
                or pc.gpu_path_now in (GpuPath.CUPY, GpuPath.NATIVE_PLUGIN, GpuPath.HOOMD_BUILTIN)
                or pc.native_forcecompute_candidate
            )
            assert ok, f"P0 compartment '{spec.name}' has no planned GPU path."


def test_requires_reference_real_compartments():
    names = set(REGISTRY.names())
    for spec in REGISTRY.all():
        for dep in spec.requires:
            assert dep in names, f"{spec.name} requires unknown '{dep}'"


# --------------------------------------------------------------------------
# Recipes
# --------------------------------------------------------------------------
def test_all_recipes_present():
    assert set(RECIPE_NAMES) <= set(list_recipes())


@pytest.mark.parametrize("recipe_name", RECIPE_NAMES)
def test_recipe_loads_and_only_names_known_compartments(recipe_name):
    recipe = load_recipe(recipe_name)
    names = set(REGISTRY.names())
    for n in recipe["enable"]:
        assert n in names, f"recipe {recipe_name} enables unknown '{n}'"
    for n in recipe.get("declare_pending", []) or []:
        assert n in names, f"recipe {recipe_name} declares unknown '{n}'"


def test_recipe_extends_chain_resolves():
    # adherent_active_spread -> adherent_passive -> suspended_round
    r = load_recipe("adherent_active_spread")
    enable = set(r["enable"])
    # inherited baseline:
    assert {"cortex", "myosin", "cytoplasm", "membrane_surface"} <= enable
    # inherited FA layer:
    assert {"fa", "rigid_ligand_coating"} <= enable
    # own active layer:
    assert {"lamellipodium", "membrane_load"} <= enable


def test_suspended_round_composes_with_adhesion_off():
    base = load_manifest("mcf7_baseline.yaml")
    recipe = load_recipe("suspended_round")
    manifest, deferred = REGISTRY.compose_manifest(recipe, base_manifest=base)
    assert deferred == []
    opt = manifest["optional_subsystems"]
    assert opt["fa"]["enabled"] is False
    assert opt["lamellipodium"]["enabled"] is False
    # baseline stays ON
    assert manifest["compartments"]["cytoplasm"]["enabled"] is True
    assert manifest["compartments"]["enclosed_volume"]["enabled"] is True


def test_adherent_passive_turns_fa_on():
    base = load_manifest("mcf7_baseline.yaml")
    recipe = load_recipe("adherent_passive")
    manifest, deferred = REGISTRY.compose_manifest(recipe, base_manifest=base)
    assert deferred == []
    assert manifest["optional_subsystems"]["fa"]["enabled"] is True
    # substrate (compliant) stays OFF — rigid ligand coating is the default
    assert manifest["optional_subsystems"]["substrate"]["enabled"] is False
    assert manifest["optional_subsystems"]["lamellipodium"]["enabled"] is False


def test_active_spread_turns_lamellipodium_on():
    base = load_manifest("mcf7_baseline.yaml")
    recipe = load_recipe("adherent_active_spread")
    manifest, deferred = REGISTRY.compose_manifest(recipe, base_manifest=base)
    opt = manifest["optional_subsystems"]
    assert opt["fa"]["enabled"] is True
    assert opt["lamellipodium"]["enabled"] is True
    assert opt["membrane_load"]["enabled"] is True
    # ventral_stress_fibers is declared pending, not enabled -> not in deferred
    # (declare_pending is documentation; only `enable` drives composition).
    assert "ventral_stress_fibers" not in deferred


def test_compose_does_not_mutate_base_manifest():
    base = load_manifest("mcf7_baseline.yaml")
    before = base["optional_subsystems"]["fa"]["enabled"]
    recipe = load_recipe("adherent_passive")
    REGISTRY.compose_manifest(recipe, base_manifest=base)
    assert base["optional_subsystems"]["fa"]["enabled"] == before  # unchanged


def test_multicell_junction_defers_cadherin():
    base = load_manifest("mcf7_baseline.yaml")
    recipe = load_recipe("multicell_junction")
    manifest, deferred = REGISTRY.compose_manifest(recipe, base_manifest=base)
    # cadherin/junc are declare_pending (not in enable) -> composes to baseline.
    assert manifest["optional_subsystems"]["fa"]["enabled"] is False
    assert deferred == []  # nothing forced into enable


def test_full_physiological_enables_only_baseline():
    base = load_manifest("mcf7_baseline.yaml")
    recipe = load_recipe("full_physiological")
    manifest, deferred = REGISTRY.compose_manifest(recipe, base_manifest=base)
    opt = manifest["optional_subsystems"]
    for name in ("fa", "substrate", "lamellipodium", "membrane_load", "erm", "turnover"):
        assert opt[name]["enabled"] is False
    assert deferred == []


def test_forcing_experimental_into_enable_raises():
    base = load_manifest("mcf7_baseline.yaml")
    recipe = {"name": "bad", "enable": list(load_recipe("suspended_round")["enable"]) + ["cadherin_junction"]}
    # cadherin_junction requires nothing extra but is EXPERIMENTAL -> strict raises.
    with pytest.raises(UnratifiedCompartmentError):
        REGISTRY.compose_manifest(recipe, base_manifest=base, strict=True)


def test_forcing_experimental_non_strict_defers():
    base = load_manifest("mcf7_baseline.yaml")
    # membrane_reservoir is still EXPERIMENTAL (osmotic_regulation / microtubules /
    # intermediate_filaments / linc graduated LIVE 2026-06-09). It requires
    # membrane_surface, which is baseline-required (suspended_round enables it),
    # so it defers cleanly (no dependency error).
    recipe = {
        "name": "bad",
        "enable": list(load_recipe("suspended_round")["enable"]) + ["membrane_reservoir"],
    }
    manifest, deferred = REGISTRY.compose_manifest(recipe, base_manifest=base, strict=False)
    assert "membrane_reservoir" in deferred


def test_dropping_baseline_compartment_raises():
    base = load_manifest("mcf7_baseline.yaml")
    enable = [n for n in load_recipe("suspended_round")["enable"] if n != "cytoplasm"]
    recipe = {"name": "bad", "enable": enable}
    with pytest.raises(BaselineDropError):
        REGISTRY.compose_manifest(recipe, base_manifest=base)


def test_missing_dependency_raises():
    base = load_manifest("mcf7_baseline.yaml")
    # membrane_load requires lamellipodium; enable it without lamellipodium.
    enable = list(load_recipe("suspended_round")["enable"]) + ["membrane_load"]
    recipe = {"name": "bad", "enable": enable}
    with pytest.raises(ValueError):
        REGISTRY.compose_manifest(recipe, base_manifest=base)


def test_registry_is_reconstructable():
    fresh = CompartmentRegistry()
    assert set(fresh.names()) == set(REGISTRY.names())


# --------------------------------------------------------------------------
# Round-trip: a composed manifest still resolves through the unchanged loader
# --------------------------------------------------------------------------
def test_composed_suspended_manifest_resolves():
    pytest.importorskip("hoomd")
    from ffn_sim.cell.manifest import resolve_baseline

    base = load_manifest("mcf7_baseline.yaml")
    recipe = load_recipe("suspended_round")
    manifest, _ = REGISTRY.compose_manifest(recipe, base_manifest=base)
    rb = resolve_baseline(manifest)
    assert rb.p_fa is None              # adhesion OFF
    assert rb.p_lamellipodium is None
    assert rb.p_cytoplasm is not None   # baseline ON
    assert rb.p_enclosed_volume is not None


def test_composed_adherent_manifest_resolves_with_fa():
    pytest.importorskip("hoomd")
    from ffn_sim.cell.manifest import resolve_baseline

    base = load_manifest("mcf7_baseline.yaml")
    recipe = load_recipe("adherent_passive")
    manifest, _ = REGISTRY.compose_manifest(recipe, base_manifest=base)
    rb = resolve_baseline(manifest)
    assert rb.p_fa is not None          # FA resolved
    assert rb.p_substrate is None       # compliant substrate OFF (rigid pin default)


# ---------------------------------------------------------------------------
# Registry <-> cortical_tension γ-denylist cross-check (registry-driven, 2026-06-09)
# ---------------------------------------------------------------------------
def test_cortical_tension_honours_registry_gamma_denylist():
    """cortical_tension excludes every registry γ-denylist prefix EXCEPT cortex_*.

    The estimator is now registry-driven: activating a new γ-contaminating
    compartment auto-excludes its load path. The single exception is any
    ``cortex_`` prefix — ``cortex_myosin_*`` IS the active-γ signal and must stay
    COUNTED (the SF-motor reuse is handled by the distinct ``sf_myosin_*`` prefix).
    """
    from ffn_sim.cortex.cortical_tension import _is_adhesion_bond_type

    deny = REGISTRY.gamma_denylist()
    assert deny  # non-empty
    for pfx in deny:
        if pfx.startswith("cortex_"):
            # cortex_myosin_ MUST remain counted (the active-γ signal).
            assert not _is_adhesion_bond_type(pfx + "backbone"), pfx
        else:
            assert _is_adhesion_bond_type(pfx + "x"), pfx
    # The cortex's own bonds are never excluded.
    assert not _is_adhesion_bond_type("cortex_bond")
    assert not _is_adhesion_bond_type("cortex_myosin_head_backbone")


def test_internal_live_recipe_composes_and_builds_without_contamination():
    """The 3 LIVE internal compartments (osmotic + MT + IF) compose in ONE cell
    with no mutual cortical-γ contamination (the activation integration capstone)."""
    from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
    ct = pytest.importorskip("ffn_sim.cortex.cortical_tension")

    base = load_manifest("mcf7_baseline.yaml")
    manifest, deferred = REGISTRY.compose_manifest(
        load_recipe("internal_live"), base_manifest=base, strict=True
    )
    assert deferred == []
    off = build_baseline_cell("mcf7_baseline.yaml", seed=1)
    on = build_baseline_cell("mcf7_baseline.yaml", manifest=manifest, seed=1)
    sn = on.simulation.state.get_snapshot()
    assert {"mt_bead", "mtoc", "if_bead"} <= set(sn.particles.types)
    # osmotic updater added (+1).
    assert len(on.simulation.operations.updaters) == len(off.simulation.operations.updaters) + 1
    off.simulation.run(0); on.simulation.run(0)
    go = ct.measure_cortical_tension(off.simulation, R_cell=off.p_cortex.R_cell,
                                     p_enclosed_volume=off.p_enclosed_volume)
    gn = ct.measure_cortical_tension(on.simulation, R_cell=on.p_cortex.R_cell,
                                     p_enclosed_volume=on.p_enclosed_volume)
    k = "gamma_soft_N_per_m" if "gamma_soft_N_per_m" in go else "gamma_soft"
    assert abs(gn[k] - go[k]) <= 1e-12 * max(1.0, abs(go[k]))
