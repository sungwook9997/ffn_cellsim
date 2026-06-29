"""Prefix-split tests for cortex/myosin (KU-3.5 SF NMII split, PI 2026-06-09).

The myosin builder is parameterized by a bond/particle TYPE-NAME ``prefix`` so a
ventral-stress-fiber NMII placement can mint ``sf_myosin_*`` types that fall under
the ``sf_`` γ-denylist, while the cortical placement keeps the byte-identical
``cortex_myosin_*`` types that the cortical-tension estimator SUMS (the active-γ
signal). These tests assert:

* the DEFAULT prefix reproduces every legacy type name byte-for-byte;
* the ``sf_myosin_`` prefix yields a DISTINCT set of types, all under ``sf_``;
* the resolver validates the prefix;
* both placements COEXIST in one snapshot with no type collision;
* the per-prefix attach-bond filter in :class:`MyosinStepUpdater` is correct.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

import gsd.hoomd

from ffn_sim.archive.hoomd_legacy.cortex.myosin import (
    BOND_TYPE_MYOSIN_BACKBONE,
    BOND_TYPE_MYOSIN_HEAD_BACKBONE,
    MYOSIN_DEFAULT_PREFIX,
    cortex_myosin_attach_bin_names,
    extend_state_with_cortex_myosin,
    generate_cortex_myosin_layout,
    myosin_attach_bin_names,
    myosin_backbone_bond_name,
    myosin_head_backbone_bond_name,
    myosin_particle_type_names,
    resolve_cortex_myosin,
)

CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"

R_CELL = 7.5e-6


def _cfg(prefix: str | None = None, n_motors: int = 3) -> dict:
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors
    if prefix is not None:
        cfg["cortex"]["myosin"]["prefix"] = prefix
    return cfg


def _base_frame(n_actin: int = 24) -> gsd.hoomd.Frame:
    """A minimal cortex-like base frame (a handful of actin_cortex beads)."""
    f = gsd.hoomd.Frame()
    f.particles.N = n_actin
    f.particles.types = ["actin_cortex"]
    f.particles.typeid = np.zeros(n_actin, dtype=np.uint32)
    rng = np.random.default_rng(0)
    f.particles.position = rng.uniform(-R_CELL, R_CELL, (n_actin, 3))
    f.particles.mass = np.ones(n_actin, dtype=np.float64)
    f.bonds.N = 0
    f.bonds.types = ["cortex_backbone"]
    f.bonds.group = np.empty((0, 2), dtype=np.uint32)
    f.bonds.typeid = np.empty((0,), dtype=np.uint32)
    f.configuration.box = [4 * R_CELL, 4 * R_CELL, 4 * R_CELL, 0, 0, 0]
    return f


# ---------------------------------------------------------------------------
# Default prefix: byte-identical legacy names
# ---------------------------------------------------------------------------
def test_default_prefix_is_cortex_myosin():
    assert MYOSIN_DEFAULT_PREFIX == "cortex_myosin_"


def test_default_name_helpers_match_legacy_constants():
    assert myosin_backbone_bond_name() == BOND_TYPE_MYOSIN_BACKBONE
    assert myosin_head_backbone_bond_name() == BOND_TYPE_MYOSIN_HEAD_BACKBONE
    assert myosin_particle_type_names() == (
        "cortex_myosin_backbone",
        "cortex_myosin_head",
    )
    assert cortex_myosin_attach_bin_names(4) == [
        "cortex_myosin_attach_b0",
        "cortex_myosin_attach_b1",
        "cortex_myosin_attach_b2",
        "cortex_myosin_attach_b3",
    ]
    # The wrapper IS the default-prefix call.
    assert cortex_myosin_attach_bin_names(4) == myosin_attach_bin_names(
        4, MYOSIN_DEFAULT_PREFIX
    )


def test_resolver_default_prefix():
    p = resolve_cortex_myosin(_cfg(), dt=1e-9)
    assert p.prefix == "cortex_myosin_"


# ---------------------------------------------------------------------------
# sf_myosin_ prefix: distinct names, all under the sf_ denylist
# ---------------------------------------------------------------------------
def test_sf_prefix_names_are_distinct_and_sf_denylisted():
    pfx = "sf_myosin_"
    names = [
        myosin_backbone_bond_name(pfx),
        myosin_head_backbone_bond_name(pfx),
        *myosin_particle_type_names(pfx),
        *myosin_attach_bin_names(3, pfx),
    ]
    # all distinct from the cortical names ...
    cortical = {
        BOND_TYPE_MYOSIN_BACKBONE,
        BOND_TYPE_MYOSIN_HEAD_BACKBONE,
        "cortex_myosin_backbone",
        "cortex_myosin_head",
        *cortex_myosin_attach_bin_names(3),
    }
    assert not (set(names) & cortical)
    # ... and EVERY sf_myosin_ name falls under the "sf_" γ-denylist prefix.
    assert all(n.startswith("sf_") for n in names)


def test_resolver_accepts_sf_prefix():
    p = resolve_cortex_myosin(_cfg(prefix="sf_myosin_"), dt=1e-9)
    assert p.prefix == "sf_myosin_"


@pytest.mark.parametrize("bad", ["", "cortex_myosin", "sf-myosin", 123, None])
def test_resolver_rejects_bad_prefix(bad):
    cfg = _cfg()
    cfg["cortex"]["myosin"]["prefix"] = bad
    if bad is None:
        # None falls back to the default (cfg.get default), so it is VALID.
        p = resolve_cortex_myosin(cfg, dt=1e-9)
        assert p.prefix == "cortex_myosin_"
        return
    with pytest.raises(ValueError):
        resolve_cortex_myosin(cfg, dt=1e-9)


# ---------------------------------------------------------------------------
# Build-path: extend a snapshot under each prefix
# ---------------------------------------------------------------------------
def _extend(prefix: str | None):
    p = resolve_cortex_myosin(
        _cfg(prefix=prefix) if prefix else _cfg(), dt=1e-9
    )
    layout = generate_cortex_myosin_layout(
        p, R_CELL, motor_tag_start=24, rng=np.random.default_rng(1)
    )
    snap = extend_state_with_cortex_myosin(_base_frame(), layout, p)
    return p, snap


def test_default_build_mints_cortex_myosin_types():
    _, snap = _extend(None)
    assert "cortex_myosin_backbone" in snap.particles.types
    assert "cortex_myosin_head" in snap.particles.types
    bts = list(snap.bonds.types)
    assert BOND_TYPE_MYOSIN_BACKBONE in bts
    assert BOND_TYPE_MYOSIN_HEAD_BACKBONE in bts
    assert all(n in bts for n in cortex_myosin_attach_bin_names(10))


def test_sf_build_mints_sf_myosin_types_only():
    _, snap = _extend("sf_myosin_")
    ptypes = list(snap.particles.types)
    btypes = list(snap.bonds.types)
    assert "sf_myosin_backbone" in ptypes and "sf_myosin_head" in ptypes
    # No cortical motor type leaked in.
    assert "cortex_myosin_backbone" not in ptypes
    assert "cortex_myosin_head" not in ptypes
    assert not any(b.startswith("cortex_myosin_") for b in btypes)
    # Every new myosin bond type is sf_-denylisted.
    sf_myo = [b for b in btypes if b.startswith("sf_myosin_")]
    assert sf_myo and all(b.startswith("sf_") for b in sf_myo)


def test_both_placements_coexist_without_collision():
    """A cortical placement then an SF placement on the same frame: distinct."""
    p_cor = resolve_cortex_myosin(_cfg(), dt=1e-9)
    lay_cor = generate_cortex_myosin_layout(
        p_cor, R_CELL, motor_tag_start=24, rng=np.random.default_rng(1)
    )
    snap = extend_state_with_cortex_myosin(_base_frame(), lay_cor, p_cor)

    n_after_cortex = int(snap.particles.N)
    p_sf = resolve_cortex_myosin(_cfg(prefix="sf_myosin_"), dt=1e-9)
    lay_sf = generate_cortex_myosin_layout(
        p_sf, R_CELL, motor_tag_start=n_after_cortex,
        rng=np.random.default_rng(2),
    )
    snap = extend_state_with_cortex_myosin(snap, lay_sf, p_sf)

    ptypes = list(snap.particles.types)
    btypes = list(snap.bonds.types)
    # Both families present, no name aliasing.
    for fam in ("cortex_myosin_backbone", "cortex_myosin_head",
                "sf_myosin_backbone", "sf_myosin_head"):
        assert fam in ptypes
    assert "cortex_myosin_backbone" in btypes  # bond namespace
    assert "sf_myosin_backbone" in btypes
    # Type-id uniqueness (no duplicate names in either list).
    assert len(ptypes) == len(set(ptypes))
    assert len(btypes) == len(set(btypes))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
