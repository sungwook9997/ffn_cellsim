"""A2 (H.7, 2026-06-05) — single physiological-baseline manifest loader.

ONE source of truth for "what a resting physiological MCF7 single cell IS":
``configs/mcf7_baseline.yaml`` declares every compartment with its in-vivo
setpoint + ON/OFF (no null-defaults, per the physiological-baseline HARD rule).
This module resolves that manifest into the subsystem dataclasses and hands them
to the A1-unified :meth:`ffn_sim.cell.cell.Cell.build`, so the whole cell is
assembled through a single entry point with the ``production_policy``
physiological-baseline gate enforced at build time.

The PI-ratified setpoints were externalised from the previously-hardcoded
``MCF7`` dict in ``scripts/mcf7_fullcell_stage1.py``.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from ffn_sim.cell.cell import Cell, CellBuildOptions
from ffn_sim.cell.cytoplasm import resolve_cytoplasm
from ffn_sim.cell.membrane_surface import resolve_membrane_surface
from ffn_sim.cell.nucleus import resolve_nucleus
from ffn_sim.common.production_policy import (
    require_full_cell_physiological_baseline,
)
from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.enclosed_volume import resolve_enclosed_volume
from ffn_sim.cortex.myosin import resolve_cortex_myosin

_CONFIG_DIR = Path(__file__).resolve().parents[1] / "configs"

# Optional subsystems declared in the manifest but PI-gated; enabling one via the
# loader is not yet ratified (Phase B / SCOPE / KU-5.x) -> raise, do not silently
# wire un-ratified physics.
_PI_GATED_OPTIONALS = ("fa", "substrate", "lamellipodium", "turnover", "erm", "membrane_load")


def _deep_merge(base: dict, override: dict | None) -> dict:
    """Recursively merge ``override`` into a deep copy of ``base``."""
    out = deepcopy(base)
    for key, val in (override or {}).items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], val)
        else:
            out[key] = deepcopy(val)
    return out


def _enabled(block: Any) -> bool:
    return isinstance(block, dict) and bool(block.get("enabled", False))


@dataclass(slots=True)
class ResolvedBaseline:
    """Resolved physiological-baseline subsystem configs, ready for Cell.build."""

    cell_type: str
    R_cell: float
    dtc: float
    p_cortex: Any
    p_myosin: Any
    p_xlinks: Any
    p_cytoplasm: Any
    p_enclosed_volume: Any
    p_nucleus: Any
    p_membrane_surface: Any
    # Phase-B / PI-gated optionals (always None at the baseline; reserved).
    p_fa: Any | None = None
    p_substrate: Any | None = None
    p_lamellipodium: Any | None = None
    p_turnover: Any | None = None
    p_erm: Any | None = None
    p_membrane: Any | None = None
    manifest: dict = field(default_factory=dict)

    def compartments(self) -> dict[str, Any]:
        """The compartment bundle for the production_policy baseline gate."""
        return {
            "p_cytoplasm": self.p_cytoplasm,
            "p_enclosed_volume": self.p_enclosed_volume,
            "p_nucleus": self.p_nucleus,
            "p_membrane_surface": self.p_membrane_surface,
        }


def load_manifest(path_or_name: str | Path) -> dict:
    """Load a manifest YAML by absolute path or by name under ``configs/``."""
    path = Path(path_or_name)
    if not path.is_absolute() and not path.exists():
        path = _CONFIG_DIR / path_or_name
    with open(path) as handle:
        return yaml.safe_load(handle)


def resolve_baseline(manifest: dict) -> ResolvedBaseline:
    """Resolve every enabled subsystem in ``manifest`` into its dataclass.

    Raises
    ------
    ValueError
        If a required baseline compartment (cytoplasm, enclosed-volume,
        nucleus, membrane-surface) is disabled — the physiological baseline is
        not optional.
    NotImplementedError
        If a PI-gated optional subsystem is enabled (wiring it through the
        manifest is not yet ratified; surface to PI).
    """
    cell_type = manifest.get("cell_type", "MCF7")
    R_cell = float(manifest["R_cell"])

    # --- cortex / myosin / crosslinkers: base H.3 config + MCF7 overrides ---
    base_cfg = load_manifest(manifest["base_cortex_config"])
    cfg = _deep_merge(base_cfg, manifest.get("cortex_overrides"))
    cfg.setdefault("cortex", {})["R_cell"] = R_cell
    p_cortex = resolve_h3_derived(cfg)
    tau_bend = p_cortex.gamma_b * p_cortex.rest_length**3 / p_cortex.bending_modulus
    dtc = 0.001 * tau_bend
    p_myosin = resolve_cortex_myosin(cfg, dt=dtc, R_cell=p_cortex.R_cell)
    p_xlinks = resolve_crosslinkers(cfg, dt=dtc)

    comp = manifest.get("compartments", {})

    cyto_b = comp.get("cytoplasm")
    if not _enabled(cyto_b):
        raise ValueError(
            "Physiological baseline requires cytoplasm enabled (measured "
            "viscosity, not water). Set compartments.cytoplasm.enabled: true."
        )
    p_cytoplasm = resolve_cytoplasm(cell_type=cyto_b.get("cell_type", cell_type))

    ev_b = comp.get("enclosed_volume")
    if not _enabled(ev_b):
        raise ValueError(
            "Physiological baseline requires enclosed_volume enabled (positive "
            "turgor pre-tensions the cortex). Set compartments.enclosed_volume.enabled: true."
        )
    p_enclosed_volume = resolve_enclosed_volume(
        {"turgor_dP0": float(ev_b["turgor_dP0"])}, R_cell=R_cell
    )

    nuc_b = comp.get("nucleus")
    if not _enabled(nuc_b):
        raise ValueError(
            "Physiological baseline requires nucleus enabled. "
            "Set compartments.nucleus.enabled: true."
        )
    p_nucleus = resolve_nucleus(
        {"E_nuc": float(nuc_b["E_nuc"]), "ratio_lamin": float(nuc_b["ratio_lamin"])},
        R_nuc=float(nuc_b["R_nuc_frac"]) * R_cell,
        n_beads=int(nuc_b["n_beads"]),
    )

    mem_b = comp.get("membrane_surface")
    if not _enabled(mem_b):
        raise ValueError(
            "Physiological baseline requires membrane_surface enabled. "
            "Set compartments.membrane_surface.enabled: true."
        )
    p_membrane_surface = resolve_membrane_surface(
        {"gamma_mem": float(mem_b["gamma_mem"])}, R_cell=R_cell
    )

    # --- PI-gated optional subsystems: declared-but-off; enabling is unratified ---
    opt = manifest.get("optional_subsystems", {})
    for name in _PI_GATED_OPTIONALS:
        if _enabled(opt.get(name)):
            raise NotImplementedError(
                f"optional subsystem '{name}' is PI-gated (Phase B / SCOPE / "
                "KU-5.x) and not yet wired through the manifest loader. Surface "
                "to PI before enabling — do not silently run un-ratified physics."
            )

    return ResolvedBaseline(
        cell_type=cell_type,
        R_cell=R_cell,
        dtc=dtc,
        p_cortex=p_cortex,
        p_myosin=p_myosin,
        p_xlinks=p_xlinks,
        p_cytoplasm=p_cytoplasm,
        p_enclosed_volume=p_enclosed_volume,
        p_nucleus=p_nucleus,
        p_membrane_surface=p_membrane_surface,
        manifest=manifest,
    )


def build_baseline_cell(
    path_or_name: str | Path = "mcf7_baseline.yaml",
    *,
    manifest: dict | None = None,
    device: Any | None = None,
    seed: int = 1,
    with_baoab: bool = True,
    allow_unpressurized_dev: bool = False,
) -> Cell:
    """Resolve a manifest and assemble the full cell via the unified Cell.build.

    Enforces the ``production_policy`` physiological-baseline gate (non-water
    cytoplasm + positive turgor) before any HOOMD state is created.
    """
    if manifest is None:
        manifest = load_manifest(path_or_name)
    rb = resolve_baseline(manifest)
    require_full_cell_physiological_baseline(
        rb.compartments(), allow_unpressurized_dev=allow_unpressurized_dev
    )
    opts = CellBuildOptions(
        with_crosslinkers=True, with_myosin=True, with_baoab=with_baoab
    )
    return Cell.build(
        rb.p_cortex,
        p_xlinks=rb.p_xlinks,
        p_myosin=rb.p_myosin,
        p_enclosed_volume=rb.p_enclosed_volume,
        p_membrane_surface=rb.p_membrane_surface,
        p_nucleus=rb.p_nucleus,
        p_cytoplasm=rb.p_cytoplasm,
        options=opts,
        device=device,
        rng=np.random.default_rng(seed),
    )
