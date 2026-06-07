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
from ffn_sim.bridge.fa import resolve_h4
from ffn_sim.cell.cytoplasm import resolve_cytoplasm
from ffn_sim.cell.lamellipodium import resolve_h5_lamellipodium
from ffn_sim.cell.membrane import resolve_membrane
from ffn_sim.cell.membrane_surface import resolve_membrane_surface
from ffn_sim.cell.nucleus import resolve_nucleus
from ffn_sim.common.production_policy import (
    require_full_cell_physiological_baseline,
)
from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.enclosed_volume import resolve_enclosed_volume
from ffn_sim.cortex.erm import resolve_erm
from ffn_sim.cortex.myosin import resolve_cortex_myosin
from ffn_sim.cortex.turnover import resolve_turnover
from ffn_sim.ecm.substrate import resolve_substrate

_CONFIG_DIR = Path(__file__).resolve().parents[1] / "configs"


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

    # --- Optional subsystems (resolve those the manifest enables) ---
    # KU-5.x membrane-load / lamellipodium PI-APPROVED 2026-06-06. FA-substrate
    # adhesion is the KU-3.5 settled operating point (does NOT require the
    # lamellipodium — cortical tension is a cortex property; lamellipodium is a
    # parallel active-spreading enhancement, PI 2026-06-06). Each optional
    # resolves from its referenced base_config (+ inline overrides); disabled
    # subsystems stay None and the build is unchanged.
    opt = manifest.get("optional_subsystems", {})

    def _opt_cfg(block: dict) -> dict:
        base = load_manifest(block["base_config"]) if block.get("base_config") else {}
        inline = {k: v for k, v in block.items() if k not in ("enabled", "base_config")}
        return _deep_merge(base, inline) if inline else base

    p_fa = p_substrate = p_lamellipodium = p_turnover = p_erm = p_membrane = None
    lam_cfg: dict = {}

    fa_b = opt.get("fa")
    if _enabled(fa_b):
        # FA auto-seeds the south-cap contact footprint (_extend_snapshot_with_fa)
        # so clutch bonds form; the SubstrateLigandPin is the rigid-dish adhesion.
        p_fa = resolve_h4(_opt_cfg(fa_b))

    sub_b = opt.get("substrate")
    if _enabled(sub_b):
        # k_sub REQUIRED (no magic default; ECM-condition stiffness). FA already
        # pins ligands rigidly — the compliant spring is an explicit ECM variant.
        p_substrate = resolve_substrate(_opt_cfg(sub_b), kT=p_cortex.kT)

    lam_b = opt.get("lamellipodium")
    if _enabled(lam_b):
        lam_cfg = _opt_cfg(lam_b)
        # The lamellipodium inherits the CELL's box, so Y_max / wave_area must
        # derive from this cell's L_box (resolve_h5 defaults Y_max=0.45*L_box,
        # wave_area=L_box^2). phase1_h5.yaml's literal values were authored for a
        # 30 um reconstitution box and would violate Y_max < L_box/2 on the
        # 7.5 um MCF7 cell; strip them so they auto-derive. (The H.7 single-cell
        # geometries place WAVEs from R_cell + the FA cap, not Y_max, anyway.)
        lam_sub = lam_cfg.get("lamellipodium", lam_cfg)
        if isinstance(lam_sub, dict):
            lam_sub.pop("Y_max", None)
            lam_sub.pop("wave_area", None)
        p_lamellipodium = resolve_h5_lamellipodium(
            lam_cfg, L_box=p_cortex.L_box, dt=dtc, kT=p_cortex.kT
        )

    mem_b = opt.get("membrane_load")
    if _enabled(mem_b):
        if p_lamellipodium is None:
            raise ValueError(
                "membrane_load requires lamellipodium enabled (KU-5.x load acts "
                "on the lamellipodial barbed ends; y_plane=Y_max, A_mem=wave_area)."
            )
        mem_cfg = _opt_cfg(mem_b) or lam_cfg
        p_membrane = resolve_membrane(
            mem_cfg, y_plane=p_lamellipodium.Y_max,
            A_mem=p_lamellipodium.wave_area, kT=p_cortex.kT,
        )

    erm_b = opt.get("erm")
    if _enabled(erm_b):
        p_erm = resolve_erm(_opt_cfg(erm_b), kT=p_cortex.kT, R_cell=p_cortex.R_cell)

    tov_b = opt.get("turnover")
    if _enabled(tov_b):
        tov_cfg = _opt_cfg(tov_b)
        # _opt_cfg strips the control `enabled`, so a turnover resolved from a
        # base_config (e.g. phase1_h3.yaml, enabled:false) would inherit
        # enabled=false and no-op. The manifest gate already said yes -> force
        # enabled True on the turnover sub-dict (mirror resolve_turnover's unwrap).
        _t = tov_cfg
        if isinstance(_t.get("cortex"), dict):
            _t = _t["cortex"]
        if isinstance(_t.get("turnover"), dict):
            _t = _t["turnover"]
        _t["enabled"] = True
        p_turnover = resolve_turnover(
            tov_cfg, dt=dtc, rest_length=p_cortex.rest_length
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
        p_fa=p_fa,
        p_substrate=p_substrate,
        p_lamellipodium=p_lamellipodium,
        p_turnover=p_turnover,
        p_erm=p_erm,
        p_membrane=p_membrane,
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
    constrained: bool = False,
    constrained_dt_safety: float = 1.0,
    equilibrate: bool = False,
    equilibrate_steps: int = 0,
    equilibrate_softstart_steps: int = 100,
    connected_mesh: bool = False,
    faithful_connected_mesh: bool = False,
    cm_z_struct: float = 3.7,
    cm_bundle_mult: int = 2,
) -> Cell:
    """Resolve a manifest and assemble the full cell via the unified Cell.build.

    Enforces the ``production_policy`` physiological-baseline gate (non-water
    cytoplasm + positive turgor) before any HOOMD state is created.

    For the GATE-B operating-point build (FA-adhered cell settled onto the
    substrate), pass ``equilibrate=True`` with a softstart/baoab budget so the
    forces ramp gently (a raw FA-adhered free run trips the BAOAB guard), and
    ``constrained=True`` to run the rigid M-SHAKE backbone (required for the
    rigid Lagrange gamma channel).
    """
    if manifest is None:
        manifest = load_manifest(path_or_name)
    rb = resolve_baseline(manifest)
    require_full_cell_physiological_baseline(
        rb.compartments(), allow_unpressurized_dev=allow_unpressurized_dev
    )
    # Lamellipodium leading-edge geometry (H.7 single-cell): the manifest's
    # lamellipodium block may declare geometry (flat_plane | basal_ring |
    # polarized_patch) + an optional polarization for the patch.
    lam_block = manifest.get("optional_subsystems", {}).get("lamellipodium", {})
    lam_geometry = lam_block.get("geometry", "flat_plane") if isinstance(lam_block, dict) else "flat_plane"
    lam_polarization = lam_block.get("polarization") if isinstance(lam_block, dict) else None

    opts = CellBuildOptions(
        with_crosslinkers=True,
        with_myosin=True,
        with_baoab=with_baoab,
        with_lamellipodium=rb.p_lamellipodium is not None,
        with_fa=rb.p_fa is not None,
        with_erm=rb.p_erm is not None,
    )
    return Cell.build(
        rb.p_cortex,
        p_xlinks=rb.p_xlinks,
        p_myosin=rb.p_myosin,
        p_enclosed_volume=rb.p_enclosed_volume,
        p_membrane_surface=rb.p_membrane_surface,
        p_nucleus=rb.p_nucleus,
        p_cytoplasm=rb.p_cytoplasm,
        p_fa=rb.p_fa,
        p_substrate=rb.p_substrate,
        p_lamellipodium=rb.p_lamellipodium,
        lamellipodium_geometry=lam_geometry,
        lamellipodium_polarization=lam_polarization,
        p_turnover=rb.p_turnover,
        p_erm=rb.p_erm,
        p_membrane=rb.p_membrane,
        constrained=constrained,
        # constrained_dt_safety < 1 shrinks the constrained step: the rigid
        # backbone removes the BACKBONE CFL but the FA integrin catch-bond /
        # LJ / turgor forces keep their own; at the bare-cortex dtc those forces
        # spike (Pereverzev |F|/F_s > 700 overflow guard) and most seeds blow up.
        constrained_dt=(rb.dtc * float(constrained_dt_safety)) if constrained else None,
        # connected_mesh seeds a PERCOLATED bridge-different-filament crosslinker
        # network (z~3.3, giant~99%) instead of the random-anchor mesh (z~1.3,
        # giant~7%, fragmented) — required for cortical tension to transmit.
        connected_mesh=connected_mesh,
        # faithful_connected_mesh=True uses the variable-length bimodal cortex (Arp2/3
        # branches + Chugh-regulated lengths) instead of the uniform fast-hybrid.
        faithful_connected_mesh=faithful_connected_mesh,
        cm_z_struct=cm_z_struct,
        cm_bundle_mult=cm_bundle_mult,
        equilibrate=equilibrate,
        equilibrate_steps=equilibrate_steps,
        equilibrate_softstart_steps=equilibrate_softstart_steps,
        options=opts,
        device=device,
        rng=np.random.default_rng(seed),
    )
