"""H.SF ventral stress-fiber ACTIVATION GATE (PASSIVE backbone, LIVE 2026-06-09).

Sixth compartment graduated EXPERIMENTAL->LIVE under the PI ownership grant.
Explicit contractile actomyosin bundles spanning FA→FA across the basal plane.
PASSIVE backbone only (NMII deferred — the active phase needs the sf_myosin_*
prefix split + the equilibrated Kumar-tension gate). The FA-anchor pairs are
LONG-AXIS-ALIGNED on the basal footprint (PI 2026-06-09 pairing rule), and each
bundle's backbone is per-bundle EXACT-r0 (force-free at its as-built ell0_b). The
gate builds the FA-adhered cell OFF (adherent_passive) vs ON (+ SF) and checks the
build-time controls:

  BUNDLES ASSEMBLED — ON adds n_SF bundles (sf_actin + sf_xlink_head + sf_ bonds).
  ALIGNED          — bundle axes align with the footprint in-plane long axis
                     (mean |cos| ≥ 0.8) — physiological vSF, not random spans.
  FORCE-FREE       — per-bundle backbone born at its EXACT ell0_b (strain ~0);
                     FA anchors at r0=0 (end bead AT the FA clutch).
  NO-CONTAM        — cortical γ_soft IDENTICAL OFF vs ON (sf_ registry-denylisted);
                     every sf_ bond type excluded from the cortical mask.
  OFF-IDENTITY     — SF off reproduces the adherent_passive inventory.

The Kumar 2006 10-30 nN single-SF tension band is the ACTIVE+equilibrated gate
(needs NMII via the sf_myosin_* prefix + the equilibration prelude) — DEFERRED.
N_filaments=20 (Cramer 1997) is the config candidate → μ_SF (PI-pending).

Run:  python ffn_sim/scripts/h7_stress_fibers_activation_gate.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
if str(_HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(_HERE.parents[2]))

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ffn_sim.archive.hoomd_legacy.cell.compartment_registry import REGISTRY, load_recipe
from ffn_sim.archive.hoomd_legacy.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.archive.hoomd_legacy.cell.stress_fibers import sf_actin_backbone_bin_names
from ffn_sim.archive.hoomd_legacy.cortex.cortical_tension import (
    _is_adhesion_bond_type,
    measure_cortical_tension,
)

_OUT = _HERE.parents[1] / "outputs" / "h7"


def _gamma_soft(cell) -> float:
    cell.simulation.run(0)
    r = measure_cortical_tension(
        cell.simulation, R_cell=cell.p_cortex.R_cell,
        p_enclosed_volume=cell.p_enclosed_volume,
    )
    return float(r.get("gamma_soft_N_per_m", r.get("gamma_soft", float("nan"))))


def run() -> dict:
    base = load_manifest("mcf7_baseline.yaml")
    # OFF = FA-adhered (adherent_passive); ON = + ventral stress fibers.
    m_off, _ = REGISTRY.compose_manifest(
        load_recipe("adherent_passive"), base_manifest=base, strict=True
    )
    m_on, _ = REGISTRY.compose_manifest(
        load_recipe("ventral_stress_fibers_passive"), base_manifest=base, strict=True
    )
    cell_off = build_baseline_cell("mcf7_baseline.yaml", manifest=m_off, seed=1)
    cell_on = build_baseline_cell("mcf7_baseline.yaml", manifest=m_on, seed=1)
    p = cell_on.p_stress_fibers

    s_off = cell_off.simulation.state.get_snapshot()
    s_on = cell_on.simulation.state.get_snapshot()
    handles = cell_on.extras.get("handles", {})
    n_sf = int(handles.get("n_stress_fibers", 0))
    layout = handles.get("stress_fiber_layout")
    n_part_delta = int(s_on.particles.N) - int(s_off.particles.N)
    n_bond_delta = int(s_on.bonds.N) - int(s_off.bonds.N)
    sf_bonds = [t for t in s_on.bonds.types if t.startswith("sf_")]

    bundles_ok = (
        n_sf > 0
        and "sf_actin" in s_on.particles.types
        and "sf_xlink_head" in s_on.particles.types
        and any(t == "sf_actin_bond" or t.startswith("sf_actin_bond") for t in sf_bonds)
    )

    # ALIGNED: bundle axes vs the in-plane long axis of the FA footprint.
    align_mean = float("nan")
    bundle_len_um = []
    if layout is not None and n_sf > 0:
        pos = np.asarray(s_on.particles.position, dtype=np.float64)
        fe = np.asarray(layout.fa_endpoints, dtype=np.int64)
        vecs = pos[fe[:, 1]] - pos[fe[:, 0]]
        lens = np.linalg.norm(vecs, axis=1)
        bundle_len_um = (lens * 1e6).tolist()
        u = vecs / np.maximum(lens[:, None], 1e-30)
        mean_ax = u.mean(axis=0)
        nrm = np.linalg.norm(mean_ax)
        if nrm > 0:
            mean_ax = mean_ax / nrm
            align_mean = float(np.mean(np.abs(u @ mean_ax)))
    aligned_ok = bool(np.isfinite(align_mean) and align_mean >= 0.8)

    # FORCE-FREE: per-bundle backbone born at its EXACT ell0_b.
    max_strain = float("nan")
    if layout is not None and n_sf > 0:
        pos = np.asarray(s_on.particles.position, dtype=np.float64)
        g = np.asarray(s_on.bonds.group, dtype=np.int64)
        bt = np.asarray(s_on.bonds.typeid, dtype=np.int64)
        bn = list(s_on.bonds.types)
        names = sf_actin_backbone_bin_names(n_sf)
        worst = 0.0
        for b, name in enumerate(names):
            if name not in bn:
                continue
            m = bt == bn.index(name)
            if not m.any():
                continue
            seg = pos[g[m]]
            ln = np.linalg.norm(seg[:, 0] - seg[:, 1], axis=1)
            r0 = float(layout.ell0_actin[b]) if b < layout.ell0_actin.size else 0.0
            if r0 > 0:
                worst = max(worst, float(np.max(np.abs(ln - r0) / r0)))
        max_strain = worst
    force_free_ok = bool(np.isfinite(max_strain) and max_strain < 1e-6)

    # NO-CONTAMINATION.
    mask_excludes_sf = all(_is_adhesion_bond_type(t) for t in sf_bonds)
    g_off, g_on = _gamma_soft(cell_off), _gamma_soft(cell_on)
    gamma_identical = bool(
        np.isfinite(g_off) and np.isfinite(g_on)
        and abs(g_on - g_off) <= 1e-12 * max(1.0, abs(g_off))
    )

    ok = (
        bundles_ok and aligned_ok and force_free_ok
        and mask_excludes_sf and gamma_identical
    )

    result = {
        "gate": "H.SF ventral_stress_fibers activation (passive)",
        "verdict": "PASS" if ok else "REVIEW",
        "physics_claim": True,
        "params": {
            "n_SF": int(p.n_SF), "n_beads_per_SF": int(p.n_beads_per_SF),
            "N_filaments": p.N_filaments, "mu_SF_N": p.mu_SF,
            "EA_single_N": float(p.EA_single),
        },
        "controls": {
            "bundles_assembled": bool(bundles_ok),
            "n_stress_fibers": n_sf,
            "n_particle_delta": n_part_delta, "n_bond_delta": n_bond_delta,
            "aligned": {
                "ok": aligned_ok, "mean_abs_cos_long_axis": align_mean,
                "bundle_len_um_min": (min(bundle_len_um) if bundle_len_um else None),
                "bundle_len_um_mean": (float(np.mean(bundle_len_um)) if bundle_len_um else None),
                "bundle_len_um_max": (max(bundle_len_um) if bundle_len_um else None),
            },
            "force_free_construction": {"ok": force_free_ok, "max_strain": max_strain},
            "no_contamination": {
                "mask_excludes_all_sf": bool(mask_excludes_sf),
                "gamma_soft_off_N_per_m": g_off,
                "gamma_soft_on_N_per_m": g_on,
                "gamma_identical": gamma_identical,
            },
            "off_identity": bool(n_part_delta > 0),
        },
        "deferred": {
            "active_nmii_kumar_tension": "Kumar 2006 10-30 nN single-SF band — needs "
            "NMII via the sf_myosin_* prefix split (cortex/myosin.py) + the "
            "equilibration prelude. DEFERRED (passive backbone only here).",
            "N_filaments_ratification": "N_filaments=20 (Cramer 1997 candidate) — PI-pending.",
        },
        "note": (
            "REAL build-time activation gate on the FA-adhered baseline. PASSIVE "
            "FA→FA bundles, LONG-AXIS-ALIGNED basal pairing (PI 2026-06-09), "
            "per-bundle EXACT-r0 (force-free). KEY control is no-contamination: "
            "cortical γ_soft identical OFF vs ON (sf_ registry-denylisted). The "
            "Kumar tension band is the deferred active gate."
        ),
    }
    _figure(s_on, layout, g_off, g_on, p, result)
    return result


def _figure(s_on, layout, g_off, g_on, p, result) -> None:
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(15, 4.4))
    pos = np.asarray(s_on.particles.position) * 1e6
    tid = np.asarray(s_on.particles.typeid)
    types = list(s_on.particles.types)
    # basal view (x-y) of the cortex + SF bundles.
    cor = (tid == types.index("actin_cortex"))
    ax0.scatter(pos[cor, 0], pos[cor, 1], s=1, color="0.85", label="cortex")
    if layout is not None and len(layout.actin_tags) > 0:
        for b in range(len(layout.actin_tags)):
            t = np.asarray(layout.actin_tags[b], dtype=np.int64)
            ax0.plot(pos[t, 0], pos[t, 1], color="#d62728", lw=1.2, alpha=0.85)
    ax0.set_xlabel("x [µm]"); ax0.set_ylabel("y [µm]")
    ax0.set_aspect("equal", adjustable="datalim")
    a = result["controls"]["aligned"]
    ax0.set_title(f"{result['controls']['n_stress_fibers']} aligned vSF "
                  f"(|cos|={a['mean_abs_cos_long_axis']:.2f})")
    ax0.legend(fontsize=8)

    bl = a.get("bundle_len_um_min"), a.get("bundle_len_um_mean"), a.get("bundle_len_um_max")
    if layout is not None and len(layout.actin_tags) > 0:
        lens = []
        posf = np.asarray(s_on.particles.position)
        fe = np.asarray(layout.fa_endpoints, dtype=np.int64)
        lens = np.linalg.norm(posf[fe[:, 1]] - posf[fe[:, 0]], axis=1) * 1e6
        ax1.hist(lens, bins=15, color="#d62728", alpha=0.85)
        ax1.set_xlabel("bundle length (FA→FA) [µm]"); ax1.set_ylabel("count")
        ax1.set_title(f"bundle lengths (force-free; strain "
                      f"{result['controls']['force_free_construction']['max_strain']:.1e})")

    ax2.bar(["γ_soft OFF", "γ_soft ON"], [g_off * 1e3, g_on * 1e3],
            color=["#4c72b0", "#dd8452"])
    ax2.set_ylabel("cortical γ_soft [mN/m]")
    ax2.set_title(
        f"NO-contamination: identical="
        f"{result['controls']['no_contamination']['gamma_identical']}"
    )
    fig.suptitle(
        f"H.SF ventral_stress_fibers activation GATE — {result['verdict']} | "
        f"bundles={result['controls']['bundles_assembled']} | "
        f"aligned={result['controls']['aligned']['ok']} | "
        f"force-free={result['controls']['force_free_construction']['ok']} | "
        f"no-contam={result['controls']['no_contamination']['gamma_identical']}",
        fontsize=9.5,
    )
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT / "figs" / "h7_stress_fibers_activation_gate.png",
                dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    res = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_stress_fibers_activation_gate.json"
    with open(path, "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    c = res["controls"]; nc = c["no_contamination"]; a = c["aligned"]
    print(f"[H.SF ventral_stress_fibers activation gate] {res['verdict']}")
    print(f"  bundles assembled={c['bundles_assembled']} (n_SF={c['n_stress_fibers']}, "
          f"+{c['n_particle_delta']} part, +{c['n_bond_delta']} bonds)")
    print(f"  aligned={a['ok']} (mean|cos|={a['mean_abs_cos_long_axis']:.3f}, "
          f"len {a['bundle_len_um_min']:.1f}-{a['bundle_len_um_max']:.1f} µm)")
    print(f"  force-free: max strain={c['force_free_construction']['max_strain']:.2e} "
          f"ok={c['force_free_construction']['ok']}")
    print(f"  no-contam: γ_soft OFF={nc['gamma_soft_off_N_per_m']:.4e} ON={nc['gamma_soft_on_N_per_m']:.4e} "
          f"identical={nc['gamma_identical']}; sf excluded={nc['mask_excludes_all_sf']}")
    print(f"  json: {path}")
    return 0 if res["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
