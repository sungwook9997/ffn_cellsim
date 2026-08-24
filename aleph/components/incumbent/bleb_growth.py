#!/usr/bin/env python3
r"""Oracle-B multi-step blebbing — emergent bleb GROWTH + pore-pressure evolution after local de-adhesion.

Extends the single-shot nucleation criterion (:mod:`ac.cell.bleb_experiment`) to the mechanical outer loop.
It drives the composed cell through many physical steps with the SAME machinery the resting driver uses
(:class:`ac.fluid.scheduler.PhysicalScheduler` + :func:`ac.cell.driver.make_inner_solve`, imported read-only —
no shared-driver edit), and reads three emergent Charras (2005/2008) signatures from the assembled physics:

  1. **patch bulge** — the bare de-adhered membrane cap moves outward under the resting pore pressure while the
     intact (tethered) membrane stays put; the readout is the baseline-subtracted outward displacement of the
     patch nodes MINUS that of the control nodes;
  2. **local ΔP non-equilibration** — the mean pore ΔP over the cap fluid cells vs. the far-field, tracked per
     step (the Charras pressure-transient signature);
  3. **emergent Laplace balance** — the cap radius ``a_patch`` vs. the measured critical radius
     ``a_crit(t) = 2 γ_mem / ΔP(t)``; the bleb is expected to expand while ``a_patch > a_crit``.

Nothing here scripts the bleb: the bulge, its expansion, and the pressure transient are outputs of the
membrane bending + tension + Biot pressure physics. This module supplies only the *perturbation*
(:mod:`ac.cell.bleb_perturbation`) and the *diagnostics*.

Two fidelity contracts distinguish this from a naïve growth loop:

* **Physiological baseline (PI 2026-06-04, HARD).** The bleb must nucleate FROM a converged resting cell, not
  from a freshly-built, un-relaxed shell. The driver first equilibrates for ``equilibrate_steps`` outer steps;
  if ANY equilibration step is rejected — i.e. resting mechanics did not converge — it HALTS and reports
  ``resting_converged=False`` / ``blocked_on="resting-convergence"`` WITHOUT de-adhering. This is the explicit
  gate on the (b) full-native resting convergence: until that lands, this run stops here by construction rather
  than measuring a bleb off an unphysical baseline.
* **Commit-safe de-adhesion.** Local ERM ablation is irreversible biology, so — exactly like the ERM Bell KMC
  and the nuclear-envelope rupture — it may mutate ``erm_bound`` only at an ACCEPTED outer physical-time
  boundary. The de-adhesion kernel is launched from inside the scheduler's ``commit_irreversible`` hook,
  predicated on the authoritative device acceptance scalar, and fires exactly once (idempotent replay is safe:
  ``1→0`` monotone). A rejected mechanical candidate never leaves a bleb behind.

Runtime: Warp-CUDA only (I0-A). The device run executes on the gbook A5000; the dev Mac can only import this
module and exercise the pure-NumPy analysis helpers (:func:`cap_membership`, :func:`patch_outward_um`,
:func:`laplace_state`, :func:`bulge_grows`) — the same functions the device path calls on its post-step host
reads — via the CPU reference gate ``tests/ac/cell/test_bleb_growth_reference.py``.

    ~/miniconda3/envs/ffn_sim/bin/python -m aleph.components.incumbent.bleb_growth \
        --n-filaments 70686 --patch-deg 30 --equilibrate-steps 8 --growth-steps 40 \
        --inner-solver contact_tournament --max-inner-retries 4 \
        --out aleph/outputs/ac/bleb/growth.npz
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.components.incumbent.assemble import PI_0_PA, CellConfig, build_cell
from aleph.components.incumbent.bleb_experiment import _cap_pressure_jump
from aleph.components.incumbent.bleb_perturbation import cap_cos_half_angle, erm_patch_deadhesion_kernel
from aleph.components.incumbent.driver import make_inner_solve
from aleph.components.fluid.scheduler import PhysicalScheduler

__all__ = [
    "cap_membership",
    "patch_outward_um",
    "laplace_state",
    "bulge_grows",
    "run_bleb_growth",
]


# ── pure analysis helpers (host-side; shared by the device path and the CPU reference gate) ────────────

def cap_membership(
    membrane_pos: npt.NDArray[np.float64],
    centroid: npt.NDArray[np.float64],
    axis: npt.NDArray[np.float64],
    cos_half: float,
) -> tuple[npt.NDArray[np.bool_], npt.NDArray[np.float64]]:
    """Return ``(in_cap, radius)`` for each membrane node relative to ``centroid`` along a unit ``axis``.

    ``in_cap`` is the boolean angular-cap membership ``dot((x−c)/|x−c|, axis) ≥ cos_half`` (the SAME test the
    device de-adhesion kernel applies), and ``radius`` is ``|x − centroid|`` in µm. Kept pure so the CPU gate
    can validate the analysis on analytic geometry without a CUDA device.
    """
    r = np.asarray(membrane_pos, np.float64) - np.asarray(centroid, np.float64)
    d = np.linalg.norm(r, axis=1)
    cos = np.where(d > 1e-12, (r @ np.asarray(axis, np.float64)) / np.maximum(d, 1e-12), -2.0)
    return cos >= float(cos_half), d


def patch_outward_um(radius_now: npt.NDArray[np.float64], base_radius_um: float,
                     mask: npt.NDArray[np.bool_]) -> float:
    """Mean baseline-subtracted outward radial displacement over the masked nodes [µm] (0 if the mask is empty)."""
    mask = np.asarray(mask, bool)
    if not mask.any():
        return 0.0
    return float(np.asarray(radius_now, np.float64)[mask].mean() - float(base_radius_um))


def laplace_state(gamma_mem_pn_um: float, dp_pa: float, a_patch_um: float) -> tuple[float, float, bool]:
    """Emergent Laplace balance for the bare cap.

    Units are FF-engine consistent: ``γ_mem`` [pN/µm], ``ΔP`` [Pa] = [pN/µm²], ``a`` [µm], so
    ``a_crit = 2 γ_mem / ΔP`` [µm] and ``ΔP_Laplace = 2 γ_mem / a_patch`` [Pa]. Returns
    ``(a_crit_um, dP_laplace_Pa, expands)`` where ``expands`` is ``a_patch > a_crit`` — the cap bulges when its
    radius exceeds the pressure-set critical radius (Charras nucleation/growth criterion). Non-positive ΔP or
    ``a_patch`` yield ``inf`` and ``expands=False``.
    """
    a_crit = (2.0 * gamma_mem_pn_um / dp_pa) if dp_pa > 0 else float("inf")
    dp_laplace = (2.0 * gamma_mem_pn_um / a_patch_um) if a_patch_um > 0 else float("inf")
    return float(a_crit), float(dp_laplace), bool(a_patch_um > a_crit)


def bulge_grows(patch_disp: npt.ArrayLike, ctrl_disp: npt.ArrayLike) -> tuple[float, bool]:
    """Final ``(bulge, grows)`` from the post-de-adhesion series.

    ``patch_disp`` / ``ctrl_disp`` are the per-growth-step outward displacements (index 0 == the de-adhesion
    step). ``bulge`` is ``patch[-1] − ctrl[-1]`` [µm]; ``grows`` requires the patch to end above the control AND
    above its own de-adhesion-step value (i.e. it expanded outward, not merely started detached).
    """
    p = np.asarray(patch_disp, np.float64)
    c = np.asarray(ctrl_disp, np.float64)
    if p.size == 0:
        return 0.0, False
    bulge = float(p[-1] - c[-1])
    grows = bool(bulge > 0.0 and p[-1] > p[0])
    return bulge, grows


# ── device driver (Warp-CUDA; gbook A5000) ─────────────────────────────────────────────────────────────

def run_bleb_growth(cfg: CellConfig, *, patch_deg: float = 30.0, equilibrate_steps: int = 8,
                    growth_steps: int = 40, dt_phys: float = 0.05, n_inner: int = 60,
                    reshape_every: int = 20, max_inner_retries: int = 0, inner_solver: str = "explicit",
                    implicit_cg_max_iterations: int = 32, implicit_line_search_steps: int = 4,
                    implicit_coarse_iterations: int = 8, rkc_stages: int = 4,
                    biot_cfl_safety: float = 0.9, axis: npt.NDArray[np.float64] | None = None,
                    out_npz: str | None = None) -> dict:
    """Equilibrate the composed cell to rest, then de-adhere a cap commit-safely and track the emergent bleb.

    Protocol: (1) equilibrate ``equilibrate_steps`` outer steps — if any is rejected the resting state did not
    converge and the run HALTS on that gate; (2) measure the patch/control baseline FROM the equilibrated
    geometry; (3) arm the commit-safe de-adhesion and advance ``growth_steps`` steps, recording the patch bulge,
    cap ΔP, and Laplace balance per step. The de-adhesion fires only at an accepted outer boundary.
    """
    t0 = time.time()
    cell = build_cell(cfg)
    dev = cell.device
    mem, grid = cell.membrane, cell.grid
    if mem is None or grid is None or cell.substrate is None or cell.domain is None:
        raise RuntimeError("bleb growth needs membrane + Biot substrate + domain (full composed cell).")
    if mem.n_erm == 0:
        raise RuntimeError("no ERM tethers on the membrane — cannot de-adhere.")

    ax = np.array([0.0, 0.0, 1.0]) if axis is None else np.asarray(axis, np.float64)
    ax = ax / np.linalg.norm(ax)
    ax_v = wp.vec3d(*ax)
    cos_half = cap_cos_half_angle(patch_deg)
    cos_v = wp.float64(cos_half)
    p_ext = float(getattr(cell.substrate, "p_ext", 0.0))
    mem_idx = mem.erm_m_d.numpy()

    # ── commit-safe de-adhesion: wrap the driver's own irreversible-commit hook so the ablation is predicated
    # on the SAME device acceptance scalar the ERM Bell KMC / nucleus rupture use. Never edits driver.py. ──
    inner_solve = make_inner_solve(
        cell, n_inner, reshape_every, max_inner_retries, inner_solver=inner_solver,
        implicit_cg_max_iterations=implicit_cg_max_iterations,
        implicit_line_search_steps=implicit_line_search_steps,
        implicit_coarse_iterations=implicit_coarse_iterations, rkc_stages=rkc_stages)
    orig_commit = inner_solve.commit_irreversible
    deadh = {"armed": False, "done": False, "centroid": wp.vec3d(0.0, 0.0, 0.0)}

    def commit_with_deadhesion(accepted_d: wp.array) -> None:
        orig_commit(accepted_d)
        if deadh["armed"] and not deadh["done"]:
            # accepted_d is the authoritative outer-acceptance latch; a rejected step is a device no-op.
            wp.launch(erm_patch_deadhesion_kernel, dim=mem.n_erm,
                      inputs=[cell.pos_d, mem.erm_m_d, mem.erm_bound_d,
                              deadh["centroid"], ax_v, cos_v, accepted_d], device=dev)

    inner_solve.commit_irreversible = commit_with_deadhesion  # type: ignore[attr-defined]

    sched = PhysicalScheduler(
        substrate=cell.substrate, domain=cell.domain, membrane_bc=cell.membrane_bc,
        inner_solve=inner_solve, osmotic_difference=lambda t: PI_0_PA)

    # ── 1. equilibrate to the resting baseline; the (b) convergence gate lives here ─────────────────────
    eq_accepted, eq_residual = [], []
    for _ in range(equilibrate_steps):
        rep = sched.outer_step(dt_phys, biot_cfl_safety=biot_cfl_safety).readback()
        eq_accepted.append(bool(rep.outer_accepted))
        eq_residual.append(float(rep.inner_residual))
        if not rep.outer_accepted:
            break
    resting_converged = bool(eq_accepted) and all(eq_accepted)
    if not resting_converged:
        report = {
            "n_filaments": cfg.n_filaments, "patch_half_angle_deg": patch_deg,
            "equilibrate_steps": equilibrate_steps, "growth_steps": growth_steps,
            "inner_solver": inner_solver, "max_inner_retries": max_inner_retries, "n_inner": n_inner,
            "resting_converged": False, "blocked_on": "resting-convergence",
            "equilibrate_accepted": eq_accepted, "equilibrate_last_residual_pN": eq_residual[-1],
            "bleb_deadhered": False, "wall_s": time.time() - t0,
        }
        if out_npz:
            _dump(out_npz, report, {})
        return report

    # ── 2. baseline geometry FROM the equilibrated (resting) state ──────────────────────────────────────
    pos_eq = cell.pos_d.numpy()
    centroid = pos_eq[: cell.n_actin].mean(axis=0)
    deadh["centroid"] = wp.vec3d(*centroid)
    r_eq = np.linalg.norm(pos_eq[mem_idx] - centroid, axis=1)
    patch, _ = cap_membership(pos_eq[mem_idx], centroid, ax, cos_half)
    ctrl = ~patch
    base_patch_r = float(r_eq[patch].mean()) if patch.any() else 0.0
    base_ctrl_r = float(r_eq[ctrl].mean())
    R_cell = float(r_eq.mean())
    a_patch = R_cell * float(np.sin(np.deg2rad(patch_deg)))       # geodesic cap radius [µm]
    bound_pre = int(mem.erm_bound_d.numpy().sum())

    # ── 3. arm the commit-safe de-adhesion and track the emergent bleb ──────────────────────────────────
    deadh["armed"] = True
    ts, patch_disp, ctrl_disp, dp_cap = [], [], [], []
    a_crit_series, expands_series, accepted_flags = [], [], []
    deadhered_step = -1
    for step in range(growth_steps):
        rep = sched.outer_step(dt_phys, biot_cfl_safety=biot_cfl_safety).readback()
        if deadh["armed"] and not deadh["done"] and rep.outer_accepted:
            deadh["done"] = True                                 # one-shot: the accepted step committed it
            deadhered_step = step
        pos = cell.pos_d.numpy()
        r = np.linalg.norm(pos[mem_idx] - centroid, axis=1)
        patch_disp.append(patch_outward_um(r, base_patch_r, patch))
        ctrl_disp.append(patch_outward_um(r, base_ctrl_r, ctrl))
        dp_mean, _, _ = _cap_pressure_jump(grid, centroid, ax, cos_half, p_ext)
        dp_cap.append(dp_mean)
        a_crit, _, expands = laplace_state(mem.gamma_mem, dp_mean, a_patch)
        a_crit_series.append(a_crit)
        expands_series.append(expands)
        ts.append(float(rep.t))
        accepted_flags.append(bool(rep.outer_accepted))

    bound_post = int(mem.erm_bound_d.numpy().sum())
    n_deadhered = bound_pre - bound_post
    bulge, grows = bulge_grows(patch_disp, ctrl_disp)
    report = {
        "n_filaments": cfg.n_filaments, "patch_half_angle_deg": patch_deg,
        "equilibrate_steps": equilibrate_steps, "growth_steps": growth_steps, "dt_phys": dt_phys,
        "n_inner": n_inner, "inner_solver": inner_solver, "max_inner_retries": max_inner_retries,
        "resting_converged": True, "blocked_on": None,
        "n_tethers_patch": int(patch.sum()), "n_tethers_control": int(ctrl.sum()),
        "bleb_deadhered": bool(deadhered_step >= 0), "deadhered_at_growth_step": deadhered_step,
        "n_tethers_deadhered": int(n_deadhered), "erm_bound_pre": bound_pre, "erm_bound_post": bound_post,
        "R_cell_um": R_cell, "a_patch_um": a_patch, "gamma_mem_pN_um": mem.gamma_mem, "k_erm_pN_um": mem.k_erm,
        "base_patch_r_um": base_patch_r, "base_ctrl_r_um": base_ctrl_r,
        "final_patch_outward_um": float(np.asarray(patch_disp)[-1]) if patch_disp else 0.0,
        "final_control_outward_um": float(np.asarray(ctrl_disp)[-1]) if ctrl_disp else 0.0,
        "final_bulge_patch_minus_control_um": bulge, "bleb_grows": grows,
        "dp_cap_first_Pa": float(dp_cap[0]) if dp_cap else 0.0,
        "dp_cap_last_Pa": float(dp_cap[-1]) if dp_cap else 0.0,
        "a_crit_last_um": float(a_crit_series[-1]) if a_crit_series else float("inf"),
        "cap_expands_last": bool(expands_series[-1]) if expands_series else False,
        "all_growth_steps_accepted": bool(all(accepted_flags)) if accepted_flags else False,
        "wall_s": time.time() - t0,
    }
    if out_npz:
        _dump(out_npz, report, {
            "t": np.array(ts), "patch_disp": np.array(patch_disp), "ctrl_disp": np.array(ctrl_disp),
            "dp_cap": np.array(dp_cap), "a_crit": np.array(a_crit_series),
            "a_patch": np.full(len(ts), a_patch), "deadhered_step": deadhered_step,
        })
    return report


def _dump(out_npz: str, report: dict, arrays: dict) -> None:
    """Persist the growth series (npz) + the report (json sidecar)."""
    np.savez_compressed(out_npz, report_json=np.array(json.dumps(report), dtype=object), **arrays)
    with open(str(out_npz).rsplit(".npz", 1)[0] + ".json", "w") as fh:
        json.dump(report, fh, indent=2, default=float)


def _print(r: dict) -> None:
    print("\n=== Oracle-B multi-step bleb growth (equilibrate -> commit-safe de-adhesion -> emergent bulge) ===",
          flush=True)
    if not r.get("resting_converged", False):
        print(f"  n_filaments={r['n_filaments']:,}  inner_solver={r['inner_solver']} "
              f"(retries={r['max_inner_retries']})", flush=True)
        print(f"  equilibrate accepted={r['equilibrate_accepted']}  "
              f"last_residual={r['equilibrate_last_residual_pN']:.4g} pN", flush=True)
        print("  ==> BLOCKED ON resting-convergence: no de-adhesion applied. "
              "Re-run once (b) full-native resting converges.", flush=True)
        print(f"      (wall={r['wall_s']:.1f}s)", flush=True)
        return
    print(f"  n_filaments={r['n_filaments']:,}  patch θ½={r['patch_half_angle_deg']}°  "
          f"solver={r['inner_solver']}(retries={r['max_inner_retries']})  "
          f"growth_steps={r['growth_steps']} (dt={r['dt_phys']}s)", flush=True)
    print(f"  resting converged; de-adhered {r['n_tethers_deadhered']} tethers at growth step "
          f"{r['deadhered_at_growth_step']} (patch tethers={r['n_tethers_patch']}, "
          f"control={r['n_tethers_control']})", flush=True)
    print(f"  final outward displacement: patch={r['final_patch_outward_um']:.3e} µm  "
          f"control={r['final_control_outward_um']:.3e} µm  ->  bulge(patch−ctrl)="
          f"{r['final_bulge_patch_minus_control_um']:.3e} µm", flush=True)
    print(f"  cap ΔP: first={r['dp_cap_first_Pa']:.3g} Pa  last={r['dp_cap_last_Pa']:.3g} Pa;  "
          f"a_patch={r['a_patch_um']:.3g} µm vs a_crit(last)={r['a_crit_last_um']:.3g} µm "
          f"(expands={r['cap_expands_last']})", flush=True)
    print(f"  ==> BLEB GROWS (patch bulges beyond control): {r['bleb_grows']}   (wall={r['wall_s']:.1f}s)",
          flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description="Oracle-B multi-step bleb growth after commit-safe ERM de-adhesion.")
    p.add_argument("--n-filaments", type=int, default=70686)
    p.add_argument("--patch-deg", type=float, default=30.0)
    p.add_argument("--equilibrate-steps", type=int, default=8, help="resting-relaxation steps before de-adhesion")
    p.add_argument("--growth-steps", type=int, default=40)
    p.add_argument("--deadhere-at", type=int, default=None,
                   help="deprecated alias — de-adhesion now fires on the first accepted growth step")
    p.add_argument("--dt-phys", type=float, default=0.05)
    p.add_argument("--n-inner", type=int, default=60)
    p.add_argument("--reshape-every", type=int, default=20)
    p.add_argument("--max-inner-retries", type=int, default=0)
    p.add_argument("--inner-solver", type=str, default="explicit")
    p.add_argument("--implicit-cg-max-iterations", type=int, default=32)
    p.add_argument("--implicit-line-search-steps", type=int, default=4)
    p.add_argument("--implicit-coarse-iterations", type=int, default=8)
    p.add_argument("--rkc-stages", type=int, default=4)
    p.add_argument("--biot-cfl-safety", type=float, default=0.9)
    p.add_argument("--no-steric", action="store_true")
    p.add_argument("--no-myosin", action="store_true")
    p.add_argument("--device", type=str, default=None)
    p.add_argument("--out", type=str, default="")
    a = p.parse_args()
    if a.deadhere_at is not None:
        print("[bleb-growth] --deadhere-at is deprecated and ignored; de-adhesion fires on the first accepted "
              "growth step (commit-safe).", flush=True)
    cfg = CellConfig(n_filaments=a.n_filaments, with_myosin=not a.no_myosin, with_steric=not a.no_steric,
                     with_membrane=True, with_pressure=True, device=a.device)
    print(f"[bleb-growth] building cell n_filaments={cfg.n_filaments}", flush=True)
    r = run_bleb_growth(
        cfg, patch_deg=a.patch_deg, equilibrate_steps=a.equilibrate_steps, growth_steps=a.growth_steps,
        dt_phys=a.dt_phys, n_inner=a.n_inner, reshape_every=a.reshape_every,
        max_inner_retries=a.max_inner_retries, inner_solver=a.inner_solver,
        implicit_cg_max_iterations=a.implicit_cg_max_iterations,
        implicit_line_search_steps=a.implicit_line_search_steps,
        implicit_coarse_iterations=a.implicit_coarse_iterations, rkc_stages=a.rkc_stages,
        biot_cfl_safety=a.biot_cfl_safety, out_npz=a.out or None)
    _print(r)
    if a.out:
        print(f"[bleb-growth] wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()
