#!/usr/bin/env python3
r"""NG-4 fluid-first PAYLOAD — drained<->undrained rate dependence + spatial poroelastic tau_p.

This is the concrete answer to "is the cytosol well-implemented, or is it a scalar-turgor balloon?".
It drives the CONSERVATIVE Biot p/mass substrate (``ac.fluid.biot_substrate``) on the live
membrane-minus-nucleus ``Domain`` with a *controlled solid dilatation* (the "prescribed dilatation"
option of NG-4) and MEASURES the two signatures a balloon cannot show:

  (b) Skempton — the UNDRAINED (fast) response is STIFFER than the DRAINED (slow) one: for the SAME total
      imposed compressive strain, a fast ramp builds a much higher LOCAL pore pressure (fluid has no time to
      redistribute) than a slow ramp (fluid drains internally as the load is applied).
  (c) spatial relaxation time tau_p ~ R^2/c_v: after a fast hemispherical load, the internal pore-pressure
      gradient RELAXES with a time constant set by the consolidation coefficient c_v = mobility/S and the
      domain size R. We verify it SCALES as 1/c_v (halve c_v -> tau doubles) — the poroelastic fingerprint.
  (a) the fluid FLOWS — grad p develops, the Darcy discharge q = -(k/mu) grad p is nonzero, and the absolute
      pore-fluid velocity v_f = v_s + q/phi differs from v_s (with a PROVISIONAL phi, LABELLED a GAP).

STERIC is irrelevant here BY CONSTRUCTION: this payload is pure fluid on the Eulerian field (no actin nodes,
no excluded volume), so the pN-scale poroelastic signal is fully isolated from the sigma_EV interpenetration
confound (~63k nodes, up to 1000 pN steric) that lives on the SOLID track.

A scalar-turgor balloon has ONE uniform pressure that adjusts instantly: no spatial gradient, no internal
flow, no rate dependence, tau=0. Any nonzero Skempton ratio + finite tau_p refutes the balloon.

Runtime: Warp-CUDA only (I0-A). Runs on the gbook A5000. Reads the fluid track modules read-only.

    ~/miniconda3/envs/ffn_sim/bin/python -m aleph.components.incumbent.ng4_payload \
        --out aleph/outputs/ac/ng4/ng4_payload.npz
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import warp as wp

from aleph.components.fluid.biot_substrate import BiotSubstrate
from aleph.components.fluid.boundary import MembraneFluxBC
from aleph.components.fluid.domain import Domain, StaticSphereMembraneProvider
from aleph.components.fluid.field_grid import FLUID, FieldGrid
from aleph.components.fluid.velocity import DarcyVelocity

__all__ = ["build_substrate", "run_ramp_relax", "ng4_demo"]

# ── physiological / I0-B1 closed anchors (engine units: um, pN, s; 1 Pa == 1 pN/um^2) ────────────────
R_CELL_UM = 7.5           # MCF7 radius
BIOT_STORAGE_S = 1.0e-4   # S = 1/M [um^2/pN]  (M ~ 1e4 Pa)                        [I0-B1 derived]
BIOT_MOBILITY = 5.0e-3    # k/mu [um^4/(pN*s)] so c_v = mobility/S = 50 um^2/s     [I0-B1 c_v anchor]
BIOT_ALPHA = 1.0          # Biot-Willis coupling                                    [PI-ratified]
L_P = 1.6e-8              # membrane hydraulic conductivity [um/(s*Pa)] (draft; ~sealed on tau_p scale)
PHI_PROVISIONAL = 0.5     # porosity for v_f = v_s + q/phi — *** I0-B1b GAP (unsourced) ***


def _cell_coords(grid: FieldGrid) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cell-centre coordinate meshes (host) for prescribing/selecting regions (out-of-hot-loop)."""
    nx, ny, nz = grid.shape
    ox, oy, oz = grid.origin
    xs = ox + np.arange(nx) * grid.dx
    ys = oy + np.arange(ny) * grid.dx
    zs = oz + np.arange(nz) * grid.dx
    return np.meshgrid(xs, ys, zs, indexing="ij")


def build_substrate(*, cv: float = 50.0, dx: float = 0.5, pad: float = 3.0, device: str | None = None):
    """Build (grid, domain, substrate, membrane_bc) for a spherical cell of radius R at consolidation cv.

    mobility is scaled to realise c_v = mobility/S = ``cv`` at the fixed storativity S (so a c_v sweep is a
    clean mobility change; S, alpha, R held). Seeds p = p_bar = 0 so the field IS the excess response.
    """
    mobility = cv * BIOT_STORAGE_S
    half = R_CELL_UM + pad
    nx = int(round(2.0 * half / dx)) + 1
    origin = (-half, -half, -half)
    grid = FieldGrid((nx, nx, nx), dx, origin, device=device)
    grid.set_pressure(np.zeros(grid.shape, np.float64))         # p_excess starts at 0
    with wp.ScopedDevice(device):
        grid.p_bar = wp.zeros(grid.shape, dtype=wp.float64)
    membrane = StaticSphereMembraneProvider(centre=(0.0, 0.0, 0.0), radius=R_CELL_UM)
    domain = Domain(grid, membrane=membrane, nucleus=None)
    domain.classify()
    substrate = BiotSubstrate(grid, mobility=mobility, storage_S=BIOT_STORAGE_S, alpha=BIOT_ALPHA)
    membrane_bc = MembraneFluxBC(grid, L_p=L_P, p_ext=0.0, sigma_refl=1.0)
    return grid, domain, substrate, membrane_bc


def _hemisphere_div_vs(grid: FieldGrid, eps_rate: float) -> np.ndarray:
    """Prescribe div(v_s) = -eps_rate (compression) in the TOP (z>0) fluid hemisphere, 0 elsewhere [1/s].

    A negative div(v_s) is a volumetric CONTRACTION of the skeleton: it squeezes pore fluid, so the
    -alpha*div(v_s) term is a POSITIVE pressure source in the compressed region.
    """
    _, _, Z = _cell_coords(grid)
    mask = grid.mask.numpy()
    div = np.zeros(grid.shape, np.float64)
    sel = (mask == FLUID) & (Z > 0.0)
    div[sel] = -eps_rate
    return div


def run_ramp_relax(
    grid, domain, substrate, membrane_bc, *, eps_rate: float, t_ramp: float, t_relax: float,
    with_membrane: bool = True, safety: float = 0.9, record_every: int = 25, flow_phi: float | None = None,
) -> dict:
    """Ramp a hemispherical compression at rate ``eps_rate`` for ``t_ramp`` s, then hold + relax ``t_relax`` s.

    Records (out-of-hot-loop) the top/bottom-hemisphere mean pressure, max |p|, total content, and the peak
    inter-hemisphere pressure difference. Returns time series + final field + measured decay diagnostics.
    """
    g = grid
    _, _, Z = _cell_coords(g)
    fluid = g.mask.numpy() == FLUID
    top = fluid & (Z > 1.0)                 # avoid the mid-plane + membrane skin
    bot = fluid & (Z < -1.0)
    div_host = _hemisphere_div_vs(g, eps_rate)

    dt = substrate.cfl_dt(safety)
    n_ramp = max(1, int(round(t_ramp / dt)))
    n_relax = max(1, int(round(t_relax / dt)))
    dt = t_ramp / n_ramp                     # exact fit so the imposed strain is deterministic

    def _snapshot(t: float) -> tuple:
        p = g.p.numpy()
        return (t, float(p[top].mean()), float(p[bot].mean()),
                float(np.abs(p[fluid]).max()), float(BIOT_STORAGE_S * p[fluid].sum() * g.dx**3))

    ts: list[tuple] = []
    # ── RAMP phase: div_vs prescribed (skeleton compressing) ─────────────────────────────────────────
    with wp.ScopedDevice(g.device):
        g.div_vs = wp.array(div_host, dtype=wp.float64)
    ts.append(_snapshot(0.0))
    t = 0.0
    for it in range(n_ramp):
        if with_membrane:
            membrane_bc.apply(0.0)           # resting osmotic balance; efflux ~ -L_p * p (near-sealed)
        substrate.step(dt)
        t += dt
        if (it + 1) % record_every == 0 or it == n_ramp - 1:
            ts.append(_snapshot(t))
    t_ramp_end = t
    eps_total = eps_rate * t_ramp_end

    # ── FLOW diagnostic at the PEAK-gradient (undrained) state == end of the ramp, BEFORE relaxation ─
    flow = _flow_diag(grid, substrate, flow_phi) if flow_phi is not None else None

    # ── RELAX phase: hold strain (div_vs = 0), let the internal gradient consolidate ─────────────────
    with wp.ScopedDevice(g.device):
        g.div_vs.zero_()
    for it in range(n_relax):
        if with_membrane:
            membrane_bc.apply(0.0)
        substrate.step(dt)
        t += dt
        if (it + 1) % record_every == 0 or it == n_relax - 1:
            ts.append(_snapshot(t))

    arr = np.array(ts)                       # (N, 5): t, p_top, p_bot, p_absmax, content
    t_arr, p_top, p_bot = arr[:, 0], arr[:, 1], arr[:, 2]
    dp = p_top - p_bot                       # inter-hemisphere pressure difference (the gradient signal)
    peak_dp = float(dp.max())
    peak_p_top = float(p_top.max())

    # ── measure tau: exponential decay of dp during the RELAX phase ──────────────────────────────────
    relax = t_arr >= t_ramp_end
    tr, dpr = t_arr[relax], dp[relax]
    tau_fit = float("nan")
    if dpr.size >= 4 and dpr[0] > 0:
        # fit toward the equilibrated floor (mean of the last few) so a nonzero asymptote doesn't bias tau
        floor = float(dpr[-3:].mean())
        y = dpr - floor
        good = y > 0.02 * (dpr[0] - floor + 1e-12)
        if good.sum() >= 3:
            slope = np.polyfit(tr[good] - tr[good][0], np.log(y[good]), 1)[0]
            if slope < 0:
                tau_fit = float(-1.0 / slope)

    return {
        "series": arr, "dt": dt, "n_ramp": n_ramp, "n_relax": n_relax,
        "eps_rate": eps_rate, "t_ramp": t_ramp_end, "eps_total": eps_total,
        "peak_dp_Pa": peak_dp, "peak_p_top_Pa": peak_p_top,
        "tau_relax_s": tau_fit, "cv": substrate.c_v, "final_p": g.p.numpy(), "flow": flow,
    }


def _flow_diag(grid, substrate, phi: float) -> dict:
    """Compute the Darcy discharge q on device and report the fluid-flow magnitudes (out-of-hot-loop)."""
    dv = DarcyVelocity(grid, mobility=substrate.mobility, phi=phi)
    dv.discharge()
    q = dv.q.numpy()                          # (nx,ny,nz,3) [um/s]
    fluid = grid.mask.numpy() == FLUID
    qmag = np.linalg.norm(q, axis=-1)[fluid]
    # v_f - v_s = q/phi (solid held during relax -> v_f magnitude ~ q/phi); grad p from q = -mobility grad p
    gradp = qmag / substrate.mobility         # |grad p| [Pa/um]
    return {"max_q_um_s": float(qmag.max()), "p95_q_um_s": float(np.percentile(qmag, 95)),
            "max_vf_minus_vs_um_s": float((qmag / phi).max()),
            "max_gradp_Pa_um": float(gradp.max()), "phi_provisional_GAP": phi}


def ng4_demo(*, dx: float = 0.5, device: str | None = None, out_npz: str | None = None,
             t_ramp_fast: float = 0.02, t_ramp_slow: float = 12.0, t_relax: float = 6.0,
             eps_total: float = 0.01, phi: float = PHI_PROVISIONAL) -> dict:
    """Run the full NG-4 payload: fast(undrained) vs slow(drained) at c_v=50, + a c_v-scaling check for tau_p."""
    t0 = time.time()
    R2_over_cv = R_CELL_UM**2 / 50.0          # nominal tau_p scale = R^2/c_v ~ 1.1 s
    report: dict = {"R2_over_cv_s": R2_over_cv, "phi_provisional_GAP": phi, "eps_total_target": eps_total}

    # ── (b) Skempton: fast vs slow, SAME total strain, c_v = 50 ──────────────────────────────────────
    g, dom, sub, mbc = build_substrate(cv=50.0, dx=dx, device=device)
    fast = run_ramp_relax(g, dom, sub, mbc, eps_rate=eps_total / t_ramp_fast,
                          t_ramp=t_ramp_fast, t_relax=t_relax, flow_phi=phi)
    fast_flow = fast["flow"]                   # flow measured AT the ramp end (peak-gradient undrained state)

    g2, dom2, sub2, mbc2 = build_substrate(cv=50.0, dx=dx, device=device)
    slow = run_ramp_relax(g2, dom2, sub2, mbc2, eps_rate=eps_total / t_ramp_slow,
                          t_ramp=t_ramp_slow, t_relax=t_relax)

    skempton_ratio_top = fast["peak_p_top_Pa"] / max(slow["peak_p_top_Pa"], 1e-12)
    skempton_ratio_dp = fast["peak_dp_Pa"] / max(slow["peak_dp_Pa"], 1e-12)
    report["skempton"] = {
        "fast_t_ramp_s": fast["t_ramp"], "slow_t_ramp_s": slow["t_ramp"],
        "fast_peak_p_top_Pa": fast["peak_p_top_Pa"], "slow_peak_p_top_Pa": slow["peak_p_top_Pa"],
        "fast_peak_dp_Pa": fast["peak_dp_Pa"], "slow_peak_dp_Pa": slow["peak_dp_Pa"],
        "undrained_over_drained_p_top": skempton_ratio_top,
        "undrained_over_drained_dp": skempton_ratio_dp,
        "eps_total_fast": fast["eps_total"], "eps_total_slow": slow["eps_total"],
    }

    # ── (c) tau_p scaling: c_v = 50 vs 25 (fast load), expect tau(25) ~ 2 * tau(50) ──────────────────
    g3, dom3, sub3, mbc3 = build_substrate(cv=25.0, dx=dx, device=device)
    half_cv = run_ramp_relax(g3, dom3, sub3, mbc3, eps_rate=eps_total / t_ramp_fast,
                             t_ramp=t_ramp_fast, t_relax=t_relax * 2.0)
    report["tau_p"] = {
        "tau_cv50_s": fast["tau_relax_s"], "tau_cv25_s": half_cv["tau_relax_s"],
        "ratio_tau25_over_tau50": (half_cv["tau_relax_s"] / fast["tau_relax_s"]
                                   if fast["tau_relax_s"] == fast["tau_relax_s"] else float("nan")),
        "expected_ratio": 2.0, "R2_over_cv50_s": R2_over_cv,
        "note": "fundamental Neumann l=1 mode has prefactor 1/(kR)^2~0.23, so tau ~ 0.23*R^2/c_v",
    }

    # ── (a) fluid FLOWS ──────────────────────────────────────────────────────────────────────────────
    report["flow_at_fast_undrained_state"] = fast_flow

    report["wall_s"] = time.time() - t0
    if out_npz:
        np.savez_compressed(
            out_npz,
            fast_series=fast["series"], slow_series=slow["series"], cv25_series=half_cv["series"],
            fast_final_p=fast["final_p"].astype(np.float32),
            slow_final_p=slow["final_p"].astype(np.float32),
            grid_shape=np.array(g.shape, np.int32), grid_dx=np.float64(g.dx),
            grid_origin=np.array(g.origin, np.float64), mask=g.mask.numpy().astype(np.int32),
            report_json=np.array(json.dumps(report), dtype=object),
        )
        with open(str(out_npz).rsplit(".npz", 1)[0] + ".json", "w") as fh:
            json.dump(report, fh, indent=2, default=float)
    return report


def _print(report: dict) -> None:
    print("\n=== NG-4 fluid-first payload ===", flush=True)
    print(f"  nominal tau_p = R^2/c_v = {report['R2_over_cv_s']:.3f} s  (R=7.5 um, c_v=50 um^2/s)", flush=True)
    sk = report["skempton"]
    print("\n(b) SKEMPTON — undrained(fast) vs drained(slow), SAME total strain "
          f"eps={report['eps_total_target']}:", flush=True)
    print(f"    fast  t_ramp={sk['fast_t_ramp_s']:.3g}s  peak p_top={sk['fast_peak_p_top_Pa']:.3f} Pa  "
          f"peak dp={sk['fast_peak_dp_Pa']:.3f} Pa", flush=True)
    print(f"    slow  t_ramp={sk['slow_t_ramp_s']:.3g}s  peak p_top={sk['slow_peak_p_top_Pa']:.3f} Pa  "
          f"peak dp={sk['slow_peak_dp_Pa']:.3f} Pa", flush=True)
    print(f"    undrained/drained  (p_top)={sk['undrained_over_drained_p_top']:.2f}x   "
          f"(dp)={sk['undrained_over_drained_dp']:.2f}x   <-- >1 == STIFFER when fast (NOT a balloon)",
          flush=True)
    tp = report["tau_p"]
    print("\n(c) tau_p spatial relaxation (fast load, then hold strain):", flush=True)
    print(f"    c_v=50 -> tau={tp['tau_cv50_s']:.3f} s ;  c_v=25 -> tau={tp['tau_cv25_s']:.3f} s ;  "
          f"ratio={tp['ratio_tau25_over_tau50']:.2f} (expect ~2 => tau ~ 1/c_v)", flush=True)
    print(f"    {tp['note']}", flush=True)
    fl = report["flow_at_fast_undrained_state"]
    print("\n(a) the fluid FLOWS (post fast ramp, solid held):", flush=True)
    print(f"    max|grad p|={fl['max_gradp_Pa_um']:.3f} Pa/um   max|q|={fl['max_q_um_s']:.3g} um/s   "
          f"max|v_f - v_s|=q/phi={fl['max_vf_minus_vs_um_s']:.3g} um/s  (phi={fl['phi_provisional_GAP']} GAP)",
          flush=True)
    print(f"\n  wall={report['wall_s']:.1f}s", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description="NG-4 fluid-first payload (Skempton + tau_p + flow).")
    p.add_argument("--dx", type=float, default=0.5)
    p.add_argument("--device", type=str, default=None, help="Warp CUDA device alias (default: current CUDA)")
    p.add_argument("--out", type=str, default="")
    p.add_argument("--t-ramp-fast", type=float, default=0.02)
    p.add_argument("--t-ramp-slow", type=float, default=12.0)
    p.add_argument("--t-relax", type=float, default=6.0)
    p.add_argument("--eps-total", type=float, default=0.01)
    p.add_argument("--phi", type=float, default=PHI_PROVISIONAL)
    a = p.parse_args()
    print(f"[ng4] building fluid substrate dx={a.dx} device={a.device}", flush=True)
    rep = ng4_demo(dx=a.dx, device=a.device, out_npz=a.out or None,
                   t_ramp_fast=a.t_ramp_fast, t_ramp_slow=a.t_ramp_slow, t_relax=a.t_relax,
                   eps_total=a.eps_total, phi=a.phi)
    _print(rep)
    if a.out:
        print(f"[ng4] wrote {a.out}", flush=True)


if __name__ == "__main__":
    main()
