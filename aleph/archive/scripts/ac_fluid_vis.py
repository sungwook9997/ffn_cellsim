"""Regenerate the fluid-spine (I1a) analytic-oracle figures — the ac/ visualization pattern-setter.

The single regeneration entry-point for the fluid-spine track (extend it as I1b/I1c oracles land).
Every figure OVERLAYS the closed-form oracle/band on the numeric result, annotates SI units, and never
truncates an axis (the professor's visualization-integrity rules). Pure numpy/matplotlib on the committed
oracle suite — runs on the dev Mac (no Warp/CUDA).

    python aleph/scripts/ac_fluid_vis.py    ->  aleph/outputs/ac/fluid-spine/figs/*.png

Uses the I0-B1 ledger closure (c_v = 50 um^2/s = mobility/S), so the figures read the real numbers.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aleph.components.fluid.consolidation_analytic import (
    terzaghi_degree_of_consolidation,
    terzaghi_excess_pressure,
)
from aleph.components.fluid.darcy_analytic import darcy_slab_flux, pore_fluid_velocity
from aleph.components.fluid.fv_reference import BiotFVReference
from aleph.components.fluid.greens_analytic import diffusion_greens_function, mean_square_radius
from aleph.components.fluid.manufactured import (
    cosine_mode,
    cosine_mode_decay_rate,
    moving_interval_content,
    moving_interval_content_rate,
)
from aleph.components.fluid.transport_analytic import (
    frap_cosine_recovery_rate,
    frap_soumpasis_recovery,
)
from aleph.components.fluid.transport_reference import RADTransportReference, rad_cfl_dt

D_C = 3.0  # um^2/s  (I0-B1c draft, crowded cytoplasm — KB-DRAFT-7-03)

C_V = 50.0        # um^2/s   (KB-3.B3.2 anchor; = mobility/S)
MOBILITY = 5.0e-3  # um^2/(Pa*s)
STORAGE_S = 1.0e-4  # 1/Pa
OUT = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "fluid-spine" / "figs"

# Terzaghi's tabulated U(T_v) — the band overlaid on every consolidation gate.
TEXTBOOK_TV_U = np.array(
    [[0.008, 0.10], [0.031, 0.20], [0.071, 0.30], [0.126, 0.40], [0.197, 0.50],
     [0.287, 0.60], [0.403, 0.70], [0.567, 0.80], [0.848, 0.90]]
)


def fig_terzaghi() -> Path:
    """Degree of consolidation U(T_v) vs the textbook band + excess-pressure profiles."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    tv = np.linspace(0.0, 1.2, 600)
    ax0.plot(tv, terzaghi_degree_of_consolidation(tv), lw=2, label="oracle series U(T_v)")
    ax0.scatter(TEXTBOOK_TV_U[:, 0], TEXTBOOK_TV_U[:, 1], c="crimson", zorder=5,
                label="Terzaghi 1943 tabulated")
    ax0.axhline(0.5, ls=":", c="gray"); ax0.axvline(0.197, ls=":", c="gray")
    ax0.annotate("U=0.5 @ T_v≈0.197", (0.197, 0.5), (0.35, 0.35),
                 arrowprops=dict(arrowstyle="->", color="gray"))
    ax0.set_xlabel("time factor  T_v = c_v t / H²  [–]"); ax0.set_ylabel("degree of consolidation U  [–]")
    ax0.set_title("Consolidation U(T_v) vs textbook band"); ax0.set_ylim(0, 1); ax0.legend(loc="lower right")

    z = np.linspace(0.0, 2.0, 121)  # z/H across the doubly-drained layer
    for tvi in (0.02, 0.05, 0.15, 0.5, 1.0):
        ax1.plot(z, terzaghi_excess_pressure(z, tvi, u0=1.0), label=f"T_v={tvi}")
    ax1.set_xlabel("normalized depth  z/H  [–]"); ax1.set_ylabel("excess pressure  u/u₀  [–]")
    ax1.set_title("Excess-pressure dissipation (drained at z/H=0,2)"); ax1.set_ylim(0, 1.05); ax1.legend()
    fig.suptitle(f"I1a Terzaghi oracle   (c_v={C_V:.0f} µm²/s  →  τ_p=R²/c_v≈{7.5**2/C_V:.1f} s @ R=7.5 µm)")
    fig.tight_layout()
    p = OUT / "i1a_terzaghi.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_greens() -> Path:
    """Diffusion heat kernel spreading + the <r²>=2 d c_v t diffusion law (numeric vs analytic)."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    x = np.linspace(-60, 60, 601)
    for t in (0.1, 0.5, 2.0):
        ax0.plot(x, diffusion_greens_function(np.abs(x), t, C_V, dim=1),
                 label=f"t={t}s  (σ={np.sqrt(2*C_V*t):.1f} µm)")
    ax0.set_xlabel("x  [µm]"); ax0.set_ylabel("G₁ᴰ(x,t)  [1/µm]")
    ax0.set_title("1-D heat kernel — pore-pressure impulse spread"); ax0.legend()

    ts = np.linspace(0.05, 3.0, 40)
    r2_analytic = mean_square_radius(1.0, C_V, 3) * ts  # = 2*d*c_v*t, linear in t
    r2_numeric = []
    for t in ts:
        sig = np.sqrt(2 * C_V * t); r = np.linspace(0, 12 * sig, 8000)
        g = diffusion_greens_function(r, t, C_V, 3); w = g * 4 * np.pi * r**2
        r2_numeric.append(np.trapezoid(r**2 * w, r) / np.trapezoid(w, r))
    ax1.plot(ts, r2_analytic, lw=2, label="analytic  ⟨r²⟩ = 2·d·c_v·t (d=3)")
    ax1.scatter(ts[::4], np.array(r2_numeric)[::4], c="crimson", zorder=5, label="numeric ∫r²G dV")
    ax1.set_xlabel("time  t  [s]"); ax1.set_ylabel("mean-square radius ⟨r²⟩  [µm²]")
    ax1.set_title("Diffusion law (reads back c_v)"); ax1.legend(loc="upper left")
    fig.suptitle("I1a Green's-function oracle — conservation + diffusive spread")
    fig.tight_layout()
    p = OUT / "i1a_greens.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_manufactured() -> Path:
    """Cosine eigenmode decay vs exp(-λt) + the moving-boundary conservation identity vs FD."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    k = 1.0; lam = cosine_mode_decay_rate([k], C_V)
    t = np.linspace(0, 0.1, 200)
    amp = np.array([cosine_mode(np.array([0.0]), ti, 1.0, [k], C_V)[0] for ti in t])
    ax0.semilogy(t, np.abs(amp), lw=2, label="cosine eigenmode |p(0,t)/p₀|")
    ax0.semilogy(t, np.exp(-lam * t), "r--", label=f"exp(-λt), λ=c_v k²={lam:.0f} /s")
    ax0.set_xlabel("time  t  [s]"); ax0.set_ylabel("amplitude  [–]  (log axis)")
    ax0.set_title("Free eigenmode decay = operator eigenvalue"); ax0.legend()

    a, b, c, w = 30.0, 1.5, -0.4, 0.7; l0, eps, om = 6.0, 0.2, 0.9
    p = lambda x, tt: (a + b * x + c * x**2) * np.exp(-w * tt)
    px = lambda x, tt: (b + 2 * c * x) * np.exp(-w * tt)
    length = lambda tt: l0 * (1 + eps * np.sin(om * tt)); lrate = lambda tt: l0 * eps * om * np.cos(om * tt)
    src = lambda x, tt: STORAGE_S * (-w * p(x, tt)) - MOBILITY * (2 * c * np.exp(-w * tt))
    ts = np.linspace(0.05, 6.0, 120); dt = 1e-4
    ana = [moving_interval_content_rate(p, px, src, length, lrate, tt, STORAGE_S, MOBILITY)["total"] for tt in ts]
    fd = [(moving_interval_content(p, length, tt + dt, STORAGE_S)
           - moving_interval_content(p, length, tt - dt, STORAGE_S)) / (2 * dt) for tt in ts]
    ax1.plot(ts, ana, lw=2, label="analytic  dΦ/dt (moving-face + flux + source)")
    ax1.scatter(ts[::6], np.array(fd)[::6], c="crimson", zorder=5, label="finite-diff of Φ(t)")
    ax1.set_xlabel("time  t  [s]"); ax1.set_ylabel("content rate  dΦ/dt  [Pa·µm/s per unit area]")
    ax1.set_title("Moving-boundary conservation identity"); ax1.legend(loc="upper right")
    fig.suptitle("I1a manufactured-solution oracle — eigenmode + moving-domain conservation")
    fig.tight_layout()
    pth = OUT / "i1a_manufactured.png"; fig.savefig(pth, dpi=130); plt.close(fig)
    return pth


def fig_convergence() -> Path:
    """I1a discrete FV reference converges to Terzaghi at 2nd order in dx (numeric vs oracle)."""
    half_h = 5.0
    t_target = 0.15  # T_v = c_v t / H^2 = 50*0.15/25 = 0.3
    tv = C_V * t_target / half_h**2
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))

    grids = (40, 80, 160)
    errors = []
    for n in grids:
        dx = 2.0 * half_h / n
        ref = BiotFVReference((n,), dx=dx, mobility=MOBILITY, storage_S=STORAGE_S,
                              dirichlet={(0, 0): 0.0, (0, 1): 0.0})
        p = np.ones(n)
        dt = ref.cfl_dt(safety=0.5)
        n_steps = int(np.ceil(t_target / dt)); dt = t_target / n_steps
        for _ in range(n_steps):
            p = ref.step(p, dt)
        z = (np.arange(n) + 0.5) * dx / half_h
        ax0.plot(z, p, lw=1.2, alpha=0.85, label=f"FV reference n={n}")
        errors.append(float(np.sqrt(np.mean((p - terzaghi_excess_pressure(z, tv)) ** 2))))
    zf = np.linspace(0, 2, 400)
    ax0.plot(zf, terzaghi_excess_pressure(zf, tv), "k--", lw=2, label="Terzaghi oracle")
    ax0.set_xlabel("normalized depth  z/H  [–]"); ax0.set_ylabel("excess pressure  u/u₀  [–]")
    ax0.set_title(f"Discrete FV vs Terzaghi (T_v={tv:.2f})"); ax0.set_ylim(0, 1.02); ax0.legend()

    dxs = np.array([2.0 * half_h / n for n in grids])
    ax1.loglog(dxs, errors, "o-", lw=2, label="FV reference L2 error")
    ax1.loglog(dxs, errors[0] * (dxs / dxs[0]) ** 2, "r--", label="2nd-order slope ∝ dx²")
    ax1.set_xlabel("grid spacing  dx  [µm]"); ax1.set_ylabel("L2 error vs Terzaghi  [–]")
    ax1.set_title("Convergence order"); ax1.legend()
    fig.suptitle("I1a discrete conservative-FV convergence (the Warp kernel's CPU twin)")
    fig.tight_layout()
    p = OUT / "i1a_convergence.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_darcy() -> Path:
    """I1b Darcy kinematics: flux linear in Δp and mobility; v_f−v_s = q/φ (φ swept, I0-B1b GAP)."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    dps = np.linspace(0, 30, 50)
    for mob, lab in ((MOBILITY, "k/µ = 5e-3 µm²/(Pa·s)"), (3 * MOBILITY, "3× mobility")):
        q = np.array([darcy_slab_flux(dp, length=5.0, mobility=mob) for dp in dps])
        ax0.plot(dps, q, lw=2, label=f"q = (k/µ)Δp/L,  {lab}")
    ax0.set_xlabel("pressure drop  Δp  [Pa]"); ax0.set_ylabel("Darcy discharge  q  [µm/s]")
    ax0.set_title("Slab flux linear in Δp and in mobility"); ax0.legend(loc="upper left")

    q = np.linspace(0, 5e-3, 50)
    for phi in (0.5, 0.7, 0.9):
        vf_minus_vs = np.array([pore_fluid_velocity(0.0, qi, phi) for qi in q])
        ax1.plot(q, vf_minus_vs, lw=2, label=f"φ={phi}  (slope 1/φ)")
    ax1.set_xlabel("Darcy discharge  q  [µm/s]"); ax1.set_ylabel("v_f − v_s = q/φ  [µm/s]")
    ax1.set_title("Pore-velocity reconstruction (φ = I0-B1b GAP, swept)"); ax1.legend(loc="upper left")
    fig.suptitle("I1b mixed q-p form — the fluid FLOWS (kinematics vs analytic)")
    fig.tight_layout()
    p = OUT / "i1b_darcy.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_transport() -> Path:
    """I1c conservative RAD: total-actin conservation, FRAP↔D_c, advection follows v_f (not q)."""
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(15, 4.2))

    # Panel 0 — total-actin conservation to machine precision.
    rng = np.random.default_rng(0)
    shape = (16, 16, 16)
    ref = RADTransportReference(shape, dx=0.4, phi=0.7, d_c=D_C)
    u = np.abs(rng.standard_normal(shape)) * 0.7
    v_f = 0.3 * rng.standard_normal((*shape, 3)); reaction = 1e-2 * rng.standard_normal(shape)
    ref.bound_monomer = 5.0
    t0 = ref.total_actin(u); dt = rad_cfl_dt(ref.dx, D_C, float(np.max(np.abs(v_f))), 3)
    resid = []
    for s in range(200):
        u = ref.step(u, dt, v_f=v_f, reaction=reaction)
        resid.append(abs(ref.total_actin(u) - t0) / abs(t0))
    ax0.semilogy(np.arange(1, 201), np.clip(resid, 1e-18, None), lw=1.5)
    ax0.axhline(1e-11, ls="--", c="crimson", label="1e-11 gate")
    ax0.set_xlabel("step"); ax0.set_ylabel("|∫φc + N_polymer − total₀| / total₀  [–]")
    ax0.set_title("Total-actin conservation (adv+diff+reaction)"); ax0.legend()

    # Panel 1 — FRAP cosine recovery reads back D_c; Soumpasis literature curve overlaid.
    n = 200; L = 20.0; dx = L / n
    r1 = RADTransportReference((n,), dx=dx, phi=0.7, d_c=D_C)
    x = (np.arange(n) + 0.5) * dx; k = np.pi / L
    u1 = 0.7 * (1 + 0.3 * np.cos(k * x)); amp0 = 0.3 * 0.7
    dt1 = rad_cfl_dt(dx, D_C, 0.0, 1); ts, amps = [], []
    for s in range(6000):
        u1 = r1.step(u1, dt1); dev = u1 - u1.mean()
        ts.append((s + 1) * dt1); amps.append(0.5 * (dev.max() - dev.min()) / amp0)
    ts = np.array(ts)
    ax1.semilogy(ts, amps, lw=2, label="FV reference cosine bleach")
    lam = frap_cosine_recovery_rate([k], D_C)
    ax1.semilogy(ts, np.exp(-lam * ts), "r--", label=f"exp(−D_c k² t), D_c={D_C} µm²/s")
    ax1.set_xlabel("time  t  [s]"); ax1.set_ylabel("bleach amplitude  [–]  (log)")
    ax1.set_title("FRAP recovery reads back D_c"); ax1.legend()

    # Panel 2 — advection front follows v_f, NOT q = φ v_f.
    n = 400; dx = 0.25; phi = 0.7; v_f = 1.0
    r2 = RADTransportReference((n,), dx=dx, phi=phi, d_c=0.0)
    x = (np.arange(n) + 0.5) * dx
    u2 = np.exp(-((x - 20.0) ** 2) / (2 * (2 * dx) ** 2)); vf = np.full((n, 1), v_f)
    dt2 = rad_cfl_dt(dx, 0.0, v_f, 1); tt, cen = [], []
    for s in range(200):
        u2 = r2.step(u2, dt2, v_f=vf)
        tt.append((s + 1) * dt2); cen.append(float(np.sum(x * u2) / np.sum(u2)))
    tt = np.array(tt)
    ax2.plot(tt, cen, lw=2, label="FV reference blob centroid")
    ax2.plot(tt, 20.0 + v_f * tt, "r--", label="x₀ + v_f·t  (correct)")
    ax2.plot(tt, 20.0 + phi * v_f * tt, ":", c="gray", label="x₀ + q·t  (q=φv_f, WRONG)")
    ax2.set_xlabel("time  t  [s]"); ax2.set_ylabel("front position  [µm]")
    ax2.set_title("Advection rides v_f, not the discharge q"); ax2.legend(loc="upper left")

    fig.suptitle("I1c conservative G-actin transport — mass conservation + FRAP + advection velocity")
    fig.tight_layout()
    p = OUT / "i1c_transport.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    figs = [fig_terzaghi(), fig_greens(), fig_manufactured(),
            fig_convergence(), fig_darcy(), fig_transport()]
    for f in figs:
        print(f"  wrote {f.relative_to(Path(__file__).resolve().parents[2])}")


if __name__ == "__main__":
    main()
