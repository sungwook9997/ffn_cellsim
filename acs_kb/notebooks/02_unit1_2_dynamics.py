"""Phase 1 Unit 1.2 — overdamped Langevin dynamics walkthrough (preliminary).

Generates:
  - 04_free_diffusion.png             — MSD vs t for free beads (sanity)
  - 05_equipartition.png              — per-link energy histogram → ½ k_BT
  - 06_preliminary_strain_stiffening.png — σ_xy(γ), G(γ), G_0 estimate

Plus a JSON summary in acs_kb/outputs/phase1/02_summary.json.

**Strain-stiffening here is preliminary.** Per PI direction (2026-05-18),
quantitative comparison to KU-1.30 #1 (G_0 ≈ 36 Pa) and #2 (γ_c ≈
0.16) is deferred until either exact intersection-node insertion or
finer bead discretization is in place.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from acs_kb.common.derived_params import load_config
from acs_kb.common.sanity_gate import (
    gate_equipartition,
    gate_unit1_2_dynamics,
)
from acs_kb.ecm.cross_links import generate_cross_links
from acs_kb.ecm.diagnostics import link_extension_energies
from acs_kb.ecm.fiber_mechanics import compute_forces
from acs_kb.ecm.fiber_network import generate_2d_fiber_network
from acs_kb.ecm.integrator import EulerMaruyama, run
from acs_kb.ecm.shear_protocol import preliminary_affine_G_curve

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "phase1"
OUT.mkdir(parents=True, exist_ok=True)
CFG_PATH = ROOT / "configs" / "phase1_unit1.yaml"


def _free_diffusion_msd(gamma_b, kT, dt, n_steps_grid, rng):
    """Run n_walkers free beads, return (t, MSD, expected)."""
    n_walkers = 4000
    pos0 = np.zeros((n_walkers, 1, 2))
    p = pos0.copy()
    forces_fn = lambda r: np.zeros_like(r)
    integ = EulerMaruyama()
    msds, ts = [], []
    cum = 0
    for n in n_steps_grid:
        steps = n - cum
        for _ in range(steps):
            p = integ.step(p, forces_fn, gamma_b, kT, dt, rng)
        cum = n
        msds.append(float(np.mean(np.sum((p - pos0) ** 2, axis=2))))
        ts.append(cum * dt)
    D = kT / gamma_b
    expected = 4.0 * D * np.array(ts)
    return np.array(ts), np.array(msds), expected


def main() -> None:
    cfg = load_config(CFG_PATH)
    ecm = cfg["ecm"]; d = ecm["derived"]; dyn = ecm["dynamics"]
    gamma_b = d["gamma_b"]; kT = ecm["kT"]; dt = dyn["dt"]
    mu = ecm["stretching_modulus"]; kappa = ecm["bending_modulus"]
    rng = np.random.default_rng(0)
    integ = EulerMaruyama()

    # ---- CFL gate ----
    cfl = gate_unit1_2_dynamics(
        dt=dt, tau_xl=d["tau_xl"], tau_stretch=d["tau_stretch"],
        tau_bend=d["tau_bend"], safety_factor=dyn["cfl_safety_factor"],
        integrator_name=dyn["integrator"],
    )
    cfl.raise_if_failed()

    # ---- Free diffusion (Fig 04) ----
    t0 = time.perf_counter()
    n_steps_grid = [20, 50, 100, 200, 400]
    ts, msds, expected = _free_diffusion_msd(gamma_b, kT, dt, n_steps_grid, rng)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(ts * 1e9, msds * 1e18, "o-", label="EulerMaruyama (4000 walkers)",
            color="steelblue")
    ax.plot(ts * 1e9, expected * 1e18, "--", label="4 D t (D=k_BT/γ_b)",
            color="darkorange")
    ax.set_xlabel("t (ns)"); ax.set_ylabel(r"⟨|Δr|²⟩ (nm²)")
    ax.set_title("Free overdamped diffusion (KU-1.26 sanity)")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(OUT / "04_free_diffusion.png", dpi=150); plt.close(fig)
    t_diff = time.perf_counter() - t0

    # ---- Equipartition (Fig 05) ----
    t0 = time.perf_counter()
    net_eq = generate_2d_fiber_network(
        L_box=ecm["L_box"], n_fibers=400, L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"], S_order=0.0, seed=3,
    )
    links_eq = generate_cross_links(net_eq, stiffness=ecm["xl_stiffness"])
    def forces_eq(p):
        return compute_forces(p, net_eq.rest_length, mu, kappa, net_eq.box_size,
                              cross_links=links_eq)
    p = run(net_eq.bead_positions, forces_eq, integ, gamma_b, kT, dt,
            net_eq.box_size, n_steps=12000, rng=rng)
    samples = []
    for _ in range(200):
        p = run(p, forces_eq, integ, gamma_b, kT, dt, net_eq.box_size,
                n_steps=50, rng=rng)
        samples.extend(link_extension_energies(p, links_eq, net_eq.box_size).tolist())
    samples = np.array(samples)
    mean_e = float(np.mean(samples))
    target = 0.5 * kT
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(samples / target, bins=60, color="steelblue", alpha=0.7, density=True)
    ax.axvline(1.0, color="red", ls="--", label="½ k_BT (equipartition)")
    ax.axvline(mean_e / target, color="darkorange", ls="-",
               label=f"measured mean = {mean_e/target:.3f}")
    ax.set_xlabel("per-link energy / ½ k_BT")
    ax.set_ylabel("density")
    ax.set_title(f"Equipartition (KU-1.26): {len(links_eq)} links × {200} time samples")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(OUT / "05_equipartition.png", dpi=150); plt.close(fig)
    eq_gate = gate_equipartition(measured_mean_link_energy=mean_e, kT=kT,
                                 n_samples=len(samples), tolerance=0.12)
    t_eq = time.perf_counter() - t0

    # ---- Preliminary strain-stiffening (Fig 06) ----
    t0 = time.perf_counter()
    net_ss = generate_2d_fiber_network(
        L_box=ecm["L_box"], n_fibers=400, L_fiber=ecm["L_fiber"],
        beads_per_fiber=ecm["beads_per_fiber"], S_order=0.0, seed=7,
    )
    links_ss = generate_cross_links(net_ss, stiffness=ecm["xl_stiffness"])
    gammas = np.array([0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.18,
                       0.22, 0.26, 0.30, 0.35])
    h_layer = ecm.get("layer_thickness_for_2d_to_3d")
    sigma_2d, G_2d, G_3d = preliminary_affine_G_curve(
        net_ss, links_ss, mu, kappa, gammas, layer_thickness=h_layer,
    )
    # Sign convention for the virial: positive γ gives negative σ_xy in
    # our derivation (b_y is the off-axis coordinate, F_y has opposite
    # sign to extension along that direction). We plot |σ_xy| to keep
    # the log axes honest and annotate that sign separately.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].loglog(gammas, np.abs(sigma_2d), "o-", color="steelblue",
                   label="|σ_xy^2D| (N/m)")
    axes[0].set_xlabel("shear strain γ")
    axes[0].set_ylabel(r"$|\sigma_{xy}^{\rm 2D}|$ (N/m)")
    axes[0].set_title("Affine 2D virial (preliminary)")
    axes[0].legend(fontsize=8)

    axes[1].semilogx(gammas, np.abs(G_3d), "o-", color="darkorange",
                     label=f"|G_affine^3D| = |σ_2D|/h, h={h_layer*1e6:.1f} μm")
    axes[1].axhspan(15, 200, color="green", alpha=0.1,
                    label="KU-1.30 #1 acceptance band (15–200 Pa)")
    axes[1].set_xlabel("shear strain γ")
    axes[1].set_ylabel(r"$|G^{\rm 3D}_{\rm affine}(\gamma)|$ (Pa)")
    axes[1].set_title("Preliminary affine stiffening (upper bound)\n"
                      "NOT publication-grade — needs Lees-Edwards + node insertion")
    axes[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "06_preliminary_strain_stiffening.png",
                                    dpi=150)
    plt.close(fig)
    t_ss = time.perf_counter() - t0

    summary = {
        "cfl_gate_report": cfl.summary(),
        "free_diffusion": {
            "t_ns": (ts * 1e9).tolist(),
            "msd_nm2": (msds * 1e18).tolist(),
            "expected_nm2": (expected * 1e18).tolist(),
            "rel_err_at_last_t": float(abs(msds[-1] - expected[-1]) / expected[-1]),
        },
        "equipartition": {
            "n_links": len(links_eq),
            "n_time_samples": 200,
            "mean_E_over_half_kT": mean_e / target,
            "gate_passed": eq_gate.passed,
            "gate_report": eq_gate.summary(),
        },
        "preliminary_strain_stiffening": {
            "gamma": gammas.tolist(),
            "sigma_xy_2D_N_per_m": sigma_2d.tolist(),
            "G_2D_N_per_m": G_2d.tolist(),
            "G_3D_Pa": G_3d.tolist(),
            "G_0_affine_2D_N_per_m": float(G_2d[0]),
            "G_0_affine_3D_Pa_using_h_2um": float(G_3d[0]),
            "layer_thickness_for_2D_to_3D_m": h_layer,
            "publication_grade": False,
            "method": "affine deformation, no Langevin relaxation, no Lees-Edwards",
            "deferred_reason": ("PI note: revisit exact intersection-node "
                                "insertion or finer bead discretization, and "
                                "implement Lees-Edwards PBC, before quantitative "
                                "KU-1.30 #1, #2 comparison."),
        },
        "timings_s": {
            "free_diffusion": t_diff,
            "equipartition": t_eq,
            "strain_stiffening": t_ss,
        },
    }
    with open(OUT / "02_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
