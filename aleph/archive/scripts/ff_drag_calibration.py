"""Calibrate + validate the bulk-η per-node drag (γ_node=6πηR/Nc) against the analytic Newtonian ground truth
σ_visc = η·ε̇, then read Zbiral 5µm/s with the calibrated drag. η_eff = F_visc/(A·ε̇); should ≈ 65.9 Pa·s.
"""
import numpy as np
from aleph.laws.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from aleph.laws.network_warp import simulate_whole_cell_compression_on_device

NU, ETA = 0.5, 65.9
def E_hertz(F, R0, d1):
    return (F / ((4/3)*np.sqrt(R0)*d1**1.5)) * (1-NU**2) if F > 0 and d1 > 0 else 0.0
def fresh(nf=1500):
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=nf, n_xl=nf, n_myo=nf//10, rng=np.random.default_rng(1))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    return cx

NF = 1500
R0 = fresh(NF).R0_mean; s = 0.03; d1 = R0*s
print(f"# drag calibration | R0={R0:.3f}µm strain=3% Nc≈{NF} eta_target={ETA} Pa·s", flush=True)

# equilibrium (fully relaxed) reference force
_, m_eq = simulate_whole_cell_compression_on_device(fresh(NF), NMIIA_MINIFIL_STALL_PN, strain=s,
                                                    v_press_um_s=5.0, dwell_steps=20000, eta_bulk_Pa_s=ETA, device="cpu")
F_eq = m_eq["F_plate_pN"]
print(f"# equilibrium F={F_eq:.1f} pN  E={E_hertz(F_eq,R0,d1):.1f} Pa", flush=True)
print("v[µm/s]  N_ramp  t_ramp[s]  F[pN]     F_visc[pN]  a_contact  eps_dot[/s]  eta_eff[Pa·s]  E[Pa]  ×249", flush=True)
for v in (50.0, 5.0, 0.5):
    _, m = simulate_whole_cell_compression_on_device(fresh(NF), NMIIA_MINIFIL_STALL_PN, strain=s,
                                                     v_press_um_s=v, dwell_steps=0, eta_bulk_Pa_s=ETA, device="cpu")
    F = m["F_plate_pN"]; Fv = max(F - F_eq, 0.0)
    a = max(m["contact_radius_um"], 1e-6); A = np.pi*a*a
    eps_dot = v / (2.0*R0)                                        # nominal compression strain rate
    eta_eff = Fv / (A*eps_dot) if (A*eps_dot) > 0 else 0.0
    print(f"{v:7.3g}  {m['n_ramp']:6d}  {m['t_ramp_s']:8.4f}  {F:8.1f}  {Fv:9.1f}  {a:8.3f}  {eps_dot:9.3f}  "
          f"{eta_eff:11.1f}  {E_hertz(F,R0,d1):7.0f}  {E_hertz(F,R0,d1)/249:5.1f}×", flush=True)
print("# done", flush=True)
