"""Physical press-speed experiment: with η=65.9 Pa·s cytoplasm drag, does pressing faster than the cell can
deform read a NON-equilibrium transient? Read force at ramp-end (dwell=0 = the press-speed reading) vs after a
long dwell (relaxed → equilibrium). Undrained/mechanical-only (no Lp, no turnover) to isolate the mechanics.
"""
import numpy as np
from ffn_sim.ff.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from ffn_sim.ff.network_warp import simulate_whole_cell_compression_on_device

NU = 0.5
def E_hertz(F, R0, d1):
    return (F / ((4/3)*np.sqrt(R0)*d1**1.5)) * (1-NU**2) if F > 0 and d1 > 0 else 0.0
def fresh(nf=2000):
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=nf, n_xl=nf, n_myo=nf//10, rng=np.random.default_rng(1))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    return cx

R0 = fresh().R0_mean; s = 0.03; d1 = R0*s
print(f"# physical press speed | R0={R0:.3f}µm strain=3% δ1={d1:.3f}µm (undrained, mechanical-only)", flush=True)
print("v_press[µm/s] dwell  N_ramp  t_ramp[s]  dt_real[s]   F[pN]    E[Pa]   ×249  read", flush=True)
for v, dwell, tag in [(5.0, 0, "ramp-end"), (5.0, 12000, "+dwell(eq)"),
                      (1.0, 0, "ramp-end"), (0.2, 0, "ramp-end"), (0.05, 0, "ramp-end")]:
    _, m = simulate_whole_cell_compression_on_device(
        fresh(), NMIIA_MINIFIL_STALL_PN, strain=s, v_press_um_s=v, dwell_steps=dwell, device="cpu")
    F = m["F_plate_pN"]; E = E_hertz(F, R0, d1)
    print(f"{v:11.3g}  {dwell:5d}  {m['n_ramp']:6d}  {m['t_ramp_s']:8.3f}  {m['dt_real_s']:.2e}  "
          f"{F:8.1f}  {E:7.0f}  {E/249:5.1f}×  {tag}", flush=True)
print("# done", flush=True)
