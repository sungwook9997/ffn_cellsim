"""Direct measurement of lateral bulge (R_eq) under parallel-plate compression, undrained vs drained.
Tests PI's hypothesis: does the cortex stretch sideways (bulge) enough to conserve volume, or does it
under-expand and lose volume (-> osmotic dP spike)? Reports the sim's OWN R_eq (not inferred).
"""
import numpy as np
from ffn_sim.ff.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from ffn_sim.ff.network_warp import simulate_whole_cell_compression_on_device

DEV = "cpu"
NF = 2000
NSTEPS = 1600

def fresh():
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=NF, n_xl=NF, n_myo=NF // 10,
                                  rng=np.random.default_rng(1))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    return cx

R0 = fresh().R0_mean
print(f"# cortex N_fil={NF}, R0_mean={R0:.3f} um, n_steps={NSTEPS}, device={DEV}", flush=True)
print(f"# constant-volume oblate prediction: R_eq_cv = sqrt(R0^3 / (R0*(1-strain)))", flush=True)
print("mode        strain  half_gap  R_eq    R_eq/R0   R_eq_cv/R0  V/V0    F_plate[nN]  dP[Pa]  gApp[mN/m]", flush=True)

for strain in (0.06, 0.12, 0.18, 0.27):
    cx = fresh()
    _, m = simulate_whole_cell_compression_on_device(
        cx, NMIIA_MINIFIL_STALL_PN, strain=strain, n_steps=NSTEPS, device=DEV)
    c = R0 * (1 - strain)
    r_cv = np.sqrt(R0 ** 3 / c) / R0
    print(f"undrained   {strain*100:4.0f}%  {m['half_gap']:7.3f}  {m['R_eq_um']:6.3f}  "
          f"{m['R_eq_um']/R0:7.4f}  {r_cv:9.4f}   {m['V_over_V0']:6.4f}  {m['F_plate_pN']/1000:9.2f}  "
          f"{m['dP_turgor_Pa']:7.0f}  {m['gamma_apparent_mN_m']:7.2f}", flush=True)

# drained limit: hold dP at the resting turgor (perfect water flux) -> tension-governed Fischer-Friedrich
for strain in (0.12, 0.27):
    cx = fresh()
    _, m = simulate_whole_cell_compression_on_device(
        cx, NMIIA_MINIFIL_STALL_PN, strain=strain, n_steps=NSTEPS, pressure_setpoint=40.0, device=DEV)
    c = R0 * (1 - strain)
    r_cv = np.sqrt(R0 ** 3 / c) / R0
    print(f"drained@40  {strain*100:4.0f}%  {m['half_gap']:7.3f}  {m['R_eq_um']:6.3f}  "
          f"{m['R_eq_um']/R0:7.4f}  {r_cv:9.4f}   {m['V_over_V0']:6.4f}  {m['F_plate_pN']/1000:9.2f}  "
          f"{m['dP_turgor_Pa']:7.0f}  {m['gamma_apparent_mN_m']:7.2f}", flush=True)
print("# done", flush=True)
