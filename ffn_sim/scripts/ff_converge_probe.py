"""Is the measured force a CONVERGED equilibrium, or an under-relaxed (fast-press) transient?
Run the same strain at increasing relaxation steps; if F keeps dropping, we are pressing/relaxing too fast
(the cell hasn't finished deforming) → the stiffness is a numerical transient, not the equilibrium.
"""
import numpy as np
from ffn_sim.ff.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from ffn_sim.ff.network_warp import simulate_whole_cell_compression_on_device

NU, R = 0.5, None
def E_hertz(F, R0, d1):
    return (F / ((4/3)*np.sqrt(R0)*d1**1.5)) * (1-NU**2) if F > 0 and d1 > 0 else 0.0

def fresh(nf=2000):
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=nf, n_xl=nf, n_myo=nf//10,
                                  rng=np.random.default_rng(1))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    return cx

R0 = fresh().R0_mean
print(f"# convergence test: R0={R0:.3f}µm, strain=3%, δ1={R0*0.03:.3f}µm", flush=True)
print("mode       n_steps   F[pN]    E_pt[Pa]   ×249", flush=True)
for mode, kw in [("undrained", {}),
                 ("drained40+Kdr", {"pressure_setpoint": 40.0, "K_drained_Pa": 300.0})]:
    prev = None
    for ns in (800, 1600, 3200, 6400, 12800):
        _, m = simulate_whole_cell_compression_on_device(
            fresh(), NMIIA_MINIFIL_STALL_PN, strain=0.03, n_steps=ns, device="cpu", **kw)
        F = m["F_plate_pN"]; E = E_hertz(F, R0, R0*0.03)
        drop = "" if prev is None else f"  (Δ={100*(F-prev)/prev:+.1f}%)"
        print(f"{mode:14s} {ns:6d}  {F:8.1f}  {E:8.1f}  {E/249:5.1f}×{drop}", flush=True)
        prev = F
print("# done", flush=True)
