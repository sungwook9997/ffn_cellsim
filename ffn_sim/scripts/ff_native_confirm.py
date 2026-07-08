"""Native-scale (A5000) confirmation of the CPU findings: (1) the small-strain undrained modulus was an
under-relaxed transient (converges DOWN toward the soft γ-floor), (2) the drained+K_drained equilibrium, and
(3) a fast physical press point with the bulk-η-calibrated drag. Run on the gbook A5000 (device=cuda:0).

Usage: python ff_native_confirm.py <NF> <mode>   mode in {probe, converge, press}
"""
import sys, time
import numpy as np
from ffn_sim.ff.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from ffn_sim.ff.network_warp import simulate_whole_cell_compression_on_device

DEV = "cuda:0"
NF = int(sys.argv[1]) if len(sys.argv) > 1 else 38000
MODE = sys.argv[2] if len(sys.argv) > 2 else "probe"
NU = 0.5
TURGOR_EVERY = 50            # native: fewer ConvexHull(Nc) calls (the cost driver)


def _build():
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=NF, n_xl=NF, n_myo=NF // 10,
                                  rng=np.random.default_rng(1))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    return cx


def E_h(F, R0, d1):
    return (F / ((4 / 3) * np.sqrt(R0) * d1 ** 1.5)) * (1 - NU ** 2) if F > 0 and d1 > 0 else 0.0


t0 = time.time()
cx = _build()
Nc, R0 = cx.net.n_nodes, cx.R0_mean
print(f"# NATIVE build NF={NF} Nc={Nc} R0={R0:.3f}µm build={time.time()-t0:.1f}s dev={DEV} turgor_every={TURGOR_EVERY}", flush=True)
s = 0.03; d1 = R0 * s

if MODE == "probe":
    t = time.time()
    _, m = simulate_whole_cell_compression_on_device(cx, NMIIA_MINIFIL_STALL_PN, strain=s, n_steps=500,
                                                     turgor_every=TURGOR_EVERY, device=DEV)
    dt = time.time() - t
    print(f"# probe n_steps=500: {dt:.1f}s = {dt/500*1000:.1f} ms/step; F={m['F_plate_pN']:.1f} pN  "
          f"E={E_h(m['F_plate_pN'],R0,d1):.1f} Pa", flush=True)
    print(f"# ETA converge (4000+8000+16000 undrained + 8000 drained ≈ 36000 steps): {36000*dt/500/60:.0f} min", flush=True)

elif MODE == "converge":
    print("mode        n_steps   F[pN]     E[Pa]    ×249   wall[s]", flush=True)
    for ns in (4000, 8000, 16000):
        cx2 = _build(); cx2.R0_mean = R0
        t = time.time()
        _, m = simulate_whole_cell_compression_on_device(cx2, NMIIA_MINIFIL_STALL_PN, strain=s, n_steps=ns,
                                                         turgor_every=TURGOR_EVERY, device=DEV)
        E = E_h(m["F_plate_pN"], R0, d1)
        print(f"undrained  {ns:6d}  {m['F_plate_pN']:8.1f}  {E:8.1f}  {E/249:5.2f}×  {time.time()-t:6.0f}", flush=True)
    cx3 = _build(); cx3.R0_mean = R0
    t = time.time()
    _, md = simulate_whole_cell_compression_on_device(cx3, NMIIA_MINIFIL_STALL_PN, strain=s, n_steps=8000,
                                                      pressure_setpoint=40.0, K_drained_Pa=300.0,
                                                      turgor_every=TURGOR_EVERY, device=DEV)
    E = E_h(md["F_plate_pN"], R0, d1)
    print(f"drained+Kd  8000  {md['F_plate_pN']:8.1f}  {E:8.1f}  {E/249:5.2f}×  {time.time()-t:6.0f}", flush=True)

elif MODE == "long":
    # push undrained convergence further to find the NATIVE elastic equilibrium (does it reach the CPU soft
    # γ-floor, or plateau stiffer at native density?)
    print("mode        n_steps   F[pN]     E[Pa]    ×249   wall[s]", flush=True)
    for ns in (32000, 64000):
        cx2 = _build(); cx2.R0_mean = R0
        t = time.time()
        _, m = simulate_whole_cell_compression_on_device(cx2, NMIIA_MINIFIL_STALL_PN, strain=s, n_steps=ns,
                                                         turgor_every=TURGOR_EVERY, device=DEV)
        E = E_h(m["F_plate_pN"], R0, d1)
        print(f"undrained  {ns:6d}  {m['F_plate_pN']:8.1f}  {E:8.1f}  {E/249:5.2f}×  {time.time()-t:6.0f}", flush=True)

elif MODE == "press":
    print("v[µm/s] N_ramp  F[pN]     E[Pa]    ×249   eta_eff  wall[s]", flush=True)
    # equilibrium reference (relaxed)
    cxe = _build(); cxe.R0_mean = R0
    _, me = simulate_whole_cell_compression_on_device(cxe, NMIIA_MINIFIL_STALL_PN, strain=s, v_press_um_s=5.0,
                                                      dwell_steps=8000, eta_bulk_Pa_s=65.9,
                                                      turgor_every=TURGOR_EVERY, device=DEV)
    Feq = me["F_plate_pN"]
    print(f"# equilibrium(dwell) F={Feq:.1f} E={E_h(Feq,R0,d1):.1f}", flush=True)
    for v in (50.0, 5.0):
        cxp = _build(); cxp.R0_mean = R0
        t = time.time()
        _, m = simulate_whole_cell_compression_on_device(cxp, NMIIA_MINIFIL_STALL_PN, strain=s, v_press_um_s=v,
                                                         dwell_steps=0, eta_bulk_Pa_s=65.9,
                                                         turgor_every=TURGOR_EVERY, device=DEV)
        Fv = max(m["F_plate_pN"] - Feq, 0.0)
        a = max(m["contact_radius_um"], 1e-6); A = np.pi * a * a; ed = v / (2 * R0)
        eta_eff = Fv / (A * ed) if A * ed > 0 else 0.0
        E = E_h(m["F_plate_pN"], R0, d1)
        print(f"{v:6.1f}  {m['n_ramp']:6d}  {m['F_plate_pN']:8.1f}  {E:8.1f}  {E/249:5.2f}×  {eta_eff:6.1f}  {time.time()-t:6.0f}", flush=True)

print("# done", flush=True)
