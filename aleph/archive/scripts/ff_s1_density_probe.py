"""CPU filament-density (NF) convergence-tail probe for S1 (NON-AUTHORITATIVE, labelled).

Shows how the bare-cortex E_fit trends with filament count toward the native NF=70686.
Coarse cortex is a dev probe only (CLAUDE.md native+full); the native GPU run is authoritative.
"""
import json
import numpy as np
from aleph.scripts.ff_s1_sphere import s1_compression_sweep, gate_s1_hertz, LOAD_PHYSIO

rows = []
for nf in [1000, 2000, 4000, 8000, 16000]:
    sw = s1_compression_sweep(
        np.linspace(0.005, 0.03, 4), n_filaments=nf, n_steps=1200,
        device="cpu", load=LOAD_PHYSIO, with_stress=False, quiet=True,
    )
    g = gate_s1_hertz(sw)
    rows.append(dict(nfil=nf, R0_um=sw["R0_um"], E_fit_Pa=g.get("E_fit_Pa"),
                     flatness=g.get("flatness"), r2=g.get("r_squared"), verdict=g.get("verdict")))
    print(f"NF={nf:6d}  E_fit={g.get('E_fit_Pa'):.1f} Pa  flat={g.get('flatness'):.3f}  {g.get('verdict')}", flush=True)

out = "aleph/outputs/mech_hier/s1_sphere/density_probe_cpu.json"
json.dump(rows, open(out, "w"), indent=2)
print(f"PROBE DONE -> {out}")
