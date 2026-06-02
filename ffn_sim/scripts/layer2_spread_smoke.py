"""L2.2 spreading-mechanism smoke: verify f_active=0 reduces to G1, and f_active>0 spreads.

Mechanism check only (reduced traction; quantitative A/A0 awaits the Omidvar-anchored nN
cohesion). Boundary: f_active=0 => A/A0 ~ 1 (stable, the G1 limit). Sign-sense: f_active>0
=> A/A0 > 1 (spreads) while staying connected (low detached fraction).
"""

from __future__ import annotations

from pathlib import Path

import yaml

from ffn_sim.spheroid.params import resolve_layer2
from ffn_sim.spheroid.spreading import run_spreading

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


def main() -> int:
    resolved = resolve_layer2(yaml.safe_load(_CFG.read_text()))
    # cohesive force scale ~ Morse max |F| = D_e*alpha/2
    f_coh = resolved.D_e * resolved.morse_alpha / 2.0
    print(f"[scale] Morse max |F| per bond = {f_coh*1e12:.2f} pN  (cohesion scale)")

    common = dict(n_cells=120, settle_steps=5_000, spread_steps=15_000, n_samples=4)

    # Sweep f_active across the cohesive-escape threshold (~2.5 pN here). Mechanism check:
    # f=0 reduces to G1 (A/A0~1); above the escape threshold A/A0 rises (spreads) and the
    # detached fraction rises (transition toward unbinding). The spread-vs-detach BALANCE
    # is what the real (Omidvar nN) cohesion sets — here we only confirm the machinery.
    # nN-scale traction (now that cohesion is nN-anchored, Iturri 2020). Sweep below the
    # measured detachment force so we see clean active-wetting (spread, stay connected).
    forces_nN = [0.0, 2.0, 4.0, 6.0]
    rows = []
    for fn in forces_nN:
        r = run_spreading(resolved, f_active=fn * 1e-9, **common)
        rows.append((fn, r["area_over_a0"][-1], r["detached_fraction_final"]))
        print(f"  f_active={fn:5.1f} nN  ->  A/A0={r['area_over_a0'][-1]:.3f}  "
              f"detached={r['detached_fraction_final']:.4f}")

    base = rows[0][1]
    ok_boundary = 0.9 <= base <= 1.15
    ok_monotone = rows[-1][1] > base + 0.05  # traction below cohesion spreads vs f=0
    ok_connected = rows[-1][2] < 0.10  # ...while staying connected (clean active-wetting)
    print(f"\n[boundary f=0 -> A/A0~1 ] {'PASS' if ok_boundary else 'FAIL'}  (A/A0={base:.3f})")
    print(f"[active-wetting: spreads ] {'PASS' if ok_monotone else 'FAIL'}  "
          f"(A/A0={rows[-1][1]:.3f} @ {forces_nN[-1]:.0f} nN)")
    print(f"[stays connected         ] {'PASS' if ok_connected else 'FAIL'}  "
          f"(detached={rows[-1][2]:.3f})")
    return 0 if (ok_boundary and ok_monotone and ok_connected) else 1


if __name__ == "__main__":
    raise SystemExit(main())
