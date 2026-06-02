"""One entry-point that regenerates ALL Layer-2 figures (project visualize convention).

Mirrors the single-cell ``h*_vis.py`` pattern: one script refreshes every Layer-2 figure
into ``outputs/layer2/figs/``. Extend as new Layer-2 units land (L2.3 A/A₀ sweep, L2.5
cadherin, L2.6 invasion). Per-driver auto-viz still applies (each smoke/driver bakes its own
figure at run end); this is the bulk-refresh entry-point.

Usage:  python -m ffn_sim.scripts.layer2_vis
"""

from __future__ import annotations

import sys

from ffn_sim.scripts import layer2_aa0_sweep, layer2_g1_smoke, layer2_spread_smoke


def main() -> int:
    print("=== [layer2_vis] G1 stable-aggregate figure ===")
    layer2_g1_smoke.main(["--n-cells", "200", "--settle", "20000", "--measure", "10000"])
    print("\n=== [layer2_vis] L2.2 motility-mechanism figure ===")
    layer2_spread_smoke.main([])
    print("\n=== [layer2_vis] L2.3 A/A0(R0) law figure (slow: R0 sweep) ===")
    layer2_aa0_sweep.main([])
    print("\n[layer2_vis] done — figures in ffn_sim/outputs/layer2/figs/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
