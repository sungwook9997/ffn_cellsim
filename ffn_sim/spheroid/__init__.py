"""Layer-2 multicellular tumor-spheroid model (center-based: 1 particle per cell).

A lightweight, broad-scale complement to the fine-grained single-cell line, built on the
SAME HOOMD-blue + BAOAB stack. Reaches the collective / multi-cell terms of the PI
experiment's spreading law ``A/A0 = a + b/R + c/R^2`` (cohesion c, A/A0 itself) that the
single-cell line cannot reach in wall-time. Isolated as code, coupled as physics to the
single-cell line via the scale-bridge (per-cell cortical tension / traction / adhesion
energy). See ``docs/LAYER2_MULTICELL_DESIGN.md``.

Modules:
- ``observables``: pure point-cloud measurement protocol (spread area, gyration radius,
  radial density profile, nearest-neighbour stats, detached fraction). No HOOMD, no
  physics parameters — written before the physics run per the Sanity-Gate rule.
"""
