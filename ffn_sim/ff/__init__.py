"""ff — Filament-FEM engine layer (Cytosim physics, NVIDIA Warp).

The going-forward fine-grained engine: a discrete cross-linked filament network
solved as a mechanical (FEM-like) system — energy assembly + implicit nonlinear
solve + autodiff — rather than explicit thermostatted MD. The reason MD is not
needed here: thermal effects that enter through the *mean* (entropic WLC
elasticity) become a deterministic constitutive term, and active/stochastic
kinetics (myosin, catch/slip bonds) enter as an event-driven kinetic layer on
top of the mechanical solve; only pure-thermal rare-event configurational search
would require sampling, and an active cell has little of that.

Physics reference = **Cytosim** (Nédélec & Foethke 2007, *New J. Phys.* 9:427)
— fibers as bending beams + implicit Brownian step + the Hand binding model —
exactly analogous to how the ``dcm/`` layer references SimuCell3D. Cytosim itself
(runnable external C++) is the independent parity oracle, so the active layer is
never validated against only itself.

Reuses ``dcm/`` machinery: discrete network forces (``dcm.network_warp``), the
matrix-free linearly-implicit / Newton solver (``dcm.dcm_warp_implicit``), and
the end-to-end differentiable loop (``dcm.dcm_warp_diff``).

Status: scaffold (2026-06-29). Build tracked as Stage 6 — first milestone is the
γ-floor prototype (WLC bending core + myosin active prestress, no explicit BAOAB),
the decisive "MD-free reproduces the same γ?" experiment.
"""
