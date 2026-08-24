"""Native (device) gates for the ac fluid spine — run on the A5000 by the lead (I0-A: CUDA only).

The fluid-spine kernels (``field_grid`` / ``biot_substrate`` / ``velocity`` / ``transport``) were authored
as SOURCE on the dev Mac, which cannot execute Warp-CUDA. Everything in this package is therefore UNPROVEN
until it runs on the gbook RTX A5000. Each runner builds a small fixed-seed case, drives the Warp-CUDA path
AND the pure-NumPy acceptance oracle (``fv_reference`` / ``transport_reference`` / ``darcy_analytic``) from
the IDENTICAL initial state + parameters, and reports the round-off parity residual.

NG-0 (``ng0_parity``): Warp <-> NumPy source-faithfulness on 16^3 — the first proof the device kernels ARE
the CPU-verified stencils. STRUCTURAL / parameter-agnostic (round-off parity, not a magnitude verdict).
"""
