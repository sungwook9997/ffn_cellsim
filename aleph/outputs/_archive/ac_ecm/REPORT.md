# ECM B2 transaction visualization

## Scope

This is an RTX A5000 execution of the Warp-CUDA ECM topology transaction. It visualizes
proposal planning, explicit rejection, accepted refinement with graph endpoint remap,
and accepted coarsening. It is not a collagen constitutive, traction, FA-chemistry, or
production-material validation.

## Sanity measurements

- Rejected maximum position change: `0 µm`.
- Source rest length before refine: `1.5000000000000002 µm`.
- Rest length after refine→coarsen: `1.5000000000000002 µm`.
- Active node/segment/bend populations: baseline `[90, 72, 54]`;
  refined `[91, 73, 55]`; coarsened `[90, 72, 54]`.
- Topology epochs: `[0, 0, 0, 1, 1, 2]`.
- Final endpoint population is conserved at 3.

## Figures

- `figs/ecm_b2_transaction.png` — common-axis six-state transaction audit.
- `figs/ecm_b2_transaction.gif` — animated step-through of the same actual CUDA states.
- `ecm_b2_transaction.json` — projected segment topology, endpoint coordinates, IDs,
  generations, populations, epochs, and proposal plan used by the figures.

## Held

The collagen modulus conflict (30–100 Pa versus 5–100 Pa) remains HELD. No force law or
damage threshold was introduced to make this visualization.

## Figures (native ECM constitutive milestone, 2026-07-25)

- `ecm_network.html` — native ECM Mikado Collagen-I network (n_fibers=220 / 3740 nodes, KERNEL_BOUND), interactive WebGL full-res, per-node |F| (log pN) over the device SoA; from a gbook `dump_state --ecm` (relaxed + 0.05 shear). Browser-checked (rendered, no JS errors) + visually inspected: crossing Mikado collagen fibers fill the R=7.5um sphere, |F| near the relaxed floor (constitutive relaxed~0 confirmed visually), no interpenetration/wrong-sign artifacts.
