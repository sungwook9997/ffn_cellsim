# Active Cell ERM Bell on/off kinetics

Status: **mechanism and CUDA transaction gates PASS; MCF7 production values remain blocked.**

## Outcome

The membrane–cortex ERM primitive now has two explicitly separated modes:

- the old membrane-energy hard threshold is retained as a diagnostic compatibility path only;
- the production-form path carries the full unilateral tensile spring force and advances Bell slip
  detachment plus capture-gated rebinding exactly once per accepted outer physical step.

The Bell path never mutates bond state during inner mechanical iterations. A rejected outer transaction leaves
`bound`, `rest`, and event counters bit-exact. A successful rebind sets `rest` to the current formation length,
so it cannot inject an unrecorded prestress. Detach/attach events and the post-step bound count remain
device-resident during the physical loop and enter the ledger only during post-step diagnostics.

## Production firewall

There are deliberately no defaults for active MCF7 ERM density, `k_on`, `k_off0`, Bell force `F0`, or capture
radius. Supplying only part of the kinetic contract fails before CUDA initialization; supplying values without
a provenance label also fails. The `erm_density_mcf7_production` latch additionally requires both an explicit
density and the complete Bell contract.

NG-3 independently checks the serialized ledger for:

- `erm_kinetics_mode == BELL_SLIP_ON_OFF`;
- nonempty kinetics provenance;
- positive `k_on`, `k_off0`, `F0`, and capture radius;
- `formation_length` rebind policy;
- a source-grounded single-ERM preload-force basis.

The final item is still absent. The historical 11.4347 pN value is derived from continuum membrane
tube-extraction energetics, not a measured single-ezrin rupture force. Consequently the old `13.102 Pa` and
`2,297 ERM` arithmetic is quarantined as diagnostic and cannot close NG-3.

## Verification

- Local full `ffn_sim/tests/ac/`: **337 PASS; 27 CUDA-only tests skipped as intended** (364 collected).
- RTX A5000 focused ERM/compartment/wiring suite: 46/46 PASS.
- RTX A5000 full `ffn_sim/tests/ac/`, split by package to avoid one-process CUDA module accumulation, plus
  the final public-solver wiring regression: **364/364 PASS**.
- Ruff on all changed Python files: PASS.
- CUDA Bell gate: full tensile load preserves Newton closure; rejected KMC is bit-exact; accepted high-hazard
  control detaches one bound tether, rebinds one in-capture tether at formation length, and leaves one
  out-of-capture tether free.

No physiological outcome or bleb timing is claimed from the test-only high-hazard kernel control.

## Figures

- `figs/fig_ac_erm_bell_contract.png` — dimensionless Bell off-rate and Poisson transition probability, with
  the unresolved MCF7 parameter and single-linker-capacity contracts named explicitly.

## Open items

1. Register and PI-ratify MCF7 active membrane–cortex ERM surface density.
2. Register a complete ERM `k_on`, `k_off0`, `F0` (or transition-state distance), and capture-radius contract.
3. Register the rebind formation-length policy or ratify a different mechanistic policy before production.
4. Replace the continuum membrane-tether force basis with a source-grounded single-ERM preload/lifetime
   contract, then rerun the full 70,686-filament physiological baseline.
