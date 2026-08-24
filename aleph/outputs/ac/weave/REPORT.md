# ac.weave (I4) — unified emergent network: analytic-gate report

Track **unified-weave** (`ac/weave`), increment **I4** — the CORE of pillar P3 ("one unified emergent
network"). Host-numpy analytic oracles are the local (Mac, CPU) gate; the `*_warp` modules are CUDA-only
kernel SOURCE (authored, never launched here). Native gates (protrusion / retrograde-flow / bundle
condensation on the full 70,686-filament cell) are the lead's, on the gbook A5000 — see `INTEGRATION.md` §3.

## What I4 delivers

`weave_cell()` + `WovenCell` assemble ONE actomyosin network from a region list, with **first-class
cross-region bonds** (dorsal-SF -> arc, cap -> LINC -> nucleus). Regions are SEEDED, not scripted: filament
orientations are isotropic-conditional on each region's declared manifold + nucleator + polarity; the
SF/cap/arc/filopodium bundles are meant to CONDENSE (the falsifiable I5 native proof, NOT claimed here). The
pre-made `ff.weave._build_bundle` sarcomeric bundle is DEMOTED to a non-authoritative control. `stress_fiber.py`
then WIRES the SF into the whole-cell loop (the biggest structural gap, G1/G2): each SF is discovered label-blind
and produces an emergent axial-prestress load path + a traction dipole (FA end -> substrate, LINC end -> nucleus),
with the SF magnitude held as an I0-B3/I0-B6 GAP.

Fidelity (honest):
- **FULL fidelity** — cortex (bit-identical delegation), lamellipodium dendritic Arp2/3 REBUILD (angle-harmonic
  branch, thermal sigma_theta, N-fixed dormant activation, flux-limited nucleation), topology-reforming
  crosslink KMC, the `walk_dir` power-stroke hand-off.
- **Honest region builders / SEEDS** — ventral SF, dorsal SF, transverse arc, perinuclear cap (LINC stubs),
  filopodium: correct manifold + isotropic-conditional seed + alpha-actinin/NMII + FA/LINC anchor sites, concat
  correctly; their bundle condensation is the I5 native proof.
- **SF wiring (`stress_fiber.py`) — the biggest structural gap CLOSED at the representation level.** Stress
  fibers were seeded but not IN the whole-cell loop (FILAMENT_SUBSYSTEM_PLAN G1/G2). `wire_stress_fibers(cell)`
  now discovers each SF LABEL-BLIND (a motorized + FA/LINC-anchored connected bundle — the `region_id` firewall
  holds) and computes its EMERGENT axial prestress from the discrete NMII head reactions accumulated along the
  fiber backbone (NOT a lumped `k_SF`), plus the traction DIPOLE it delivers — the FA end into the substrate loop
  (G1/G2), the LINC end into the I7 nucleus loop. The emergent taxonomy (ventral/dorsal SF vs cap SF vs
  non-contractile filopodium) falls out of anchor-kind + contractility (INV-4/INV-5), never a label. SF MAGNITUDE
  is a GAP (report-not-tune).
- **Deferred** — the myosin-turnover -> local-density SF condensation TRIGGER (build-plan 6.1) and full
  growing-N nucleation stay I5/I9+; the cap -> nucleus tether's nucleus-side endpoints are wired at I7
  (ac.nucleus).

## Gates (all CPU-green, 53 tests)

| Gate | Result |
|---|---|
| **Gate-1** `weave_cell([CORTEX])` all-OFF -> bit-identical to `ff.weave.weave(CORTEX)` | **PASS** — every array (pos/xl/myo) max\|diff\| = 0; identical gamma by construction. Delegation is population-independent -> CPU gate runs a small cortex; full-70,686 parity is the native gate. |
| Branch-angle distribution matches theta0 +- sigma_theta (3A.b) | **PASS** — measured mean 70.4 deg, SD 8.9 deg vs theta0=70 deg, sigma_theta=9 deg (Faessler 2020); the angle FLUCTUATES (not a rigid 72 deg). |
| Angle-harmonic FD-gradient sign arbiter + Newton 3rd law | **PASS** — analytic force = -dU/dtheta (FD), sum force = 0, restores toward theta0. |
| Equipartition sigma_theta = sqrt(kT/k_theta) round-trip | **PASS** — the ARP23 anchor (0.173 pN.um/rad^2 <-> sigma=9 deg @310 K) inverts. |
| Cross-region bond concat correctness | **PASS** — declared bonds resolve to global endpoints in the two named regions, within reach. |
| Single-crosslinker-channel double-count guard | **PASS** — both channels (kmc_reattach + r0_creep) -> raise; I4 selects `kmc_reattach`. |
| Crosslink KMC detailed balance + topology reform + conservation | **PASS** — bound fraction -> k_on/(k_on+k_off); partners hop (~80% re-partnered) while count is conserved. |
| Flux-limited nucleation monotone in monomer + NPF; N-fixed activation | **PASS** — rich monomer activates more dormant daughters; c=0 or npf=0 -> no nucleation. |
| walk_dir from barbed-end polarity; passive default; bipolar contraction sign | **PASS** — bound head walk_dir -> barbed end; free head zero (passive); anti-parallel bipolar -> contractile. |
| **SF prestress sign-sense** (antiparallel FA<->FA -> contractile closed dipole; residual->0) | **PASS** — organized sarcomere: uniform positive axial prestress, `closed_dipole=True`; a one-end (dorsal) SF is open (residual = dipole). |
| **SF wiring firewall** (invariant to a `region_id` permutation) | **PASS** — `wire_stress_fibers` discovers SFs label-blind (motorized + anchored bundle); scrambling `region_id` leaves the ledger + dipoles bit-identical. |
| **SF ledger partition** (unique IDs; cortex excluded, no double-count) | **PASS** — bundles partition the anchored-motorized set (`N_unique == N_summed`); the un-anchored cortex is disjoint from every SF bundle. |
| **SF magnitude GAP guard** (no `f_head` -> `tension=None`, GAP status) | **PASS** — the magnitude-independent SHAPE is populated; the magnitude stays null until I0-B3/I0-B6 close (never tuned). |

## Figures

- `figs/i4_branch_angle.png` — dendritic Arp2/3 branch-angle distribution (measured, N=4000) OVERLAID with the
  theta0 +- sigma_theta Boltzmann oracle (prop sin(theta) exp(-U/kT)); theta0=70 deg line + sigma_theta band.
  Measured mean 70.4 deg / SD 8.9 deg matches — the branch angle FLUCTUATES thermally, NOT a rigid 72 deg
  (CLAUDE.md worked-example).
- `figs/i4_cortex_parity.png` — Gate-1: cortex node cloud (weave_cell vs ff.weave, perfectly overlapping) +
  per-array max\|diff\| bars (all 0e+00) — the bit-identical regression anchor -> identical gamma.
- `figs/i4_dendritic_net.png` — the rebuilt lamellipodium: mother filaments (red) + Arp2/3 daughters (blue)
  branching at ~70 deg junctions (black) — a genuine dendritic mother/daughter TREE, not the ratchet proxy.
- `figs/i4_crosslink_kmc.png` — topology-reforming crosslink KMC: bound fraction tracks
  k_on/(k_on+k_off)=0.714 (detailed balance) while ~80% of bonds re-partner over time (the topology reforms
  -> bundles condense at I5).
- `figs/i4_stress_fiber.png` — the SF prestress load path (magnitude-independent SHAPE). LEFT: an organized
  antiparallel FA<->FA sarcomere carries ~uniform positive axial prestress (condensed SF) while the isotropic-
  conditional ventral-SF SEED (S=0.14) is disordered / mixed-sign (uncondensed). RIGHT: the organized bundle is a
  CLOSED contractile dipole (balance residual -> 0); the seed is OPEN (residual ~ dipole) — the closed dipole is
  what CONDENSES at I5, not the seed. Caption states MAGNITUDE IS GAP (I0-B3 f_stall + I0-B6 head density);
  shape-only, report-not-tune.

Regenerate: `PYTHONPATH=. python ffn_sim/scripts/ac_weave_vis.py`.

## I0-B4 parameter GAPs (surfaced to PI — `params_i0b4.yaml`)

GROUNDED: theta0=70 deg, sigma_theta=9 deg (Faessler 2020), k_theta=0.173 pN.um/rad^2 (DERIVED by
equipartition, not tuned). GAP — PI (native protrusion/bundle gate INVALID until closed): Arp2/3
branch-nucleation rate `k_arp0`, monomer half-saturation `K_m`, capping rate `k_cap`, NPF areal density;
filopodium fascin bundle stiffness + inter-filament spacing + bundling angle + filaments-per-bundle +
tip-nucleator density. Consumed read-only: I1c monomer field (D_c/c0), I3 hand API (k_xb/f_stall/v0/kappa at
I0-B3).

**SF magnitude GAP (`stress_fiber.py`, report-not-tune).** The single-SF force scale (~5-6 nN) is NOT closable
from the coarse seed: it is `f_head` (I0-B3 per-head NMII stall) × the engaged-head DENSITY per cross-section
(I0-B6 minifilaments-per-bundle) — and the density is geometrically not back-solvable to hit a 10-30 nN network
band (build-plan §4C, [[project-sf-nmii-forcescale-result]]). So `wire_stress_fibers` defaults to `f_head_pN=None`
-> the magnitude-independent SHAPE only (`tension=None`, GAP status). The native SF-traction magnitude gate is
INVALID until PI closes f_stall + head density; the load-path SHAPE + firewall + ledger gates hold regardless.

## Honest scope / deepest risk

Emergence may NOT condense bundles in an overdamped quasi-static solve (build-plan 6.1): the
myosin-turnover -> local-density SF trigger is still absent, so the emergence claim is genuinely falsifiable —
it is I5's native proof, not this increment's claim. I4 delivers the buildable network + the regression + the
branch-angle gate. No I0-B magnitude was tuned to hit any outcome.
