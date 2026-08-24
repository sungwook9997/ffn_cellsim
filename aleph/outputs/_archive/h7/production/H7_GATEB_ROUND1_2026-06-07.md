# H.7 GATE-B + A/A0 — first full-scale GPU production (round 1, 2026-06-07)

gbook RTX A5000, full ×40 mesoscopic scale (1000 effective filaments). GATE-B ~6 min,
A/A0 (both geometries) ~3.5 min. **PROVISIONAL — no KU-3.5 conclusion (GATE-B contract).**

## GATE-B emergent cortical tension (3 channels, reported SEPARATELY)
Operating point: FA-adhered (2026 clutches), turgor Π₀ = 133 Pa, full compartment stack.

| channel | value (mN/m) | reading |
|---|---|---|
| `γ_active` (soft MOP, cortical-only) | **0.000** | active cortical tension FLOORS at zero — 10538 cortical bonds counted, 2026 adhesion excluded |
| `γ_rigid` (M-SHAKE Lagrange) | **0.391** | in band, but likely turgor-transmitted structural backbone tension (+ possible M-SHAKE shunt) |
| `γ_passive` (turgor Young-Laplace) | **0.499** | **CIRCULAR** — identity with the B3 turgor setpoint (133 Pa = 2·0.50/R); NOT a measurement |
| `γ_structural` (active+rigid) | 0.391 | turgor NEVER folded in (B3/GATE-B contract) |
| band overlay | [0.35, 0.65] | rounded/de-adhered non-MCF7 proxy — gate-contract review pending PI |

**Honest read (provisional):** in the full cell the in-band total is carried by **structural
(rigid backbone, ~turgor pre-tension) + passive turgor**, while **active myosin tension still
floors near zero** — consistent with the project's γ-floor reframe (active = generation-limited).
Caveat (direction-review #4): µs-scale warmup (2000 CFL steps ≈ ms physical) cannot reach the
seconds-scale myosin/turnover timescale, so `γ_active`≈0 is the GENERATION floor, not a settled
value. Turnover was OFF. `γ_passive` is circular with the chosen turgor; `γ_rigid` needs an
unconstrained-soft cross-check to bound the M-SHAKE shunt.

## A/A0 spreading (both single-cell lamellipodium geometries, membrane OFF, single seed)
| geometry | A0 (µm²) | A/A0 final | polarization | reading |
|---|---|---|---|---|
| basal_ring (isotropic) | 17.47 | **1.140** | 0.23 | footprint grows ~radially-symmetric (+14%) |
| polarized_patch (migrating) | 25.06 | **1.008** | 0.90 | extends directionally; less hull-AREA gain (+0.8%) |

Fig: `outputs/h7/figs/h7_spreading_compare_full.png`. Provisional: single seed, short run,
membrane OFF (capping uniform; the +ŷ CappingUpdater normal is dormant here, step-2 fix flagged).

## Next experiments (from direction-review + novelty-analysis)
1. **Probe γ_active vs run-length / physical time** (does it rise, or is it truly floored?) +
   report physical time, not steps.
2. **5-seed ensemble** for γ_rigid reproducibility + A/A0.
3. **Turnover ON** (sets the tension-vs-length maximum; Chugh 2017) — and the novelty-defining
   experiment: reproduce Chugh's non-monotonic tension-vs-filament-length maximum EMERGENTLY, untuned.
4. **Cross-check γ_rigid** vs an unconstrained soft-backbone run (M-SHAKE shunt bound).
5. **De-circularize γ_passive** (independent osmotic datum) OR keep labeling it a setpoint.
6. **Gate-contract**: the [0.35,0.65] band is a non-MCF7 proxy — surface to PI for the right target.
