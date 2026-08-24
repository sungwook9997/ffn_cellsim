# emergence track — REPORT (I5 detector-core)

Track: **emergence** (I5 detector-core + synthetic oracles) · branch `ac/emergence` · base `ac/new-engine`.
Status: detector-core + synthetic-config oracle suite committed + visualized (**44 host-numpy tests**, all
green on the dev Mac). The native emergence PROOF and the 4-toggle ablation harness need I1c + I3 + I4 hooks
that do not exist yet — they are **specified, not built**, in `ffn_sim/ac/emergence/INTEGRATION.md`.

## What this delivers

A **label-blind** nematic / bundle-condensation detector (`ffn_sim/ac/emergence/`) coded against the FF
fiber-array geometry contract (`pos` / `fiber_offsets` / `n_fibers` — predates I4). It reuses the validated
Q-tensor order parameter `S = λ_max((3⟨n⊗n⟩−I)/2)` (bit-parity with
`ff/architecture_metrics.parallel_order_parameter`) and turns the **vacuous `bundle_count = n_fibers`** into
a real, LOCALIZING condensation measure with a **pre-registered effect-size rule** against the finite-N
isotropic null. No function reads a region / type / construction label (the anti-coupling firewall,
structurally enforced by `test_detector_contract.py`).

Modules: `nematic` (validated S) · `null_model` (finite-N null band, analytic `E[Σλ²]=3/(2N)` anchor + the
5σ effect-size rule) · `condensation` (local nematic+density fields, KD-tree spatial clustering, bundle
localization) · `detector` (`EmergenceDetector.detect()` structured report) · `synthetic` (oracle fixtures).

## Figures

Regenerate all: `python ffn_sim/scripts/ac_emergence_vis.py` (pure numpy/scipy/matplotlib on the synthetic
oracle configs — no Warp/CUDA).

- **`figs/i5_null_vs_aligned.png`** — the two nematic corners. LEFT: the finite-N isotropic null band μ±σ
  shrinking as `1/√N` (per-realisation isotropic draws overlaid, landing in the band), with the aligned
  bundle `S=1` line above — so "isotropic → S~0" is a quantitative band, not a hand-wave. RIGHT: the analytic
  cross-check — Monte-Carlo `⟨Σλ²⟩` sits exactly on the closed form `3/(2N)` (log-log). The null is anchored
  to a derivation, not taken on faith.
- **`figs/i5_localization.png`** — a planted bundle CONDENSING out of an isotropic background. LEFT: the
  label-blind local-order field `S_local` — the dense aligned patch lights up against the dark isotropic
  background. RIGHT: the recovered members (red rings) overlap the planted-truth bundle, and the recovered
  centre (★) sits on the planted centre (✕): recall=1.00, precision=0.87, centre error 0.22 µm.
- **`figs/i5_invariance.png`** — the firewall property. LEFT: `|S − S_ref|` under 40 random rotations + 40
  relabelings stays at ~1e-15, far under the 1e-12 machine-precision floor (S is rotation- and
  label-invariant). RIGHT: the recovered bundle centre is rotation-EQUIVARIANT (centre error ~0 across
  rotations).
- **`figs/i5_effect_size.png`** — the pre-registered falsifiability rule. LEFT: effect size `z` vs aligned
  fraction crosses the pre-registered `Z_crit=5` (5σ) bar as order rises (dose-response). RIGHT: specificity —
  60 independent isotropic draws all score `z<Z_crit` (no false emergence) while aligned draws sit at z≈62.

All figures overlay the closed-form null / oracle / pre-registered bar on the measured detector output,
annotate axes (S, z, fractions are dimensionless `[–]`; positions in µm), do not truncate axes, and show
ensembles as per-realisation markers + mean/band (the visualization-integrity rules).

## Gates owned (CPU, host numpy — all green)

`tests/ac/emergence/` (44 tests): isotropic → S in the null band (S~0); aligned → S=1 exactly (to the
validated formula's 1e-12 normalization floor); **antiparallel bundle still S~1** (nematic, not polar — the
SF graded-polarity regime); rotational + permutation/label invariance of S; monotonicity in aligned fraction;
MC null `⟨Σλ²⟩` == analytic `3/(2N)`; `1/√N` null-bias scaling; effect-size specificity (isotropic z<5) +
sensitivity (aligned z≫5); planted-bundle localization (recall/precision/centre); localization
rotation-equivariance + permutation-invariance; NEIGHBORHOOD_SCALE k-plateau grid-invariance; crowding ≠ order
(dense isotropic patch NOT flagged); label-blindness of the whole public API (no region/type argument).

## Deferred → native (spec in `INTEGRATION.md`, not built this session)

The 4 ablation-control toggle hooks (myosin-off / unanchored / KMC-off / angle-off) and the native emergence
PROOF (needs live I1c monomer field + I3 head-resolved NMII + I4 unified weave). **Falsifiability contract
(hard-truth #1 / §5):** a flat native result → FINDING to PI + a SEEDED-labeled scaffold fallback, NEVER
re-tuned to force condensation.

## Native-gate viz spec (lead runs on gbook)

Interactive 3-D cell-morphology **HTML** overlaying, on the LIVE cell geometry, the per-fiber `S_local`
condensation field + the recovered bundle members/centres/axes (colour by `S_local`, ring the recovered
members) — real geometry, peel/slab/cut views, full-res (no downsampling), browser-verified
(`browser_check.py`). Extend `ac_emergence_vis.py` with the ablation per-realisation traces (S/bundle
persistence rises vs all 4 controls) as those hooks land.
