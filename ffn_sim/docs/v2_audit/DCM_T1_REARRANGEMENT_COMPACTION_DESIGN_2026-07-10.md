# DCM (b) — T1 rearrangement-driven compaction: junction-gated interfacial-tension minimisation (DESIGN)

**Date:** 2026-07-10 · **Owner:** Lead session · **Branch:** `dcm/main` · **Status:** DESIGN (PI chose (b) full build after the S4 G4 FAIL)
**Prereq:** `DCM_CADHERIN_CLUSTER_REARRANGEMENT_DESIGN_2026-07-09.md` §6c (S4 verdict: aggregate-σ compaction is NOT
maturation-rate-limited — the driver bypasses the junction). Machinery map: this session's Explore recon.

---

## 0. Why we are here (the S4 finding that forces (b))

S1–S3 built a correct cadherin junction mechanism (load-sharing cluster → emergent lifetime; contact-age maturation →
junctions lock over τ_mature; force ∝ engaged count `m`). S4 measured the payoff at native scale and it **FAILED G4**:
the compaction RATE is not τ_mature-gated (τ=30 s ≈ τ=600 s trajectories). Root cause: the compaction DRIVER is the
**aggregate-σ Foty liquid-drop** (`aggregate_laplace_kernel`) — an external, inward **radial** Laplace pressure on the
free faces toward the global centroid. It densifies a loose aggregate by **radial gap-closing** (cells translate inward
together), which needs **no neighbour exchange and no bond breaking** → junction state never enters the rate. The first
~20 s of every S4 run is identical across all maturation settings = pure drag-limited gap-closing.

## 1. Target physics — active-foam rearrangement (T1), rate-limited by junction remodelling

Real spheroid compaction (24–48 h) proceeds by **cells rearranging** — neighbour exchange (T1 transitions) — driven by
**differential adhesion** (Steinberg DAH: cells maximise low-energy cell–cell contact, minimise high-energy cell–medium
free surface). The tissue-surface-tension σ is the *macroscopic* summary of this; its *microscopic* origin is the
cadherin adhesion differential `w_adh`. Crucially, rearrangement is **tangential** (cells slide past each other to
re-pack / fill voids), and tangential sliding requires the cadherin junction at the trailing interface to **remodel**
(break + reform). So the rearrangement rate — hence the compaction rate — is set by **junction turnover**: matured
(long-lived, high-`m`) junctions resist sliding → slow rearrangement → slow compaction.

## 2. The key realisation — the rate lever is ALREADY in the model; only the DRIVER is wrong

From the machinery map, two facts make (b) mostly a **re-driving**, not a new force:

1. **The cadherin cluster force already provides maturation-gated tangential friction.** `cadherin_bond_force_cluster_kernel`
   applies `F = m·k_trans·(L−r0)` per bond — the engaged count `m` scales the stiffness. A matured junction has high `m`
   (→ stiffer tangential spring) AND long lifetime (→ the bond persists through more of a sliding displacement before it
   turns over). So matured junctions ALREADY resist tangential sliding more, and nascent ones less. The friction-gating
   is present — it just never gets exercised, because the current driver (radial σ) produces almost no tangential sliding.
2. **A tangential, rearrangement-demanding driver already exists**: `differential_surface_tension_kernel` (DAH/foam) —
   contact faces carry reduced tension `γ_eff = γ_free − w_adh/2`, free faces keep `γ_free`. The differential drives cells
   to **grow cell–cell contact area** (Young–Dupré spreading), which is intrinsically **tangential** (a cell spreads over
   its neighbours). This is the microscopic σ, and it densifies by rearrangement, not radial translation.

**So (b) = replace/subordinate the radial `aggregate_laplace_kernel` with the tangential differential-tension driver, and
let the existing cluster+maturation cadherin friction gate the resulting rearrangement.** The compaction rate should then
track τ_mature — the G4 that σ failed. This is the fine-grained option (the driver is the actual microscopic adhesion
differential, not a lumped global radial field), and it reuses validated kernels.

## 3. Mechanism — concrete

- **Driver**: `differential_surface_tension_kernel` ON (DAH), `aggregate_laplace_kernel` OFF (or demoted to a weak
  envelope-rounding term, not the densification driver). `w_adh` = the cadherin adhesion energy (lit-anchored, KB-4.x /
  the existing `w_cs`); free-vs-contact face classification is the existing hash-grid-over-face-centroids test.
- **Gate**: cadherin `--cad-cluster --cad-mature` (S1–S3). The bond force `m·k_trans` provides tangential friction; the
  contact-age maturation raises `m`→`n_mature` and the lifetime over τ_mature → sliding friction rises over τ_mature.
- **Start**: **confluent (touching) init** (`build_confluent_cluster` / `builder=confluent`), so there is no radial gap to
  close — the ONLY way to densify/round is tangential rearrangement, which is exactly the junction-gated process. (A loose
  start reintroduces the fast radial gap-closing that masks the effect — the S4 confound.)
- **(if needed) active T1 contraction**: if DAH + friction alone does not eliminate internal voids (memory §2e: differential
  tension facets/rounds but may not globally densify), add an **active junctional contraction that shrinks high-energy
  interfaces** (drives T1), realised via `f_contract` **gated by maturation state** (a matured junction contracts its
  interface toward disappearance at a rate set by its remodelling) — new code, only if S1 shows DAH is insufficient.

## 4. Validation gates (write BEFORE implementing; do not loosen; PI-authored contracts)

- **G1 — DAH drives compaction from a confluent start**: with `diff_tension` ON + `aggregate_tension` OFF, a confluent
  aggregate densifies/rounds (porosity ↓ and/or asphericity ↓) via rearrangement — proving the tangential driver works
  without the radial σ field.
- **G2 — the rate is junction-turnover-gated (THE payoff, the S4 G4 that σ failed)**: sweep τ_mature (e.g. 30 / 300 / 600 s)
  → the compaction half-time **tracks τ_mature** (τ↑ → slower), NOT the drag/σ timescale. Equivalently, freezing junction
  turnover (very long lifetime / high n_nascent) stalls rearrangement; fast turnover speeds it.
- **G3 — it is REARRANGEMENT, not gap-closing**: measure neighbour-exchange (cadherin bond partner turnover / cell
  neighbour-set change), not just porosity — the densification must come with actual T1 events, and their rate gates it.
- **G4 — σ-independent**: the compaction does not require `aggregate_laplace_kernel` (DAH alone drives it). If σ is kept,
  it is a minor envelope-rounding term, not the densifier.
- **G5 — back-compat / physiological baseline**: all existing runs unaffected (diff_tension already default-off); w_adh,
  γ_free at lit values; no gate tuned to a compaction target.
- **G6 — NATIVE (A5000)**: N=400+ full-compartment, stable, τ_mature-gated rate.

## 5. Staging (isolate-first, measure before building new force)

- **(b)S1 — the DECISIVE isolate experiment (mostly config, ~no new force)**: confluent N=64/400, `diff_tension` ON,
  `aggregate_tension` OFF, `--cad-cluster --cad-mature`, sweep τ_mature ∈ {30, 300, 600 s}. Measure porosity(t)/asphericity(t)
  + a neighbour-exchange counter. **Decision:** if the compaction rate tracks τ_mature (G2) → the mechanism was there;
  (b) is done at the driver level → harden + native + G3/G6. If DAH does NOT densify or is NOT gated → build (b)S2.
- **(b)S2 — active maturation-gated junctional T1 contraction (only if S1 insufficient)**: `f_contract` (or a new interface-
  shrink force) gated by per-junction maturation so matured junctions drive/limit T1 at the remodelling rate. New kernel,
  isolate-validated, gated default-off.
- **(b)S3 — G2/G3 payoff at native scale (A5000)** + a neighbour-exchange metric committed.
- **(b)S4 — long run toward the min–hr regime + morphology HTML viewer** (per the visualize rule).

## 6. Risks / open questions (surface to PI)

1. **DAH may round/facet without globally densifying** (memory §2e: pairwise/differential levers didn't eliminate voids).
   If so, G1 fails → (b)S2 active T1 contraction is required. (b)S1 measures this cheaply first.
2. **A confluent start is "already compact"** (porosity low) — need a confluent-but-frustrated init (packing defects/voids to
   relax by rearrangement), else there is nothing to compact. Use a confluent init with imperfect packing (Voronoi + defects)
   or a mildly loose-then-confluent state; define in (b)S1.
3. **Tangential friction magnitude**: whether `m·k_trans` friction is strong enough to actually gate the rate (vs cells
   sliding freely regardless) is empirical — (b)S1/G2 tests it; do not tune k_trans to force gating (magic-number rule).
4. **w_adh anchor**: the adhesion differential must be the lit cadherin adhesion energy, not fitted to a compaction target.

## 7. One-line summary
Re-drive compaction with the **tangential differential-adhesion (DAH) interfacial-tension** driver instead of the radial
aggregate-σ field, from a **confluent start**, so densification proceeds by **cadherin-gated tangential rearrangement** —
the maturation-scaled bond friction (`m·k_trans`, already built) then throttles the rate to ~τ_mature, the junction-rate-
limited compaction the S4 σ-driver could not produce; add active maturation-gated T1 contraction only if DAH alone
does not eliminate voids.
