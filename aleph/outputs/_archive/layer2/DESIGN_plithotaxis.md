# Layer-2 §E DESIGN — coherent collective traction (SPP plithotaxis)

**Status:** DESIGN (PI-ratified mechanism 2026-06-04: "Full SPP plithotaxis" + "deep-research ξ·τ → SE first").
Numeric anchors marked ⏳ PENDING-RESEARCH are filled from the in-flight deep-research pull
(velocity correlation length ξ, persistence time τ, plithotaxis stress magnitudes, front velocity,
CIL, SPP↔ξ/τ mapping) before any production run. **No value is fitted to the PI A/A₀ (overlay-only).**

---

## 1. The diagnosis this mechanism answers (why it is needed, data-backed)

§B2/§C/§D eliminated scale, measurement-definition, and cohesion as the ~5–9× magnitude-gap suspect.
The decisive piece is the `daxis_yield` table (turnover/ductile cohesion, N₀=4000, R₀≈151 µm, 2 seeds):

| basal-crawl f | A/A₀ core | A/A₀ raw | eject |
|---|---|---|---|
| 0 | 1.35 / 1.30 | 1.45 | no |
| 1.6 nN (whole-cell) | 1.31 / 1.34 | 1.44 | no |
| **9.4 nN (protrusion)** | **1.34 / 1.29** | **1.48** | no |

**The crawl force does essentially nothing (f=0 → 1.35; f=9.4 nN → 1.34).** The A/A₀≈1.3 is almost
entirely proliferation; the crawl adds ~0. Cause: `substrate_crawl.substrate_crawl_forces` applies a
**radially-symmetric outward body force, uniform on every basal cell, recomputed each epoch from the
instantaneous centroid**. That is an *isotropic internal pressure* — cohesion (surface tension)
balances it exactly → a static force-balance equilibrium with **net expansion ≈ 0**. Raising f only
moves the equilibrium (stronger pressure ↔ stronger cohesion), it does not spread.

The missing ingredient is **coordination, not force magnitude**: (i) **persistence** (the current field
has *zero* — reset radially every epoch), (ii) **spatial correlation / alignment** (neighbours moving
coherently → distributed intercellular stress, the Trepat-2009 "tug of war"), (iii) **free-edge
polarization** (CIL: boundary cells crawl persistently into free space → a continuously advancing
front, not an equilibrium). A symmetric radial body force has none of these; SPP plithotaxis adds all three.

---

## 2. Mechanism — self-propelled-particle (SPP) plithotaxis on the center-based CBM

Each cell `i` carries an **in-plane polarization unit vector** `p_i = (cosθ_i, sinθ_i)` (persistent
cell state, the active-matter standard rendering of collective migration). The active traction is
applied the *faithful D-axis way* (`substrate_crawl` lineage): **substrate-reacted** (drag = clutch
γ_cell, not cohesion-competing), **in-plane** (z-force ≡ 0, cells crawl across the dish), on **basal**
cells:

```
F_active,i = f_active · p̂_i · w_basal,i          (z-component ≡ 0)
```

`f_active` = the anchored self-propulsion value, bracket {1.6, 5.0, 9.4} nN (1.6 net-whole-cell, **5.0 =
MCF10A v_m·γ_cell — the principled breast-epithelial anchor**, 9.4 protrusion; with the yield/ductile
cohesion none ejects, all admissible). `w_basal,i` = the same half-cosine basal membership weight as
`substrate_crawl` (geometry, not a constant).

**Per-epoch polarization update — Smeets-2016 CIL-SPP core** (cadence Δt = `epoch_steps·dt`, identical
to the existing crawl recompute — host-side numpy + one cKDTree, same as proliferation already does; the
per-*step* force is a fixed `SettableForce` vector → GPU-resident, integrator untouched, BAOAB freeze
intact). The continuous Smeets form `θ̇_i = −f_cil·(θ_i − θ_i^free) + √(2D_r)·η_i` discretised over Δt:

```
θ_i ← θ_i + Δθ_CIL + Δθ_noise   [ + Δθ_align ]      (bracketed term default-OFF, see §3)

Δθ_CIL   = f_cil · w_edge,i · Δt · sin(θ_free,i − θ_i)   θ_free,i = arg(p_i_pos − weighted-centroid_{N(i)})
Δθ_noise = sqrt(2 · D_r · Δt) · ζ          ζ ~ N(0,1)                     rotational diffusion (persistence τ=1/D_r)
Δθ_align = J · Δt · sin(θ̄_i − θ_i)        θ̄_i = arg( Σ_{j∈N(i)} p_j )    Vicsek (OPTIONAL, J=0 default)
```

- `N(i)` = neighbours within the cohesion interaction range (reuse the cKDTree the epoch loop builds).
- `θ_free,i` = direction away from the **weighted mean position of contacting neighbours** (toward open
  dish) — Smeets' exact CIL form; the *local* centroid (not global) is what makes a finger/front not a balloon.
- `w_edge,i` = edge-ness ∈ [0,1] from local coordination, e.g. `clip(1 − n_neighbours/n_bulk, 0, 1)`
  (geometry: a fully-surrounded bulk cell has w_edge≈0 → only diffuses; a free-boundary cell has w_edge≈1
  → driven outward). CIL is the literature mechanism; the edge detection is geometry. Smeets gates CIL on
  contact; we additionally weight by free-edge openness so bulk cells (no free space) are not spuriously driven.

**Why this spreads where the radial field did not:** persistence (D_r, τ=20 min) keeps a cell crawling in
one direction instead of being re-randomised every epoch; the catch+turnover **cohesion** elastically
couples neighbours so their persistent motions become *correlated* over ~ξ (Garcia 2015: this is where the
~200 µm coherence comes from — no Vicsek term needed); CIL keeps the free edge persistently advancing into
open space (the directed front). The sheet does **not** reach a static cohesion-balanced equilibrium — it
is a continuously migrating, ductile-cohesion (`yield_remodel`) sheet whose footprint grows over the 60 h
window → A/A₀ can rise toward the PI magnitude. (Contrast: the §D radial field is instantaneous, zero-
persistence, and globally symmetric → isotropic pressure → static equilibrium → net spread ≈ 0.)

---

## 3. Anchoring — ALL inputs measured; ξ & stress are EMERGENT targets (no-magic-number block)

✅ RESOLVED 2026-06-04 by deep-research `wf_e0df20b3-470` (anchors + citation-integrity flags in
`outputs/tag_kb/SE_REGISTRATION_CANDIDATES_2026-06-04_collective-migration.md`). Every model INPUT is a
measured value; ξ and intercellular stress are EMERGENT overlay-validation targets (not fitted); the
Vicsek term is DROPPED on literature grounds. No value is tuned to A/A₀.

| param | meaning | anchored value | source (system) | role |
|---|---|---|---|---|
| `D_r` | rotational diffusion | **0.05 min⁻¹ = 8.33e-4 s⁻¹** (τ=1/D_r=20 min) | Smeets 2016 (**MCF10A** MSD fit) | persistence INPUT |
| `f_cil` | CIL repolarization rate | **0.1 min⁻¹ = 1.67e-3 s⁻¹** | Smeets 2016 | free-edge outward drive INPUT |
| `f_active` | per-cell self-propulsion | bracket **{1.6, 5.0, 9.4} nN** | motility_bridge + **5.0 = v_m(1µm/min,MCF10A Smeets)·γ_cell** | 3 arms; 5.0 nN = principled breast-epi anchor, < cohesion 6.5 nN |
| `J` (Vicsek) | neighbour alignment | **OMITTED (default-off)** | Garcia 2015 / Henkes 2020 | correlations EMERGE from persistence+cohesion; imposing J = unanchored tuned term |
| ξ | velocity correlation length | **~200 µm TARGET (overlay)** | Petitjean 2010 (MDCK) | emergent self-consistency diagnostic, NOT fitted |
| stress | intercellular tension | **>300 Pa tensile, edge→in build-up (overlay)** | Tambe 2011 | emergent stress-field validation |
| ψ | f_cil/(2 D_r) | **= 1.0** (reproduces Smeets MCF10A) ✓ | Smeets 2016 | derived consistency check (0 free params) |

`τ = 1/D_r` = white-angular-noise definitional relation (Bi/Manning 2016; Szabó 2006). ⚠️ `τ=1/D_r` is the
*polarity-vector* autocorrelation (= single-cell persistence only in the dilute limit; renormalised by
neighbours in a dense monolayer) — but `D_r` is the correct SPP INPUT and Smeets fit it directly from
MCF10A MSD, so we feed `D_r`, not a monolayer τ. ⚠️ ξ has 3 inequivalent definitions (Petitjean PIV /
Angelini ξh / Garcia ξvv) — the emergent-ξ estimator uses the Petitjean PIV velocity-correlation form to
match the ~200 µm anchor.

**⭐ Literature refinement of the ratified mechanism (surfaced to PI).** PI ratified "Full SPP plithotaxis"
incl. Vicsek neighbour alignment. Deep-research finds (Garcia 2015 PNAS + Henkes/Marchetti 2020 Nat Commun,
both 3-0 verified) that the ~200 µm velocity correlations **EMERGE from per-cell traction persistence +
cohesive elastic coupling WITHOUT any explicit Vicsek alignment**, and the closest breast-epithelial-anchored
model (Smeets 2016, MCF10A) has **no Vicsek term** — it is exactly `persistence(D_r) + CIL(f_cil) + cohesion`.
Imposing a Vicsek `J` would (a) add an UNANCHORED tunable (ξ emerges without it ⇒ J would be fit to ξ,
anti-magic-number) and (b) double-count the coordination the cohesion already supplies. **Decision
(literature-first):** the production mechanism is the **Smeets-2016 CIL-SPP core** (persistence + CIL,
correlations emergent from the catch+turnover cohesion). The Vicsek term is still *implemented but
default-OFF* (PI's "full" machinery exists and is togglable); the first bracket measures the emergent ξ — if
it lands ~200 µm without alignment (as Garcia predicts), the drop is confirmed. The core code is identical
either way, so this is not a blocker.

---

## 4. Integration plan (additive, freeze-safe, default-off)

1. **New module** `spheroid/plithotaxis.py` — `PolarizationField` (stateful: holds `θ`, RNG; `.update(positions, neighbour_index, dt_epoch)`; `.forces(positions, f_active, z_substrate, basal_band)` → (N,3) in-plane substrate-reacted force). Pure mechanism; full Sanity-Gate docstring.
2. **`proliferation.run_growth_pooled`** — add `crawl_mode='plithotaxis'` (alongside `'edge'`/`'basal'`). It (a) constructs the `PolarizationField` once, (b) each epoch calls `.update(...)` using the cKDTree the loop already builds, then `.forces(...)` → `edge_force.set_vectors(...)`. `crawl_mode='basal'`/`'edge'` paths unchanged; default stays `'edge'` (regression-green). Requires a substrate (in-plane crawl).
3. **Driver** `scripts/layer2_plithotaxis.py` — bracket {f0, whole-cell 1.6 nN, protrusion 9.4 nN} with `cohesion='catch', yield_remodel=True`, substrate on, at R₀≈153 µm (N₀=4000, ≥2 seeds), GPU. Records `aa0_core`, `aa0_raw`, ejected, **+ emergent velocity-correlation-length ξ and front velocity** (overlay-check vs the anchored measured values — these are *diagnostics*, the model is not fitted to them). Auto-viz at run end (production-driver-auto-viz rule).
4. **Tests** `tests/test_plithotaxis.py` — Sanity Gate: f_active=0 ⇒ zero force (reduces to D-axis baseline); z-component ≡ 0; persistence ⇒ ⟨cosΔθ⟩ decays with the right τ; alignment ⇒ correlation length increases with J; CIL ⇒ free-edge cells polarize outward; pure-noise (J=β=0) ⇒ isotropic diffusion, no net spread (recovers the §D "radial does nothing" null as a *limit*, confirming the new spread comes from coordination).

**Hard-rule compliance:** mechanistic (not lumped) ✓ PI-ratified ✓; `integrator/baoab*.py` untouched (force is a `SettableForce`, per-step application only) ✓; per-epoch host update is the same pattern as the existing crawl/proliferation (GPU-main respected — no new per-step cpu_local_snapshot) ✓; overlay-only (ξ/τ/stress anchored to measurement, A/A₀ never fitted) ✓; physiological baseline (runs on the full catch+substrate+turnover stack at its setpoints) ✓; default-off for backward-compat, production turns it ON at the anchored value ✓.

---

## 5. Sanity Gate (to be the module docstring)

- **Dimensional:** θ [rad]; D_r [1/s]; J,β [1/s]; Δt [s]; f_active [N]; force [N]. `sqrt(2 D_r Δt)` dimensionless·rad ✓.
- **Boundary:** f_active=0 ⇒ zero force (= D-axis baseline). J=β=0 ⇒ pure rotational diffusion ⇒ isotropic, COM diffusive not ballistic, net spread→0 (recovers §D null). β large, J=0, no noise ⇒ pure CIL radial-outward (≈ the old `substrate_crawl`, as a *limit*).
- **Sign/sense:** larger β ⇒ stronger outward front ⇒ more spread; larger D_r (smaller τ) ⇒ shorter persistence ⇒ less spread; larger J ⇒ longer ξ ⇒ more coherent fronts.
- **Conservation:** z-force ≡ 0 (substrate reacts z); in-plane net force is *not* required zero (the substrate is an external reservoir — the front is genuinely propelled, unlike the zero-mean ActiveMotility); count conserved.
- **Numerical:** per-epoch angular step `|Δθ|` must stay < O(1) rad for the discrete update to track the SDE — assert `J·Δt, β·Δt, sqrt(2 D_r Δt) ≲ 0.3`; if `epoch_steps` makes Δt too coarse, sub-step the angle update (host-side, cheap) — does NOT touch the integrator.
- **Measurement consistency:** emergent ξ (from the simulated velocity field) and v_front are reported against the anchored measured ξ, v_front as an *overlay self-consistency check*, never fitted.

---

## 6. Open / next

- ⏳ Fill §3 anchors from deep-research (run ID wf_e0df20b3-470); write SE_REGISTRATION_CANDIDATES.
- Then implement §4 (1→2→3→4), Mac smoke, gbook GPU bracket, REPORT §E + figure, Notion Dev Log.
- D2 σ-bridge absolute value to be re-confirmed once the main session finalizes KU-3.5 (soft coupling, not a prerequisite for driving — per PI).
