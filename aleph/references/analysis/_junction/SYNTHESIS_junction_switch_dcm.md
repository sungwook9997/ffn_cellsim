---
id: SYNTHESIS_junction_switch_dcm
topic: J1+J2 synthesis → concrete DCM implementation plan for the bulk-pressure-driven cadherin→integrin junction switch + pressure-release outward spreading
kind: implementation-plan (junction track, DCM tier)
year: 2026
ffn_relevance: High (DCM-tier mechanism design, gates existing modules)
parents: [J1_bulk_pressure_unjamming, J2_cadherin_integrin_switch]
ffn_themes: [junction, DCM, bulk-pressure, jamming-unjamming, cadherin-integrin-switch, FA-clutch, spreading, mechanism-design]
entities: [DcmTurgorForce, CellCellAdhesion, DcmSubstrateForce, FaClutchForce, FaCatchSlipUpdater, CadherinTransJunctionUpdater, IntegrinBondUpdater, ProliferationUpdater]
modules_gated: [cell/dcm.py, cell/dcm_prolif.py, cell/dcm_ecm.py, junction/cadherin.py, bridge/integrin_bonds.py]
keywords: [local pressure gate, contact stress, neighbour crowding, cadherin weakening, integrin strengthening, unjamming, pressure release, outward spreading, shape index oracle]
tags: ["#junction", "#DCM", "#bulk-pressure", "#jamming-unjamming", "#cadherin-integrin-switch", "#mechanism-design", "#implementation-plan"]
has_transferable_params: true
status: design (PI sign-off required before any new ValidationGate; new module is additive/default-off until production config turns it ON at the physiological band)
---

# [SYNTHESIS] Bulk-pressure-driven cadherin→integrin junction switch + pressure-release outward spreading — concrete DCM implementation

**Tags:** #junction #DCM #bulk-pressure #jamming-unjamming #cadherin-integrin-switch #mechanism-design #implementation-plan

> Combines **J1** (a compacting aggregate builds 0.1–20 kPa solid stress; the
> proliferation/motility gate engages ~0.5 kPa and saturates ~5 kPa; unjamming is a
> sharp shape-index threshold; release re-fluidizes and re-licenses spreading) with
> **J2** (cadherin AJ and integrin FA are two clutch populations competing for one
> finite actin/vinculin pool; rising tension/stiffness reallocates the budget toward
> integrin while cadherin softens; both arm at a shared ~5 pN gate; no single sharp
> nN "switch tension" exists — it is ratio/stiffness-set) into a **runtime DCM
> mechanism**: a per-cell **local-pressure gate** that, above a threshold band,
> **weakens cell-cell cadherin adhesion** and **strengthens cell-ECM integrin/FA
> clutch**, unjamming the cluster so the released pressure drives outward spreading.

---

## 0. Design constraints carried in from CLAUDE.md

- **No abstractions (PI 2026-05-19).** The switch is NOT a phenomenological scalar
  flag per cell. It is driven by a **mechanically measured** per-cell local pressure
  (already computed inside `DcmTurgorForce`) and a **measured contact stress** (already
  computed inside `CellCellAdhesion`), and it acts by **retuning the existing explicit
  clutch populations** (cadherin adhesion well + FA clutch stiffness/binding), not by a
  lumped order parameter. J2's "shared finite actin/vinculin pool" competition is the
  fine-grained ideal; in the DCM tier (where cadherin and integrin are adhesion-energy
  wells / clutch springs, not particle-resolved adaptors) we implement the **reallocation
  of one conserved adhesion budget** — see §4 budget-conservation — which is the
  DCM-faithful image of the same competition.
- **Physiological operating point (PI 2026-06-04, HARD).** The gate trips at a
  PHYSIOLOGICAL pressure band (~0.5 kPa onset, ~5 kPa saturation — Dolega 2021 eLife
  e63258; Delarue 2014 Biophys J 107:1821), measured FROM the resting baseline (turgor
  `turgor_dP0` already 133 Pa in `ResolvedDCM`). The new updater is additive/default-off
  for regression cleanliness, but the **production config MUST turn it ON** at the band.
- **Vertex/SPV are oracles only.** Shape index `q*≈3.81` (2D) / `SI≈5.4` (3D) and the
  Han core/periphery split (5.84 jammed / 6.6 unjammed; Han 2021 iScience 24:103252) are
  **acceptance gates measured off the GSD**, never imposed as a force law.
- **No empirical magic numbers.** Every band value below is anchored to a cited
  literature number with units; the dimensionless gate softness `n_hill` and budget
  exponents are grid-invariant shape parameters, not gate-fitting knobs (Magic-Number Block).

---

## 1. Mechanism summary (one paragraph)

As a DCM aggregate compacts, each cell's enclosed volume `Vc` drops below `V0` and its
neighbours crowd it. `DcmTurgorForce` already turns that compression into a per-cell
osmotic over-pressure `ΔP_compress_c = K_vol·(V0−Vc)/V0`; `CellCellAdhesion` already
computes the repulsive contact forces between cells. We sum these into a **per-cell local
pressure** `P_c` (Pa). When `P_c` rises through the **0.5 → 5 kPa band** (J1), a new
batched `JunctionSwitchUpdater` (a `hoomd.custom.Action`, same cadence convention as
`ProliferationUpdater` / `FaCatchSlipUpdater`) applies a smooth Hill gate `g(P_c) ∈ [0,1]`
that (a) **weakens the cadherin cell-cell adhesion** for that cell's nodes (scale the
`CellCellAdhesion` well depth `W_cc`/`F_adh0`, and for the explicit junction raise the
trans-dimer off-rate / drop `k_on`), and simultaneously (b) **strengthens the cell-ECM
integrin/FA clutch** (raise `DcmSubstrateForce` adhesion `W_cs`/`k_well`, and on the ECM
build raise `FaClutchForce.k_fa` + `FaCatchSlipUpdater` rebind `k_on`/`capture`). The
two moves are coupled by a **conserved adhesion budget** (J2: one finite actin/vinculin
pool — cadherin and integrin draw on the same currency), so weakening one funds
strengthening the other rather than adding free adhesion. The cluster's effective
adhesion-to-cortical-tension ratio rises on the ECM side and falls on the cohesion side →
the cohort crosses the shape-index rigidity line → **unjams** → the locally released
pressure (the cell can now lose neighbours and spread its footprint) **drives outward
spreading**. Because compression-arrest is fully reversible (Delarue 2014), when `P_c`
falls back below the band the gate relaxes and adhesion reverts — release → re-spread is
the gate running in reverse.

---

## 2. The per-cell local pressure signal `P_c` (the gate input)

Two mechanically-real contributions, both ALREADY computed in existing modules — we only
expose and combine them (no new physics, no new magic number):

**(2a) Compressive volumetric pressure** — read directly off `DcmTurgorForce`. Inside
`DcmTurgorForce.set_forces` (`cell/dcm.py:139`) the exact per-cell volume is

    Vc = (1/6) Σ_faces v0·(v1×v2)              # divergence theorem, already there
    ΔP_compress_c = max(0, K_vol·(V0 − Vc)/V0) # the compaction over-pressure (Pa)

Expose `self.dP_compress` (an `(n_cells,)` array, the `K_vol·(V0−Vc)/V0` term, clamped
≥0) as a cached attribute after each `set_forces`. With `K_vol = 1e3 Pa` (`ResolvedDCM`)
a 0.5 % volume compaction = 5 Pa; to reach the ~0.5 kPa onset a cell must be compressed
~50 % OR `K_vol` must sit at the physiological osmotic bulk modulus. **Use the J1
physiological `K_vol`**: the osmotic/poroelastic bulk modulus of a crowded cell aggregate
is ~1–10 kPa per unit volumetric strain (consistent with Stylianopoulos 2012 PNAS solid
stress 0.37–19 kPa over the compaction strains observed; Delarue 2014 cell-cell distance
−20 % at 10 kPa → bulk modulus ~10 kPa/0.5 ≈ ... order kPa). `ResolvedProlifDCM` already
runs `K_vol = 5e3 Pa`, which puts a 10 % compaction at 500 Pa — exactly the J1 onset.
**Adopt `K_vol = 5e3 Pa` as the DCM physiological default for the switch runs** (already
the prolif default; promote it in `ResolvedDCM`/`ResolvedDcmEcm` for switch production).

**(2b) Contact / crowding stress** — read off `CellCellAdhesion`. Inside its
`set_forces` (`cell/dcm_prolif.py:153`) the per-pair repulsive force `fmag[rep]` is the
contact load between two different cells. Accumulate the inter-cell **repulsive** force
magnitude per cell and divide by that cell's contact area to get a contact stress:

    σ_contact_c = (Σ_{j≠c} |F_rep,cj|) / A_contact_c       # Pa

where `A_contact_c ≈ (n_contact_nodes_c)·area_per_node`. This is the supracellular
crowding term (neighbour caging) that J1 §3 ties to the shape index. For `cell/dcm.py`
(per-cell-type LJ, no `CellCellAdhesion`) the same number is the sum of the repulsive
branch of the diff-cell LJ; simplest is to compute `σ_contact_c` from neighbour count ×
mean overlap in the new updater (see §5) so it works for BOTH dcm.py and dcm_prolif.py.

**Combine:**

    P_c = ΔP_compress_c + σ_contact_c                     # Pa, per cell

`P_c` is the gate input. It is mechanically measured every batch from the live snapshot —
no per-cell state variable is integrated, so it is regression-clean and ParticleSorter-safe
(read tag-ordered like every other DCM force).

---

## 3. The threshold band + the switch rule

**Band (J1, cited):**

| Quantity | Value | Source |
|---|---|---|
| Gate ONSET `P_lo` | **500 Pa** | Dolega 2021 eLife e63258 (Πd onset); Delarue 2014 Biophys J 107:1821 |
| Gate SATURATION `P_hi` | **5 000 Pa** | Dolega 2021 (migration −50 %, doubling 36→68 h saturates); Montel 2011 PRL 107:188102 / Delarue 2014 (5–10 kPa drastic arrest) |
| Resting baseline (no switch) | **~40–133 Pa** | project turgor `turgor_dP0`; J1 §1 |
| Shape-index unjamming oracle (2D) | **q\* ≈ 3.81** | Bi 2015 Nat Phys 11:1074; Bi 2016 PRX 6:021011; Park 2015 Nat Mater 14:1040 |
| Shape-index unjamming oracle (3D) | **SI ≈ 5.4** (core 5.84 jammed / periphery 6.6 unjammed) | Merkel & Manning 2018 NJP 20:022002; Han 2021 iScience 24:103252 |

**Smooth gate (no hard cliff — J2: the switch is graded/ratio-set, NOT a sharp nN
threshold).** A Hill function in log-pressure between `P_lo` and `P_hi`:

    g(P_c) = P_c^h / (P_c^h + P_mid^h),   P_mid = sqrt(P_lo·P_hi) ≈ 1581 Pa,   h = n_hill (≈ 4)

`g = 0` below the band (jammed, cadherin-dominant, no spreading), `g → 1` above it
(unjammed, integrin-dominant, spreading). `n_hill ≈ 4` is a grid-invariant softness, not
a gate-fit (set it once from the Dolega dose-response steepness; surface to PI if it ends
up being tuned to make a gate pass — Magic-Number Block).

**Switch rule (exact updater action, per cell c, every batch):**

    g = g(P_c)                                              # ∈ [0,1]
    # (a) WEAKEN cadherin cell-cell adhesion:
    W_cc_eff_c   = W_cc0   · (1 − f_cad · g)                # f_cad ≈ 0.8 (down to 20 %)
    # (b) STRENGTHEN cell-ECM integrin/FA clutch:
    W_cs_eff_c   = W_cs0   · (1 + f_int · g)                # f_int ≈ 1.5 (up to 2.5×)
    k_fa_eff_c   = k_fa0   · (1 + f_int · g)                # FA clutch stiffens (J2: stiff→FA)
    # budget coupling (J2 shared pool): the integrin gain is FUNDED by the cadherin loss
    # (see §4) so total adhesion currency is conserved, not created.

Sign-sense (Sanity Gate): higher `P_c` → larger `g` → **less** cadherin cohesion + **more**
ECM grip → adhesion-to-cortical-tension ratio rises on the substrate side → shape index
climbs above `q*` → unjamming → footprint grows (J1 §4 mechanistic chain). At `P_c ≤ P_lo`,
`g→0`, all `*_eff = *0` (baseline) → reversible (Delarue 2014). ✔

---

## 4. Budget conservation (the J2 shared-actin competition, DCM-faithful)

J2's load-bearing mechanism is that cadherin and integrin **compete for one finite
actin/vinculin pool** (Barcelona-Estaje 2024 Nat Commun 15:8824) — engaging cadherin
*siphons* actin from FA and vice versa; cell-pair intercellular force is a **constant
fraction 0.47±0.07 of total ECM traction** (Maruthamuthu 2011 PNAS 108:4708), i.e. they
co-scale around a homeostatic budget, not bond-for-bond. In the DCM tier we honour this by
making the integrin gain **funded by** the cadherin loss rather than free:

    B_total = W_cc0·A_cc + W_cs0·A_cs                       # conserved per-cell adhesion budget (J)
    ΔB_cad  = f_cad · g · W_cc0 · A_cc                      # released by weakening cadherin
    W_cs_eff_c = W_cs0 + (ΔB_cad / A_cs)·route_eff          # poured into integrin (route_eff ≤ 1)

`route_eff ≤ 1` is the fraction of freed cadherin budget that the FA side can actually
re-anchor this batch (capped by available basal contact / ligand density `ligand_density`)
— the rest is lost (not all freed actin reaches a productive FA), matching Barcelona-Estaje
2024's "more, weaker bonds siphon actin away, capping force." Keep `f_int` as the *cap*
(§3) and `route_eff` as the *realised* transfer; in the simplest first cut set
`f_int·g ≡ ΔB_cad/(W_cs0·A_cs)` so the two are literally the same quantity (full budget
conservation), which removes `f_int` as an independent magic number. **Validation oracle
(never runtime):** at steady state the emergent intercellular/traction ratio should track
~0.47 (Maruthamuthu 2011) across the gate sweep — a clean acceptance test, exactly the J2
prescription.

---

## 5. Exactly which modules to gate (and how)

### NEW module: `cell/junction_switch.py` — `JunctionSwitchUpdater(hoomd.custom.Action)`
One batched updater (cadence `switch_batch_steps`, default 200, same convention as
`ProliferationUpdater` and `FaCatchSlipUpdater`; fired via
`hoomd.update.CustomUpdater(action=…, trigger=hoomd.trigger.Periodic(switch_batch_steps))`).
Each `act(timestep)`:
1. read tag-ordered positions from `sim.state.cpu_local_snapshot` (ParticleSorter-safe);
2. read `turgor.dP_compress` (cached by `DcmTurgorForce`, §2a);
3. compute `σ_contact_c` from neighbour overlap (§2b) — works for both dcm.py and
   dcm_prolif.py, so it does not depend on `CellCellAdhesion` existing;
4. `P_c = dP_compress_c + σ_contact_c`; `g_c = Hill(P_c)`;
5. mutate the gated parameters in place on the existing force objects it holds references
   to (constructor takes `turgor`, `adhesion`, `substrate`, and optionally `fa_force`,
   `fa_updater`, `cadherin_updater`). It mutates ONLY force/updater parameters — never
   particle positions (BAOAB stays the sole stepper, CLAUDE.md hard rule);
6. store `g_c`, `P_c` as public arrays for the vis/oracle layer.

Resolved-params dataclass `ResolvedJunctionSwitch` (additive, `enabled: bool = False`):
`P_lo=500.0`, `P_hi=5000.0`, `n_hill=4.0`, `f_cad=0.8`, `route_eff=1.0`,
`switch_batch_steps=200`. All band values cited in §3.

### `cell/dcm.py`
- **`DcmTurgorForce`** — cache `self.dP_compress` (the `K_vol·(V0−Vc)/V0` array, clamped
  ≥0) at the end of `set_forces` (one line; it is already computed at `dcm.py:153`). This
  is the §2a pressure tap. Promote the switch-run default `K_vol = 5e3 Pa` (§2a).
- **`DcmSubstrateForce`** — make `W_cs`/`k_well` **per-cell mutable**: change the scalar
  `self.W_cs` / `self.k_well` to `(n_cells,)` arrays indexed by `cell_of_tag`, so the
  updater can raise integrin adhesion only for switched cells. (Currently node-uniform;
  the switch needs per-cell granularity — the substrate force already iterates nodes, so
  index `k_well[cell_of_tag[node]]`.) This is the **integrin/cell-ECM strengthen** lever.
- **diff-cell LJ `epsilon=W_cc`** — for the per-cell-type build, the updater rescales the
  diff-cell LJ `epsilon` per cell-type-pair `(ti,tj)` via `lj.params[(ti,tj)]` (HOOMD
  allows live param edits). This is the **cadherin/cell-cell weaken** lever for dcm.py.

### `cell/dcm_prolif.py`  (primary switch host — recommended first implementation)
- **`CellCellAdhesion`** — (i) accumulate per-cell repulsive contact force → expose
  `self.sigma_contact` `(N_max,)` for §2b; (ii) make `self.W_cc` / `self.F_adh0`
  **per-cell-pair scalable** via a `self.W_cc_scale` `(N_max,)` array (default 1.0) that
  the updater sets to `(1 − f_cad·g_c)`; multiply `F_adh0` by `0.5·(scale_i+scale_j)` per
  pair. This is the **cadherin weaken** lever (well depth `W_cc` ↓). The existing pair
  loop (`dcm_prolif.py:182-193`) already has `acell[ii]`,`acell[jj]` in hand to index it.
- **`DcmSubstrateForce`** — per-cell `W_cs`/`k_well` as above (integrin strengthen).
- **`ProliferationUpdater`** — couple, don't fight: a switched (unjammed, `g→1`) rim cell
  should keep dividing; a jammed core (`g→0`, contact-inhibited) stays quiescent — already
  consistent with J1 (core jammed/arrested, periphery unjammed/spreading). No change
  needed beyond optionally reading `g_c` to bias `p_div` (out of scope for v1).

### `cell/dcm_ecm.py`  (the explicit catch-slip FA clutch host)
- **`FaClutchForce.k_fa`** — let the updater scale `self.k_fa` (or a per-bond multiplier
  keyed by basal cell) by `(1 + f_int·g_c)`: the FA clutch **stiffens** under pressure
  (J2: stiff/high-tension → integrin wins; Kong 2009 JCB 185:1275 catch-bond reinforces
  10–30 pN). This is the mechanistic integrin-strengthen for the ECM build.
- **`FaCatchSlipUpdater`** — raise rebind `k_on` and/or `capture` for switched cells so
  the FA clutch population **grows** (more engaged integrin) when `g→1`. The Pereverzev
  `k_off(F)` catch-slip law is untouched (it stays the explicit bond physics — CLAUDE.md);
  the switch only modulates the engagement rate, the J2-faithful "more FA under load."
- The DCM↔ECM `DcmTurgorForce` instance supplies `dP_compress` exactly as §2a.

### `junction/cadherin.py`  (explicit E-cadherin trans-dimer catch-bond, when used)
- **`CadherinTransJunctionUpdater`** — when a cell's `g_c` is high, (i) drop its `k_on`
  (fewer new trans-dimers form) and/or (ii) add a `g`-scaled multiplier to the sampled
  `effective_k_off(F)` (faster rupture) for that cell's cadherins. This is the
  particle-resolved cadherin-weaken: the explicit Rakshit-2012 sliding-rebinding catch-slip
  law (`validation.cadherin_sliding_rebinding`) stays the bond physics; the switch only
  shifts the bind/unbind balance — the J2 "softened cadherin" without inventing a ΔG. Do
  NOT alter `k_trans` (the dimer stiffness is a fixed mechanical constant).

### `bridge/integrin_bonds.py`  (single-cell H.4 FA clutch, when on the same run)
- **`IntegrinBondUpdater`** — symmetric to the cadherin updater: raise rebind `k_on` /
  `capture_radius_R_FA` for the switched cell so integrin engagement grows under pressure.
  The Pereverzev catch-slip `k_off(F)` is untouched (the explicit bond physics, KU-2.5).
  This keeps the single-cell FA track consistent with the DCM-tier switch if both run.

---

## 6. Band values (the numeric ledger to encode)

| Symbol | Value | Meaning | Source |
|---|---|---|---|
| `P_lo` | **500 Pa** | gate onset (cadherin starts to yield, integrin starts to win) | Dolega 2021 eLife e63258; Delarue 2014 BpJ 107:1821 |
| `P_hi` | **5 000 Pa** | gate saturation (full unjamming / integrin dominance) | Dolega 2021; Montel 2011 PRL 107:188102 |
| `P_mid` | **≈1 581 Pa** (=√(P_lo·P_hi)) | Hill half-max | derived (grid-invariant) |
| `n_hill` | **≈4** | gate steepness (dose-response) | Dolega 2021 dose-response; grid-invariant shape param |
| `K_vol` (switch runs) | **5×10³ Pa** | osmotic bulk modulus → 10 % compaction = 500 Pa = onset | `ResolvedProlifDCM` default; consistent w/ Stylianopoulos 2012 PNAS 0.37–19 kPa |
| `f_cad` | **0.8** | max cadherin weakening (W_cc → 20 %) | J2: E-cad softens/loses in EMT (Canel 2013 JCS 126:393; Figueiredo 2021) |
| `f_int` / `route_eff` | **budget-conserved** (≡ ΔB_cad transfer) | integrin gain funded by cadherin loss | J2 shared-pool (Barcelona-Estaje 2024 Nat Commun 15:8824); 0.47 ratio (Maruthamuthu 2011 PNAS 108:4708) |
| `turgor_dP0` | **133 Pa** | resting baseline (no switch) | project; J1 §1 |
| `q*` (2D) / `SI` (3D) | **3.81 / 5.4** | unjamming oracle (measured, not imposed) | Bi 2015/2016; Merkel-Manning 2018; Han 2021 |

---

## 7. Validation oracles (acceptance gates — measured off the GSD, NEVER runtime)

1. **Shape-index crossing.** Measure 2D `q = ⟨P/√A⟩` (monolayer / dcm_prolif) or 3D
   `SI = surface/volume^(2/3)` (dcm.py 3d cluster) per cell off the trajectory. PASS:
   pre-switch cohort `q < 3.81` (jammed), post-switch / released-periphery `q > 3.81`
   (unjammed). Han target: core ~5.84 vs periphery ~6.6 (3D). Same observable family as
   the existing L2 hull/shape metrics + Roffay `ΔP=σ(1/R+1/R')` oracle.
2. **Pressure→spreading dose-response.** Footprint `A(t)/A0` (already sampled by
   `h7_dcm_spheroid_spread.py` / `h7_dcm_prolif_spread.py`) vs imposed/emergent `P_c`:
   monotone rise across 0.5→5 kPa, saturating ~5 kPa (Dolega: migration −50 % at 5 kPa).
3. **Reversibility.** Lower confinement (or `K_vol`/external stress) so `P_c` drops back
   below `P_lo`; verify `g→0`, `W_cc`/`W_cs` revert, `q` crosses back, footprint dynamics
   resume — Delarue 2014 "completely reversible"; Han 2020 stress-release volume flip.
4. **Cadherin/integrin force ratio.** Emergent intercellular / total-ECM-traction ratio
   tracks **0.47±0.07** across the gate sweep (Maruthamuthu 2011) — the J2 homeostasis test.
5. **Figures (CLAUDE.md visualize-at-closeout):** extend `scripts/h1_h2_vis.py` (or a
   `h7_junction_switch_vis.py` entry point) → `outputs/h7/figs/junction_switch/`:
   (a) `P_c` map over the cluster, (b) `g_c` gate field, (c) `q`/SI vs `P_c` with the 3.81/5.4
   line overlaid, (d) `A/A0` vs `P_c` dose-response with the 0.5/5 kPa band overlaid,
   (e) reversibility hysteresis loop. Per-realisation thin lines + ensemble mean.

---

## 8. Sanity Gate (record in the new module docstring per CLAUDE.md)

- **Dimensional:** `K_vol·(V0−Vc)/V0` [Pa]; `σ_contact = ΣF/A` [N/m²=Pa]; `g` [–];
  `W_cc`,`W_cs` [J]; `k_fa` [N/m]; `Δt_batch = switch_batch_steps·dt` [s]. ✔
- **Boundary:** `enabled=False` → strict no-op (all `*_eff = *0`, snapshot unchanged).
  `P_c ≤ P_lo` → `g=0` → baseline. `P_c ≥ P_hi` → `g=1` → full switch (W_cc→20 %,
  budget fully routed). Single cell / no contacts → `σ_contact=0`. ✔
- **Conservation:** the switch mutates force *parameters*, not positions → injects no
  net momentum; the budget coupling (§4) conserves total adhesion currency (no free
  adhesion created). ✔
- **Sign/sense:** `dP_compress`,`σ_contact` clamped ≥0; `g` monotone↑ in `P_c`; cadherin
  ↓ and integrin ↑ with `g` (the unjamming direction). ✔
- **Numerical:** read tag-ordered (ParticleSorter-safe); RNG-free (deterministic gate);
  per-cell arrays sized to `n_cells`/`N_max`; no double-stepping. ✔
- **Measurement consistency:** the `P_c` that drives the gate is the SAME `K_vol·(V0−Vc)`
  pressure HOOMD applies via `DcmTurgorForce` + the SAME repulsive force `CellCellAdhesion`
  applies — no second pressure model. ✔

---

## 9. Implementation order (recommended)

1. `DcmTurgorForce.dP_compress` cache (1 line, `cell/dcm.py`) + `CellCellAdhesion.sigma_contact`
   accumulation (`cell/dcm_prolif.py`) — the pressure taps. STATIC test: known compression → known `P_c`.
2. Per-cell `W_cs`/`k_well` in `DcmSubstrateForce` + per-cell `W_cc_scale` in `CellCellAdhesion`
   (default 1.0 → bit-identical to current behaviour; regression gate).
3. `cell/junction_switch.py` + `ResolvedJunctionSwitch` (`enabled=False` default).
4. Wire into `build_prolif_simulation` (and `build_dcm_simulation`, `build_dcm_ecm_*`) as an
   additive, default-off updater; production config turns it ON at the §6 band.
5. `junction/cadherin.py` + `bridge/integrin_bonds.py` `k_on`/off-rate hooks (explicit-bond builds).
6. Vis + the §7 oracles. PI sign-off before any new `ValidationGate` row (gate-contract change).

---

## 10. Caveats

- **DCM-tier honesty.** In the DCM tier cadherin/integrin are adhesion-energy wells /
  clutch springs, not particle-resolved talin/α-catenin/vinculin. The §4 "shared budget"
  is the DCM-faithful image of J2's actin-pool competition, NOT the literal pool — the
  literal pool lives in the single-cell H.4 FA track + `junction/cadherin.py`. Flag this
  explicitly so the switch is not over-claimed as the fine-grained competition.
- **No single nN switch tension exists** (J2): encoded as the graded Hill gate + budget
  competition, never a constant threshold force.
- **2D vs 3D shape index** (J1): dcm_prolif monolayer → 3.81; dcm.py 3d cluster → 5.4. Do
  not mix conventions in the oracle.
- **`K_vol` is the load-bearing pressure-scale choice.** At `K_vol=1e3 Pa` the onset needs
  ~50 % compaction (unphysical); at `5e3 Pa` ~10 % hits 500 Pa (physiological). This is a
  physiological-operating-point call, surface to PI if the aggregate's emergent compaction
  strain does not put `P_c` in the 0.5–5 kPa band at realistic densities.
- **PI sign-off** required for the new ValidationGate (shape-index crossing / dose-response /
  0.47-ratio) before any production claim; the module itself is additive/default-off until then.

---

## 11. Primary citations (verify `source_audit` verdict OK before deliverable use)

- Dolega ME et al. *eLife* 10:e63258 (2021). doi:10.7554/eLife.63258 — onset 500 Pa, sat 5 kPa, migration −50 %.
- Delarue M et al. *Biophys J* 107(8):1821 (2014). doi:10.1016/j.bpj.2014.08.031 — reversible compression arrest; −20 % cell-cell dist at 10 kPa.
- Montel F et al. *Phys Rev Lett* 107:188102 (2011). doi:10.1103/PhysRevLett.107.188102 — 5–10 kPa drastic arrest.
- Stylianopoulos T et al. *PNAS* 109(38):15101 (2012). doi:10.1073/pnas.1213353109 — endogenous solid stress 0.37–19 kPa.
- Bi D et al. *Nat Phys* 11:1074 (2015) doi:10.1038/nphys3471; *Phys Rev X* 6:021011 (2016) doi:10.1103/PhysRevX.6.021011 — q*≈3.81.
- Park JA et al. *Nat Mater* 14:1040 (2015). doi:10.1038/nmat4357 — q*≈3.81 in human bronchial epithelium.
- Merkel M, Manning ML. *New J Phys* 20:022002 (2018). doi:10.1088/1367-2630/aaaa13 — 3D SI≈5.4.
- Han YL et al. *Nat Phys* 16:101 (2020) doi:10.1038/s41567-019-0680-8; *iScience* 24:103252 (2021) doi:10.1016/j.isci.2021.103252 — core 5.84 / periphery 6.6; stress-release flip.
- Barcelona-Estaje D et al. *Nat Commun* 15:8824 (2024). doi:10.1038/s41467-024-53107-6 — cadherin/integrin shared-actin clutch competition.
- Maruthamuthu V et al. *PNAS* 108:4708 (2011). doi:10.1073/pnas.1011123108 — intercellular = 0.47×traction.
- Canel M et al. *J Cell Sci* 126:393 (2013) doi:10.1242/jcs.100115; Figueiredo J et al. *Gastric Cancer* 25:124 (2021) doi:10.1007/s10120-021-01239-9 — EMT cadherin→integrin switch.
- Kong F et al. *J Cell Biol* 185:1275 (2009). doi:10.1083/jcb.200810002 — integrin α5β1 catch-bond 10–30 pN.
- In-tree modules gated: `cell/dcm.py`, `cell/dcm_prolif.py`, `cell/dcm_ecm.py`, `junction/cadherin.py`, `bridge/integrin_bonds.py`.
