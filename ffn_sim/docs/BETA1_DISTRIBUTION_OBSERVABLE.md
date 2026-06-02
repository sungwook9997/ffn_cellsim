# β1-distribution observable — design (PI-exp Layer-1 single-cell readout)

> **STATUS: DESIGN + standalone observable — 2026-06-02. Prep, NOT wired.** Parallel
> session (`AUTONOMOUS_LOG_2026-06-02.md`). Ships a tested, import-isolated geometry
> module `bridge/clutch_spatial.py`. Completes the PI-exp Layer-1 single-cell toolkit
> alongside `LIGAND_IDENTITY_FA_SPEC.md` (ligand identity) + the existing traction /
> KU-3.5 readout. Read-only post-hoc analysis — touches NO runtime, NO Lead file.

## The experiment observable

The PI study reads **integrin β1 immunofluorescence** spatial pattern, which tracks the
three ligand conditions (PI-exp map §Layer-1 mapping table):

| PI condition | β1 IF pattern | platform analog |
|---|---|---|
| Bare | weak / **diffuse** | low engaged-clutch count, area-proportional spread |
| Pre | **peripheral** | engaged clutch ringed at the contact edge |
| Lam4 | **uniform** | engaged clutch evenly spread, no edge bias |

This is a **spatial distribution of the engaged integrin clutch** — exactly the data
the platform's FA/clutch already carries (engaged-integrin tags + positions). The IF
data is **qualitative** (diffuse/peripheral/uniform), so the observable validates the
**pattern / ordering**, never an absolute number (overlay-only, literature-first).

## Metrics (over engaged-clutch positions projected to the substrate plane, en-face)

Let `{r_i}` be the xy positions (substrate at z≈0 ⇒ the xy projection is the en-face IF
view) of the **engaged** integrins; `c` their centroid; `r̃_i = |r_i − c|`.

1. **Edge-localization fraction `f_edge`** = fraction with `r̃_i > ρ · r_max`
   (`r_max = max r̃_i`, default `ρ = 0.7`). **Reference:** a uniform areal disk gives
   `f_edge_uniform = 1 − ρ²` (= 0.51 at ρ=0.7). `f_edge ≫ 0.51` ⇒ **peripheral** (Pre);
   `f_edge ≈ 0.51` ⇒ **uniform** (Lam4); `f_edge ≪ 0.51` ⇒ central. The single
   cleanest peripheral-vs-uniform separator.
2. **Angular uniformity `angular_cv`** = CV (std/mean) of engaged counts over
   `n_sectors` equal angular bins about `c`. **Low** ⇒ uniform (Lam4); **high** ⇒
   patchy / clustered. Distinguishes uniform from clumped at equal `f_edge`.
3. **Engaged count `n_engaged`** = overall β1 IF intensity proxy. Bare (weak) ⇒ low;
   Pre/Lam4 ⇒ higher.
4. **(documented extension, not in the tested module)** Clark-Evans NN index
   `R = d̄_obs / d̄_CSR`, `d̄_CSR = 0.5/√λ`, `λ = n/A` — clustering vs CSR vs dispersed.
   Needs an edge-correction + an area estimate (convex hull); specified here, deferred
   to keep the shipped module robust + trivially testable. Add when a richer clustering
   readout is needed.

## Phenotype acceptance (qualitative ordering — overlay targets, NOT fits)

| Condition | n_engaged | f_edge (vs 0.51) | angular_cv |
|---|---|---|---|
| Bare (diffuse) | low | ≈ uniform | moderate (sparse noise) |
| Pre (peripheral) | higher | **> uniform** (edge ring) | elevated (ring ≠ filled) |
| Lam4 (uniform) | higher | ≈ uniform | **low** |

Validation = reproduce the **ordering** `f_edge(Pre) > f_edge(Lam4) ≈ f_edge_uniform`
and `angular_cv(Lam4) < angular_cv(Pre)`; compare PI IF qualitatively at the end.

## The shipped module — `bridge/clutch_spatial.py`

```python
from ffn_sim.bridge.clutch_spatial import beta1_distribution_metrics
m = beta1_distribution_metrics(positions_xyz, engaged_mask, center=None,
                               edge_rho=0.7, n_sectors=12)
# -> {n_engaged, f_edge, f_edge_uniform_ref, angular_cv}
```

Pure numpy geometry (no HOOMD, import-isolated). Operates on a positions array + an
engaged boolean mask + an optional center (defaults to the engaged centroid) — so it is
**layout-agnostic**: the Lead wires it to the real engaged-integrin tags
(`bridge/fa.py` clutch range) at analysis time without this module knowing the tag
layout. Verified GREEN in isolation (`tests/test_clutch_spatial.py`).

## Wiring (Lead analysis step — NOT done here)

A post-hoc analysis over the production GSD / a snapshot: slice the engaged-integrin
tags + xy positions, call `beta1_distribution_metrics`, emit a figure per condition
(per the CLAUDE.md visualize rule) overlaying the three conditions' `f_edge` / `cv`.
No runtime change; a measurement, not a force.

## Open items for PI / Lead

- [ ] Confirm the en-face (xy) projection is the right IF analog (vs the full basal
      contact patch if the cell is not flat on the substrate).
- [ ] Set `ρ`, `n_sectors` conventions (defaults 0.7 / 12 are reasonable, not tuned).
- [ ] Decide if the Clark-Evans / Ripley extension is needed for the cohesion story.
