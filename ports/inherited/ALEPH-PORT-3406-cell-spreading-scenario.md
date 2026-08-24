# ALEPH-PORT-3406 — cell spreading on a 2D ECM, and the treadmill it must be distinguishable from

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3406` |
| Lane | `f70de564` scenario lane |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Amended | `2026-07-31` — §4d, §4e, §7a, §8a, §9a and a rewritten §14.11, **written before the substrate-adhesion code**, same rule |
| Port class | `RE-DERIVED` |
| Exists because | the PI asked to see a single cell **move**, not a cell running in place |

---

## 0. The instruction this entry is written against

The PI, in Korean, on the difference between two things a cell model can do:

> 크롤이 런닝머신 위에서 되는 것이라면 그거 말고 Spreading처럼 단일세포도 움직이는 거 그거를
> 보아야하고

*If the crawling happens on a treadmill, that is not what I want — I want the kind where the single
cell actually moves, like spreading.*

That is a precise technical complaint and it sets the observables. A cell model can produce
retrograde flow, traction and adhesion turnover with its **centre of mass exactly where it started**;
the substrate scrolls past and nothing translocates. Traction is therefore **not** the primary
measurement here. The primary measurements are:

1. **contact area against the substrate, versus descent step** — the classic spreading observable,
   which rises and saturates;
2. **centre-of-mass displacement of the cell, in the lab frame**;
3. traction, third, as the thing that produced them.

And the treadmill is built as a **control that must fail**: a world in which the clutch tensions are
of the same order while the lab-frame footprint and the lab-frame centre of mass do not move.

The PI's second instruction, same session, governs §12:

> 단순히 데이터들만을 믿지 말고 시각화도 같이 시켜서 시각화로도 판단해야돼

*Do not trust the numbers alone — render it and judge by eye.*

## 1. Aleph API

```python
from aleph.scenarios.spread import (
    SpreadingConfig,          # every scenario knob, with its provenance recorded per field
    MaterialPointCoupling,    # site-indexed connector -> full-array connector adapter
    build_spreading_world,    # -> (CellWorld, SpreadingHandles)
    run_spreading,            # -> SpreadingRun: the observable time series
    contact_area_um2,         # the footprint, against a LAB-FRAME plane
    contact_footprint_radius_um,
    centre_of_mass_um,        # vertex-area-weighted, not a plain node mean
    clutch_traction_pn,
    render_frames,            # pure-stdlib PNG; no third-party plotting dependency
)
```

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` |
| Source path / symbol | **none named, none read** |
| Read from | **neither.** No file under `/Users/sw1/ffn_cellsim` was opened by this lane. |

`PLAN.md` §7 records that the reference project's whole-cell solve has a declaration and **zero
production tasks**. There is no spreading result there to port, and that absence is also the warning
this entry is written under: a world that assembles and produces a curve is not a world that spread.

## 3. Why source-derived porting beats clean-room

**It does not, and nothing was ported.** Every force law used here already exists in this repository
and is gradient-checked in its own owner's controls: the Helfrich membrane, the elastic cortex, the
osmotic envelope, the collagen network, the lamellipodial branched network, and the focal-adhesion
series joint. This module composes them and measures. It **defines no force law of its own**, and
`test_the_scenario_defines_no_force_law` is the control that keeps it that way.

## 4. What this is

A scenario, in three configurations, over one composition.

**The composition** — six owners, four couplings, assembled through `aleph.scenarios.world`:

| Owner | Module | Role here |
|---|---|---|
| `membrane` | `vertical/membrane.py` (via `build_vertical`) | the cell surface whose footprint is measured |
| `cortex` | `vertical/cortex.py` (via `build_vertical`) | the surface that resists flattening |
| `turgor` | `vertical/pressure.py` (via `build_vertical`) | the osmotic envelope, a field owner loading the cortex |
| `ecm` | `vertical/ecm_network.py` — `CollagenNetwork` | **the substrate.** A planar woven collagen sheet, perimeter-clamped |
| `lamellipodium` | `vertical/lamellipodium.py` | the protrusive population at the basal rim; the clutch's actin side |
| `focal_adhesion` | `vertical/focal_adhesion.py` | the clutch graph; geometry-less, supplies the series compliance and the gate |

| Coupling | Class | Registry row |
|---|---|---|
| membrane↔cortex | `CortexMembraneContact` | `membrane_cortex_contact` |
| membrane↔cortex | `ErmTether` | `membrane_erm_cortex` |
| lamellipodium↔cortex | `CrosslinkConnector` | `lamellipodium_cortex_seam` |
| lamellipodium↔ecm | `SeriesJointConnector` | `lamellipodium_nascent_fa` |
| membrane↔ecm | `SubstratePlaneCoupling` | **none — see §4e** |

**Five couplings since 2026-07-31, and the fifth has no registry row.** `SubstratePlaneCoupling` is
the area-proportional adhesion energy and the non-penetration constraint §14.11 named as missing;
§4d is its design and §4e records that it is defined in the scenario module because this lane's
allocation forbids `aleph/vertical/**`. It is deliberately **not** given the name of a registry row,
for the reason `ecm_network.py` gives about `FarFieldClamp`: a class named after a registry row reads
as that row being implemented.

**Four couplings, and this table said five when it was written.** `membrane_ecm_contact` was in the
composition and has been taken out; the entry is corrected here rather than left to disagree with the
code, per the rule that a ledger section the code proves wrong is fixed and said so. The reason is in
§14.9 and it is a real property of that connector rather than a convenience: it is a **line-of-centres**
non-adhesive contact over a **fixed** site pairing, its own docstring records the line-of-centres
choice, and on a surface that spreads laterally the pairing goes stale — after a micron of spreading
the line joining a membrane node to the collagen point it was paired with is nearly horizontal, so the
non-penetration element pushes the membrane *sideways*. That force is plausible, quiet, and wrong. It
is left out and the penetration it would have prevented is **measured** instead, by
`lowest_membrane_height_um`, which is reported in every sample of every run.

**The load path that makes the cell spread**, stated once so it can be checked rather than assumed:

```
ecm ligand --[integrin clutch, series compliance]-- lamellipodial tip
            --[lamellipodium_cortex_seam]-- cortex --[contact + tether]-- membrane
```

Adhesion is the **only** thing in this world that lowers energy by putting cell surface near the
substrate. It is a real spring between two owners' material points, not a prescribed displacement.

**The three configurations:**

| Name | Difference | What it is for |
|---|---|---|
| `SYMMETRIC` | clutches on the whole basal rim, substrate clamped | area rises and saturates; the centre of mass descends |
| `POLARISED` | clutches on one half of the rim only, substrate clamped | the centre of mass also moves **laterally** — the single cell moves |
| `TREADMILL` | identical to `POLARISED` except `clamp_boundary=False`, and the free sheet is given the larger descent rate | **the control that must fail** |

`TREADMILL` is built against `POLARISED` and not against `SYMMETRIC`, because a treadmill is a
statement about **translocation** and symmetric spreading has none to lose. The free sheet is also
given the larger `mobility_um_per_pn_s`, which `relax.py` defines as a **descent rate and not a drag
coefficient** — that is exactly the right knob, because "the belt is the compliant body, the cell is
not" is a statement about which side gives first, and in a quasi-static descent that is what a
relative descent rate says. It changes no equilibrium; the landscape is identical and only the path
down it differs. `UNSOURCED`, applied in one mode only, and measured: without it the free sheet also
relaxes the cell, and the treadmill's centre of mass moves *more* than the anchored one's on the same
budget rather than less.

`ecm_network.build_planar_collagen_sheet` states the treadmill mechanism in its own docstring:
"Without it the sheet is a free body and every traction a cell exerts on it merely translates it."
That is the treadmill in one sentence, and it is why the control is built by removing the clamp
rather than by breaking the cell. Nothing about the cell's machinery differs between `SYMMETRIC` and
`TREADMILL`: same clutches, same stiffnesses, same initial separations, therefore the same tensions
at step zero. What differs is whether the traction has anything to pull against.

## 4d. The two elements added on 2026-07-31, and why they had to be new

The first landing of this scenario did not spread. §14.11 named two causes and neither was a tuning
problem, so this amendment adds the two missing elements. **This section is written before their
code**, and it is the design they are held to.

### 4d.1 An area-proportional adhesion energy

```text
E_adh = - W * sum_t  A_t * phi(h_t)
```

over the membrane's triangles, with `A_t` the triangle's area [µm²], `h_t` the height of its
centroid above the substrate plane [µm], and `W` the **work of adhesion** [pN/µm² = pN·µm/µm²].
`phi` is a C¹ plateau-then-smoothstep contact weight: exactly `1` for `h <= h_plateau`, exactly `0`
for `h >= h_range`, and `1 - 3u² + 2u³` in between with `u = (h - h_plateau)/(h_range - h_plateau)`.

Three properties, each of which is the reason for a design choice rather than decoration:

* **It is area-proportional, which is the whole point.** A ring of clutches rewards *pulling*; this
  rewards *covering*. `dE/dA_contact = -W` is a driving force for spreading that exists at every
  configuration, which is the thing §14.11 said was absent.
* **The plateau is exactly flat**, so `phi = 1` and `phi' = 0` throughout the contact zone. That
  makes `W` mean exactly what Young–Dupré needs it to mean — energy per unit area of adhered
  surface, with no residual dependence on how far above the plane the adhered surface happens to
  sit — and it puts the equilibrium standoff under the control of the non-penetration element alone.
* **The gradient is exact and analytic.** `dA_t/dx_a = (1/2) n_hat x (c - b)` and cyclically, and
  `dh_t/dz_v = 1/3`. Both are checked term-by-term against a central finite difference.

### 4d.2 A non-penetration constraint, reusing the element that already exists

`aleph/vertical/connectors.py` already ships compressive-only unilateral contact with an
**exact-zero** inactive branch and a `SOFT_CORE` law that diverges before the surfaces can cross —
`ContactLaw.SOFT_CORE` records, in its own docstring, that a linear unilateral contact presents a
barrier of only `f²/2k` and a monotone descent will punch through it one site at a time. The
substrate contact **reuses `UnilateralSpring` unchanged** rather than writing a second one. Nothing
about the law is re-derived here; what is new is only the geometry it is applied over.

The geometry is the one thing that makes this different from `membrane_ecm_contact`, which §4 threw
out. That connector pairs sites **by a fixed index** over a line of centres, and §14.9 records the
consequence: after a micron of spreading the line joining a membrane node to the collagen point it
was paired with is nearly horizontal, so the non-penetration element starts pushing the membrane
sideways. Here there is no pairing at all. The separation is the **signed height** `h = z - z_plane`
and the direction is exactly `+z` by construction, recomputed from the live configuration at every
evaluation. A surface that slides a micron sideways changes nothing about which direction is up.

Signed, and that is load-bearing: `relax.py` warns that under a line-of-centres law "penetration
reads as separation on the far side". A signed height cannot do that. `h <= 0` is **refused**
(`SubstrateContactError`) rather than evaluated, and `SOFT_CORE` diverges before a descent can get
there.

### 4d.3 Where the plane is, and what the substrate is idealised as

The plane is at the **mean height of the ECM's nodes**, read live from the position block handed to
the connector. So:

* the substrate is idealised as a **rigid plane carried by the sheet's mean height**. It is not the
  sheet's local surface. Stated as a limitation in §14.14, not buried;
* the reaction is real. `E` depends on the ECM only through `z_plane`, so
  `sum_i dE/dz_i^(membrane) = -dE/dz_plane`, and distributing `-dE/dz_plane` equally over the ECM's
  nodes closes force **exactly** — every component, to float64. `newton_pair=True` is a measured
  claim;
* `central_force=False` is declared, and that is honest rather than convenient. A rigid plane
  carries a couple: the load lands under the cell and the reaction is spread over the whole sheet,
  so the moment about the pair centroid does not vanish. The moment residual is **measured and
  reported** by `test_the_substrate_coupling_closes_force_exactly_and_declares_its_couple` rather
  than skipped silently;
* it keeps the treadmill real. A lab-frame plane would flatten the cell identically whether or not
  the dish was clamped, which would delete the discrimination this entry exists for. Because the
  plane follows the sheet, a **free** sheet given the larger descent rate is pulled up to the cell
  instead of the cell being pulled down to it — which is the treadmill, arriving as a consequence of
  the physics rather than as a special case in the code.

### 4d.4 What is switched off with the adhesion

`adhesion_engaged=False` now releases **every** adhesive element in the world: the clutches *and*
`W`. That is what the control it exists for needs — "with no adhesion of any kind, the cell does not
spread" — and splitting it into two flags would have let the negative control pass while half the
adhesion was still on. The non-penetration is **not** released by it: a cell with no adhesion still
does not sink through the dish, and a flag that removed the floor along with the glue would be
testing two things at once.

## 4e. This module now defines two force laws, and that is a finding

`test_the_scenario_defines_no_force_law` was a true statement about the first landing and is no
longer one. `substrate_adhesion_energy_and_forces` and `substrate_repulsion_energy_and_forces` are
force laws, they are named by the house convention so that the guard **sees** them rather than being
routed around, and the control is tightened to `test_the_scenario_defines_exactly_the_two_substrate_laws_it_declares`,
which asserts the set of `*_energy_and_forces` functions in this module is **exactly** those two. An
undeclared third one still fails it.

They are here because this lane's allocation is `aleph/scenarios/spread.py` and forbids
`aleph/vertical/**`, which is where an owner-published element belongs. **The proposal that they be
moved into an owner is `docs/decisions/PROPOSAL-substrate-adhesion-belongs-to-an-owner.md`**; it is
not this lane's decision to take.

## 4a. Mutation testing

Four mutants applied to `aleph/scenarios/spread.py`, the **whole** control file run against each,
then restored and verified byte-identical. Baseline: `25 passed in 305.30s`.

| Mutant | Result | Tests killed | Survived |
|---|---|---|---|
| **M1** — swap the scatter weights in `MaterialPointCoupling._scatter`: force still closes, the load lands at the other end of the segment | `8 failed, 17 passed` | 8, led by `test_the_adapter_conserves_force_and_moment` | no |
| **M2** — `clamp_substrate` returns `True` in every mode, so the treadmill's dish is anchored | `2 failed, 23 passed` | 2: `test_the_treadmill_is_distinguishable_from_spreading`, `test_traction_alone_cannot_tell_the_treadmill_apart` | no |
| **M3** — `MaterialPointCoupling.spring` returns `None`, so `CellWorld` treats the adhesion as purely dissipative and drops it from the descent | `6 failed, 19 passed` | 6 | no |
| **M4** — `centre_of_mass_um` returns a plain node mean | `1 failed, 24 passed` | 1: `test_the_centre_of_mass_is_area_weighted_and_not_a_node_mean` | no |

**All four killed.** M1 is the one worth naming: it conserves force exactly and is invisible to every
closure check, and what caught it is the moment assertion — plus, downstream, five physics controls,
because a load applied at the wrong end of a lamellipodial segment changes what the cell does.

**A measurement defect in the first harness, recorded because it is the exact trap the house brief
warns about.** The first run passed `-q` to `pytest`; `pyproject.toml` already sets
`addopts = "-q"`, so the run got `-qq` and the `N passed` summary line vanished entirely. The harness
then fell through to printing the last 300 characters, which happened to be two `FAILED` lines, and
reported "2 killed" for M1 — against the true 8. The numbers above are from a re-run without `-q`,
with each mutant's full output kept at `mut_M*.txt`.

## 4b. What was measured

Every number here is from `runs/spreading/`, produced by `python -m aleph.scenarios.spread`-style
driving at `chunks=8, steps_per_chunk=600`. The abscissa is **accepted descent steps**.

| | `SYMMETRIC` | `POLARISED` | `TREADMILL` | no adhesion |
|---|---|---|---|---|
| accepted steps | 4,255 | 4,243 | 4,219 | 2,971 (converged) |
| traction at step 0 [pN] | 3063.84 | **1044.12** | **1044.12** | 0.0 |
| contact area, first -> last [µm²] | 2.509 -> 5.629 | 2.509 -> 3.350 | 2.509 -> 2.513 | 2.509 -> 2.508 |
| footprint radius, last [µm] | 0.004 | **0.474** | 0.005 | 0.000 |
| centre of mass displacement [µm] | 0.411 | **0.198** | **0.0031** | 0.00003 |
| lateral displacement [µm] | 0.016 | **0.149** | **0.0024** | 0.000 |
| dish (perimeter) displacement [µm] | **0.0** | **0.0** | **1.407** | 0.0 |

The treadmill row is the result this entry exists for. `POLARISED` and `TREADMILL` start from the
**same** traction to the last bit, and end with the cell 65x further from where it started in one and
the dish 1.4 µm from where it started in the other. `POLARISED`'s contact area also does the thing
the classic assay curve does — it rises and **saturates**: 2.509, 2.755, 3.056, 3.227, 3.314, 3.350,
3.361, 3.359, 3.350.

## 5. Units, domains, singular cases, invariants

Length µm, force pN, energy pN·µm. **There is no time axis and the module must not pretend there is
one.** `aleph/vertical/relax.py` states plainly that its `dt` is a step size and not a time and that
its mobility is a descent rate and not a drag coefficient. So the abscissa of every series produced
here is **accepted descent steps**, spelled `step` in the artefacts and never `t` or `time_s`, and
`test_the_series_abscissa_is_steps_and_not_seconds` is the control that pins it.

**Singular cases, each refused rather than clamped:**

| Case | Response |
|---|---|
| a clutch whose two ends are coincident | `AdhesionCouplingError` from the joint — the axis is undefined |
| a lamellipodial site outside its filament | `ValueError` from the owner |
| an ECM ligand outside its fibre | `ValueError` from the owner |
| a cell placed so no clutch can reach the sheet | `SpreadingConfigError` — a world with no adhesion spreads by nothing and would report a flat curve that looks like a physics result |
| a substrate with neither clamp nor far-field | permitted **only** in `TREADMILL`, and named as such |
| a contact plane taken from the ECM's live position | refused; the footprint plane is the lab frame |

**Invariants:**

* **I1 — the footprint plane is fixed in the lab frame.** Measuring contact against the *current*
  ECM position is exactly the mistake that hides a treadmill: if the substrate rises to meet the
  cell, a relative measurement reports contact that the dish never saw. Both are computed;
  `contact_area_um2` takes the lab plane and `relative` is reported beside it as a diagnostic.
* **I2 — the centre of mass is vertex-area-weighted.** A plain node mean over an icosphere is
  *nearly* the centroid, and "nearly" drifts as the mesh deforms anisotropically — which is precisely
  what spreading does to it.
* **I3 — force and moment close through the adapter.** `MaterialPointCoupling` gathers with
  `(1-w, w)` and scatters with the same `(1-w, w)`, so the applied point is unchanged and the moment
  about any origin is preserved. Force closure alone permits a load to land at the wrong place.
* **I4 — traction is a sum of magnitudes; translocation needs a resultant.** They are reported
  separately and never conflated. `PLAN.md` §2.5 records a resultant-based tolerance rejecting 40,000
  consecutive steps of a correct relaxation for the mirror-image reason.

## 6. Source evidence class and known retractions

Source evidence class: **none** — nothing was read. Aleph evidence rung: `ANALYTIC_ORACLE` for the
composition (it reproduces owners whose gradients are checked), `BLOCKED` for any quantitative
claim about a real cell. See §14.

## 7. Independent oracle or derivation

Three oracles that do not depend on any physical constant in this scenario:

1. **Adhesion off ⇒ no spreading.** With every clutch released (`BindingState.UNBOUND`), the series
   stiffness is **exactly** `0.0` by the owner's own gate, so the world reduces to the verified
   vertical plus a non-adhesive contact. The contact area must not grow. This separates "the cell
   spread" from "the cell fell", which no amount of looking at a rising curve does.
2. **Rigid translation invariance.** Translating the cell and the sheet together by a constant vector
   changes no energy and no force. A footprint measured against a hard-coded `z = 0` would fail this
   and would be measuring the coordinate system.
3. **The treadmill partition.** Cell displacement plus substrate displacement is a conserved
   statement about who moved; with the sheet clamped the substrate's share is exactly zero, and with
   it free the substrate's share is the majority. The ratio is the oracle, not either number.

## 7a. The independent oracle for the adhesion energy: Young–Dupré

**This is the check with no free constant in it, and it is the reason §4d.1 is an energy rather than
a force.** At equilibrium a drop of surface tension `sigma` wetting a solid satisfies Young's
relation `gamma_SV = gamma_SL + sigma cos(theta)`, and Dupré's definition of the work of adhesion
`W = gamma_SV + sigma - gamma_SL` gives

```text
W = sigma * (1 + cos theta)
```

with `theta` the contact angle measured **inside** the cell. The two limits are what pin the
convention and both must work:

| Limit | Prediction | Why it is the right test |
|---|---|---|
| `W -> 0` | `cos theta -> -1`, `theta -> 180°` | a sphere touching the dish at a point |
| `W >= 2 sigma` | `cos theta -> 1`, `theta -> 0°` | complete wetting; the drop spreads without bound |

**A note on the sign convention, because the instruction this work was done under wrote
`W = sigma(1 - cos theta)` and then gave those two limits.** The two statements are not consistent
with each other: `sigma(1 - cos theta)` sends `W -> 0` to `theta -> 0`, which is complete wetting and
not a point contact. The **limits** are unambiguous physics and they select `W = sigma(1 + cos theta)`,
which is also the standard form. `theta` here is the angle inside the cell; the other convention is
the same relation read with `theta' = 180° - theta`. Both numbers are reported by
`WettingMeasurement` so a reader can check rather than having to trust this paragraph.

**How `theta` is measured, without using `W`.** Not by fitting a circle and not by finding the
contact line, both of which are unreliable on a 162-node icosphere. For a spherical cap of height
`H` and volume `V`,

```text
V = (pi/3) H^3 (2 + cos theta)/(1 - cos theta)      =>      cos theta = (k - 2)/(k + 1),  k = 3V/(pi H^3)
```

so the contact angle follows from two quantities a coarse mesh reports accurately: the **enclosed
volume**, which is exact for a closed triangulation and which `membrane.enclosed_volume_um3` already
computes, and the **apex height above the substrate plane**. Sanity: a sphere resting on a point has
`H = 2R` and `V = pi H^3/6`, so `k = 1/2` and `cos theta = -1`; a hemisphere has `H = R` and
`V = (2/3) pi H^3`, so `k = 2` and `cos theta = 0`. `k < 1/2` is **refused** — it names a shape with
less volume than any cap of that height, which is not a cap at all.

`sigma` is measured too, in the sense that matters: the oracle world is built with **no cortex**, so
the only surface energy in it is the membrane's own `sigma * A` term (`membrane.MembraneEnergy`
carries it as exactly that), and `sigma` is the material constant that term is built from. The
recovered `W` is then compared with the `W` that was put into the energy, and neither number was used
to obtain the other.

**The oracle world is a sessile drop and nothing else**: a Helfrich membrane, an osmotic envelope in
`VOLUME_PENALTY` so the volume is fixed, a clamped collagen sheet, and the one substrate coupling.
No cortex, no lamellipodium, no clutches, no ECM ligands. It runs the **shipped** connector, not a
copy, which is the only reason a result on it says anything about the scenario.

## 7b. What Young–Dupré measured

Every row is from `runs/spreading/young_dupre.json`, produced by the **shipped**
`SubstratePlaneCoupling` on the sessile-drop world at `sigma = 3.0` pN/µm, level-1 icosphere
(42 nodes), 8 chunks of 600 descent steps. `W_in` appears in **none** of `theta`, `H` or `A_wet`.

| `W` in | predicted `theta` | measured `theta` | `W` recovered | error |
|---|---|---|---|---|
| 0.0 | 180.0° | 169.97° | 0.046 | the `W -> 0` limit; see below |
| 1.5 | 120.0° | 120.12° | 1.4945 | **−0.4%** |
| 3.0 | 90.0° | 87.18° | 3.1477 | **+4.9%** |
| 4.5 | 60.0° | 63.84° | 4.3225 | **−3.9%** |
| 6.0 | 0.0° (complete) | 45.96° | 5.0855 | −15.2%; see below |

**The `W -> 0` limit** converged in 109 accepted steps at 169.97°, covering 1.67 µm² of dish. That is
a sphere touching at a point, to the accuracy a 42-node mesh resting on a 0.15 µm non-penetration
standoff can represent one; a true point contact is not representable on a mesh and the residual
10° is that standoff, not a defect in the energy.

**The `W >= 2 sigma` limit is complete wetting, which means spreading without bound**, so it has no
equilibrium angle to be measured at and the −15.2% row is not an error in the recovery — it is a
drop that has not stopped. What it does have is the right *behaviour*: it covers 117.75 µm² of dish,
more than the drop's entire initial surface area of 105 µm², against 22.87 µm² for the partially
wetting `W = 1.5` case on the same budget. `test_the_complete_wetting_limit_spreads_far_further_than_the_partial_one`
asserts the comparison rather than a number, because pinning a number there would be pinning the
step budget.

**Convergence.** The middle rows are not converged either, and the bias is a descent that is still
running rather than a systematic error. Measured at level 2, `W = 3.0`, by running the same drop for
four budgets:

| accepted steps | `theta` | `W` recovered | error |
|---|---|---|---|
| 7,037 | 93.90° | 2.796 | −6.8% |
| 17,586 | 88.81° | 3.063 | +2.1% |
| 35,168 | 88.44° | 3.082 | +2.7% |
| 63,300 | 88.09° | 3.100 | **+3.3%** |

So the recovery settles at about **+3%** on a 162-node mesh, and the residual is a discretisation
bias rather than an unfinished relaxation. Nothing was tuned to produce these numbers: `sigma` is the
membrane's material card, `theta` is two extrema of the vertex set, and `W` is what the energy was
built from.

## 8. Positive control

`tests/scenarios/test_spread.py`:

* `test_the_contact_area_grows_and_the_centre_of_mass_descends`
* `test_the_polarised_configuration_moves_the_centre_of_mass_laterally`
* `test_the_composition_holds_the_six_owners_and_five_couplings_it_claims`
* `test_the_adhesion_load_path_is_closed_at_construction`
* `test_the_coverage_verdict_is_complete_after_the_scenario_runs`
* `test_the_adapter_conserves_force_and_moment`
* `test_the_centre_of_mass_is_area_weighted_and_not_a_node_mean`
* `test_traction_is_reported_as_two_numbers_and_never_one`
* `test_the_series_abscissa_is_steps_and_not_seconds`
* ~~`test_the_scenario_defines_no_force_law`~~ — replaced by
  `test_the_scenario_defines_exactly_the_two_substrate_laws_it_declares`; the ownership debt is
  explicit and bounded while `PROPOSAL-substrate-adhesion-belongs-to-an-owner.md` is open
* `test_the_frames_are_written_and_are_real_pngs`
* ~~`test_the_symmetric_configuration_lifts_its_own_centre`~~ — **replaced 2026-07-31** by
  `test_the_symmetric_configuration_goes_from_round_toward_pancake`. The old control pinned the
  §14.11 shortfall and its own docstring instructed that it be replaced by the spreading assertion
  it was standing in for once a support element was added. One was; it is.

Added 2026-07-31, for the two substrate elements (§4d):

* `test_the_area_adhesion_gradient_matches_a_central_difference` — per term, order > 1.6, energy
  asserted non-zero **first**, because a term returning zero energy and zero force agrees with its
  own finite difference perfectly
* `test_the_non_penetration_gradient_matches_a_central_difference`
* `test_the_contact_weight_is_exactly_one_then_exactly_zero` — the plateau and the cut-off are the
  numbers, not small numbers
* `test_force_closes_exactly_and_the_couple_is_measured_rather_than_claimed_away`
* `test_translating_the_whole_scene_changes_no_substrate_energy_or_force`
* `test_covered_area_is_not_grown_by_crumpling_the_membrane` — §14.15
* `test_the_symmetric_configuration_goes_from_round_toward_pancake`
* `test_the_picture_has_a_control_and_the_footprint_survives_it` — §12
* the Young–Dupré battery: `test_the_contact_angle_of_an_undeformed_sphere_is_exactly_180_degrees`,
  `test_the_contact_angle_cannot_be_told_the_work_of_adhesion`,
  `test_young_dupre_recovers_the_work_of_adhesion_the_measurement_never_saw`,
  `test_the_zero_adhesion_limit_is_a_sphere_touching_the_dish`,
  `test_the_complete_wetting_limit_spreads_far_further_than_the_partial_one`,
  `test_the_oracle_world_has_no_cortex_to_hide_a_tension_in`

## 9. Deliberately failing negative control

`tests/scenarios/test_spread.py`:

* `test_the_treadmill_is_distinguishable_from_spreading` — **the control the PI asked for.** Clutch
  tension of the same order in both; lab-frame footprint growth and lab-frame centre-of-mass
  displacement an order of magnitude apart; and the substrate's own displacement the opposite way
  round. If this test cannot separate them, the scenario has not shown what was asked.
* `test_with_every_clutch_released_the_cell_does_not_spread` — adhesion is what spreads it.
* `test_a_world_whose_clutches_cannot_reach_the_sheet_is_refused`
* `test_traction_alone_cannot_tell_the_treadmill_apart` — the *measurement* negative control. It
  asserts that the step-zero traction is **identical** in the two worlds and that it decays
  monotonically in both, so a traction-first reading reports two cells doing the same thing and one
  of them never went anywhere. Without it, putting contact area and centre of mass ahead of traction
  would look like a preference rather than a requirement.
* `test_a_swapped_scatter_weight_conserves_force_and_moves_the_load` — the negative control for the
  adapter's moment test, so that its tolerance is not being carried by force closure.
* `test_an_interpolation_weight_outside_the_segment_is_refused`
* `test_translating_the_whole_scene_does_not_change_the_footprint`

Added 2026-07-31, each driving a **shipped** deliberate-break flag on `SubstratePlaneCoupling` and
asserting in **both** directions — the broken branch misbehaves and the healthy branch is exactly
inert where it should be:

* `test_withholding_the_substrate_reaction_breaks_force_closure` — without it `newton_pair=True` is
  an unchecked claim
* `test_omitting_the_area_adhesion_leaves_a_contact_that_cannot_spread` — the missing-term control.
  The world still assembles, still relaxes, still conserves energy perfectly
* `test_omitting_the_non_penetration_removes_the_floor_the_cell_flattens_against` — asserted on the
  **difference** between the two worlds rather than on the total force, because at the initial
  configuration the adhesion pulls the basal nodes down harder than the contact pushes them up, so a
  test on the sign of the total would have been a test on which of two terms was larger
* `test_flipping_the_adhesion_sign_pushes_the_cell_off_the_dish`
* `test_a_coupling_with_both_branches_switched_off_is_refused`
* `test_a_surface_at_or_below_the_substrate_plane_is_refused`,
  `test_a_negative_work_of_adhesion_is_refused`, `test_a_zero_width_adhesion_ramp_is_refused`,
  `test_a_plateau_below_the_contact_standoff_is_refused_at_composition_time`,
  `test_a_prolate_body_has_no_contact_angle`

## 10. Numerical and precision envelope

float64 throughout. The descent is explicit and the clutch is the stiffest element, so the step size
is set by it and is discovered by rejection. Relaxation is run in **chunks** so observables can be
sampled; each chunk is a real `relax_to_equilibrium` call over the same world, and the chunk boundary
is a sampling point and not a physical event. Convergence within a chunk is legitimate and is
recorded per chunk rather than assumed.

## 11. Production-backend residency and transfer

CPU / `NumpyBackend` only. **No GPU work, and this lane holds no authorization.**

## 12. Visualisation

Required by the PI's instruction quoted in §0 and therefore part of this entry, not an extra.

**matplotlib is not a dependency of this project and is not installed in the `aleph` environment** —
`pyproject.toml` line 10 declares `dependencies = ["numpy>=1.26", "scipy>=1.11"]`, and
`import matplotlib` raises `ModuleNotFoundError` under `/Users/sw1/miniconda3/envs/aleph/bin/python`.
Reported here as a finding rather than worked around by adding a dependency this lane has no
authority to add.

So the renderer is **pure standard library plus numpy**: it rasterises into an RGB array and writes
a PNG with `zlib` and `binascii.crc32`. It has no import that can fail, so frames are produced rather
than skipped. Each frame is a two-panel figure — a side elevation (x–z) and a plan view (x–y) — with
the substrate fibres, the membrane surface, the lamellipodial filaments, the engaged clutches, and
the measured footprint drawn. Frames are also tiled into one contact sheet so the whole time series
can be judged in a single image.

`aleph/viz/` was read and **not** modified: `scene.py`, `frames.py`, `capture.py` and
`publisher.py` are a streaming capture/replay transport for a live buffer feed, not a still
renderer for a composed scenario, and there is no entry point in them that takes a `CellWorld`.

## 13. Acceptance

Accepted when: the six owners and five couplings are all present and the coverage verdict is
complete; the contact area rises and saturates and the centre of mass moves; the treadmill is
separated by the assertions in §9; the mutants are killed; and the frames are written and looked at.

**Measured after the substrate-law amendment:** `tests/scenarios/test_spread.py` collects 48
controls, and the full 5,032-test repository run on 2026-08-03 passed all 48. The cell now goes from
round toward pancake, §14.10 and §14.11 are discharged, the Young–Dupré reduced-world controls pass,
and the no-adhesion and treadmill branches remain discriminating negative controls.

Status remains `PROPOSED`, not `ACCEPTED`, for three current reasons rather than the superseded
§14.11 shortfall: the four published finite-budget runs do not meet their final residual threshold;
the area adhesion and plane contact still live in the scenario package pending
`PROPOSAL-substrate-adhesion-belongs-to-an-owner.md`; and all material magnitudes remain
`UNSOURCED`. The encoding is accepted at `ANALYTIC_ORACLE`; no quantitative real-cell claim is.

## 14. Honest limits — what this entry does NOT establish

Stated as absent rather than left to be inferred.

1. **There is no time in this scenario.** The abscissa is descent steps. Nothing here is evidence
   about spreading *rate*, about the ~minutes timescale of real spreading, or about any comparison
   to an experimental area-versus-time curve. `relax.py` says why in its own first paragraph.
2. **No rate constant is supplied, invented, or defaulted.** `lamellipodium.rate_prior()` refuses
   every one of the four it owes evidence on, and this scenario does not route around it. Where a
   protrusion increment is applied it is a **prescribed kinematic input to the scenario**, declared
   on the config field, labelled `UNSOURCED`, and reported in the run artefact. A cell that spreads
   because it was told how fast to spread is a demonstration, and §4 says which of the two this is:
   **the spreading here is adhesion-energy-driven and needs no prescribed velocity at all.** The
   prescribed-elongation path exists as an option and is **off** in all three shipped
   configurations, so that no result in `runs/spreading/` depends on it.
3. **Clutch kinetics are absent, not implied.** No clutch binds or releases during a run. Binding
   would need on/off rates and a force-dependent lifetime, none of which this project holds. So the
   number of engaged adhesions is constant, and nothing here is evidence about adhesion turnover,
   about the maturation sequence, or about the catch-bond.
4. **No actin flow.** There is no retrograde flow field in this scenario, so the traction reported
   is a clutch tension and not a flow-coupled traction. The treadmill control is therefore a control
   on the *observables* — footprint and centre of mass — and not a reproduction of a treadmilling
   cell's internal kinematics.
5. **`nmii` is not in the world.** Contractility is absent, so nothing here is evidence about the
   contractile ring, about stress-fibre tension during spreading, or about the retraction phase that
   follows spreading in a real cell. Absent, not implied.
6. **The substrate is a woven sheet with ~10² nodes**, its stiffness is the collagen card's default,
   and no stiffness sweep was run. Nothing here is evidence about rigidity sensing.
7. **The membrane is an icosphere at subdivision level 1 or 2.** `PLAN.md` §2.5 records that step
   count depends on tessellation as much as on resolution. No mesh-family invariance study was run
   for *this* scenario; the one in `runs/mesh_family_invariance` is about the vertical.
8. **A converged descent is a local minimum near where it started.** `relax.py` says it has nothing
   to say about any other one.
10. ~~**There is no non-penetration element between the membrane and the substrate.**~~
   **Discharged 2026-07-31.** There is one now, and §4d.2 records that it is `UnilateralSpring`
   under `SOFT_CORE` applied along the plane normal rather than a second implementation of the same
   idea. The reason `membrane_ecm_contact` was thrown out — a fixed site pairing whose line of
   centres turns horizontal after a micron of spreading — does not apply to it: there is no pairing,
   the separation is the signed height, and the line of action is `+z` by construction.
   `lowest_membrane_height_um` is still reported in every sample, and it is now a **check** on the
   contact rather than a measurement of an absence.
11. ~~**The cell does not flatten, and the projected outline area falls rather than rises.**~~
   **Discharged 2026-07-31, and this is what the amendment was for.** The old text is kept below
   because it is the measurement the fix is judged against, not because it is still true.

   > Measured: `spread_area_um2` goes 26.16 -> 24.60 µm² over 1,353 accepted steps in `SYMMETRIC`
   > and 26.16 -> 25.67 in `POLARISED`. [...] Two causes, both named rather than guessed at: the
   > adhesion is a **ring of discrete clutches at the rim** and not an area-proportional adhesion
   > energy over the basal surface, and there is **no non-penetration element** [...] so a downward
   > pull descends the cell instead of flattening it against something.

   Both elements were added (§4d). Measured now, `SYMMETRIC`, same 3 x 500-step budget:

   | | before | after |
   |---|---|---|
   | contact angle | 180.0° -> 180.0° (never a cap) | **180.0° -> 116.7°** |
   | projected outline `spread_area_um2` [µm²] | 26.16 -> **24.60** | 26.16 -> **39.59** |
   | proximity footprint `contact_area_um2` [µm²] | 2.509 -> 5.629 | 2.509 -> **13.196** |
   | covered area `wetted_area_um2` [µm²] | — | 2.53 -> **39.58** |
   | lowest membrane node [µm] | 0.100 -> 0.247 (**lifting away**) | 0.100 -> **0.15** (on the standoff) |

   `test_the_symmetric_configuration_lifts_its_own_centre` has been **replaced** by
   `test_the_symmetric_configuration_goes_from_round_toward_pancake`, on the old test's own written
   instruction: *"If that is because a support element was added, this test has done its job and
   should be replaced by the spreading assertion it was standing in for."*
14. **The substrate is a rigid plane at the ECM's mean height, not the sheet's local surface.**
   §4d.3 gives the design; this is what it costs. Vertical force closes **exactly** — the identity
   is the chain rule, not an assertion — but the reaction is distributed uniformly over the dish
   while the load lands under the cell, so the pair carries a **couple**. `central_force=False` is
   declared on the `AdjointPair` for that reason, which means the ledger does not check moment
   closure on this connector; the residual is measured instead by
   `test_force_closes_exactly_and_the_couple_is_measured_rather_than_claimed_away`, which asserts
   the one part of it that is *not* an idealisation — that there is no torque about the plane
   normal. A local surface coupling against `CollagenNetwork`'s own geometry would remove the
   couple and reintroduce the staleness problem §14.9 records; which of the two is wanted is a
   modelling decision and is in
   `docs/decisions/PROPOSAL-substrate-adhesion-belongs-to-an-owner.md`.
15. **The adhesion energy is proportional to the area of substrate *covered*, not to membrane area.**
   The distinction is not cosmetic and it was found the expensive way. With the triangle's own area
   the scenario does not spread: it rewards making membrane area near the dish, and the cheapest way
   for a coarse mesh to do that is to **crumple**. Measured over 1,790 accepted steps at `W = 6`:
   membrane area rose 105.3 -> 111.7 µm² while the projected outline **fell** 26.16 -> 23.64 µm² and
   the greatest in-plane radius fell 3.012 -> 2.825 µm. `test_covered_area_is_not_grown_by_crumpling_the_membrane`
   pins it.
16. **The Young–Dupré recovery is checked on a reduced world, and it does not transfer numerically
   to the composed cell.** The oracle drop has no cortex, so its `sigma` is the membrane's material
   constant. The composed cell carries an elastic cortex whose rest lengths never turn over, and the
   effective tension the adhesion works against is several times `sigma` — which is why the shipped
   `W` is 20 pN/µm² while `2 sigma` is 6. **No effective tension was fitted and none is claimed**:
   what the scenario reports is the contact angle it reached, and §7a's recovery is a statement
   about the *element*, not about the composed cell.
12. **The proximity footprint is quantised by the mesh.** On a level-1 icosphere of radius 3 µm the
   node spacing is ~1.6 µm and adjacent basal nodes differ in height by ~0.43 µm, so a 0.2 µm
   contact band admits at most the single lowest node and the footprint jumps between one node's
   lumped area and exactly zero. At level 2 it oscillates between 1.15 µm² and 0 on consecutive
   samples. Read it with `lowest_membrane_height_um` beside it, never alone.
13. **The membrane tension and the two cortex moduli are below `build_vertical`'s defaults and are
   `UNSOURCED`.** They are a scenario choice, recorded on the config fields with the measurement
   that forced them: at the vertical's 30 pN/µm the footprint stays exactly `0.0000` µm² for 1,581
   accepted steps. A low tension stands in for the membrane reservoir a spreading cell draws area
   from, and this model has no reservoir. Nothing here is evidence about tension during spreading.
