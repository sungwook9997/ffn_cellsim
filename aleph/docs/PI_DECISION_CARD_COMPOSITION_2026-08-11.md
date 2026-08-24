# PI decision card — everything that blocks the composed cell, in one place

**Why this exists.** Wiring one connector stops at one determination, every time. On 2026-08-11 that
happened four times in a row inside a single session: a mobility, a head-exclusivity rule, a bond
selection rule, an α-actinin rate law. Each stop was correct — the rails did their job — but
discovering them one at a time makes every step feel like a new wall.

So this enumerates them **all**, from the code that actually refuses, not from memory. Decide them
together and the wiring after is mechanical.

**Read the split first.** Three of these are not PI decisions at all, and one of them is what stopped
the session's last attempt. Sorting them was the point of writing this down.

---

## A. NOT a PI decision — implementation, data already sourced

These need code. Surfacing them to PI would waste a decision on a value that exists.

| # | What | Every constant it needs | Blocks |
|---|---|---|---|
| A1 | **α-actinin bind/unbind delegate** | `ALPHA_ACTININ`: `k_on` 10 /s, `k_off0` 0.066 /s, Bell `f0`, `capture_radius_um` 0.06 µm, `link_k` 4.6e5 pN/µm — Ferrer 2008, PI-approved 2026-06-30 | `sf_cortex_transient`. Topology + specs already landed (`15f86c9c`); the builder refuses to bind without kinetics because a transient crosslink that cannot unbind is a permanent weld. |
| A2 | **Per-connector force accumulators** | none | Attribution on the one cut that carries two power ports (`cortex`/`membrane`). Today a shared cut scores both edges jointly and `attributable_cuts()` honestly excludes it. |
| A3 | **`ac_composed_world_dump.py` raises on the two new mechanical groups** | none | A dump path. Standing since 2026-08-10 with **no lane owning it**. |

**A1 is the immediate one.** It is the whole distance between "one connector wired" and "two".

---

## B. PI must choose a RULE — the datum exists, the modelling choice does not

| # | Decision | Why it cannot be defaulted | Blocks |
|---|---|---|---|
| B1 | **NMII head exclusivity** — may one head bind both a cortex and an SF target? | Each MOTOR connector owns its own `bound_d` over the *same* 8,840 heads. Binding a head twice double-counts its force — the co-location trap the charter names. Options: (a) a shared claim array, (b) partition minifilaments geometrically, (c) one motor edge at a time. | `nmii_sf_motor`, and every further MOTOR edge (`nmii_lamellipodium_motor`, `nmii_filopodium_motor`) |
| B2 | **ECM candidate → bond selection** | `mikado_topology` produces `SegmentContactCandidates` and states the separation deliberately: *"candidates are not crosslink bonds… this prevents a topology builder from silently creating static chemistry."* Something must select and associate them under accepted-step kinetics. `ecm_library` carries `k_xl_pN_um` and `xl_contact_um`, so the chemistry is there; the **selection rule** is not. | `ecm_crosslink` |
| B3 | **Per-component overdamped mobility γ** | Engine components own `position_d` and `force_d` and have **nothing that integrates one from the other** — only the frozen incumbent driver has an integrator. Any per-component relaxation needs γ, and γ is a derived physical quantity per component, not a wiring detail. `laws/relax.py` derives one for a `FiberNetwork` (`U.fiber_point_drag`); the other components have none. | every component owning private arrays; the partitioned coupled solve |
| B4 | **Membrane resolution** — is subdiv 8 required, or does the bending length justify 7? | subdiv 6/7/8 costs 6.6 / 24.2 / **93.8 ms** per inner iteration, and the membrane is ~71% of the step at native. "native" means physiological density × geometry, not a subdivision number someone typed; the right criterion is how many points resolve ℓ = √(κ/σ). | up to 4× of the step cost, and the whole speed budget |

---

## C. PI must source a DATUM, or ratify the gap — the code already refuses

Every one of these raises rather than inventing a value. They are honest blocks, not bugs.

| # | Missing datum | Where it refuses | Blocks |
|---|---|---|---|
| C1 | **LINC nesprin/plectin on/off-rate cards** | `SourceGatedLincKinetics.propose_candidate` — *"no validated on/off-rate card exists… Surface to PI before binding"* | `actin_cap_linc`, `mt_nucleus_linc`, `if_nucleus_linc` (3 edges) |
| C2 | **Keratin K8/K18 `EA` / `x_max`** | `intermediate_filament_rig` — paywalled source, `NotImplementedError`, no default | every `intermediate_filament` edge (3) |
| C3 | **MT catastrophe / rescue Poisson rates** | `microtubule_rig` — UNSET PI-GAP slot | every `microtubule` edge (4) |
| C4 | **Fascin crosslink stiffness** | `protrusion.FASCIN_K_CROSSLINK_PN_PER_UM = None` — *"NOT FOUND — PI-GAP slot, no default"* | `filopodium` edges (5) |
| C5 | **`mu_medium`** — medium viscosity at 37 °C | `medium_exterior` — REQUIRED with no default; the KB's 14 viscosity rows are all intracellular or membrane | `membrane_medium_traction`; already on the dashboard |
| C6 | **SF inventory** — stress fibers per cell, NMII per SF | no KB axis exists (TAG returns 0 rows) | how many `sf_cortex_transient` joints is physiological. **Today's 15 is placement-derived, not sourced.** |

---

## What each answer buys, measured

The composed native world binds **3 of 14 components** and **2 of 38 connectors**. The 2026-08-11
interface residual measured the consequence: where a connector IS wired the adjoint closes to
round-off (95.0 pN present, cancellation 2.47e-15), and the two other measured cuts carry **47.5 pN
with nothing crossing them** — absent connectors, not a solver failure.

**33 of the 37 unbound edges are blocked by a missing COMPONENT, not by their own connector code.**
Component coverage is the lever: `cytosol` alone gates 9 edges, `membrane` 7, `focal_adhesion` 6,
`nucleus` 5, `lamellipodium` 5, `filopodium` 5, `microtubule` 4, `intermediate_filament` 3.

And the interior fluid column — `cytosol` + `membrane` + `nucleus`, the largest single block — is the
**Card-5 seam PI already approved on 2026-08-11**. That approval is not yet spent.

## The acceptance test — ⚠ CORRECTED 2026-08-11 (the first version was wrong)

**The criterion this card first carried was broken, and it was already failing on a correctly wired
connector when it was written.** It said: *a newly bound connector must move its cut's cancellation
from ~1 to ~0.* It cannot, because `measure_interface_residual` reduced **whole-body** force arrays —
so a component sitting on two edges carries BOTH into every cut it appears in, and the other edge's
force shows up as an uncancelled residual.

`sf_cortex_transient` was bound, accumulating, and conserving, and its cut still read
`cancellation = 0.99999999998`. `cortex|nmii` passed only because it was, that day, the only edge
carrying force. Every edge landed from here sits on a body with more than one edge, so the original
criterion would have declared **every future correct wiring a failure**.

The arithmetic that settled it: `47.5077854660593 = 95.0155709321186 / 2` to 2e-15 relative. The
"unpaired" force was the motor's per-side half, inherited by every cut touching cortex or nmii.

**The corrected criterion — LEAVE-ONE-IN.** Zero the watched force arrays, re-run ONLY the connectors
crossing that cut, then score. The reading then belongs to the connector rather than to the bodies.
Records carry `isolated: true` when taken this way; a non-isolated reading is a property of the BODIES
and may not be quoted as a connector result. `reaction_pn` and `traction_pn` are recorded separately,
because a residual and a scale alone cannot say which side carries the load — and inferring it by
arithmetic is exactly how the wrong root cause above was reached.

Measured under the corrected criterion (job 98, `outputs/ac/interface_residual/composed_native_isolated.json`),
**both bound power ports now score and both adjoints close**:

| cut | A [pN] | B [pN] | residual [pN] | cancellation |
|---|---|---|---|---|
| `cortex\|nmii` (`nmii_cortex_motor`) | 4.7508e+01 | 4.7508e+01 | 1.656e-13 | **1.74e-15** |
| `cortex\|sf_arc` (`sf_cortex_transient`) | 4.5140e-10 | 4.5140e-10 | 5.780e-26 | **6.40e-17** |

The motor's two halves are 47.5077854 each — exactly `95.0155709 / 2`, which is what the whole-body
reading had been smearing across every cut touching either body. Isolation recovers it as a pair.

**Scope.** `sf_cortex_transient`'s magnitude is small because its rest length is the separation at
bind, so the joints are unstrained. Both rows show the scatter CONSERVES; neither shows the edge
carries a physically meaningful load, and no magnitude here is quotable (`k_axial`, `k_xb` are PI-GAPs
and all 30 steps are force-accepted).

Not reducible to the fan-out's first suggestion (reduce `joint_force_on_a` twice): that is vacuous,
since `_segment_joint_force_kernel` scatters `+w·f` to A and `−w·f` to B and cancellation is guaranteed
by construction. Leave-one-in asks the question that can fail — did this connector put force into
**both** arrays.

⚠ **Two edge classes still have NO acceptance test**: `SERIES_HALF` (all six FA edges — `balance_cuts()`
has no key containing `focal_adhesion`) and `RESERVOIR_PORT` (the far side does not move, so nothing
cancels). Wiring either and choosing the evidence afterwards is the shape the charter forbids. **New
PI item B5 below.**

## Recommended order

1. **A1** — no decision needed, and it is the difference between 2 wired connectors and 3.
2. **B3 + the Card-5 approval** — the coupled solve needs γ, and the interior column is the biggest
   block and is already approved.
3. **B1, B2** — the two remaining rules, each unblocking a family of edges.
4. **C1–C6** — source or ratify as gaps; each is a family that stays honestly blocked until then.
5. **B4** — pure speed, and it should be argued from the bending length rather than from cost.

---

## ⚠ Added 2026-08-11 after a 12-agent audit — items the first draft missed

**New PI rule (B5): what closes a `SERIES_HALF` or `RESERVOIR_PORT` wiring?** The cancellation gate
structurally cannot see either. Either ratify leave-one-in per connector as the instrument for series
edges BEFORE the FA series is wired, or state that such a landing is reported "bound, unscoreable".

**New PI data (C7–C10), none of which appear above:**

| # | Missing datum | Where it refuses |
|---|---|---|
| C7 | `k_hub_pn_per_um` — the MTOC anchor stiffness. **No constant, no derivation and no PI-GAP declaration exists anywhere in the tree**; the only numbers are a test literal and an archived substitution. Harder than the MT catastrophe/rescue rates, which are `has_events=False` and not needed to bind. | `microtubule_rig.py:1053-1054` |
| C8 | The LINC stiffness / rest / stiffening magnitudes — separate from C1's rate cards. | `linc_connector.py:770-772`, `:1394-1400` |
| C9 | May `INTEGRIN_A5B1` (α5β1–fibronectin, Kong 2009) stand as proxy for the **α2β1–collagen** clutch the architecture requires — and if so, ratified *with* the unreconciled 7-vs-30 pN F\* discrepancy? | `hand_kmc.py:150,163`; `fa_clutch_warp.py:13-17` |
| C10 | Π₀ — three live values, none MCF7 (40 / 133 / ~72 Pa, a ~3.3× spread on the resting baseline), plus `L_p`. | `assemble.py:140-155` |

**Correction to C2 above:** keratin `EA`/`x_max` blocks the WLC regime and the three IF *edges*. It does
**not** block binding the `intermediate_filament` component, which this card implied.

**Three latent defects found, none yet fired — all class A, no decision needed:**

1. **The composed ECM is an unanchored free-floating volume.** The caller asserts
   `BoundaryAnchorMode.FAR_FIELD_DIRICHLET` (`ac_gate_b_composed_native.py:198`) and `_ecm_run` never
   applies the anchor. `ECMWorld.accumulate_mechanics` is its only caller and is constructed only in
   tests. A live physics defect in a running configuration.
2. **Binding LINC today would be a silent permanent weld.** `LincRigConnectorAdapter` names its hook
   `propose_candidate` (`linc_connector.py:725-727`); `CellTransaction` scans for `propose_events`
   (`transaction.py:78,141`). The proposal would be skipped and `SourceGatedLincKinetics`' refusal
   would **never fire** — the charter's forbidden case, passing quietly. Fix before anyone touches C1.
3. **`accumulate_rig` does not exist.** Three rigs call it on `ContractJointConnector`
   (`microtubule_rig.py:467`, `intermediate_filament_rig.py:497`, `stress_fiber.py:353`); the class
   exposes `accumulate` / `accumulate_contact` / `accumulate_actor` / `accumulate_motor`. Guaranteed
   `AttributeError` the day `mt_cortex_capture` binds.

**And a trap for whoever wires the membrane:** `build_native_surface_owner` hardcodes
`_NoOpSurfaceMechanics` (`erm_cortex_slice.py:151`). Binding it as-is yields a bound component carrying
zero force — the job-94 failure exactly. Pass a real delegate.
