# Worked example — the V1 vertical as an L1 plan

| Field | Value |
|---|---|
| Status | `SPEC` — the acceptance example `aleph/represent/plan.py` must reproduce. That module now exists, and §4's tie-break disagreement is a live test in `tests/represent/test_plan.py`. |
| Written | 2026-07-30 12:47 KST |
| Scope | `ALEPH-DQ-101` V1: membrane–cortex–pressure, the smallest slice with a closed-form oracle |
| Purpose | An abstract language cannot be wrong. A worked plan can. |

## 0. Why this document exists

`ALEPH-TN-LANGUAGE-SPEC.md` states rules. Rules are unfalsifiable until something is derived from
them, so this document runs the language once, by hand, over the smallest real slice of the census —
and every step below is a claim that the implementation must reproduce or contradict. **If the
implementation disagrees with this document, one of the two is wrong and the disagreement is the
finding.** Reconciling by quietly editing this file would waste the only test the design has.

## 1. The slice

Owners (from the landed census), with scope tags:

| Owner | Scope | Why it is in V1 |
|---|---|---|
| `membrane` | `E` | the surface whose fluctuation spectrum is the `ALEPH-DQ-102` observation |
| `cortex` | `E` | the load-bearing shell under it, and the osmotic envelope |
| `cytosol` | `E` | carries the pressure that loads the two above |
| `extracellular_medium` | `E` | the only owner that can load the six rigid-body modes |

Connectors between them:

| Connector | Family | Endpoints | Note |
|---|---|---|---|
| `membrane_erm_cortex` | KINETIC | membrane ↔ cortex | tensile only; owns bound/unbound state |
| `membrane_cortex_contact` | CONTINUOUS | membrane ↔ cortex | the compressive path; **not** the same load path as the tether |
| `membrane_cytosol_boundary` | CONTINUOUS | membrane ↔ cytosol | pressure load onto the surface |
| `membrane_medium_traction` | CONTINUOUS | membrane ↔ extracellular_medium | exterior Stokes resistance |

The tensile/compressive pair is the trap the earlier audit found and it is the reason this slice is
the acceptance example: a plan that treats the two as one connector will look tidier and be wrong.

## 2. Blocks (spec §3)

| Owner | Blocks | Derived from |
|---|---|---|
| `membrane` | `MESH_SURFACE`, `EVENT_STREAM` | vertex/element state; topology events |
| `cortex` | `OWNER_PARTICLES`, `EVENT_STREAM` | disjoint filament population; binding events |
| `cytosol` | `CONTINUUM_FIELD` | pressure and flow coefficients with units and a domain |
| `extracellular_medium` | `CONTINUUM_FIELD` | committed surface reference, traction, banked dissipated work |
| all four connectors | `CONNECTOR` | endpoint refs, force/work/energy/dissipation, kinetic state |
| the fluctuation protocol | `OBSERVATION` | raw output plus calibration and analysis lineage |

Non-state entries that still register: temperature and osmotic condition are `B` →
`ContextConditioningAxis` in the context manifest, never state axes. Unresolved solvent drag is `H` →
a `HomogenizedLawSlot` with its declared inadequacy. Organelle mechanics is `X` → an `ExcludedRegion`
carrying its re-entry condition.

**The check with teeth:** `kB*T` may not enter as a state index. It is context, its evidence class is
`ASSUMED`, and no entry point may supply it silently. A plan that indexes temperature has violated
invariant I3 and must be refused.

## 3. Candidate cuts (spec §4)

| Cut | Crossing connectors | Size | Legal? |
|---|---|---|---|
| `{membrane} \| {cortex, cytosol, medium}` | erm_cortex, cortex_contact, cytosol_boundary, medium_traction | 4 | legal, but the widest — the membrane is the most connected owner |
| `{membrane, cortex} \| {cytosol, medium}` | cytosol_boundary, medium_traction | 2 | **legal and narrowest** |
| `{membrane, cytosol} \| {cortex, medium}` | erm_cortex, cortex_contact, medium_traction | 3 | legal, and worse — it separates the two owners that share two load paths |
| `{cortex} \| rest`, splitting a composite group | — | — | **refused** `CUT_SPLITS_COMPOSITE` if the tensile/compressive pair is registered as a composite group |

The narrowest cut puts the two mechanically-paired surfaces on one side and the two fluids on the
other. That is not a spatial statement: the medium and the cytosol are on opposite sides of the
membrane in space, and they sit together here because neither shares a connector with the other.
**This is the concrete demonstration that the bond structure follows the connector graph.** A plan
built on spatial adjacency would have grouped cytosol with cortex and paid for it in bond width.

Every crossing connector must be `implemented` (rule C4). At the time of writing, `implemented=False`
is the default in the connector registry, so a plan over this cut refuses today with
`CUT_CROSSES_UNIMPLEMENTED` — which is correct behaviour and worth keeping visible rather than
special-casing.

## 4. Axis order (spec §5)

Query family: the membrane fluctuation spectrum, so the root set is `{membrane}`.

1. BFS from `membrane` over `E`-scope owners: depth 0 `membrane`; depth 1 `cortex`, `cytosol`,
   `extracellular_medium`.
2. Depth-1 tie broken by **coupling descending** (`TieBreak.STRONG_COUPLING_ADJACENT`, PI-decided
   2026-07-30 ~14:56): `cortex` is reached by two connectors, `cytosol` and `extracellular_medium` by
   one each → `cortex`, then `cytosol`, then `extracellular_medium`.

   *This is where working the example by hand paid for itself.* The draft rule sorted cut size
   ascending, which exiled the doubly-connected `cortex` to the far end of the chain while its own
   stated rationale — keep weakly coupled owners at the ends — described the opposite. Since the bond
   dimension across a cut is set by the correlation crossing it, separating membrane from cortex would
   have made every intermediate cut carry the two load paths between them. The rule now matches the
   reason, and `WEAK_COUPLING_ADJACENT` remains available for the measurement that could overturn it.
3. Within `membrane`: `MESH_SURFACE` before `EVENT_STREAM` (fixed block order).
4. Within `MESH_SURFACE`: `SECTOR`, then `ENTITY` (vertices), then `SPATIAL` (extent 3), then `FIELD`,
   then `MODE`.
5. The membrane's sector-dependent axis leads, so a topology jump re-plans a prefix rather than the
   whole chain.

## 5. Families (spec §6)

| Block | Planned family | Why, and what it must beat |
|---|---|---|
| membrane `MESH_SURFACE` | `SPECTRAL_BASIS` | the Helfrich mode spectrum is the observable; the basis is the physics, not a compression |
| cortex `OWNER_PARTICLES` | `DENSE` | local TT on particles is admissible only after measured qualification, and there is nothing to measure against |
| cytosol `CONTINUUM_FIELD` | `TT_MPS` planned, `SVD_POD` as the baseline it must beat | a compression that cannot beat a linear subspace has bought complexity with nothing |
| connectors | `TYPED_EDGE` with the exact `+f/−f` layer | action–reaction is architectural, never a penalty |
| events | `DENSE` | the stream is short in V1; a marked-process encoder is L2 |

Every entry above says *planned*. None says *works*. Qualification is L2 and sealed.

## 6. What this example refuses, today, on purpose

| Question | Answer now |
|---|---|
| What bond dimension does the narrow cut need? | `BondDimensionUnmeasured` — it is the interaction rank, and the rank is measured |
| What is the truncation error of the cytosol TT? | `RepresentationSealed` — no accepted native trace exists |
| Is `SPECTRAL_BASIS` good enough for the membrane? | unanswerable; the observation-channel tolerance slot is empty |
| Can we train the SIM reader on this slice? | no — the lane gate needs an accepted trace digest and none exists |

Four refusals out of four questions is the expected result at this stage, and a version of this
document with fewer refusals would mean someone had started guessing.

## 7. Two findings this exercise produced, before any code ran

1. **Cut size versus query relevance at equal BFS depth is undecided** (§4 step 2). The spec asserted
   cut size without noticing that the two criteria disagree on the very first realistic example.
2. **`implemented=False` makes every V1 cut refuse today.** That is correct, and it means the
   representation lane's first green integration test will arrive with the runtime rather than with
   the language — worth knowing now rather than discovering it as a mysterious failure later.
</content>
