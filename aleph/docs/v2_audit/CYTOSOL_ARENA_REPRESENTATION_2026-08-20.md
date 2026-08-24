# What the cytosol IS in the arena — the question `arena.py` left open, answered with evidence

**Status:** DECISION REQUEST to Lead. Session `world/cytosol` (ComponentRole `FLUID_VOLUME`),
2026-08-20. **Selects no constant and sets no threshold.** It answers one representation question,
reports what the answer costs, and records one KB row that does not hold for this cell.

**Nothing was built.** `aleph/world/build/cytosol.py` does not exist yet, deliberately: the answer
below needs one line of `arena.py`, which this session does not own, and guessing it is the exact
failure `arena.py` warns about — *"guessing would bake an answer into the layout."*

---

## 1. The answer, in one line

**The cytosol is a GRID, not a node population, and the arena needs an eighth `Kind` —
`GRID_CELL` — which `arena.py` already presupposes it has.**

## 2. Why it is not a node population — three independent pieces, none of them a preference

1. **The port source is Eulerian and says so.** `components/fluid/field_grid.py:9-11` —
   *"the Eulerian Biot/RAD PDE fields live on this structured, conservative finite-volume grid"* —
   and the same paragraph separates it from the Warp `HashGrid` used for filament neighbour search:
   *"a neighbour structure, NOT the PDE field grid. They may share origin/cell size but are distinct
   allocations."* The distinction the arena would collapse is one the port source already made
   explicitly.
2. **The component contract requires rank 3.** `engine/fluid_core.py:215-243`
   (`FluidVolumeStateOwner`) validates `pressure_d` and `mask_d` with `ndim=3`, float64 and int32,
   and rejects anything else. A node population is rank 1. This is not a formatting detail — the
   Biot update kernel indexes `p[i, j, k]` with a 6-point masked Laplacian
   (`biot_substrate.py:80`), so the stencil IS the data structure.
3. **The six `*_cytosol_transfer` connectors do not join node to node.** `ImmersedPorousTransfer`
   (`engine/immersed_transfer.py:141-152`) scatters `-alpha * V_node * grad(p)` onto the SOLID's own
   force array through a Peskin-4 interpolation. The cytosol end of that edge is **a grid stencil
   around a position**, not an ID. There is no cytosol node to bond to, and inventing one would be a
   lumped stand-in for an interpolation — the thing the charter's architectural principle forbids.

So the six `NOT_BOUND` connectors are not waiting for a node population. They are waiting for a field.

## 3. Why `GRID_CELL` and not something else — the arena already wrote this down

`arena.py:71-83` says the seven kinds are *"closed by argument rather than by convenience"* and then
justifies rejecting a species pool like this (`arena.py:77-78`):

> a species pool is a named channel on a ``GRID_CELL`` field plus one conserved scalar per consuming
> range.

**`GRID_CELL` is not in the enum.** The set that is closed by argument has an argument naming an
eighth member. This session is not proposing a new primitive; it is reporting that one was assumed
and never added.

## 4. What adding it costs — almost nothing, because the precedent is already there

The arena allocates device arrays **only for `NODE`** (`arena.py:170-181`). `SEGMENT`, `ANGLE3` and
`FACE` are pure bookkeeping — capacity, claim, partition invariant — and the builder owns the
kind-specific arrays: `build/cortex.py:429-440` owns `seg_node (S,2)`, `seg_rest_um (S,)`,
`angle_idx (A,3)`; `build/membrane.py:577-589` owns `face_idx`, `hinge_idx`, `area0`.

`GRID_CELL` fits that shape exactly:

| | arena | `build/cytosol.py` |
|---|---|---|
| capacity `nx*ny*nz`, one claim, `assert_partitioned` | ✅ | |
| `p`, `p_bar`, `s_water`, `s_membrane`, `s_total`, `div_vs`, `mask` (`wp.array3d`) | | ✅ |
| `dx`, `origin`, shape | | ✅ |

The 1-D/3-D mismatch is not a problem and is not new: a `Claim` counts `S` segments while the array
is `(S, 2)`. Here the claim counts `nx*ny*nz` cells while the arrays are `(nx, ny, nz)`. The
partition and byte-span invariants hold on the flat count either way.

**One arena change is needed and it is one enum member.** No new allocation branch, no new array set.

## 5. The option that was considered and is rejected on a measurement

*Keep the grid outside the arena entirely, as the incumbent does.* Rejected because **PHASE 1 exists
to measure `exact_peak_gpu_bytes`** — the plan calls it *"a number no session has ever taken"* — and
§7 below shows the fluid grid is plausibly the **largest single allocation in the cell**. A PHASE 1
measurement taken with the biggest allocation outside the arena is a measurement of the smaller half
of the cell, reported under a whole-cell name. That is the failure mode the plan was rewritten to
avoid.

## 6. What the "density" is for a FLUID_VOLUME — and the honest answer is that it is not one

There is no per-cell count to source. A grid has a **resolution**, and `dx` is a discretisation
choice, not a physiological density. But it is bounded on both sides by sourced physics:

| bound | value | KB row | status |
|---|---|---|---|
| **floor** — below the pore/mesh size a continuum pore pressure is not defined | ξ ≈ 14–40 nm | `KB-DRAFT-3.B-27`, `KB-DRAFT-3.B-18` | **draft** |
| **ceiling** — the poroelastic diffusion length over one outer step, √(D_p·dt_phys) = **1.58 µm** at D_p = 50 µm²/s, dt_phys = 0.05 s | D_p = 40–60 µm²/s | `KB-3.B3.2` | **verified** |

The incumbent's value sits inside that band and was not chosen from it:
`components/incumbent/assemble.py:333` — `dx_um: float = 0.5  # conservative field-grid resolution
(numerical sizing; not literature)`. That is the same class of value as the cortex's `seg_um = 0.5`,
which this phase already replaced with a sourced band.

**Request:** `dx_um` is a declared axis with no default, and `build/cytosol.py` refuses without it —
the same refusal every other builder makes. **Which value inside the band is a PI decision, not a
derivation**, and §7 is why it is not a small one.

## 7. ⚠ The finding Lead needs before deciding: `dx` enters cost as **dx⁻⁵**

Cells scale as dx⁻³. The explicit consolidation CFL is `safety·S·dx²/(2·d·mobility)`
(`field_grid.py:125-127`), i.e. `0.9·dx²/(6·D_p)` — so **subcycles per outer step scale as dx⁻²**.

**Reproducible from code**, not from a scratch script: `build/cytosol.cytosol_counts()` computes
every column with no card and no `Kind`, and `python -m aleph.world.build.cytosol` prints it. Box
derived as R_cell + one OUTSIDE cell layer, D_p = 50 µm²/s, dt_phys = 0.05 s, 60 B/cell (7×float64 +
int32):

| dx [µm] | shape | cells | GB | CFL dt [s] | subcycles/step | relative fluid work |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 (incumbent, *"not literature"*) | 33³ | 35,937 | 0.002 | 7.5e-04 | 67 | 1× |
| 0.25 | 63³ | 250,047 | 0.015 | 1.9e-04 | 267 | 28× |
| 0.2 | 79³ | 493,039 | 0.030 | 1.2e-04 | 417 | 85× |
| 0.1 | 153³ | 3,581,577 | 0.215 | 3.0e-05 | 1,667 | **2,480×** |
| 0.05 (the cortex's sourced floor) | 303³ | 27,818,127 | 1.669 | 7.5e-06 | 6,667 | **77,027×** |

⚠ An earlier draft of this table used a fixed 16 µm cube and read 33.1 M / 1.985 GB / 91,586× at
dx = 0.05. The box is now derived from `R_cell` rather than typed, which is smaller and is the
number the build would actually claim. **The conclusions are unchanged and the table above is the
one to quote.**

Three consequences, none of them this session's to resolve:

* **`WORLD_PORT_PLAN` §1's node budget has no cytosol row at all.** At dx = 0.05 µm the grid is
  **27.8 M cells against 4.86 M solid nodes** — 5.7× the entire rest of the cell in element count and
  3.9× in bytes (solid ≈ 0.43 GB at 88 B/node). The 10 M capacity figure was never asked to hold a
  field.
* **It runs against §3's only lever.** `k_xl` is worth ~6,500×; moving the grid from 0.5 to 0.05 µm
  costs ~77,000×. Whatever the `k_xl` decision buys, the `dx` decision can spend several times over,
  and the two have never been priced on the same page.
* **CLAUDE.md already said this and it is now quantified.** *"the cost is the fluid grid and the
  moving membrane boundary."*

## 8. KB-3.21 checked as instructed — **verified, and it does not hold for this cell**

`KB-3.21` is `status: verified`, `confidence: High`, and does carry `Phase-1 c_G~25 uM V_cyto~4200
um^3`. Citation audit: `PollardBorisy2003_Cell` **OK**, `MogilnerOster2003_BJ` **OK**,
`Carlsson2018_JPhysCondMat` **`NO_DOI_FOUND`** (the audit's candidate match is an unrelated T4
lysozyme paper).

But the arena's own sourced geometry gives a different number:

| | value |
|---|---:|
| `laws/cell_geometry` MCF7 `R_cell` | 7.50 µm |
| `V_cell` | 1,767 µm³ |
| nucleus at sourced `nc_ratio` 0.68 → `R_nuc` 5.10 µm, `V_nuc` | 556 µm³ |
| **V_cell − V_nuc (upper bound on V_cyto)** | **1,211 µm³** |
| KB-3.21 `V_cyto` | 4,200 µm³ |
| **ratio** | **3.47×** |

4,200 µm³ is a sphere of radius **10.0 µm**. The arena builds a 7.5 µm cell. The row is not wrong —
it is **scoped to a different cell**, which is precisely `bond.SourceClass.UNRATIFIED_PROXY`:
*"a value can be perfectly SOURCED and still be wrong for the run … the SCOPE is what was never
recorded."*

**Consequence:** V_cyto cannot be used as this build's acceptance check. The check available is the
one the build produces itself — `Σ_FLUID dV` from the mask, against `V_cell − V_nuc` from the same
geometry the membrane and envelope builders already use. That is self-consistent and needs no
external row. **KB-3.21's V_cyto is not quoted by this population.**

## 9. Out of scope, flagged not solved

* **`BondFamily` cannot express an immersed transfer.** `bond.py:189` stores `node_i`/`node_j` as
  global NODE indices and `component_pairs` resolves them with `arena.population_of(Kind.NODE, i)`,
  hardcoded. A cytosol endpoint is a GRID_CELL. With the kind added, the fix is a kind argument —
  but it is PHASE 3 and it belongs to whoever owns those six family modules.
* **Nucleoplasm — CLOSED 2026-08-20 by session F, against this memo's leaning. One grid,
  mask-indexed properties, one `dx_um`.** This memo had left it open and pointed at KB-3.31's
  nucleoplasm η ≈ 52 Pa·s as a reason the resolutions *might* differ. That inference was loose and F
  refuted it; the refutation is verified here against source, not relayed:
  - **The CFL does not see bulk η.** Darcy mobility is `k/µ` with µ the PORE fluid viscosity
    (`biot_substrate.py:5`), and KB-3.31's own text instructs *"Wire pore viscosity into the fluid
    channel … poroelastic, not one lumped viscosity"* — pore viscosity is 2–8× water in both
    compartments. Quoting the bulk rheology at a resolution decision is precisely the lumping that
    row warns against.
  - **The quantity that would split `dx` is absent.** `D_p ~ E ξ²/µ`, and `Moeendarbary2013_NatMater`
    p9 says verbatim: *"measurements were acquired in several locations in the cytoplasm **avoiding
    the nucleus**."* So `KB-3.B3.2`'s 41 ± 11 µm²/s is cytoplasm-only and `ξ_nuc` is not in the
    corpus. A second `dx` would rest on a number that does not exist.
  - **Different properties ≠ different grid.** `biot_pmass_update_kernel` already reads the mask per
    cell; making `mobility`/`storage_S` mask-indexed gives the nucleoplasm its own properties on one
    grid. ⚠ Not free, and F corrected their own claim here: `biot_substrate.py:62-64` holds every
    non-FLUID cell (`mask[i,j,k] != _FLUID → p_new = p; return`), so **the nucleus today is a
    no-flux hole with no pressure DOF**, not a fluid with other properties. Promoting it is a
    prerequisite — one central branch and a 2-element coefficient lookup, orders of magnitude below a
    second grid.
  - **Reversal is asymmetric.** Adding a second `GRID_CELL` claim later is what the arena's range
    model made cheap; breaking "a flat range indexes one grid" now, on a guess, is the expensive
    direction. The landed `CytosolField` carries `shape`/`dx_um`/`origin_um` per instance, so a
    second grid stays representable and nothing yet assumes `n_live(GRID_CELL)` is one grid.
  - ⚠ **Reopen trigger, with the sign reversed from this memo's guess:** a sourced `ξ_nuc` or
    `D_p_nuc`. If it arrives, the candidate is a **coarser** nucleus, not a finer one — lamin
    meshwork 0.4 µm and chromatin domains 0.6 µm (F, `Stephens2017_MBoC` Table 1, `13dfcff7`) against
    the cortex's 0.05 µm, 8–12× the other way. That proposition — *"`dx` is set by the finest
    structure exchanging water"* — is itself unratified.
* **`fluid_core.py` in `aleph/world/`** — named in this session's scope as a possible second file.
  Not needed under this design: the grid arrays live with the builder, exactly as `seg_node` lives
  with the cortex builder. One file.

## 10. LANDED — the refusing builder, at Lead's direction

`aleph/world/build/cytosol.py` is in the tree and **refuses**, which is its product for now. It does
not touch `arena.py`: the enum line is Lead's to add after approval, and a line added before it would
read as the approval.

* `grid_cell_kind()` raises **`ArenaKindGapError`** — a named type, the sibling of
  `families.ConnectorGapError`, carrying the whole argument and pointing at this document. Never an
  `AttributeError`: a reader opening the file must see *waiting*, not *broken*.
* All three axes in `CYTOSOL_AXES` carry `value: None`. `dx_um` has no default and the build refuses
  without one; a `dx` below the sourced continuum floor is refused with both ends of the band named,
  because a finer grid there is not a finer discretisation of the same field — it is a different one.
* `cytosol_counts()` is the §7 table, **as code**, with no card and no `Kind`. That is what makes the
  `dx` decision arithmetic instead of taste, and it is available to the PI now.
* The build path is written and its four Warp kernels **compile** (syntax only — nothing was
  launched and no number was produced; there is no CUDA on the dev machine and a CPU number would not
  be a result). It has never executed.

**The classifier was ported, not invented.** Classification is against the surfaces that were BUILT,
not against ideal spheres — `build/cortex.py`'s own rule. The algorithm is
`wp.mesh_query_point_sign_winding_number` over a `wp.Mesh` of the live arena vertices, with an
explicit float32 surface-tie epsilon, taken from `components/incumbent/live_mesh_domain.py:110-125`.
⚠ Worth knowing: the fluid track's own `components/fluid/domain.py` ships only
`StaticSphereMembraneProvider` / `StaticSphereNucleusMaskProvider`, which its docstring labels a
**placeholder**. The real mesh classifier lives in `incumbent/`, so this is a port past that
placeholder rather than a port of it.

Proof when the `Kind` lands: the population stands, `assert_partitioned` passes, and `Σ_FLUID dV` is
reported against `volume0_um3(membrane) − volume0_um3(envelope)` — both measured from what those
builders built. **No physics is claimed**: fields are allocated and left at zero, no Biot step, no
pressure, no flux.

`world/fluid_core.py` is not needed under this design — the grid arrays live with the builder, the
way `seg_node` lives with the cortex builder. One file.

**Still blocked on two things, and only two:** (1) `Kind.GRID_CELL`, (2) the PI's `dx_um`.
