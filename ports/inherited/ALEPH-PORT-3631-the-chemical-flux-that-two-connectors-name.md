# ALEPH-PORT-3631 — the chemical flux two connectors name and nothing carries

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3631` |
| Lane | `46143f30` Lane W2 (lead) |
| Status | `PROPOSED` |
| Written | `2026-08-05` — **before the code**, per `CLAUDE.md` §3 |
| Port class | `NOT A PORT` — nothing is taken from `/Users/sw1/ffn_cellsim`; the transport law is written from the control-volume mesh this owner already has |
| Aleph target | `aleph/vertical/cytosol.py` (append-only: one phase's worth of law and its state evolution), `tests/vertical/test_cytosol_species_transport.py` (new) |
| Exists because | `CytosolField.species_content` is allocated, committed and rolled back, and **nothing advects or diffuses it**. Two registry rows name that channel and cannot be wired without it. |
| PI instruction | Implement chemical flux **without inventing parameters** — "chemical_flux 그냥 파라미터 없이 구현, 이것도 나중에 스윕에 추가될 것임". The parameters join the sweep later; this entry adds the *law*, not numbers. **Coordinates, per `CLAUDE.md` §2 rule 2** (added 2026-08-05 after an audit found this quote uncited): transcript `~/.claude/projects/-Users-sw1-Project-Aleph/46143f30-d1ab-4f8a-9c66-b63204a3b7a4.jsonl`, record uuid `fbfdbd98-b70a-4f5f-9b70-1cc026aae84d`, timestamp `2026-08-04T16:02:48.971Z`, delivered as a `/goal` command. A session that did not receive it can open that record and check. |

---

## 1. What is missing, in the owner's own words

`CytosolField.species_content`:

> State key `transported_species_content`. Allocated, committed, rolled back — and **nothing advects
> or diffuses it**. The G-actin channel that two of the nine connectors name does not exist.
> Registering the slot and leaving it inert is the honest encoding.

That is exactly right and it is why two rows are withdrawn. `connectors_transfer`'s own guard
refuses to wire half of them:

> Both protrusion contracts declare two couplings on one edge, and the monomer half cannot be
> written today […] **Wiring the drag half alone supplies exactly the infinite reservoir both
> contracts name as the error to avoid.**

An infinite reservoir is a cell whose actin never runs out: a lamellipodium would extend forever at
its free polymerisation rate, and every protrusion number the engine produced would be an artefact
of monomer that was never accounted.

## 2. No new parameters, and that is a constraint on the law rather than a wish

The PI's instruction is that this entry adds no numbers. It does not need to, and that is worth
stating as a fact about what was already built rather than as a convenience:

| what the law needs | where it already is |
|---|---|
| diffusivity `D` | `CytosolField.poroelastic_diffusivity_um2_per_s` — derived from the card's mobility, storage and Biot coefficient; **3.983 µm²/s** on this fixture, not a constant anyone typed |
| pore fluid velocity | `internal_face_flux_um_per_s()`, the Darcy flux this owner already solves |
| face topology and areas | `faces.lo_dof`, `faces.hi_dof`, `faces.axis`, `face_area_um2` |
| control volume | `node_volume`/`dx³`, already used by the Biot step |

So the entry supplies **structure**: which terms exist and how they are discretised. Every magnitude
comes from state the owner already computes. `ALEPH-PORT-1803` forbids a viscosity magnitude in
code; this entry does not add one and does not need one.

## 3. The law

Content `c_i` per control volume [µm³ of monomer, or whatever unit the sweep later fixes]:

    dc_i/dt = - Σ_faces( J_adv + J_dif )·A_f  +  S_i

with, per internal face `f` between `lo` and `hi`:

    J_dif = -D (c_hi/V_hi - c_lo/V_lo) / dx        Fick, on CONCENTRATION not content
    J_adv =  q_f · (c_lo/V_lo if q_f > 0 else c_hi/V_hi)      first-order upwind

**Fick's law is on concentration, not content**, and getting that wrong is the first thing a reader
should check. Content is extensive; two cells holding the same content at different volumes are at
different concentrations and diffusion is driven by the difference of the intensive quantity. On a
uniform grid every `V` is `dx³` and the distinction vanishes numerically — which is exactly why it
must be written correctly now, before a non-uniform mesh makes it wrong silently.

**Upwind, not central**, for advection. A central difference on an advection-dominated cell produces
oscillations that go negative, and a negative monomer content is not a small error: it is a cell that
has anti-actin in it, and every downstream rate reads it as real. Upwind is first-order and diffusive
and that cost is stated rather than hidden.

## 4. Conservation, which is the whole check

The scheme is a **flux form**: every internal face contributes `+F` to one cell and `−F` to the
other, by construction. Therefore

    Σ_i c_i  is conserved exactly, to round-off, when no boundary flux and no source

and that is asserted as an equality rather than a tolerance, because the two contributions are the
same number with opposite signs. It is the same argument that makes a connector's Newton pair
meaningful, applied to a scalar.

**Boundary faces are a source or a sink, never silent.** A face on the membrane or the nuclear
envelope either carries a declared flux or carries none; a scheme that quietly let content leave
through an unmodelled boundary would look conservative on the interior and lose mass at the edge.

## 5. Stability, and a refusal rather than a fudge

Explicit in time, so the step is bounded:

    dt ≤ dx² / (6D)          diffusion
    dt ≤ dx / max|u|         advection (CFL)

A step exceeding either is **refused by name**, not clamped and not silently sub-stepped. A clamped
step is a different simulation from the one the caller asked for, and returning it as though it were
the requested one is how a trajectory stops meaning what its `dt` says. The refusal names which
bound was crossed and by how much, so the caller can pick a `dt` rather than guess.

On this fixture `dx = 0.5 µm` and `D = 3.983 µm²/s`, so the diffusive bound is `dt ≤ 1.05e-2 s`,
which is far above the relaxation's `1e-4`. **The bound is therefore not expected to bite here**, and
saying so now means a future failure is a change in the physics rather than a surprise.

## 6. What this does NOT do

- **It does not wire the two protrusion rows.** It supplies the channel they need; the connectors
  themselves need an elongation rate that `filopodium` currently refuses as unsourced, and a
  conversion from µm of F-actin to a monomer count that `lamellipodium` cannot do. Those are the
  next entry's, and the guard that withholds the builders stays exactly as it is until they exist.
- **It does not model actin chemistry.** There is one transported scalar with no reactions, no
  binding, no nucleotide state. Naming it "G-actin" would be a claim; it is a transported species
  and the sweep will decide what it is.
- **It does not source or sink anything yet.** `S_i` is present in the discretisation and zero
  everywhere, because a source with no polymerisation attached to it is a number nobody chose.

## 7. Units, domains, singular cases, invariants

Units: length µm, time s, diffusivity µm²/s, face flux µm/s, area µm², volume µm³. Content is in
whatever unit the sweep later fixes; the law is linear in it, so the choice does not enter here.

**Singular cases:**

| case | treatment |
|---|---|
| a face between a fluid cell and a non-fluid one | not an internal face; carries no diffusive term |
| zero diffusivity (`decouple_solid=True` gives `D = 0` hard) | the diffusive term vanishes exactly; advection survives, which is correct |
| zero face flux | the upwind branch is taken on the `q > 0` side by convention and both give the same answer, since the coefficient is zero |
| negative content arriving from a source | refused: content is a count and a negative one is not a small error |

**Invariants:**

| invariant | how it holds |
|---|---|
| total content conserved with no boundary flux and no source | flux form: `+F` and `−F` are the same float |
| content stays non-negative under the CFL bound | upwind + the step bound; asserted, and the bound is refused rather than clamped |
| `D = 0` reduces to pure advection | the term is multiplied by `D`, not branched on it |
| a uniform field does not move | zero concentration gradient and divergence-free flux both give exactly zero |

## 8. Positive control, negative control

**Positive** — `test_a_uniform_field_does_not_move`: a constant concentration with any flux field
must stay constant to round-off. It is the discretisation's own consistency: any sign or index error
in the face loop breaks it immediately, and it needs no reference solution.

**Positive** — `test_total_content_is_conserved_to_round_off`: sum before and after, `==` rather than
`approx`, over many steps.

**Negative controls**, shipped as deliberate breaks on the field so they can be driven:

| flag | what it breaks | what must catch it |
|---|---|---|
| `advect_with_central_difference` | second-order but oscillatory | content goes negative on a sharp front |
| `diffuse_on_content_not_concentration` | Fick on the extensive quantity | wrong on a non-uniform mesh; **identical on this one**, and that is recorded as a declared blind spot rather than left to look like a passing test |
| `leak_one_side_of_each_face` | the flux form's `+F/−F` symmetry | total content stops being conserved |

The second is the interesting one and it is named as blind on the uniform fixture **before** anyone
runs it, so the count of caught mutants is not quietly inflated by a mutant that cannot differ.

## 9. Numerical envelope

Host `float64`. The conservation claim is exact (`==`), because the two face contributions are the
same value negated. Nothing else here is a precision question; there is no reference implementation
to grade ULP against, and that is stated rather than implied — this is a **new law**, not a port, so
it has no parity oracle and its evidence is the invariants in §7 plus the controls in §8.

## 10. Residency

Host NumPy, evolving with the rest of the cytosol. The face loop is a gather/scatter over
`faces.lo_dof` and `faces.hi_dof`, which is exactly the shape `backend.scatter_add` takes, so it is
written to move to the device with the rest of this owner and not as a special case.

## 11. Source identity, why not clean-room, discarded prose

No source repository, path, symbol or commit: nothing was ported. Upwind finite volume on a
control-volume mesh is textbook and was written here against this owner's own face table. No foreign
prose exists to discard.

## 12. Acceptance, reviewer, rollback

Status `PROPOSED`. Reviewer: the PI. Rollback: the phase does nothing if the law is not scheduled,
and `species_content` returns to being inert state — which is where it was, and was honestly labelled
as such.

## 13. Open

**What is the species?** This entry transports one scalar and does not name it. The two protrusion
contracts want G-actin monomer; a real cell also transports ATP, calcium and every diffusing
regulator. Whether this channel is *the* monomer channel or the first of several is a modelling scope
question and it is not settled here.
