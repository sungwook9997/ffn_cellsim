# Closed-Loop ECM Gate — Phase A Evidence Ledger

**Status**: **NON-LOCKING evidence ledger only.** Records existing
tests/docstrings/commits that already satisfy the items the
phased-plan brief
(`docs/v2_closed_loop_ecm_gate_dependency_brief.md`) calls "Phase A
items provable now under existing ECM-OL preflight". The ledger
makes no design decision, commits to no constitutive law, and
holds no opinion on the closed-loop gate items beyond reading off
the current state.

**Source ledger**: implementation-work, drafted while
design-discussion `id=1223` is in flight on the broader gate
phased plan. Per Codex impl `id=1225`, this Phase A note is
explicitly readback/evidence-only and is **superseded** as soon as
the design-discussion lock artifact (target:
`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`) lands with
a different phase split or different evidence requirements.

**What this note is not**:
- Not a Sanity Gate document. Phase A items are pre-existing and
  already passed their own gates (ECM-OL preflight commits were
  reviewed individually).
- Not a closed-loop activation. No FA→ECM, no ECM bias, no
  constitutive law decision is recorded here.
- Not a commitment that Phase A is sufficient. The full gate has
  six items; this note covers only the two that are already
  trivially or near-trivially provable under the current open-loop
  preflight state.

---

## Item 3 — No response when traction is zero

### Closed-loop gate text (verbatim, roadmap §"Closed-Loop ECM Gate" item 3)
> No response when traction is zero.

### Open-loop side: trivially satisfied

Test:
- `tests/test_v2_ecm_ol_harness.py::test_zero_traction_scenario_no_accumulation_change`
  - Uses `_zero_traction_factory((4, 4))` over 5 steps.
  - Asserts `run.final_ecm.accumulated_traction_nNs_per_um2.max() == 0.0`.
  - This is the open-loop preflight's direct realisation of "zero
    traction in ⇒ zero accumulator change out".

Code that ensures this:
- `acs/v2/dynamics/ecm_open_loop.py::accumulate_prescribed_traction`
  computes `accumulated_new = accumulated + traction · dt`. With
  `traction == 0`, `accumulated_new == accumulated` exactly to
  float64 round-off (in fact bit-exact since the addition is
  `x + 0`).

### Closed-loop side: pending

Once a closed-loop response field is added (e.g., a
`stiffness_kpa` reduction under cumulative traction), the same
zero-traction test must be re-run on the response field. That test
will be added at the time the response law lands; the design
decision on the response law is a hard blocker enumerated in the
phased-plan brief.

### Status

Open-loop side ✓. Closed-loop side blocked on Hard Blocker 1
(constitutive-law direction).

---

## Item 6 — Rule 10 unit-chain proof for traction-density storage and comparisons

### Closed-loop gate text (verbatim, roadmap §"Closed-Loop ECM Gate" item 6)
> Rule 10 unit-chain proof for traction-density storage and
> comparisons: `nN`, `um2`, `nN*s/um2`, and kPa / force-per-area
> comparisons before any ECM force or remodeling update.

### Open-loop side: storage proven, comparisons N/A

The current ECM-OL preflight performs **storage** updates only;
it does no force or remodeling comparison across mismatched units.
Each preflight function performs exactly one Rule 10 unit-chain
reduction:

| Function | Unit chain | Storage field |
|---|---|---|
| `accumulate_prescribed_traction` | `[nN/μm²] · [s] = [nN·s/μm²]` | `accumulated_traction_nNs_per_um2` |
| `apply_prescribed_stiffness_rate` | `[kPa/s] · [s] = [kPa]` | `stiffness_kpa` |
| `apply_prescribed_density_rate` | `[1/s] · [s] = [dimensionless]` | `ligand_density`, `fiber_density` |
| `apply_prescribed_orientation_rate` | `[1/s] · [s] = [dimensionless]` | `orientation_tensor` |

Source of truth: `acs/v2/dynamics/ecm_open_loop.py` module
docstring (post commit `22b9456`) §1, lifted verbatim from each
function's individual unit reduction in 6.4-open-A/B/C commits
(`7dc1767`, `e4a71c9`, `8b2c1ab`).

### kPa ↔ nN/μm² comparison: not yet performed (closed-loop concern)

The Hard Rule 10 reminder in the gate explicitly cites "kPa /
force-per-area comparisons before any ECM force or remodeling
update". The open-loop preflight performs no such comparison. The
two units are dimensionally identical:
`kPa = kN/m² = 10³ N/m²`, and `nN/μm² = 10⁻⁹ N / (10⁻⁶ m)²
= 10³ N/m²`, so **1 kPa = 1 nN/μm² exactly**. (The first commit
of this note carried a `10⁻³` factor that was wrong by three
orders of magnitude — flagged by Codex review id=1234 and
re-derived above.) Even with the factor of 1, the project rule
still mandates that the comparison be written with both quantities
reduced to the same unit at the call site, not implicitly at the
storage layer.

This comparison only arises when a closed-loop FA→ECM coupling
function feeds per-FA `nN` traction into the ECM as a per-cell
`nN/μm²` (force-per-area) deposition, or when an ECM stiffness
(`kPa`) field scales an FA `traction_scale_nN`. Neither exists
yet; both are gated behind the phased-plan brief's Phase D
scaffolding.

### Status

Open-loop storage side ✓. Closed-loop force-per-area comparison
side pending the closed-loop coupling functions, each of which
must record the unit-chain proof in its own Sanity Gate.

---

## What is deliberately out of this evidence note

- Item 1 (response monotonicity): only the trivial accumulator
  monotonicity is provable now; Item 1 in the closed-loop sense
  awaits a response-field definition. Not recorded here as
  "satisfied" because the gate text implies a response on a field
  beyond the accumulator.
- Item 2 (saturation): no current open-loop equivalent.
  Saturation form is a hard blocker.
- Item 4 (bounded feedback): no closed-loop loop exists.
- Item 5 (sensitivity sweep): no current sweep test exists. Phase
  C will add one; not in scope for this Phase A evidence ledger.

These items are **explicitly not satisfied** by this note. The
phased-plan brief covers the resolution path for each.

---

## Caveat on supersession

If the design-discussion lock artifact changes:
- the phase split (e.g., merges Phase A into another phase, or
  reassigns Item 6 to a later phase),
- the evidence requirements (e.g., demands a per-function Sanity
  Gate proof per item rather than a docstring reference),
- the gate items themselves (e.g., adds a 7th item, drops Item 5,
  rewords any item),

then this note is superseded by the lock artifact and should be
either rewritten to match the new phase split or deleted in favour
of the lock artifact's per-phase evidence sections.

This note's commitment is solely "as of `923d0b2`, the open-loop
side of Items 3 and 6 has the evidence cited above". Any further
interpretation requires the design-discussion lock.

---

## References

- Phased-plan brief (in flight to design-discussion):
  `docs/v2_closed_loop_ecm_gate_dependency_brief.md` (`923d0b2`)
- ECM-OL preflight module:
  `acs/v2/dynamics/ecm_open_loop.py` (`22b9456` post)
- ECM-OL harness:
  `acs/v2/ecm_open_loop_harness.py`,
  `tests/test_v2_ecm_ol_harness.py` (`04ee5a7`)
- Forward roadmap:
  `docs/v2_phase1_forward_roadmap.md` lines 63–77
- Hard Rule 10 (CLAUDE.md): dimensional comparison verification.
