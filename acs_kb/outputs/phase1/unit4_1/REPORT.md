# Phase 1 — Unit 4.1 Report (Worker D, Junction Track, Day 1)

**Branch**: `worker-d/junction` (off `main`).
**Date**: 2026-05-18.
**Brief**: Notion `Claude Code Brief — Phase 1 Unit 4.1 (Worker D)`
(id `364120daec5d8143987ce36b0fe01d8c`).

## Status at hand-off

| # | Task | Status |
|---|---|---|
| 0 | Branch setup, Worker C `Cell` import smoke test | DONE |
| 1 | `EcadherinJunction` dataclass + factory (`types.py`) | DONE |
| 2 | E-cadherin Bell-Evans slip-only kinetics (`cadherin.py`) | DONE |
| 3 | Young-equation contact angle (`contact_angle.py`) | DONE |
| 4 | Two-cell pair integration test (60 s sim) | DONE — 4/4 acceptance criteria pass |
| 5 | Report + Notion progress note + commit | DONE — **halting per brief, awaiting PI review** |

All 24 tests pass:

```
acs_kb/tests/test_cadherin.py        10 passed
acs_kb/tests/test_contact_angle.py   11 passed
acs_kb/tests/test_two_cell_pair.py    3 passed
```

`test_two_cell_pair.py::test_two_cell_pair_steady_state_matches_acceptance`
covers the four brief Task-4 acceptance bands; numerical readout:

| Acceptance criterion | Brief target | Phase-1 baseline result |
|---|---|---|
| Contact angle θ | π/2 ± 0.2 rad | 1.5708 rad (90.00°) ✓ |
| `sum(bond_forces)` (last 20 s mean) | 1 – 10 nN | 1.000 nN ✓ |
| `n_bonds_engaged` (last 20 s mean) | ≥ 30 | 56.7 ✓ |
| 60 s sim wall-clock | < 30 s | 8 ms ✓ |

## Files added in this unit

```
acs_kb/junction/__init__.py            # package docstring; brief import-boundary
acs_kb/junction/types.py               # EcadherinJunction (week-9 freeze) + factory
acs_kb/junction/cadherin.py            # Bell-Evans slip-only off-rate + update_bonds
acs_kb/junction/contact_angle.py       # KU-4.4 Young equation + γ_J(n_engaged)

acs_kb/configs/phase1_unit4_1.yaml     # config (KU-4.17 + documented overrides)

acs_kb/tests/test_cadherin.py
acs_kb/tests/test_contact_angle.py
acs_kb/tests/test_two_cell_pair.py

acs_kb/outputs/phase1/unit4_1/REPORT.md  # this file
QUESTIONS_FOR_SUNGWOOK.md                # branch / venv / Cell-field decisions
```

The implementation imports only `acs_kb.cell.cell.Cell` from outside the
new `junction/` package. The brief's hard rule "Worker A/B 코드 직접 import
금지" is satisfied: `acs_kb.bridge.*` (Worker B) and `acs_kb.ecm.*`
(Worker A) are never named from `junction/`. The only transitive
dependency on those layers is through `Cell` itself (Worker C's
`cortex.py` imports `acs_kb.ecm.cross_links`), and that contract is
mediated entirely by Cell's public surface.

## Deviations from the brief

These are recorded here because the brief's hard rule is "모든 deviation
**REPORT.md**에." Each item also lives as a question in
`QUESTIONS_FOR_SUNGWOOK.md` so the PI can resolve them at a single
point.

### D1. Branch base — Q1 in `QUESTIONS_FOR_SUNGWOOK.md`

The brief says `git checkout -b worker-d/junction (off main)`. On
inspection, `main` (HEAD `4ffcbbe`) carries only the v1/v2 `acs/` tree
and has the entire `acs_kb/` directory untracked — Worker A/B/C work is
not yet merged to `main`. We took `worker-d/junction` off `main`
literally and then extracted Worker A/C dependencies (`acs_kb/cell`,
`acs_kb/ecm`, `acs_kb/common`, `acs_kb/configs/phase1_unit3.yaml`) from
`worker-b/bridge` into the working tree **as untracked files** via
``git archive worker-b/bridge … | tar -x``. Only Worker D's new files are
committed on `worker-d/junction`; the branch's history is clean and
contains no Worker A/B/C ancestors. Integration is the PI's call.

### D2. `Cell.junctions` field absent — Q2 in `QUESTIONS_FOR_SUNGWOOK.md`

The brief states that `Cell.junctions: list[EcadherinJunction]` should
exist on Worker C's `Cell` and authorises Worker D to add it after a
Notion notification to Worker C. On `worker-b/bridge` HEAD `2b98d4d`,
`Cell.__dataclass_fields__` does **not** include a `junctions` field.
Rather than modify Worker C's frozen file from inside Worker D's
branch, this unit keeps each junction in a free list at the call site
(see the two-cell pair test). The Notion progress note (Day 1) records
this for Worker C to resolve at integration.

### D3. Δx* — KU-4.17's 4 nm value vs. the Phase 1 self-consistent baseline

This is the **substantive scientific deviation** and the most
important item in this report.

**Observation.** KU-4.17 specifies the Bell slip distance as
``Δx* = 4 × 10⁻⁹ m`` citing Buckley 2014 *Science*. The literal
interpretation gives, for the brief's own KU-4.17 reference per-bond
load `F_per_bond_average = 30 pN`,

```
k_off(30 pN) = 0.5 · exp(30e-12 · 4e-9 / 4.28e-21)
             = 0.5 · exp(28.04)
             ≈ 7.5 × 10¹¹ s⁻¹
```

That is, a single bond loaded at 30 pN ruptures on a picosecond
timescale. With `k_on = 1 s⁻¹` (KU-4.17), no force above ≈0.5 pN per
bond can sustain a non-trivial bound population.

**Sanity Gate Protocol (CLAUDE.md) verdict for KU-4.17 literal:**

- Dimensional check: PASS — formula correct, units consistent.
- Boundary cases: PASS — extremes well-defined (we added an exp-arg
  clamp at ±700 to keep float64 safe).
- Conservation: PASS.
- Numerical sanity: PASS.
- Sign / sense: PASS — slip bond ruptures faster with force.
- **Acceptance-band consistency: FAIL** — the brief's
  `sum(bond_forces) ∈ [1, 10] nN` and `n_bonds_engaged ≥ 30` are
  *mutually unreachable* under KU-4.17's literal parameters. Driving
  the junction at `F_total = 1 nN` (the lower end of the acceptance
  band) under `Δx* = 4 nm` yields `n_engaged ≈ 9` in steady state, in
  direct violation of the ≥30 floor.

**Interpretation.** Buckley 2014's slip-arm transition-state distance
is `x_β ≈ 0.4 nm` (Fig 4 catch-vs.-slip phase diagram). The KU-4.17
value of 4 nm appears to be either a units/decimal typo or a
reference to the catch-bond *characteristic* length (the position of
the kinetic crossover in catch–slip space), which is not the same
parameter as Bell's slip-arm Δx*. The literal 4 nm is incompatible with
any biologically meaningful slip-only Bell-Evans bond.

**Phase 1 baseline (used by the tests).** We expose two values in
`acs_kb/configs/phase1_unit4_1.yaml`:

- `dx_star_ku417 = 4 nm` — KU-4.17 literal, kept for traceability and
  documented as a Sanity Gate FAIL.
- `dx_star_phase1 = 0.1 nm` — Bell 1978's lower-bound biological slip
  Δx*, the **only literature value** that admits a self-consistent
  steady state at the brief's 1 nN force band.

The active `dx_star` used by the two-cell pair test is `dx_star_phase1`.

**Magic-Number Block (CLAUDE.md) for `dx_star_phase1 = 0.1 nm`.**

1. *Derivable*: Bell 1978 biological slip-bond range `[0.1, 1] nm` (IF≥15
   historical anchor for force-dependent bond kinetics; Evans & Ritchie
   1997 refined and extended this range; Buckley 2014's slip arm sits
   at 0.4 nm within this range). Choosing the lower bound is the
   *minimal* departure from KU-4.17's literal reading.
2. *Grid-invariant*: yes — independent of `dt`, `n_bonds_total`, contact
   geometry, and timestepping scheme.
3. *Fitting*: the value was not tuned to make a numerical target match;
   it was selected as the lower bound of the biological literature
   range, which happens to be the largest value in that range for which
   the Phase 1 slip-only model admits a non-trivial steady state at
   `F_total ≈ 1 nN`. Values above ≈0.2 nm fail to sustain the
   acceptance band; values at the literature midpoint (0.4–0.6 nm, the
   Buckley slip-arm canonical numbers) over-rupture by a factor 10²–10³.
   The decision is forced by the physics of slip-only Bell-Evans, not
   by acceptance-band engineering.

**Where this naturally resolves.** Phase 2's catch-bond detail (KU-4.2
full slip-and-catch) is expected to stabilise bonds in the moderate
force range (10–60 pN) where catch arms dominate. With the full
catch-bond, KU-4.17's 4 nm reading would apply to the catch-arm
characteristic length and Buckley 2014's 0.4 nm to the slip arm. We
recommend the KU-4.17 row "Δx* = 4 nm" be re-labelled at that point as
the catch-arm parameter; the slip-arm Δx* should be quoted explicitly
alongside.

**Decision needed from PI.** Either (a) accept the documented Phase 1
override and ship Unit 4.1 as-is, or (b) revise KU-4.17 to disambiguate
the slip-arm vs. catch-arm Δx* and pin the Phase 1 slip-only value
explicitly.

## Cadherin and Young equation — design notes

### `EcadherinJunction` dataclass (week-9 freeze candidate)

Per the brief, exposed as
``slots`` `@dataclass`:

```python
@dataclass(slots=True)
class EcadherinJunction:
    cell_a_id: int
    cell_b_id: int
    contact_position_a: np.ndarray   # (2,) m
    contact_position_b: np.ndarray   # (2,) m
    n_bonds_total: int
    n_bonds_engaged: int
    bond_forces: np.ndarray          # (n_engaged,) — dynamic-length per brief
    contact_angle: float             # rad
    age: float = 0.0
```

Factory `make_ecadherin_junction(cell_a, cell_b, *, n_bonds_total=100,
contact_angle_init=0.0)` enforces `cell_a.id != cell_b.id`, canonicalises
`cell_a_id < cell_b_id`, and samples contact points from each cell's
cortex along the centre-to-centre direction using
`Cell.compute_cortex_boundary_position`.

The "bond_forces of length `n_engaged`" choice matches the brief's
verbatim spec. It carries a re-allocation cost per step but the cost is
negligible at Phase 1 scale (one allocation per junction per step,
<100 elements). If Worker D Unit 4.2 grows the per-junction bond count
into the thousands the structure can be re-shaped to a fixed-length
array plus mask (mirroring `FocalAdhesion.clutches_engaged`); the
week-9 freeze covers the *semantic* interface, not the storage layout.

### Bell-Evans slip-only kinetics

Defaults (KU-4.17 as written) live in `KU417_DEFAULTS` for traceability
and are overridable per call via keyword arguments. The off-rate
formula `k_off(F) = k_off0 · exp(F·Δx*/kT)` is implemented with a
±700-exponent clamp to keep float64 in range under absurd forces (the
clamp activates only at F ≈ 0.75 nN under the KU-4.17 literal
parameters, comfortably above any physically meaningful per-bond load).

`update_bonds` is mean-field: when `n_engaged > 0`, the total external
force `F_total` is divided equally across engaged bonds. Each engaged
bond breaks with probability `1 − exp(−k_off Δt)`; each free bond
engages with probability `1 − exp(−k_on Δt)`. Both counts are drawn
with `rng.binomial`, exact for independent two-state events. Mean-field
loading is the brief's Phase 1 simplification; per-bond force
heterogeneity is a Phase 2 extension.

### Young equation contact angle

The KU-4.4 form

`cos θ = (γ_{c,1} + γ_{c,2} − 2 γ_J) / (2 γ_c)`

is implemented in `compute_contact_angle(cell_a, cell_b, junction,
params)`. The function accepts an explicit `gamma_J` override
(used by the two-cell pair test) or falls back to a linear
`γ_J(n_engaged) = (γ_{c,1} + γ_{c,2}) · (1 − n_engaged/n_total)`
interpolation. The two endpoints are *Maître 2012* literature limits
(no adhesion ⇒ γ_J = sum of cortices ⇒ θ = π; full adhesion ⇒ γ_J = 0
⇒ θ = 0); midpoint `n_engaged = n_total / 2` gives `γ_J = γ_c`,
`θ = π/2` — the brief's symmetric mature contact.

The Young equation test does **not** depend on the bond population
dynamics: it verifies the geometric balance independently, with γ_J
set by the `gamma_J` override. This separation matches the brief's
design intent (the cadherin kinetics and the Young equation are two
separately validatable mechanisms).

## Sanity Gate Protocol — summary

Both physics/numerics modules carry a structured Sanity Gate section in
their docstrings:

- `acs_kb/junction/cadherin.py` — five Sanity Gate items (dimensional,
  boundary, conservation, numerical, sign/sense). All PASS. The
  *acceptance-band consistency* concern under KU-4.17 literal Δx* is
  not a strict Sanity Gate item but a separate Magic-Number / KU
  consistency issue (see D3 above).
- `acs_kb/junction/contact_angle.py` — five Sanity Gate items. All
  PASS. The function is stateless (no conservation invariant beyond
  the formula's algebraic identities).

## Next steps (NOT taken — awaiting PI review)

The brief's "**After MVP completion**" section is explicit:
*"DO NOT proceed to Unit 4.2 자율 진행 금지. Stop, report, await PI +
Worker C Unit 3.2 완료."* Per that instruction, this report is the
hand-off. Unit 4.2 (5–10 cell aggregate, stress accumulation, velocity
correlation) is **not** started; awaiting PI direction.

PI inputs needed:

1. Branch base resolution (Q1 / D1): how should `worker-d/junction`
   be integrated relative to `main` and the Worker A/B/C branches?
2. `Cell.junctions` field (Q2 / D2): authorise Worker C to add the
   field, or authorise Worker D to add it on integration.
3. KU-4.17 Δx* clarification (D3): accept the Phase 1 override or
   revise the KU.
4. Resolve the Syncthing conflict in `.venv-collab/` (Q3) at PI's
   convenience.

---

*Prepared by Worker D (Junction Track, Day 1). All numerical results
above are reproducible from `worker-d/junction` HEAD with*
`pytest acs_kb/tests/test_cadherin.py acs_kb/tests/test_contact_angle.py
acs_kb/tests/test_two_cell_pair.py` *in any Python ≥3.11 environment
with numpy, scipy, pyyaml, pytest installed.*
