# ALEPH-PORT-3655 — the crossover-repeat floor was on the wrong object

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3655` |
| Lane | `a40e55a2 S-OBSERVE handover` |
| Status | `PROPOSED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Corrects | **`ALEPH-PORT-3650` §4/§5 I2 and its `assert_cortex_is_sourced` filament check.** That entry is this session's own, landed 2026-08-08 as `6a57b2d`. |
| Found by | the PI, uuid `5603243b-f9c6-4981-9d8d-c56044d9716d`, `2026-08-08T11:17:58.703Z`: *"필라멘트 길이 분포 역시 데이터들이 있었던것 같은데 항상 같을 수 없으니"* |

---

## 1. Aleph API

`aleph/vertical/cortex_sourcing.py`:

- `assert_cortex_is_sourced` **stops refusing on filament length.** The crosslinker checks stay.
- `minimum_filament_length_um()` is **retired** as a construction floor and replaced by
  `can_host_a_crosslinker(filament_length_um) -> bool` — the same 36 nm, asked of the right object.
- `ACTIN_CROSSOVER_REPEAT_UM` is kept unchanged. **The measurement was never wrong; its placement
  was.**

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** Nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | **not applicable** |
| Source path | **not applicable** |
| Source symbol(s) | **not applicable**; §1 names the Aleph symbols |
| Read from | Fritzsche et al. 2016 *Sci Adv* 2:e1501337, **read in full** — the same paper `-3652` sourced the mean from, and the same section |
| Working tree == commit? | **not applicable** |

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. Stated rather than left blank.

## 4. Physical or mathematical law represented

**What `-3650` got right:** α-actinin's bound length is 36 nm, that is also the F-actin two-start
crossover repeat, and it is also the minimal spacing between crosslinkers along a filament. So a
filament shorter than one repeat **cannot present a crosslinker's binding geometry**.

**What `-3650` got wrong:** it turned that into a floor on **whether a filament may be built**.
Filament length is a **distribution**, not a value — Fritzsche §Results, in the same section the
120 nm came from: *"The actin turnover processes considered above do imply an **exponential
distribution** for the length of formin- and Arp2/3-nucleated filaments"*, with Fig 3A plotting
`P(L)` for both cell lines. Under `P(L) = (1/L̄)e^{−L/L̄}`:

| population | `L̄` | `P(L < 36 nm)` | `P(L < 10 nm)` |
|---|---:|---:|---:|
| Arp2/3, HeLa | 120 nm | **25.9 %** | 8.0 % |
| Arp2/3, M2 | 60 nm | **45.1 %** | 15.4 % |
| formin, HeLa | 1,200 nm | 3.0 % | 0.8 % |

**A quarter to a half of the real population is shorter than one crossover repeat.** Those filaments
are in cells. What they cannot do is host a crosslinker.

**So the constraint moves to the crosslink, and there it is already enforced.** A crosslink exists
iff two filaments come within the crosslinker's reach; a filament too short to present a binding
site simply never satisfies it. **The correct rule needs no floor at all** — the reach criterion
`ALEPH-PORT-3650` itself installed does the work, and `can_host_a_crosslinker` exists to *report*
the fact rather than to gate a build.

## 5. Units, domains, singular cases, invariants

- µm. `ACTIN_CROSSOVER_REPEAT_UM` unchanged at 0.036.
- **I1.** `assert_cortex_is_sourced` accepts any positive filament length. It still refuses a
  missing crosslinker and two disagreeing reaches.
- **I2 (replaces `-3650`'s I2).** `can_host_a_crosslinker(L) == (L >= 0.036)`, and it **gates
  nothing** — it is a predicate a caller may report on, not a build-time refusal.
- **I3.** Under the sourced distribution the fraction of a real population that cannot host a
  crosslinker is `1 − e^{−0.036/L̄}`, which is 25.9 % at `L̄ = 120 nm`. A construction that refuses
  those filaments cannot represent a cortex.

## 6. Source evidence class and known retractions

Peer-reviewed, open access, single-molecule imaging with stochastic simulation, the Arp2/3 figure
independently consistent with electron tomography. **This entry is itself a retraction**, of
`-3650`'s §5 I2 and of the filament branch of its `assert_cortex_is_sourced`.

## 7. Independent oracle or derivation

**Arithmetic on the sourced distribution, with no code in the loop.** `1 − e^{−0.036/0.120}` =
0.2592, and `1 − e^{−0.036/0.060}` = 0.4512. Any construction sampling `P(L)` and passing the old
guard would have to reject that fraction of its own draws — so the old guard and the sourced
distribution cannot both stand, and the distribution is the measurement.

**And a second one that needs no distribution at all:** `-3652` set the *mean* at 120 nm and the old
floor at 36 nm. A guard whose threshold is 30 % of the mean of an exponential variable rejects a
quarter of it by construction. The two entries were inconsistent from the moment `-3652` landed and
nobody noticed for a day.

## 8. Positive control

`tests/vertical/test_cortex_sourcing.py`:

- `test_a_short_filament_no_longer_refuses_the_cortex` — 10 nm passes with a sourced crosslinker,
  where `-3650` refused it.
- `test_the_crosslinker_checks_are_unchanged` — no crosslinker still refuses; two disagreeing
  reaches still refuse. **The correction removes one check and not the entry.**
- `test_can_host_a_crosslinker_is_the_crossover_repeat` — read from the band, not a literal.
- `test_the_predicate_gates_nothing` — calling `assert_cortex_is_sourced` with a filament below the
  repeat returns `None`, so the predicate is reporting and not refusing.

## 9. Deliberately failing negative control

- `test_the_old_floor_would_reject_a_quarter_of_a_real_population` — computes
  `1 − e^{−0.036/L̄}` for both sourced means and asserts 25.9 % and 45.1 %, so the reason for this
  correction is pinned and cannot quietly stop being true.
- `test_a_filament_of_zero_length_is_still_refused` — removing the floor must not remove the domain:
  a non-positive length has no tangent and no axial law.
- `test_can_host_a_crosslinker_can_return_both` — vacuity control on the predicate.

## 10. Numerical and precision envelope

Exact float64 comparisons against declared constants; the distribution fractions are asserted to
three decimals, which is far inside `numpy`'s `exp`.

## 11. Production-backend residency and transfer

**None.** A build-time host check being removed, and a pure predicate added.

## 12. Comments and docstrings to discard

Nothing copied. `-3650`'s own refusal message is rewritten rather than deleted, because it argued
the floor and that argument is what is being corrected.

## 13. Acceptance

`PROPOSED`. `AUDITED` when §8 and §9 pass. **`-3650` must be amended in the same commit** — an entry
whose §5 I2 has been retracted may not read as though it still holds.

## 14. Honest limits

- **This removes a check; it does not add the sampling.** The construction still builds one length
  for every filament. Drawing from `P(L)` is the next entry, and until it lands the *mean* is what
  is built and the distribution is only what justifies not refusing short ones.
- **`can_host_a_crosslinker` is reported by nobody yet.** It is the right home for the constraint and
  currently an unused predicate; wiring it into a build report is a separate, small entry.
- **It does not revisit `-3650`'s other checks.** The missing-crosslinker refusal and the
  disagreeing-reaches refusal stand and are re-controlled here.
- **The error was mine and it was live for one day**, on `main`, behind controls that passed. It was
  found by the PI and not by this session or by any guard.
