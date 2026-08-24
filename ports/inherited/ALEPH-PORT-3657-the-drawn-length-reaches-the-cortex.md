# ALEPH-PORT-3657 — the drawn length reaches the cortex

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3657` |
| Lane | `b4fad06b` (construction layer, picked up from the closed `a40e55a2`) |
| Status | `PROPOSED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Completes | `ALEPH-PORT-3656` §13 — the sampler exists and is wired to nothing; this is the wiring **and** the before/after that `-3656` made a precondition of `ACCEPTED` |
| Does not re-source | the mean is `-3652`'s and the distribution is `-3656`'s. This entry adds no new number |

---

## 1. Aleph API

`aleph/vertical/cortex_surface_coupling.py`:

- `build_radial_cortex_network(..., filament_length_law=FilamentLengthLaw.UNIFORM,
  filament_length_seed=20260809)`.
- The law is consumed **only on the areal-density branch**. On the retained one-filament-per-vertex
  branch a non-`UNIFORM` law is **refused by name**, not ignored — see §5, I5.

`aleph/vertical/assembly.py`:

- `build_vertical(..., cortex_filament_length_law=..., cortex_filament_length_seed=...)`, forwarded
  verbatim. Both names carry the `cortex_` prefix, so `build_whole_cell(cortex={...})` reaches them
  through `_validated_cortex_parameters`'s existing allowlist with no change to that allowlist.

**What this does NOT do.** It does not change a default, it does not move the filament mean, and it
does not touch `assert_cortex_is_sourced`. `UNIFORM` remains the default at every level of the call
chain, so a caller who asks for nothing gets the construction that every committed cortex number in
this repository was measured on.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** Nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | **not applicable** — no provider commit is behind this |
| Source path and symbol | **not applicable** — no `ffn_sim` path and no provider symbol was read |
| Read from | the distribution is `ALEPH-PORT-3656`'s, sourced there from Fritzsche M., Erlenkämper C., Moeendarbary E., Charras G., Kruse K. (2016) *Actin kinetics shapes cortical network structure and mechanics*, **Sci Adv 2:e1501337** |
| Working tree == commit? | **not applicable** |

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. Stated rather than left blank.

## 4. Physical or mathematical law represented

No new law. The law is `-3656`'s `P(L) = (1/L̄)·exp(−L/L̄)`; this entry makes the cortex builder
consume a **per-filament** half-length array where it previously broadcast one scalar.

**What the wiring changes physically, and it is not cosmetic.** Two filaments form a crosslink when
their segments come within the crosslinker's reach `c`. Two segments of half-lengths `h_a`, `h_b`
whose centres are `d` apart can only do so if

```
d  <=  c + h_a + h_b
```

— the collinear worst case, and an upper bound on the real segment–segment criterion. On a surface
of areal density `ρ` the expected number of partners a filament has is therefore proportional to
`ρ · E[(c + h_a + h_b)²]`, i.e. to the **second** moment of a length distribution whose **first**
moment the uniform construction already gets right.

## 5. Units, domains, singular cases, invariants

µm throughout. The law is a `FilamentLengthLaw`; the seed is an int.

- **I1 — `UNIFORM` is bit-identical to the pre-change builder.** `np.array_equal` on node
  positions, and the same crosslink count, ids and rest lengths. A new option that silently moves
  the old path is the defect, not the feature.
- **I2 — determinism.** The same `(count, mean, seed)` gives the same cortex to the last bit. The
  seed is a **named argument**, following `placement_seed`'s precedent in the same module.
- **I3 — the length seed is independent of `placement_seed`.** Changing one must not move the other's
  output, or the two randomisations are one randomisation wearing two names.
- **I4 — the drawn population's mean half-length is the scalar it replaces**, to the sample standard
  error. The wiring may not move the first moment; that is `-3652`'s number and this entry does not
  re-source it.
- **I5 — a non-`UNIFORM` law on the one-filament-per-vertex branch is REFUSED.** That branch is the
  **retained reference control** (`-3656` §1: *"moving what a reference references is worse than
  leaving it unsupported"*). Silently ignoring an argument the caller supplied is the same class of
  defect as the geometric reach fallback that `-3650` removed: the call site reads as if a law were
  applied and no message says otherwise.

## 6. Source evidence class and known retractions

Peer-reviewed, open access. No retraction known. This entry adds **no** evidence class of its own —
it is a wiring and a measurement, and every number it reports is a property of this tree, marked as
such and never as physiology.

## 7. Independent oracle or derivation

**Written before the code runs, so the outcome can refute it.** All three follow from §4 with no
appeal to the engine.

**(a) — WRITTEN, THEN REFUTED BY ITS OWN MEASUREMENT. The retraction is left in place with the
claim it replaces, because a prediction that is quietly deleted after it fails is not a prediction.**

*As written, before the code ran:* with `S = h_a + h_b` and `r = c + S`, the reachable area goes as
`E[r²]`; for `EXPONENTIAL`, `S` is Erlang-2 with the same mean and `Var(S) = 2h̄²`, so
`E[r²]/(c+2h̄)² = 1 + 2h̄²/(c+2h̄)² = 1.296` at `c = 0.036`, `2h̄ = 0.120`. **≈ +30 % crosslinks at a
fixed mean.** *"A measured excess above 30 %, or any deficit, refutes this section rather than the
code."*

*Measured, 8 length seeds, ρ = 100 /µm², 1 µm shell, α-actinin:*

| | uniform | exponential |
|---|---:|---:|
| crosslinks | 1,703 | 1,651.5 ± 72.8 (sd) |
| ratio to uniform | 1 | **0.9698 ± 0.0428**, range 0.924–1.038 |

**A deficit, and one indistinguishable from no change at all. §7(a) as written is refuted.**

*Why it was wrong, which is the part worth keeping.* `d ≤ c + h_a + h_b` is the **collinear** case:
a one-dimensional statement about the axial direction. Squaring it to make an area counts the
length a second time in a direction where a rod has **no extent**, and that double count is what
manufactured the variance term. The correct object is the two-dimensional excluded area of two rods
of half-lengths `h_a`, `h_b` at angle `θ` with capture distance `c`:

```
A_excl(θ)  =  4·h_a·h_b·|sin θ|  +  4·c·(h_a + h_b)  +  π·c²
```

Every term is **linear**, or **bilinear in two independent draws**. The expected crosslink count is
`(ρ²A/2)·⟨A_excl⟩`, and with `E[h_a·h_b] = E[h_a]·E[h_b]` for independent lengths, that is a
function of the **first moment alone**. **No distribution at a fixed mean may move it**, which is
what was measured.

*The general lesson, and it is the one this lane was warned about.* A quantity built from a
**product of two independent** draws is mean-preserving; a quantity built from the **square of one**
is not. Writing the pair reach as `(h_a + h_b)²` silently converted the first into the second. The
closed row's trap (i) — *"a single number was trusted twice"* — recurred here as a single number
squared once too often, and it was caught by the control rather than by the author.

**(b) The fraction of filaments that reach NOTHING must RISE. — CONFIRMED, and it survives (a)'s
retraction untouched**, because it rests on the mode and not on the second moment. An exponential's
mode is at zero, so the population gains filaments whose own half-length contributes nothing to any
pair's reach; a filament at `h_a → 0` keeps only the `4c·h_b + πc²` terms of `A_excl` and loses the
bilinear one entirely.

*Measured, same 8 seeds:* uniform **0.0000** — on a Fibonacci lattice at a uniform length, every
single filament reaches a partner — against exponential **5.05 % – 8.37 %**.

**This is the whole result, and it is sharper than the one that was predicted.** The distribution
does **not** change how many crosslinks there are. It changes **which filaments have them**: the same
1,700 crosslinks redistribute onto the long tail, and 7 % of the population is left touching
nothing at all. A measurement that reported only a crosslink count would have seen **no effect
whatsoever** and concluded the law was cosmetic.

**(c) The mean is untouched.** The measured mean half-length of the drawn population equals the
scalar it replaced to within `h̄/√N` — 0.57 % at the native-density count. This is I4 and it is what
distinguishes "the shape changed" from "the cortex got more actin".

## 8. Positive control

`tests/vertical/test_filament_length_law.py` (the file `-3656` created; these are added to it):

- `test_uniform_through_the_builder_is_bit_identical` — I1, on node positions, crosslink count and
  crosslink rest lengths, against a build that does not pass the argument at all.
- `test_the_drawn_lengths_reach_the_built_filaments` — the built network's per-filament rest lengths
  are the drawn array, not the scalar. Without this the other controls could pass against a builder
  that accepted the argument and dropped it.
- `test_the_law_is_deterministic_in_its_seed_through_the_builder` — I2.
- `test_the_length_seed_is_independent_of_the_placement_seed` — I3.
- `test_the_drawn_mean_is_the_scalar_it_replaced` — I4.
- `test_the_crosslink_count_is_a_first_moment_and_does_not_move` — §7(a) **as corrected**. A
  two-sided ±15 % band around the uniform count, which the 8 measured seeds fit inside with margin
  and which **the refuted +30 % would fail.** The retracted claim is kept as a thing the control can
  catch, not deleted.
- `test_the_exponential_leaves_more_filaments_unreached` — §7(b). 0.0000 → 0.0505–0.0837.

## 9. Deliberately failing negative control

- `test_a_length_law_on_the_one_per_vertex_control_is_refused` — I5. The refusal must name the
  argument and the branch, so a caller can tell an ignored argument from an inapplicable one.
- `test_an_unknown_length_law_is_refused_by_the_builder` — the refusal comes from
  `draw_filament_half_lengths_um` and reaches the caller rather than being swallowed.
- `test_the_law_argument_is_not_vacuous` — `UNIFORM` and `EXPONENTIAL` at the same seed produce
  different node positions. Without it every control above would pass against a builder that
  ignored the law entirely, which is exactly the failure I5 exists to prevent elsewhere.

## 10. Numerical and precision envelope

float64 throughout. `UNIFORM` is compared with `np.array_equal` — an exact comparison, because the
claim is bit-identity and not agreement. The crosslink-count comparisons in §8 are **inequalities
with a predicted bound**, not tolerances, because §7 predicts a direction and a ceiling rather than
a value.

**One ordering fact that is load-bearing.** `_closest_approach_crosslinks` returns pairs in ascending
`(index_a, index_b)` order and that order fixes the sequence the crosslink force accumulates in;
floating-point addition is not associative. Per-filament lengths change **which** pairs are within
reach, so the crosslink set changes — but the ordering rule is untouched, so a rebuild at the same
seed is reproducible to the last bit.

## 11. Production-backend residency and transfer

**No new residency.** Per-filament lengths change node *positions* and therefore the crosslink
count; array **sizes** move with that count, which the existing host/device parity controls already
cover. No dtype, no layout and no kernel changes. The device kernels read positions and a crosslink
table and neither gains a field.

## 12. Comments and docstrings to discard

Nothing copied. `build_radial_cortex_network`'s docstring gains one paragraph stating that the
half-length is a population mean once a law is given, and the `ρ_crit = 1/(2·half)²` note in
`resolve_cross_filament_overlaps` is qualified: under a law that is a mean, not a bound, so the
overlap threshold becomes a statement about the average filament.

## 13. Acceptance

`AUDITED` — §8 and §9 pass, `tests/vertical/test_filament_length_law.py` **22 passed**. The
before/after on a real cortex is measured and written up:
`docs/results/2026-08-09-filament-length-law/` (31,165 filaments, 5 builds, `measured.json`).

**§7's predictions, checked against it in writing, including the one that failed:**

| | predicted | measured | |
|---|---|---|---|
| (a) crosslink count | **+30 %** | **0.9990 ± 0.0032** | **REFUTED**, and the derivation is corrected above |
| (b) filaments reaching nothing | rises | **0.0000 % → 6.47–6.72 %** | confirmed |
| (c) mean unmoved | within 0.57 % | 0.11915–0.12016 vs 0.120 | confirmed |

**One result nobody predicted, and it is the largest:** the **giant component falls from 100 % to
79 %**. The same ~42,800 crosslinks, redistributed onto the long tail, leave a fifth of the cortex
off the main network. That is not in §7 at all — it was found because the measurement reported three
numbers instead of one.

**`ACCEPTED` is withheld and the reason is not a missing measurement.** Whether a cortex that is not
singly connected is *correct* is an empirical question this repository has read no source for, and
`ACCEPTED` on this entry would be read as endorsing it. The wiring is audited; the modelling
consequence is reported. See the results README §4.

**The default does not move in this entry.** `-3652` §13 said a default may not move inside the
entry that sources its value; the same rule applies to the entry that wires a law. `EXPONENTIAL`
becoming the default is a separate decision with its own measured blast radius, and it is behind
`ALEPH-PORT-3652`'s `cortex_filament_half_length_um` question, which is itself behind `PLAN.md`
item #2.

## 14. Honest limits

- **One population, still.** Arp2/3 at 120 nm; the formin population at 1,200 nm is not drawn. A
  two-population mixture is the honest object and is not built here.
- **No length–position correlation.** Real cortical filaments turn over locally, so length and
  neighbourhood are not independent. This wiring draws them independently and says so.
- **The card is still a scalar.** `build_vertical` derives `CortexCard` from
  `filament_length_um = 2 × cortex_filament_half_length_um`, which under a law is the population
  **mean**. That derivation is already flagged in its own comment as a discretisation scaling and
  not a biological claim; a per-filament card is not built here and would be a different entry.
- **The retained one-per-vertex control keeps a scalar**, by refusal rather than by omission (I5).
- **The overlap-resolution threshold `ρ_crit = 1/(2·half)²` becomes an average.** Under a law some
  filaments are far longer than the mean and overlap below the nominal threshold. The sweep count
  default is `0`, so nothing silently changes; it is recorded because a caller who turns sweeps on
  is now resolving a population the threshold no longer describes exactly.
