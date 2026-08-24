# ALEPH-PORT-3656 — the filament length is drawn, not assigned

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3656` |
| Lane | `a40e55a2 S-OBSERVE handover` |
| Status | `PROPOSED` |
| Written | `2026-08-09` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |
| Completes | `ALEPH-PORT-3652` (the mean) and `ALEPH-PORT-3655` (why a floor on it was wrong) |
| Asked for by | the PI, uuid `5603243b-f9c6-4981-9d8d-c56044d9716d` — *"필라멘트 길이 분포 역시 데이터들이 있었던것 같은데 항상 같을 수 없으니"* |

---

## 1. Aleph API

`aleph/vertical/cortex_sourcing.py`:

- `FilamentLengthLaw` — `UNIFORM` (every filament the same length, the historical behaviour) and
  `EXPONENTIAL` (the sourced law).
- `draw_filament_half_lengths_um(count, *, law, mean_full_length_um, seed) -> np.ndarray`.

**The construction is NOT wired to it in this entry, and §1 said otherwise until it was corrected.**
`cortex_surface_coupling.py` still applies one scalar half-length to every filament. Threading a
per-filament array through the areal-density path changes crosslink counts, connectivity and node
positions on a real cortex, and §13 already required a measured before/after for exactly that — so
the wiring belongs with the measurement and not with the sampler.

**And when it is wired, the retained one-per-vertex control keeps a scalar**: it is a reference
construction, and moving what a reference references is worse than leaving it unsupported.

**Default is `UNIFORM`, so nothing moves unless a caller asks.** That is deliberate and it is the
lesson of `ALEPH-PORT-3652` §13, which this session wrote and then broke within the hour.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none.** Nothing read from `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | **not applicable** — no provider commit is behind this |
| Source path and symbol | **not applicable** — no `ffn_sim` path and no provider symbol was read |
| Read from | Fritzsche M., Erlenkämper C., Moeendarbary E., Charras G., Kruse K. (2016) *Actin kinetics shapes cortical network structure and mechanics*, **Sci Adv 2:e1501337**, read in full |
| Working tree == commit? | **not applicable** |

## 3. Why source-derived porting beats clean-room

It does not; nothing is source-derived. Stated rather than left blank.

## 4. Physical or mathematical law represented

§Results: *"The actin turnover processes considered above do imply an **exponential distribution**
for the length of formin- and Arp2/3-nucleated filaments."* Fig 3A plots `P(L)` for HeLa and M2,
finite-cortex simulations against analytic approximations.

```
P(L) = (1/L̄) · exp(−L / L̄)          L̄ = 0.120 µm (Arp2/3, HeLa) or 0.060 µm (M2)
```

**The law is not decorative.** An exponential has its mode at zero: the *most common* filament is
short, and the mean is carried by a tail. A construction that gives every filament `L̄` has the
right first moment and **the wrong network** — it has no short filaments that fail to reach a
neighbour and no long ones that reach several.

## 5. Units, domains, singular cases, invariants

- µm. `mean_full_length_um > 0`, `count >= 1`, `seed` an int.
- **I1.** Drawn lengths are strictly positive. An exponential can return arbitrarily small values
  and **zero is excluded**, because a zero-length filament has no tangent — `-3655` kept that domain
  when it removed the floor.
- **I2.** Determinism: the same `(count, mean, seed)` gives a bit-identical array. The seed is a
  **named argument**, following `placement_seed`'s precedent in the same module, not a hidden
  default.
- **I3.** `UNIFORM` reproduces the historical construction **bit-identically**.
- **I4.** The sample mean converges to `L̄`; at `count = 31,165` (the native density) the standard
  error is `L̄/√N ≈ 0.57 %`.

## 6. Source evidence class and known retractions

Peer-reviewed, open access; the Arp2/3 figure is stated in the paper to be independently consistent
with electron tomography. No retraction known. **This entry does not re-source the mean** — that is
`-3652`, and this one adds only the shape around it.

## 7. Independent oracle or derivation

**Three, none of which needs the engine.**

1. **The drawn population's mean is the sourced mean**, to the standard error at the count used.
2. **The fraction below one crossover repeat is `1 − e^{−0.036/L̄}` = 25.9 % at `L̄` = 120 nm** —
   the same arithmetic `-3655` used to retract the floor, now produced by the sampler rather than
   asserted about it. **If the sampler and that number disagree, one of the two entries is wrong.**
3. **`UNIFORM` against the pre-change construction, bit for bit.** A new option that silently moves
   the old path is the defect, not the feature.

## 8. Positive control

`tests/vertical/test_filament_length_law.py`:

- `test_uniform_reproduces_a_constant_length` — every entry equal, and equal to `L̄/2`.
- `test_the_exponential_mean_converges_to_the_sourced_mean` — 31,165 draws, within 3 standard
  errors.
- `test_the_draw_is_deterministic_in_its_seed` — same seed bit-identical, different seed different.
- `test_the_fraction_below_one_crossover_repeat_matches_ALEPH_PORT_3655` — the sampler reproduces
  25.9 % to within sampling error, tying the two entries together.

## 9. Deliberately failing negative control

- `test_no_drawn_length_is_zero` — I1 driven over a large sample; a zero would be a filament with no
  tangent, which is the domain `-3655` kept.
- `test_a_non_positive_mean_is_refused` and `test_a_zero_count_is_refused`.
- `test_uniform_and_exponential_differ` — a vacuity control: if the law argument changed nothing the
  three controls above would pass against a sampler that ignored it.
- `test_the_exponential_is_not_a_narrow_band` — asserts the drawn spread is at least half the mean,
  because an exponential whose samples cluster at `L̄` is not an exponential and would pass the
  mean control.

## 10. Numerical and precision envelope

`numpy.random.Generator.exponential`, float64. Comparisons on the mean use 3 standard errors, which
is a stated tolerance rather than an arbitrary one. `UNIFORM` is compared with `np.array_equal`.

## 11. Production-backend residency and transfer

**None at the sampler.** Downstream, per-filament lengths change node **positions** and therefore
crosslink counts — array *sizes* move with the crosslink count, which the existing host/device
parity controls already cover, and dtypes and layout do not change.

## 12. Comments and docstrings to discard

Nothing copied. The paper is quoted verbatim in §4 and cited.

## 13. Acceptance

`PROPOSED`. `AUDITED` when §8 and §9 pass. **`ACCEPTED` requires a measured before/after on a real
cortex** — crosslink count, connectivity, and the fraction of filaments that reach nothing — and
that is a separate step, for the same reason `-3652` did not move a default inside the entry that
sourced its value.

## 14. Honest limits

- **`EXPONENTIAL` is not the default.** Nothing changes until a caller asks. The blast radius is
  therefore zero and also the benefit is zero until someone opts in; §13 is what closes that.
- **One population, still.** The cortex is Arp2/3 *and* formin, 120 nm and 1,200 nm, and this draws
  from one of them. A two-population mixture is the honest object and is not built here.
- **The retained one-per-vertex control is untouched**, so it keeps a construction the source does
  not support. It is a reference, and moving a reference silently is worse.
- **No length correlates with position.** Real cortical filaments are turned over locally, so length
  and neighbourhood are not independent. This draw treats them as independent and says so.
