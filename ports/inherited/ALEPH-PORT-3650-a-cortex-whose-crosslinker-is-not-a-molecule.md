# ALEPH-PORT-3650 — a cortex whose crosslinker is not a molecule must refuse to build

> **PARTIALLY RETRACTED 2026-08-09 by `ALEPH-PORT-3655`. Read that entry before acting on this
> one.** §5 **I2** and the filament-length branch of `assert_cortex_is_sourced` are **withdrawn**:
> they made "shorter than one F-actin crossover repeat" a condition on whether a filament may be
> **built**, and filament length is an exponential distribution in which **25.9 %** of a real
> Arp2/3 population in HeLa — **45.1 %** in melanoma M2 — falls below that. A guard there refuses
> every real cortex.
>
> **The measurement is not retracted.** 36 nm is right, and a filament shorter than it genuinely
> cannot present a crosslinker's binding geometry. It belongs on the **crosslink**, where the reach
> criterion this entry installed already enforces it, and it is now reported by
> `can_host_a_crosslinker` rather than gating a build.
>
> Everything else here stands and is re-controlled by `-3655` §8: the missing-crosslinker refusal,
> the two-disagreeing-reaches refusal, and the density-independence oracle.
>
> **Found by the PI**, not by this session and not by a guard, one day after this entry landed —
> uuid `5603243b-f9c6-4981-9d8d-c56044d9716d`.

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3650` |
| Lane | `a40e55a2 S-OBSERVE handover, scope extended by PI reassignment` |
| Status | `PROPOSED` |
| Written | `2026-08-08` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |

---

## 1. Aleph API

`aleph/vertical/cortex_sourcing.py`, a new module:

- `ACTIN_CROSSOVER_REPEAT_UM: LiteratureRange` — the F-actin two-start helix crossover repeat.
- `ALPHA_ACTININ_LENGTH_UM: LiteratureRange`, `FILAMIN_A_LENGTH_UM: LiteratureRange`.
- `UnsourcedCortexError(ValueError)`.
- `assert_cortex_is_sourced(*, crosslinker, filament_half_length_um, crosslink_capture_um) -> None`
  — refuses, by name and with the number, a cortex that (a) supplies no crosslinker molecule, or
  (b) carries a filament shorter than one crossover repeat, or (c) takes its capture radius from a
  geometric fallback rather than from a molecule.

Called from `aleph/scenarios/whole_cell.py::build_whole_cell`, which is where the silent fallback
actually happens.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | **none — no provider file was read for this entry.** The bounds come from primary literature, cited below |
| Source commit | **not applicable** — nothing was read from `/Users/sw1/ffn_cellsim` (READ ONLY) for this port |
| Source path | **not applicable**; the nearest provider artefact is `ffn_sim/ff/hand_kmc.py:129-132,144-148`, and it is **quoted for refutation** below rather than ported |
| Source symbol(s) | **not applicable**; the Aleph symbols this entry approves are named in §1 |
| Read from | **primary literature**, not a revision of the provider — see the two citations below |
| Working tree == commit? | **not applicable** — no provider revision is cited, so there is none to have drifted |

> The template's last two rows exist because a ledger citing a commit it never read is citing a
> revision that may not exist. The honest answer here is that **no provider revision is cited at
> all**, and writing "not applicable" six times is the way to say so without inventing a sha.

**Not a port of provider code.** Two of the three bounds are re-derived from primary structural
literature, read in this session and recorded in `docs/design/SOURCED_CONSTRUCTION_VALUES.md`:

- Meyer R. & Aebi U. (1990) *J Cell Biol* — α-actinin free 35 nm (negatively stained) / 37 nm
  (rotary-shadowed), **36 nm bound to F-actin**, *"roughly coinciding with the crossover repeat of
  the two-stranded F-actin helix (i.e., 36 nm)"*, minimal crosslinker spacing along a filament
  ≈36 nm.
- Suphamungmee W. et al. (2012) *J Mol Biol* — *"Flexible, 160-nm-long FLNa molecules are
  tail-to-tail dimers"*.

The third input, `Crosslinker.capture_radius_um = 0.060` for **both** presets, is
`ALEPH-PORT-3640`'s transcription of `ffn_sim/ff/hand_kmc.py`. **A finding this entry records
rather than repairs:** that entry's provenance string sources `link_k` to Ferrer 2008 PNAS and
**cites nothing for the capture radius**, and 60 nm matches neither α-actinin's 36 nm nor
filamin's 160 nm — it is one lumped number for two molecules that differ by 4.4×. Routed to
`-3640`, not edited here.

## 3. Why source-derived porting beats clean-room

It does not, and nothing is source-derived here. The provider is read only to record what
`-3640` already transcribed; both usable bounds come from primary papers. This section exists to
say that explicitly rather than leave it blank.

## 4. Physical or mathematical law represented

Not a force law — a **domain of definition**. A crosslinker is a molecule with a contour length ε.
A crosslink exists **iff** the two filaments' closest approach is within ε; beyond it there is no
bond, only a distance. And a filament shorter than one F-actin crossover repeat cannot present the
binding geometry a crosslinker requires, because the crosslinker's own bound length *is* that
repeat and its minimal spacing along a filament is likewise that repeat.

## 5. Units, domains, singular cases, invariants

- `capture_radius_um`, `filament_half_length_um`, all lengths: **µm**.
- Domain: `capture_radius_um ∈ (0, 1)` µm, already enforced by `Crosslinker.__post_init__`.
- Singular case: **no crosslinker at all** — today `cortex_crosslinker=None` silently selects
  `crosslink_capture_from_spacing_um`, i.e. `2.2 × mean nearest-neighbour spacing`. That is the
  case this entry converts from a fallback into a refusal.
- **I1.** A reach may not depend on the population's density. Measured today it does: the median
  crosslink rest is `0.53469 µm` at ρ=8 and `0.15376 µm` at ρ=100 — **1.512× and 1.538× the mean
  spacing `1/√ρ`**, across a 12.5× change in ρ.
- **I2.** `2 × filament_half_length_um ≥ ACTIN_CROSSOVER_REPEAT_UM`. Measured today:
  `2 × 0.005 = 0.010 µm` against `0.036 µm` — **3.6× under**.

## 6. Source evidence class and known retractions

Primary structural-biology measurements, peer-reviewed, EM and cryo-EM. No retraction known for
either. `-3640`'s 60 nm is transcription-of-a-transcription and is marked in §2 as carrying no
primary citation for that field; this entry does **not** use it as a bound, only reports it.

## 7. Independent oracle or derivation

The oracle is **dimensional and independent of the tree**: a reach that is a molecule's property
must not move when the density moves. So build the same cortex at two densities and compare the
reach. Aleph's present construction fails this by 1.512 vs 1.538 × spacing (i.e. the *reach*
tracks `1/√ρ` exactly), while a molecular reach is constant by construction. No provider number is
consulted.

## 8. Positive control

`tests/vertical/test_cortex_sourcing.py`:

- `test_a_sourced_crosslinker_is_accepted` — `crosslinker("ALPHA_ACTININ")` with a filament at or
  above one crossover repeat passes `assert_cortex_is_sourced` and returns `None`.
- `test_the_reach_does_not_move_with_density` — the accepted crosslinker's `capture_radius_um` is
  bit-identical at two densities, where the geometric fallback differs by `√(100/8) = 3.54×`.

## 9. Deliberately failing negative control

- `test_no_crosslinker_is_refused` — `crosslinker=None` raises `UnsourcedCortexError` naming the
  fallback it would otherwise have taken.
- `test_a_filament_shorter_than_one_crossover_repeat_is_refused` — the shipped default
  `filament_half_length_um=0.005` raises, and the message carries both `0.010` and `0.036`.
- `test_the_guard_can_pass_and_therefore_can_fail` — a vacuity control: the same call with sourced
  arguments does **not** raise, so a green suite is not a guard that never runs.

## 10. Numerical and precision envelope

Exact comparisons only — `2 * half_length >= repeat` on float64, no tolerance, because both sides
are declared constants and a tolerance here would be a place for a value to hide. `LiteratureRange`
bands are inclusive.

## 11. Production-backend residency and transfer

**None.** This is a build-time guard on host scalars; it runs once per `build_whole_cell` and
transfers nothing to a device. No kernel, no Warp array, no residency implication.

## 12. Comments and docstrings to discard

Nothing is copied, so nothing is discarded. One provider-facing string is **quoted for refutation**
in §2 (`-3640`'s provenance) and is cited, not inherited.

## 13. Acceptance

`PROPOSED`. It becomes `AUDITED` when §8 and §9 pass and the failure count §4 predicts is measured
and published. It may not become `ACCEPTED` while `V-4` — the cortical actin filament length — is
still a named blank in `SOURCED_CONSTRUCTION_VALUES.md`, because this entry can then reject 10 nm
but cannot say what to use instead.

## 14. Honest limits

- **This refuses; it does not repair.** After it lands, `build_whole_cell` at its shipped defaults
  raises. That is the pre-registered step 2 of
  `PROPOSAL-a-bond-has-no-formation-criterion-and-no-rupture-criterion.md` §5, and the failures it
  produces are the map for steps 3 and 4 — not a regression to be silenced.
- **It does not set the filament length.** It supplies a floor (one crossover repeat) that rejects
  the shipped 10 nm. Choosing a value needs V-4 and a source this session did not read.
- **It does not change `rest_length_um`.** `cortex_surface_coupling.py:1407,1433` still take the
  measured separation. That is `-3651`, and it is deliberately separate so that this entry's
  refusal is measured before another change moves the same numbers.
- **`-3640`'s 60 nm capture radius is reported, not corrected.** Two molecules that differ by 4.4×
  cannot share one reach, but repairing that is a change to `-3640`'s transcription and belongs
  with its pin.
