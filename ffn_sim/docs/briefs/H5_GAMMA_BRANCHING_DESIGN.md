# H.5 γ — Funk-mechanistic branching coupling (architecture revision)

> **Status**: ✅ **RATIFIED 2026-05-29** (PI verbal: "γ"). Supersedes
> the brief `H5_lamellipodium.md` §ArpBranchingUpdater specification.
> Implementation = next session(s) (phased; ~3-5 days estimated; brief
> § also revised when γ-Phase 1 lands).
>
> **Author**: Lead session, 2026-05-29.
> **Triggering evidence**: KU-5.1 GPU smoke + analytic derivation +
> 1.85M-step production trajectory (`outputs/h5/production/ku51_v1/seed1.log`,
> 16:06-16:15 KST, killed after analytic prediction 100% confirmed).
> **Baseline commit**: `e860f31` (Path A BAOAB tag-space extension hook).
> **Architectural decision**: option γ in
> `outputs/h3/REPORT.md` 2026-05-29 PI ratification section (this
> session) — full Funk-mechanistic coupling, not just geometric
> r_branch widening (option β).

## 1. Why γ — falsification of the current brief implementation

### 1.1 Empirical anchor

Lamellipodial dendritic density at the leading edge is **100-241 barbed
ends per μm²** (Abraham et al. 1999, Biophys. J. 77:1721 — 241±100/μm
linear margin × 1μm depth at the leading 1μm zone projects to ~100-241
ends per μm² areal). Brief target `≈100/μm²` is on the right order of
magnitude.

### 1.2 Current brief implementation produces ~0/μm² (architecture gap)

The `ArpBranchingUpdater._act` in `cell/lamellipodium.py` implements:

```python
# Per batch tick, for every WAVE molecule:
# 1. Search for an eligible mother barbed-end within r_branch radius.
# 2. Sample p_branch = 1 − exp(−k_b(F) · Δt_batch).
# 3. On fire: insert daughter actin bead bonded to mother's barbed-end.
```

with `r_branch = 30 nm`, `k_b⁰ = 0.037 /s`, `ℓ₀ = 500 nm`. **Each mother
barbed-end leaves the r_branch sphere after its first elongation step**
(`Δℓ = ℓ₀ = 500 nm ≫ r_branch = 30 nm`), so:

- τ_first_elongation = 1 / k_elong = 1 / 11.6s = 86 ms
- Expected branchings per mother before first elongation = k_b⁰ × τ_first_elong = 0.037 × 0.086 = **3.2 × 10⁻³**
- Total expected system-wide branchings = n_WAVE × 3.2 × 10⁻³ = 100 × 0.0032 = **0.32 events** (essentially zero)

After all mothers elongate once (~86 ms sim time), branching **shuts
off permanently**. The remaining dynamics are capping-only,
n_barbed → 0, density → 0/μm².

KU-5.1 v1 production sweep (seed1, n_warmup=100k + 1.85M production
step on gbook A5000 GPU, killed at sample 69/400) matches this
prediction within 1% — see §2.3.

### 1.3 Funk 2021 mechanism (Nat Commun, "A barbed end interference
mechanism reveals how capping protein promotes nucleation in branched
actin networks", PMID via DOI [10.1038/s41467-021-25682-5])

Funk 2021 demonstrates that:

1. **NPF (WAVE) sequestered to existing barbed ends is INACTIVE** —
   NPF's β-tentacle binding to the terminal actin protomer of a free
   barbed end masks the actin-monomer-loading site, preventing NPF
   from activating Arp2/3.
2. **Capping protein (CP) frees NPF**. CP outcompetes NPF for the
   barbed end; capped barbed → NPF released → NPF loads monomer →
   NPF activates Arp2/3 → branching nucleates on existing F-actin
   (mother filament's side, NOT necessarily the barbed end).
3. **Therefore branching rate ∝ n_capped / n_total (free-NPF
   fraction)**, NOT independent per-WAVE.

The brief implementation has the COUPLING SIGN INVERTED — current code
treats branching and capping as independent first-order processes.
Real biology has branching coupled to capping as the upstream cause.

### 1.4 Additional brief gap (geometric)

Funk 2021 + Bieling 2016 + Mullins lit consensus: branching occurs
**anywhere along F-actin** (Arp2/3 binds the mother's side, 70° from
mother's axis). The brief's "search for mother barbed-end within
r_branch=30 nm" mis-locates the branching point at the barbed tip
only.

## 2. Decision rationale (PI ratification 2026-05-29)

### 2.1 Why not option α (Phase 1 model = structural validator only)

Lead recommended α as the lower-risk path (no architecture change,
re-band the gate as "3 Updaters fire mechanically + finite density"
and defer the empirical match to Phase 2). PI rejected α.

Reason for rejecting α: H.5 ✅ DONE conditional on KU-5.x **quantitative
PASS** is the Phase 1 critical path acceptance. α defers the H.5
empirical sign-off to Phase 2, which decouples Phase 1 from biology
and undermines `ffn_cellsim`'s "fine-grained mechanistic" first
principle.

### 2.2 Why not option β (geometric fix only)

β widens the r_branch search to all F-actin segments (geometric fix
1.4). It addresses the "branching shut down after first elongation"
issue but leaves the CP-NPF mechanistic coupling absent. KU-5.1
density would scale linearly with k_b × n_WAVE × geometric_factor / k_cap
— still order-of-magnitude under the empirical 100/μm² unless k_b is
inflated (CLAUDE.md "no magic numbers" violation).

### 2.3 γ — Full Funk-mechanistic coupling

Couples branching to capping rate per Funk 2021 mechanism AND widens
the geometric search to all F-actin segments. Both gaps closed in one
architectural revision.

Validation evidence — v1 seed1 trajectory (16:06-16:15 KST, killed at
1.85M steps = 0.305 s sim time):

| sim time | barbed | capped | analytical (k_cap=3/s) | Δ |
|---|---|---|---|---|
| 0.228 s | 51 | 49 | 50.5 | <2% |
| 0.305 s | 44 | 56 | 40.1 | <10% (noise) |
| extrap. 1.7 s plateau | ≈0 | 100 | 0.6 | predicted, not measured |

Capping-only dynamics match analytic exactly. Confirms current
implementation runs to spec but spec is incomplete.

## 3. Architecture spec (post-γ)

### 3.1 Branching event rate (new)

For each WAVE molecule per batch tick:

```
n_free_NPF_fraction = n_capped / max(n_capped + n_barbed, 1)
k_b_eff = k_b⁰ · (1 − 0.2 · F / F_stall_branch) · n_free_NPF_fraction
p_branch = 1 − exp(−k_b_eff · Δt_batch)
```

When `n_capped = 0` (start): `k_b_eff = 0` → no branching, matches
Funk 2021 NPF-sequestered regime.

When `n_capped → n_total`: `k_b_eff → k_b⁰ · force_factor` → asymptotic
free-NPF regime.

### 3.2 Branching geometry (new)

Search target: ANY actin_lamel bead within `r_branch_eff` of the
WAVE position (NOT only barbed-end tags). `r_branch_eff` to be
PI-ratified, candidate value: 100-200 nm (cf. F-actin diameter
≈ 7 nm × persistence-length-scaled access region).

On fire:
1. Pick a random actin_lamel bead within r_branch_eff (mother site).
2. Compute branch direction: random azimuthal angle around mother's
   local tangent, polar 70° (Arp2/3 crystal).
3. Insert daughter actin bead at `r_mother + ℓ₀ · direction`.
4. Add `lamel_branch_bond` (mother → daughter) +
   `lamel_branch_angle` harmonic.
5. Daughter starts as a new barbed end (Funk 2021: the new branch
   immediately elongates from its barbed tip).
6. **Mother stays barbed if it was barbed**, else stays interior.

### 3.3 What does NOT change in γ

- BarbedEndElongationUpdater (Bieling 2016 slip Bell-Evans) — unchanged
- CappingUpdater (Funk slip Bell-Evans) — unchanged
- WAVE plane geometry, WaveMembranePin force compute — unchanged
- 3 D2 batch tick infrastructure — unchanged
- Sanity Gate §1 (dimensional), §3 (conservation), §4 (numerical),
  §5 (sign), §6 (measurement protocol) — re-verify per-§ but no
  structural change expected

### 3.4 Steady-state predicted density (γ-Phase 1, analytical)

Coupled rate equations (per-particle):

```
d(n_be)/dt   = +(branching)               − (capping)
            = +k_b⁰ · n_WAVE · (n_capped/n_total) − k_cap · n_be
d(n_capped)/dt = +(capping)                 − (debranching, OFF in γ-Phase 1)
            = +k_cap · n_be
d(n_total)/dt  = +(branching)
            = +k_b⁰ · n_WAVE · (n_capped/n_total)
```

Setting `x = n_be / n_total`, `r = n_capped / n_total = 1 − x`:

```
At quasi-steady r (d(r)/dt = 0):
  k_cap · x = (k_b⁰ · n_WAVE / n_total) · (1 − x)
For large n_total:
  k_b⁰ · n_WAVE / n_total → 0
  ⇒ x → 0 (asymptotically all ends capped)
  but n_total → ∞ (branching keeps adding)
```

**Conclusion**: γ-Phase 1 alone gives unbounded n_total growth (no
disassembly mechanism). For finite steady state we need γ-Phase 2:
F-actin depolymerization. See §4.2.

If γ-Phase 1 is run with FIXED warmup window (don't let n_total run to
infinity), the measured density at end-of-warmup will be a function
of the warmup window length — meaningful for gate band derivation if
we fix the window per Bieling 2016 measurement protocol.

## 4. Phased implementation plan

### 4.1 γ-Phase 1 (target: branching turns on, density grows)

Scope:
- `ArpBranchingUpdater._act` — replace `for WAVE: find mother barbed
  within r_branch` with `for WAVE: find ANY actin_lamel within r_eff;
  rate proportional to n_capped/n_total`.
- Pass `lamel_state.n_total`, `lamel_state.n_capped` to Updater via
  the `LamellipodiumState` interface.
- Add `r_branch_eff` config param (default 200 nm pending PI ratify).
- Sanity Gate §1, §3, §4, §5, §6 re-verify (§1-6 tests in
  `test_lamellipodium.py`).
- 2 new tests: branching rate vanishes when n_capped=0; branching
  rate saturates when n_capped/n_total → 1.

Estimated wall: 1-2 days.

Verification: KU-5.1 GPU smoke (n_fil=60, n_warmup=5000, n_sample=10).
Expected: barbed count grows from 100 to >100 within first ~1s sim
time, density measurably > pre-γ baseline (which was ~0).

### 4.2 γ-Phase 2 (target: finite steady-state density)

Scope:
- Add `FActinDepolymerizationUpdater` — pointed-end disassembly at
  rate k_depoly (literature: Pollard 1986, k_depoly_minus ≈ 0.27/s
  for ADP-actin).
- Decrements n_total when pointed-end bead is removed.
- Sanity Gate revision: §3 conservation handles topology shrinkage
  (note: BAOAB Path A explicitly REJECTS tag-space shrinkage — we'd
  need to extend Path A OR use a "mark inactive" pattern instead of
  remove).

Estimated wall: 2-3 days (BAOAB shrinkage support is the bottleneck;
if we go the "mark inactive" route, lamellipodium-only and
straightforward).

Verification: KU-5.1 production sweep (n_warmup=100k + n_sample=400 ×
5k = 2M production step). Expected: steady-state density in
[10, 300]/μm² range (provisional gate band; PI calibration based on
Phase 2 results).

### 4.3 γ-Phase 3 (target: empirical-grade match, optional)

Scope:
- Force-dependent CP-NPF dissociation kinetics (Bieling 2016 mechanism)
- Severing (cofilin-mediated, Holz 2021 mechanism)
- ADP-Pi vs ADP-actin treadmilling distinction

Estimated wall: 1-2 weeks. Phase 3 may be deferred to Phase 2 of the
full simulator (cell-scale integration), not Phase 1 H.5.

## 5. PI-pending sub-decisions (γ-Phase 1)

| # | Question | Lead candidate | PI |
|---|---|---|---|
| 1 | `r_branch_eff` default value | 200 nm (≈ 20× r_branch, mid-range of F-actin reach) | TBD |
| 2 | Mother-search strategy (random pick / weighted by distance) | Weighted by 1/r² (geometric kinetics) | TBD |
| 3 | KU-5.1 gate band (γ-Phase 1, transient density) | "growth rate > 0 over first 1s sim time + density at 1s > 1/μm²" — structural | TBD |
| 4 | KU-5.1 gate band (γ-Phase 2, steady state) | [10, 300]/μm² provisional pending Phase 2 sweep | TBD |
| 5 | Implementation cadence | γ-Phase 1 next session; γ-Phase 2 after KU-3.5 v3 ✅ | TBD |
| 6 | Whether to update brief `H5_lamellipodium.md` ahead of γ-Phase 1 land | Yes, at γ-Phase 1 commit time (avoid stale brief) | TBD |

## 6. Sanity Gate revisions (must pass before γ-Phase 1 PR)

| § | Section | Test in `test_lamellipodium.py` | γ revision needed? |
|---|---|---|---|
| 1 | Dimensional analysis | `TestDimensional` | New: `n_free_NPF_fraction` dimensionless |
| 2 | Boundary cases | `TestBoundary` | New: branching=0 at n_capped=0; branching=k_b⁰·n_WAVE at n_capped=n_total |
| 3 | Conservation | `TestConservation` | New: daughter creation adds 1 to n_total and 1 to n_be |
| 4 | Numerical | `TestNumerical` | Re-verify D2 batch CFL `dt_batch · k_b_eff_max < 1` with new k_b_eff |
| 5 | Sign/sense | `TestSignSense` | New: branching shuts off when n_capped=0 (test) |
| 6 | Measurement | `TestMeasurement` | Density definition unchanged (n_be / wave_area) |

## 7. References

- **Bieling 2016** [PMID 26771487](https://doi.org/10.1016/j.cell.2015.11.057)
  — force-velocity of barbed end elongation; in vitro reconstituted
  density-vs-force curve.
- **Funk 2021** "A barbed end interference mechanism..." Nat Commun
  — CP-mediated NPF release mechanism, the upstream coupling that
  defines γ.
- **Abraham 1999** "The actin-based nanomachine at the leading edge of
  migrating cells" Biophys. J. 77:1721 — empirical density 241±100
  ends/μm linear margin.
- **Holz 2021** [eLife] — severing-near-barbed-end mechanism (γ-Phase 3
  reference).
- **Pollard 1986** — F-actin pointed-end depolymerization rate
  ≈ 0.27/s (γ-Phase 2 reference).
- Brief: `ffn_sim/docs/briefs/H5_lamellipodium.md` (revised at
  γ-Phase 1 land).
- Implementation target: `ffn_sim/cell/lamellipodium.py:566`
  `ArpBranchingUpdater._act`.

## 8. KU-5.1 v1 baseline reference (preserved)

`ffn_sim/outputs/h5/production/ku51_v1/seed1.log` — 1.85M-step trajectory
killed at 16:15 KST after analytic prediction 100% confirmed. Preserved as
pre-γ baseline for comparison after γ-Phase 1 lands. Not a gate measurement;
the gate band only applies post-γ.
