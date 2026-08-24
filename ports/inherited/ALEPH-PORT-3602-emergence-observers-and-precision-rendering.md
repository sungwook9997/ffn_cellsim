# ALEPH-PORT-3602 — Q-tensor nematic order, a finite-N isotropic null band, and a condensation detector, as observers

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3602` |
| Lane | `46143f30` Track V (visualisation) |
| Status | `PROPOSED` |
| Written | `2026-07-31` — **before the code**, per `PLAN.md` §0.2.5 |
| Port class | `SOURCE-DERIVED` |
| Exists because | Aleph's renderer can show *what* the step did and cannot say whether the structure it draws is **ordered** or is a picture of noise. At finite `N` an isotropic population does not score zero, so a bare threshold on an order parameter reports order in noise. This entry adds the order parameter, the null band that makes it a measurement, and the renderer that draws the result. |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
from aleph.viz.emergence import (
    ALIGNMENT_COSINE_MIN,
    MIN_NEIGHBOURS,
    NEIGHBOURHOOD_SCALE,
    ORDER_Z_CRITICAL,
    CondensedBundle,
    EmergenceObservation,
    LocalFields,
    NullBand,
    condensed_bundles,
    isotropic_band_asymptote,
    isotropic_eigenvalue_scale_sq,
    isotropic_null_band,
    local_fields,
    nematic_order,
    nematic_order_and_director,
    neighbourhood_radius,
    observe_emergence,
    order_excess_z,
    q_tensor,
    sample_isotropic_directors,
    segment_directors,
    segment_lengths,
    segment_midpoints,
)

from aleph.viz.render import (
    ISOLATION_STYLE,
    DrawnScene,
    IsolationStyle,
    OwnerElements,
    draw_scene,
    render_scene_svg,
)

from aleph.viz.figure import (
    EXIT_OK,
    EXIT_REFUSED,
    ArtefactError,
    figure_from_artefact,
    main,
    read_artefact,
)
```

The one command (`V2`):

```bash
/Users/sw1/miniconda3/envs/aleph/bin/python -m aleph.viz.figure <artefact> -o <figure.svg>
```

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e587`) |
| Source path | `ffn_sim/ac/emergence/nematic.py`, `null_model.py`, `condensation.py`, `detector.py` |
| Source symbol(s) | `q_tensor`, `nematic_order`, `nematic_order_and_director`, `fiber_axes`, `fiber_centroids`; `nematic_null_scale_sq`, `sample_isotropic_axes`, `nematic_null_band`, `effect_size`, `is_ordered`, `EFFECT_SIZE_Z_CRIT`, `NullBand`; `neighborhood_radius`, `compute_local_fields`, `localize_bundles`, `LocalFields`, `Bundle`; `EmergenceDetector`, `EmergenceReport`, `detect` |
| Read from | **working tree** |
| Working tree == commit? | **yes** — `git diff be0e5876 -- ffn_sim/ac/emergence/` is empty, and `git status --porcelain` on that directory is empty. Verified 2026-07-31, both commands run read-only. |

Nothing in that repository was written, moved, or modified. The four files were read once, in full, to
understand the method; no file was copied and no line survives into `aleph/`.

## 3. Why source-derived porting beats clean-room

**One idea crossed the boundary, and it is not the mathematics.** It is the discipline:

> An order parameter measured on a finite population must be scored against **the distribution that
> population would produce if it were disordered**, not against zero — and the detector that uses it
> must be an observer that no force channel can read.

That is a hard-won empirical insight rather than a convention any author reaches by default. The
tempting design — "call it ordered when `S > 0.3`" — is what a first implementation looks like, and it
is wrong in a way that does not announce itself: for `N = 16` segment directors drawn **isotropically**
the expected order parameter is already `E[S] ≈ 0.21`, so that threshold fires on pure noise roughly a
third of the time. A detector built that way reports emergence in a run with no structure in it, and
every downstream artefact inherits the claim. Knowing in advance that the null band is the load-bearing
part — and that it must be cross-checked against a closed form rather than trusted as a Monte-Carlo
number — is worth citing a source for.

**What was NOT taken, and was derived here instead:** every quantitative statement. §4 derives the
Q-tensor normalisation, the exact eigenvalue-scale identity, and — going past the source, which states
plainly that the band "has no simple closed form" and Monte-Carlos it — **closed forms for the null
band's own mean and standard deviation**. Those are new in Aleph and are the reason this port is worth
more than the thing it was derived from; see §7.

An honest statement of the alternative: the mathematics in §4 is standard liquid-crystal and
random-matrix material and clean-room re-derivation was entirely feasible — as the derivation there
demonstrates, since it was done. The citation is for the *discipline*, and it is recorded rather than
quietly absorbed.

## 4. Physical or mathematical law represented

Restated in Aleph's vocabulary and derived here.

### 4.1 The nematic order parameter

A filament population is a set of segments. Segment `k` runs between two node positions and has a unit
**director** `n_k` (µm/µm, dimensionless). A director is defined only up to sign — a filament has no
head — so any honest order measure must be invariant under `n -> -n`. That rules out the polar mean
`|<n>|`, which reports zero for a perfectly aligned antiparallel bundle. The lowest-order sign-blind
object is the outer product `n (x) n`, so the measure is built from its population mean.

Define the second-moment tensor and the traceless **Q-tensor**:

```
M = < n (x) n >            (3x3, symmetric, trace 1)
Q = (3/2) M - (1/2) I      (3x3, symmetric, traceless)
S = lambda_max(Q)          (the scalar order parameter, dimensionless)
```

The normalisation is fixed by two anchors and nothing else:

- **Perfect alignment.** All `n_k = u`: `M = u (x) u`, whose eigenvalues are `(1, 0, 0)`. Then `Q` has
  eigenvalues `(1, -1/2, -1/2)` and `S = 1` **exactly**.
- **Perfect isotropy.** `M = I/3` by symmetry, so `Q = 0` and `S = 0`.

`S` is bounded in `[-1/2, 1]`: `Q` is traceless, so its largest eigenvalue is at least `0` only when the
three coincide, and the minimum of `lambda_max` over traceless symmetric tensors is `-1/2` (attained by
`(-1/2, -1/2, 1)` reordered, i.e. a perfectly *oblate* distribution confined to a plane). The
**director** is the unit eigenvector of `lambda_max`; its sign is arbitrary and is fixed by convention
so a figure is reproducible.

### 4.2 The finite-N null band — why zero is the wrong reference

`S = 0` holds only in the infinite-population limit. For `N` directors drawn independently and
uniformly on the sphere, `M` is a *sample* mean and fluctuates about `I/3`, so `lambda_max(Q) > 0`
almost surely. The size of that bias is the whole question.

**An exact identity, valid at every `N`.** Let `t = n (x) n - I/3`, the single-sample traceless part.
Its Frobenius norm is **deterministic**, not random:

```
|t|_F^2 = tr[(n (x) n)^2] - (2/3) tr[n (x) n] + (1/9) tr(I)
        = 1 - 2/3 + 1/3
        = 2/3
```

using `(n (x) n)^2 = n (x) n` for a unit `n`. The `t_k` are independent and zero-mean, so

```
E |M - I/3|_F^2 = (1/N^2) * N * (2/3) = 2/(3N)
```

and since `Q = (3/2)(M - I/3)`,

```
E [ sum_i lambda_i^2 ] = E |Q|_F^2 = (9/4)(2/(3N)) = 3 / (2N).        (*)
```

`(*)` is **exact for every `N >= 1`**, with no appeal to a limit. It is checkable by hand at `N = 1`,
where a single director gives `Q` eigenvalues `(1, -1/2, -1/2)` and `sum lambda^2 = 1 + 1/4 + 1/4 =
3/2`, matching `3/(2*1)`. This is the cross-check that keeps the Monte-Carlo honest.

**Closed forms for the band itself.** `(*)` bounds the eigenvalue *scale* but not `E[S]`. Going
further: the covariance of `t` is, using `E[n_i n_j n_k n_l] = (d_ij d_kl + d_ik d_jl + d_il d_jk)/15`,

```
Cov(t_ij, t_kl) = (1/15)(d_ik d_jl + d_il d_jk) - (2/45) d_ij d_kl
                = (2/15) [ (d_ik d_jl + d_il d_jk)/2 - d_ij d_kl / 3 ]
```

which is exactly the isotropic traceless-symmetric form. So by the central limit theorem `Q` converges
to the traceless Gaussian orthogonal ensemble with variance parameter

```
v = (9/4)(1/N)(2/15) = 3 / (10 N)
```

whose eigenvalue density is `p(l) ~ prod_{i<j} |l_i - l_j| * exp(-sum l_i^2 / 2v) * delta(sum l_i)`.
Parameterise the traceless plane in polar coordinates, `l_k = R sqrt(2/3) cos(psi - 2 pi k / 3)`, which
makes `sum l_k = 0` automatic and `sum l_k^2 = R^2`, and `psi` the true polar angle. The Vandermonde
factor becomes `R^3 |sin 3 psi| / sqrt(2)` (via the cubic discriminant, using
`prod_k cos(psi - 2 pi k/3) = cos(3 psi)/4`), so with the area element `R dR dpsi`,

```
p(R, psi) ~ R^4 exp(-R^2 / 2v) * |sin 3 psi|
```

— and **`R` and `psi` separate**. `R = sqrt(v) * chi_5`, and on the fundamental domain
`psi in [-pi/3, pi/3]` the largest eigenvalue is `R sqrt(2/3) cos psi`. Both integrals are elementary:

```
E[chi_5] = 8 sqrt(2) / (3 sqrt(pi))
E[cos psi]   = (9/8) / (4/3) = 27/32
E[cos^2 psi] = (29/30) / (4/3) = 29/40
```

giving, after simplification, the two closed forms this entry contributes:

```
mu_inf(N)    = 9 / (2 sqrt(10 pi N))                  ~= 0.802866 / sqrt(N)
sigma_inf(N) = sqrt( (29 pi - 81) / (40 pi N) )       ~= 0.283590 / sqrt(N)
```

Both are **asymptotic** (they invoke the CLT), unlike `(*)` which is exact. §10 states the measured
approach rate. The `1/sqrt(N)` scaling is the substantive claim: **the null bias falls only as the
square root of the population**, which is why a fixed threshold cannot be right at two population sizes.

### 4.3 Condensation

A bundle is a patch that is simultaneously **orientationally ordered** and **spatially dense**; density
alone is a crowd, order alone is a texture. So two per-segment fields are measured over a neighbourhood
whose radius is set by the configuration's own median nearest-neighbour spacing (making it
grid-invariant): the local order `S_local`, and the neighbour count. A segment is flagged when its local
order beats the **local** null band at its own neighbour count — the band is a function of `N`, so a
neighbourhood of 6 and one of 40 are held to different bars, which is the entire point — and when its
own director lies within the alignment cone of the local director, which distinguishes a *member* of a
bundle from a segment merely adjacent to one. Flagged segments are then grouped by spatial connectivity.

**This is not a law.** It is a *measurement convention*: it defines what Aleph will call a bundle. It
generates no force, enters no residual, and no accepted step can read it. §11 states that structurally.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| node positions | µm | 1e-6 m | finite |
| segment director `n` | dimensionless | — | unit sphere |
| `Q` | dimensionless | — | symmetric, traceless |
| `S`, `mu`, `sigma` | dimensionless | — | `S in [-1/2, 1]` |
| `z` (order excess) | dimensionless | — | finite; `sigma > 0` required |
| neighbourhood radius | µm | 1e-6 m | `> 0` |
| density | count | — | integer `>= 1` |

Singular and boundary cases, each with the behaviour Aleph requires:

- **A zero-length segment** (coincident endpoints) has no director. Refused with a typed error naming
  the segment index — *not* silently regularised by an epsilon, because a degenerate segment in a
  filament population is a defect in the population and a renderer that hides it is the failure this
  package exists to prevent.
- **Fewer than 2 segments**: `q_tensor` refuses. A single director gives `S = 1` trivially and is never
  evidence of order.
- **A neighbourhood below `MIN_NEIGHBOURS`**: local order is `nan`, never `0.0`. A missing measurement
  and a measurement of zero are different, and `0.0` would be read as "isotropic here".
- **`sigma = 0`** in a band (degenerate Monte-Carlo): `order_excess_z` refuses rather than dividing.
- **An owner with no readable geometry**: listed in `DrawnScene.undrawn` with a reason. Never skipped.

Invariants that must hold, each with the test that asserts it:

- **I1.** `S = 1` exactly for perfectly aligned directors —
  `tests/viz/test_emergence.py::test_aligned_directors_score_exactly_one`
- **I2.** `S` is invariant under a global rotation (eigenvalues are) and under `n -> -n` on any subset
  (the measure is nematic) — `tests/viz/test_emergence.py::test_order_is_rotation_and_polarity_invariant`
- **I3.** `S` is invariant to segment relabelling —
  `tests/viz/test_emergence.py::test_order_is_permutation_invariant`
- **I4.** `E[sum lambda^2] = 3/(2N)` exactly, at every `N` —
  `tests/viz/test_emergence.py::test_the_exact_eigenvalue_scale_holds_at_every_population`
- **I5.** The closed-form band agrees with the Monte-Carlo band, and the residual falls as `1/sqrt(N)` —
  `tests/viz/test_emergence.py::test_the_closed_form_band_and_the_monte_carlo_band_agree`
- **I6.** Every component of a scene is either drawn or named undrawn — never neither —
  `tests/viz/test_render.py::test_every_component_is_either_drawn_or_named_undrawn`
- **I7.** The five isolation states are pairwise distinguishable in the emitted SVG —
  `tests/viz/test_render.py::test_the_five_isolation_states_are_pairwise_distinguishable`

## 6. Source evidence class and known retractions

**Where I looked:** the source's `STATE.md` §(c) not-quotable blocklist (all 17 active entries, read in
full), its `STATE_NONQUOTABLE.md` link targets by name, the four module docstrings' own self-declared
gate claims, and the filesystem for the control files those docstrings name.

- **No entry on the 17-row blocklist covers the emergence corpus.** The blocklist's subjects are motor
  tensions, gate results, particle counts, crawl speeds, dose responses and cortex-motor magnitudes.
  Nothing nematic, condensation- or null-model-shaped appears on it.
- **This is moot for this entry in any case: no magnitude from that repository is quoted here or in the
  code.** Every number in §4 and §10 was derived or measured in Aleph. The check is recorded because
  the discipline requires it to be recorded, not because it was load-bearing.
- **The named controls exist.** Each module docstring cites a control file; all four resolve to real
  files under `ffn_sim/tests/ac/emergence/`. I record this because my first check looked at the wrong
  path prefix and briefly suggested they did not exist — the correction is the finding, and a
  "the source names tests that do not exist" claim was one path-check away from entering this ledger
  falsely.
- **Reachability: not dead code.** `ac/emergence` is imported by `ac/weave/woven_cell.py`,
  `ac/weave/stress_fiber.py` and a visualisation script, in addition to its own tests.
- **Whether those tests exercise the law or only the plumbing: NOT ASSESSED, deliberately.** I did not
  read them, so that Aleph's controls in §8/§9 are written against §4's derivation rather than against
  someone else's idea of what to assert. Carried as a limit in §14.

## 7. Independent oracle or derivation

Four oracles, none of which is the source:

1. **An exact closed-form identity, checkable by hand.** `E[sum lambda_i^2] = 3/(2N)` from §4.2, derived
   from the deterministic `|t|_F^2 = 2/3`. It holds at every `N`, is verifiable analytically at `N = 1`
   (`3/2`), and the Monte-Carlo sampler must reproduce it at every population size tested.
2. **Closed forms for the band's mean and standard deviation** — `mu_inf`, `sigma_inf` of §4.2 — derived
   independently via the traceless-GOE eigenvalue density. **The source states these have no simple
   closed form and Monte-Carlos them.** They are therefore a genuinely independent oracle for the
   sampler, and the direction of the improvement is the thing worth having: the sampler is now
   falsifiable against algebra instead of being the only estimate available.
3. **Exact algebraic anchors.** `S = 1` for aligned directors and `S = 0` for the exactly-isotropic
   tensor `M = I/3`, both from the normalisation in §4.1, both asserted to machine precision.
4. **Symmetry.** Rotation invariance, polarity invariance and permutation invariance are properties of
   the construction and hold to machine precision or the construction is wrong.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/viz/test_emergence.py::test_aligned_directors_score_exactly_one` | Perfectly aligned directors give `S == 1.0` and the director recovers the alignment axis to `1e-15` |
| Positive | `tests/viz/test_emergence.py::test_the_exact_eigenvalue_scale_holds_at_every_population` | Sampled `E[sum lambda^2]` matches `3/(2N)` within Monte-Carlo error at `N in {1,2,4,16,64,256}`; at `N = 1` it is exact to `1e-15` |
| Positive | `tests/viz/test_emergence.py::test_the_closed_form_band_and_the_monte_carlo_band_agree` | `mu_MC/mu_inf -> 1` and the relative residual falls with a measured order of `~0.5` in `log2 N` |
| Positive | `tests/viz/test_emergence.py::test_a_planted_bundle_is_localised_against_an_isotropic_background` | A planted aligned bundle is recovered: its centroid within one neighbourhood radius, recall `>= 0.8`, and background contamination `<= 0.2` |
| Positive | `tests/viz/test_render.py::test_every_component_is_either_drawn_or_named_undrawn` | Drawn owners and `undrawn` owners **partition** the scene's component set exactly |
| Positive | `tests/viz/test_render.py::test_the_five_isolation_states_are_pairwise_distinguishable` | All five states differ pairwise in the emitted SVG on stroke, dash and marker |
| Positive | `tests/viz/test_figure_cli.py::test_one_command_turns_a_run_artefact_into_a_figure` | A single `main([...])` call on a real artefact exits `EXIT_OK` and writes parseable SVG naming every owner |

## 9. Deliberately failing negative control

Each of these asserts a **wrong** input or a **broken** guard is rejected, and each fails if the
rejection stops working. Three carry their own vacuity guard — a scanner that cannot fire is the
recorded failure mode this project has hit four times.

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/viz/test_emergence.py::test_isotropic_directors_do_not_fire_the_detector` | **Specificity.** Across many independent isotropic draws the detector stays silent; the false-positive rate is bounded well under the `z` threshold's nominal rate |
| Negative (must fail) | `tests/viz/test_emergence.py::test_a_bare_threshold_fires_on_noise_where_the_null_band_stays_silent` | **The load-bearing control.** On the *same* isotropic sample a fixed threshold `S > 0.3` fires while the null-band detector does not — the null band is doing work, not decorating |
| Negative (must fail) | `tests/viz/test_emergence.py::test_a_dense_but_isotropic_patch_is_not_called_a_bundle` | Crowding without alignment is not condensation |
| Negative (must fail) | `tests/viz/test_emergence.py::test_a_degenerate_segment_is_refused_rather_than_regularised` | A zero-length segment raises, naming the index, instead of being epsilon-smoothed into a fake director |
| Negative (must fail) | `tests/viz/test_emergence.py::test_no_runtime_or_vertical_module_imports_the_observer` | **The observer rule.** No module under `aleph/vertical/`, `aleph/runtime/`, `aleph/state/` or `aleph/scenarios/` imports `aleph.viz.emergence` |
| Negative (must fail) | `tests/viz/test_emergence.py::test_the_observation_carries_no_force_field` | No field of the observation is force-shaped; **and** the scanner is proved live by `test_the_force_field_scanner_detects_a_planted_force_field` |
| Negative (must fail) | `tests/viz/test_emergence.py::test_the_observer_does_not_write_to_its_inputs` | Given read-only arrays the observer still runs, so it cannot be writing back into producer state |
| Negative (must fail) | `tests/viz/test_render.py::test_not_scheduled_never_renders_like_isolated` | The one state that is a defect is distinct from the one that is a decision, on every visual channel |
| Negative (must fail) | `tests/viz/test_render.py::test_an_unreadable_owner_is_listed_undrawn_rather_than_skipped` | An owner whose arrays cannot be read appears in `undrawn` with a reason; **and** `test_the_undrawn_partition_detects_a_planted_silent_skip` proves the partition check can fail |
| Negative (must fail) | `tests/viz/test_render.py::test_engaged_zero_is_refused_when_the_force_is_merely_small` | `ENGAGED_ZERO` requires an **exact** zero; a small non-zero force is rejected, because the exactness is the evidence |
| Negative (must fail) | `tests/viz/test_figure_cli.py::test_the_cli_refuses_an_artefact_it_cannot_read_and_exits_two` | A refusal exits `2`, never `0` — a tool that refuses and exits `0` lets a pipeline treat "nothing was drawn" as success |

## 10. Numerical and precision envelope

- **Working and accumulation precision:** float64 throughout. `q_tensor` accumulates the outer-product
  mean in float64; the eigen-decomposition is `numpy.linalg.eigh` on a symmetrised `3x3`.
- **Conditioning.** `Q` is formed as `(3 M - I)/2` from a mean of outer products of unit vectors, so
  every entry is `O(1)` and there is no subtraction of large nearly-equal quantities. `eigh` on a
  symmetric `3x3` is backward-stable; the eigenvalue error is `O(eps)` and never `O(eps/gap)` because
  eigenvalues, not eigenvectors, carry the answer. **The director is the exception**: when the two
  largest eigenvalues are nearly degenerate (a nearly oblate distribution) the director is
  ill-conditioned even though `S` is not. That is a property of the physics — a planar-isotropic
  population genuinely has no unique axis — so the director is reported and never asserted tightly.
- **Tolerances and why each number.**
  - `S = 1` for aligned directors: `1e-15`, i.e. a few `eps`. Anything looser would pass a wrong
    normalisation such as `(3M - I)/3`.
  - `3/(2N)` at `N = 1`: `1e-15`, because it is exact algebra and not a sampled quantity.
  - `3/(2N)` sampled: `4 sigma/sqrt(n_draws)` at the measured spread, not a hand-picked percentage —
    a fixed `1%` would be vacuous at large `n_draws` and flaky at small.
  - Rotation/permutation/polarity invariance: `1e-14`. These are exact identities perturbed only by
    the eigensolver.
- **Range of validity, and what happens outside it.** The **exact** identity `(*)` holds at every
  `N >= 1`. The **closed-form band** is asymptotic; measured relative excess of the Monte-Carlo mean
  over `mu_inf` is `+7.6%` at `N = 4`, `+3.6%` at `N = 16`, `+1.8%` at `N = 64`, `+0.83%` at `N = 256`
  and `+0.35%` at `N = 4096` — i.e. `O(1/sqrt(N))`, measured, not assumed.
  **Consequence, and it is a design decision rather than a caveat:** the *detector* uses the
  **Monte-Carlo** band, never the asymptote. Using the asymptote at small `N` would understate `mu` by
  `7.6%` and inflate every `z`, which biases the detector toward *firing* — the exact failure mode this
  entry exists to prevent. The closed form's job is to falsify the sampler, not to replace it.
- **Refusal, not silent degradation:** every domain violation in §5 raises a typed error naming the
  offending index or field.

## 11. Production-backend residency and transfer

**This code never runs on the production backend, and that is the point.**

`aleph/viz/emergence.py` is a host-side observer over host-side NumPy arrays. It is not called from any
accepted step, it contributes to no residual, it returns no force, torque, stress or energy, and no
module under `aleph/vertical/`, `aleph/runtime/`, `aleph/state/` or `aleph/scenarios/` may import it —
asserted structurally by `test_no_runtime_or_vertical_module_imports_the_observer`, which walks those
trees with `ast` rather than trusting a convention.

When the CUDA path lands, the observer reads whatever the viewer already holds: it takes plain arrays
and never a device handle, so it cannot originate a device-to-host copy. The existing
`aleph/viz/buffers.py` contract (no authoritative reach, no writeable view) is unchanged and unweakened
by this entry, and `tests/viz/test_import_boundaries.py` — which globs `aleph/viz/*.py` — picks up all
three new modules with no edit to the guard.

The renderer and the CLI are likewise host-only and write text.

## 12. Comments and docstrings to discard

Every line of source prose is discarded. Nothing from that repository's vocabulary survives. Discarded
specifically, and what replaces it:

| Discarded | Replaced by |
|---|---|
| Module paths (`ff/architecture_metrics`, `ac/emergence/...`) and the "validated form"/"parity with" framing | Aleph's own derivation in §4; the code cites this ledger id and nothing else |
| Gate names, `Sanity Gate` blocks, `Magic-Number Block` headers, `hard-truth #1`, `§5 falsifiability contract`, `I5 anti-coupling firewall`, `STRESS_FIBER_TARGET`, `SEEDED`, `native run`, `PI` escalation language | Aleph's own control table (§8/§9) and the observer rule stated in Aleph's terms |
| Identifiers carrying provider vocabulary (`fiber_*`, `bundle_count`, `EFFECT_SIZE_Z_CRIT`, `NEIGHBORHOOD_SCALE`) | Aleph's segment/µm vocabulary: `segment_directors`, `ORDER_Z_CRITICAL`, `NEIGHBOURHOOD_SCALE` |
| The claim that the band has no closed form | §4.2 derives one, and §10 states the regime where it is and is not the right thing to use |
| Their normalisation comment and their `AXIS_EPS` epsilon-regularisation of degenerate axes | §4.1's two-anchor derivation; and Aleph **refuses** a degenerate segment instead of regularising it (§5) |

`tests/ports/test_port_discipline.py::test_no_provider_vocabulary_leaks_into_the_package` scans
`aleph/**` for provider tokens and is the mechanical half of this section.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | **Not accepted.** `PROPOSED`. `tests/viz`: **277 passed** on 2026-07-31 at the tolerances in §10; the closed forms of §4.2 are corroborated by Monte-Carlo across `N in [4, 4096]` at the residuals recorded there. **Mutation study: 6 planted defects, 6 killed** — see below. |
| Reviewer | **None. Agent-proposed and unratified** — written by session `46143f30` Track V. No PI review has occurred and this entry carries no `decided_by` field. |
| Rollback | Delete `aleph/viz/{emergence,render,figure}.py`, their three test files, and the appended names in `aleph/viz/__init__.py`. **Nothing breaks**: no module outside `aleph/viz/` imports any of them, which is enforced by §9's import-direction control. That the rollback is inert is the strongest available statement that the observer rule holds. |

### 13.1 Mutation study — the controls were made to fail on purpose

Run 2026-07-31 with `PYTHONDONTWRITEBYTECODE=1` (a same-byte-length mutant with a fast edit/revert
leaves a stale `.pyc` that CPython treats as valid, which produced phantom kills for another agent on
this repository tonight). Each mutant was reverted with `git checkout --` and the tree verified clean.

| # | Planted defect | Killed by | Result |
|---|---|---|---|
| M1 | `order_excess_z` scores against **zero** instead of the band mean — i.e. the null band deleted and the naive reference restored | `test_a_bare_threshold_fires_on_noise_where_the_null_band_stays_silent` | **KILLED**, and it is the **only** test that fails. The other 29 pass. |
| M2 | `Q = (3M - I)/3` instead of `(3M - I)/2` | `test_aligned_directors_score_exactly_one` | **KILLED** |
| M3 | Degenerate segment epsilon-regularised instead of refused | `test_a_degenerate_segment_is_refused_rather_than_regularised` | **KILLED** |
| M4 | `draw_scene` silently skips an owner with no positions | `test_every_component_is_either_drawn_or_named_undrawn` | **KILLED** |
| M5 | `NOT_SCHEDULED` given `ISOLATED`'s exact style | `test_not_scheduled_never_renders_like_isolated` | **KILLED** |
| M6 | CLI refusal returns `EXIT_OK` instead of `EXIT_REFUSED` | `test_the_cli_refuses_an_artefact_it_cannot_read_and_exits_two` | **KILLED** |

**M1 is the one worth reading.** The claim made in §9 was that the bare-threshold control is the only
control that would notice the null band being replaced by a constant. That was a prediction, and it is
what happened: 29 of 30 emergence controls — including every invariance, every closed-form check and
the planted-bundle recovery — pass happily with the null band removed. A test suite without that one
control would have reported a healthy detector that fires on noise.

## 14. Honest limits

What this entry does **not** establish:

1. **A number computed from a fixture is not a measurement of a cell.** Every `S`, `z` and bundle count
   produced by the controls comes from synthetic directors constructed by the test. `UNVERIFIED` for any
   claim about real cytoskeletal organisation. Nothing here has been run against a converged whole-cell
   trajectory, because none exists yet in this repository.
2. **The closed-form band is asymptotic and is not claimed otherwise.** `mu_inf` and `sigma_inf` are CLT
   limits with a measured `O(1/sqrt(N))` approach (§10). No finite-`N` error bound is *proved*; the
   convergence order is **measured** across `N in [4, 4096]`, which is weaker than a proof and is stated
   as such. The exact identity `(*)` is the only band-related result carried as exact.
3. **`ORDER_Z_CRITICAL` is a statistical convention, not a physical constant.** It is a significance
   threshold declared in advance. It is not derived from Aleph physics and must never be adjusted so
   that a run passes.
4. **The neighbourhood scale and alignment cone are conventions.** They define what Aleph *calls* a
   bundle (§4.3). Their stability is checked across a scale sweep, but no claim is made that another
   defensible convention would recover the same bundles.
5. **Whether the source's own controls exercise their law is NOT ASSESSED** — deliberately not read, so
   that §8/§9 are independent (§6).
6. **No GPU.** Nothing here has run on a device; the authorization record is expired and an agent may
   never write one. The CUDA-residency statement in §11 is a design constraint that is *structurally*
   enforced (the observer takes arrays, not handles) and is `UNVERIFIED` as a device measurement.
7. **The renderer is not a viewer process.** It emits SVG text. There is no interactivity, no camera
   model beyond a fixed axis-aligned projection, and no claim that the projection is the best one for
   any particular structure.
8. **Occlusion is not resolved.** Elements are drawn in a documented order with depth-derived opacity;
   this is a legibility aid and not a correct hidden-line computation. A dense population will overdraw,
   and the element counts in the legend — not the picture — are the quantitative record.
