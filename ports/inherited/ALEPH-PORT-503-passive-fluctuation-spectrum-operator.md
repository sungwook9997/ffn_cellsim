# ALEPH-PORT-503 — passive thermal fluctuation spectrum operator

| Field | Value |
|---|---|
| Lane | L5 (observation layer) |
| Target | `aleph/observe/fluctuation.py` |
| Reference consulted | **NONE.** No `ffn_cellsim` module was read for this operator. |
| Source | `ALEPH-DQ-102` (first observation modality), manuscript §8.1 for the manifest coupling |
| Port class | **ORIGINAL** |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. What the operator is

The first concrete `ObservationOperator`. Given a time series of membrane vertex displacements, it
projects onto a mode basis, forms the mode power spectrum with a per-mode sample count, and reports
the trusted wavevector range. It is the raw-observable side of `ALEPH-DQ-102`; the `(kappa, sigma)`
fit that turns a spectrum into parameters belongs to `aleph/infer/` (L8) and to the analytic oracle
in `validation/analytic/helfrich.py` (L4), which this module must not and does not import.

## 1a. Why not clean-room

It **is** clean-room. Nothing was ported and no `ffn_cellsim` module was read for this operator, so
the question "why does source-derived porting beat clean-room here" has the default answer: it does
not. The projection, the spectrum, and both trusted-range bounds are derived from scratch in §2, and
the source-identity fields of the template are therefore not applicable rather than unrecorded.

## 2. Derivation

### 2.1 Projection

Write the normal displacement field on `V` vertices as `h_v(t)` [µm]. A mode basis supplies
`Phi` of shape `(V, M)` with columns `phi_n`, each labelled by a wavevector magnitude `q_n`
[1/µm], plus per-vertex areas `A_v` [µm²]. The basis is orthonormal under the **area-weighted**
inner product, which is the discretisation of the continuum `L²` inner product on the surface:

    <f, g> = sum_v A_v f_v g_v  /  A_total

so that `<phi_m, phi_n> = delta_mn`. Mode amplitudes are then

    a_n(t) = sum_v A_v phi_n(v) h_v(t) / A_total          [µm]

The area weighting is not cosmetic. On a non-uniform mesh an unweighted dot product weights dense
regions more heavily, which biases the spectrum toward whatever the mesher happened to refine — a
mesh artefact that would be read as physics.

Units: `a_n` is in µm, matching the L2 lane's agreed convention (`q` in 1/µm, areas in µm², mode
amplitudes in µm).

### 2.2 Power spectrum

    P(q_n) = < |a_n(t)|^2 >_t                              [µm²]

The mean is subtracted per mode first. Not subtracting it measures the *static* shape of the
configuration rather than its fluctuation, and for a mode with a non-zero equilibrium offset that
offset dominates the variance. `ddof = 1` is used, because the mean was estimated from the same
samples.

Per-mode sample counts are reported two ways:

- `n_samples[n]` — finite samples actually used.
- `n_effective[n]` — the independent-sample count from `aleph.observe.stationarity`, since
  consecutive frames of a relaxing membrane are strongly correlated and `n_samples` would overstate
  the precision of `P(q_n)` by `sqrt(2 * tau_samples)`. This is the ALEPH-PORT-501 correction applied
  where it actually bites.

The per-mode standard error uses `n_effective`, never `n_samples`.

### 2.3 The trusted wavevector range, stated honestly

A discrete mesh cannot report a spectrum at arbitrary `q`. Two bounds, both measured from the mesh:

**Upper — the mesh Nyquist limit.** With mean vertex spacing `a` [µm], the shortest resolvable
wavelength is `2a` (two samples per period), so

    q_max = pi / a

Above this, a mode is aliased: the sampled field cannot distinguish it from a longer-wavelength one,
and any power reported there is the power of some other mode wearing its label. This is the bound
the brief singles out, and it is the one most easily violated by accident, because a mode basis is
usually generated from a formula that will happily produce as many modes as it is asked for.

**Lower — the finite box.** The longest wavelength that fits in a patch of total area `A_total` is
`L = sqrt(A_total)`, so

    q_min = 2*pi / L

Below this there is no full period inside the domain; the "mode" is a gradient across the patch,
its variance is dominated by the boundary condition rather than by the fluctuation spectrum, and one
sample of it per frame is one sample of a quantity with no defined mean inside the window.

Modes outside `[q_min, q_max]` are **excluded from the reported spectrum and listed by index in
`excluded_modes`, with the reason**. They are not silently dropped and they are not reported with a
warning attached to a number a reader will use anyway.

`trusted_wavevector_range` travels in the result. A result that has a spectrum but no range is a
spectrum someone will extrapolate.

### 2.4 Coordination with `aleph.state.spectral`

L2 owns the mode basis on the real manifold. As of writing, `aleph/state/spectral.py` does not exist.
Rather than block or guess at its API, this module defines a **`ModeBasis` Protocol** — the concept
and the units, not the implementation — with the members `basis_matrix`, `wavevectors_inv_um`,
`vertex_areas_um2`, `total_area_um2` and `mean_vertex_spacing_um`. Any L2 object with those members
satisfies it structurally; nothing needs to be imported in either direction, so the two lanes cannot
break each other.

`PlaneWaveModeBasis` ships here as the concrete reference implementation for a flat periodic patch.
That is not a stand-in: a flat periodic patch **is** the Helfrich geometry of `ALEPH-DQ-102`, the
geometry the closed-form oracle is written for. When L2's curved-manifold basis lands it should
satisfy the same Protocol and reproduce this one in the flat limit — a cross-lane test worth writing
once both exist. **Stub declared:** see §5.

## 3. Why `raw_observable` and `apparent_quantity` are separate here

`raw_observable()` returns the mode power spectrum — the thing an experiment measures, after its own
projection and averaging, with its own sample counts and trusted range. `apparent_quantity()` reports
a single scalar from it, and the default is the **equipartition check**: the mean mode energy in
units of `kB*T/2`.

That choice is deliberate and it is the negative-control-friendly one. `ALEPH-DQ-102`'s argument for
this modality is that on the exact taut-string case equipartition is an *identity*, not a limit — an
estimator that cannot recover `<E_n> = kB*T/2` there is simply broken, with no fitting and no
threshold involved. So the apparent quantity is chosen to be the thing that can be wrong
unambiguously. The `(kappa, sigma)` fit is a different, harder question and lives elsewhere.

## 4. Controls shipped

All under `tests/observe/test_fluctuation.py`.

| Control | Test | Asserts |
|---|---|---|
| Positive | `test_plane_wave_basis_is_orthonormal_under_area_weighting` | `<phi_m, phi_n> = delta_mn` to 1e-12 |
| Positive | `test_plane_wave_basis_satisfies_the_mode_basis_protocol` | the §2.4 seam holds structurally |
| Positive | `test_basis_areas_sum_to_the_patch_area` | quadrature weights are areas |
| Positive | `test_basis_wavevectors_are_ascending_and_positive` | ordering |
| Positive | `test_injected_mode_variances_are_recovered` | known `P(q)` recovered within sampling error |
| Positive | `test_single_mode_excitation_appears_in_exactly_one_bin` | no leakage between modes |
| Positive | `test_parseval_mean_square_displacement_matches_the_spectrum_sum` | `sum_n P_n = <h^2>_A` to 1e-9 |
| Positive | `test_static_offset_does_not_enter_the_spectrum` | an offset 1000x the signal changes nothing |
| Positive | `test_trusted_range_is_the_derived_formula` | `q_min = 2 pi / sqrt(A)`, `q_max = pi / a` |
| Positive | `test_trusted_range_carries_the_provenance_of_both_bounds` | the two bounds fail differently and say so |
| Positive | `test_nyquist_scales_inversely_with_vertex_spacing` | refinement moves the upper bound and not the lower |
| Positive | `test_equipartition_identity_is_recovered_on_the_exact_case` | apparent quantity `= 1.0 +/- 0.02` in `kB T / 2` units |
| Positive | `test_effective_sample_count_is_the_frame_count_for_independent_frames` | uncorrelated limit |
| Positive | `test_effective_sample_count_is_below_the_raw_count_for_correlated_frames` | ALEPH-PORT-501 applied: `n_eff < 0.2 N`, error bars grow accordingly |
| Positive | `test_apparent_quantity_records_the_full_analysis_chain` | the §3 separation is visible in the output |
| Positive | `test_identifiability_declares_the_amplitude_degeneracy` | the degeneracy is declared, not discovered |
| Positive | `test_the_result_carries_the_manifest_hash` | the join key travels |
| Positive | `test_the_operator_satisfies_the_observation_operator_protocol` | structural conformance |
| **Negative** | `test_modes_beyond_mesh_nyquist_are_excluded_with_a_reason` | the honesty requirement, enforced; every basis mode is either reported or listed as excluded |
| **Negative** | `test_a_basis_entirely_beyond_nyquist_is_refused` | nothing left to report -> typed refusal, not an empty spectrum |
| **Negative** | `test_single_frame_is_refused` | one snapshot has no ensemble |
| **Negative** | `test_a_one_dimensional_array_is_refused` | a single frame cannot be mistaken for a series |
| **Negative** | `test_non_finite_displacement_is_refused` | no silent NaN spectrum |
| **Negative** | `test_shape_mismatch_between_series_and_basis_is_refused` | typed refusal |
| **Negative** | `test_missing_owner_and_missing_array_are_refused` | refusal precedes arithmetic |
| **Negative** | `test_out_of_domain_temperature_is_refused` | applicability enforced |
| **Negative** | `test_apparent_quantity_without_stiffnesses_is_refused_not_guessed` | the operator must not invent the numbers that make the ratio come out at 1 |
| **Negative** | `test_unweighted_projection_biases_a_nonuniform_mesh` | demonstrates the §2.1 failure the area weighting prevents |
| **Negative** | `test_result_without_a_trusted_range_cannot_be_constructed` | structural |
| **Negative** | `test_result_rejects_ragged_arrays` | structural |
| **Negative** | `test_basis_rejects_impossible_requests` | four impossible geometries, including over-requesting modes |
| **Negative** | `test_range_helpers_reject_impossible_geometry` | zero, negative, NaN, infinite spacing and area |

### A note on what the mesh-Nyquist control taught us

The first `PlaneWaveModeBasis` built integer wavevectors up to `n//2`, the discrete transform's
Nyquist index. On a **cell-centred** grid that index is degenerate: `cos(pi*(j + 1/2))` samples to
exactly zero at every vertex, so the mode is not representable and the basis construction failed.
The honest default is therefore `n//2 - 1`, one short of the Nyquist index, and `include_beyond_
nyquist=True` builds up to `n - 1` — producing columns that are genuine **aliases** of lower-`q`
ones. Such a basis is deliberately *not* orthonormal, and that is the aliasing being demonstrated
rather than a construction bug. It is documented on the class.

## 4a. Units, domains, singular cases, invariants

**Units**, matching L2's agreed convention throughout: `q` in 1/µm, areas in µm², mode amplitudes
and vertex displacements in µm, and therefore mode power `P(q)` and its standard error in µm².
`tau_int_s` is in seconds, inherited from `ObservationContext.sample_interval_s`. `n_samples` is an
integer count and `n_effective` a float count, both dimensionless. The apparent quantity is
dimensionless by construction — a mode energy in units of `kB*T/2`.

**Domain.** A `(frames, vertices)` array of finite normal displacements with at least 4 frames, a
vertex count matching the basis, a temperature inside the declared applicability band, and a
positive sampling interval. On the wavevector axis the domain is the **measured** band
`[2*pi/sqrt(A_total), pi/a]` — not a caller preference, and not negotiable.

**Singular cases.**

- Every basis mode outside the trusted band: `Refusal(BEYOND_RESOLUTION_LIMIT)`. Returning an empty
  spectrum would look like a measurement that found nothing.
- A single frame, or a 1-D array: `Refusal(INSUFFICIENT_SAMPLES)` / `Refusal(INCONSISTENT_INPUT)`.
  One snapshot has a shape, not a spectrum.
- A mode that samples to zero on the grid (the cell-centred Nyquist index, where
  `cos(pi*(j+1/2)) = 0` at every vertex): skipped during basis construction as not representable,
  with the requested-count check then refusing to pad with duplicates.
- Missing stiffnesses or `kB*T`: `Refusal(UNSUPPORTED_QUANTITY)` rather than a guess. The operator
  must not be able to choose the numbers that make the equipartition ratio come out at 1.
- Any non-finite displacement: refused, never masked. Dropping frames would change the sampling
  interval that every correlation time here is expressed in.

**Invariants.**

- The basis is orthonormal under the area-weighted inner product to 1e-12 (whenever
  `include_beyond_nyquist` is off; with it on the basis is deliberately *not* orthonormal, and that
  aliasing is the thing being demonstrated).
- Mode accounting is closed: `len(reported) + len(excluded) == len(basis modes)`. Nothing vanishes.
- Every reported `q` satisfies `q_min <= q <= q_max`, and every exclusion carries a reason naming
  which bound it violated.
- `sum_n P_n` equals the area-weighted mean square displacement to 1e-9 (Parseval's identity on the
  spanned subspace).
- `n_effective <= n_samples` per mode, inherited from ALEPH-PORT-501's floor.
- A `SpectrumResult` cannot be constructed without a `TrustedRange` or without a manifest hash, and
  its arrays must all have the same width. All three are enforced in `__post_init__`.
- Adding a constant to every displacement changes no reported power (1e-9).

## 4b. Numerical and precision envelope, and production-backend residency

**Numerical.** float64 throughout. The projection is a single `(T, V) @ (V, M)` matrix product with
area weights folded in; the spectrum is `numpy.var(..., ddof=1)` per mode. `ddof=1` because the mean
was estimated from the same samples. The per-mode relative error is `sqrt(2/(m-1))` for a variance
from `m` independent Gaussian samples, with `m = n_effective` and never the frame count.

Tolerances in the controls are statistical rather than numerical and are stated per test: 8% on
variance recovery at 4,000–6,000 frames, 2% on the equipartition identity at 40,000 frames. The
exact-arithmetic assertions — orthonormality, Parseval, offset invariance — are held to 1e-9 to
1e-12, which is round-off at these magnitudes. Basis-mode degeneracy is detected at a norm below
1e-12.

**Residency.** Host-side, CPU only, dense `numpy`. The transfer boundary is explicit and is the
reason this is affordable: the operator consumes a *recorded* `(frames, vertices)` displacement
series that the runtime has already reduced and copied to host, not live device state. The
projection is `O(T V M)` and the correlation-time estimate `O(M T log T)`, both trivial at the
frame counts a passive measurement produces. A device implementation would buy nothing and would
couple the observer to the backend, which the `ModeBasis` and `AcceptedState` Protocols exist to
prevent.

## 5. Stubs and cross-lane dependencies — declared

| Dependency | State at writing | What this module did |
|---|---|---|
| `aleph.state.spectral` (L2) | **DOES NOT EXIST** | Defined the `ModeBasis` Protocol locally with the agreed units; ships `PlaneWaveModeBasis` as the concrete flat-patch implementation. No import in either direction. |
| `aleph.units` (L1) | Partially landed — `aleph/units/__init__.py` imports `aleph.units.guard`, which was not on disk | **Not imported at all.** Units are carried as documented strings on the result objects and asserted in tests. When L1 settles, the `units` string fields are the obvious attachment point for a `Dimension`. |
| `validation.analytic.helfrich` (L4) | Not consulted | Forbidden by the firewall and by design. This module produces a spectrum; judging it against Helfrich is L4's job and must stay on the other side of the wall. |

## 6. Compliance with PLAN §0.2

1. No source text carried — no source read. ✔
2. Projection, spectrum, and both range bounds derived here (§2). ✔
3. New prose. ✔
4. Positive and failing negative controls (§4). ✔
5. Entry written before the code. ✔
6. No `ffn_cellsim` import; no `validation` import. ✔
