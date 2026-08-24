# ALEPH-PORT-203 — spherical mode basis and trusted wavevector band

| Field | Value |
|---|---|
| Lane | L2 (geometry and topology primitives) |
| Target | `aleph/state/spectral.py` |
| Reference consulted | **None.** No `ffn_cellsim` file was read for this module. |
| Port class | **CLEAN ROOM** — nothing ported; entry filed so the record has no gap |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. Why an entry exists for un-ported code

PLAN §0.2 requires a ledger entry before ported code lands. This module is not ported: there is no
counterpart in the reference tree, and none was consulted. The entry is filed anyway so that a later
reader can tell "clean room" from "nobody wrote it down", which are the same absence in the ledger
and very different facts.

## 2. Decisions taken, with the alternative that was rejected

### 2.1 Real spherical harmonics, not a flat-patch 2-D FFT

The task allowed either. Spherical was chosen:

1. Real spherical harmonics are the exact eigenbasis of the Laplace-Beltrami operator on the sphere
   the mesh actually is. A 2-D FFT is the eigenbasis of a periodic plane, so the patch route needs
   either a boundary artifact or a taper, and both contaminate the long-wavelength end where the
   tension term `sigma q^2` is read off.
2. The icosphere has no regular grid. An FFT route would have to resample onto a Cartesian grid
   first, and that interpolation is an uncontrolled low-pass filter sitting between the geometry and
   the observation — biasing the short-wavelength end, which is where `kappa q^4` is read off. The
   harmonic route evaluates the basis at the vertices themselves, so the resampling error is zero
   rather than small.
3. The quasi-spherical (Milner-Safran) form of the Helfrich spectrum is already a statement about
   `<|u_lm|^2>`.

The cost is that the geometry assumption becomes explicit, so the module *enforces* it:
`assert_quasi_spherical` refuses a mesh whose sphericity departs from 1 by more than 5%. A strongly
deformed membrane gets a typed refusal, not a number.

### 2.2 Degree to wavevector

`-Delta_S Y_lm = l(l+1)/R^2 Y_lm` against the plane-wave eigenvalue `q^2` gives
**`q_l = sqrt(l (l+1)) / R`** in 1/um, with `R` the area-equivalent radius `sqrt(A / 4 pi)` so that a
mildly deformed membrane degrades gracefully rather than needing a hand-supplied nominal radius.

### 2.3 The Nyquist argument, stated in full

The field is sampled at the vertices, whose mean spacing is the mean edge length `d`. Two samples
per wavelength is the information floor, so `lambda_min = 2d` and `q_nyquist = pi / d`. Above that
nothing is resolved; content there aliases *down into the band below*, which is why exceeding the
limit corrupts the coefficients that are resolved rather than merely failing to add new ones.

Two samples per wavelength assumes a perfect grid, a bandlimited signal, and an anti-alias filter.
None of the three holds: the icosphere spacing is irregular (min/max edge ratio 0.837 at level 5),
there is no filter, and a discrete operator on the mesh loses accuracy well before the floor. The
default `safety_factor = 0.5` therefore sets the *trusted* ceiling at `q_max = 0.5 * pi / d`, i.e. at
least four vertices per wavelength. `trusted_wavevector_range` returns the hard Nyquist value and the
trusted value as separate fields so the choice is visible rather than baked in.

A second, independent ceiling: a real expansion to degree `L` has `(L+1)^2` coefficients and the mesh
supplies `V` samples, so `L <= sqrt(V) - 1`. The right comparison is against the **hard Nyquist**
degree, since both are then information-theoretic limits rather than one of them being a policy
choice. On an icosphere at level `n`, hard Nyquist gives `l ~ 2.61 * 2^n` and the count gives
`~3.16 * 2^n`. **Two arguments with nothing in common agree in scaling and differ by about 20%**
— measured ratio 1.18 at level 3 (20 vs 24) and 1.19 at level 4 (41 vs 49). That is a real check on
both. The *trusted* ceiling then sits a further factor `safety_factor = 0.5` below the Nyquist one
(9 at level 3, 20 at level 4), so in practice the sampling limit binds; the intersection with the
count limit is taken regardless, because on a mesh coarser or more irregular than an icosphere it
need not.

`q_min` is not a mesh property at all. On a closed vesicle degree 0 is a volume change and degree 1 a
rigid translation; neither bends the membrane, and the `(l-1)` factor of the quasi-spherical Helfrich
spectrum kills them. The first deforming mode is `l = 2`, so `q_min = sqrt(6) / R`, and refinement
never lowers it — only a larger vesicle does.

### 2.4 Least squares over quadrature

Quadrature (`a_lm = sum_v f_v Y_lm(v) A_v / R^2`) is one matrix-vector product and is exact for
perfectly uniform sampling. The icosphere's vertex areas are not uniform, and the measured leakage
into other coefficients when decomposing a pure `Y_lm` is **1.2e-2 to 1.4e-2** at level 3. The
area-weighted least-squares solve returns the same field with leakage **~2e-15**, i.e. exact for any
field inside the span regardless of sampling irregularity. Least squares is the default; quadrature
stays available and is documented with its cost.

### 2.5 Normalised Legendre recurrence, not `factorial * lpmv`

The textbook form needs `sqrt((l-m)!/(l+m)!) * P_l^m`. For moderate `l` the factorial ratio underflows
to zero while `lpmv` overflows, so the product silently becomes `0 * inf`. A three-term recurrence
carrying the normalisation inside every step is used instead, and its orthonormality is verified
numerically against a Gauss-Legendre times uniform-phi quadrature: **max |G - I| = 8.0e-14** at
`L = 8`.

## 3. Firewall

`aleph/state/spectral.py` imports `numpy` and `aleph.state.manifold`. It does **not** import
anything under `validation/`, and must not: the runtime side may never depend on the oracle side.
Compatibility with `validation/analytic/helfrich.py` is by convention only, and the convention is
written into the module docstring as a unit table — `q` in 1/um, areas in um^2, mode amplitudes in um,
`power_per_mode` in um^2 as the mean over the `2l+1` degenerate orders (not the sum).

## 4. Controls shipped

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/state/test_spectral.py::test_basis_is_orthonormal_on_the_sphere` | `max abs(G - I) < 1e-10` under quadrature |
| Positive | `...::test_basis_reproduces_the_first_harmonics_in_closed_form` | `Y_00`, `Y_10`, `Y_11`, `Y_20` match their analytic expressions |
| Positive | `...::test_least_squares_recovers_a_pure_harmonic_exactly` | coefficient 1.0, leakage < 1e-12 |
| Positive | `...::test_quadrature_recovers_a_pure_harmonic_only_approximately` | recovers it, with leakage ~1e-2 — the documented cost |
| Positive | `...::test_degree_and_wavevector_are_mutual_inverses` | round trip to 1e-12 |
| Positive | `...::test_trusted_band_matches_the_nyquist_derivation` | `q_nyquist == pi / mean_edge_length` |
| Positive | `...::test_sampling_and_dof_limits_agree_within_about_twenty_percent` | both `~2^n`, agreeing within 25% |
| Positive | `...::test_in_band_spectrum_agrees_across_resolutions` | same field, three meshes, in-band error falls with refinement |
| **Negative** | `...::test_power_above_the_trusted_band_is_wrong` | out-of-band error is >10x the in-band error (measured 31x) |
| **Negative** | `...::test_non_spherical_surface_is_refused` | `NotQuasiSphericalError` on an ellipsoid |
| **Negative** | `...::test_underdetermined_degree_is_refused` | more coefficients than vertices raises |
| **Negative** | `...::test_bad_safety_factor_is_refused` | outside `(0, 1]` raises |
| **Negative** | `...::test_mesh_too_coarse_for_degree_two_is_refused` | no trusted band raises rather than returning an empty one |

## 5. Compliance with PLAN §0.2

Items 1-3 are vacuous (nothing was read, so nothing could be carried). Items 4-6 hold: controls
shipped, entry written before the code, no `ffn_cellsim` import.
