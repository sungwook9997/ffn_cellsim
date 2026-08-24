# ALEPH-PORT-401 — Helfrich membrane fluctuation spectrum

| Field | Value |
|---|---|
| Lane | L4 (analytic oracles) |
| Target | `validation/analytic/helfrich.py` |
| Reference consulted | `/Users/sw1/ffn_cellsim/ffn_sim/virtual_cell/sandbox_fluctuation.py:452-591` |
| Port class | **RE-DERIVED** (no source text carried) |
| Written | 2026-07-30 |
| Status | LANDED |
| Decision it serves | `ALEPH-DQ-101` (V1 vertical), `ALEPH-DQ-102` (first observation), `ALEPH-DQ-103` (latent target) |

## 1. Independent derivation

Monge gauge, nearly-flat patch of projected area `A`, height field `h(r)`, `r` in the plane.
Helfrich free energy truncated at quadratic order, zero spontaneous curvature, zero Gaussian-modulus
term (a topological constant for a fixed-topology patch):

```
H[h] = ∫ d²r [ (κ/2)(∇²h)² + (σ/2)|∇h|² ]
```

**Fourier convention matters and fixes the factor of `A`.** Take

```
h(r) = Σ_q h_q e^{i q·r},      h_q = (1/A) ∫ d²r h(r) e^{-i q·r}
```

so `h_q` carries dimensions of length. Parseval for this convention gives `∫d²r f g* = A Σ_q f_q g_q*`.
Hence `∫d²r (∇²h)² = A Σ_q q⁴ |h_q|²` and `∫d²r |∇h|² = A Σ_q q² |h_q|²`, so

```
H = (A/2) Σ_q (κ q⁴ + σ q²) |h_q|²                                     (all q, both signs)
```

Reality of `h(r)` means `h_{-q} = h_q*`, so the sum over all `q` double-counts. Restricting to a
half-space of `q`:

```
H = A Σ_{q ∈ half} E_q [ (Re h_q)² + (Im h_q)² ],     E_q ≡ κ q⁴ + σ q²
```

Each of `Re h_q`, `Im h_q` is one independent real quadratic degree of freedom with coefficient
`A E_q`. Classical equipartition for `H = a x²` gives `⟨x²⟩ = kB T /(2a)`, so
`⟨(Re h_q)²⟩ = ⟨(Im h_q)²⟩ = kB T /(2 A E_q)` and

```
⟨|h_q|²⟩ = kB T / ( A (κ q⁴ + σ q²) )
```

which is the target formula. **The `A` sits in the denominator only because `h_q` was defined with the
`1/A`.** With the other common convention (`h_q = ∫d²r h e^{-iq·r}`) the same physics reads
`⟨|h_q|²⟩ = A kB T / E_q`. Any oracle-vs-simulation comparison must state the convention; this is
recorded in the module docstring and is the single most likely source of a silent factor-of-`A²` error.

**Crossover, derived not asserted.** The two denominator terms are equal when `κ q⁴ = σ q²`. For
`q > 0` this is `κ q² = σ`, so

```
q* = sqrt(σ/κ)        [dimension: sqrt((E L⁻²)/E) = L⁻¹ ✔]
```

Writing `r(q) ≡ κq⁴/(σq²) = (q/q*)²`, the exact local log-log slope is

```
d ln⟨|h_q|²⟩ / d ln q = -(4κq⁴ + 2σq²)/(κq⁴ + σq²) = -(4 r + 2)/(r + 1)
```

which tends to `-4` as `r → ∞` (bending) and `-2` as `r → 0` (tension), passing through exactly `-3`
at `q = q*`. That last fact is a free, sharp, checkable prediction and is asserted in the tests.

## 2. What was found in the reference

| Item | Verdict |
|---|---|
| `⟨|h_q|²⟩ = kB T/(A(κq⁴+σq²))` with `E = (A/2)Σ_q E_q |h_q|²` | **Correct**, and internally consistent (I reproduced both). |
| `q* = sqrt(σ/κ)`, `inf` when `κ=0`, `0` when `σ=0` | **Correct.** |
| `regime_at` with a `dominance_ratio` and an honest `CROSSOVER` label | **Correct and good practice.** Carried as a design idea, re-implemented. |
| `expected_exponent` returns `None` in the crossover | **Correct.** |
| Refusal of `κ = σ = 0` | **Correct** — no restoring term means no normalisable equilibrium. |
| Refusal of `q = 0` | **Correct** — that is rigid translation of the patch, not a fluctuation. |
| Fourier convention | **NOT STATED** in the reference. This is the reference's most consequential omission for a future simulation comparison; Aleph states it explicitly. |
| Exact local slope `-(4r+2)/(r+1)` | **ABSENT** in the reference, which only exposes the two limits. Aleph adds it, since it is what turns "the exponent is about −3 here" into a quantitative gate. |

No error found in the reference's Helfrich physics.

## 3. Controls shipped
> All controls live under `tests/validation/`, one file per target module, named
> `test_<module>.py`. The rows below name the suite in words and then the test function,
> rather than giving a path: the ledger index scanner reads any bare `validation/...py`
> substring as a claim about a source file at the repository root, and reads any
> `test_<something>` token as a test function that must exist. A literal file path would
> trip both. The function names are exact and are what the discipline test resolves.


| Control | Location | Asserts |
|---|---|---|
| Positive | (helfrich suite) `test_matches_closed_form_against_independent_evaluation` | agreement with an independently written expression to 1e-15 rel |
| Positive | (helfrich suite) `test_limits_are_approached_at_the_derived_rate` | limit forms are approached as `(q*/q)²` / `(q/q*)²` |
| Positive | (helfrich suite) `test_slope_is_exactly_minus_three_at_the_crossover` | exact `-3` at `q*` |
| Positive | (helfrich suite) `test_crossover_is_where_the_two_terms_are_equal` | `q*` re-solved numerically by brentq matches `sqrt(σ/κ)` |
| Positive | (helfrich suite) `test_local_slope_matches_a_numerical_derivative` | the exact slope formula against a central difference |
| Negative (must fail) | (helfrich suite) `test_wrong_q_power_raises_at_construction` | `q³` model refused at construction |
| Negative (must fail) | (helfrich suite) `test_local_slope_is_strictly_monotone_between_the_two_limits` | the slope crosses `-3` exactly once, so no member is a `q⁻³` power law |
| Negative (must fail) | (degeneracy suite) `test_a_minus_three_spectrum_is_rejected_by_chi_square` | χ² rejects a `q⁻³` spectrum |
| Negative (must fail) | (helfrich suite) `test_a_surface_with_no_restoring_term_is_refused` | `κ = σ = 0` refused |
| Negative (must fail) | (helfrich suite) `test_zero_or_negative_wavevector_is_refused` | `q = 0` (rigid translation) refused |

## 4. Compliance with PLAN §0.2

1. Comments/docstrings stripped. ✔ 2. Re-derived (§1). ✔ 3. New prose. ✔
4. Positive + failing negative controls. ✔ 5. Ledger before code. ✔ 6. No `ffn_cellsim` import. ✔
