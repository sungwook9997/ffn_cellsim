# ALEPH-PORT-400 — analytic dimension algebra and physical constants

| Field | Value |
|---|---|
| Lane | L4 (analytic oracles) |
| Target | `validation/analytic/dimensions.py`, `validation/analytic/constants.py` |
| Reference consulted | `/Users/sw1/ffn_cellsim/ffn_sim/virtual_cell/sandbox_fluctuation.py:241-331` |
| Port class | **RE-DERIVED** (no source text carried) |
| Written | 2026-07-30 |
| Status | LANDED |

## 1. What the reference does

Defines a frozen `Dimension` 4-tuple over (length, mass, time, temperature) with `*`, `/`, `**`,
then `helfrich_variance_dimension(bending_power, tension_power)` which composes
`kB*T / (A * (kappa*q**b + sigma*q**t))` symbolically and raises when the two denominator terms are
incommensurate or when the result is not an area.

## 2. Independent re-derivation

Base assignments, taken from the definitions of the quantities, not from the reference:

- `kappa` is the coefficient of `(1/2)(nabla^2 h)^2` integrated over area. `[kappa] * L^-4 * L^2 = E`,
  so `[kappa] = E = M L^2 T^-2`.
- `sigma` is the coefficient of `(1/2)|grad h|^2` integrated over area. `[sigma] * L^-2 * L^2 = E`
  wait — `[sigma] * (L/L)^2 * L^2 = E` gives `[sigma] = E L^-2`, i.e. force per length. Both readings
  agree: `E L^-2 = (M L^2 T^-2) L^-2 = M T^-2 = (M L T^-2)/L = F/L`.
- `[q] = L^-1`, `[A] = L^2`, `[kB T] = E`.

Commensurability of the denominator:
`[kappa q^b] = E L^-b`, `[sigma q^t] = E L^-2 L^-t = E L^-(2+t)`.
Equal iff **`b = t + 2`**.

Result must be an area:
`[kB T / (A kappa q^b)] = E / (L^2 · E L^-b) = L^(b-2)`, which equals `L^2` iff **`b = 4`**.

The two constraints together pin `(b, t) = (4, 2)` uniquely. `b = 3` fails the first test
(`E L^-3` vs `E L^-4`); `(b, t) = (5, 3)` passes the first and fails the second (`L^3`, not an area).

**Verdict on the reference: the code is correct.** Its docstring claim that any pair other than
`(4,2)` "up to a common shift" fails commensurability is *imprecise* — a common shift `(4+s, 2+s)`
stays commensurate and is caught by the second (area) test, not the first. Both tests are needed and
the reference does run both. Aleph's version documents which test catches which mistake.

## 3. Constants

`kB = 1.380649e-23 J/K` is **exact by definition** of the kelvin in the 2019 SI redefinition
(CODATA/BIPM), not a measured value with an uncertainty. Recorded as exact with that provenance.
`T = 310.15 K` (37 degC) is a *modelling choice*, not a constant, and is labelled as such.

## 4. Controls shipped
> All controls live under `tests/validation/`, one file per target module, named
> `test_<module>.py`. The rows below name the suite in words and then the test function,
> rather than giving a path: the ledger index scanner reads any bare `validation/...py`
> substring as a claim about a source file at the repository root, and reads any
> `test_<something>` token as a test function that must exist. A literal file path would
> trip both. The function names are exact and are what the discipline test resolves.


| Control | Location | Asserts |
|---|---|---|
| Positive | (dimensions suite) `test_helfrich_powers_4_2_compose_to_an_area` | `(4,2)` yields `L^2` |
| Negative (must fail) | (dimensions suite) `test_wrong_q_power_raises` | `(3,2)`, `(4,3)`, `(5,3)`, `(2,0)` all raise |
| Negative (must fail) | (dimensions suite) `test_only_four_two_survives_a_brute_force_sweep` | a brute-force sweep of 169 power pairs leaves exactly `(4,2)` |
| Negative (must fail) | (helfrich suite) `test_wrong_q_power_raises_at_construction` | `HelfrichSpectrum(bending_power=3)` raises |
| Negative (must fail) | (dimensions suite) `test_the_two_tests_are_not_redundant` | `(5,3)` passes commensurability and is caught only by the area test |

## 5. Compliance with PLAN §0.2

1. All original comments/docstrings stripped — nothing carried. ✔
2. Law re-derived independently (§2 above). ✔
3. New prose written for Aleph. ✔
4. Positive + failing negative control. ✔
5. This entry written before the code. ✔
6. Passes with `ffn_cellsim` absent from `sys.path` (no import of it anywhere). ✔

## 6. Deviation from the reference

- Aleph parameterises `HelfrichSpectrum` by `bending_power`/`tension_power` so the dimensional guard
  is reachable from a constructor call, making "someone writes q^3" a literally testable event
  rather than a hypothetical.
- Aleph attempts `aleph.units.dimension` first (L1) and falls back to the local algebra. See the
  `TODO(L1)` in `validation/analytic/dimensions.py`.
