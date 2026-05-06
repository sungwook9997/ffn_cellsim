# V2 Adhesion-Network Layer: Cadherin/Integrin Adhesion-Program Dynamics — Locked Design

**Date**: 2026-05-05 KST
**Authors**: Claude + Codex design-discussion (3-round adversarial lock,
PI id=809 + id=1985 aggressive debate posture; PI id=2003 framing
invitation + id=2025/2035 (α) confirmation explicitly absorbed)
**Source unit**: design-discussion `topic=v2-layer-3-adhesion-network-framework`
+ `topic=v2-layer-3-adhesion-network` (consolidated), MCP id 1989–2035
**PI directives**:
- `id=1985`: 병행 진행 + 건강하고 활발한 토론 + framework 우선 + MCF7/epithelial
  generalization (260313 fitting NOT)
- `id=2003`: int-β1 ↔ E-cad 완전 대립 가정 NOT (framing invitation)
- `id=2025`: (α) 충분
- `id=2035`: (α) 채택하고 발전시키는 방향
**Tier**: **A-tier full debate** (NOT B-tier compressed). PI `id=1985`
explicitly requested "건강하고 활발한 토론"; Codex 13 catches across
2 review rounds + Y17 acceptance trace.
**Naming disambiguation (Y1+Y16)**: This unit is the **v1 Layer 3
(adhesion network φ-ODE) reformulation in v2 cell-resolved
architecture**. The v2 roadmap's `Layer 3` (per
`docs/v2/00_project_vision_v2.md`) is **ECM fiber network** — a separate
unit. PI's "Layer 3" reference in `id=1980` etc. means the
adhesion-network unit per CLAUDE.md original Layer 3 definition.
This lock uses "adhesion-network layer" wording to avoid the v2
roadmap conflict.
**PI ratify status**: full delegation per PI id=939/1008/1985/2025/2035.
impl-work uses this for the adhesion-network Sanity Gate doc + code
entry.

---

## 0. Scope + PI directive integration (Y16 central anchor)

### What this unit IS

- v2 reformulation of v1 Layer 3 adhesion-network state contract +
  first ODE
- Cell-level **two-program independent state** representing cadherin
  (cell-cell adhesion) vs integrin (cell-substrate adhesion) program
  preferences
- Reciprocal-coupling ODE with **exact affine analytical update** (no
  CFL, no explicit Euler tolerance games)
- Pure read-only diagnostic: state evolves over time, but does NOT
  feed back into FA rates / traction / motility (separate future
  A-tier units)
- Generalization target: **MCF7 and similar epithelial cell spheroids**
  (PI `id=1985`); literature provenance only, NO PI 260313 data fitting

### Wording boundary (Y4 + Y16, central anchor — PI `id=2003` directive)

> "This framework absorbs PI `id=2003` directive: **integrin-β1 and
> E-cadherin are NOT a strict-opposition or conserved-binary-switch
> relationship**. The two programs are **partly competing axes, not
> exclusive**. Coexistence (both high, mixed transitional state, both
> low) MUST remain representable. Reciprocal cross terms in the rate
> law are **signal-dependent coupling**, NOT enforced antagonism."

This wording is locked verbatim in §0 + module docstring + function
docstring + return-class docstring per Y11 sister discipline.

### What this unit is NOT

- ❌ NOT strict opposition / single φ canonical state (Y3)
- ❌ NOT `cadherin_program + integrin_program = 1` conserved (Y16
  forbidden — PI directive)
- ❌ NOT "E-cad ↔ Intβ1 transition" wording (Y16 — implies strict
  binary)
- ❌ NOT FA `bound_fraction` / `maturity` field reuse as program state
  (Y2 — distinct from FA mechanical state; circular definition risk)
- ❌ NOT magic-blend signal weights (Y7)
- ❌ NOT explicit Euler ODE (Y9 — exact analytical only)
- ❌ NOT default rate values (Y6 — caller must supply)
- ❌ NOT silent epsilon in derived φ undefined case (Y11 — explicit
  `None`)
- ❌ NOT silent zero-coercion for missing junction `cadherin_proxy`
  (Y13 — fail-closed)
- ❌ NOT ULA/pV4D4 phenotype as core framework axis (Y4 — cross-check
  only, future unit)
- ❌ NOT FA rate / traction / motility feedback (Y5 — first unit
  diagnostic-only)
- ❌ NOT v2 roadmap Layer 3 (ECM fiber network) — separate unit (Y1)
- ❌ NOT PI 260313 data fitting (Y8 + Hard Rule 1)
- ❌ NOT EMT / spatial polarity / signaling pathway / leader-follower
  framing (PI `id=2025` confirmed (α) sufficient; future expansion
  paths)

---

## 1. Final Lock Summary

### Constants + naming

```python
# No literature-anchored constants in this lock body.
# Rates are caller-supplied per Y6 (no defaults to avoid magic-number
# leakage). Test fixtures may use simple values like 1e-4.
# Production/science configs must carry external provenance (Cho 2020
# + general epithelial spheroid adhesion remodeling literature).
```

### Typed dataclasses

```python
import numpy as np
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True, slots=True)
class CellAdhesionNetworkState:
    """Cell-level adhesion-program state (v2 first lock).

    Upstream cell-program signal favoring cadherin (cell-cell) vs
    integrin (cell-substrate) adhesion mode. **Distinct from FA
    mechanical state (FocalAdhesionState.bound_fraction/maturity) and
    junction mechanical state (JunctionState.cadherin_proxy)** — those
    are downstream effectors / input signals, NOT the
    adhesion-network program itself.

    PI `id=2003` directive: NOT strict opposition; NOT conserved sum.
    Coexistence (both high, both low, mixed) is representable.
    """
    cell_id: str
    cadherin_program: float           # ∈ [0, 1]; cell-program preference for cadherin mode
    integrin_program: float           # ∈ [0, 1]; cell-program preference for integrin mode

    def __post_init__(self) -> None:
        if not self.cell_id:
            raise ValueError("cell_id must be non-empty")
        for name, v in [("cadherin_program", self.cadherin_program),
                        ("integrin_program", self.integrin_program)]:
            if isinstance(v, bool):
                raise ValueError(f"{name} must be float (not bool — Python bool ⊂ int trap), got {type(v).__name__}")
            if not np.isfinite(v):
                raise ValueError(f"{name} must be finite, got {v}")
            if not (0.0 <= v <= 1.0):
                raise ValueError(f"{name} must be in [0, 1], got {v}")

    @property
    def phi_integrin_dominance(self) -> Optional[float]:
        """Derived diagnostic. None when both programs are 0 (Y11).

        NOT a stored field; computed property only. Branch-defined,
        no epsilon hiding. Returns None for both-zero (undefined
        dominance), explicit ratio otherwise.
        """
        total = self.integrin_program + self.cadherin_program
        if total == 0.0:
            return None
        return self.integrin_program / total


@dataclass(frozen=True, slots=True)
class AdhesionNetworkSignals:
    """Realized engagement signals input to the ODE (Y10 + Y14).

    Built from FA/Junction objects via separate signal builders
    (compute_fa_engagement_signal, compute_junction_engagement_signal).
    Caller cannot pass raw counts ambiguously; signals are pre-normalized
    to [0, 1].
    """
    fa_engagement_signal: float       # ∈ [0, 1]
    junction_engagement_signal: float  # ∈ [0, 1]

    def __post_init__(self) -> None:
        for name, v in [("fa_engagement_signal", self.fa_engagement_signal),
                        ("junction_engagement_signal", self.junction_engagement_signal)]:
            if isinstance(v, bool):
                raise ValueError(f"{name} must be float (not bool), got {type(v).__name__}")
            if not np.isfinite(v):
                raise ValueError(f"{name} must be finite, got {v}")
            if not (0.0 <= v <= 1.0):
                raise ValueError(f"{name} must be in [0, 1], got {v}")


@dataclass(frozen=True, slots=True)
class AdhesionNetworkRates:
    """All rates required, no defaults (Y6 sister with HB#4-active Y3
    k_active discipline). Production must carry external provenance."""
    k_int_on: float                   # integrin assembly rate (1/s)
    k_int_off: float                  # integrin disassembly rate (1/s)
    k_cad_on: float                   # cadherin assembly rate (1/s)
    k_cad_off: float                  # cadherin disassembly rate (1/s)

    def __post_init__(self) -> None:
        for name, v in [("k_int_on", self.k_int_on), ("k_int_off", self.k_int_off),
                        ("k_cad_on", self.k_cad_on), ("k_cad_off", self.k_cad_off)]:
            if isinstance(v, bool):
                raise ValueError(f"{name} must be float (not bool — Python bool ⊂ int trap), got {type(v).__name__}")
            if not np.isfinite(v):
                raise ValueError(f"{name} must be finite, got {v}")
            if not v > 0.0:
                raise ValueError(f"{name} must be strictly positive, got {v}")


@dataclass(frozen=True, slots=True)
class AdhesionNetworkStepDiagnostics:
    """Typed diagnostics for one adhesion-network step (Y12)."""
    integrin_steady_state: Optional[float]      # a_int/(a_int+b_int), None if both zero
    cadherin_steady_state: Optional[float]      # a_cad/(a_cad+b_cad), None if both zero
    integrin_relaxation_rate_per_s: float       # a_int + b_int
    cadherin_relaxation_rate_per_s: float       # a_cad + b_cad
    fa_engagement_signal: float                  # echoed for reproducibility
    junction_engagement_signal: float            # echoed for reproducibility
    dt_s: float                                  # echoed for run traceability (Codex id=2034 FYI)


@dataclass(frozen=True, slots=True)
class AdhesionNetworkStepResult:
    """Single-step result (Y12)."""
    state: CellAdhesionNetworkState
    diagnostics: AdhesionNetworkStepDiagnostics
```

### Step function

```python
def step_adhesion_network_state(
    state: CellAdhesionNetworkState,
    signals: AdhesionNetworkSignals,
    rates: AdhesionNetworkRates,
    dt_s: float,
) -> AdhesionNetworkStepResult:
    """One adhesion-network state step using exact affine update.

    ODE form (PI `id=2003` directive: signal-dependent coupling, NOT
    strict opposition):

        dI/dt = k_int_on · S_fa       · (1 - I) - k_int_off · S_junction · I
        dC/dt = k_cad_on · S_junction · (1 - C) - k_cad_off · S_fa       · C

    Exact analytical solution (Y9/Y17, sister with HB#1+#2 expm1 pattern):

        a = k_on  · signal_on
        b = k_off · signal_off
        if a + b == 0: x_new = x   (exact no-op for zero signals)
        else:
            x_ss = a / (a + b)
            x_new = x_ss + (x - x_ss) * exp(-(a + b) · dt_s)

    First unit is **state-evolving diagnostic-only** (Y5):
    - State evolves over time
    - No mutation of FA rates, traction, junction forces, or contour motility
    - Future feedback into FA rates is a separate A-tier unit
    """
    _validate_dt_s(dt_s)

    # Integrin axis (signal_on = S_fa, signal_off = S_junction)
    integrin_new, i_ss, i_rate = _exact_affine_step_with_diagnostics(
        x=state.integrin_program,
        signal_on=signals.fa_engagement_signal,
        signal_off=signals.junction_engagement_signal,
        k_on=rates.k_int_on,
        k_off=rates.k_int_off,
        dt_s=dt_s,
    )

    # Cadherin axis (signal_on = S_junction, signal_off = S_fa)
    cadherin_new, c_ss, c_rate = _exact_affine_step_with_diagnostics(
        x=state.cadherin_program,
        signal_on=signals.junction_engagement_signal,
        signal_off=signals.fa_engagement_signal,
        k_on=rates.k_cad_on,
        k_off=rates.k_cad_off,
        dt_s=dt_s,
    )

    new_state = CellAdhesionNetworkState(
        cell_id=state.cell_id,
        cadherin_program=cadherin_new,
        integrin_program=integrin_new,
    )

    diagnostics = AdhesionNetworkStepDiagnostics(
        integrin_steady_state=i_ss,
        cadherin_steady_state=c_ss,
        integrin_relaxation_rate_per_s=i_rate,
        cadherin_relaxation_rate_per_s=c_rate,
        fa_engagement_signal=signals.fa_engagement_signal,
        junction_engagement_signal=signals.junction_engagement_signal,
        dt_s=dt_s,
    )

    return AdhesionNetworkStepResult(state=new_state, diagnostics=diagnostics)


def _exact_affine_step_with_diagnostics(
    x: float, signal_on: float, signal_off: float,
    k_on: float, k_off: float, dt_s: float,
) -> tuple[float, Optional[float], float]:
    """Y9/Y17 exact analytical solution to dx/dt = k_on·signal_on·(1-x) - k_off·signal_off·x.

    Returns (x_new, x_steady_state_or_None, relaxation_rate).
    """
    a = k_on * signal_on
    b = k_off * signal_off
    rate = a + b
    if rate == 0.0:
        return x, None, 0.0          # exact no-op (Y9 boundary)
    x_ss = a / rate
    x_new = x_ss + (x - x_ss) * float(np.exp(-rate * dt_s))
    return x_new, x_ss, rate


def _validate_dt_s(dt_s: float) -> None:
    """Y14 dt_s validation. dt_s=0 valid no-op (analytical solution makes 0 trivial)."""
    if isinstance(dt_s, bool):
        raise ValueError(f"dt_s must be float (not bool — Python bool ⊂ int trap), got {type(dt_s).__name__}")
    if not np.isfinite(dt_s):
        raise ValueError(f"dt_s must be finite, got {dt_s}")
    if dt_s < 0.0:
        raise ValueError(f"dt_s must be non-negative, got {dt_s}")
    # dt_s == 0 OK → analytical solution gives x_new = x (exact no-op)
```

### Signal builders (Y7 + Y13 fail-closed)

```python
def compute_fa_engagement_signal(adhesions: tuple[FocalAdhesionState, ...]) -> float:
    """Y7 normalized FA signal: mean of bound_fraction · maturity over all FAs.

    Empty FA list → 0.0 (no integrin engagement signal). Branch-defined
    normalization (Y7), no magic blend weights.
    """
    if not adhesions:
        return 0.0
    products = np.array(
        [fa.bound_fraction * fa.maturity for fa in adhesions],
        dtype=np.float64,
    )
    return float(products.mean())


def compute_junction_engagement_signal(junctions: tuple[JunctionState, ...]) -> float:
    """Y7 normalized junction signal: mean of cadherin_proxy over all junctions.

    Y13 fail-closed: every junction must have non-None cadherin_proxy.
    Silent coercion to 0 is forbidden; caller must filter or populate.

    Empty junction list (single-cell context) → 0.0 (no cadherin
    engagement signal).
    """
    if not junctions:
        return 0.0
    proxies = []
    for j in junctions:
        if j.cadherin_proxy is None:
            raise ValueError(
                f"Junction with cadherin_proxy=None cannot contribute to engagement signal. "
                f"Caller must populate cadherin_proxy or filter junctions before passing."
            )
        proxies.append(float(j.cadherin_proxy))
    return float(np.array(proxies, dtype=np.float64).mean())
```

---

## 2. Reasoned-acceptance trace (Y1–Y17)

The lock converged after 3 rounds (id 1989–2035) under PI directive
id=1985 aggressive debate posture, with PI `id=2003` framing
invitation explicitly absorbed and (α) confirmed by PI `id=2025/2035`.

**Y1 (round 1, accepted Codex C1 = naming disambiguation)**: Lock title
uses "adhesion-network layer" wording to avoid v2 roadmap Layer 3
(ECM fiber network) conflict. v1 Layer 3 reformulation.

**Y2 (round 1, accepted Codex C2 = state ownership distinct from FA)**:
`AdhesionNetworkState` cell-program state distinct from
`FocalAdhesionState.bound_fraction/maturity` (FA mechanical state) and
`JunctionState.cadherin_proxy` (junction mechanical state). The latter
two are downstream effectors / input signals, NOT canonical Layer 3
state; conflation creates circular definition.

**Y3 (round 1, accepted Codex C3 = two-axis bounded ODE, NOT single
canonical φ)**: Both `cadherin_program` and `integrin_program` are
[0,1] bounded, can independently be high/low/mixed. Single φ canonical
state would lose information; derived as diagnostic only.

**Y4 (round 1, accepted Codex C4 — CRITICAL biology direction
correction)**: Original round 1 had ULA/pV4D4 direction REVERSED.
Correct: ULA → cadherin-dominant (no substrate), pV4D4 → integrin-FA
rises (substrate engagement). Verified per `docs/03_adhesion_dynamics.md`.
Without Codex catch, framework would have locked with reversed
biology. ULA/pV4D4 is cross-check only, not core framework axis.

**Y5 (round 1, accepted Codex C5 = state-evolving diagnostic-only)**:
First unit evolves state but does NOT feed back into FA / motility
/ junction dynamics. Future feedback into FA rates is a separate
A-tier unit.

**Y6 (round 1, accepted Codex C6 = no default rates)**: All rate
constants in `AdhesionNetworkRates` required, no defaults. Sister with
HB#4-active Y3 (k_active no-default discipline). Test fixtures may use
simple values; production must carry external provenance.

**Y7 (round 1, accepted Codex C7 = branch-defined normalization, no
magic blend weights)**: `compute_fa_engagement_signal` and
`compute_junction_engagement_signal` use mean of declared per-object
quantity; NO weighted blend like `S = w_a · S_fa + w_b · S_junction`
(magic-number trap).

**Y8 (round 1, accepted Codex C8 = literature first, PI overlay
later)**: Cho 2020 + general epithelial spheroid literature for rate
constants provenance. PI 260313 data NOT used for fitting (Hard Rule
1). Overlay is future visualization unit.

**Y9 (round 2, accepted Codex C12 = exact affine update)**: ODE form
`dx/dt = k_on·signal_on·(1-x) - k_off·signal_off·x` solved analytically
as `x_new = x_ss + (x - x_ss)·exp(-(a+b)·dt_s)`. No CFL, no explicit
Euler tolerance games. Sister with HB#1+#2 expm1 pattern.

**Y10 (round 2, accepted Codex C10 = typed `AdhesionNetworkSignals`,
NOT raw FA/Junction objects)**: Core ODE consumes typed signals
contract; signal builders are separate helpers. Avoids tight coupling
of ODE to optional FA/Junction schema details.

**Y11 (round 2, accepted Codex C9 = derived φ `Optional[float]` with
`None` for both-zero)**: `phi_integrin_dominance` returns `None` when
`integrin_program + cadherin_program == 0` (undefined dominance).
Returning 0.0 would falsely imply E-cad dominance in downstream plots.
Branch-defined, no epsilon hiding.

**Y12 (round 2, accepted Codex shape proposal = result + diagnostics)**:
`AdhesionNetworkStepResult` carries new state + typed
`AdhesionNetworkStepDiagnostics` (steady states, relaxation rates,
echoed signals + dt_s for traceability — Codex `id=2034` FYI).

**Y13 (round 3, accepted Codex C11 = junction signal fail-closed on
`cadherin_proxy is None`)**: `compute_junction_engagement_signal`
raises `ValueError` if any junction has `cadherin_proxy=None`. Caller
must populate or filter; silent zero-coercion forbidden.

**Y14 (round 3, accepted Codex C13 = `dt_s` + `AdhesionNetworkSignals`
validation)**: `dt_s` must be float (not bool — Python bool ⊂ int
trap), finite, non-negative. `dt_s = 0` valid no-op (analytical
solution gives `x_new = x` exactly). Signals dataclass validates
finite + [0,1] + non-bool.

**Y15 (round 3, accepted Codex round 2 #4 + PI `id=2003` = coexistence
test)**: Test 13 with both signals high + symmetric rates →
both programs converge to ≈0.5 (coexistence), NOT winner-take-all.
Forward guard against future framework drift back to strict
opposition. PI directive `id=2003` integration test.

**Y16 (round 3, accepted Codex round 2 wording requirements + PI
`id=2003`)**: Lock §0 + module/function/return-class docstrings
contain PI directive verbatim ("partly competing, not exclusive
axes"). Forbidden list includes `cadherin_program + integrin_program
= 1` enforcement, "E-cad ↔ Intβ1 transition" wording (implies binary),
winner-take-all enforcement.

**Y17 (round 3, accepted Codex C12 reaffirm = lock body specifies
exact affine)**: Lock pseudocode includes the analytical update form
verbatim. Test 14 verifies arbitrary-large `dt_s` keeps both programs
in [0, 1] (boundedness is analytical, not approximate).

**Codex `id=2034` final fixture corrections (incorporated into test
catalog §4)**:
- Test 11 (FA-only "cadherin↓"): initialize `cadherin_program > 0`
  so decrease is observable
- Test 12 (junction-only "integrin↓"): initialize `integrin_program
  > 0` so decrease is observable

---

## 3. Sanity Gate (impl writes in module docstring or sibling sanity doc)

1. **Units**:
   - `cadherin_program`, `integrin_program`: dimensionless ∈ [0, 1]
   - `fa_engagement_signal`, `junction_engagement_signal`:
     dimensionless ∈ [0, 1]
   - `k_int_on`, `k_int_off`, `k_cad_on`, `k_cad_off`: 1/s
   - `dt_s`: s, ≥ 0
   - `a, b` (rate products): 1/s
   - `(a + b) · dt_s`: dimensionless (argument to exp)
2. **Boundary**:
   - Both signals zero: `a + b = 0` → `x_new = x` exact no-op (Y9)
   - Both programs zero: derived `phi = None` (Y11)
   - `dt_s = 0`: every axis no-op via analytical solution
   - Single-cell context (no junctions): `S_junction = 0`; cadherin
     program decays to 0 if `S_fa > 0` (FA pushes integrin up + cadherin
     down via cross term)
   - Adversarial large `dt_s`: bounded in [0, 1] by analytical exp (Y17)
3. **Conservation**:
   - This unit is state-evolving (programs change), but NOT
     conservation-constrained (cadherin + integrin sum is NOT
     constrained to 1, per PI `id=2003` directive Y16)
   - Per-step boundedness `[0, 1]` is analytical (convex combination of
     `x` and `x_ss`)
   - No mutation of FA / junction / motility (Y5)
4. **Numerical**:
   - `np.float64` enforced throughout
   - Exact analytical update (Y9/Y17) — no CFL stability gate, no
     tolerance tuning
   - `np.exp(-(a+b)·dt_s)` numerically stable for any non-negative
     argument (does NOT need expm1 since we're computing decay factor,
     not difference)
5. **Sign**:
   - `a, b ≥ 0` always (rates positive, signals non-negative)
   - `x_ss = a / (a + b) ∈ [0, 1]` (convex combination ratio)
   - `x_new ∈ [0, 1]` if `x ∈ [0, 1]` (convex combination of two values
     in [0, 1])
6. **Measurement-protocol** (Hard Rule 11, central anchor):
   - **PI `id=2003` directive integration (Y15+Y16)**: ODE represents
     two axes that are **partly competing, not strict opposition**.
     Coexistence (both high) is the diagnostic case test 13 enforces;
     winner-take-all (one near 1, other near 0) under symmetric high
     signals would indicate framework drift back to strict-opposition
     and is forbidden by test 13.
   - **Y1 layer naming**: this unit is NOT v2 roadmap Layer 3 (ECM
     fiber network); naming separation prevents future audit confusion.
   - **Y4 biology direction (ULA/pV4D4)**: cross-check only, not
     framework axis. Direction (ULA→cadherin, pV4D4→integrin) verified
     per `docs/03_adhesion_dynamics.md`.
   - **Y8 literature provenance**: rate constants must cite Cho 2020 +
     general epithelial spheroid literature. PI 260313 data is NOT
     used for fitting (Hard Rule 1).
   - **Y2 state ownership**: `cadherin_program` / `integrin_program`
     are upstream cell-program signals, NOT FA mechanical state
     (`bound_fraction`, `maturity`) or junction mechanical state
     (`JunctionState.cadherin_proxy`). Conflation creates circular
     definition with HB#4 multipliers and FA dynamics.

---

## 4. Test Catalog (16 tests, A-tier full)

### State validation (5)

1. `test_cell_adhesion_network_state_cadherin_in_unit_interval` —
   `cadherin_program ∉ [0, 1]` raises `ValueError`
2. `test_cell_adhesion_network_state_integrin_in_unit_interval` —
   `integrin_program ∉ [0, 1]` raises `ValueError`
3. `test_cell_adhesion_network_state_bool_program_rejected` — Python
   `bool ⊂ int` trap, both programs reject `True`/`False`
4. `test_cell_adhesion_network_state_empty_cell_id_raises`
5. `test_cell_adhesion_network_state_sum_constraint_not_enforced`
   (Y16 PI directive forward guard) — `cadherin=0.7, integrin=0.7`
   (sum=1.4) is valid; `__post_init__` does NOT reject

### Derived φ diagnostic (2)

6. `test_phi_integrin_dominance_both_zero_returns_none` (Y11) —
   `cadherin=0, integrin=0` → `phi_integrin_dominance is None`
7. `test_phi_integrin_dominance_general_returns_ratio` —
   `cadherin=0.3, integrin=0.7` → `phi == 0.7 / 1.0 == 0.7`

### Signals + Rates validation (2)

8. `test_adhesion_network_signals_validation` (Y14) — bool rejected,
   finite, [0, 1] for both fields
9. `test_adhesion_network_rates_validation` (Y6) — all 4 rates
   required, finite positive, bool rejected

### ODE evolution (5)

10. `test_step_adhesion_network_state_both_signals_zero_no_op` (Y9
    boundary) — `signals = (0, 0)` → `state` unchanged exactly;
    diagnostics `integrin_steady_state is None` + `cadherin_steady_state
    is None` + relaxation rates 0
11. `test_step_adhesion_network_state_fa_only_integrin_up_cadherin_down`
    (Codex `id=2034` fixture correction) — initial `cadherin=0.7,
    integrin=0.0`, `S_fa=1.0, S_junction=0.0`, equal rates; assert
    `integrin` increases AND `cadherin` decreases (cross term active)
12. `test_step_adhesion_network_state_junction_only_cadherin_up_integrin_down`
    (Codex `id=2034` fixture correction) — initial `integrin=0.7,
    cadherin=0.0`, `S_fa=0.0, S_junction=1.0`, equal rates; assert
    `cadherin` increases AND `integrin` decreases
13. **`test_step_adhesion_network_state_coexistence_under_symmetric_high_signals`**
    (Y15 — PI directive forward guard) — both signals=1.0, all 4
    rates equal, run 100 steps; assert both programs converge to
    ≈0.5 (NOT winner-take-all). Critical PI `id=2003` integration test.
14. `test_step_adhesion_network_state_large_dt_stays_bounded` (Y17) —
    `dt_s = 1e6`; assert both programs ∈ [0, 1] after step (analytical
    boundedness, not approximate)

### Junction signal fail-closed + dt validation (2)

15. `test_compute_junction_engagement_signal_none_cadherin_proxy_raises`
    (Y13) — junction with `cadherin_proxy=None` → `ValueError` (no
    silent zero coercion)
16. `test_step_adhesion_network_state_dt_validation` (Y14) — bool
    `dt_s` rejected, `dt_s=-1.0` rejected, `dt_s=inf` rejected,
    `dt_s=0.0` valid no-op

### Exports + meta-test (1, bundled)

17. (Bundled into existing tests via import statements) Exports
    `CellAdhesionNetworkState`, `AdhesionNetworkSignals`,
    `AdhesionNetworkRates`, `AdhesionNetworkStepResult`,
    `AdhesionNetworkStepDiagnostics`, `step_adhesion_network_state`,
    `compute_fa_engagement_signal`,
    `compute_junction_engagement_signal` through both
    `acs.v2.adhesion_network_state` (state classes) +
    `acs.v2.dynamics.adhesion_network_dynamics` (functions) +
    `acs.v2.dynamics` (re-exports) + `acs.v2` (mirror)

(16 tests catalogued; test 17 export coverage verified via successful
import in test fixtures.)

---

## 5. Files

- `docs/v2/v2_adhesion_network_layer_locked.md` (this file, source of truth)
- `acs/v2/adhesion_network_state.py` (NEW):
  - `CellAdhesionNetworkState` dataclass (Y2 + Y11 derived φ property)
  - `AdhesionNetworkSignals` dataclass (Y10 + Y14)
  - `AdhesionNetworkRates` dataclass (Y6)
  - `AdhesionNetworkStepDiagnostics` dataclass (Y12 + Codex id=2034 dt_s)
  - `AdhesionNetworkStepResult` dataclass (Y12)
- `acs/v2/dynamics/adhesion_network_dynamics.py` (NEW):
  - `step_adhesion_network_state` function (Y9 exact affine + Y3
    reciprocal)
  - `compute_fa_engagement_signal` helper (Y7)
  - `compute_junction_engagement_signal` helper (Y7 + Y13 fail-closed)
  - Private `_exact_affine_step_with_diagnostics` (Y9)
  - Private `_validate_dt_s` (Y14)
- `acs/v2/dynamics/__init__.py` + `acs/v2/__init__.py` — new exports
- `tests/v2/test_v2_adhesion_network_dynamics.py` (NEW, 16 tests per §4)

---

## 6. References

- PI directives:
  - `id=1985`: 병행 진행 + 건강한 토론 + framework 우선 + MCF7 일반화
  - `id=2003`: int-β1 ↔ E-cad 완전 대립 NOT (framing invitation)
  - `id=2025` + `id=2035`: (α) 채택 confirmation
- Sister locks (v2 architecture upstream):
  - `docs/v2/v2_p1_derivation_locked.md`
  - `docs/v2/v2_phase_e_v2_composition_step_2_locked.md`
  - `docs/v2/v2_hard_blocker_4_active_locked.md` (Y6 sister:
    no-default rate constant discipline)
  - `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`
  - `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md`
    (Y9 sister: analytical exact update via expm1 pattern)
  - `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2/v2_phase_f_minimal_cell_motility_pilot_locked.md` (Y4
    sister: HB#4 multipliers diagnostic-only first-unit pattern)
- v1 Layer 3 reference (reformulation source):
  `docs/03_adhesion_dynamics.md` (Y4 biology direction verification:
  ULA → cadherin-dominant, pV4D4 → integrin-FA dominant)
- Naming disambiguation source (v2 roadmap Layer 3 = ECM fiber network):
  `docs/v2/00_project_vision_v2.md` (Y1)
- Literature provenance (Y8 — for impl, not in this lock body):
  - Cho et al 2020 (E-cad ↔ Int-β1 dynamics in epithelial spheroids)
  - General epithelial cell adhesion remodeling literature (TBD per
    impl, MCF7-relevant)
- Hard Rule 1 (no PI data fitting): CLAUDE.md
- Hard Rule 11 (sim-experiment measurement matching): CLAUDE.md +
  memory `hard_rule_11_wording_boundary_meta_test.md`
- Magic-Number Block: CLAUDE.md (Y6 no-default rate discipline)
- Aggressive debate posture: memory
  `feedback_aggressive_design_debate.md` (PI id=809) + this lock
  (PI id=1985 explicit "건강하고 활발한 토론")

---

## 7. Cross-room dispatch

This file is the design-team input to implementation-work for:

1. impl Claude or Codex implements
   `acs/v2/adhesion_network_state.py` +
   `acs/v2/dynamics/adhesion_network_dynamics.py` + tests per §1
   pseudocode + §4 16-test catalog. **Skip Sanity Gate doc** (B-tier
   compressed-style despite A-tier debate) — A-tier debate already
   done in design-discussion.
2. impl-side single review pass on commit covering:
   - Y3 ODE form correctness (reciprocal cross terms, NOT strict
     opposition)
   - Y9/Y17 exact affine update (analytical, not Euler)
   - Y6 no-default rates (caller must supply, bool rejected)
   - Y10 typed signals contract (NOT raw FA/Junction)
   - Y11 derived φ `None` for both-zero
   - Y13 junction signal fail-closed on `cadherin_proxy=None`
   - Y15 coexistence test 13 — PI directive forward guard
   - Y16 wording: lock title + module docstring contain "partly
     competing, not exclusive" + NOT "E-cad ↔ Intβ1 transition"
   - Y2 state distinct from FA `bound_fraction`/`maturity`
3. On PASS: adhesion-network layer first lock milestone reached.
4. After milestone: design-discussion idle on adhesion-network layer.
   Next forward units (per PI directives + Codex id=2027 expansion
   path):
   - **Future Layer 3 expansion (β/γ/δ/ε)**: separate A-tier full
     debate cycles; only opened if PI/data demand it
   - **Adhesion-network → FA rate feedback** (separate A-tier unit,
     Y5 deferred): currently diagnostic-only; second unit can
     introduce actual feedback into FA assembly/disassembly rates
   - **PI 260313 overlay** (B-tier visualization): cross-check
     trajectories against PI experimental data (overlay only, NOT
     fitting per Hard Rule 1)

A-tier debate cycle MCP id 1989–2035 (rounds 1-3 + Codex SEAL ack
+ PI confirmation).
