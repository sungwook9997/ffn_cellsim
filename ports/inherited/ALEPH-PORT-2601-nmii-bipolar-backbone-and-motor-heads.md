# ALEPH-PORT-2601 — NMII bipolar minifilament backbone with individual, head-resolved motor crossbridges

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-2601` |
| Lane | `L26 nmii` |
| Status | `PROPOSED` |
| Written | `2026-07-30` (before the code, per PLAN §0.2.5) |
| Port class | `RE-DERIVED` |

---

## 1. Aleph API

The exact public surface this entry authorises. Nothing outside this list is covered.

```python
from aleph.vertical.nmii import (
    # errors
    NmiiError, NmiiParameterError, NmiiGeometryError, NmiiKineticsError, NmiiLedgerError,
    # parameter carrier and the two laws
    DetachmentLaw, StrokeModel, MotorKinetics, motor_parameter, bell_force_scale_from_temperature,
    # state
    STATE_BLOCKS, ActinBindingSites, RigidActinTrack, NmiiPopulation,
    # coupling
    NmiiMotorConnector, MOTOR_CONNECTOR_NAMES,
    # accounting
    MotorPowerBudget, assert_power_budget_closes,
    # measurement and checks
    IsokineticSlidingAssay, SlidingMeasurement,
    assert_state_blocks_match_contract,
    assert_bipolar_backbone,
    assert_load_dependent_detachment,
)
```

The Aleph target file is `aleph/vertical/nmii.py`. Its controls live in
`tests/vertical/test_nmii.py`.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/motor/{hand,hill_fv_analytic,bell_kinetics_analytic,ensemble_stall_analytic,powerstroke_analytic,force_budget_ledger,minifilament_topology,segment_motor}.py`, `ffn_sim/ac/motor/params_i0b3.yaml`, `ffn_sim/ff/myosin_linear.py` |
| Source symbol(s) | Read at the module level for structure and parameter inventory only. No symbol was lifted; see §3. |
| Read from | working tree — **stated explicitly** |
| Working tree == commit? | `yes` — verified with `git diff be0e5876 -- ffn_sim/ac/motor/ ffn_sim/ff/myosin_linear.py`, which is empty |

The last two rows matter here: at `be0e5876` the reference working tree carries uncommitted
changes elsewhere, so "read from the working tree" and "read from the commit" are not
interchangeable statements about that repository in general. For these paths specifically they
coincide, and that was checked rather than assumed.

## 3. Why source-derived porting beats clean-room

**It does not.** Nothing was ported. This entry records a re-derivation and exists so that the
reading is on the record.

The reference's motor tree was read for *inventory and technique*: which quantities a head-resolved
myosin model needs to carry, which analytic special cases are worth owning as oracles, and how a
force budget for an active element is arranged. That is prior art about *what to check*, not about
*how to compute*. Every law in the Aleph module is derived in §4 below from the two-state
crossbridge master equation, and the derivation is reproduced in the module docstring in Aleph's own
vocabulary and Aleph's own units.

The three candidate reasons for a source-derived port were considered and each fails:

- *An empirical failure mode found by a real run.* None applies. The failure mode this module is
  organised around — a motor ensemble that cannot stall because detachment is load-independent — is
  derivable on paper and is shipped here as a negative control, not inherited as a scar.
- *A case analysis that is easy to get subtly wrong.* The branch structure here is three states
  (detached / attached-mid-stroke / attached-post-stroke) and two event types. It is small enough to
  write from the derivation.
- *A numerically conditioned form.* One genuinely delicate expression exists — the mean-field bound
  fraction and mean strain, which in their textbook form involve `exp(z) * E_1(z)` and underflow for
  small sliding velocity. Aleph does **not** inherit a conditioning trick for it: the integrals are
  re-derived in a shifted variable (§4) in which no exponential overflows at any velocity, which is
  a better form than a guarded special-function call, and is Aleph's own.

## 4. Physical or mathematical law represented

All of the following is derived here. Aleph units throughout: length µm, force pN, energy pN·µm,
time s, stiffness pN/µm, rate 1/s.

### 4.1 One head

A head is a linear elastic crossbridge of stiffness `k_xb` [pN/µm] anchored on the minifilament
backbone and, when attached, bound to one material coordinate of one actin filament. Let `q` be the
barbed-end-ward offset of the bound site from the head's anchor, and `q*(φ)` the value of `q` at
which the crossbridge is unstrained given lever-arm phase `φ ∈ [0,1]`. The working stroke moves the
unstrained point by the stroke displacement `d` [µm] in the pointed-end-ward sense:

    q*(φ) = q*(0) − d·φ

Define the head strain `x ≡ q − q*(φ)`. Then the force the head applies to the bound actin site is

    f_site = −k_xb · x · t̂          (t̂ = unit actin tangent, barbed-end-ward)

and the reaction on the backbone anchor is `+k_xb·x·t̂`. Positive `x` therefore pulls the actin
toward its own pointed end, which is the direction myosin drives it. Over a step,

    x_new = x_old + d·Δφ + Δs,      Δs ≡ t̂ · (Δr_site − Δr_anchor)

so shortening at sliding velocity `v` (site moving pointed-end-ward relative to the anchor, `Δs =
−v·Δt`) gives the strain equation

    dx/dt = d·φ̇ − v                                                     (L1)

### 4.2 Load-dependent detachment — the law and why it is Bell-type

Detachment is a bond rupture across a single transition state, so its rate carries an Arrhenius
factor in the mechanical work done on the transition state, `k = k_0 · exp(±F·δ/k_B T)`. That is the
Bell form, and it is chosen because it is the *only* one-parameter form that follows from putting a
force into an activation barrier; anything else would be a fit.

The sign is the content. For myosin II the rate-limiting detachment step is slowed by a load that
resists the stroke, which is what lets the motor bear force. Written in the strain variable:

    k_det(x) = k_det_0 · exp(−x / x_bell),        x_bell ≡ f_bell / k_xb                (L2)

with `f_bell` [pN] the Bell force scale and `x_bell` [µm] the equivalent Bell strain. A head
strained in its productive sense (`x > 0`) hangs on; a head dragged past its neutral point (`x < 0`),
which is a head that has become a brake, lets go fast. `f_bell = k_B T / δ` if one prefers to
parameterise by a transition-state distance; Aleph provides that conversion as a helper and refuses
to apply it with an implicit temperature, per PLAN §6 (`k_B·T` inherits `ASSUMED` from the
temperature, so no entry point may default to it).

Attachment is load-free at rate `k_att` [1/s] per detached head, scaled by that head's activation
state, and places the head unstrained (`x = 0`, `φ = 0`) — weak binding precedes the stroke.

### 4.3 The ensemble force–velocity relation is *derived*, not imposed

No Hill hyperbola is assumed anywhere. Let `n(x)` be the steady-state density of attached heads per
head in the ensemble at strain `x`, with the stroke fast compared with detachment so that heads are
injected at `x = d`. Steady shortening at `v > 0` in (L1) gives, for `x < d`,

    v · dn/dx = k_det(x) · n(x)

with injection flux `v·n(d⁻) = k_att·p_detached`. Hence

    n(x) = (k_att·p_d / v) · exp( −(1/v) ∫_x^d k_det(x′) dx′ )                          (L3)

With (L2), substituting `y = exp(−x/x_bell)` and then `s = λ(y − y_d)` where `λ = k_det_0·x_bell/v`
and `y_d = exp(−d/x_bell)`, and writing `z ≡ λ·y_d = k_det(d)·x_bell/v`:

    p_b / p_d = (k_att / k_det(d)) · z·A₀(z)
    ⟨x⟩       = d − x_bell · Φ(z)

    A₀(z) = ∫₀^∞ e^{−s} / (z + s) ds
    Φ(z)  = [ ∫₀^∞ e^{−s} · log1p(s/z) / (z + s) ds ] / A₀(z)                           (L4)

Both integrals are bounded, positive, and free of any growing exponential at every `v > 0`, which is
the conditioning claim of §3. The ensemble force on the actin follows as

    F(v) = N_heads · k_xb · p_b(v) · ⟨x⟩(v)                                              (L5)

**Stall force.** As `v → 0⁺`, `z → ∞`, `z·A₀(z) → 1` and `Φ(z) → 1/z → 0`, so

    F_stall = N · k_xb · d · ρ_duty,      ρ_duty = k_att / (k_att + k_det(d))             (L6)

`ρ_duty` is the isometric duty ratio and `k_xb·d·ρ_duty` is the **single-head stall force** [pN].
(L6) is the exact statement of ensemble additivity at stall and is what the ensemble control
measures.

**Unloaded sliding velocity.** `v_0` is the root of `⟨x⟩ = 0`, i.e. of `Φ(z) = d/x_bell`. `Φ` is
strictly decreasing, so the root is unique; `v_0 = k_det(d)·x_bell / z_root` [µm/s].

**Shape.** (L4)–(L5) are not a Hill hyperbola. They are **convex** over `[0, v_0]`, as Hill's
hyperbola is, so Aleph reports the best-fit Hill curvature `a/F_stall` as a **measured** number
rather than accepting it as an input. Writing the hyperbola with the stall force and the zero-force
velocity substituted in leaves a single shape parameter,

    F(v)/F_stall = kappa * (1 - v/v_0) / (v/v_0 + kappa)

so `kappa` is a statement about curvature alone and cannot absorb an error in either endpoint. A
three-parameter Hill fit will match almost any decreasing curve, and reporting that it matched would
be reporting nothing.

This is the one place where the conventional choice (assume Hill, fit `a` and `b`) is deliberately
not made. The hyperbola is an empirical summary of whole-muscle data; adopting it at the head level
would also specify the same degree of freedom twice, since a head that carries an elastic crossbridge
already has its velocity–force relation determined by (L1) and (L2).

**Measured, 2026-07-30, at the fixture parameters:** the derived curve fits `kappa = 0.19299` with a
worst point deviation of `5.8e-02`, and the curve measured from the stepped kernel fits
`kappa = 0.19285` with `6.0e-02`. The two curvatures agree to `7e-4` relative. That the derived law
lands near Hill's classic value without being asked to is a result, not an input — and the sign took
a measurement to establish: the first draft of this entry asserted concavity and was wrong.

### 4.4 Why load-dependent detachment is load-bearing, stated as a derivation

Set `k_det(x) = k_det(d)` — the same *isometric* rate, load dependence deleted. Then (L3) integrates
in closed form and

    p_b = k_att / (k_att + k_det(d))        — independent of v
    ⟨x⟩ = d − v / k_det(d)                  — F(v) exactly linear, unbounded below

The stall force (L6) and the isometric duty ratio are **unchanged**. What is lost is:

1. the velocity dependence of the bound fraction — `dp_b/dv ≡ 0`, exactly;
2. the convexity of `F(v)` — `d²F/dv² ≡ 0`, exactly, because the relation is then the straight line
   joining the stall force to the zero-force velocity (measured deviation from that line:
   `1.1e-16`);
3. boundedness of the strain — a head can be dragged to arbitrarily negative strain and keeps
   pushing, so the ensemble resists being pulled linearly and without limit.

Together those three say the broken ensemble is a **dashpot with an offset**, not a motor. Items 1
and 2 are exactly zero for the broken model and strictly non-zero for the correct one, which is why
they and not the stall force are what the negative control keys on.

### 4.5 The active-power ledger

Per attached head over one step, with `x̄ = (x_old + x_new)/2`:

    W  = −k_xb·x̄·Δs                          endpoint work done by the head on its endpoints
    dU = ½k_xb(x_new² − x_old²) = k_xb·x̄·(d·Δφ + Δs)
    A  = k_xb·x̄·d·Δφ                          ATP-derived stroke work
    D  = ½k_xb·x_detach²                       stored strain released when a head lets go

Then `W + dU + D − A = 0` identically, term by term, with `A` computed from the stroke increment and
the strain — **never** from the residual. This is the identity the module asserts to round-off, and
it is the reason `A` is a first-class field rather than whatever is left over.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| crossbridge stiffness `k_xb` | pN/µm | N/m | `> 0`, finite |
| working stroke `d` | µm | m | `> 0`, finite |
| attachment rate `k_att` | 1/s | 1/s | `> 0`, finite |
| unloaded detachment rate `k_det_0` | 1/s | 1/s | `> 0`, finite |
| Bell force scale `f_bell` | pN | N | `> 0`, finite |
| stroke rate | 1/s | 1/s | `> 0`, finite |
| head strain `x` | µm | m | unbounded sign; magnitudes ≲ a few `x_bell` |
| head force | pN | N | `k_xb·x` |
| sliding velocity `v` | µm/s | m/s | `≥ 0` for the derived F–v; `0` handled by the limit (L6) |
| stall force | pN | N | `> 0` |
| active power `A/Δt` | pN·µm/s | W | `≥ 0` while strokes are being executed |
| contractile dipole | pN·µm | N·m | sign is the observable |

Singular and boundary cases, each with the behaviour Aleph requires:

- `v = 0` exactly: (L4) has `z = ∞`. The module takes the analytic limit (L6) rather than dividing
  by zero, and the limit is checked against `v → 0⁺` by convergence, not asserted.
- `Δx = 0` over a step (no sliding, no stroke): the exact interval hazard `H` in the head update has
  a `0/0`; the module uses the `Δx → 0` limit `H = k_det(x)·Δt` and a test drives that branch.
- Zero heads, zero minifilaments, zero binding sites: legal. Every aggregate is `0.0` and no
  division by a population count occurs.
- A minifilament with **no bound heads at all**: legal and load-bearing for
  `nmii_cytosol_transfer` — it still has a backbone that feels drag, so the drag quadrature must
  include it. Asserted.
- Two heads proposing the same actin material coordinate in one candidate: refused by the
  reservation set; the second head does not bind. Site exclusion is nmii's own binding topology, not
  a write into the actin owner's arrays.
- Non-unit actin tangent: refused (`NmiiGeometryError`). The tangent *is* the polarity; a
  non-normalised one silently rescales every head force.

Invariants that must hold, each with the test that asserts it:

- **I1.** The module's state blocks cover exactly the eleven `owned_state` keys of the registered
  `nmii` contract, with no key claimed twice.
  `tests/vertical/test_nmii.py::test_state_blocks_match_the_registered_contract`
- **I2.** Ensemble additivity at stall: `F_stall(N) = N · f_stall(1)`.
  `tests/vertical/test_nmii.py::test_ensemble_stall_is_n_times_single_head_stall`
- **I3.** The measured steady-state bound fraction equals `ρ_duty` implied by the rates.
  `tests/vertical/test_nmii.py::test_measured_duty_ratio_matches_the_rate_ratio`
- **I4.** The measured force–velocity curve matches the derived form (L4)–(L5) in *shape*, not only
  in monotonicity. `tests/vertical/test_nmii.py::test_force_velocity_curve_matches_the_derived_form`
- **I5.** `W + dU + D − A = 0` to round-off, with `A` independently computed.
  `tests/vertical/test_nmii.py::test_active_power_ledger_closes_with_independent_atp_term`
- **I6.** A rejected candidate leaves every head array and the RNG stream position bit-identical,
  *and the candidate demonstrably changed attachment state before being rejected*.
  `tests/vertical/test_nmii.py::test_rejected_candidate_restores_every_head_and_the_rng_exactly`
- **I7.** Newton's third law holds per head and in the resultant, with the two sides accumulated
  independently. `tests/vertical/test_nmii.py::test_motor_connector_is_adjoint_per_head`
- **I8.** The backbone is bipolar: two head groups on opposite sides of the backbone centre, and the
  contractile dipole is positive while the net force on the backbone vanishes.
  `tests/vertical/test_nmii.py::test_bipolar_backbone_is_contractile_not_a_transporter`

## 6. Source evidence class and known retractions

Where I looked: the motor directory's parameter file with its per-parameter evidence fields, its
ensemble-stall and Bell-kinetics modules, its force-budget module, and the audit record already in
this repository at `PLAN.md` §1.1 and §7 — which was produced by *running* the reference rather than
reading it. A subagent audit of the remaining motor modules was dispatched and died on a platform
error before returning, so the reading below is only as broad as it says it is: I did not read
`segment_motor.py`, `minifilament_warp.py`, `backbone_warp.py`, `two_filament_reference.py` or
`resting_setpoint.py`, and this entry makes no claim about them.

**Its own evidence labels are the most useful thing in the tree, and they are honest.** The parameter
file classifies nearly every magnitude the model needs as an open gap awaiting the PI, with the value
held at null rather than filled in: the per-head stall force, the head count per side, the unloaded
velocity, the force–velocity curvature, the crossbridge stiffness, the attachment rate, the backbone
length and bead count, the working-stroke displacement and the ATP free energy are all null. Two
competing literature claims are recorded side by side for several of them rather than one being
silently chosen. Three positions are recorded for the force–velocity curvature — a muscle hyperbola,
a non-muscle single-molecule hyperbola, and a ratified position that the non-muscle law is linear —
and the conflict is left open. The file states as a hard rule that no value may be tuned to pass a
gate.

Two consequences for Aleph, and they point in opposite directions:

- Nothing was available to inherit even if the policy allowed it. **No parameter value was taken.**
  Aleph's module ships no numeric default at all — see §14.
- The one number the file does carry as a "physical constant" is the thermal energy at body
  temperature. It is not a physical constant: it embeds a temperature convention, and PLAN §6 already
  records the cost of that exact confusion elsewhere in the reference. Aleph does not have a
  `k_B·T` default anywhere, and this module's Bell scale is a caller-supplied force.

**One defect found by reading, and it is the defect this lane's positive control exists to avoid.**
The ensemble-stall module offers a mean-field route and a "stochastic" route and states that their
agreement demonstrates the ensemble force is assembled from heads rather than imposed. It does not.
The stochastic route draws a bound-head count from a binomial distribution whose success probability
is the *same* engaged fraction the mean-field route uses, and multiplies it by the *same* imposed
per-head stall force. Its mean is the mean-field value identically, by construction — the comparison
is a binomial mean against its own parameter. Neither route contains a crossbridge, a strain, a
working stroke, or time; the engaged fraction is a Bell factor evaluated at one imposed force rather
than an average over a strain distribution. So the quantity described as emergent is the bound
*count*, and the per-head force it multiplies is a free parameter that the file itself records as an
open gap with two conflicting candidate values.

That is why Aleph's additivity control (§8) measures the single-head stall force from an independent
long single-head run and compares it against an independently measured ensemble, and why the
single-head stall force in §4 is *derived* as `k_xb·d·ρ_duty` from three parameters that have
separate meanings, rather than being a parameter of its own.

**A second, structural divergence, recorded as a disagreement rather than as a defect.** The
reference imposes a Hill hyperbola on each individual head's stepping velocity and separately applies
a Bell *slip* factor in the absolute value of the head force. Aleph does neither: the force–velocity
relation is derived for the ensemble from the head kinetics (§4.3), and the Bell factor is a catch
form in the *signed* strain (§4.2). Two reasons, both physical. Imposing a velocity–force law on a
head that also carries an elastic crossbridge specifies the same degree of freedom twice. And taking
the absolute value of the force gives a head that has been dragged backwards into a brake exactly the
same detachment rate as a head pulling productively — which removes precisely the asymmetry that lets
an ensemble shed its brakes and hold isometric load. Aleph's negative control (§9) keys on that
asymmetry, so this module would fail a check built on `|f|`.

Reachability and evidence class: the motor tree is reachable from the reference's production path and
carries live tests. That is a statement about the code, not about the numbers. On Aleph's ladder every
magnitude in it is `INHERITED_UNVERIFIED` at best and most are absent, and its own file says the
native magnitude gate is invalid until the gaps close. Aleph inherits none of it.

## 7. Independent oracle or derivation

Three routes, mutually independent, none of them the reference:

1. **The derivation of §4**, evaluated as (L4)–(L6) by quadrature in the shifted variable. This is
   the mean-field prediction and it lives in `aleph/vertical/nmii.py` as the model's own statement
   about itself.
2. **An exact event-driven renewal estimator**, written in the test from the *analytic* inverse
   cumulative hazard of (L2) — for a head attached at strain `x₀` and sliding at constant `v`, the
   detachment time solves `τ = (x_bell/v)·log1p(E·v·e^{x₀/x_bell}/(k_det_0·x_bell))` for a unit
   exponential `E`. It has **no timestep** and therefore no discretisation error; it differs from (1)
   only by Monte-Carlo noise, and it shares no code with (1).
3. **The production head arrays themselves**, stepped through real `Transaction` candidates against
   a rigid sliding actin track. This is the path that has to be right, and (1) and (2) are what judge
   it.

Two exact identities are additionally available and are used as oracles in their own right: the
energy identity of §4.5, which must close to round-off at any parameter values, and the `v → 0` limit
(L6), which must be reached by (1) as `v` decreases and must equal the isometric duty ratio computed
from the rates alone.

## 8. Positive control

| Control | Location | Asserts |
|---|---|---|
| Positive | `tests/vertical/test_nmii.py::test_ensemble_stall_is_n_times_single_head_stall` | Ensemble additivity at stall: `F_stall(64)/(64·f_stall(1))` within the Monte-Carlo band of the two independent measurements, with the band computed from the event counts rather than chosen |
| Positive | `tests/vertical/test_nmii.py::test_force_velocity_curve_matches_the_derived_form` | The measured curve agrees with (L4)–(L5) point by point, and the *shape* agreement is asserted through the fitted Hill curvature rather than through monotonicity alone |
| Positive | `tests/vertical/test_nmii.py::test_measured_duty_ratio_matches_the_rate_ratio` | Measured bound fraction equals `k_att/(k_att+k_det(d))` |
| Positive | `tests/vertical/test_nmii.py::test_active_power_ledger_closes_with_independent_atp_term` | `\|W + dU + D − A\|` relative to the sum of term magnitudes is at round-off, over a multi-step run with attachments and detachments actually occurring |
| Positive | `tests/vertical/test_nmii.py::test_mean_field_quadrature_agrees_with_the_event_driven_renewal_estimator` | Routes (1) and (2) of §7 agree, which is what licenses using (1) to judge (3) |

## 9. Deliberately failing negative control

| Control | Location | Asserts |
|---|---|---|
| Negative (must fail) | `tests/vertical/test_nmii.py::test_load_independent_detachment_is_caught` | With `DetachmentLaw.LOAD_INDEPENDENT` — the same isometric rate, load dependence deleted — `assert_load_dependent_detachment` refuses. The two signatures it reads are exactly zero for the broken model: `dp_b/dv = 0` and `d²F/dv² = 0`. The test asserts the refusal is raised, names both measured signatures, and separately asserts that the *stall force and isometric duty ratio are unchanged*, so the control cannot be passing for the wrong reason |
| Negative (must fail) | `tests/vertical/test_nmii.py::test_the_load_dependence_check_passes_on_the_correct_law` | The paired half of the above: the same check must **not** refuse the Bell law, or it is a check that refuses everything |
| Negative (must fail) | `tests/vertical/test_nmii.py::test_unipolar_assembly_is_caught_being_a_transporter` | A backbone whose two head groups sit on the same side of its centre is detected: the contractile dipole collapses and a net backbone force appears. Bipolarity is what makes the element contractile rather than a transporter, so this is the control for that claim |
| Negative (must fail) | `tests/vertical/test_nmii.py::test_kinetics_committed_during_a_rejected_candidate_is_caught` | Deliberately commits binding during a rejected candidate and shows the integrity assertion fails. Without this the rejected-step control could pass by never having bound anything |
| Negative (must fail) | `tests/vertical/test_nmii.py::test_unsourced_parameter_is_refused_where_a_source_is_required` | A parameter carried at `ASSUMED` cannot be used where the census policy requires a source; the refusal names the parameter |

## 10. Numerical and precision envelope

Working and accumulation precision: float64 throughout, on CPU. No float32 path exists in this
module and none is planned before the backend question is settled.

- **The energy identity (§4.5).** Asserted at `1e-12` relative to the sum of the term magnitudes,
  not to the residual. That number is chosen because the identity is algebraically exact: the only
  error is the round-off of a few thousand `fsum`-accumulated products, whose relative floor is a few
  ULP times the term count. `1e-12` is roughly three decades above that floor and roughly four
  decades below any modelling error, so it is a real gate rather than a rubber stamp. A looser
  tolerance would stop distinguishing "the books close" from "the books nearly close", which is the
  whole content of the check.
- **The mean-field quadrature (L4).** `scipy.integrate.quad` on `[0, ∞)` with absolute and relative
  tolerances at `1e-12`; the integrands are smooth, positive and exponentially damped, so the
  reported quadrature error is ~`1e-14` relative across `z ∈ [1e-3, 1e12]`. Outside that range the
  module refuses rather than degrading: `z > 1e12` means the velocity is so far below stall that the
  `v = 0` limit (L6) is the correct call and is used instead.
- **The Monte-Carlo comparisons.** These are the loose ones and their tolerances are *derived from
  the event counts*, not chosen: the standard error of a bound-fraction estimate from `n`
  attachment episodes is `sqrt(p(1−p)/n)`, and every stochastic assertion is stated at a stated
  number of standard errors of that quantity, with the event count reported alongside. A stochastic
  control asserted at a hand-picked absolute tolerance is a control that will be widened the first
  time it flakes, so the tolerance is computed instead.
- **Conditioning as written.** The one hazard is `exp(−x/x_bell)` for large negative `x`, which
  overflows for `x < −700·x_bell`. Strains that deep are unreachable under the derived dynamics
  (the double-exponential tail of (L3) cuts them off), but "unreachable under correct dynamics" is
  not a guarantee, so the detachment rate is evaluated through a guarded exponent that refuses above
  a stated exponent bound rather than returning `inf`.

## 11. Production-backend residency and transfer

Host, CPU, numpy float64, for now — and that is a decision with a reason, not a placeholder. Every
array in this module is host-resident and every operation goes through the `Backend` protocol, so
the numerical definitions live in `NumpyBackend` and a device backend substitutes for it rather than
reimplementing it.

What would move to a device and what would not, when a device is authorised: the per-head force
evaluation, the strain update and the per-head hazard evaluation are element-wise over the head
arrays and are the natural device kernels. The reservation set that enforces one head per actin
material coordinate is a scatter with conflict detection and is the part that needs care. The
mean-field quadrature (L4) is a diagnostic and a prediction, evaluated a handful of times per run on
the host, and must **never** be inside a step. No host round-trip is required per step by anything
in this module.

Zero GPU work was run for this entry, and no GPU authorization exists or was sought — PLAN §0.1.

## 12. Comments and docstrings to discard

Nothing to discard, because nothing was copied: no comment, docstring, identifier, parameter value
or file layout crossed the boundary. What is being deliberately *not* carried across is listed so a
reviewer can check for its absence:

- the reference's package namespace, module paths and file names;
- its gate labels and status vocabulary, and any claim inherited from them;
- its parameter file's values and their secondhand attributions — Aleph ships no numeric default at
  all, so there is nothing for an unread citation to attach to;
- its naming for the head/backbone data layout, including its choice of which quantities are
  per-head and which are per-minifilament;
- its arrangement of the force budget, and in particular any route by which an active-work term is
  obtained as a residual of the other three.

What replaces them: the derivation of §4 restated in the module docstring in Aleph's units, the
eleven state-block names taken from Aleph's own registered `nmii` contract in
`aleph/state/census_actomyosin.py`, and comments written for this module's situation — most of them
explaining why a check is arranged the way it is rather than what a line does.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | Not yet accepted. All 28 controls in `tests/vertical/test_nmii.py` pass at the tolerances of §10 (2026-07-30). Measured: ensemble stall ratio `1.00102 ± 0.00418` at N=128 and `1.00567 ± 0.00207` at N=512, per-head stall force varying by 2.3% across N ∈ {8, 32, 128}; duty ratio `0.384722 ± 0.001988` against `0.380767` implied by the rates, ratio `1.01039` over 59,880 attachment events; worst per-step power-ledger residual `6.4e-17` relative; force–velocity agreeing with the derived form at all eight velocities within the stated error bars, Hill curvature `0.19285` measured against `0.19299` derived. But every parameter the module runs on is `ASSUMED` (§14), so no quantitative claim about NMII is supported and `ACCEPTED` would overstate what exists |
| Reviewer | Agent-proposed, unratified. Lane L26, 2026-07-30. No PI review has occurred |
| Rollback | Delete `aleph/vertical/nmii.py` and `tests/vertical/test_nmii.py` and set this entry to `REJECTED`. Nothing else in the tree imports either file, so nothing breaks except that the five connectors listed in §14 stay blocked and the seed census keeps its only active element unimplemented |

## 14. Honest limits

What this entry does **not** establish:

- **No parameter is sourced.** All six physical parameters are supplied by the caller with an
  explicit `EvidenceClass`, and every value used in the controls is `ASSUMED` with a source string
  that says so. Consequently the module establishes nothing about the magnitude of NMII
  contractility, the real duty ratio, the real stall force, or the real unloaded velocity. It
  establishes that the *machinery* is self-consistent and that the additivity, duty-ratio, F–V and
  energy identities hold, which is a statement about the model and not about myosin.
- **`k_B·T` is not used implicitly anywhere**, per PLAN §6. The Bell force scale is a caller-supplied
  parameter; the temperature route is offered as an explicit helper that requires the temperature and
  a transition-state distance to be passed, and inherits `ASSUMED` from the temperature.
- **The mean-field prediction assumes a fast stroke.** (L3) injects heads at `x = d`. With
  `StrokeModel.FIRST_ORDER` and a stroke rate comparable to the detachment rate, the measured stall
  sits below (L6) by a bias of order `k_det/k_stroke`. The controls that compare against (L6) use
  `StrokeModel.INSTANT`, and the size of the finite-stroke-rate bias is measured and reported rather
  than assumed small.
- **The mean-field prediction is a mean field.** It ignores the correlation between heads that share
  one backbone and one bound filament. In the isokinetic assay that correlation is absent by
  construction, which is why the assay is the right instrument for comparing against it — and it is
  also why agreement in the assay is *not* evidence about a compliant, many-filament network.
- **Nothing here is a whole-cell result, or evidence about one.** No connector to `sf_arc`, `cortex`,
  `lamellipodium` or `filopodium` is wired into an assembly by this entry. What lands is one generic
  `NmiiMotorConnector` parameterised by target owner, plus the binding-site interface the four
  `[K]` motor connectors and the `[C]` `nmii_cytosol_transfer` drag quadrature need in order to be
  written. Those five remain unbuilt.
- **The actin-side interface is Aleph-internal and provisional.** `ActinBindingSites` is defined in
  this module because the four actin owners are being written in parallel lanes and depending on
  their APIs would couple this lane to theirs. Reconciling the two is future work and may change the
  shape of `ActinBindingSites`.
- **The runtime's `AdjointPair` energy channel cannot carry this connector's work.** It computes
  endpoint work as `F·Δx` on the reported *resultants*, and a correctly contractile bipolar motor
  has a vanishing actin-side force resultant by symmetry — so the resultant-based work is not the
  distributed work, and no work-equivalent reduction exists at the point where the resultant goes to
  zero. The motor's power ledger is therefore a separate first-class object,
  `MotorPowerBudget`, and the pair reported into `ForceWorkLedger` carries forces only. This is a
  shape limitation of the runtime contract, recorded here rather than worked around; consequently
  `BalanceVerdict.total_active_work` does **not** see this connector's ATP term.
- **No load-history or catch-slip crossover.** (L2) is a pure catch bond in the strain variable over
  the whole range. Real myosin II detachment is non-monotone in load at large stretch. The module
  refuses above a stated exponent bound rather than extrapolating, and the non-monotone regime is
  simply not modelled.
- **Minifilament assembly and disassembly are not modelled.** The contract classes them as topology
  events on this owner; the module carries the topology arrays and the accepted-topology commit path
  but proposes no assembly events.
- **Site assignment does not scale, and this is a known cost rather than a hidden one.** The
  exclusion rule that stops two heads sharing one actin material coordinate is enforced by a linear
  scan over free sites per attaching head, so binding is `O(H·S)` in the worst case. That is fine at
  the sizes the controls run and it is wrong at R5 scale. The correct structure is a free-list per
  target with a device-side scatter and conflict detection, which is named in §11 as the part that
  needs care. Recorded here because a performance trap that nobody wrote down becomes a correctness
  compromise later, when somebody quietly relaxes the exclusion to make it fast.
- **Isoform, activation and spatial-population labels exist and carry no sourced kinetics.**
  Activation gates the attachment rate multiplicatively, which is a modelling choice with no
  per-isoform evidence behind it, and the module says so at the field.
