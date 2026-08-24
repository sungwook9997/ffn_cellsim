# Port ledger — ALEPH-PORT-304 · force / moment / work ledger

| Field | Value |
|---|---|
| Aleph target | `aleph/runtime/ledger.py` |
| Lane | L3 runtime |
| Written | 2026-07-30 (before the code landed, per PLAN §0.2.5) |
| Status | LANDED |
| Evidence class | ANALYTIC for the closure relations (exact for a linear spring + dashpot under the midpoint rule); STRUCTURAL for the reporting surface. |

## What was read in `ffn_cellsim` (read-only)

- `ffn_sim/ac/engine/ledger.py` (603 lines) — `GlobalCellLedger`, its two separately sourced force
  channels (`reaction_resultant`, `traction_resultant`), `assemble_balance`, and
  `derive_balance_tolerance` / `assemble_balance_tolerance_ratio`.

## What is genuinely good there, and was carried as an *idea*

Three ideas, restated in Aleph's own words and re-implemented from scratch:

1. **Two channels must be independently sourced.** The reference module prose is exactly right that
   computing one channel from the other (`traction := −reaction`) converts a check into a tautology.
   Aleph carries this by having each connector report `force_a` and `force_b` from separate
   accumulations, and by documenting the requirement on `AdjointPair` rather than in a module header.
2. **The tolerance is not a magic number.** The reference derives it from a float64 summation bound
   scaled by the ledger's own resultant magnitudes. Aleph's `Tolerance` is the same shape of idea in
   simpler form — `absolute + relative * scale`, with the scale built from the magnitudes of the terms
   that entered the residual, per channel, per connector.
3. **The gate must not decide acceptance.** The ledger reports; the transaction decides.

## What was rejected

- The reference collapses to `balance_ok_d`, a single device flag. Aleph returns a `BalanceVerdict`
  with a `ConnectorResidual` per connector. A global sum can be small while two connectors are wrong
  in opposite directions — which is precisely the situation where a whole-cell gate reports health it
  does not have.
- Warp/CUDA residency. All of this is float64 on the host, per PLAN §0.1.

## Physics re-derived independently

Let a connector apply step-average forces `F_a`, `F_b` at points `r_a`, `r_b`, with endpoint
displacements `dx_a`, `dx_b` and declared couples `M_a`, `M_b`.

**Force closure** (Newton III for a complete momentum pair): `F_a + F_b = 0`.

**Moment closure** about the centroid `c = (r_a + r_b)/2`:
`(r_a − c) × F_a + (r_b − c) × F_b + M_a + M_b = 0`. Substituting force closure gives
`(r_a − r_b) × F_a + M_a + M_b`, which vanishes iff the transmitted force is central. This is the
two-force-member theorem and is **not** implied by force closure — it catches a force of the right
magnitude and sign pointing the wrong way.

**Energy closure.** Let `W = F_a·dx_a + F_b·dx_b` be the work the connector does *on* its endpoints,
`dU` the change in its internal stored energy, `D ≥ 0` its dissipation, `A` the energy drawn from a
chemical reservoir (0 if passive). Conservation for the connector as a subsystem:

    W + dU + D − A = 0

**The audit note this ledger was written to honour:** for a passive connector this says
`W = −(dU + D)`, which is *not generally zero*. A relaxing spring does positive work on its endpoints
and pays out of `dU`; a dashpot does negative work and puts it into `D`. A test asserting `W = 0`
would pass only for a connector that did nothing — the failure mode this whole lane exists to prevent.
So the tests assert the closure relation *and separately assert that `W ≠ 0`*.

**Why the midpoint rule closes exactly for the positive control.** For a linear spring of stiffness
`k` and rest length `L0` extending from `l1` to `l2`, the step-average force `−k((l1+l2)/2 − L0)` times
the displacement `(l2 − l1)` equals `−½k[(l2−L0)² − (l1−L0)²] = −dU` identically. For a dashpot of
coefficient `c` at constant relative velocity `v`, `W = −c|v|²dt = −D` identically. The positive
control therefore closes to machine precision, not to a fudged tolerance.

## Controls

- Positive control: `test_spring_dashpot_pair_closes_energy_books` (residual at float64 round-off, and
  `endpoint_work` asserted **non-zero** so the test cannot pass vacuously);
  `test_active_connector_keeps_its_power_channel_separate`.
- Negative (must fail): `test_broken_third_law_fails_force_closure`,
  `test_non_central_force_fails_moment_closure`,
  `test_missing_dissipation_fails_energy_closure`,
  `test_negative_dissipation_is_rejected`,
  `test_duplicate_connector_contribution_is_refused`.

## Provenance and envelope

- **Source repository / commit read:** `ffn_cellsim` at commit `be0e5876`, read-only.
- **Discarded prose:** all original comments and docstrings discarded. The three ideas listed above
  were re-stated from the physics, not translated; the reference's Sanity Gate block, its
  ECM-far-field framing and its device-residency notes were all dropped as inapplicable.
- **Units, domains, invariants:** forces pN, points and displacements um, moments pN.um, energies
  pN.um. Invariants: `dissipation >= 0` (one-sided, checked separately from the balance so a sign
  error cannot hide in it); one connector may contribute at most once per candidate; the three energy
  channels are never netted before the check.
- **Numerical and precision envelope:** float64; aggregates use `math.fsum` so a long connector list
  does not accumulate summation error into the totals. The positive control closes to float64
  round-off *identically* rather than approximately, because the midpoint force rule is exact for a
  linear spring and for a dashpot at constant relative velocity.
- **Residency:** host.

## Independence

No `ffn_cellsim` import. No shared identifiers. The reference's ECM far-field cut is not modelled here
at all; Aleph's ledger is per-connector by construction.
