# Before building: are the parameters enough, and are the detectors wired?

**Status:** AUDIT. PI-requested 2026-08-20 as a parallel stream to the arena port —
*"전부 구축하기 전에 파라미터는 충분한지 … 우리 검출기들은 충분하고 그들에 대한 세팅은 다 잘 되어 있는
건지 병렬로"*. Measured, not surveyed. Decides nothing.

---

## 1. Detectors — 5 of 7 are wired to nothing

`aleph/observe/` holds **4,888 lines** of observation stack. Call sites counted across
`aleph/scripts/**` and `aleph/engine/**` (archive excluded):

| module | lines | call sites | what it is |
|---|---:|---:|---|
| `stationarity` | 1,030 | **3** | whether a series supports quoting a mean |
| `operator` | 493 | **3** | the observation seam: accepted state in, apparent quantity out |
| `fluctuation` | 904 | **0** | passive membrane thermal fluctuation spectrum |
| `tether` | 682 | **0** | **the AFM membrane-tether forward operator** |
| `kinematics` | 728 | **0** | spreading rate, volume flux, PIV speed |
| `traction` | 619 | **0** | traction force per cell, from TFM |
| `manifest` | 432 | **0** | what must be declared before a measurement means anything |

⚠ **3,365 lines — five of the seven — are called from nowhere.**

⚠ ~~`tether` is the one that matters most right now.~~ **CORRECTED 2026-08-20 by session C, and the
correction is worse than the original claim.**

`tether.py` reads a **bead-pulling tether plateau**. The AFM apparatus this project actually has —
`outputs/ac/afm_cortical_tension_readiness/` — is **spherical indentation**. They are different
experiments, so `tether` was never the forward operator for the PI's 3-minute protocol.

**And there is no indentation forward operator anywhere in `aleph/observe/`** (session C, grep). So the
protocol the entire 6,514x throughput budget is sized for does not have a detector that is *unwired* —
**it has no detector at all.** This audit pointed at the wrong module and thereby understated the gap.

⚠ Session C found a second thing in the same pass: **today's tether force is an echo of three
constants.** `gamma_mem` is a constant plateau (`compartments.py:579`, `K_A` off by default), `kappa` is
constant, `GAMMA_MCA_PN_UM` is a literal — so the scatter is zero and `standard_error` came out
**0.0**, which is a claim of infinite precision. The `ApparentQuantity` docstring already forbade it and
nothing enforced it; C's guard now rejects `standard_error == 0.0` across **all seven operators**, and
the driver reports `INPUT_ECHO_NOT_A_MEASUREMENT`.

⚠ **`manifest` being unwired is a second-order defect.** Its own docstring: *"A number without its
protocol is not a measurement, it is a rumour ... Thirty-four fields in six families. None of them has
a default."* It distinguishes `NOT_APPLICABLE` (the protocol has no such quantity — an assertion) from
`NOT_RECORDED` (it existed and nobody wrote it down — a confession). Nothing calls it, so nothing is
currently forced to declare either.

**This is not a request to wire all five.** It is the statement that the detector side of the port is
not "already done and just needs porting" — for five of seven there is nothing downstream to port TO.

---

## 2. Parameters — the gate is green and the coverage is not

`verify_params.py --gate`:

    declaration gate OK - 43 declared constants, none drifted
    coverage ratchet    - 71 KU-tagged runtime constants in aleph/configs/**
                          declared in params_manifest.yaml: 1
                          UNDECLARED (ratchet debt): 70   <-- may only SHRINK
                          42/43 declared rows point OUTSIDE aleph/configs/**

⚠ **The declaration gate and the coverage ratchet are measuring different things, and only the first is
green.** 43 constants are declared and none has drifted — but **42 of those 43 point at retired oracle
configs that no runtime path reads**, so the gate is clean over a set that is almost entirely not the
runtime. Of the 71 constants a runtime actually reads, **1 is declared and 70 are debt.**

## 3. PI-GAPs that will surface during PHASE 1

These refuse to default, so the build stops rather than inventing. Found in
`components/{incumbent,motor,solid}/`:

| where | constant | note carried in the source |
|---|---|---|
| `assemble.py:164-166` | `NMII_N_SIDE` = 10, `NMII_L_BB_UM` = 0.301, `NMII_HEAD_OFFSET_UM` = 0.200 | all three labelled PI GAP; `N_side` is the 10-vs-28/30 conflict |
| `driver.py:1592,1597` | resting bound-myosin fraction, per-head force | PI-GAP, **no default** |
| `compartments.py:102-113` | nuclear `r_eq_um` 5.0, `aspect` 1.0, lamina params, `eps_rupture` 0.50, `k_linc` 1.0e2 | every one *"I0-B2 PI GAP, TEST values"*; `k_linc` carries *"8 pN is a TENSION not a stiffness"* |
| `params_i0b3.yaml` | `F_stall_head` 0.5 vs 2.0, `N_side` 10 vs 28–30 | CONFLICTED; the 2026-08-20 source audit finds **neither end of the force bracket is a measurement** |
| ERM | `erm_density_per_um2` | `ABSENT_FROM_CONTRACT_GRAPH`; the fallback is one tether per mesh vertex |

**None of these blocks PHASE 1's geometry**, except where a density sets a node count — which is exactly
where the build should refuse. They block PHASE 2 and PHASE 3, and they are the reason PHASE 3 is a
decision queue rather than a coding queue.

---

## 4. What this changes about the plan

1. **PHASE 3's 28 connectors are gated on PI decisions, not on implementation time.** Each family must
   answer `BondCount`'s four questions, and the answers for ERM, NMII and LINC are PI-GAPs today.
   Parallel sessions can build the scaffolding; they cannot invent the counts.
2. **A detector stream belongs beside the port, not after it.** `tether` in particular: the throughput
   target is sized from an AFM protocol whose forward operator is wired to nothing.
3. **The params coverage ratchet is measuring the wrong set** — 42 of 43 declared rows are outside the
   runtime. Worth surfacing to the PI as its own item; it is not a blocker for PHASE 1.
