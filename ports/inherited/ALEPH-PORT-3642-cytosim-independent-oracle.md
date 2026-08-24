# ALEPH-PORT-3642 — an oracle that fails differently

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3642` |
| Lane | `46143f30` Lane W2 (cortex) |
| Status | `PROPOSED` |
| Written | `2026-08-06` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **EXTERNAL ORACLE, driven as a process.** No Cytosim source enters Aleph. The binary is executed and its output parsed. |
| Third party | Cytosim, Nédélec & Foethke — `https://gitlab.com/f-nedelec/cytosim`, **GPL v3**, built locally at commit `3e9a5fa` |
| Aleph target | `validation/independent/` — **never** `aleph/**` |
| Authorisation | `validation/independent/CITATION-cytosim-build.md`, PI 2026-08-06T11:34:11.600Z |
| Exists because | Every oracle this project owns is computed on the same discretisation as the thing it judges, and today it met that limit twice. |

---

## 0. Why an oracle that shares nothing

Aleph's internal oracles are good and they have a fixed blind spot: **they cannot catch a defect
they share.** Two from today:

* The bending sum omitted the two end half-segments, so it was `(n−2)/(n−1)` of the continuum and
  **2× low at the cortex's `n = 3`**. Every internal check used the same interior-triple sum, so none
  could see it. It surfaced because the *provider* had recorded the factor against Cytosim.
* The affine estimator reports **14,270 pN/µm** for a shell that is 2,436 disconnected fragments. It
  is blind to connectivity, and no internal cross-check reveals that, because every internal check
  runs on the same graph.

An independent implementation is the only instrument that fails differently. That is the whole
argument for spending a build on it.

## 1. Licence, and what it forbids

Cytosim is **GPL v3**. Checked at the clone rather than assumed.

| | |
|---|---|
| build and run it | fine |
| **drive the binary as a subprocess and compare outputs** | fine — no linking, no derived work |
| copy any Cytosim source into `aleph/**` or `validation/**` | **would make Aleph a derived work.** Forbidden here independently of the porting rule. |

So the subprocess form is not a convenience, it is the only admissible one, and this entry records
that the alternative was considered and is closed.

## 2. The first comparison, and why it is NOT the provider's

`ffn_sim/ff/cytosim_parity.py` used the **bending energy of a fixed circular arc**, and found a real
defect that way. Repeating it here would be worth little: `ALEPH-PORT-3639` has already anchored
Aleph's bending energy to the closed form `κL/2R²`, so an arc-energy comparison would mostly confirm
a term that already has an analytic check.

And one of this lane's parallel surveys measured why that observable is weak:

> the arc bending **energy** is only **0.79% sensitive at `n = 129`** to a change that alters the
> end-node bending **force by 50%** at every `n ≥ 5`, and by 100% at the cortex's `n = 3`.

**So the first comparison is the end-node bending force**, not the energy:

* it is **~63× more sensitive** to the thing `ALEPH-PORT-3639` changed;
* it is exactly where that port changed behaviour, and Aleph has **no closed form** for the discrete
  end-node force — the continuum has no "end node";
* it is one filament, so both codes are set up in a few lines and nothing about crosslinks,
  placement or density enters.

The configuration: a single fibre of length `L` held in a prescribed circular arc of radius `R`, with
`rigidity` = κ and `segmentation` = `L/(n−1)`, at `n` = 3, 5, 9, 17, 33. Cytosim's units are
**seconds, micrometres and piconewtons** — identical to Aleph's, so nothing is converted. Compare the
force on the **end node** and on its neighbour.

**A gate written by the survey, not by this lane, and it matters**: a naive "L-bend" configuration
gives a ratio of exactly 1.0 when the kink sits at an interior node and 1.5 only when it sits at
node 1 or `n−2`. **The comparison must state where the curvature is**, or it measures nothing.

## 3. What a disagreement means, and what it does not

**Cytosim is not authority.** If the two disagree, that is a finding and which side is wrong is a
separate measurement. The provider's own record is the precedent: it found Cytosim matching the
continuum and its own engine 17% low, *surfaced it as a fidelity choice*, and did not silently
change its law.

Three outcomes and what each licenses:

| | |
|---|---|
| they agree | Aleph's end-node force has an independent check for the first time. `ALEPH-PORT-3639` gains external support it does not currently have. |
| they disagree and Aleph reproduces the continuum in a limit Cytosim also reproduces | Cytosim's end treatment differs; record it, do not change Aleph. |
| they disagree and Aleph does **not** reproduce a limit both should | a defect in Aleph, found by the only instrument that could. |

## 4. The gates — written before the code

| # | Gate | Threshold |
|---|---|---|
| **E-1** | **The harness degrades.** With no binary present and `CYTOSIM_SIM` unset, every entry point returns a clear skip and the committed suite stays green. | `tests/validation` green on a machine that has never built Cytosim |
| **E-2** | **The firewall holds.** `aleph/**` does not import `validation/**`, and `validation/**` does not import Cytosim source. | `tests/firewall` green |
| **E-3** | **The same configuration reaches both codes.** The generated `.cym` and Aleph's builder agree on `L`, `R`, `κ`, `n` and the node positions to 1e-12 µm. | positions bitwise-comparable after unit-identity |
| **E-4** | **The end-node force comparison runs and reports a number.** At `n` = 3, 5, 9, 17, 33. | a table, not a verdict |
| **E-5** | **Aleph's end-corrected and interior-only quadratures are BOTH compared.** | both rows present |
| **E-6** | **The disagreement is reported, not resolved.** No Aleph law changes in this port. | `git diff` touches nothing under `aleph/` |

**E-4 carries no threshold and E-6 forbids acting on it.** This entry buys a measurement, not a
correction; whatever it finds gets its own ledger. Writing a threshold now would be choosing the
answer before the instrument exists, which is the error this project has recorded five times today.

## 5. What is NOT claimed

- **Not that Cytosim is correct.** It is *independent*, which is a different and weaker property,
  and it is the only one being relied on.
- **Not that this validates the cortex.** One observable on one filament. The areal modulus, the
  connectivity and the network response remain unchecked against anything external.
- **Not that the units claim is verified.** That Cytosim's configuration units are s/µm/pN is taken
  from the provider's parity record and from Cytosim's own documentation; this lane has not derived
  it, and E-3 is what would catch it being wrong.
- **Not a comparison of solvers.** Aleph relaxes; here both codes are given a *prescribed*
  configuration and asked only for forces. Nothing is integrated on either side.

---

## Amendment, 2026-08-07 — the harness exists and **E-4's observable is not readable out of Cytosim**

### What was there before today

The authorisation and the binary; **zero lines of harness**. `validation/independent/` held exactly
one file, `CITATION-cytosim-build.md`, and a grep for `cytosim` across `aleph/`, `tests/` and
`validation/` returned nothing. E-1 was vacuously true because no entry point existed to degrade.

And one gap that looked like a working install and was not: **only `sim` had been built.**
`report` — the one tool that emits per-vertex forces — is a separate cmake target and
`build/bin/` contained a single binary. `sim info` worked, so nothing said anything was missing.
`make report` in the Cytosim build tree fixes it and now `build/bin/report` exists.

### What exists now

`validation/independent/cytosim.py` and `tests/validation/test_cytosim_oracle.py`.

| # | outcome |
|---|---|
| **E-1** the harness degrades | **PASS** — 5 controls, all green with `CYTOSIM_BIN` pointed at an empty directory. A missing binary returns a `CytosimUnavailable` **value**, never an exception, and a build missing `report` says *"make report"* by name |
| **E-2** the firewall holds | **PASS for this port** — `test_runtime_never_imports_the_analytic_oracles` green; `aleph/**` imports `validation/**` nowhere. `tests/firewall` is still red for `aleph/scenarios/ecm_provider.py`, which is **not this row's file** and is reported, not patched |
| **E-3** both codes get the same configuration | **Aleph half done** — the arc is gated to `1e-14` on radius and `1e-12` on chord equality. The Cytosim half is blocked by E-4 |
| **E-4** the end-node force comparison | **BLOCKED, and the blocker is Cytosim's** |
| **E-5** both quadratures compared | not reached |
| **E-6** no Aleph law changes | **held** — `git diff` touches nothing under `aleph/` |

### E-4: `report fiber:force` returns exactly zero for a curved free fibre

Measured, on the build at commit `3e9a5fa`:

| configuration | rigidity | reported force |
|---|---:|---|
| 5-point circular arc, `L` = 1 µm, `R` = 2 µm | 0.075 pN·µm² | `+0 +0 +0` at every vertex |
| 5-point **90° kink** | **100** pN·µm² | `+0 +0 +0` at every vertex |
| the same arc at `precision=12` | 0.075 | `+0 +0 +0`, and `tension = 0` |

**It is not a formatting artefact** — `precision=12` shows the same zeros — and it is not a broken
configuration: the reported vertex positions are demonstrably *curved* (consecutive segment
directions turn), so Cytosim did receive and keep the shape.

`Simul::computeForces` (`src/sim/simul_solve.cc:313`) calls `Meca::calculateForces`, whose own
comment reads `// vFOR <- external forces`. **A free fibre with no confinement, no attachment and no
steric has no external force**, and its bending enters the solve as a matrix operator rather than as
an assembled force vector. So `fiber:force` appears to report what is *applied to* a fibre, not what
its own curvature does.

**`ALEPH-PORT-3642` §2 chose that observable because it is 63× more sensitive than the arc energy.
That reasoning stands. What does not stand is the assumption that Cytosim will hand it over.**

### Three things that would unblock it, none of them attempted here

1. **Confine or attach the fibre**, so the bending force appears as the reaction that holds the arc.
   Changes the comparison — both codes would then need the same constraint — and that is a design
   question, not a fix.
2. **Find the right report.** `fiber:tension` and `fiber:confine_force` exist; whether either
   exposes the internal bending is unread.
3. **Fall back to the arc energy**, which `-3642` §2 explicitly rejected as 63× less sensitive.
   Worth having as a *weak* check rather than none.

### What is NOT claimed

- **Not that Cytosim cannot do this.** It is that *this route* returns zeros, measured on three
  configurations. Somebody who knows Cytosim's reporting surface may have this in one line.
- **Not that the harness is wasted.** It builds the configuration, drives both binaries, parses the
  report **by column name** rather than by position, and degrades cleanly. When a readable force is
  found, only the reader changes.
- **Not that Aleph's end-node force is unverified because of this.** It is unverified *by an
  independent implementation*, which is what this port was for and what it has not delivered.

---

## Amendment, 2026-08-07 — **E-4 and E-5 delivered.** An implementation sharing nothing grades `-3639`

The end-node force is unreadable (previous amendment). `report fiber:energy` is readable, and it is
enough for the one thing this port existed to settle.

### The measurement

One fibre held in a circular arc, `L` = 1 µm, `R` = 2 µm, `κ` = 0.075 pN·µm².
Closed form `E = κL/2R²` = **0.009375 pN·µm**.

| `n` | **Cytosim** / continuum | Aleph `INTERIOR_ONLY` / continuum | `(n−2)/(n−1)` | **`END_CORRECTED` / Cytosim** |
|---:|---:|---:|---:|---:|
| **3** | **0.99723** | **0.49870** | **0.50000** | **1.000164** |
| 5 | 0.99970 | 0.74951 | 0.75000 | 0.999647 |
| 9 | 0.99877 | 0.87486 | 0.87500 | 1.001065 |
| 17 | 1.00256 | 0.93746 | 0.93750 | 0.997405 |
| 33 | 1.00447 | 0.96874 | 0.96875 | 0.995543 |

### What it says, and it is the first of its kind in this repository

**1. The defect `ALEPH-PORT-3639` found is real, and Cytosim does not have it.** Aleph's
`INTERIOR_ONLY` sum reproduces `(n−2)/(n−1)` to **four decimal places at every `n`** — it is exactly
half the continuum at the cortex's `n` = 3. Cytosim sits within **0.5% of the continuum at every
`n`**, including at `n` = 3 where Aleph's default is 2× out.

**2. The correction is right.** `END_CORRECTED` agrees with Cytosim to **0.016% at `n` = 3**, which
is the mesh the cortex actually runs on, and to better than 0.5% everywhere.

`-3639` was found because the *provider* had recorded the factor against Cytosim, and this lane
adopted the correction on an internal closed-form check. **Until now nothing outside this repository
had graded it.** Now something has, and it agrees.

### Where the agreement is weakest, stated rather than smoothed

Cytosim drifts **above** the continuum as `n` rises — 1.00256 at 17, 1.00447 at 33 — so
`END_CORRECTED/Cytosim` falls to 0.9955 there. That is Cytosim's own discretisation and this entry
does not explain it; `-3642` §3 says a disagreement is a finding and which side is wrong is a
separate measurement. It is under 0.5%, it is at meshes the cortex does not use, and it moves the
wrong way to be Aleph's.

**No Aleph law changed.** E-6 holds: `git diff` touches nothing under `aleph/`.

### The gates now

| # | outcome |
|---|---|
| **E-1** the harness degrades | **PASS** — 5 controls green with no binary |
| **E-2** the firewall holds | **PASS for this port**; `tests/firewall` red only for `aleph/scenarios/ecm_provider.py`, reported not patched |
| **E-3** the same configuration reaches both | **PASS** — the arc is gated to 1e-14 on radius, and both codes are given the same 17-digit vertex list |
| **E-4** the comparison runs and reports a number | **PASS, on the ENERGY not the end-node force.** The observable changed and the reason is measured, not preferred |
| **E-5** both Aleph quadratures compared | **PASS** — both columns above |
| **E-6** the disagreement is reported, not resolved | **PASS** — nothing under `aleph/` changed |

**E-4's threshold was "a table, not a verdict", and the table is above.** The verdict it happens to
support — that `-3639` was right — is the strongest evidence this project has for any of its laws,
because it is the only one that did not come from inside.
