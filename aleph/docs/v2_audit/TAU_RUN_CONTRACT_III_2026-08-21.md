# Contract III — the instrument changes, and it changes because it was measured bit-identical

**2026-08-21 07:35 KST, before any step of run III executes.** Contracts
[I](TAU_RUN_CONTRACT_2026-08-21.md) and [II](TAU_RUN_CONTRACT_II_2026-08-21.md) stand as the record.
**§4 — the read order, the trailing windows, the resolution criterion — is carried here unchanged for
the third time.** It has not been edited once, under any of the three sizings, and that is the point
of keeping it in its own section.

## 1. What changed, and the evidence for it

Contract II kept the host estimator because an A/B rejected the device path at **1.54e-01**. The cause
was mine: `gamma_from_resting_readback`'s `centre` defaults to the actin centroid and I had passed the
**origin**. Measured, the centroid is **2.05e-02 µm** from origin — twenty nanometres — and a centre
offset does not perturb the sum, it **flips which elements cross a plane**. Discrete, not round-off.

Reducing the centroid on the device (`wp.utils.array_sum` over 4.2 M vec3d, **0.474 ms**) rather than
reading 101 MB back, and passing it:

```
device with ORIGIN    vs host    worst rel  1.539205e-01
device with CENTROID  vs host    worst rel  5.129989e-16      ← 12 steps, seed 7, per step
```

**Bit-identical to machine precision, step by step**, against a host-vs-host repeat of **2.17e-16**
that establishes the engine is deterministic for a fixed seed. And the step cost:

| | s/step | |
|---|---:|---|
| host estimator | 6.4204 | contract II's instrument |
| **device kernel** | **0.2805** | **22.9×** |

⚠ **This is an instrument change and it is why this is a new contract rather than an edit.** The
quantity is the same to 5.13e-16; the code computing it is not.

## 2. What the change buys, which is why it is worth making

| | host | device |
|---|---|---|
| 1 seed × 5,250 steps (contract II) | 9.4 h | 0.41 h |
| **3 seeds × 25,000 steps** | 133 h — impossible | **5.84 h** |

⚠ **Contract II could answer τ on one seed, at a 50-τ bar of T\* ≤ 104, or answer nothing.**
**Run III can answer τ over a window 4.8× longer AND supply the seed scatter D-2 asks for** — which
`(e) 1` calls *"the single most likely way to write the amendment and still be measuring nothing"* if
it is taken from within-run variance instead.

## 3. Run III

**Three seeds — 1, 2, 3 — each 30,000 steps**, run sequentially so an interrupted budget leaves whole
seeds rather than three partial ones. γ every step, device estimator, heartbeat every 30 s. Same
thermostat configuration throughout: cortex bound, `k_axial` 88000 pN/µm, κ 0.07 pN·µm², mobility
1e-06 µm/pN·s ⚠ UNSOURCED, T = 310 K, **full cell, no `--core-only`**.

**7.36 h against ~9.5 h of grant at launch.**

### 3.1 ⚠ The τ ceiling this budget buys, computed BEFORE the run

§4.3 needs a discarded transient plus an analysed window. At the declared bars that is
`7τ + 50τ ≤ N`:

| steps / seed | total | **covers τ up to** |
|---:|---:|---:|
| 25,000 | 6.13 h | 439 steps |
| **30,000** | **7.36 h** | **526 steps** |
| 35,000 | 8.58 h | 614 steps |

⚠ **Run I's last measured τ was 188 steps and it was a LOWER BOUND on an unresolved estimate**, from a
ladder (1.139 → 2.575 → 9.384 s) that grew while `n_eff` stayed pinned at 4.16 — which means τ was
tracking the budget, not settling at 188. **So τ > 526 is a live possibility and it would end all three
seeds REFUSED.**

**30,000 rather than 25,000 buys 87 more steps of τ headroom for 1.2 h**, and the insurance that
matters here is **steps, not seeds**. ⚠ **Seeds are not tradeable**: three is D-2's minimum and
`observe_gamma_seed_scatter` refuses fewer (`min_replicates=3`), because a scatter over two replicates
has one degree of freedom and **produces a number without being a measurement.**

⚠ **And if τ exceeds 526, that is the result, not a failure to be re-budgeted.** *"Stationary γ is not
reachable in this regime within this budget"* is a decision-grade statement, and it is precisely what
run I could not deliver.

⚠ **The 50-τ bar now passes if T\* ≤ 500 steps** — against contract II's 104 and contract I's 195, and
a last measured lower bound of 188. **This is the first sizing under which the bar is not on a knife
edge**, and that is a consequence of the instrument, not of a threshold being moved.

## 4. Carried unchanged, third time

**4.1** Read order: `min_tau_reliability` → `max_tau_relative_error` → `required_window_samples`.
**4.2** Windows are TRAILING — last 1/8, 1/4, 1/2, whole.
**4.3** RESOLVED when τ(last half) and τ(whole) agree within their combined Madras–Sokal 1σ.

## 5. ⚠ The control this contract owes, declared before the run rather than after

**The bit-identity is measured over 12 steps. Run III is 25,000.** Nothing observed says it degrades —
the engine is deterministic and the per-step agreement is at machine precision with no trend — but
**12 is not 25,000**, and a per-step identity that held for twelve steps is not a proof it holds for
twenty-five thousand.

> **A host-γ control of seed 1, 250 steps, is run AFTER the seeds and compared against run III's first
> 250 values. It costs 27 minutes of the remaining margin. If it does not reproduce them to ~1e-15,
> run III's series are withdrawn, not reinterpreted.**

⚠ Declared here so it cannot be skipped once the numbers look reasonable, and so its failure condition
is written before there is anything to be reluctant about.

## 5.1 ⚠ THE CELL RUN III STEPS HAS NO MOTORS, and nothing said so until a sequence was exported

**`world_phase4_native.py` never builds NMII.** It stands `build_all` (cortex, membrane, envelope) plus
`build_remaining_populations` (the nine), which is **eleven populations**. `world_phase1_native.py`
builds a twelfth when `--nmii-heads-per-side` is passed. **The driver that steps the cell has no such
flag.**

So: **run III's τ is the correlation time of a γ series from a cell with no myosin minifilaments.**

⚠ **That is consistent with the resting path and it is not obviously wrong** — the resting γ has two
families, actin and crosslink, and no motor term. But **"FULL NATIVE" appears throughout this contract
and in `STATE.md`, and it reads as "everything".** The counts are distinguishable if a reader compares
`4,558,554 / 11` against PHASE 1's `4,591,262 / 12`; **nothing anywhere says which population the
difference is, or that the missing one is the contractile machinery.**

⚠ Found by exporting a position sequence for the viewer and reading its population list: eleven, with
`nmii` absent. **Not by reading either driver**, both of which state their builds plainly — the same
route as the lamina, where a consumer failing was the first thing to ask the question.

**Nothing about run III changes.** The scope was always this; only its statement was missing. It is
recorded here so that no reading of run III's τ can quietly become a statement about a contractile
cortex.

## 6. What may not be done with run III

⚠ **No magnitude.** Mobility is UNSOURCED, so τ in seconds is not a physical time; the criterion uses τ
in **steps**. ⚠ **No γ value is quotable** — `STATE.md` (c) 3 and (c) 17 stand and the coefficients are
declared test points. ⚠ **No stationarity verdict except by §4.3.** ⚠ **A truncated run is analysed as a
shorter run**, never padded, never restarted-and-concatenated — and runs I and II's frames are
discarded rather than joined to III's.
