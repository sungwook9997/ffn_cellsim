# ALEPH-PORT-3637 — the cortex relaxation, made device-resident

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3637` |
| Lane | `46143f30` Lane W2 (lead) |
| Status | `PROPOSED` |
| Written | `2026-08-06` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **RE-DERIVATION onto kernels.** Aleph's own three NumPy laws, transcribed to warp kernels and driven from a resident state. Nothing is read from the provider for the laws themselves. |
| Provider | the *pattern* only: `ffn_sim/ff/network_warp.py`, `ff/relax.py`, `ff/implicit_ff.py` — read for the residency discipline, not for constitutive content |
| Aleph target | `aleph/runtime/cortex_relax_kernels.py` (new), `aleph/runtime/resident_cortex.py` (new) |
| Depends on | `ALEPH-PORT-3609` (`residency.py`) for the invariant this reuses; `ALEPH-PORT-3636` for the shell that makes it necessary |
| Exists because | The measurement that decides whether this engine's cortex is a network cannot be finished on a CPU. |

---

## 0. The case for the device, measured rather than preferred

This lane spent yesterday arguing the opposite and was wrong about the scale, so the case is stated
as numbers.

`docs/results/2026-08-06-overlap-free-native-cortex` §4: the cortex's real modulus is the virial
**after** the strained shell relaxes, and it is 34–66× below the affine bound the project has
reported all along. Measured on the workstation's CPU:

| ρ | nodes | 4,000 descent steps | peak residual force |
|---:|---:|---:|---:|
| 8 | 7,479 | 41 s | 49.4 pN |
| 24 | 22,440 | 154 s | 57.7 pN |
| 48 | 44,877 | 361 s | 40.0 pN |
| **100 (sourced)** | **93,495** | **~750 s, extrapolated** | — |

**None of those converged.** 4,000 steps leaves 40–58 pN of residual tangential force, so every
relaxed modulus reported so far is an upper bound. Driving the residual to a scale-relative 1e-6
needs one to two orders of magnitude more steps, which at ρ = 100 is **hours per strain amplitude**,
and a modulus needs the amplitude swept.

That is the first thing this engine has needed that a laptop cannot finish, and it arrived by
measurement two hours after the overlap-free shell made ρ = 100 buildable at all. The engine cost
audit's §5.7 condition — *cortex ≥ 120,000 nodes* — is 93,495 nodes away from met on the count and
already met on the wall clock, because it counted one force evaluation and this is `O(10⁵)` of them.

## 1. What goes on the device

Every array `internal_forces()` touches is already flat, which is why this is a port and not a
redesign:

| term | arrays | shape |
|---|---|---|
| axial | `segments`, `segment_rest_um`, modulus | `(S,2) i64`, `(S,) f64`, `(S,) f64` |
| bending | `triples`, `voronoi_length_um`, `kappa` | `(T,3) i64`, `(T,) f64`, `(T,) f64` |
| crosslink | `ia, wa, ib, wb, k, rest` from `crosslink_arrays()` | `(X,2) i64`, `(X,) f64` ×4 |
| steric | `steric_pairs` | `(P,2) i64` |

**The steric term is not ported, and that is a decision with a measurement behind it.**
`ALEPH-PORT-3636`'s resolution settles every cross-filament pair just *outside* `r_c`, so at the
relaxed shell the steric force is **exactly zero** — measured, 0 pairs at every density from 8 to
100. Porting a term that evaluates to zero would buy nothing and would hide the fact that it is
zero. The driver therefore **refuses to start** if the shell has any steric pair, rather than
silently omitting a force that has become non-zero. That refusal is gate G-4 below.

## 2. The invariant, reused rather than restated

`ALEPH-PORT-3609` established, for `ResidentMembraneState`:

> Between `push_positions` and `pull_positions`, no array whose element count scales with the mesh
> crosses the host boundary.

This port makes the same claim for `ResidentCortexRelaxation` between `push` and `pull`. The descent
reads back **one float64 per iteration** — the peak tangential force, which sets the step — and
nothing else. At ρ = 100 a mesh-sized crossing would be 280,485 elements; the budget is 1.

**The adaptive step is why there is a readback at all**, and the alternative was considered: a fixed
step needs no scalar and no sync, but it is either unstable at the start or uselessly slow at the
end, and this measurement is about the residual reaching zero. One float64 per iteration is inside
the invariant as `ALEPH-PORT-3609` wrote it, and it is named here so nobody has to infer it.

## 3. The arithmetic decision

**Float64 end to end**, unlike `ALEPH-PORT-3609`'s integrator, which takes forces in the ratified
float32 compute channel. The reason is specific to what is being measured: the descent runs `O(10⁵)`
steps and the quantity of interest is a *difference of virials* at a strain of 1e-3. In float32 the
virial's own round-off is `~1e-7` relative, which is 1e-4 of the signal — and the residual has to
fall six orders below the initial force for the modulus to be converged. This is a measurement whose
answer lives in the bits float32 does not have.

That is a departure from the ratified channel and it is stated here, not buried: **this port is a
measurement instrument, not a step of the production pipeline**, and nothing in `aleph/vertical/`
calls it.

## 4. The gates — written before the code

| # | Gate | Threshold |
|---|---|---|
| G-1 | **Parity with the NumPy laws.** Assembled force, resident vs `internal_forces()`, on the same shell. | `≤ 4` ULP per component |
| G-2 | **Energy parity.** Assembled energy, same comparison. | `≤ 4` ULP |
| G-3 | **Residency.** Mesh-sized host crossings between `push` and `pull` over 100 iterations. | exactly `0` |
| G-4 | **The steric refusal fires.** A shell with a non-zero steric pair list is rejected, not silently relaxed without it. | typed refusal |
| G-5 | **The relaxed modulus agrees with the CPU** at ρ = 8, both stopped at 4,000 steps. | within 1% |
| G-6 | **It converges where the CPU could not.** Peak residual at ρ = 100. | `< 1e-3` pN |

**G-1 at 4 ULP and not at 0.** The scatter is an atomic add and CUDA does not order atomics, so a
node receiving contributions from several segments accumulates them in a nondeterministic order.
`ALEPH-PORT-3601` measured that cost and it is the reason a bitwise gate would be a gate on the
scheduler rather than on the physics. The integrator in `ALEPH-PORT-3609` *is* graded at 0.0 ULP,
because it has no scatter.

**G-6 can fail honestly.** If the residual stalls above 1e-3 the descent is the wrong solver for
this problem — steepest descent on a near-singular Hessian is famously slow — and the finding would
be that the cortex relaxation needs a conjugate-gradient or an implicit step, not a faster device.
That is a result and would be reported as one.

## 5. What is NOT claimed

- **No performance claim is made here**, and none of the gates is a speed. Wall-clock will be
  reported as a measurement with its card and its allocation, never as a ratio without one.
- **Not a production path.** `aleph/vertical/` does not call this and must not; the cortex's
  registered step is unchanged.
- **Not a general resident cortex.** It relaxes tangentially with the radius pinned. Any other
  boundary condition is a different instrument.

## 6. The allocation this runs under

RTX 4090 on `sungwook@100.110.26.26`, through `gpu-submit`, under the authorization cited in
`docs/design/GPU_POLICY.md` with its transcript coordinates. Nothing is stored on the host: code up
per job, results back by rsync, scratch removed.

---

## Amendment, 2026-08-06, after the first run: two gates were miswritten and one law was

Recorded as an amendment rather than by editing §4, because a gate that is quietly retuned to what
the code produced is not a gate. Both changes below are stated with the measurement that forced
them, and the original thresholds stay visible above.

### A1. G-1 and G-2 at 4 ULP were unachievable by construction

Measured, ρ = 2, level 2, shell strained by 1e-3 (`scripts/_ws_resident_cortex.py cpu 2.0`):

```
  FAIL  G-1  force parity with the NumPy laws (<= 4 ULP)     64 ULP over 5,607 components
  FAIL  G-2  energy parity (<= 4 ULP)                        28 ULP,  2.29454946 pN.um
```

The ledger's stated reason for 4 rather than 0 was *"the scatter is an atomic add and CUDA does not
order atomics."* **That reason is wrong here, and the test that would have caught it was run second
instead of first.** Summing the same per-element energies in five different orders — index, reverse,
sorted ascending, sorted descending, and exactly-rounded `math.fsum` — spans **2 ULP**. The device
sits **27 ULP** from the exactly-rounded sum. Reordering cannot produce that, so associativity is
not the cause.

What is: **writing the same formula two ways.** With every other operation identical,

| the squared length, written as | axial energy [pN·µm] | gap |
|---|---|---:|
| `np.linalg.norm(d, axis=1)` | 2.1359999999994437 | — |
| `np.sqrt(dx*dx + dy*dy + dz*dz)` | 2.1359999999994437 | **0 ULP** |
| `np.sqrt(np.einsum("ij,ij->i", d, d))` | 2.1359999999994326 | **25 ULP** |

NumPy disagrees with *itself* by 25 ULP depending on how the dot product is accumulated, which is
the same order as the device's 27. `wp.length(v)` is `sqrt(dot(v, v))` with warp's own accumulation.

So G-1 and G-2 as written required the transcription to reproduce not the formula but **NumPy's
specific dot-product accumulation order** — a gate on BLAS, not on physics. They are recorded as
**FAILED as written**, and the corrected criterion is:

> **G-1′ / G-2′.** Assembled force and energy agree with the NumPy laws to a **relative** tolerance
> of `1e-13`, which is two orders above the ~25-ULP floor that writing the same expression
> differently already costs, and eleven orders below anything the physics is sensitive to.

### A2. The first descent could not converge, and G-6 caught it

G-6 failed on the first run: 400,000 iterations at ρ = 2, residual flat at **27.28 pN**.

The ledger anticipated a stall and named the wrong reason — *"steepest descent on a near-singular
Hessian is famously slow"*. The real cause was in the step rule. It read `step = c / peak`, which
makes the **peak node's displacement** the constant `c`: the descent never slows as it approaches
the minimum and instead oscillates around it at fixed amplitude. It was not a slow solver, it was a
solver that structurally cannot converge.

**This is a defect in the instrument and it is also a correction owed on already-published numbers.**
`docs/results/2026-08-06-overlap-free-native-cortex` §4 used the same rule on the CPU and reported
relaxed moduli with a 40–58 pN residual, described as "not converged". They were not converging.

Replaced by a mobility with energy backtracking — `x + µF`, µ halved when a step raises the energy
and grown 10% when it lowers it, with the rollback a **device-to-device copy** so the residency
invariant is untouched. Measured at ρ = 2: **76,926 iterations to 9.99e-4 pN in 18.3 s**, converged.

The relaxed modulus barely moved (327.5 pN/µm against the oscillation's plateau), so the published
*conclusion* survives; the *method* did not, and the numbers are re-measured rather than kept.

### A3. Gate status after the amendment

| # | as written | outcome |
|---|---|---|
| G-1 | ≤ 4 ULP | **FAILED** — 64 ULP. Superseded by G-1′ (relative 1e-13). |
| G-2 | ≤ 4 ULP | **FAILED** — 28 ULP. Superseded by G-2′. |
| G-3 | 0 mesh-sized crossings | **PASS** — 0 bulk, 144,557 scalar readbacks over 76,926 iterations. |
| G-4 | typed refusal | **PASS** — a shell with 11 steric pairs is refused. |
| G-5 | within 1% of CPU | see the results directory |
| G-6 | residual < 1e-3 pN | **FAILED, then PASS after A2** — 9.99e-4 pN. |
