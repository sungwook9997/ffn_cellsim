# ALEPH-PORT-3638 — the cortical microstructure: one resonance, three defects

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3638` |
| Lane | `46143f30` Lane W2 (lead) |
| Status | `PROPOSED` |
| Written | `2026-08-06` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **§A is an internal FIX**, not a port: nothing is read from the provider for it. **§B is a RE-DERIVATION** of the provider's arrangement *shape* against Aleph's own structures. §C and §D are internal. |
| Provider | `/Users/sw1/ffn_cellsim/ffn_sim/ff/cortex_assembly.py` — read for §B's shape only, READ ONLY, untouched |
| Aleph target | `aleph/vertical/cortex_surface_coupling.py`, `aleph/vertical/cortex_filaments.py` |
| Depends on | `ALEPH-PORT-3636` (the overlap-free shell), `ALEPH-PORT-3637` (the device relaxation that measures it) |
| Exists because | The PI, 2026-08-06 14:08 KST: the crosslink and the filament arrangement have to be right before the motor will work. They are not right, and the reason is one line of arithmetic. |

---

## 0. The finding, in one sentence

**Position, depth and orientation are all assigned by the same golden ratio, and the positions are a
Fibonacci lattice — on which spatial neighbours sit at Fibonacci index offsets, where
`frac(F_k·φ) → 0` by definition.** So the sequence that is a superb low-discrepancy sequence *in the
index* is catastrophically correlated *in space*, and every property it assigns is handed to
neighbours nearly unchanged.

```
position     fibonacci_sphere_directions(i)          spatial neighbours: index offset = F_k
depth        T · frac(0.5 + i·φ)                     cortex_surface_coupling.py:575
orientation  2π · frac(i·φ)                          cortex_surface_coupling.py:872
```

| F_k | 34 | 55 | 89 | 144 | 233 | 377 |
|---|---:|---:|---:|---:|---:|---:|
| `frac(F_k·φ)` | 0.0132 | 0.0081 | 0.0050 | 0.0031 | 0.0019 | 0.0012 |

Measured index offsets between the six nearest spatial neighbours at ρ = 100: **34, 55, 89, 144,
233, 377.** Exactly the Fibonacci numbers.

The builder's own docstring (`:868-871`) claims the opposite:

> *"Orientation is decorrelated from position by **the same irrational rotation** that placed the
> centres."*

**"The same irrational rotation" is the cause, not the cure.**

## 1. What it costs, measured

| | neighbours | random pairs | ratio |
|---|---:|---:|---:|
| **depth** difference (ρ = 100) | **1.41 nm** | 54.30 nm | **38.6×** more alike |
| **orientation** difference | **4.63°** | 90.75° | **19.6×** more aligned |

### 1a. The shell has no local thickness

`cortex_shell_thickness_um = 0.2` was added to give the cortex radial neighbours — `cortex.py:22`
says *"A shell with no thickness is not a shell… it has no radial neighbours, so nothing in it
resists bending or shear **as a shell**."*

| | global radius span | **local span (6 nearest neighbours)** | actin diameter |
|---|---:|---:|---:|
| ρ = 8 | 199.8 nm ✓ | **4.96 nm** | 7.0 nm |
| ρ = 100 | 203.0 nm ✓ | **1.41 nm** | 7.0 nm |

**Globally 200 nm, locally 1.4–5 nm.** It is a corrugated sheet, not a filled layer, and the
parameter added to create radial neighbours does not create them. The global span is right, so
nothing checks it.

### 1b. The arrangement is not isotropic, and its director is a construction axis

| | measured | isotropic |
|---|---:|---:|
| global nematic order `S` (largest eigenvalue of `⟨(3/2)t⊗t − I/2⟩`) | **0.2167** | 0 |
| director | **[0, −1, 0]** | none |
| `⟨\|t·ê_θ\|⟩` | 0.7676 | 0.6366 |
| `⟨\|t·ê_φ\|⟩` | 0.4718 | 0.6366 |

Not a cap artefact: `S` = 0.2377 on `|z| > 0.9` and 0.1871 on `|z| < 0.5`.

**This is what the non-affinity field was showing.** `docs/results/2026-08-06-native-density-on-the-4090`
§5 recorded a correlation of −0.479 between the non-affinity and latitude and called it "the
lattice" without saying which property of it. It is these two, and they are one.

### 1c. Every published cortex number carries both

`S = 0.2167` and a 1.4 nm local thickness are in every modulus this lane measured. The *relative*
statements (affine ≫ relaxed; exponent 1 → 0.35) are measured on the same lattice on both sides and
survive. **The absolute moduli do not**, and §5 already says so.

## 2. What density does, and does not, rescue — measured before designing

`scripts/_ws_crosslink_study.py` on the 4090 (job 26), crosslinker reach swept as a **material**
parameter, every row relaxed to a 1e-3 pN residual with the translation mode projected out:

| ρ | reach | X/F | z | link [nm] | **components** | `K_A` relaxed |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 35 (α-actinin) | 0.02 | 0.05 | 31 | **2,436** | 0.019 |
| 8 | 160 (filamin) | 0.70 | 1.40 | 106 | **746** | 0.091 |
| 8 | 778 (today) | 10.09 | 20.19 | **462** | 1 | ~412 |
| 24 | 60 | 0.65 | 1.30 | 31 | **2,645** | 0.057 |
| 24 | 100 | 1.32 | 2.64 | 57 | **1** | 0.978 |
| **24** | **160 (filamin)** | **2.78** | **5.56** | **97** | **1** | **9.94** |
| 24 | 778 (today) | 13.01 | 26.02 | 258 | 1 | 615 |

**Density rescues it.** At ρ = 8 filamin leaves 746 fragments; at ρ = 24 it percolates, `z = 5.56`
clears the isostatic threshold of ~5, and **`X/F = 2.78` lands inside the physiological 3–6.** `K_A`
climbs 0.091 → 9.94 for a 3× density — an exponent near 4.3, the signature of a percolation
transition rather than of a spring count.

**So the crosslinker does not need to be long. The population needs to be dense.** That is the PI's
claim — fewer filaments is not the problem, the joining is — with the correction that at ρ = 8 there
is no joining to be had at any physical reach.

Rows for ρ = 48 and ρ = 100 were still running when this was written and are **not** used below.

## 3. §A — break the resonance. The change

Replace the golden-ratio sequences that assign **depth** and **orientation** with a sequence that
does not resonate with the Fibonacci lattice. The **positions** stay a Fibonacci lattice: it is the
right way to put near-uniform points on a sphere and nothing here is wrong with it.

**The replacement: the van der Corput radical inverse in base 2.** `vdc(i)` reverses the bits of `i`
about the binary point. Its correlation structure is governed by powers of two, and the Fibonacci
numbers are not powers of two — `F_k` and `2^m` coincide only at 1, 2 and 8 — so the offsets that
make neighbours neighbours carry no information about the assigned value.

```python
depth_i        = T · vdc2(i)                 replaces  T · frac(0.5 + i·φ)
orientation_i  = 2π · vdc2(i ⊕ MASK)         replaces  2π · frac(i·φ)
```

The two need **different** sequences from each other as well, or depth and orientation become
functions of one another; the XOR mask is what separates them, and it is a bit pattern rather than
a second irrational so the two are provably independent in the bits rather than approximately so.

**This is opt-in.** `_shell_depths` and the orientation assignment take a mode argument whose default
reproduces today's sequence bit-for-bit, because every committed cortex number was measured on it.

## 4. §B, §C, §D — declared here, specified after §A is measured

Written now so the programme is visible, and deliberately **not** specified in detail before §A's
result is in — §A alone may move `S` and the local thickness enough to change what the others need.

* **§B — the arrangement becomes a declared material.** `orientation ∈ {isotropic, circumferential,
  meridional}` and a nematic order `S ∈ [0,1]`, with `t = normalize(S·director + (1−S)·random)`, plus
  a `nematic_order()` diagnostic so "isotropic" is a measurement and not a docstring. This is the
  provider's shape (`ff/cortex_assembly.py:104-152`), re-derived. **No physiological `S` is adopted
  until it is sourced**; the default stays whatever §A produces.
* **§C — the filament node count stops being 3.** `Filament(..., 3, ...)` is hard-coded at
  `cortex_surface_coupling.py:855` and `:884`; a filament is two 150 nm segments with **one** bending
  triple. The provider uses a 75 nm mesh. Every downstream dependence on `node_count == 3` must be
  enumerated before this moves.
* **§D — the crosslinker becomes a molecule.** A card with a contour length (fascin ~5 nm,
  α-actinin ~35 nm, filamin ~160 nm) and a stiffness, both sourced; a crosslink exists **iff** the
  closest approach is within that length, and does not exist otherwise. Today
  `crosslink_capture_from_spacing_um` returns `2.2 × spacing` — a geometric consequence of the
  placement rule, which is the same class of error `filament_areal_density_per_um2` was introduced
  to fix for the filament *count*.

## 5. The gates for §A — written before the code

| # | Gate | Threshold |
|---|---|---|
| **A-1** | **Depth decorrelates.** Mean radius difference to the 6 nearest spatial neighbours, ρ = 100, T = 0.2. | `≥ 40 nm` (random pairs give 54.3; today gives 1.41) |
| **A-2** | **Orientation decorrelates.** Mean angle difference to the 6 nearest spatial neighbours. | `≥ 70°` (uniform gives 90.75; today gives 4.63) |
| **A-3** | **The arrangement becomes isotropic.** Global nematic order `S`. | `≤ 0.05` (today 0.2167) |
| **A-4** | **The shell acquires local thickness.** Local radius span per solid-angle cell, against the 203 nm global span. | `≥ 100 nm` (today 31 nm) |
| **A-5** | **The retained control is untouched.** Default arguments reproduce `K_A = 1081.819088` at ρ = 0.5, level 2. | bitwise |
| **A-6** | **Determinism.** Two builds with identical arguments give identical positions. | bitwise |
| **A-7** | **The non-affinity field stops tracking latitude.** Correlation between `\|u\|` and `\|cos(latitude)\|` at ρ = 100. | `\|corr\| ≤ 0.15` (today −0.479) |

**A-3 and A-7 can fail honestly, and that is why they are here.**

If `S` falls to ~0 and the latitude correlation **persists**, then the anisotropy is not the
orientation sequence and something else in the construction carries a polar axis — a result, and one
that would invalidate this ledger's §0 rather than confirming it.

If `S` falls and the latitude correlation goes with it, then §1c is confirmed and **every absolute
cortex modulus in the repository has to be re-measured**, which is hours on the 4090 and is planned
rather than discovered.

## 6. What is NOT claimed

- **Not that §A fixes the modulus.** It removes an artefact. Whether `K_A` then lands in the
  physiological band is §2's question and depends on density and on §D.
- **Not that the Fibonacci lattice is wrong for positions.** It is the right construction; the error
  is reusing its own irrational for the properties assigned on it.
- **Not that ρ = 24 with filamin is the answer.** `K_A` relaxed there is 9.94 pN/µm against a band of
  100–4,000, i.e. still an order low. ρ = 48 and ρ = 100 were unmeasured when this was written.
- **Not that van der Corput is the only fix.** Any sequence whose correlation structure is
  incommensurate with the Fibonacci offsets would do; base 2 is chosen because the non-resonance is
  provable from the bits rather than argued from a plot.
- **No physiological nematic order is adopted anywhere in this entry.** §B names the mechanism and
  explicitly defers the value to a source.

## 7. Reproducing §0 and §1

```bash
cd /Users/sw1/Project_Aleph
OPENBLAS_NUM_THREADS=1 /Users/sw1/miniconda3/envs/aleph/bin/python \
    scripts/_ws_lattice_resonance.py            # written with the §A code, not before it
MEM_GB=48 CPUS=16 scripts/wsrun.sh xlink2 4090-1 04:00:00 \
    scripts/_ws_crosslink_study.py cuda 8,24,48,100 0.035,0.06,0.10,0.16,0.25,0.40,auto 400000
```

---

## Amendment, 2026-08-06, after §A ran: A-2's threshold was unreachable, and the frame carries residue

Recorded as an amendment rather than by editing §5. **This is the fourth time this lane has written a
threshold against a number computed by a different estimator** (`ALEPH-PORT-3637` A1 records the
first three), and the pattern is now the finding: a gate written from a remembered number rather
than from a run of the estimator it will be tested with.

### Measured, ρ = 8, level 2, `scripts/_ws_lattice_resonance.py`

| quantity | `GOLDEN_RATIO` | `VAN_DER_CORPUT` | independent draws |
|---|---:|---:|---:|
| neighbour depth difference [nm] | 4.9586 | **82.8273** | 65.70 |
| neighbour angle difference [deg] | 6.7107 | **48.9038** | **56.94** |
| global nematic order `S` | 0.2167 | **0.0022** | 0 |
| local radius span [nm] | 32.2615 | **186.4516** | 199.86 *(global span)* |

| gate | outcome |
|---|---|
| A-1 depth ≥ 40 nm | **PASS** — 82.83 |
| A-2 angle ≥ 70° | **FAILED** — 48.90. See below. |
| A-3 `S` ≤ 0.05 | **PASS** — 0.0022, from 0.2167 |
| A-4 local span ≥ 100 nm | **PASS** — 186.45 of a 199.86 global, i.e. 93% of nominal, from 16% |
| A-5 retained control | **PASS** — `K_A` 1081.819088, bitwise |
| A-6 determinism | **PASS** — bitwise |
| A-7 latitude correlation | not run here; needs a relaxation |

### A-2: the threshold could not be met by anything

70° was taken from a measurement of the **scalar angle parameter**, `|((θᵢ − θⱼ + π) mod 2π) − π|`,
whose uniform mean is 90° and which measured 90.75° for random pairs. The gate script measures a
different quantity: the angle **between tangent vectors**, taken through `|cos|` so that `t` and
`−t` are the same director, which is bounded by 90° and — because tangents on a sphere live in
different tangent planes — gives **56.94° for independent draws**. Nothing can reach 70° on this
estimator.

**A-2 is recorded FAILED as written.** The corrected criterion, stated separately rather than
substituted:

> **A-2′.** The neighbour angle difference is at least **0.80×** what independent draws give on the
> **same** estimator, measured in the same run. Here: 48.90 / 56.94 = **0.859**, which passes.

### RETRACTED, an hour later: there is no residue and `_tangent_frame` is not a defect

The paragraph that stood here claimed a "second defect" — that 48.90° against 56.94° left a 14%
neighbour-alignment residue, and that `_tangent_frame`'s fixed helper axis was the cause. **Both
halves are wrong, and the error is the same one A-2 already committed: the wrong baseline.**

Neighbours on a sphere **share a tangent plane**. Two independent in-plane angles differ by 45° on
average. The 56.94° figure is the mean over *cross-sphere random pairs*, whose tangent planes are
different, and that difference inflates the angle — it is not what neighbours would give under
independence.

Measured with the control built correctly — each filament's angle recovered in **its own** frame,
the angles permuted, the tangents rebuilt in the original frames:

| angle sequence | neighbour angle | ÷ control |
|---|---:|---:|
| `GOLDEN_RATIO` | 6.71° | **0.149×** |
| `VAN_DER_CORPUT` | **48.90°** | **1.089×** |
| permuted control (independent angles, same frames) | **44.89°** | 1.000 |
| 2-D uniform, theory | 45.00° | — |

`VAN_DER_CORPUT` is **above** independence, not below it. That is what a low-discrepancy sequence
does — it stratifies — and it means the orientation correlation is gone rather than reduced.

And `_tangent_frame` cannot cause the effect that was attributed to it: `cos(θ)·e₁ + sin(θ)·e₂` with
θ uniform gives a uniform direction in the plane **whatever orthonormal basis (e₁, e₂) is**. A smooth
frame is not a correlated one.

**Three baselines were used for one gate before the right one was found** — 90.75° (signed angle
parameter), 56.94° (cross-sphere pairs), 57.45° (permuted tangent *vectors*, which hands a
north-pole filament a south-pole tangent), and finally 44.89° (permuted *angles* in their own
frames). The gate script now measures its control **in the same run** rather than carrying a
remembered number, which is the actual fix and is worth more than the threshold it repairs.

### A-2 status

**FAILED as written** (≥ 70° is unreachable: the control itself is 44.89°). Superseded by:

> **A-2′.** The neighbour angle difference is at least **0.90×** the permuted-angle control measured
> in the same run. Measured: 48.90 / 44.89 = **1.089**, PASS. `GOLDEN_RATIO` gives 0.149.
