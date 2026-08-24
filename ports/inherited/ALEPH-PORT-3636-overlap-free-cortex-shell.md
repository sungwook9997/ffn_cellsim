# ALEPH-PORT-3636 — the overlap-free cortical shell, and the wall that has a closed form

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3636` |
| Lane | `46143f30` Lane W2 (lead) |
| Status | `PROPOSED` |
| Written | `2026-08-06` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **RE-DERIVATION.** No file is transcribed and nothing is imported at runtime. The provider's *algorithm* is read, its geometric target is adopted, and the implementation is written against Aleph's own `Filament` / `node_start` layout. |
| Provider | `/Users/sw1/ffn_cellsim/ffn_sim/ff/cortex_assembly.py` (466 lines), READ ONLY, untouched |
| Aleph target | `aleph/vertical/cortex_surface_coupling.py` — one new placement mode, one new relaxation function, one new diagnostic |
| Exists because | Aleph's cortex cannot be built above ρ ≈ 8 filaments/µm² against a sourced 100, and the reason is a construction defect this lane measured today. |

---

## 0. The claim this port is answering

`docs/results/2026-08-06-cortex-n-dependence` established that Aleph's areal modulus `K_A` is
proportional to the filament count rather than converging — that the engine has measured *one
filament's stiffness times how many there are*. It could not go further because the builder breaks
down above ρ ≈ 8/µm².

That break has now been measured rather than assumed, and it is not what the earlier note supposed.

## 1. The wall, measured

`scripts/_ws_density_wall.py`, run on the workstation (64-thread CPU, `~/wp/bin/python`):

| ρ [/µm²] | filaments | nodes | min cross-filament separation [µm] | WCA σ [µm] | pairs < r_c | `K_A` | build [s] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 2,493 | 7,479 | **0.034551** | 0.00700 | **0** | 1.397e+04 | 2.4 |
| 10 | 3,117 | 9,351 | **0.001963** | 0.00700 | 30 | 1.757e+04 | 3.2 |
| 12 | 3,740 | 11,220 | **0.000836** | 0.00700 | 18 | 2.919e+06 | 4.1 |
| 16 | 4,986 | 14,958 | — | 0.00700 | 18 | 2.771e+04 | 5.6 |
| 24 | 7,480 | 22,440 | — | 0.00700 | 11 | 4.156e+04 | 8.2 |
| 32 | 9,973 | 29,919 | — | 0.00700 | 16 | 1.264e+06 | 11.6 |
| 48 | 14,959 | 44,877 | — | 0.00700 | 537 | 8.990e+08 | 19.0 |
| 64 | 19,946 | 59,838 | — | 0.00700 | 777 | 2.126e+09 | 27.7 |
| 100 | 31,165 | 93,495 | — | 0.00700 | 1,290 | 1.135e+15 | 49.1 |

(The separation column is exact brute force and was computed only where the pairwise matrix fits;
`—` is not measured.)

**Three things follow, and the first two contradict what this lane wrote yesterday.**

1. **The build does not fail.** ρ = 100 builds 31,165 filaments and 93,495 nodes in 49 s on one
   core. There is no ceiling in the constructor. The earlier note said "the golden-angle builder's
   interpenetration ceiling is at 8"; what is at 8 is the last density whose `K_A` is finite.
2. **`pairs` is a neighbour list, not an overlap count.** At ρ = 24 there are 11 pairs and
   `K_A`/filament is **5.5555** — indistinguishable from the clean 5.6045 at ρ = 8. Eleven pairs
   inside `r_c` cost nothing. What costs is *how deep* the worst one is.
3. **The wall has a closed form.** Filaments are tangent segments of length `2·half` on one sphere,
   with orientations drawn from a golden-angle sequence *deliberately decorrelated from position*
   (so the population is isotropic — the docstring says so). Nothing prevents two neighbours from
   crossing. Mean centre spacing is `1/√ρ`, so crossings must begin when

   ```
   2 · half_length_um  ≥  1/√ρ        ⟹        ρ_crit = 1 / (2·half)² = 11.1 /µm²
   ```

   at `half = 0.15 µm`. Measured onset: between 8 and 10. **The wall was predictable from the
   construction and no one predicted it.**

Once a pair is inside the core the WCA `r⁻¹²` does the rest: at ρ = 12 the closest pair sits at
0.836 nm against σ = 7 nm, and `(7/0.836)¹² = 1.2e11`. That single pair is the 2.9e6 modulus.

## 2. Radial thickness alone does not fix it — measured, not assumed

`build_radial_cortex_network`'s `shell_thickness_um` defaults to **0.0**: every node at exactly one
radius. A real cortex is ~0.2 µm thick, which is ~28 filament diameters of depth, so the obvious
hypothesis is that the wall is a dimensional error. `scripts/_ws_shell_thickness.py` tested it:

| ρ | T = 0 | T = 0.05 | T = 0.10 | T = 0.20 |
|---:|---:|---:|---:|---:|
| 8 | 0 pairs, `K_A` 1.397e4 | 0, 1.411e4 | 0, 1.425e4 | 0, 1.454e4 |
| 12 | 18, **2.919e6** | 11, 3.234e4 | 13, 2.195e4 | 10, 2.171e4 |
| 24 | 11, 4.156e4 | 12, 4.186e4 | 11, 4.229e4 | 11, 4.315e4 |
| 48 | 537, **8.990e8** | 573, 1.069e8 | 590, 8.345e8 | 602, 3.297e7 |
| 100 | 1,290, **1.135e15** | 1,315, 2.638e9 | 1,383, 1.927e8 | 1,544, 8.223e8 |

**Thickness is worth seven orders of magnitude at ρ = 100 and is still four orders short.** The pair
*count* rises slightly with thickness rather than falling. Dispersing radially breaks the exact
radial coincidence of crossing nodes — which is why the worst pair gets much less bad — but it does
not separate them to contact, because the depth assignment knows nothing about which filaments are
laterally adjacent.

**This is the finding that decides the port**: the missing half is not thickness, it is a
relaxation. Aleph has no build-time relaxation of any kind.

## 3. What the provider does, and what is adopted

`ffn_sim/ff/cortex_assembly.py` carries `overlap_free_shell(...)` and, in its own config,
`overlap_free_cortex: bool = True  # radial ~0.2µm thickness + WCA relaxation`. Two things are worth
recording:

* **The provider reached the same diagnosis and needed both halves too.** Their default is thickness
  *and* relaxation, not thickness. This lane arrived at the same split by measurement, which is the
  only reason it is stated here as a finding rather than as deference.
* **σ_EV = 0.007 µm on both sides.** Aleph's `f_actin_cortex_card` sets `steric_diameter_um = 2·r`
  with the F-actin radius; the provider hard-codes 7 nm. They agree, so the geometric target is
  the same object in both projects and the comparison in §5 is meaningful.

The algorithm adopted, stated as geometry rather than as code:

> Repeat until no cross-filament node pair is closer than `r_c = 2^(1/6)·σ`: find every such pair,
> and move each of its two nodes apart along their separation direction by half the shortfall to
> `1.01·r_c`. Same-filament neighbours are never pushed — they are the bonded arc.

**The target is param-free.** `r_c` is fixed by σ, and σ is sourced. There is no tuned constant, no
force gate, and no fitted relaxation rate — which is the property that makes this a construction
step rather than a knob. The `1.01` settles just past the cutoff instead of exactly on it.

The provider also offers a second mode (`radial_span`) that separates crossings out-of-plane and
spreads the displacement over a cosine bump along each filament, to avoid kinking a fine-meshed
filament and blowing up its bending force. **That mode is not ported here.** Aleph's cortical
filaments are three nodes long; there is no arc to kink and no span to spread over. Porting it would
be porting a solution to a problem this engine does not have. It is named here so that the omission
is a decision on the record rather than an oversight, and so the day Aleph's filaments get a fine
mesh, this entry says where to look.

## 4. What is NOT copied

- No line of `cortex_assembly.py` is transcribed. Aleph stores filaments as `Filament(filament_id,
  node_start, node_count, polarity)` against one `positions` array; the provider uses a
  `fiber_offsets` prefix-sum. The node-ownership map is derived from Aleph's own structure.
- The provider's `rng`-drawn per-fiber inward offset is **not** adopted. Aleph's `_shell_depths` is
  a golden-ratio low-discrepancy sequence and is deterministic; this project's builds are
  reproducible without a seed and that property is not traded away for a closer match.
- `count_interpenetrations` is re-derived as a diagnostic returning what Aleph's tests need (the
  minimum separation in units of σ), not the provider's dict.
- No import of `ffn_sim` is added. This is a re-derivation, unlike `ALEPH-PORT-3635`, because the
  algorithm is fifteen lines of geometry against a data layout Aleph does not share.

## 5. The gate — what would make this port a failure

Stated before the code, so it cannot be adjusted to what the code produces.

| # | Gate | Threshold |
|---|---|---|
| G-A | **Zero interpenetration at native density.** Cross-filament pairs closer than `r_c` at ρ = 100, T = 0.2. | exactly `0` |
| G-B | **`K_A` becomes finite and per-filament.** `K_A`/filament at ρ = 100 against the clean value 5.6045 measured at ρ = 8. | within 25% |
| G-C | **The relaxation does not move filaments off the shell.** Max node displacement. | `< σ = 0.007 µm` |
| G-D | **The old behaviour is retained as a control.** Default arguments reproduce today's numbers bit for bit. | `K_A` = 1081.819088 at ρ = 0.5, level 2 |
| G-E | **Determinism.** Two builds with identical arguments give identical positions. | bitwise |

**G-B is the one that matters and it is the one that can fail honestly.** If `K_A`/filament stays
near 5.6 once the overlaps are gone, then the cortex is *still* a sum over non-interacting filaments
even at native density, and the N-dependence finding stands unweakened — the engine will have proved
that its cortex is not a network at any density it can reach. If instead it drops, the filaments
have started to constrain each other and the network is real. **Both outcomes are results.** This
port is built to tell them apart, not to produce the second.

## 6. Cost, and why the GPU is the next entry and not this one

The relaxation is a KD-tree query plus a scatter-add per sweep. At ρ = 100 that is 93,495 nodes and
~60 sweeps, entirely on the CPU, once per build. It is not the hot loop, so **this port is not a GPU
port** and does not need one.

What it *unlocks* is the hot loop: a cortex that can be built at 93,495 nodes is a cortex whose
per-step steric and bending evaluation is 12× the largest Aleph has ever stepped. That is the
device-residency question, and it belongs to its own entry.

## 7. Reproducing §1 and §2

```bash
# on the workstation, results returned to the Mac; nothing is left on the shared host
MEM_GB=24 CPUS=8 wsrun.sh wall  cpu 00:40:00 scripts/_ws_density_wall.py
MEM_GB=32 CPUS=8 wsrun.sh thick cpu 01:00:00 scripts/_ws_shell_thickness.py
```
