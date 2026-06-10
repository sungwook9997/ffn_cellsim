# M-band targeting (loop24d) — A/B verification

**Date:** 2026-06-10  **Branch:** h7/full-cell-integration  **Device:** CPU (M1 Max), preliminary
**Verdict:** ❌ NO IMPROVEMENT → **preserved but not adopted** (default `mband=off`)

## What was tested

loop24d (uncommitted WIP) added *sarcomeric M-band targeting* to the cortex/SF
minifilament placer (`generate_cortex_myosin_layout`). Hypothesis: placing
minifilaments at the sarcomere **M-bands** (antiparallel-overlap zones) instead of at
random actin beads fixes "single-sided pulls" (a minifilament whose one side finds no
antiparallel partner contributes random-sign force) and should raise the **coherent
contractile traction** the ventral SF delivers to its FA anchors.

A/B harness: added `--mband-mode {on,off}` to `h7_ventral_sf_traction.py`
(threaded as `mband=` into `build_sf_sim`; `off` = pre-loop24d random greedy placement).
Each run is the existing same-seed paired force-ON−OFF coherent-traction differential
(thermal noise cancels exactly; +ve = contractile). M-band code only fires for layouts
that expose `m_band_x`, i.e. `--sarcomeric`.

**Config (matches the reference `h7_sf_2c_sarcomeric_gpu` geometry):**
`--sarcomeric --n-filaments 6 --n-sarcomeres 3 --n-motors 16
 --equilibrate 120000 --contract 120000 --n-samples 20`, seeds 1/2/3, modes on/off.
(Reference single-seed GPU run was 150k/150k/30samp → +131 ± 8 pN, CONTRACTILE.)

## Result

| seed | ON (M-band) [pN] | OFF (random) [pN] | paired Δ (on−off) |
|---|---|---|---|
| 1 | −1.94 ± 25.51 (eng 43) | **+101.66 ± 6.79** (eng 41) | −103.60 |
| 2 | +88.98 ± 7.59 (eng 40) | −285.20 ± 33.12 (eng 43) | +374.17 |
| 3 | −416.57 ± 19.38 (eng 49) | −333.12 ± 18.51 (eng 48) | −83.45 |
| **ensemble (n=3)** | **−109.84 ± 155.59** | **−172.22 ± 137.64** | **+62.38 ± 156.0 → NS** |

Figure: `mband_ab.png`. Summary: `AB_SUMMARY.json`. Raw: `ab_s{1,2,3}_{on,off}.json`.

## Interpretation

1. **No detectable improvement.** The paired differential (+62 ± 156 pN, n=3) is far
   below 2·SEM — M-band targeting does **not** measurably raise the coherent traction.
   On seed 1 it *destroys* a clean random-placement contractile signal
   (random +102 → M-band ~0, SEM 3.8× larger). → **not adopted**; code preserved with
   `mband` defaulting to **off** (the validated random path), opt-in via `--mband-mode on`.

2. **Bigger flag (separate from M-band): the per-SF coherent traction is seed-unstable.**
   Both arms swing in *sign* across realizations (random: +102, −285, −333; M-band: −2,
   +89, −417), spread ±150 pN. The single-seed reference (+131 ± 8) is therefore **not
   robust** — it was one favourable realization. This undermines the per-SF magnitude read
   itself, not just the M-band variant. This is a **magnitude-characterization** matter →
   **HALT → PI** (do not self-resolve): needs either a seed-ensemble at full GPU statistics
   (150k/150k/30samp × many seeds) or a measurement-protocol fix before any per-SF traction
   number is quoted.

## Caveats

- CPU preliminary at slightly reduced statistics (120k/120k/20 vs reference 150k/150k/30);
  3 seeds. Conclusion (1) "no improvement" is robust to this; the absolute magnitudes and
  the instability finding (2) should be reproduced at full GPU statistics before quoting.
