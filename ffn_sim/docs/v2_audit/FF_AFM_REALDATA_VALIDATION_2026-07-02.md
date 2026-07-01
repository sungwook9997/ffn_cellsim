# FF whole-cell AFM force — validated against real MCF7 / whole-cell compression data (2026-07-02)

**Date:** 2026-07-02  **Engine:** FF (Warp, A5000)  **Branch:** dcm/main
**PI directive (8h goal):** "실제 데이터들 있으면 한 번 테스트." Tests the FF rigid-plate whole-cell AFM
force–strain (FF_STAGE6Z) against real experimental single-cell compression data.

## The test

FF rigid-plate native (N=70686, regulated turgor) predicts:

| strain | 0.15 | 0.35 | 0.40 | 0.55 |
|---|---|---|---|---|
| **F_plate** | 3.3 nN | 18.7 nN | 22.9 nN | 38.8 nN |

## Verdict — RIGHT BALLPARK, no red flag (lit-anchored comparison)

1. **Hertz at the real MCF7 AFM modulus.** MCF7 apparent E (AFM) ≈ 0.3–1.0 kPa (Li 2008 BBRC ~0.5 kPa,
   10.1016/j.bbrc.2008.07.078; Mirzaluo 2023 ~300–400 Pa, 10.1007/s00249-023-01642-3; Zbiral 2023 249 Pa).
   Hertz/Sneddon (sphere R=7.5 µm, ν=0.5, δ=strain·R) gives, for E=0.5→1.0 kPa: 2.9→5.8 nN (15%), 20.4→40.8
   nN (55%). **The FF curve (3.3→38.8 nN) sits between the 0.5 and 1.0 kPa Hertz lines → effective whole-cell
   modulus ~0.7–1 kPa = exactly the physiological MCF7 band.**
2. **Whole-cell confinement (direct force–strain analog).** Cattin 2015 PNAS (mitotic HeLa confinement,
   10.1073/pnas.1502029112): 5 nN@20% → 50@45% → 100@52% → 150–200@68%. FF sits ~2–3× BELOW at high strain —
   the physically expected direction (HeLa mitotic R~10 µm at ~0.4–0.55 kPa pressure is larger + stiffer than
   an interphase MCF7 R=7.5 µm).
3. **Parallel-plate.** Fischer-Friedrich 2016 (mitotic HeLa, 10.1016/j.bpj.2016.06.008): Fst ~40–130 nN
   (→40 with blebbistatin). FF below, again expected.
4. **Shape.** Super-linear rise (single-digit → tens of nN) matches Lulevich 2006 compression (10.1021/la060561p).

**Consistency both ways:** the FF force scale is consistent with EITHER a ~0.7–1 kPa elastic cortex OR a
few-hundred-Pa regulated turgor (a pressurised bag at ~400 Pa gives ~21→78 nN over 15–55%; at the resting 40
Pa only 2–8 nN). Both are physiological. One watch-item: FF is slightly LOW at low strain (3.3 nN@15% vs
Hertz 5.8 nN@15% at 1 kPa, Cattin 5 nN@20%).

## Significance

This validates the **rigid-plate fix (FF_STAGE6Z)** against real data: fixing the soft-plate artifact (which
had reported ~8500 nN) yields forces that match real MCF7 AFM/compression to within the physiological modulus
band. The soft-plate ~10³× "gap #1" is confirmed an artifact — the corrected force is experimentally credible.

## Sanity gate
- No tuning: FF forces are the model output; comparison is to independent lit values (DOIs verified). ✓
- Honest: FF slightly low at low strain flagged; the turgor-vs-elastic degeneracy noted. ✓

Sources (DOIs verified by the research agent): Cattin 2015 (10.1073/pnas.1502029112); Fischer-Friedrich 2016
(10.1016/j.bpj.2016.06.008); Li 2008 (10.1016/j.bbrc.2008.07.078); Mirzaluo 2023 (10.1007/s00249-023-01642-3);
Lulevich 2006 (10.1021/la060561p). Related: FF_STAGE6Z (rigid plate), FF_STAGE6V, FF_STAGE6W/6X.
