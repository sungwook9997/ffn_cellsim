# FF active movement — piece 3/5: FA-clutch substrate traction (integrin catch-slip, GPU) 2026-07-02

**Engine:** FF (Warp, µm·pN·s) **Branch:** dcm/main **Scope: 3 of 5.** Piece-1 protrusion, piece-2 polarization,
piece-3 = TRACTION: how the cell grips the substrate/ECM and pulls (the molecular-clutch model,
Mitchison-Kirschner 1988; Chan-Odde 2008). De-adhesion is EMERGENT from load-dependent catch-slip bond rupture,
never a binary latch ([[feedback-junction-switch-fine-grained]]).

## What this is

`ff/fa_clutch_warp.py` — the GPU-resident port of the host-side `fa_anchor.py` clutch: two Warp kernels keeping
the whole cell GPU-native.
- **`clutch_spring_kernel`** — an engaged clutch is a Hookean spring k_int·(L−rest)·û pulling an actin END node
  toward its FIXED substrate anchor (reaction goes to the immovable ECM = traction on the cell).
- **`clutch_catchslip_kmc_kernel`** — one KMC tick per clutch: a bound clutch on load F=k_int·|L−rest| detaches
  with 1−exp(−τ·off_rate(F)); off_rate = Pereverzev two-pathway `k_catch0·exp(−F·x_c/kT) + k_slip0·exp(+F·x_s/kT)`
  (catch-then-slip); a detached clutch within the capture radius re-attaches with 1−exp(−τ·k_on).

## Validation (analytic ground truth first — oracle-is-crosscheck)

`tests/ff/test_fa_clutch.py` (3/3):
- **Spring force** = k_int·(L−rest) toward the anchor, exact; detached clutch exerts 0.
- **Catch-slip KMC** — the Warp kernel's empirical off-rate (bound-fraction decay, k_on=0 ensemble) == the
  analytic Pereverzev off-rate to **rel <5%** across F=0–30 pN (measured 0.0–0.2%).
- **Catch-then-slip signature** — off-rate is minimal at the analytic peak F* = kT/(x_c+x_s)·ln[(k_c0 x_c)/(k_s0
  x_s)]: off(0)=0.90 > off(F*)=0.78 < off(30)=1.37 (the bond STRENGTHENS with load up to F*, then slips).

## Demo — emergent traction (Chan-Odde motor-clutch)

`scripts/ff_clutch_traction.py` (`outputs/ff/figs/clutch_traction.png`): retrograde actin flow (v_retro=0.02
µm/s) drags the actin against 200 integrin clutches on a fixed substrate; each stretches → loads → catch-slip
ruptures → re-engages stress-free. Emergent **traction ≈ 1.48 nN (200 clutches, ~15 pN mean per bound clutch,
bound fraction 0.49)**, fluctuating from the stochastic (desynchronized) catch-slip rupture/rebinding — NOT a
latch. Order-of-magnitude cross-check: scales to tens of nN at the ~thousands-of-adhesions whole-cell level
(Gil-Redondo MCF7 ~102 nN; a cross-check, not a gate — H.7 SF-array is HALTED).

## ⚠️ SURFACED TO PI — F* = 7 pN vs KB-claimed 30 pN (integrin catch-slip peak)

The recorded INTEGRIN_A5B1 Pereverzev params (k_catch0=0.4/s, x_catch=0.61 nm, k_slip0=0.5/s, x_slip=0.14 nm;
hand_kmc, KB-2.5, "Kong 2009") give a catch-slip peak **F* = 6.99 pN**. But Kong et al. 2009 (α5β1-fibronectin)
report the catch bond strengthening up to ~**30 pN**, and the KB claim says 30 pN. **These are inconsistent.**
Per the no-magic-number rule I did NOT retune x_catch/x_slip to hit 30 pN (that would be fitting the params to a
target). The kernel is validated against the analytic off-rate AS RECORDED (7 pN). **PI decision needed:** are
the recorded Pereverzev x-values wrong (should reproduce Kong's 30 pN peak), or is the KB-2.5 "30 pN" claim the
error? This is a KB-content question (ModelContract / Parameter), gate-contract change → PI sign-off.

## Honest scope + next
- Kernel + analytic validation + a 1-D retrograde-flow traction demo. **The full cell integration** — FA clutches
  on the basal cortex/filopodium end nodes to a substrate, traction map, catch-slip de-adhesion during
  spreading/migration, rendered as CELL MORPHOLOGY on a substrate — is the next increment (native run).
- Constants PI-gated (INTEGRIN_A5B1 k_int, k_on, catch-slip x-values + the F* discrepancy above).

## Files
- `ff/fa_clutch_warp.py`, `tests/ff/test_fa_clutch.py`, `scripts/ff_clutch_traction.py`,
  `outputs/ff/figs/clutch_traction.{png,json}`.

Related: [[project-ff-active-movement-pieces]], fa_anchor.py (host scaffolding), FF_POLYMERIZATION_PIECE1,
FF_POLARIZATION_PIECE2, oracle-is-crosscheck.
