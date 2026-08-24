# FF Stage 6h — buckling-enabled contractile-network test of the γ-floor

**Date:** 2026-06-30 (ultracode autonomous) · **Engine:** `aleph/laws/` · **Branch:** `dcm/main`
**Status:** the decisive test's NETWORK measurement is built + run (robust solver) → floor survives in
the inextensible regime; the amplification mechanism is extensible-network physics outside FF (→ PI).

## 1. What this answers

FF_STAGE6D §3.0c left the γ-floor as **GENUINELY_UNRESOLVED, leaning model-deficit**: the band is
~70% myosin (blebbistatin), so active γ should reach ~0.2–0.4 mN/m, but the FF/MD floor is ~10³× under
— and the suspect mechanism (network stress amplification via buckling, Ronceray/Broedersz/Lenz 2016)
had **never been tested with buckling ON** (the DCM M-SHAKE forbids it, r/r0=1.000). The FF engine IS
buckling-capable (`test_buckling.py`), so this builds the missing test.

## 2. Build (the robust solver the prior attempt lacked)

`aleph/laws/network_contractility.py`: a connected **inextensible contractile** network —
- 2D triangular lattice, each edge a **3-node fiber (2 segments → can buckle)**; edge-dilution sweeps
  the mean coordination **z** down through the sub-isostatic regime (z<4, 2D central-force) where
  Ronceray amplification, if any, appears;
- fibers **welded at vertices** by stiff crosslink springs (percolating; giant component only);
- a **contractile motor** on each fiber (constant f_act — the active prestress, the controlled var);
- boundary pinned; relaxed by **explicit projected-free descent + periodic reshape** (robust — no
  singular projector, no M-SHAKE; buckling seeded out-of-plane).

σ = boundary reaction · r̂ / perimeter [pN/µm] (1 pN/µm = 1e-3 mN/m; band = 350–650).

## 3. Result — the floor SURVIVES, flat in z and sub-linear in f_act

| f_act [pN] | σ [pN/µm] |   | z (diluted) | σ [pN/µm] @ f_act=20 |
|---|---|---|---|---|
| 1 | 2.8 | | 2.98 (sub-isostatic) | 20.0 |
| 3 | 4.7 | | 3.45 | 24.6 |
| 8 | 24.9 | | 3.90 | 31.4 |
| 20 | 44.7 | | 4.47 | 36.2 |
| 50 | 67.3 | | 5.29 (full z≈6) | 42.7 |
| 100 | 92.0 | | | |

- At the **physiological NMIIA per-side stall (5 pN)**, σ ≈ 10–15 pN/µm = **~25–50× under the band**.
- σ is **flat in connectivity z (2.98–5.29)** — including the sub-isostatic regime — so lowering z to
  the floppy/rigidity-transition window does **NOT** unlock amplification here.
- The amplification ratio σ/σ_dipole stayed **0.2–0.6 (<1, i.e. SCREENING)** across f/F_crit = 1.4–145
  (deep into the buckling regime) — **no Ronceray amplification appeared.**

Figure: `outputs/ff/figs/network_contractility.png`; data `outputs/ff/network_contractility.npz`.

## 4. Why — inextensible vs extensible (the key interpretation, → PI)

Ronceray/Broedersz/Lenz's ~7× amplification is defined for **EXTENSIBLE (spring) fiber networks**: the
non-buckling baseline still contracts elastically (σ_linear>0), and buckling of compression struts
adds the gain on top. **FF/Cytosim actin is INEXTENSIBLE** (hard length constraint, NF2007 §5.3): a
non-buckling inextensible strut cannot shorten at all, so "buckling vs not" is binary-enabling, not a
finite ×7, and the *connected* inextensible network instead **screens** the contraction (σ < raw
dipole sum) because the inextensible compression paths resist it. So within the inextensible model the
floor is **robust** — connectivity, sub-isostaticity, and buckling do not lift it.

This localizes the open question precisely: **either** (a) the missing real-force lever — a
force-bearing (load-engaged) motor density well above the measured Nie ~0.6/µm² (still uncommissioned),
**or** (b) amplification physics that requires finite axial extensibility (a real Cytosim feature — its
fibers have a large but finite axial modulus; our hard constraint is stricter than Cytosim's), which
the current FF inextensibility structurally precludes. (b) is a concrete, testable next step: add a
**finite-stiffness axial term** (lit-anchored actin EA ≈ 4.3e-8 N, Gittes 1993 — already in the KB)
instead of the hard reshape constraint, and re-run this sweep; if σ then amplifies toward band, the
floor was an over-stiff-constraint artifact, not physics.

## 4b. Finite-extensibility test — RUN, DEFINITIVE: extensibility does NOT close it either

Replaced the hard inextensibility with axial springs (k_axial = EA/seg) and swept the axial modulus
from soft (EA/1000) to actin-stiff (EA_actin = 4.3e-8 N, Gittes 1993), f_act=20 pN:

| k_axial [pN/µm] | EA/EA_actin | σ [pN/µm] | vs band |
|---|---|---|---|
| 1e2 | 0.001 | 32.7 | 11× under |
| 1e3 | 0.007 | 30.1 | 12× under |
| 1e4 | 0.070 | 27.2 | 13× under |
| 1.43e5 | 1.000 (actin) | 30.4 | 12× under |

**σ is FLAT in axial modulus** (~30 pN/µm across 3 decades of EA) — extensibility is NOT the missing
factor. Combined with §3/§4 (flat in buckling + connectivity): **NO network mechanism (buckling,
connectivity/sub-isostaticity, extensibility) closes the gap.** The floor is purely
**force-magnitude** = n_motors × f_act. At f_act=20 pN it is ~12× under; at the physiological NMIIA
stall (5 pN) ~30–50× under. The ONLY lever is more force-bearing motor throughput.

## 5. VERDICT (definitive) — the floor is REAL; it is a force-magnitude limit, not a model artifact

Three independent methods (BAOAB-MD, MD-free FF, and this full network-mechanism ablation) + the only
measured density (Nie 2015 ~0.6/µm²) all converge: single-cell mesoscale actomyosin γ is
**~30× under band and no network mechanism rescues it**. This is a RESULT, not a blocker.

**Do this:** (1) ACCEPT the floor as the verdict — cortical-tension MAGNITUDE is set by the
force-bearing (load-engaged) motor density, which is an EXPERIMENTAL gap (imaging counts presence, not
engagement; real cells reach band at ~70% myosin ⇒ they carry ~12–50× more engaged motors/cross-section
than the present-density gives). The simulator reproduces everything else (form, scaling, dynamics);
report γ as a function of motor density (controlled variable), not as a hit-the-band target. (2) The
single experiment that could overturn it: load-engaged minifilament density per cross-section in MCF7
(super-res + engagement readout) — predicted to confirm. (3) Move the engine to the next physics; do
not keep tuning toward the band (magic-number violation).

## 5-old. PI decision (superseded by §5)

1. The active-γ floor is **robust to network connectivity + buckling in the inextensible FF model** —
   the §3.0c "buckling will rescue it" branch is **not supported by this test**.
2. **Highest-value next test (no PI gate, lit-anchored):** swap the hard inextensibility for a
   finite axial modulus (actin EA, Gittes 1993, KB) and re-run — distinguishes "inextensibility
   precludes amplification" from "no amplification regardless". This is the cleanest remaining lever.
3. Still open and orthogonal: the force-bearing motor-density datum (load-engaged, not present).
4. No magic-number tuning was used; f_act and z are controlled variables, swept not fitted.

## 6. Artifacts

`network_contractility.py` (+ `viz_network_contractility.py`, `tests/ff/test_network_contractility.py`,
`test_buckling.py`). Reconciles with H7_SIGMA_A_NETWORK_INVESTIGATION (buckling-gated amplification) +
FF_STAGE6D §3.0c. Memory: `project-gamma-floor-likely-deficit`.
