# SF traction seed-instability cross-check — independent estimator vs coherent FA-traction

**Date:** 2026-06-11  **Host:** Mac (mp14, CPU dev).  **Scale:** loop24d PRODUCTION config
(sarcomeric `--n-filaments 6` → 36 total filaments / 12 FA anchors, n_motors=16, n_sarcomeres=3,
equil=120k, contract=120k, n_samples=20, seeds 1–3).  **Status:** DONE — boot autonomous step 2a.

## Question (boot autonomous step 2a)
The per-SF coherent FA-traction differential (`h7_ventral_sf_traction.py:anchor_traction`, myosin
force ON−OFF) is seed-unstable: the single-seed "+131 pN" did not reproduce (loop24d: paired
Δ +62±156 pN, NS; sign-flips even across the 20 time-samples of one run). Is that instability
**physics** (the realised bond-tension STATE genuinely fluctuates — the myosin signal is below the
per-realisation noise floor) or an **estimator artefact** (the specific left/right anchor-bond split
in `anchor_traction` fabricates it)? The new Slater radial stress(r) estimator
(`ffn_sim/cortex/radial_stress.py`, `b546c98`) gives an INDEPENDENT projection of the same state.

## Method
The loop24d SF states were never checkpointed (only scalar observables are in `mband_ab/*.json`),
so the same-seed paired ON/OFF design was **re-run** (`scripts/h7_sf_traction_crosscheck.py`),
computing on every matched sample THREE readouts on the SAME states:

1. **coherent** — the existing observable (anchor bonds, left/right inward split).
2. **midplane** — INDEPENDENT: net axial tension `Σ T·|û·axis|` across the bundle mid
   cross-section over ALL backbone bonds (no left/right split, different bonds, different spatial
   location). Force balance ⇒ it reads the SAME physical tension the anchors feel, so if it tracks
   `coherent` sample-by-sample the instability lives in the STATE, not the aggregation.
3. **σ(r)** — the new radial estimator on the SF bundle (≈1D along the axis → spherical shells
   sample the axial-tension profile).

**Config-match (critical).** The replication MUST use `--n-filaments 6` (loop24d
`mband_ab/run_ab.sh`: 6 per-cross-section → 36 total filaments / 12 anchors). A first run at the
script's old default 36 gave a 6× bundle (coh|abs| ~3200 vs ~860 pN) and a spurious "robust −160 pN"
— flagged and corrected before any conclusion. The matched run below reproduces loop24d's
mband-OFF arm to the decimal, validating the harness.

## Results
γ here is force [pN]; myosin force ON − OFF, same-seed paired.

| seed | coherent ON−OFF [pN] | midplane ON−OFF [pN] | loop24d off_diff [pN] |
|---|---|---|---|
| 1 | **+101.7 ± 6.8** | +116.4 ± 35.1 | +101.66 |
| 2 | **−285.2 ± 33.1** | −793.6 ± 96.9 | −285.20 |
| 3 | **−333.1 ± 18.5** | −322.0 ± 47.7 | −333.12 |
| ensemble | **−172.2 ± 112.4** (noise-dominated) | **−333.1 ± 214.6** (noise-dominated) | −172±138 |

- **coherent reproduces loop24d's off arm exactly** (+101.7/−285.2/−333.1 vs +101.66/−285.20/−333.12)
  → the re-run is a faithful replication.
- **The independent midplane estimator flips sign on the SAME seeds** (+,−,−), and per-sample
  **Pearson r = 0.694** with coherent → the two read the same fluctuating state.
- **BOTH ensemble differentials are noise-dominated** (|mean| < 2·sem): coherent −172±112,
  midplane −333±215. midplane is noisier (more crossing bonds, sensitive to the diverged config)
  but tracks sign.

## Verdict — BOOT HALT CORROBORATED
The seed/sample instability is **NOT an artefact of the coherent left/right split**: an independent
estimator (different bonds, different location, no split) shows the **same sign-flipping,
noise-dominated** behaviour and tracks the coherent one (corr 0.69). The instability is a
**physical property of the realised state** — myosin force ON and OFF trajectories diverge under a
weak contractile perturbation, so the ON−OFF differential is a difference of two diverged ~hundreds-
of-pN configs. **No robust per-SF traction number exists by either readout → it is not fixable by a
better estimator.** loop24b "generation/engagement-limited" stands, now on an independent estimator;
the boot HALT (full GPU seed-ensemble or a measurement-protocol change before any per-SF traction is
quoted) is the right call.

## Figures
- `sf_traction_crosscheck_full.png` —
  (1) per-sample coherent vs independent midplane differential (Pearson r=0.69, identity line);
  (2) per-seed ON−OFF differential, both estimators, with SEM — sign flips together across seeds;
  (3) Slater σ(r) on the SF bundle (seed 1, ON vs OFF) — the new `radial_stress` tool on live state.

## Next
- This used the **scalar** σ(r) tool qualitatively (panel 3); the decisive cross-check is the
  midplane vs coherent agreement. The remaining boot HALT item (per-SF seed-instability resolution)
  is unchanged: a full GPU seed-ensemble (gbook) or a protocol that does NOT difference two diverged
  trajectories (e.g. read the engaged-myosin force directly) — **PI-scoped**, magnitude
  characterisation.
- Independent of this: autonomous step 1 (cm-z-struct density sweep, gbook GPU) remains open.
