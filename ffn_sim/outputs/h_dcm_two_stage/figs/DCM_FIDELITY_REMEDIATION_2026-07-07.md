# DCM engine fidelity — audit verdict + remediation plan (2026-07-07)

## Verdict (4-agent adversarial audit + synthesis)

**The DCM engine is TRUSTWORTHY infrastructure. The N=2000 confluent run is NOT YET thesis-grade data.**
Trust the pipeline; do not quote its numbers as measured MCF7 mechanobiology yet.

### Genuinely realistic (defensible today)
- Junction **connectivity**: 11 bonded-neighbors/cell (lit 3D packing ~13-14, surface-corrected to ~11 for a finite 2000-cell spheroid), 0/2000 orphaned, ≤1-bond/node invariant holds.
- Cadherin **per-bond mechanics**: active bonds ~28 pN = Rakshit catch-bond optimum (KB-4.2). Emergent form+break churn (~2.7 s lifetime), not a latch.
- **Numerics**: V/V0=1.000, CFL bounded, displacement decaying (settling, not diverging), GPU-native at N=2000/324k nodes.
- The G2 gate **caught** the interpenetration rather than hiding it = engine integrity.

### Red flags (ranked) → remediation

| # | Red flag | Fix | Status |
|---|---|---|---|
| 1 | `fmag` "stress" is a **solver residual** (∝ v, →0 at convergence), does NOT concentrate at junctions | compute real **virial/Cauchy σ=(1/V)Σf⊗r** (cadherin+cortex+turgor+nucleus+contact); rename fmag "residual force" | ▶ DEMO'd (cadherin tension concentrates at junctions ∞× vs fmag 0.65×); save-hook TODO |
| 2 | Physical time **0.12 s** = mechanics, not mechanobiology (5 orders short of min-hr) | implicit large-dt + long run + active remodeling | TODO (deep) |
| 3 | A/A0=V/V0=1.000 is a **construction tautology** (confluent asserts confluence) | run from **dispersed** init → measure emergent compaction | TODO |
| 4 | **G2 interpenetration** 19%R median (membranes cross) | non-overlapping **INSET** + **IPC** log-barrier contact | ✅ knobs wired **+ VALIDATED**: N=200 INSET=0.05+IPC → deep-overlap 37.6%→**2.4%** (15×), median 19%→**6.2%R**, stable |
| 5 | Not converged: residual v 0.82 µm/s ≈ 50× real cell speed | run to residual ≪ 1 bond force, v→0 | TODO (ties to #2) |
| 6 | Cortex γ 10-20× soft (contested) **+** nucleus E_nuc 8-12× stiff | **E_nuc 4700→399 Pa** (audit#19); **γ = SWEEP** (PI 2026-07-07) | ✅ env knobs wired (ENUC=399, GAMMA sweep) |
| 7 | Junction force budget ~1-2 orders low (single-molecule, not bundled) | **bundle_n = N_cad = 100** (ρ_cad·A_junction, KB-4.1/4.11/4.17) (PI 2026-07-07) | ✅ env knob wired (BUNDLE=100) |

## PI decisions (2026-07-07)
- **γ (#6b): SWEEP** 5e-4 ↔ 1e-2 N/m — measure shape/compaction sensitivity, compare to experiment. (Reviewers split: band-center-defensible vs suspended-MCF7-floor → resolve by sweep, not a point pick.)
- **bundle_n (#7): derived N_cad = 100** = ρ_cad(100-500/µm², KB-4.1 Buckley2014) × A_junction(~1µm²), conservative low end; KB-4.17 ratified Phase-1 value. → per-junction force N_cad·⟨F_bond⟩ = 1-10 nN (KB-4.11).

## Code changes landed (this session)
`_gbook_fullcomp_assembly.py` env knobs (all default to the fixed/physiological value):
`ENUC=399` · `BUNDLE=100` · `GAMMA=5e-4`(sweep) · `INSET`(>0.01 non-overlap) · `IPC=1`(log-barrier).

## Next execution (sequenced)
1. **Virial save hook** in `dcm_warp_decohesion.record()` — per-node σ from cadherin+cortex-edge+turgor
   (exact laws known), self-validated by reconstruction vs the sim's `force_d`. New npz keys + honest
   viz channels (real stress; fmag renamed "residual force").
2. **Small validation** (N=200, fast): confirm ENUC/BUNDLE/INSET/IPC behave; check interpenetration
   drops with INSET+IPC; check per-interface cadherin force lands in KB-4.11 1-10 nN band with BUNDLE=100.
3. **γ sweep** at a tractable N (400) across 5e-4 … 1e-2 → shape/compaction sensitivity vs experiment.
4. **Thesis-grade run**: dispersed init + implicit large-dt + long time + converged, with all fixes on.

## Honest claim envelope (today)
CAN claim: "DCM is a stable GPU-native volume-conserving deformable-cell engine with mechanistically
emergent, physically-realistic cadherin catch-bond junction networks at N=2000." CANNOT claim (until
fixed): any stress/tension NUMBER from fmag; any compaction/sorting/rounding result; any quantitative
MCF7 tension magnitude.
