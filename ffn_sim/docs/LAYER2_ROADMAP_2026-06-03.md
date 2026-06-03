# Layer-2 spheroid line — remaining roadmap (2026-06-03)

> Recorded for the next session to decide + continue from. Branch `layer2/spheroid-cbm`.
> Current state: **L2.0–L2.6 DONE; G3 PASS (the PI law A/A₀=a+b/R+c/R² emerges, r²=0.98)**.
> Boot: `project-layer2-spheroid-line` memory + `outputs/layer2/REPORT.md` (headline) +
> `docs/LAYER2_ANCHORS_2026-06-02.md` (provenance, incl. the catch-bond SI inferences).

## Where we are (done)

- L2.1 CBM + G1; L2.2 edge-traction; L2.3 A/A₀ sweep (cohesion-locked finding); L2.4
  contact-inhibited proliferation + G4; L2.4a connected-core observable; **L2.4b leak-free
  pooled growth** (fixed the HOOMD per-rebuild leak, 7 GB→240 MB).
- **L2.5 E-cadherin catch-bond cohesion** (faithful Rakshit-2012 sliding-rebinding) → resists
  proliferation fragmentation (catch 0/5 vs morse 1/5) → **G3 PASS r²=0.98** (A/B: morse
  r²=0.80 FAIL on the same path → catch is necessary).
- L2.6 substrate wall (MCF7 = 3D cap, correct low-invasion phenotype); ligand axis (passive
  adhesion → mild separation); **active traction** validated as the spreading driver
  (A/A₀ b-coefficient grows with traction = mechanism→a/b/c).
- 58 layer-2 tests green; 12 figures; gates G1/G3/G4/G5/G6.

## Remaining — prioritized, with dependencies

### A. Make it experiment-comparable (the scientific completion)

- ~~**A1. Ligand mechanism (Bare/Pre/Lam4) — the real driver**~~ ✅ **DONE 2026-06-03**
  (REPORT §A1; `spheroid/ligand_traction.py` + `scripts/layer2_ligand_traction_sweep.py`;
  fig `fig_layer2_ligand_traction_conditions.png`; +8 tests, 96 green).
  Mapped each condition → ACTIVE edge-traction `f = T_ref·density·(φ·F_s)` from the
  `bridge/ligand_species.py` clutch kinetics: **col-I clutch 1.0 vs laminin 0.61× (measured
  ordering, not guessed)**; Bare<Pre = flagged density axis. Resolved Bare 1.50 / Pre 2.50 /
  Lam4 1.53 nN (all in the B1 ≤3 nN band). **Result: the active mechanism SEPARATES the
  conditions — Δ(A/A₀)≈0.14 at mid-R₀, ~4.6× the L2.6 passive (~0.03); ordering Pre>Lam4≳Bare
  tracks resolved traction.** ⚠️ Honest limits → feed A2/A3: per-condition a/b/c are
  UNDER-DETERMINED (4 R₀, 1 dof — curves+separation are robust, the a/b/c split is not; need
  A3); absolute T_ref unanchored (relative ordering is the science); density axis pending the
  collaborator pV4D4 datum; **single-cell↔collective laminin split** — Lam4 single-cell traction
  is low (weak laminin clutch) so it does NOT out-spread via single-cell traction; the poster's
  *collective* Lam4 enhancement (if shown) is that split, the A2 question (do NOT engineer it).
- ~~**A2. PI poster overlay (the payoff).**~~ ✅ **DONE 2026-06-03** (REPORT §A2;
  `scripts/layer2_pi_overlay.py`; fig `fig_layer2_pi_overlay.png`; PI CSVs gitignored/local).
  PI provided `references/260313_{Bare,Pre,Lam4}.csv` (59 spheroids total, R₀ 87–419 µm, ~82 h).
  **⭐ HEADLINE: the novel law's SHAPE is reproduced** — PI A/A₀ DECREASES with R₀ in all 3
  conditions (corr −0.81…−0.91, b>0) = the platform's emergent 1/R. Three honest gaps, each →
  a next step: (1) **ordering SPLIT** — PI collective Lam4>Pre>Bare vs model single-cell
  Pre>Lam4≳Bare = the single-cell↔collective laminin split CONFIRMED → Lam4's collective
  enhancement is NOT single-cell traction (→ a collective mechanism: uniform-β1 proliferation /
  cohesion-modulation lifting the c-penalty = next hypothesis, "D/A4"); (2) **magnitude ~4–5×**
  (PI raw-segmented area + 82 h vs platform connected-core + 60 h → recoverable: raw-area readout
  + longer time); (3) **R₀ range no overlap** (PI 87–419 µm ≫ model 40–78 µm → native-N GPU, B2).
- **A3. Statistics.** Current fit = 5 R₀ × 3 seeds × 2 doublings; the ligand fits used only
  3 R₀ (degenerate r²=1). More R₀ / seeds / biological time → a/b/c with real error bars
  (needs B2 for scale).

- ~~**A4. Collective-Lam4 mechanism (uniform β1 → traction localization).**~~ ✅ **DONE
  2026-06-03** (REPORT §A4; `ligand_traction.py` β1→Lp axis + `scripts/layer2_a4_uniform_beta1.py`;
  fig `fig_layer2_a4_uniform_beta1.png`; +2 tests, 98 green). Tested the A2-implicated collective
  hypothesis: Lam4's measured "uniform β1" → uniform (Lp≫spheroid) traction (vs Bare/Pre edge,
  Lp 11 µm), magnitude unchanged. **Direction CONFIRMED**: uniform β1 engages the interior
  (∝volume) → flips Lam4's size-dependence (corr −0.99→+0.95 = small-size penalty removed/
  reversed, the documented "scale-independent" Lam4 phenotype emergent) and lifts Lam4 above Pre
  at large R₀ (crossover ≈67 µm) — the collective resolution of the single-cell↔collective split
  (Lam4 enhancement = uniform-β1 collective engagement, NOT single-cell clutch traction). ⚠️ But
  full-uniform OVERSHOOTS: PI Lam4 still decreases (corr −0.91) while full-uniform reverses the
  slope; + large seed variance (±0.6–1.2, near cohesion-destabilization). → **A4′ next: partial
  β1 uniformity (intermediate Lp sweep, overlay-only) + more seeds (B2/GPU)** = the physical Lam4.

### B. Numerics / production robustness (needed to run A at scale)

- ~~**B1. Box-sizing + substrate-wall fix.**~~ ✅ **DONE 2026-06-03** (REPORT §B1).
  Root cause was the pool box being sized for a *3D ball* (`r0·N^(1/3)`) while a substrate-
  confined spheroid *wets into a √N 2D disk* → footprint exceeded the box → PBC out-of-bounds.
  Fixed with `cbm.pool_cluster_radius` (3D-pack vs 2D-wetting + spread safety) in both pool
  builders. ALSO surfaced a genuine model limit: f_traction ≳ cohesion (~6.5 nN) physically
  detaches an edge cell (instant in the overdamped large-dt CBM); the old 6 nN "success" was a
  too-small-box PBC artifact. Added a graceful **ejection guard** (`ejected=True`, no crash, no
  silent truncation). In-regime (≤3 nN) clean; G3 headline unchanged; +4 tests (88 green).
  **Refines A1:** anchor the ligand→traction map into the stable ≤~3 nN band (above = a
  detachment regime needing GPU sub-stepped bonds, D3 — not a numerics patch). The substrate
  Morse-wall repulsion spike was a *downstream* symptom of the xy runaway (now caught by the
  guard), so the wall potential is left untouched — changing it touches physics and would need
  PI sign-off; flag if a future in-regime case ever penetrates the wall.
- **B2. GPU port.** Port `run_growth_pooled` GPU-resident (project hard rule: RTX A5000+ is the
  mandatory production minimum) for native-N + long-time sweeps.

### C. Project integration (couple the two lines as physics)

- **C1. Cross-line consistency seam.** A test asserting the single-cell line and Layer-2 share
  ONE MCF7 anchor set (traction / cohesion / size scales agree — the scale-bridge).
- **C2. PI-ratification queue.**
  - The catch-bond SI inferences (Bell exponent POSITIVE/slip; sliding term uses k₋₁ not k₊₁) —
    both physics/conservation-forced and validated by reproducing the ~30 pN catch peak, but
    flagged for PI sign-off (rationale in `docs/LAYER2_ANCHORS_2026-06-02.md`).
  - Fix `PI_EXP_VALIDATION_MAP.md` (Buckley Δx*=4 nm misattribution; retired pN D_e seed).
  - `ffn/foundation` push approval (Layer-2 work is on `layer2/spheroid-cbm`, not pushed).

### D. Roadmap extensions (after A–C)

- **D1. L2.7 invasion** — 3D-Mikado ECM; cells invading a fibrous matrix (dimensionality, the
  physical meaning of the c-term, invasion phenotype).
- **D2. Surface-tension validation** — emergent aggregate surface tension vs a proxy
  (MCF10DCIS ~21 mN/m, Nagle 2022; proxy-flagged — no MCF7 datum).
- **D3. Catch-bond fidelity upgrade** — the current cohesion is the adiabatic effective-force
  limit (valid: bond kinetics ~0.1 s ≪ CBM dt ~1.7 s). Explicit stochastic bonds would need a
  sub-stepped bond integrator (the 4-OOM timescale gap) — a GPU-era option, not required.

## Dependency graph / recommended next step

```
B1 ✅ ─► A1 ✅ ─► A2 ✅ (SHAPE matches; 3 gaps) ─► A4 ✅ (uniform-β1 = right knob; full overshoots) ─┬─► A4′ partial-uniformity Lp sweep (+ more seeds)
                                                                                                     ├─► magnitude: raw-area + longer time
                                                                                                     └─► B2 native-N GPU ─► A3 robust a/b/c (PI R₀ range)
C (integration) and D (extensions) run in parallel / after.
```

**Recommended next single step (A4 done — collective mechanism validated in direction):**
- **A4′ — partial β1 uniformity (intermediate Lp sweep).** A4 showed full-uniform is the right
  knob but overshoots (slope sign flip + variance). Sweep Lp between edge (11 µm) and uniform to
  find the partial uniformity that lifts Lam4's magnitude / flattens the penalty *without*
  reversing the sign — the physical Lam4. Overlay-only; needs more seeds (variance) → pairs with B2.
- **Magnitude (quick) — raw-area + longer biological time** readout (close the §A2 ~4–5× gap).
- **B2 → A3 — native-N GPU** for the PI R₀ range (87–419 µm) + robust per-condition a/b/c.
Ask PI which to prioritise. (The platform now reproduces the law SHAPE (A2) and the collective
ligand mechanism in direction (A4); the remaining gaps are scale/statistics = the GPU axis B2.)
