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
- **A2. PI poster overlay (the payoff).** Overlay the PI's measured A/A₀(R) on the emergent
  curves (OVERLAY-ONLY, never fit — hard rule). Does the platform's emergent a/b/c match the
  experiment's a/b/c (qualitatively / quantitatively)? = the validation comparison.
- **A3. Statistics.** Current fit = 5 R₀ × 3 seeds × 2 doublings; the ligand fits used only
  3 R₀ (degenerate r²=1). More R₀ / seeds / biological time → a/b/c with real error bars
  (needs B2 for scale).

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
B1 ✅DONE ──► A1 ✅DONE (ligand→traction, 3 separated curves) ──► A2 (PI overlay) = experiment reproduced
                                                                  A3 (more R₀/seeds, needs B2/GPU) ──┘
C (integration) and D (extensions) run in parallel / after.
```

**Recommended next single step (A1 done):** **A2 — overlay the PI poster A/A₀(R) on the three
emergent curves** (OVERLAY-ONLY, never fit). Does the platform's emergent ordering/shape match
the experiment's Bare/Pre/Lam4? Key A2 question carried from A1: the **single-cell↔collective
laminin split** — A1's single-cell traction puts Lam4 ≈ Bare; if the poster shows Lam4 spreading
*most* (collective enhancement), that gap is the finding (it implicates a collective mechanism —
cohesion modulation / uniform-β1 proliferation — not single-cell traction). A2 needs the PI
poster numbers (ask PI). **A3** (more R₀/seeds for robust per-condition a/b/c error bars) needs
B2/GPU and runs alongside.
