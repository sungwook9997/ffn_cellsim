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

- **A1. Ligand mechanism (Bare/Pre/Lam4) — the real driver** ⭐ TOP PRIORITY.
  Passive substrate adhesion alone barely separates the conditions (L2.6 finding); the actual
  driver is ACTIVE traction (validated: traction↑→A/A₀↑). Anchor ligand→traction via the
  per-species integrin catch-slip in `bridge/ligand_species.py` (col-I vs laminin-111 → engaged
  clutch count → traction magnitude) so the three conditions yield three distinct a/b/c curves
  emergently. This is the experiment's actual variable. ⚠️ Anchor the relative-traction ordering
  carefully (laminin notes say LOWER traction than FN; don't guess the Bare/Pre/Lam4 order — use
  the measured relative tractions / clutch kinetics). Depends on B1.
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
B1 (box/eject-guard) ✅DONE ──► A1 (ligand→traction anchor, ≤~3 nN band) ──► A2 (PI overlay) = experiment reproduced
                                 A3 (more stats, needs B2/GPU) ──┘
C (integration) and D (extensions) run in parallel / after.
```

**Recommended next single step (B1 done):** A1 — anchor the Bare/Pre/Lam4 ligand conditions to
per-species integrin catch-slip → edge-traction *within the stable ≤~3 nN band B1 established*,
yielding three distinct a/b/c curves emergently. Then A2 (PI overlay) = the next headline
("the platform reproduces the PI Bare/Pre/Lam4 conditions mechanistically").
