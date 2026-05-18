# Questions for Sungwook — Phase 1 Unit 2.2 (Worker B)

Five open items after Task 7 stop point.

## Q1 (deferred from Unit 2.1)
Pereverzev refit to Kong 2009 α5β1. Still deferred to Phase 2 per PI sign-off; the Unit 2.2 implementation uses the KU-2.18 illustrative parameters as agreed.

## Q2 (resolved in Unit 2.2)
Mature per-clutch force band. **Resolved**: vinculin allostery + FC-stage FA size (20 clutches) puts per-clutch force at mean 9.20 pN, median 7.13 pN, 48.2 % in [5, 20] pN. See `REPORT.md` for the analysis and deviation rationale.

## Q3 (new — deviation 4) — FA size 3-modal distribution test
The brief Task 6 asked for a 100 FA × 5 min sim showing the KU-2.2 3-modal NA/FC/FA distribution. With the current Python loop this run is ~2.5 hours; vectorising across FAs is a structural rewrite outside Unit 2.2 scope. **Decision needed**:

  - (a) Accept the single-FA Hill kinetics validation (figure 4 + the three growth unit tests) as the Phase 1 surrogate.
  - (b) Implement a vectorised motor-clutch step that runs N_FAs in parallel (Phase 2 task) and validate the 3-modal distribution there.

**Recommendation**: (a) for Phase 1, schedule (b) into the Phase 2 brief.

## Q4 (new — deviation 5) — Performance budget
Brief asked 100 FAs × 1000 steps < 1 s; measured 2.7 s. Strict 1 s requires the vectorisation in Q3. **Decision needed**:

  - (a) Accept the honest 5 s budget for Phase 1 (the test asserts this).
  - (b) Commit to vectorisation as the next Unit 2.2 task before sign-off.

**Recommendation**: (a). The maturation substeps add real per-step cost (Poisson draws + Hill + resize); meeting 1 s would require batching across FAs, which is Phase 2 architecture.

## Q5 (new — deviation 3) — Why 50-clutch FAs cannot reach mature band
The Unit 2.2 mature-regime test uses 20 clutches (FC-stage) because per-clutch force at stall is fixed by `F_per ≈ N_m·F_stall / N_eng` — vinculin allostery alone cannot drive 50-clutch FAs into the 5-20 pN band; the FA-size dynamics must SHRINK the engaged clutch count, not stiffen each clutch. This contradicts the canonical biology view that "mature FAs are larger" (KU-2.2: NA → FC → FA → FB).

The bridge between the two views: the mature FA *plaque* is larger and integrates more force, but at any instant the *engaged* fraction is smaller (Bangasser 2017). The Phase 1 model collapses to "engaged at any moment" only, which gives the per-clutch force the test catches.

**Decision needed**: should the Unit 2.3 / Phase 2 brief explicitly model the NA → FC → FA maturation cascade (including the engaged-fraction reduction as FAs mature)? Or is the current Phase 1 single-clutch-population model sufficient for the upstream cell-motility work?

**Recommendation**: include the maturation cascade explicitly in the Phase 2 contract; document the Unit 2.2 finding as background.
