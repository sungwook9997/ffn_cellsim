# WCA contact-graph Schwarz — CPU spectral reference (predictor for Open-work #1)

Status: **host predictor, not a native pass.** This is a pure-NumPy reference (`ac/cell/wca_schwarz_spectral_reference.py`)
that answers, before any device implementation cost, the single question left open by the contact-cluster
Schwarz audit in `REPORT.md`:

> *Open numerical work item 1 — "Extend the proved pair Schwarz space into a multi-fiber star that
> simultaneously contains each WCA edge and its distinct crosslink neighbors … merely adding crosslink pairs as
> independent overlaps was experimentally falsified."*

The harness predicts **the star subdomain will help**: it removes ~71–78 % of the coupled residual that the
retained rigid-fiber block leaves — *provided the star stays local*. It does not claim full-native convergence;
it says the device star block is worth implementing and states exactly which native run validates it.

## Method (faithful to the tournament, not an idealized preconditioner)

A reduced multi-fiber contact network reproduces the documented plateau structure: parallel fibers in a bundle
with a stiff radial WCA contact chain (fiber `f`↔`f+1` at their centre nodes) and Hookean crosslinks from a
contacting fiber to a **distinct** non-adjacent fiber (`f`↔`f+2`), plus fiber bending, NF2007 inextensibility,
and the `sqrt(eps64)·k_max` regularizer. The local tangents are dense mirrors of the production
`contact_schwarz._wca_tangent` / `_central_tangent` (WCA `radial·uuᵀ`, central `transverse·I + (k−transverse)·uuᵀ`),
`k_max = 2 828 099.29 pN/µm`, crosslink `8.2e5 pN/µm` (Ferrer). Nothing is tuned to a convergence outcome.

The metric mirrors how `contact_schwarz` actually uses the correction — as **one residual-monotone tournament
direction**, not a symmetric PCG preconditioner. For the documented plateau residual (rigidly translating a
WCA-contacting fiber along the contact normal, which stretches both its WCA contact and its crosslink to the
distinct fiber), we apply the raw additive-Schwarz direction `d = P M P r` and report the best line-searched
reduction `1 − ‖r − α*Ad‖/‖r‖`. An idealized `cond(MA)` lens was deliberately dropped: naive additive Schwarz
over large overlapping subdomains is a poor *symmetric* preconditioner regardless of subdomain choice, so it
would misrepresent the tournament method.

The pair path is cross-checked against the production `contact_schwarz.additive_pair_schwarz_oracle`
(`test_pair_path_matches_contact_schwarz_oracle`), so the harness operator is the same one the device solves.

## Result — the star removes the coupled mode the pair/rigid blocks leave

Coupled-residual remainder after one Schwarz direction (lower is better), across network sizes:

| n_fibers | star fiber-coverage | pair remainder | rigid remainder | **star remainder** | star / rigid |
|---:|---:|---:|---:|---:|---:|
| 8 | 0.75 | 11.69 % | 11.71 % | **2.52 %** | 0.215 |
| 12 | 0.50 | 12.09 % | 12.10 % | **3.56 %** | 0.294 |
| 16 | 0.38 | 12.23 % | 12.27 % | **3.56 %** | 0.290 |

The pair and rigid blocks leave ≈12 % of the coupled residual (they omit the crosslink-neighbor fiber, so no
single block spans the coupled cluster). The star block — the WCA edge's two fibers **plus** the distinct
crosslink-neighbor fiber — cuts that to ≈2.5–3.6 %, a **~3.4–4.6× smaller remainder**, and the advantage is
stable as the network grows (star coverage falls 0.75→0.38 while the reduction holds). This is the direct
spectral analogue of the report's contact-cluster hypothesis: the stuck mode is the coupled WCA + distinct
crosslink motion, and only a subdomain that contains all three fibers can absorb it in one solve.

### Load-bearing caveat — the star must stay LOCAL

At `n_fibers = 6` (or the transitive expansion where a star's crosslink reach spans the whole neighborhood,
coverage → 1.0) the star *loses* its advantage (remainder rose to ~96 %): a subdomain that covers the entire
local system degenerates and its single full-block direction is no better than a global solve. This is the same
failure mode as the falsified "independently-overlapped crosslink" variant. **Implementation consequence:** the
device star must be bounded to each WCA edge's *immediate* crosslink neighbors — do not transitively grow it.

## Residual-dominance decomposition

Projected condition number `cond(P(aI+K)P)` with each force family removed (12-fiber network):

| configuration | cond | reading |
|---|---:|---|
| full | 2.74e8 | baseline |
| − WCA | 7.27e7 | WCA is a **primary** driver (removal cuts cond to 0.27×, ~3.7×) |
| − crosslink | 2.65e8 | minor |
| − bending | 2.73e8 | negligible |
| − regularizer | 6.42e18 | the `sqrt(eps64)·k_max` floor is load-bearing for rank |

The WCA-removal ratio is stable across sizes (0.25× at 8 fibers, 0.27× at 12 and 16). This corroborates the
report's WCA-on vs WCA-off plateau gap (2204 vs 21 pN): the stiff WCA curvature is the dominant conditioning
source, and the regularizer is essential to keep the projected operator full-rank.

## What this does NOT establish, and the native run that validates it

- It is a reduced host model, not the 70,686-fiber state. The *magnitude* of the full-native plateau drop is
  not predicted — only that the star direction captures a mode the pair/rigid blocks structurally cannot.
- **Validating native run (gbook A5000):** implement the bounded local star block in `contact_schwarz` (each
  WCA edge's two fibers + immediate crosslink-neighbor fibers, degree-weighted, projected, tournament-accepted),
  then re-run the existing full-native 200-iteration residual-monotone tournament (`native70686_wca_*`). The
  prediction is a plateau **below** the retained `18.2513 pN` rigid-cluster plateau, toward the `0.210670 pN`
  force gate. Do not weaken `sqrt(eps64)`, raise the PCG budget, or lower the population to authorize a step;
  keep the explicit fallback. If the native drop does not appear, re-audit the star locality bound above.

## Validation

- 9/9 host tests (`tests/ac/cell/test_wca_schwarz_spectral_reference.py`), CPU-only (pure NumPy): projector
  symmetric-idempotent, projected operator SPD on the free subspace, WCA/central tangents vs closed forms,
  force-cap zeroing, all local blocks SPD (Cauchy interlacing ≥ reg > 0), **pair path == production
  `additive_pair_schwarz_oracle`**, the star-beats-rigid verdict (`remainder_star < 0.5·remainder_rigid`), and
  the WCA/regularizer dominance ranking.

## Figures

- `figs/fig_ac_wca_schwarz_spectral_reference.png` — (left) coupled-residual remainder per subdomain family
  across network sizes, showing star ≈¼ of pair/rigid; (right) residual-dominance condition numbers on a log
  axis with WCA as the driver and the regularizer holding rank.
