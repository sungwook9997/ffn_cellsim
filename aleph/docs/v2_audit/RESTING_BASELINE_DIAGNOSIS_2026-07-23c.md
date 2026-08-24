# Resting-baseline diagnosis (2026-07-23c) — the cortex form-finding / force-density physics fix is REFUTED at full native; a force-free ERM decouples cortex tension from the membrane, and the resting baseline lacks a conservative cortical-tension mechanism (PI decision)

Branch `codex/ff-ac-codex`. This session executes the PI-approved SHIFT FROM SOLVER TO PHYSICS: after three
preconditioner families were exhausted (23a/23b), the plan was to construct the cortex at a clean Laplace hoop
pre-tension by a force-density / form-finding method so the resting membrane turgor is balanced by cortical
tension (`T = ΔP·R/2`) through the force-free-rest ERM — WITHOUT the counterproductive ERM preload. **Measured
at full native (70,686 filaments / 494,802 actin nodes / subdiv 6, A5000), the approach is DECISIVELY REFUTED
on two independent grounds**, and a pure-NumPy oracle pins the mechanism. The gate stays **OPEN (no
loosening)**; the honest finding is a physics/model-completeness gap that needs a PI decision, not a solver.

## TL;DR (the two refutations + the number that settles it)

| native run (subdiv 6, 70,686, radial ERM, GS-tournament, n_inner=1600, NO ERM preload) | projected residual (gate metric) | converged | raw max\|F\| after solve |
|---|---|---|---|
| NO cortex pretension | **0.4785 pN** | False | 0.7766 (membrane) |
| force-density cortex pretension (T=140 pN/µm) | **0.4840 pN** | False | 159.88 (cortex crosslinks) |

The projected residual is **IDENTICAL** with and without the cortex pretension (0.478 vs 0.484). The cortex
pretension is **mechanically inert for the gate**, and neither run reaches 0.21. **Gate OPEN.**

## 1. The construction that was built and tested (force-density, per-edge — the task's exact prescription)

Set each cortex crosslink's force density `q = t_e/L_e` to the uniform surface value consistent with
`T = ΔP·R/2`. Implemented as the per-edge rest length `r0_e = L_e·(1 − q*/k_e)` with a UNIFORM force density
`q*` (NOT the prior uniform strain, which yields non-uniform q because `k_xl ∈ {4.6e5, 8.2e5}`). Everything is
DERIVED from the sourced ΔP=40 Pa:

- `T_cortex = ΔP·R/2 − γ_mem = 40·7.5/2 − 10 = 140 pN/µm` (γ_mem = LIT-anchored bilayer tension, KB-3.B1.1).
- `q* = 2·A_cortex·T_cortex / Σ_e L_e²` (`T = q/(2A)·Σ L²`). Native: `q* = 479.09 pN/µm`, `t_e = q*·L_e` median
  36.9 pN. **Verified on the real network**: realized `T = (1/2A)·Σ t_e L_e = 140.000 pN/µm` exact; per-edge
  `q_e = 479.093` uniform (min=med=max). The force-density MATH is correct — the failure is physical, not
  arithmetic.

Native structure (why the method cannot work here): `n_xl = 70,686` (ONE crosslink per filament), crosslink
chords are SHORT (median 77 nm, i.e. nearly tangential), only **24.7 %** of actin nodes carry any crosslink,
and only **17.5 %** of the ERM-attached cortex nodes carry one. The crosslink graph is a sparse set of
near-isolated stiff springs, not a balanced tension mesh.

## 2. Refutation A — the pretension raises the RAW cortex force to 160 pN, but is PROJECTED-INERT for the gate

Build-only nodal-force decomposition (raw max\|F\| by region, pN, radial ERM, NO ERM preload):

| state | raw max\|F\| | membrane (radial mean) | actin cortex | ERM-cortex nodes |
|---|---|---|---|---|
| force-free build | 0.7766 | 0.7766 (0.6442) | 0.0356 | 0.0305 |
| + force-density pretension (T=140) | **159.88** | 0.7766 (0.6442) | **159.88** | **122.80** |
| + uniform-strain pretension (same T̄) | 132.14 | 0.7766 (0.6442) | 132.14 | — |

Two things happen, both fatal to the approach:

1. **The membrane residual is UNCHANGED (0.7766).** A force-free-rest ERM transmits nothing, so the cortex
   pretension NEVER reaches the membrane at build. The membrane still carries the full turgor.
2. **The cortex RAW force explodes 0.036 → 159.9 pN** (≈ the crosslink tension `t_e`): each pre-tensioned
   crosslink is an isolated point-couple. BUT this 160 pN lies almost entirely ALONG the actin backbone and is
   removed by the NF2007 inextensibility projector, so it does NOT enter the gate's PROJECTED residual. After
   the full 1600-iteration GS-tournament solve the projected residual is **0.4840 pN — statistically identical
   to the no-pretension floor 0.4785 pN**. The cortex pretension buys the gate NOTHING; it only injects a large
   raw crosslink imbalance that the projector discards.

⇒ **Cortex pretension is the wrong lever.** Because the ERM is force-free, cortex tension (balanced or not) is
decoupled from the membrane; the gate residual is membrane-turgor-limited and does not depend on the cortex
tension state at all.

## 3. Refutation B — the pure-NumPy oracle: edge-spring force density is not conservative; only an area-gradient balances

`formfound_cpu_oracle.py` (host NumPy, no CUDA) compares three tension representations on the SAME geometry, so
the "balanced nodal forces BY CONSTRUCTION" premise can be tested directly:

| mechanism | tangential residual (max) | radial resultant | verdict |
|---|---|---|---|
| **(0) area-gradient `γ·∇A`** on a regular icosphere (= the membrane's `membrane_area_kernel`, `−γ·½·n̂×e_opp`) | **0.015–0.24 pN** | = Young-Laplace 2T/R·A_node (ratio 0.995–0.999) | **BALANCED** |
| (1) **edge-spring uniform force density** on the SAME regular icosphere | **20.5–30.9 pN** | = Young-Laplace | NOT balanced |
| (2) edge-spring on a **sparse cortex-like graph** (27 % nodes, ~1 edge/node) | **250–1270 pN** | — | totally unbalanced |

A clean Laplace hoop tension with balanced nodal forces requires a CONSERVATIVE area-gradient (cotangent /
mean-curvature-normal) tension on a triangulated shell — the gradient of a scalar energy `γ·Area`, tangentially
balanced by construction. An edge-spring force density is NON-conservative: even on a perfectly regular
icosphere it leaves ~30 pN/node tangential (100–1000× the area-gradient), and on the sparse cortex crosslink
graph each pre-tensioned edge is an isolated point-couple. **The plasma membrane HAS the conservative shell
(γ_mem area-tension; native tangential 0.0022 pN). The actin cortex does NOT — it is a sparse Hookean crosslink
network. "Set each crosslink's force-density to the uniform surface value … balanced nodal forces BY
CONSTRUCTION" cannot hold for this cortex.** (`membrane_area_kernel` was read to confirm it is exactly the
area-gradient `−γ·½·n̂×e_opp`; the oracle's `area_gradient_forces` replicates it.)

## 4. The projected residual is pinned at the membrane-turgor floor (~0.48) — confirming 22c/23a/23b exhaustion

Both native runs plateau at ~0.48 pN, `converged=False`, `maxdisp ≈ 2.6e-7 µm` (the solver barely moves — the
residual lives on the soft-membrane mode). This matches the established structural blocker (22c): an
equilibrium exists at subdiv 6 but the soft-membrane mode (`k_erm=4600`) needs a ~0.14 nm drift that the global
`dt_mu = 0.1/kmax`, throttled by the stiff crosslink `k=8.2e5` (a 174× ratio), cannot take. Solver floors, all
NO-preload, all > the 0.21 gate: `erm_jacobi_pure_tournament` 0.5062 (23a), `erm_tournament` 0.5705 (23a),
`erm_gauss_seidel_tournament` **0.4785** (this session, radial ERM, full 1600-budget — the best to date, still
2.3× the gate). No preconditioner family coordinates the soft-membrane / stiff-crosslink split below ~0.48.

## 5. The ΔP–γ finding (Card 3, PI-gated — surfaced, NOT changed)

The membrane residual is EXACTLY the applied turgor, node by node:
`residual_mean = (ΔP − 2γ_mem/R)·A_node = (40 − 2.667)·0.017256 = 0.6442 pN` (measured 0.6442; max 0.7766 at the
12 larger pentagonal-vertex tributaries; `A_node = 4πR²/40962`). It scales linearly with ΔP:
- at the running **40 Pa** (Fischer-Friedrich HeLa proxy): mean 0.644, max 0.78 — **the measured value**;
- at the **MCF7 72 Pa** (Card 3: Hosseini γ_eff=0.27 mN/m → ΔP=2γ/R): mean would be 1.20, max 1.44.

So the residual is CONSISTENT with the 40 Pa proxy and does NOT independently point to 72 Pa — it simply *is*
the applied ΔP. The physics is upstream of the ΔP value: **whichever ΔP the PI-authored pressure-consistency
gate lands, the cortex must carry a matching RESTING cortical tension `γ_cortex = ΔP·R/2 − γ_mem` = 140 pN/µm
(0.14 mN/m) at 40 Pa, or 260 pN/µm (0.26 mN/m) at 72 Pa — and the model has NO mechanism to carry it at the
unbound-myosin resting t0, nor a way to transmit it through a force-free ERM.** Do NOT hand-swap 40→72 Pa; it
does not fix the missing mechanism, it only rescales the same unheld residual.

## 6. Root cause + next step (PI decision — model-completeness, not numerics)

**Root.** The resting turgor (0.644 pN/node = 40 Pa · A_node, on the membrane) can only be balanced by cortical
tension TRANSMITTED THROUGH THE ERM. Two coupled requirements, neither met at the current resting build:
1. the ERM must carry force (a force-free-rest ERM transmits nothing → the membrane residual is untouched by
   any cortex construction — §2); and
2. the cortex must supply a CONSERVATIVE tension that reacts it with balanced nodes — which the sparse Hookean
   crosslink graph cannot (§3), and whose physiological source (myosin) is OFF at the unbound-myosin t0.
Passive crosslink pre-straining is refuted; ERM preload is counterproductive (23a) and forbidden; the solver
cannot reach the soft-membrane equilibrium (§4). The resting baseline is missing its cortical-tension source
AND a force-transmitting membrane↔cortex coupling.

**Candidate resolutions (each a PI call, each with a tradeoff):**
1. **Resting bound-myosin tension (fine-grained-faithful).** Bind a resting NMII population so active cortical
   tension holds the turgor — the physiological source. Tradeoff: motor magnitudes are GAP'd; a "resting bound
   fraction" must be PI-defined; it makes the tension the model is meant to VALIDATE part of the baseline. It
   is the natural home for the Card-3 pressure-consistency gate (the resting bound fraction is exactly what
   sets `γ_cortex = ΔP·R/2 − γ_mem`), and its tension is carried on the ACTIN network (not lumped).
2. **Conservative cortex area-tension shell.** Give the cortex a triangulated `γ_cortex·∇A` shell (like the
   membrane's) so it carries balanced hoop tension. Tradeoff: a continuum surface-tension proxy — a LUMPED
   mechanism the fine-grained principle forbids; makes the crosslink cortex mechanically redundant at rest.
3. **Rebuild the membrane↔cortex ERM as a dense, force-transmitting coupling** (params_i0b1.yaml
   `realization_gap`: "the ERM linkage does not transmit force — off-radial, one tether/node"). Radial pairing
   (tested here) lifts the median alignment 0.92→1 but is inert while the ERM is force-free; physiological ERM
   density (~235/µm² at subdiv 7) would stiffen the coupling — but this still needs (1) or (2) to give the
   cortex something to react with, and ERM density is a PI-GAP.

**Recommendation.** Surface to PI that the resting Laplace balance is a PHYSICS / model-completeness gap: the
membrane turgor cannot be balanced by cortex form-finding through a force-free ERM (proven inert, §2), and the
cortex cannot carry balanced hoop tension by crosslink pre-straining (proven, §3). Option (1) — a PI-defined
resting bound-myosin setpoint carried on the actin network, coupled to the Card-3 pressure-consistency gate —
is the fine-grained-faithful path. Do NOT add a lumped cortex area-tension, hand-swap ΔP, or loosen the gate
without PI sign-off.

## 7. Reproduce
```
# gbook A5000 (both 1600-iter runs: 1502 s each, concurrent)
PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim ~/miniconda3/envs/ffn_sim/bin/python \
  scratch_ac_pretension_solve.py --n-filaments 70686 --subdiv 6 --pretension {none,forcedensity} --radial \
  --solver erm_gauss_seidel_tournament --n-inner 1600
# structure probe:  scratch_ac_cortex_probe.py --n-filaments 70686 --subdiv 6
# CPU oracle (dev Mac, no CUDA):  formfound_cpu_oracle.py  (area-gradient vs edge-spring vs sparse)
```
Data: the TL;DR table + §2 build-only decomposition are full-native measurements; scratch probes are a
diagnostic harness (not committed). No code path was changed: the force-density pretension was applied as a
one-off host construction on the built `r0xl` array, and — being refuted — is NOT wired into the production
driver.
