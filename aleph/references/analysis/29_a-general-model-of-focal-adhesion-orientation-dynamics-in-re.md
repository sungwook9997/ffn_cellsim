---
id: 29_a-general-model-of-focal-adhesion-orientation-dynamics-in-re
paper_n: 29
title: "A general model of focal adhesion orientation dynamics in response to static and cyclic stretch"
authors: "Rumi De (single author)"
year: "2018"
venue: "Communications Biology (Nature)"
doi: "10.1038/s42003-018-0084-9"
paper_type: continuum-model
ffn_relevance: Medium
ffn_themes: [FA/clutch, validation-oracle, parameter-source, ECM/collagen]
entities: [focal-adhesion, integrin, ligand-receptor-bond, actin, stress-fiber, catch-bond, slip-bond, substrate, extracellular-matrix]
methods: [monte-carlo, master-equation, gillespie-algorithm, stochastic-simulation]
measurables: [bond-lifetime, adhesion-cluster-size, bond-off-rate, bond-on-rate, fa-orientation-angle, threshold-strain]
keywords: [focal-adhesion, catch-bond, slip-bond, master-equation, gillespie, cyclic-stretch, mechanosensing, adhesion-cluster-stability, two-pathway-model, ligand-receptor-rebinding, fa-reorientation]
tags: ["#fa-clutch", "#catch-bond", "#adhesion-cluster", "#validation-oracle", "#parameter-source", "#mechanotransduction", "#stochastic-model"]
has_transferable_params: true
---

# [29] A general model of focal adhesion orientation dynamics in response to static and cyclic stretch

**Tags:** #fa-clutch #catch-bond #adhesion-cluster #validation-oracle #parameter-source #mechanotransduction #stochastic-model

| Field | Value |
|---|---|
| Authors | Rumi De (IISER Kolkata, single author) |
| Year / Venue | 2018 / Communications Biology (Nature) |
| DOI / ID | 10.1038/s42003-018-0084-9 |
| Type | continuum-model (stochastic master-equation / Gillespie Monte-Carlo of an adhesion bond cluster) |
| Pages | 7 |
| ffn_cellsim relevance | Medium — catch-slip bond cluster oracle for the FA/clutch unit; uses the same Pereverzev two-pathway off-rate ffn_cellsim adopts, but it is a lumped cluster (not particle-resolved) and the headline question (whole-cell/FA reorientation under stretch) is outside the single-cell mechanistic scope |

## 1. Summary
De presents a minimal stochastic theoretical model of a focal-adhesion (FA) cluster as a parallel array of Hookean ligand-receptor springs (stiffness k_b) coupling an actin stress fiber (via a cellular spring k_c) to a stretched substrate. Bond dissociation follows a force-dependent catch-slip two-pathway rate (Pereverzev-type, k_off = k_slip·exp(+f_b/f0) + k_catch·exp(−f_b/f0)); rebinding (k_on) is made stretch-rate-dependent through a novel ligand-receptor "contact-time" argument: under cyclic strain the ligand on the substrate moves away from the receptor, so k_on = γ(N−n)·exp(−v_ω²/v0²) with v_ω = l_R·ε0·ω. The cluster occupancy n(t) evolves under a one-step master equation, solved by Gillespie Monte-Carlo. The model reproduces, within one unified framework, the puzzling experimental observation that FAs orient **parallel** to a static/quasi-static stretch but **perpendicular** to a fast cyclic stretch, plus optimal-strain stability (catch-bond strengthening up to a threshold then disassembly), cell-type-dependent orientation (via intrinsic rebinding time T_R), and a threshold stretch magnitude/frequency for reorientation.

## 2. Problem & motivation
FAs are the micron-scale integrin clusters that mechanically link the cell cytoskeleton to the ECM and are central to mechanosensing. Experiments show FAs reorient differently under static vs. cyclic substrate stretch (parallel vs. perpendicular), and prior theories addressed only the perpendicular/cyclic case or the whole-cell/stress-fiber response, leaving no single theory explaining both regimes. The goal is a unified, mechanistically-motivated stochastic model that predicts FA growth, stability, and orientation across the static→cyclic-stretch spectrum from bond-level kinetics alone.

## 3. Methods / model
- **Model class**: coarse-grained stochastic adhesion-cluster model (NOT particle-resolved). One FA = N parallel Hookean ligand-receptor bonds (stiffness k_b) in series with a single cellular spring k_c representing the stress fiber/cell.
- **Elasticity / force balance**: geometric constraint u_b + u_c = Lε (L = FA distance from cell center, ε = strain). Force balance k_c·u_c = n·k_b·u_b. Single-bond force f_b = k_b·u_b = Lε·k_b·k_c/(k_c + n·k_b). Effective strain on a cluster at angle θ: ε_a = ε0·cos²θ.
- **Dissociation (catch-slip, Eq 2)**: k_off = k_slip·exp(+f_b/f0) + k_catch·exp(−f_b/f0) — slip term promoted by force, catch term opposed by force; cited from Pereverzev two-pathway model (ref 39) and Thomas catch-bond review (ref 15).
- **Association (Eq 3)**: k_on = γ(N−n)·exp(−v_ω²/v0²), with intrinsic binding rate v0 = l_R/t_R, reaction radius l_R, contact time t_R, and cyclic displacement rate v_ω = l_R·ε0·ω. Static limit (ω=0): k_on = γ(N−n).
- **Time evolution (Eq 1)**: one-step master equation dP_n/dt = K_on·P_{n−1} + K_off·P_{n+1} − (K_on + K_off)·P_n, with K_off = n·k_off.
- **Numerics**: all parameters non-dimensionalized (τ = k0·t with k0 = spontaneous off-rate; K_s = k_slip/k0, K_c = k_catch/k0, Γ = γ/k0, T_R = k0·t_R). Solved by Monte-Carlo / Gillespie's exact stochastic algorithm (ref 42). Cluster starts all-closed, runs ≥1 million events per stable trajectory; statistics over 500 trajectories.
- **Scales**: cluster size N = 200 bonds; FA distance L = 20 μm; reaction radius l_R ~ 1 nm; single-bond/cell stiffness ratio k_b/k_c = 5; f0 ~ piconewton scale.

## 4. Key results (quantitative)
- Static stretch: cluster grows with strain, peaks at an optimal value, becomes **unstable above ~20% strain** (Fig 2a: 1% and 10% stable, 20% unstable; N starts at 200). Mean closed-bond count rises with Γ (Fig 2b: Γ=1 vs Γ=2). (p.4)
- Parallel direction (θ=0, ε_a = ε0) is most stable under static stretch; perpendicular (θ=90°, ε_a=0) least stable → FAs align parallel to static stretch (Fig 2, p.4-5).
- Cyclic fast stretch (T_ω ≪ T_R): cluster stability maximal near θ=90° (perpendicular); cluster grows toward perpendicular orientation; stronger effect at higher ω (Fig 3a: ω=10 vs ω=5; ~110–140 bonds vs orientation). (p.5)
- Frequency sweep at 10% strain, T_R=10: as ω drops from 1 → 0.01, the most-stable orientation shifts from perpendicular back toward parallel (Fig 3b). (p.5)
- Intrinsic rebinding time T_R controls cell-type behavior: smaller T_R → stable over a wide range of θ (less reorientation); T_R = 10,5,2,1 sweep at ω=10, 10% strain (Fig 4). (p.5)
- Threshold stretch: reorientation only above a threshold ε0; threshold shifts to lower ε0 with increasing ω (since T_ω ∝ 1/ω) (Fig 5: ω = 1–10 sweep, T_R=10). (p.5-6)
- Fixed scaled parameters used: Γ = 2, K_s = 0.10, K_c = 120 (following ref 39 Pereverzev), k_b/k_c = 5 (p.4).

## 5. Parameters & constants of interest
| Quantity | Value + units | Source in paper |
|---|---|---|
| Catch-slip off-rate form | k_off = k_slip·e^{+f_b/f0} + k_catch·e^{−f_b/f0} | Eq 2, p.3 (cites Pereverzev ref 39, Thomas ref 15) |
| Molecular force scale f0 | ~ piconewton (order of magnitude) | p.3 |
| Scaled catch rate K_c = k_catch/k0 | 120 (following Pereverzev ref 39) | p.4 |
| Scaled slip rate K_s = k_slip/k0 | 0.10 | p.4 |
| Scaled binding rate Γ = γ/k0 | 2 | p.4 |
| Cluster size N | 200 bonds | p.4 |
| FA distance from cell center L | 20 μm | p.4 |
| Reaction radius l_R | ~1 nm | p.4 |
| Bond/cell stiffness ratio k_b/k_c | 5 | p.4 |
| Static instability threshold strain | ~20% (cluster disassembles) | Fig 2a, p.4 |
| Single-bond force | f_b = Lε·k_b·k_c/(k_c + n·k_b) | p.3 |
| Effective strain at angle θ | ε_a = ε0·cos²θ | p.4 |
| Association (rate-dep.) | k_on = γ(N−n)·e^{−v_ω²/v0²}, v_ω = l_R·ε0·ω | Eq 3, p.3 |
| Master equation | one-step birth-death, Gillespie MC, ≥1e6 events, 500 traj | Eq 1, p.3-4 |

## 6. Relevance to ffn_cellsim  <-- MOST IMPORTANT
This paper sits squarely on the **FA/clutch theme** (the bridge/junction units; H.4 FA-clutch track) and is useful in two concrete ways:

(a) **Validation oracle** — ffn_cellsim's FA/clutch will be built from explicit integrin clutches with **Bell-Evans/catch-slip force-dependent off-rates**, and this model gives a closed-form *adhesion-cluster stability* result (mean closed-bond count vs. strain; optimal strain then disassembly; threshold magnitude/frequency for reorientation). Those are exactly the kind of emergent cluster-level observables ffn_cellsim could regress its particle-resolved clutch ensemble against in a validation test (cluster grows with tensile load up to an optimum, then catastrophically disassembles). The one-step master-equation / Erdmann-Schwarz (ref 25) "stability of adhesion clusters under constant force" lineage this paper extends is a canonical clutch oracle.

(b) **Parameter source / mechanism reference** — it uses the **Pereverzev two-pathway catch-slip model (ref 39)** as its off-rate, the same family ffn_cellsim's CLAUDE.md lists (Pereverzev is named alongside Bell-Evans/Hill as a v1 oracle). The non-dimensional rate constants (K_c=120, K_s=0.10, Γ=2, f0~pN, N=200, l_R~1 nm) are concrete anchors for sanity-checking a clutch parameterization, and the catch+slip exponential split is directly portable as an acceptance oracle (never as runtime — runtime is the particle/bond dynamics).

**Caveats on scope**: The *headline contribution* — FA **orientation** (parallel vs. perpendicular) under static vs. cyclic substrate stretch — is a substrate-stretch / whole-cell-reorientation phenomenon that ffn_cellsim does not currently model (no externally cyclically-stretched substrate; the project's ECM is collagen-I + clutches, validated against MCF7-spheroid traction/cohesion). So the orientation result is **context, not a direct oracle**. The model is also explicitly lumped: one cellular spring stands in for the whole stress-fiber/cytoskeleton and all bonds are assumed to share the same force — the opposite of ffn_cellsim's per-filament/per-clutch fine-graining. Net: **Medium** relevance — a solid catch-slip cluster oracle and parameter cross-check for the FA/clutch unit, but the stretch-reorientation framing is tangential to the single-cell mechanistic simulator.

## 7. Limitations & caveats
- **Mean-field force**: "all bonds in the adhesion cluster experience the same elastic force or deformation" (p.3) — no spatial bond-force heterogeneity except stochastic; ffn_cellsim resolves this explicitly, so the cluster is a coarse projection.
- **Single cellular spring** lumps the entire stress fiber/cytoskeleton viscoelasticity into one k_c; author concedes stress-fiber viscoelasticity and actin-myosin contraction (refs 43,44) are omitted and "outside the present theory" (p.6).
- **No explicit substrate/ECM mechanics, no actin polymerization, no myosin** — purely bond-kinetics + linear-elastic springs.
- **Substrate-stretch experiment** (cyclic strain on the FA's anchoring substrate) is the validation target; this is a different boundary condition than ffn_cellsim's traction-on-compliant-ECM setup.
- Parameters are non-dimensionalized; absolute SI values (k_slip, k_catch in s⁻¹, k_b in pN/nm) are not given — only ratios and scales.

## 8. Key figures / tables
- **Fig 1 (p.2-3)**: schematic — SF adhered via two FAs, ligand-receptor bonds as Hookean springs; FA oriented at angle θ to stretch. Defines the geometry (u_b+u_c=Lε, ε_a=ε0cos²θ).
- **Fig 2 (p.4)**: static stretch — (a) n(t) for 1%/10% (stable) vs 20% (unstable) strain; (b) mean n vs strain, peak-then-fall, larger for Γ=2. The optimal-strain / disassembly oracle.
- **Fig 3 (p.5)**: cyclic stretch — (a) stability vs θ peaks at perpendicular for fast stretch (ω=10,5); (b) frequency sweep ω=1→0.01 shifts most-stable angle back to parallel.
- **Fig 5 (p.6)**: stability vs strain ε0 at θ=0 for ω=1–10 — threshold strain shifts lower with higher ω.

## 9. Notable quotes / citable claims
- "the dissociation rate k_off of the closed bond is proposed to demonstrate the catch behavior as k_off = k_slip·e^{f_b/f0} + k_catch·e^{−f_b/f0}" (Eq 2, p.3) — the catch-slip oracle form.
- "Tensile force, up to an optimal value, is found to strengthen and reinforce the molecular bonds; and bonds' lifetime decreases with further increase in force... These type of force-strengthening bonds are called catch bonds and are believed to play a crucial role in stabilizing the FA cluster." (p.3)
- "the sole consideration of the stretch-dependent association and dissociation processes of adhesion clusters could successfully predict the orientational response of FAs." (Discussion, p.6)
- "it explains the puzzling observations of parallel orientation of focal adhesions under static stretch and nearly perpendicular orientation in response to fast varying cyclic stretch." (Abstract, p.1)
