# H.7 GATE-B — consolidated findings (overnight 2026-06-07)

Full physiological MCF7 cell, full ×40 scale, gbook A5000, rounds 1-5. **PROVISIONAL —
this is the honest GATE-B read, NOT a ratified KU-3.5 conclusion (which is a PI call).
No gate-loosening: the band is NOT chased.**

## The three γ channels, measured cleanly (separate, turgor never folded in)

| channel | clean value | n | reading |
|---|---|---|---|
| **γ_active** (myosin, soft MOP, cortical-only) | **0.0001 ± ~0 mN/m** | 8 (5 soft + 3 turnover) | FLOORS at ~0; turnover-ON does NOT raise it |
| **γ_rigid** (free pressurized cortex, no FA, M-SHAKE) | **0.063 ± 0.016 mN/m** | 3 | turgor transmits only ~0.06 to the backbone — BELOW band |
| **γ_passive** (turgor Young-Laplace = ΔP·R/2) | **0.499 mN/m** | — | CIRCULAR identity with the B3 turgor setpoint (133 Pa); not emergent |
| ~~γ_rigid (FA-adhered)~~ | ~~0.391~~ | 1 (unstable) | **CONTAMINATED** — FA integrin r0=0 overload artifact, discard |

band overlay [0.35, 0.65] mN/m = a rounded/de-adhered non-MCF7 proxy (gate-contract review
pending PI).

## Conclusion (honest): no clean emergent channel reaches the band
- **Active**: γ_active ≈ 0 across 8 seeds, constrained AND unconstrained, turnover ON and OFF.
  The active myosin cortical tension floors at ~0 at the accessible operating point.
- **Structural**: the free turgor-pressurized cortex carries only γ_rigid ≈ 0.06 mN/m — turgor
  IMPLIES 0.50 (Young-Laplace) but the mesh transmits only ~1/8 of it to the backbone
  (a TRANSMISSION gap; the rest is balanced by the enclosed-volume / membrane terms).
- **Passive**: γ_passive 0.50 is a setpoint identity (we set turgor = 2·0.50/R), not a measurement.
- The only "in-band" number we ever saw (0.391) was the FA-adhered constrained run, now shown to
  be an **integrin-overload artifact** (the dynamic integrin–ligand bond is wired r0=0 → born at
  50–1500 pN; fix spec'd in H7_FA_INTEGRIN_OVERLOAD_FIX_2026-06-07.md).

⇒ In the full physiological cell, the EMERGENT cortical tension (active + clean structural)
is ≈ 0.06 mN/m — ~6-10× BELOW the band. This robustly reproduces + sharpens the project's
γ-floor reframe in the complete cell.

## Why: a TIMESCALE GAP (the central method finding)
The tension-generating/-remodelling processes are SECONDS-scale — myosin contraction builds over
~seconds; actin turnover τ½ = 10 s (Chugh 2017). The fine-grained constrained-MD step is µs, so a
feasible run (1e3–1e4 steps) reaches only ~ms. The cortex never reaches its active-tension steady
state. Confirmed directly: turnover-ON (the one reviewed mechanism that could un-floor active γ)
left γ_active at 0.0001 — because τ½=10 s ≫ the run. This is not a missing mechanism; it is the
seconds-vs-µs scale separation. (Cross-checked: γ_active ≈ 0 in BOTH constrained and unconstrained
→ real generation/timescale limitation, NOT an M-SHAKE shunt artifact.)

This is exactly the gap the conventional active-gel / active-shell continuum (Salbreux/Jülicher/
Joanny) sidesteps by working at the seconds+ continuum scale directly — the novelty-analysis's
sharpest threat ("emergent γ is a promise, not a result"). The honest status: emergent active γ
is NOT yet a result at accessible fine-grained timescales.

## Path forward (for PI — options, not chased autonomously)
1. **Accelerated dynamics** (raise the myosin v0 / turnover rates) to reach steady state in feasible
   steps — BUT this breaks Hill-validity (must read γ at F/F_stall ≤ 1) and is an accelerated-knob,
   so it needs an explicit validity gate. (The KU-3.5 memory already flags v0_accel.)
2. **Multiscale seam**: run the fine-grained cell to a quasi-steady local state, hand the
   coarse-grained stress to an active-gel continuum for the seconds-scale relaxation (a γ-seam, the
   inverse of the planned Layer-2 seam).
3. **Re-target the gate**: the [0.35,0.65] band is a non-MCF7 rounded-cell proxy; a spread adherent
   MCF7 is traction/stress-fibre dominated — the right observable may be traction, not cortical γ.
4. **De-circularize γ_passive** (independent osmotic datum) or keep labeling it a setpoint.

## ×40 coarse-graining convergence (round-6 — novelty-defending grid study)
On the STABLE paths, across the sanctioned ×40 mesoscopic range:
- free-cortex **γ_rigid**: nf 500 → 0.0593, nf 1000 → 0.0546, nf 2000 → 0.0858 mN/m — grid-stable
  at ~0.06-0.09 (within the n=3 seed spread 0.063±0.016), no divergence ⇒ **scale-robust, well
  below band**. The structural floor is NOT a coarse-graining artifact.
- **γ_soft**: nf 500 → 0.0000, nf 2000 → 0.0001 mN/m ⇒ the active floor is **scale-independent**.
Both channels are mesh-converged: the floor finding holds across resolution, not just at one scale.

## Verification corrections (2026-06-07, after PI component-check) — IMPORTANT
A PI-requested verification of the membrane / crosslinker / measurement found the overnight
"timescale gap" framing was incomplete and three things had to be CONTROLLED:
- **Cortex connectivity (controlled → RULED OUT as the cause).** The GATE-B baseline had used the
  random-anchor crosslinker mesh (fragmented: z~1.3, giant~7%), NOT the percolated connected_mesh
  (z~3.3, giant~99%). Re-test with connected_mesh: γ_rigid 0.04, γ_soft 0.0000 — the floor HOLDS.
  So fragmentation is NOT the floor cause (control clears the confound). connected_mesh is now a
  GATE-B lever (`--connected-mesh`); full-scale confirmation in round-7.
- **Myosin barely engages (a real contributor).** Only ~12 of ~2000 myosin heads bind the cortex
  after 200 steps → almost no active force is applied (binding-throughput limited, a known KU-3.5
  issue). This + the timescale gap (not cortex connectivity) are the operative active-floor causes.
- **The membrane tension is UNMEASURED.** membrane_surface (H.8 Young-Laplace shell, γ_mem=0.10
  mN/m) IS wired (a FORCE, not particles → invisible in the particle viz) but is NOT in the 3 γ
  channels (which are bond-MOP + turgor). So the reported "cortical tension" omits a real ~0.10
  mN/m membrane component. Flagged; to be added to the report as a 4th component.
- **Band re-anchored to a REAL MCF7 datum.** Hosseini 2020 (Adv Sci, AFM/FF lineage): MCF-7
  interphase γ ≈ 0.27 mN/m (IQR 0.18-0.40). The emergent ~0.06 is ~4-5× below this real MCF7
  number (tighter than the rounded-proxy 6-10×). No spread-adherent MCF7 γ exists → the adherent
  operating point's observable may be traction. (docs/v2_audit/H7_MCF7_CORTICAL_TENSION_DATUM.)

Net: the floor is REAL and robust to the cortex-percolation control, but it is the timescale gap
PLUS myosin under-binding (not connectivity), measured against an incomplete (membrane-omitting)
estimator and re-anchored to MCF7 γ≈0.27. The multiscale active-gel seam remains the path, and its
σ_a input must use connected_mesh + a myosin-engagement check.

## Active-gel seam M1 diagnostic (2026-06-07) — the floor is GENERATION-LIMITED, not timescale
The PI-directed M1 diagnostic (cortex/active_gel_seam.py + scripts/h7_active_gel_seam.py;
fig h7_active_gel_seam_M1.png) relaxes the FG active-stress drive over the seconds timescale the
MD can't reach (active-Maxwell, τ=14.4 s) and computes the steady tension three ways:
- γ_ss realized (FG soft) ≈ 0.0000 | capacity (190 engaged heads) ≈ 0.0001 | **CEILING (ALL 2000
  heads at the 0.5 pN stall × 700 nm dipole) ≈ 0.0010 mN/m** — **353× BELOW band**.
- **VERDICT: GENERATION-LIMITED.** Even the absolute ceiling (every head at stall) is 353× short,
  so the seam (which only adds TIME) CANNOT close the floor. This supersedes the overnight
  "timescale gap" framing: the dominant cause is that the myosin force-dipole DENSITY is far too
  small, not that the MD ran out of time.
- **Root cause surfaced: myosin density = 0.141/µm² vs the Salbreux 3/µm² target (21× low)**, plus
  a low per-head stall (0.5 pN vs NMII ~1-4 pN). Fixing both (×21 density, ×~8 force) lifts the
  ceiling to ~0.17 mN/m — approaching the Hosseini MCF7 datum 0.27 — so the generation deficit is
  largely a PARAMETERIZATION shortfall, fixable, OR re-target the observable to traction.
⇒ Per the PI contingency: M1 is low because of GENERATION (not time) → M2 (1-D shell/flow, also a
relaxer) would NOT help; the actionable forks are (i) fix the myosin generation params to literature
(density 0.14→3/µm², stall 0.5→~3 pN) and re-run, or (ii) re-target the gate to traction/spreading.

## What IS solid (deliverables)
- 3-channel γ estimator with the turgor-separation discipline (cortical_tension.py).
- The full physiological cell assembles + runs at full ×40 scale (cortex+xlink+myosin+nucleus+
  turgor+cytoplasm+membrane+FA), GPU.
- Two single-cell lamellipodium spreading geometries + the A/A0 comparison + cell-movement
  animation (h7_spreading_movement.gif): basal_ring isotropic +24%, polarized_patch directional.
- Stability map: unconstrained-soft is robust (5/5); constrained+FA is integrin-overload-fragile
  (fix spec'd); free-cortex constrained is stable (3/3).
