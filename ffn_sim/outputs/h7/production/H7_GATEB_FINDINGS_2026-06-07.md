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

## What IS solid (deliverables)
- 3-channel γ estimator with the turgor-separation discipline (cortical_tension.py).
- The full physiological cell assembles + runs at full ×40 scale (cortex+xlink+myosin+nucleus+
  turgor+cytoplasm+membrane+FA), GPU.
- Two single-cell lamellipodium spreading geometries + the A/A0 comparison + cell-movement
  animation (h7_spreading_movement.gif): basal_ring isotropic +24%, polarized_patch directional.
- Stability map: unconstrained-soft is robust (5/5); constrained+FA is integrin-overload-fragile
  (fix spec'd); free-cortex constrained is stable (3/3).
