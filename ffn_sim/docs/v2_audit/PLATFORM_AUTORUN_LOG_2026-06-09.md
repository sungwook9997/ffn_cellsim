# Platform Autorun Log — 2026-06-09

Autonomous compartment-platform deepening session (CWD `ffn_cellsim-platform`,
branch `h7/compartment-platform`). Sibling session works Gate-B cortical-γ in
`ffn_cellsim` (`h7/full-cell-integration`) — its files/branch untouched here.

**Mission:** safely deepen the 8 default-OFF compartments — only work possible
WITHOUT PI parameter ratification or shared-file edits. Smoke harnesses are
PLUMBING checks (never band pass/fail). None-gated constants use DETAILED_PLANS
candidate values, labelled SMOKE-ONLY. Shared/PI-gated work → PLATFORM_PI_QUEUE.md.

| time | item | commit | result |
|---|---|---|---|
| boot | read AGENTS/CLAUDE + 4 boot docs + registry + detailed plans; env ffn_sim, hoomd 7.0.1; clean tree @ fd00309 | — | OK |
| 02:00 | smoke infra: scripts/compartment_smoke/_smoke_common.py (BAOAB stepping sim + cortex stub + SMOKE watermark/fig/json helpers) | (this commit) | OK |
| 02:00 | smoke#1 microtubules: enabled aster+stub, BAOAB step ×3000 @ dt=0.5·dt_cfl; backbone len=l0 (rms dev 0.03%), rod stable, finite | acfa577 | PLUMBING_OK |
| 02:10 | smoke infra: +to_hoomd_snapshot round-trip (bare-gsd None-group hardening for non-gsd-hardened extenders) | (this commit) | OK |
| 02:10 | smoke#2 intermediate_filaments: linear cage (if_backbone + 13 per-r0-bin if_crosslink on one shared Harmonic) ×3000; construction 53.7 kT (force-free ✓ = r0=0 audit-fix holds live), rms dev 0.34% | 172fcbc | PLUMBING_OK |
| 02:20 | smoke#3 linc: perinuclear-cap topology (acceptors seeded at r0 → force-free), 80 linc_nesprin bridges (no zero-bonds trap), construction 0.0 kT, len=r0=50nm, rms tension 4.2 pN (thermal); k_linc=1e-2 SMOKE-ONLY. capture<R_cell-R_nuc warning EXPECTED for cap topology | 98f1b07 | PLUMBING_OK |
| 02:35 | smoke#4 osmotic_regulation: live EnclosedVolumePressure (turgor 133 Pa) + setpoint updater (object-identity), 60 ticks, reads ev.last_pressure + mutates ev.p.V0; hyperosmotic Δc=-100 → V0 RVD-sign down (-0.02%, run 2.1ms ≪ τ_RVD=3.2s honest slope-only); Lp=1e-12 SMOKE-ONLY | f5548b8 | PLUMBING_OK |
| 02:55 | ⭐BUG FOUND+FIXED: cadherin.extend_snapshot_with_cadherins did NOT carry the per-particle `image` array (handled mass/charge/diameter/velocity, omitted image) → Snapshot.from_gsd_frame crash (broadcast (n,3)→(N,3)) on any host frame with image. Fixed in junction/cadherin.py + regression test test_snapshot_extender_extends_velocity_and_image_arrays (32 pass). First smoke catch. | (this commit) | FIX |
| 02:55 | smoke#5 cadherin_junction: 2-cell interface (A/B facing cadherins), Rakshit catch-bond binder forms 10 trans-dimers ALL A↔B (0 intra-cell ✓), elastic force finite; membrane-anchored drag SMOKE-ONLY (cadherins not free-3D-diffusing); n_cad=40 layout not Iturri-223. NOT GATE-J | (this commit) | PLUMBING_OK |
