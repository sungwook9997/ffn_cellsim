# H.7 — σ_a / network-amplification investigation (2026-06-07)

A 4-agent Workflow (engagement-throttle code diagnosis + network-amplification literature +
realized-network-σ_a measurement design → synthesis) testing the PI-chosen hypothesis: does the
independent-dipole σ_a formula UNDERCOUNT the real percolated cortex's active stress (~1000 Pa)?

## Headline: the formula-undercount hypothesis is REFUTED — the gap is GENERATION, not the formula
γ_soft (method-of-planes, cortical_tension.py:213-379) sums the REAL per-bond tension T=k·(L−r0) of
every cortical bond crossing the cut-plane — it reads the ACTUAL deformed-network state. If a
contracting motor stretches a crosslink three filaments away, that tension is already counted. So
**γ_soft is the realized network-transmitted tension, NOT a sum of isolated dipoles** — it already
contains any percolation/amplification the connected mesh produces. The M1 independent-dipole
formula (σ_a ≈ 100 Pa) is an UPPER BOUND on what isolated engaged motors could inject; γ_soft is the
REALIZED value. Decisive fact: realized γ_soft ~1e-4 mN/m is ~1000× BELOW even the dipole ceiling
(0.10). A measure that already captures amplification, reading 3 orders UNDER the unamplified
estimate, cannot be suffering "formula undercount". The deficit is UPSTREAM of transmission.

Two serial GENERATION causes (both confirmed by prior diagnostics):
- **(a) engagement** — too few heads bind in short runs (~170/5600). KINETIC, not geometric: p_bind
  = 1−exp(−k_on·batch_dt) caught mid-transient (equilibrium needs ~20 ms; runs don't reach it). The
  k_off0 re-anchor 10→0.35/s (53ca7df) makes naive equilibrium 0.993, so binding is no longer the
  RATE limiter — the lever is kinetic pre-equilibration + run-length, NOT loosening any gate.
- **(b) contraction development** — s_grip ≈ 0: heads bind + load to ~0.2-0.44 F_stall but DON'T
  WALK (s_grip 0.012 nm vs ℓ₀ 500 nm), so bond extensions never grow → T stays tiny. The
  "FORCE-AGGREGATION" wall (track1_gsoft_verify_VERDICT.json: γ_soft decoupled from bound frac,
  corr=−0.10, saturates at binding onset).

## Engagement throttle — gates are PHYSICAL, only kinetics is the (legit) lever
- 824 nm head-to-actin placement: STALE/already FIXED (current built state median 200-232 nm,
  actin-aware path myosin.py:509-550). REFUTED — do not re-investigate.
- gate4 bipolar pairing (~50% cut, myosin.py:938-978): PHYSICAL (Lenz 2012/Stam 2017/Murrell 2012
  antiparallel-contractility) — KEEP.
- capture_perp 210nm, MAX_HEADS_PER_BEAD=3, degree≤6: PHYSICAL (head reach; HOOMD 7-exclusion nlist
  limit); gate5_degree_blocked=0 → not shedding heads. KEEP.
- **Gate-chasing risk flagged HARD**: widening capture/degree/bipolar to push engagement adds
  non-contractile pairs without real stress — forbidden. Only kinetic pre-equilibration is legit.
- **LATENT BUG (myosin.py:544)**: `fil_idx = bead_choice // beads_per_filament` assumes uniform
  bead-blocks — WRONG for the variable-length bimodal cortex (cortex.py:1024). Dormant under the
  current uniform fast-hybrid path; fix before the faithful bimodal path feeds myosin.

## Network amplification IS real (~10×) but GATED by buckling — gate B behind gate A
Literature (verified): a crosslinked/percolated network produces macroscopic active stress FAR
above the dipole sum, via (1) a percolation/connectivity THRESHOLD (Alvarado/Koenderink 2017;
Bendix 2008; Ennomani 2016) and (2) nonlinear-fiber AMPLIFICATION via buckling (Ronceray-Broedersz-
Lenz 2016 PNAS "Fiber networks amplify active stress", ~10×; Murrell-Gardel 2012). So 100 Pa →
~1000 Pa is quantitatively plausible. **CRITICAL CAVEAT**: BOTH require the network to actually
BUCKLE/CONDENSE above threshold. The constrained M-SHAKE backbone STRUCTURALLY FORBIDS buckling
(r/r0=1.000; per-filament load ~2 pN < F_Euler ~6.9 pN; ACTIN_ARCHITECTURE_NOTES §16,§★). So gate B
(buckling) sits behind gate A (generation); even a fully-engaged walking network may floor at ~the
dipole ceiling unless M-SHAKE is relaxed (a PI integrator-freeze sign-off).

## The decisive experiment (the plan)
1. [cheap] confirm k_off0=0.35 loaded + D2 CFL now permits batch_dt ~2.8e-3 s (×2000 — makes a long
   contraction run affordable); re-run h3_ku35_estimator_audit synthetic-shell check (trust the estimators).
2. [code] construction-time binding PRE-EQUILIBRATION (seed bound heads at k_on/(k_on+k_off0)
   occupancy among bipolar-eligible, force-free) — the legit engagement lever (physiological-baseline rule).
3. [code] promote the IK whole-shell virial as a 4th channel γ_soft_ik (h3_ku35_estimator_audit.py:89-95;
   ~12× lower variance, validated <0.5% vs MOP) + HOOMD pressure_tensor cross-check + σ_αβ tensor.
4. [run] extend track1_gsoft_verify.py: connected_mesh + pre-equilibrated binding → run to s_grip→ℓ₀
   (≳5-10 s sim, Hill F/F_stall≤1 filter); plot γ + s_grip vs time. CONFIRM = γ climbs toward the band
   as s_grip→ℓ₀; REFUTE = γ floors at the dipole ceiling with the network driven+organized.
5. [contingent, PI sign-off] if γ floors at the ceiling → test gate B: relax M-SHAKE so filaments
   buckle/condense (success = r/r0<1 AND γ above the dipole ceiling).
6. [hygiene] fix the latent bimodal bug (myosin.py:544); mark stale binding-diagnosis docs superseded.

**Honest bottom line**: likely PARTIALLY closes (γ off 1e-4 → ~0.10 dipole ceiling via the
generation fix) then re-confirms the floor at the ceiling unless gate B (buckling, M-SHAKE relax)
is enabled. Reaching the band needs the ~10× buckling amplification M-SHAKE currently forbids. The
active γ_soft/γ_soft_ik channel at F/F_stall≤1 is the ONLY test — never fold in γ_rigid/γ_passive.

## 2D crosslink connectivity test (PI ask; scripts/crosslink_2d_test.py, fig crosslink_2d_test.png)
Validates the connectivity mechanism in a clean 2D plane at ×40 mesoscale: at ~1880 nm filament
spacing the DYNAMIC 60 nm crosslinker reach gives z=0.20, giant 1.5% (181 components) → FRAGMENTED;
the BRIDGE √(A/n) reach (the connected-mesh recipe) gives z=9.4, giant 99% (2 components) →
PERCOLATED. ⇒ the ×40 coarse-graining puts filaments far beyond the physical 60 nm crosslinker
reach, so the connected-mesh long-reach bridge seeding is a coarse-graining-compensation device (the
crosslinker analogue of the mesoscale myosin force-scaling). The dynamic short-reach binding CANNOT
percolate the coarse-grained cortex — this is why connected_mesh is required.
