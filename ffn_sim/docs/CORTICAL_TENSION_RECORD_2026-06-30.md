---
kb_record:
  topic: KU-3.5-cortical-tension
  claim: KB-3.5
  gate: VG-H3-KU35-cortex-tension
  status: authoritative
  primary: true
  authoritative_as_of: 2026-06-30
  aliases: [KU-3.5, KU3.5, cortical tension, cortical-tension floor, g_soft, gamma-floor, active gamma floor]
  supersedes:
    - CORTICAL_TENSION_RECORD_2026-06-04
  conclusion: >-
    The KU-3.5 active cortical-tension (g_soft) floor is CONFIRMED a force-MAGNITUDE
    (generation/aggregation) limit, and the network mechanisms that could have
    rescued it are now RULED OUT. Three independent methods converge: BAOAB-MD
    g_soft, the MD-free FF mechanical solve, and a full network-mechanism ablation
    in the FF engine — all give active actomyosin gamma ~100-2300x under the
    0.35-0.65 mN/m band at the physiological NMIIA stall. The single measured
    cortical NMII minifilament areal density (Nie et al. 2015, ~0.6/um2, HeLa) is
    ~30x below the naive band requirement, so the measured density CONFIRMS the
    floor (it does not open it). The FF buckling-enabled contractile-network test
    (Stage 6h, 2026-06-30) shows the floor is FLAT in connectivity z (2.98-5.29,
    incl. sub-isostatic) and in axial extensibility (EA swept 3 decades incl. actin
    EA, Gittes 1993) and shows NO Ronceray amplification (sigma/sigma_dipole 0.2-0.6,
    screening) across f/F_crit=1.4-145 with buckling active. Reason: Ronceray's ~7x
    amplification is an EXTENSIBLE-network effect; FF/Cytosim hard inextensibility
    structurally cannot exceed it, and even finite EA does not lift it. So NO network
    mechanism (buckling, connectivity/sub-isostaticity, extensibility, and per the
    earlier work binding kinetics + FA-anchoring + compliant backbone) closes the
    gap. The g_rigid ~0.57 mN/m "myosin-independent passive" channel is reclassified
    as a turgor pressure-PARTNER double-book (turgor_dP0=133 Pa is band-implied/tuned,
    no sourced row; Young-Laplace Delta P*R/2 is the pressure's partner, not an
    independent passive cortical tension); the genuine blebbistatin-insensitive
    passive floor is ~0.04 mN/m (~9-12%), and the band is ~50-90% myosin-dependent
    (central ~70%; Fischer-Friedrich 2016 -68%, Tinevez 2009, Warmt 2021). VERDICT:
    the active-gamma floor is a RESULT, not a blocker — cortical-tension MAGNITUDE is
    set by the force-bearing (load-engaged) motor density, which is an EXPERIMENTAL
    gap (imaging counts presence, not engagement). Report gamma as a function of
    motor density (controlled variable), never tuned to the band.
---

# KU-3.5 cortical tension — authoritative record (2026-06-30)

> **⚠️ MAGNITUDE UPDATE 2026-07-01 (native scale — see `FF_STAGE6O_NATIVE_GAMMA_2026-07-01`).** The
> γ_active magnitudes below were measured at the ×40 mesoscale (cortex N=1000); the ×40 was a CPU
> constraint, RETIRED now that FF is GPU-native. γ_active is N-DEPENDENT and ×40 UNDER-reported it ~5×.
> At the native ~38000-filament cortex (A5000): **γ_active ≈ 6.2e-4 mN/m, floor ~530–570× under band**
> (N-converged; re-verified robust to buckling/turnover/crosslink-stiffness at native). The CONCLUSION
> (force-magnitude floor, force-generation/engaged-density limited) is UNCHANGED — the floor is real,
> not a coarse-graining artifact; only the number updates (the gap narrows ~2500×→~530× but persists).

Supersedes `CORTICAL_TENSION_RECORD_2026-06-04` (which left g_soft OPEN and made no ruling on the
network mechanisms). This record CONFIRMS and EXTENDS the 2026-06-04 force-generation diagnosis with
the rulings the earlier chain lacked.

## What is settled (this session, ultracode + Codex, all cross-checked vs project files)

1. **Force-magnitude limit — confirmed by 3 independent methods.** BAOAB-MD g_soft, MD-free FF, and
   the FF network-mechanism ablation all give active actomyosin γ ~100–2300× under band.
2. **Network mechanisms RULED OUT (the new rulings):**
   - buckling — FF buckles (`test_buckling.py`) yet the contractile network σ stays floored;
   - connectivity / sub-isostaticity — σ flat across z=2.98–5.29 (Stage 6h);
   - extensibility — σ flat across axial EA swept 3 decades incl. actin EA (§6h §4b);
   - (earlier) binding kinetics, FA-anchoring, compliant backbone — all shown not to lift g_soft.
   No Ronceray amplification appears (σ/σ_dipole 0.2–0.6 = screening); the ~7× is an EXTENSIBLE-network
   effect that FF/Cytosim hard inextensibility structurally precludes.
3. **The measured density confirms the floor.** Nie 2015 ~0.6 minifil/µm² (HeLa) is ~30× below the
   naive band requirement; using the real datum drives γ deeper under band, not toward it. The
   runtime "3/µm² = Salbreux 2012" was a confirmed misattribution (re-anchored to Nie 2015).
4. **g_rigid / γ_passive reclassified.** turgor_dP0=133 Pa is band-implied/tuned (no sourced row,
   PARAM_AUDIT 2026-06-25); ΔP·R/2 is the pressure's partner, not an independent passive cortical
   tension → the "passive carries the in-band tension" reading is a double-book. Real passive floor
   ~0.04 mN/m; band is ~70% myosin (blebbistatin).

## The one remaining lever (experimental, not simulation)

The **force-bearing (load-engaged) NMII minifilament density per cross-section** — imaging counts
*presence*, not *engagement*; real cells reach band (~70% myosin) so carry ~12–50× more engaged
motors/cross-section than the present-density gives. This datum does not exist; commissioning it
(super-res + engagement readout, MCF7) is the only way to overturn the floor — predicted to confirm.

## Provenance

FF Stages 6d–6h (git da1812a→92b7a1d, 2026-06-30): `docs/v2_audit/FF_STAGE6D_GAMMA_FLOOR_2026-06-29.md`,
`FF_STAGE6H_BUCKLING_NETWORK_2026-06-30.md`, `FF_CYTOSIM_PARITY_2026-06-30.md`,
`H7_CORTICAL_MYOSIN_DENSITY_DATUM_2026-06-07.md`, `H7_SIGMA_A_NETWORK_INVESTIGATION_2026-06-07.md`.
Code: `ff/gamma_floor.py`, `ff/network_contractility.py`, `ff/cytosim_parity.py`, `ff/buckling`/tests.
No magic-number tuning: f_act, density, connectivity z, EA are all controlled variables (swept) or
measured values. **Pending PI ratification of this supersession + promotion into the Notion
Contract-Graph (SoT) authoritative_record.**
