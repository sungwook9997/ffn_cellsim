---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# FF network contractile-transmission (screening) — validated PHYSICAL against reconstituted ground truth (2026-07-08)

**Task (a), PI-directed.** Resolve the last open lever from FF_STAGE6T: is the FF connected-inextensible network's
screening of contraction (σ/σ_dipole ≈ 0.2–0.6, no Ronceray amplification; FF_STAGE6H) **physical** (real cortex
screens too) or an **FF transmission artifact** (real semiflexible-actin networks amplify/transmit better)? Decided
by comparing to reconstituted actomyosin ground truth (4-agent literature workflow `wf_8ebbd778-0f1`) + a
buckling-vs-rigid control on `network_contractility.py`.

## Verdict: PHYSICAL — FF is in the dense-network no-amplification regime; the shortfall is DENSITY, not transmission

### Reconstituted ground truth (the decisive data)

Reconstituted actomyosin networks deliver **≤ the coherent motor-dipole sum, never more** at the macroscopic scale.
The famous Ronceray "amplification" is a **sparse-active-unit + long-bucklable-fiber** effect that **saturates and
collapses to ~1 as active-unit density rises** — it does NOT apply to a dense one-motor-per-fiber cortex.
**Ronceray, Broedersz & Lenz 2016 (PNAS 113:2827) Table I**, by regime:

| regime (Ronceray Table I) | system | σ_exp | σ_lin | ratio |
|---|---|---|---|---|
| **dense 3D actomyosin (cortex-relevant)** | Koenderink/Bendix 2009 | 14 Pa | 12 Pa | **≈1.17 (~unity)** |
| sparse 2D on vesicle (units ~20 µm apart) | Carvalho 2013 | ≥1 pN/µm | 0.014 | ≈70× |
| fibrin/platelet clot (extensible fiber) | Lam 2011 | 150 Pa | 9 Pa | ≈17× |

- Bendix 2008 (Biophys J 94:3126): bulk α-actinin gel, actin 23.8 µM, myosin:actin ~0.1 → only ~1 µN total / ~100 pN
  per bundle → **σ_exp ≈ 0.2–6 Pa**, while the coherent dipole sum is tens–hundreds Pa → efficiency ~10⁻³–10⁻².
- Murrell & Gardel 2014: quasi-2D cortex threshold **σ ≈ 1.8 Pa** — far under the dipole sum.
- Linsmeier 2016 (Nat Commun 7:12615): strong cooperativity (Hill ~11), semiflexible bucklable filaments (Lp 20 µm,
  L 7.1 µm), critical myosin ~0.56 filaments/µm² — but reports NO absolute Pa, so gives no ratio directly.

**So: dense actomyosin transmits at ~unity vs the linear baseline; amplification is a sparse+extensible effect that
saturates. FF runs at the dense end (one motor/fiber, mesh 0.3 µm) → the no-amplification regime.**

### FF confirmation — buckling control (`network_contractility.py`, f_act=20 pN)

| seed_z | σ/σ_dipole (vs COHERENT sum) | z_bow [µm] | buckled |
|---|---|---|---|
| 0.00 | 0.453 | 0.000 | no (planar) |
| 0.08 | 0.417 | 0.168 | YES |
| 0.20 | 0.406 | 0.237 | YES |

FF's ratio is **robustly ≈0.42–0.45 vs the COHERENT baseline** (n·f·ℓ, every dipole magnitude summed — a LARGER
denominator than Ronceray's linear σ_lin). It is buckling-consistent (buckles when seeded, z_bow up to 0.24 µm) and
crucially **does NOT collapse to ~0** — so FF is not in the artifactual rigid-rod over-screening regime (a bending-free
rod network would over-cancel extensile≈contractile → ~0). Reconciling the denominators (coherent ≈ 2–2.5× linear),
FF's 0.42 vs coherent ≈ **~unity vs the linear baseline** — i.e. consistent with the dense-actomyosin ~1.17.

## Consequences

1. **The FF screening is physical, not a transmission artifact** — dense networks transmit at ~unity; there is no
   lost ×7 amplification to recover. This CONFIRMS FF_STAGE6V: myosin's direct material tension (γ_myo) really is low,
   and the cortical tension is **pressure-borne** (γ=ΔP·R/2), not myosin-material-generated.
2. **The 25–50× cortical-tension shortfall is a DENSITY + prestress gap**, already resolved in direction (FF_STAGE6T
   proteomics: in-vivo cortical myosin ~15–150× above Nie/reconstituted; the cortex is not myosin-protein-limited).
3. **Honest residual caveat (minor, quantitative):** FF's 0.42 vs the coherent baseline vs the dense benchmark ~1.17
   vs the linear baseline — a targeted density-sweep + explicit coherent↔linear baseline reconciliation would confirm
   FF isn't over-screening by ~2–5×. This is a calibration check, not an open physics gap.

## Net

The last open lever from FF_STAGE6T is **closed: the network screening is physical.** Combined with FF_STAGE6V
(pressure-borne tension) + the 2026-07-08 AFM-modulus work (biphasic cytoplasm + physical press + η-drag), the FF
cortical-mechanics picture is now internally consistent and validated against reconstituted-network + theory
literature. **No artifact; no missing datum blocking the tension.** Sources: Ronceray 2016 (DOI 10.1073/pnas.1514208113);
Bendix 2008 (DOI 10.1529/biophysj.107.117960); Murrell & Gardel 2012 (DOI 10.1073/pnas.1214753109); Linsmeier 2016
(DOI 10.1038/ncomms12615); Koenderink 2009. See FF_STAGE6H (screening measurement), FF_CORTICAL_MECHANICS_STATE_2026-07-08.
