# FF ECM library — overnight autonomous execution plan (2026-07-11 night → 2026-07-12)

PI is asleep; full autonomous operation authorized. Cadence: self-paced ~1 h check-ins (session loop, not cloud
— the work monitors session-bound gbook GPU jobs). Standing rules for tonight (per PI + CLAUDE.md + memory
[[project-dcm-8h-autonomous-mandate]]):

- **Never block on a PI decision.** If a step needs a PI call, make the best literature-anchored choice, mark it
  **PROVISIONAL (PI to ratify)** in the doc + STATUS, and continue. Hard rules = document + continue, NOT halt.
- **No tuning to an outcome.** Every constant derived / literature-anchored; validation bands are acceptance
  oracles, never fit targets. Surface any magic-number trigger, don't silently tune.
- **NATIVE + full compartments is the only authoritative basis.** Coarse = labelled dev smoke; every finding is
  reconfirmed native (gbook A5000) before it is reported as a result.
- **Always compare against the KB/literature + verify figures visually** (real render read-back), per the session
  standing instructions.
- **Commit each landed unit** with honest framing; refresh STATUS + figures at each milestone.

## GPU pipeline (serialized on the single A5000 — one job at a time)

Chained on gbook so the GPU keeps working between my check-ins (a waiter launches the next job when the prior
PID exits). Log per job under `~/ff_scratch/*.log`; I poll by log file each loop.

| # | job | status | authoritative output |
|---|---|---|---|
| A | **Native #2** — Path (a) pa_gel traction + ECM-deformation vs stiffness (nf=38000, 3 E) | RUNNING | `ecm_stiffness_ecmnet.json` (native), the clean ECM-deform curve |
| B | **Native #6** — contact-guidance A_F(S) + R_σ(S) on aligned collagen (5 S × 3 seeds) | QUEUED (after A) | `contact_guidance_cg_native.json` + figure |
| C | **Native #6b** — static traction-REMODELING on aligned collagen (the S6-extension half of NEAR #6): does the cell reorient/recruit aligned fibers? (ecm_force reaction already wired) | TO BUILD then queue | remodel displacement field vs S |
| D | **MID — long-range stress propagation** |σ|(r)~r^−n from a contractile inclusion (KB-1.10; needs the #3 virial tensor, DONE) | TO BUILD | σ(r) log-log slope vs lit n≈1 |

## Work queue (each item: build → coarse smoke (verify plumbing) → native → validate vs KB → commit → figure)

### Tonight (autonomous, non-gated)
1. **Land Native #2** when job A finishes: fold the authoritative pa_gel curve into README/STATUS, mark ROADMAP
   NEAR #2 ✅ DONE (native-confirmed), commit.
2. **Land Native #6** when job B finishes: the A_F(S) + R_σ(S) result vs the Ray >3× / Szulczewski 35× bands.
   Expected: R_σ monotone with S (grid-invariant, robust); A_F resolves above the few-clutch noise with the
   3-seed ensemble → report whichever the native data supports (traction anisotropy present or matrix-stress-only),
   honestly. Mark ROADMAP NEAR #6 accordingly, commit + figure.
3. **Build + run Job C (traction-remodeling on aligned collagen)** — the remodeling half of NEAR #6. The cell's
   `ecm_force` reaction already remodels the matrix (validated in S6); measure the fiber reorientation / recruitment
   as a function of the pre-existing alignment S (contact guidance ⊗ remodeling). Coarse smoke → native → commit.
4. **Build + run Job D (long-range stress propagation)** — a contractile inclusion (a small contracting cell or a
   point contraction) in the collagen library; measure |σ|(r) via `ecm_material_stress` on radial shells; fit the
   decay exponent vs the literature |σ|~r^−1 (KB-1.10, 3D elastic). Oracle-overlay only. Coarse → native → commit.

### When the queue drains — WRITE THE NEXT PLAN (per PI) and continue
Draft the next roadmap tranche and keep executing. Candidate next items (autonomous unless marked PI-gated):
- **Per-material native full-extent annotated HTML viewers** (NEAR #4 remaining half) — interactive 3D morphology
  viewers for each of the 6 materials at native extent (per the viz-integrity + no-downsample rules).
- **Kim BD acceptance-oracle cross-validation** for collagen G, K(γ), N1 (oracle-only, per PI-2026-06-30). Build
  the Cytosim/Kim-BD oracle comparison as a validation cross-check (NOT a runtime mechanism, NOT KB-registered).
- **ECMRemodeler topological bond-mutation operator** (MID) — emergent plasticity/MMP/LOX. *Needs the matrix Bell-
  rate datum → PROVISIONAL: use the α-actinin k_off0=0.066 s⁻¹ (Ferrer 2008, already in KB) as the transient-
  crosslink rate, covalent = k_off→0 limit; document as PI-to-ratify.*
- **⟨z⟩(c) LOX crosslink-density law → predicted c²** (PI DECISION). *PROVISIONAL: implement a literature-anchored
  ⟨z⟩(c) growth (LOX density ∝ c, capped above the rigidity threshold) as an OPTION (default OFF), show it yields
  c≈2, mark PI-to-ratify — never an n-fit.*
- **ν-faithful continuum gel** (PI-gated) — add a bending/multibody term to break the Cauchy ν=1/4. *PROVISIONAL
  option, default OFF, documented.*

## KB actions (out-of-band, PI-gated apply — draft only tonight)
- Verify the NEAR #6 validation-band DOIs' `source_audit` verdict (Ray 2017, Szulczewski 2021, Riching 2014, Han
  2016, Niraula 2025) before they land in a deliverable; register missing SourceEvidence rows as DRAFT.
- Draft (do NOT auto-create) the NEAR #6 contact-guidance KnowledgeClaim: FF-measured A_F(S)/R_σ(S) as the KB
  number, Ray/Szulczewski as acceptance oracles (per PI-2026-06-30 oracle-is-crosscheck).
- Run `make kb-check` at each commit; halt + surface only if a gate reports DRIFT.

## Loop protocol tonight
Each ~1 h wake: (1) poll the gbook GPU pipeline logs; (2) if a native job finished → validate vs KB + commit +
figure + advance the queue; (3) if a job hung/failed → diagnose (GPU util + log age), relaunch or fall back to the
next job, document; (4) if the queue drained → build the next item or extend this plan; (5) re-arm the 1 h loop.
Never idle-wait passively — always advance a buildable no-GPU item (KB verification, next driver, docs) while the
GPU is busy. Surface nothing to the PI unless a hard-rule DRIFT or an unrecoverable blocker appears.

---

## MID execution tranche (2026-07-12, autonomous — appended per PI "계획서 계속 확장")

All 6 NEAR done + native. Now executing MID. GPU pipeline (A5000) still serialized; each item follows
build → coarse smoke → native → validate-vs-KB → figure + interactive HTML → commit. Never block on a PI
decision (choose the best literature-anchored option, mark PROVISIONAL, continue).

### In flight
- **Job D — long-range stress propagation** (KB-1.10): native 3-seed running. FINDING (single-seed native, 3-seed
  refining): aligned collagen (S≥0.59) channels contractile stress ~10× farther than isotropic; exponent n drops
  6→4.4 with alignment. Absolute n≈4-6 is STEEPER than elastic (r^-2..-3) — the athermal sub-isostatic Mikado
  localizes stress (short-range), the SAME stretch-dominated limitation as the c-scaling + N1 findings; the true
  fibrous long-range (Notbohm n≈1) needs the PI-gated nonlinear/bending physics. TO DO: land 3-seed n(S)±sd +
  stress-field interactive HTML (deformation around the inclusion).
- **Fibrin contact-guidance generalization** (extends NEAR #6): native run in flight — does R_σ(S) rise like
  collagen? (early: A_F~1, R_σ low at S=0; awaiting S≥0.59).

### Next MID items (priority order, autonomous)
1. **#6 active-motility follow-up** — the honest next step from the NEAR #6 finding. A POLARIZED / actively-
   contracting cell (myosin-biased or a protrusion along +x) on aligned collagen: does the ACTIVE traction A_F
   exceed 1 (Ray >3×) when the resting cell's did not? Reuses the crawl/motility layer + the ecm-network coupling.
   *This converts the R_σ matrix signal into the biological traction guidance — the key open question.*
2. **Interactive stress-field + deformation viewers** for Job D (per PI viz directive) — the collagen network
   colored by |displacement| around the contractile inclusion, isotropic vs aligned (the channeling made visible).
3. **Long-range r^-n with a REAL contractile cell** (not a prescribed inclusion) — embed the FF cell (cortex+
   myosin) in collagen, let its own contraction load the matrix, measure σ(r). Ties Job D to the cell engine.
4. **Kim BD acceptance-oracle cross-validation** for collagen G / K(γ) / N1 (oracle-only, per PI-2026-06-30) —
   a Cytosim/Kim-BD comparison as an independent cross-check of the emergent moduli. NOT a runtime mechanism.
5. **ECMRemodeler topological bond-mutation** (MMP/LOX plasticity) — PROVISIONAL rate = α-actinin k_off0=0.066/s
   (Ferrer 2008, in KB), covalent = k_off→0; PI-to-ratify.
6. **⟨z⟩(c) LOX law → predicted c²** — PROVISIONAL literature-anchored growth option (default OFF), show it yields
   c≈2, never an n-fit; PI-to-ratify. (Closes the documented c-scaling gap as an OPTION.)

### PI decisions parked (auto-decided PROVISIONAL, listed for morning ratification)
- collagen concentration for sweeps = spec.ref_conc 1.5 mg/mL (physiological; not tuned).
- sensing cell MICROTUBULES-OFF (consistent with the NEAR #1/#2/#6 sensing family).
- new KnowledgeClaims (contact-guidance A_F/R_σ; stress-propagation n(S)) DRAFTED, not auto-created (PI-authored).

### Finding (2026-07-12): the #6 ACTIVE-motility follow-up is BLOCKED (documented, not abandoned)
The crawl driver (`ff_crawl_on_substrate.py`) ALREADY supports a polarized cell on aligned ECM (`ecm_align_s`, `phat`, `myo_rear_bias`) — but a NATIVE crawling cell is ~500× slow (the grid-drag Σγ∝Nc blocker, `--bulk-drag` unstable = OPEN; memory [[project-ff-crawl-mechanism]]). So the full active-crawl-on-aligned-collagen test is gated by the SAME Stokes-drag migration unblocker that gates the FAR vision — it is NOT a quick MID item. Deferred to FAR (post-unblocker). INSTEAD generalized Job D (stress propagation) to FIBRIN to confirm the aligned-channeling is a general fibrous-network property.
