# PI decisions — FF γ-floor closeout + KB/gate hygiene (2026-06-30)

Everything that did NOT need a PI decision is already done (see §0). This doc lays out the items that
DO need a PI call, each with: what it is, why it's blocked on you, the recommended choice, and (where
relevant) the exact change drafted so approval is one step.

## 0. Already done (no PI decision needed) — for the record

- FF engine 6b-iv→6h built + validated; `tests/ff` 61 pass + 3 Cytosim-parity (skip w/o binary).
- γ-floor characterized + CLOSED as a result (3 methods + network-mechanism ablation converge).
- Magic-numbers I introduced this session FIXED: the "3/µm² = Salbreux 2012" mis-citation (→ Nie 2015),
  the γ_passive turgor double-book (reclassified in `gamma_floor.py`).
- KB authoritative record promoted: `CORTICAL_TENSION_RECORD_2026-06-30` supersedes 2026-06-04;
  supersession.py rebuilt → `tag_query` returns the CURRENT verdict; old record back-linked
  (`superseded_by`). kb-check disk gates GREEN.
- Memory: `project-gamma-floor-likely-deficit`, `feedback-oracle-is-crosscheck-not-truth`.

---

## DECISION 1 — gate-track the untracked magic-numbers (turgor_dP0, cortex myosin density)

**What.** `verify_params` only audits constants registered in `outputs/tag_kb/params_manifest.yaml`.
`turgor_dP0 = 133 Pa` is **not in the manifest at all** (`grep -c turgor` = 0) — so a band-implied /
tuned constant (flagged "tuned, no sourced row" in PARAM_AUDIT_SIMUCELL3D_2026-06-25) sits in the
gate's blind spot. Same for the cortex NMII areal density (was mis-cited to Salbreux; the only datum is
Nie 2015 ~0.6/µm²).

**Why it needs you.** Registering them is gate-TIGHTENING (good, catches the magic-number going
forward), but they have no clean source → the gate would go **UNSOURCED (red)**. Resolving a red gate
is a PI call: (a) accept it as a *controlled-variable / tuned-with-justification* declaration, or (b)
source it. AND `turgor_dP0` is a **DCM constant on the shared gate**, while the **/loop dcm session is
live** — reding the shared gate would break its closeout, so this must be coordinated, not done
unilaterally.

**Recommended.** Register both as `declared: tuned-controlled-variable` with the PARAM_AUDIT / this
record as provenance (honest: they are NOT literature-sourced; γ/turgor are swept, not fitted). This
keeps the gate honest without falsely claiming a citation. Coordinate timing with the /loop dcm
session.

**Drafted manifest rows (apply on approval):**
```yaml
  - id: dcm-turgor-dp0
    desc: "Pa MCF7 baseline osmotic turgor — BAND-IMPLIED/TUNED (no sourced row; PARAM_AUDIT 2026-06-25)"
    config: ffn_sim/dcm/geometry.py
    key: ResolvedDCM.turgor_dP0
    value: 133.0
    ku: KU-3.5
    declared: tuned-controlled-variable   # NOT literature-sourced — see CORTICAL_TENSION_RECORD_2026-06-30
  - id: cortex-myosin-areal-density
    desc: "minifilaments/µm² cortical NMII areal density — only direct datum Nie 2015 (HeLa, non-MCF7)"
    config: ffn_sim/ff/gamma_floor.py
    key: NIE2015_DENSITY_UM2
    value: 0.625
    ku: KU-3.5
    citation_key: Nie2015_Cytoskeleton   # needs a SourceEvidence row (Decision 2)
    declared: source-unverified
```

---

## DECISION 2 — promote to the Notion Contract-Graph (SoT) + harvest the FF runs

**What.** The SoT is Notion (the duckdb/obsidian are read layers I CAN regenerate; the SoT I cannot
edit without you). To make the authoritative chain durable beyond this repo's docs, the KU-3.5 verdict
+ the Nie 2015 SourceEvidence row + the FF 6d–6h RunResults need to enter Notion.

**Why it needs you.** Per CLAUDE.md, Contract-Graph edits + `harvest_ops --apply` + new
SourceEvidence/ValidationGate rows are **PI-authored / PI-gated** (never auto-created).

**Recommended.** Approve: (a) register `Nie2015_Cytoskeleton` (PMID 25641802) as a SourceEvidence row
(real paper, verifiable); (b) update the KU-3.5 ModelContract/DecisionLedger with the 2026-06-30
verdict; (c) `harvest_ops --apply` for the FF 6d–6h runs. I'll prepare the harvest manifest (dry-run)
on request.

---

## DECISION 3 — the next FF unit (roadmap)

**What.** The FF engine is built + validated and the γ-floor (its decisive first experiment,
ENGINE.md §4) is CLOSED. The next FF target is a roadmap call I should not make unilaterally.

**Options.**
- **(3a) Layer-2 single-cell coupling** — route the magnitude question to the fine-grained line
  (consistent with the `layer2-spheroid-line` memory: FORM reproduced, magnitude → fine-grained).
- **(3b) Cytosim dynamic / Hand-kinetics parity** — extend the oracle to the motor/relaxation
  dynamics (needs Cytosim trajectory-output config; deferred, drag-model calibration).
- **(3c) WLC (Marko–Siggia) constitutive** — only if a flexible-filament consumer appears (actin is
  inextensible; currently no consumer — low priority).
- **(3d) Use FF for what it captures** — form / scaling / dynamics deliverables, γ vs motor density.

**Recommended.** (3a) or (3d) — they use the validated engine for results now; (3b) is a dedicated
parity build; (3c) is speculative without a consumer.

---

## DECISION 4 (optional, wet-lab) — commission the floor-overturning measurement

The single experiment that could overturn the γ-floor: **force-bearing (load-engaged) NMII
minifilament density per cross-section in MCF7** (super-res + an engagement readout). Imaging to date
counts *presence*, not *engagement*; real cells reach band at ~70% myosin ⇒ ~12–50× more engaged
motors/cross-section than the present-density (Nie ~0.6/µm²) gives. **Predicted to confirm the floor.**
Only pursue if you want an empirical overturn rather than accepting the result.
