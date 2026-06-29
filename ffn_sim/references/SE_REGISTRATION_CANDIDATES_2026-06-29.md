# SourceEvidence registration candidates — 2026-06-29

Two engine-reference papers added to `references/` during the FF/DCM two-layer restructure.
Per the corpus workflow these are **candidates only** — they enter BM25 on the next ingest,
but a Notion SourceEvidence row is a manual PI-gated step (no auto-register). Verify DOIs
before citing in any deliverable (citation-integrity rule).

---

## 1. Cytosim — `ff/` (Filament-FEM) physics reference

- **Citation:** Nédélec F. & Foethke D. (2007). *Collective Langevin dynamics of flexible
  cytoskeletal fibers.* New J. Phys. **9**:427.
- **DOI:** `10.1088/1367-2630/9/11/427`  (open access; confident)
- **Code:** https://gitlab.com/f-nedelec/cytosim  (GPLv3)
- **Type:** Method / model (overdamped implicit Langevin dynamics of fibers represented by points).
- **Why:** the going-forward physics basis of the `ff/` Filament-FEM engine — fibers as bending
  beams, motors/crosslinkers as Hands, the implicit large-timestep solver. Cytosim itself is the
  intended **independent parity oracle** for `ff/` (mirrors `dcm/` ← SimuCell3D). See `ff/ENGINE.md`.
- **PDF status:** ✅ in `references/` — **`Nedelec_2007_New_J._Phys._9_427.pdf`** (the published
  IOP NJP version, PI-added 2026-06-29 — authoritative for citation). The arXiv preprint
  `0903.5178` is also present as a mirror. Grounding (Stage 6b) unblocked.
- **Suggested KB linkage:** new ModelContract for the FF fiber-mechanics engine; KnowledgeClaims for
  the discrete bending operator + the implicit integration scheme.

## 2. CellSim3D — `dcm/` (Deformable Cell Model) reference

- **Citation:** Madhikar P., Åström J., Westerholm J. & Karttunen M. (2018). *CellSim3D: GPU
  accelerated software for simulations of cellular growth and division in three dimensions.*
  Comput. Phys. Commun. **232**:206–213.
- **DOI:** `10.1016/j.cpc.2018.05.024`  (⚠️ VERIFY — derived from PII S0010465518302091, not confirmed)
- **PII:** `S0010465518302091`  (definite — from the PDF filename)
- **Code:** https://github.com/SoftSimu/CellSim3D  (C++/CUDA/Python)
- **Type:** Method / model (GPU deformable-cell colony growth + division).
- **Why:** a deformable-cell-model reference for the `dcm/` layer, registered alongside SimuCell3D
  (per PI 2026-06-29). GPU (CUDA) elastic-shell cells + internal pressure + intercellular forces +
  **division** (cf. our `dcm_cleave`), fixed-topology-on-GPU (cf. our node-pool). Useful as a
  cross-check for dcm physics and a CUDA implementation reference for the Warp port. NOT a filament
  model (no cytoskeleton) — does not serve `ff/`. See `dcm/ENGINE.md` §Physics references.
- **PDF status:** ✅ in `references/` (`1-s2.0-S0010465518302091-main.pdf`, added 2026-06-29).
- **Suggested KB linkage:** attach to the existing DCM ModelContract as a secondary
  deformable-cell-model reference (SimuCell3D primary).
