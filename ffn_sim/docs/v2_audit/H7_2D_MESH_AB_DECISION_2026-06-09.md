# H.7 — "2D surface mesh to optimize the cortex" — (a) vs (b) adversarial study

**Date** 2026-06-09 · **Branch** `h7/full-cell-integration` · **Method** 5-agent workflow
(`wf_c27dc37d-e34`): advocate-A ∥ advocate-B → cross-adversarial audit ∥ → synthesis.
PI directive being interpreted: long-ago "use a 2D surface mesh to OPTIMIZE the cortex
network." Two readings studied:
- **(a)** wire the existing geometry manifold (`surface_manifold.py`) into the cortex as a
  SEARCH/MEMORY/COORDINATE accelerator ONLY (patch-local geodesic k-ring replacing per-tick
  global cKDTree; shared deformable normal frame; no mechanics, mesh-edges ≠ filaments).
- **(b)** reorganize the cortex NETWORK itself onto the 2D surface (connectivity + transmission
  in the surface metric).

## Verdict (both audits + synthesis agree)

**Option A SURVIVES** (moderate; no fatal flaw) — but it is ORTHOGONAL to the active-γ
transmission problem and its headline speedup is false on CPU:
- Fidelity-compliant: the sanctioned Layer-2 use (zero mechanics, no γ DOF; VG-6 candidate-set
  identity + VG-1 resolution-invariance PASS). NOT the rejected mesh-as-cortex.
- **CPU search is 7–30× SLOWER** than scipy's C cKDTree at the ×40 mesoscale (measured,
  `manifold_search_benchmark.py`); the index build itself is 25–140 ms > the 1.2–5 ms tree. The
  "~50–250×" speedup is unrealized; any win is conditional on GPU-residency (cupy).
- Real value = a **shared deformable coordinate frame** for a spreading/non-spherical cell
  (genuine, independent of speed) — not speed, not γ.
- NEW correctness hole the audit surfaced: lazy-refresh of a stale bead→patch map silently
  drops **21–56 % of true within-reach binding pairs** once beads drift ~1 patch (~1.2–2 µm) on
  a deforming cortex → that DOES change explicit connectivity. The current VG suite only proves
  identity on a STATIC sphere at refresh time → a refresh-cadence-vs-drift gate must be added.

**Option B does NOT survive** (major; `survives_audit=false`):
- Closes **0×** of the 30–80× γ gap — touches none of the soft scalar couplings nor the
  ~1.3-seg load-path.
- Its sole physics-touching claim (geodesic-vs-chord binding metric) is **illusory**:
  `kring_for_reach` is a broad-phase OVER-selector whose narrow phase still filters by Euclidean
  reach → the IDENTICAL chord candidate set (this is what VG-6 asserts). The geodesic
  differentiator only appears via an unjustified narrow-phase change with ~0 benefit at
  R=7.5 µm / 60 nm reach (chord ≈ arc to <0.01 %).
- A 4–7 day full-cortex rebuild, one tripwire (resolution-invariance) from the
  H7_CORTEX_AS_MESH-rejected coarse-graining, with deeper re-meshing/GPU-sync costs.
- "Correct substrate" framing = the refuted topology-blaming pattern in new dress.
- **If any spatial substrate is wanted, A strictly dominates B** on every axis. Decide
  **A-or-nothing, not A-vs-B**.

## The crux — NEITHER option fixes transmission

The active-γ transmission floor magnitude is set by THREE soft scalar stiffnesses + the
physiological load-path, verified at the literal level against live code:
- `configs/phase1_h3.yaml:122-123` `k_intra = k_attach = 1.0e-7 N/m` (crosslink, fiber↔fiber).
- `configs/phase1_h3.yaml:331,333` `k_head_spring = k_head_actin = 1.0e-6 N/m` (cross-bridge).
- `configs/phase1_h4.yaml:98` canonical motor `head_spring_k = 1.0e-3 N/m` → the **1000×
  myosin-head unit slip is real and in-config**.
- `crosslinkers.py:514` literally documents "k_intra = 0.1 pN/μm is so soft that any [stretch is
  force-free]".
- Load-path ~1.3 seg is a CONSEQUENCE of physiological dense crosslinking (z≈3.3, L/lc≈6) that a
  faithful rebuild REPRODUCES (H7_ARCHITECTURE_LOADPATH).

A coordinate/search/geodesic reorganization changes WHERE candidates are found and HOW FAST — it
changes none of these scalars and does not lengthen the load-path, so γ_soft (a sum over explicit
Cartesian bond tensions) is unchanged by construction. **The actual fix is SEPARATE and must
happen regardless of A/B**: the PI-gated stiffness unit-slip contract correction (k_head_actin
1e-6→~1e-3; crosslink k re-anchor) + the O2 myosin-overlap accumulation mechanism
(H7_MYOSIN_OVERLAP_MECHANISM_DESIGN). ~3-literal-class edit, entirely independent of 2D-vs-3D.

## Composability

A and B are NOT composable (mutually exclusive network-coordinate stance), but **A ⊂ B's
machinery** — A is the reversible first step B would also require. The stiffness/transmission fix
is fully ORTHOGONAL and composable with either or neither.

## Three independent decisions for PI

1. **TRANSMISSION (urgent, separate, do regardless):** the stiffness unit-slip contract
   correction + O2 myosin-overlap is the ONLY thing that touches the 30–80× γ gap. Surface as a
   brief-contract change with the in-config unit-slip evidence NOW; do not let a "cortex
   optimization" narrative displace it.
2. **DIRECTIVE INTERPRETATION:** A-or-nothing, not A-vs-B. B's differentiator is illusory and it
   buys 0× γ + ~0 search for rebuild + fidelity-slip cost.
3. **WHETHER TO DO A:** if yes, scope it honestly as **spatiality/coordinate frame, not speed,
   GPU-resident**, add the refresh-cadence-vs-drift correctness gate, and defer the k_conf
   confinement spring (magic-number, PI-gated). Tell PI plainly it is NOT a γ fix.

Artifacts: workflow `wf_c27dc37d-e34` (full advocate cases + audits in the run transcript).
