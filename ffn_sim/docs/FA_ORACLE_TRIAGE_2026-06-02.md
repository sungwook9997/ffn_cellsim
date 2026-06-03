# FA literature triage — 3 papers → oracle / observable / constant candidates (2026-06-02)

> **STATUS: TRIAGE scaffold — NOT a contract, NOT ratified bands.** Classifies three
> PI-supplied focal-adhesion papers against (a) the current production wall and (b) the
> FA-integration α-restart + PI-exp validation roadmap. Every item below enters the
> platform **only** as an acceptance oracle, a literature constant, or an observable
> target — **never a runtime mechanism** (CLAUDE.md inversion rule). Author: Lead session
> on direct PI request. Companion to
> [`briefs/H4_FA_INTEGRATION_DESIGN.md`](briefs/H4_FA_INTEGRATION_DESIGN.md),
> [`PI_EXP_VALIDATION_MAP.md`](PI_EXP_VALIDATION_MAP.md),
> [`BETA1_DISTRIBUTION_OBSERVABLE.md`](BETA1_DISTRIBUTION_OBSERVABLE.md).

## Architectural principle (why these are oracles, not mechanisms)

All three are **lumped / node-level models or pure measurement**:
- Kim–Kamm–Asada 2012 — node-level ODE (integrin *bundle* per node, stress fiber = a
  variable-stiffness spring + a myosin-sliding ODE, friction = a lumped rupture term).
- Rumi De 2018 — **mean-field** master equation (all bonds in a cluster share one force).
- Zhuo 2018 (PROM) — experimental imaging; the photonic-crystal physics is irrelevant,
  only the FA quantification transfers.

Per CLAUDE.md (PI 2026-05-19, hard rule) the platform **inverts** such paper-models: they
become validation **oracles** / literature constants, never the runtime. Two of the items
below merely **corroborate anchors we already use** (Pereverzev `k_off`; Kanchanawong
λ = 30 nm).

## Relation to the current production wall (force isotropy) — honest scope

The current γ-floor wall is **cortex-internal actomyosin coherence** (commit `c3fafed`:
random minifilament orientation cancels in the great-circle tension sum; a coherent field
gives γ = 0.30 mN/m in-band). Its proper literature is **active-gel / active-nematic
cortical tension** (Salbreux–Jülicher, Mayer/Grill, Fischer-Friedrich, Hannezo) — **none
of these three FA papers address it directly.** Their relevance to the wall is *indirect /
conceptual* (recorded in §4); the direct value is to the FA roadmap (§1–§3).

---

## 1. Kim, Kamm & Asada 2012, *Integr. Biol.* 4:1386–1397 — integrated FA↔SF↔nucleus traction ORACLE

DOI 10.1039/c2ib20159c. Integrated migration model: integrin clutch + stress fibers +
myosin + nucleus on a 3D curved ECM surface. **Closest published analog to the
ffn_cellsim target architecture.** Kamm/Asada (MIT) group; PDF retrieved via KAIST
(likely PI-adjacent). **Role: acceptance oracle + literature constants — NOT runtime.**

### Maps to — `H4_FA_INTEGRATION_DESIGN.md` build sequence

The design's **S2 `fa_actin_clutch`** stage ("NEW topology: cortex/stress-fiber actin →
engaged integrin … *no current module provides it*") **is exactly this paper's load
path** (SF from integrin node → nucleus node). **S5** (close the clutch loop → γ) and the
**S7 membrane / H.9 nucleus** coupling are its force-balance.

### Extractable literature constants (all carry primary provenance in the paper)

| Quantity | Value | Paper basis | Our attach point | Flag |
|---|---|---|---|---|
| FA force law | `F_FA = n_b·k_LR·(L_b−λ)` (eq 6) | — | S2 `fa_actin_clutch` force form (oracle) | oracle |
| Single bond stiffness k_LR | **≈ 1 pN/nm** | Dembo 1994 (ref 31) | clutch spring `κ_c` cross-check (cf. Bangasser 0.8 pN/nm) | oracle |
| Integrin equilibrium length λ | **30 nm** | Kanchanawong 2010 *Nature* (ref 32) | **already our anchor** — corroboration | ✓ existing |
| SF stiffness | `k_SF = E_SF·A_SF/L`; **E_SF = 230 kPa**, A_SF r = 250 nm | Deguchi 2006 (ref 33); ref 34 | stress-fiber bond modulus (S2/S5) | oracle |
| Non-muscle myosin-II sliding v_m | **10 nm/s** | refs 35–37 | cross-check vs our grip-walk / Hill v_u (Bangasser 120 nm/s) | oracle |
| Critical bond-formation height h_c | **300 nm** | current work | FA capture geometry (S1/S2 talin/clutch capture distance) | oracle |
| Friction (bond-rupture dissipation) C_c | 0.001 N·s/m | refs 28–30 | drag oracle (lumped — NOT our per-bond runtime) | proxy |
| SF assembly / disassembly time | ~180 s / ~1 s | refs 38–44 | turnover timescale (H4 co-req "actin turnover") | oracle |

### Validation curves (oracle-only comparison)

- **Cancer-cell migration speed vs conduit width** — reproduces Irimia–Toner data,
  max speed at ~20–30 µm width / ~10 µm cylinder, **R² = 0.77** (Fig 4). → future
  confined-migration validation; pairs with [`MIKADO_3D_DESIGN.md`](MIKADO_3D_DESIGN.md)
  confined-channel preset.
- **Nucleus aspect ratio under confinement** — **1.91 (<25 µm lumen) vs 1.37 (>25 µm)**
  (Fig 1E). → H.9 nucleus-deformation oracle (currently deferred integration).

### Relevance to the current wall (§4 link)

Its SF traction is **large when SFs coherently align with the load axis, near-zero
(cancellation) when misaligned** — the closest *published* analog of the γ
isotropy-cancellation, and it shows alignment is **geometry/anchoring-emergent**, not
imposed. Conceptual cross-reference for the `FFN_MYOSIN_ALIGN` diagnostic (§4).

---

## 2. Rumi De 2018, *Commun. Biol.* 1:81 — cluster-level validation of our EXISTING Pereverzev clutch + orientation observable

DOI 10.1038/s42003-018-0084-9. Catch-bond FA cluster master equation; FA orientation
under static vs cyclic stretch. **Role: acceptance oracle (cluster level) + a future
observable — NOT runtime.**

### Key fact — its `k_off` IS our runtime mechanism

`k_off = k_slip·exp(f_b/f_0) + k_catch·exp(−f_b/f_0)` (eq 2) is the **Pereverzev
two-pathway catch-slip** model (its ref 39 = Pereverzev 2005) that the platform
**already runs**: [`bridge/integrin_bonds.py`](../bridge/integrin_bonds.py)
`::pereverzev_k_off`. This paper is therefore an **independent application of our own
off-rate law at the N-bond cluster level.**

| Item | Value / form | Our attach point | Flag |
|---|---|---|---|
| Catch-slip `k_off` | eq 2 = Pereverzev | **already runtime** (per-bond, finer than this mean-field) | ✓ existing |
| Cluster master equation | eq 1, Gillespie SSA (ref 42); N=200, Ks=0.10, Kc=120, Γ=2 | **cluster-stability oracle** for KU-2.x: growth→optimum→threshold-disassembly vs strain (Fig 2b) | oracle |
| Rate-dependent `k_on` | `γ(N−n)·exp(−v_ω²/v_0²)` (eq 3); v_ω = l_R·ε₀·ω | **NEW vs platform** — future **cyclic-loading** option only; not needed for static FA/cortex | future |
| FA orientation | parallel (static) / perpendicular (fast cyclic); threshold stretch; cell-type via T_R | future FA-orientation observable; the §4 alignment motif | observable |

### Caveat

Mean-field lumped (one shared force per cluster) → **oracle only**; our per-bond runtime
is already finer-grained. The rate-dependent `k_on` is a genuinely new mechanism but
applies to **dynamic (cyclic) substrate loading**, which the current static cortex/FA
program does not exercise — flag as a deferred extension, not a now-need.

---

## 3. Zhuo et al. 2018, *Light Sci. Appl.* 7:9 (PROM) — FA size band + periphery "ring" OBSERVABLE TARGET

DOI 10.1038/s41377-018-0001-5. Label-free FA imaging. **Only the FA quantification
transfers** (the photonic-crystal method is irrelevant). Cells = **mHAT9a dental
epithelial stem cells on fibronectin — NOT cancer/MCF7** ⇒ qualitative / order-of-mag
only. **Role: observable target + a size-distribution reference.**

| Quantity | Value | Our attach point | Flag |
|---|---|---|---|
| Focal complex (nascent) | **< 0.2 µm²** | FA size distribution / `A_baseline` calibration | qual |
| Mature FA | **1–10 µm²** | same | qual |
| Typical FA cross-section | **0.2–1.0 µm²** | KU-2.17 FA-growth size check ([`bridge/fa_growth.py`](../bridge/fa_growth.py); independent of the existing **Stricker 2013** size ref) | qual |
| Spatial: FA "ring" at cell periphery (vinculin-colocalized; **Edge > Inner**, Fig 4c/5c) | edge-localized | **direct target for `f_edge`** in [`bridge/clutch_spatial.py`](../bridge/clutch_spatial.py) (`beta1_distribution_metrics`; uniform-disk ref = 1−ρ²) ↔ PI-exp Bare(diffuse)/Pre(peripheral)/Lam4(uniform) | observable |

### Note

KU-2.17 = the Walcott–Sun 2010 FA-growth Hill oracle (`test_ku217_fa_growth.py`); the
runtime is the emergent `FAGrowthMonitor` (N_engaged / F_per_FA), **no Hill ODE at
runtime.** Zhuo supplies the **area band** the emergent `A_emergent` should land in, and
the **periphery pattern** the `clutch_spatial` β1-observable should reproduce — both
*pattern/ordering* checks, never absolute fits (esp. given the non-cancer cell type).

---

## 4. The current wall — what these papers DO and DON'T do

- **DON'T**: solve the cortex-internal force-isotropy wall (commit `c3fafed`,
  [`KU35_FLOOR_ROOT_CAUSE`](v2_audit/KU35_FLOOR_ROOT_CAUSE_2026-05-31.md)). That needs
  active-gel / active-nematic cortical-tension literature (separate sweep, not run).
- **DO (indirect, conceptual)**:
  1. **Both** Kim–Kamm–Asada and Rumi De demonstrate the exact motif the wall hinges on —
     a population of force-bearing elements yields large **net** force only when
     **coherently oriented**; misalignment cancels — independently confirming the
     cancellation is a *real, general phenomenon*, not a platform artifact.
  2. Both show orientation is **selected by an anisotropic field / geometry / anchoring**,
     not imposed — reinforcing the **α-restart bet** that closing the FA/substrate/membrane
     boundary loop ([`H4_FA_INTEGRATION_DESIGN.md`](briefs/H4_FA_INTEGRATION_DESIGN.md))
     supplies the boundary the free shell lacks.
  3. Neither resolves **spontaneous symmetry-breaking in a round cortex** (no preferred
     axis without a polarity cue) — that remains the active-gel open question.

---

## 5. Open items (candidate → PI ratification)

- [ ] Ratify whether Kim–Kamm–Asada enters as a named **S2/S5 traction oracle** in the H4
      design (alongside Chan-Odde) + its constants as cross-checks (k_LR, E_SF, v_m, h_c).
- [ ] Ratify Rumi De as the **cluster-stability oracle** for the KU-2.x catch-bond cluster
      (growth→optimum→disassembly vs strain). Rate-dependent `k_on` → defer (cyclic only).
- [ ] Ratify Zhuo FA size band (0.2–1.0 µm²) + periphery "ring" as **observable targets**
      for `fa_growth` (KU-2.17) and `clutch_spatial` `f_edge` — qualitative, non-cancer.
- [ ] Reverse cross-links (one-liners the Lead can add at integration): pointer to this
      triage from `H4_FA_INTEGRATION_DESIGN.md` §5 and `PI_EXP_VALIDATION_MAP.md` Layer-1.
- [ ] (separate decision) active-gel cortical-tension sweep for the wall itself — offered,
      not run.
