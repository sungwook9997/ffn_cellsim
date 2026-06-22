# Phase-C Warp DCM — Force-Law Correctness Audit (2026-06-23)

Adversarial physics-correctness pass over the force laws and their composition in the
Phase-C Warp DCM engine. Scope: `dcm_neighbor_warp.py`, `dcm_turgor_warp.py`,
`dcm_cadherin_host.py`, `dcm_ecm_clutch_host.py`, `dcm_lamellipodium_host.py`,
`dcm_substrate_warp.py`, `dcm_contact_warp.py`, `dcm_cohesion_warp.py`, and their
composition in `dcm_warp_decohesion.py::run_decohesion` (the device step loop) plus
the implicit operator `stiff_force_into` / `device_cg`.

Each candidate was pushed against the actual convention / a possible derivation before
being reported. SIGN traces and double-count claims survived refutation; the magic-number
items are flagged with their best in-repo anchor so PI can rule.

## Convention baseline (verified, NOT bugs)

These were checked and are CORRECT — recorded so they are not re-litigated:

- **Turgor sign.** `icosphere_mesh` forces OUTWARD winding (`dcm.py:57-64`,
  `n·(v_a−centroid)>0`). `cross=(v1−v0)×(v2−v0)=2A·n̂_out`; `dP=dP0+K_vol·(V0−V)/V` so a
  compressed cell (V<V0) has dP>0 and `ff=cross·dP/6` points outward → expands. ✓ Matches
  the numpy ground truth `DcmTurgorForce` (`dcm.py:148-155`).
- **Divergence-theorem volume** `V=(1/6)Σ v0·(v1×v2)` (`dcm_turgor_warp.py:47`) is the exact
  signed enclosed volume for outward-wound closed meshes. ✓ Internal turgor forces sum to ~0
  (each face's `ff` is split equally to its 3 nodes; a closed mesh's face-normal pressure is a
  self-equilibrated load up to the V-change work). ✓
- **Contact repulsion / adhesion direction** (`contact_grid_kernel`, `dcm_neighbor_warp.py:249-265`).
  `r_vec=node−closestpoint`, `sign=r_vec·n̂`. Repulsion (`sign<0`): node penetrated, `r_vec`
  points inward, `amp>0`, `fvec=r_vec·amp` inward, `node −= fvec` → node pushed outward. ✓
  Adhesion (`sign>0`): node outside, `r_vec` outward, `fvec` outward, `node −= fvec` → node
  pulled toward face. ✓ Newton-3 face reaction `atomic_add(ia, bary·fvec)` with node getting
  `−fvec` is momentum-clean (bary sums to 1). ✓
- **Cohesion tent**, **edge bond**, **cadherin bond** (`cadherin_bond_force_kernel:457-463`,
  attractive-only `L>r0`, equal-and-opposite), **ECM clutch** (one-sided to rigid dish, by
  design), **gravity** (`fz_node=−Δρ·g·v<0`, −z toward dish ✓), **bending** (`−k·Δ²r`),
  **surface tension** (`−γ·∂A/∂r`, area-reducing ✓) — all sign-correct.
- **Wetting** xy-gradient sign and the substrate z-well sign were checked against
  `dcm_substrate_warp` reference and are correct.

## Findings

| # | Issue | file:line | Rule violated | Severity | Confidence | Fix |
|---|-------|-----------|---------------|----------|------------|-----|
| 1 | **Triple-counted cell-cell adhesion** under `--cadherin --coupling`: the node-NODE cohesion tent (`cohesion_grid_kernel`, omega=`coh_adh`), the node-FACE bilinear adhesion (`contact_grid_kernel`, adh=`coh_adh`), AND the explicit cadherin trans-dimer bonds are all attractive and all active simultaneously. The cadherin docstring states it is "the SOLE cell-cell adhesion (run adhesion-OFF…)", and SIMUCELL3D_INTEGRATION §351 explicitly forbids it: "cadherin 모듈과 double-count 금지 — repulsion은 항상, adhesion은 선택". | `dcm_warp_decohesion.py:410` (`coh_adh = adh_strength if coupling`), fed to cohesion kernel `:578-582`/`:572-576` AND contact kernel `:606-610`/`:600-604`; cadherin `:659-664` | De-cohesion must EMERGE from cadherin rupture, not be one of three additive adhesions; no-double-count | High | High | In cadherin mode, route `coh_adh` to the contact (node-face) kernel ONLY, and pass omega=0 to the cohesion (node-node) kernel; OR (faithful) keep cadherin as the only attraction and use `--coupling` for node-face *non-adhesive* interface continuity only. At minimum, gate the node-node tent off whenever `cad is not None`. |
| 2 | **Double-counted excluded-volume repulsion**: every production step runs BOTH the node-node cohesion repulsion (`d<r_contact`, `fmag=rep·A·(r_contact−d)`) AND the node-face contact repulsion (`sign<0, min_d<c_rep`, `amp=rep·area`). A node approaching another cell's surface receives both penalties with the SAME `rep_strength`. The node-face penalty is the SimuCell3D-correct one (face-plane, shell-non-interpenetration); the node-node center-line repulsion is the superseded model (`dcm_face_contact.py:3` "Replaces the node-NODE tent contact"). | cohesion repulsion `dcm_neighbor_warp.py:191-192`; contact repulsion `:249-250` — both launched unconditionally in `step_once` (`:578`, `:606`) | No lumped/duplicated mechanism; mechanism should be single-sourced | High | Medium-High | Decide one repulsion source. Intended design (per `dcm_face_contact` docstring) is node-FACE only. Run `cohesion_grid_kernel` with `rep=0` (adhesion-tent only) and let `contact_grid_kernel` own excluded volume, OR document and justify the additive stacking with a magnitude budget. Currently `r_contact=0.30·edge` and `c_rep=0.30·edge` coincide, so the two repulsions overlap in range and stack. |
| 3 | **Deep-penetration tunnelling: zero restoring force past `c_rep`.** Contact repulsion fires only for `min_d < c_rep` (`=0.30·mean_edge`). A node that crosses deeper than `c_rep` into another cell gets `amp=0` — NO restoring force (the penalty vanishes exactly where it is most needed). Same in the cohesion node-node branch (`d<r_contact`). Faithful to the numpy reference (`dcm_face_contact.py:14`), so not a port regression, but a genuine physics hazard: the penalty is non-monotone in depth. | `dcm_neighbor_warp.py:249` (`sign<z and min_d<c_rep`); `:191` (`d<r_contact`) | Conservation / non-penetration (penalty must be monotone restoring) | Medium | High | The D8 `_vaxpy_active_capped` (cap=`c_rep`) + A3 substepping (`cfl_limit`) are the existing mitigations — but they are OPT-IN (`pen_cap` default True only on implicit path; `cfl_limit` default 0=off). Recommend: make the penalty monotone past `c_rep` (clamp `min_d` to `c_rep` inside the repel branch so amp keeps growing), OR always-on CFL substepping. The grid query radius (`con_q=c_adh+0.7·l_max`) is generous enough that a deep node is still FOUND — the cutoff is in the force law, not the neighbour search, so deepening it is safe. |
| 4 | **`k_vol = 7.73e5` driver default vs `ResolvedDCM.K_vol = 1e3`.** The osmotic bulk modulus used in production (`7.73e5`) is 773× the resolved-param default. It IS derivable (COMPARTMENT_ACTIVATION §370: Π_osm = R_gas·T·c_phys ≈ 7.74e5 Pa at 310 K) — so NOT a fabricated magic number. BUT it is the *osmometer restoring stiffness*, which that same doc (§394) flags is a DIFFERENT mode from the constant-target law the code integrates (`τ_Kvol = 1.7e3–1.9e4 s` ≫ physiological τ_RVD = 3–32 s). The two K_vol values (1e3 in the dataclass, 7.73e5 in the driver) disagree by 773× with no single sourced reconciliation in the code. | `dcm_warp_decohesion.py:203` (`k_vol=7.73e5`) vs `dcm.py:82` (`K_vol=1e3`); hardcoded `7.73e5` also in `dcm_warp_implicit.py:300,361,431` | No-magic-number / physiological-baseline (one sourced value) | Medium | Medium | Pin ONE K_vol with the derivation comment at the definition site and have both paths import it. If 7.73e5 (Π_osm) is the intended production value, update `ResolvedDCM.K_vol` and delete the bare `7.73e5` literals in `dcm_warp_implicit`. Surface the τ_Kvol vs τ_RVD mode mismatch (§394) to PI — it is a known open item, not closed here. |
| 5 | **`rep_strength = 2e8` default is ~5× stiffer than the MCF7 literature anchor (ξ≈4e7).** MCF7_PARAMETER_COLLECTION §85 derives ξ≈4×10⁷ Pa/m (ξ̄=0.48) from the MCF7 non-dimensional ratio; the `--rep-strength` help text itself says "soft lit value ~4e7 (MCF7 ξ̄)… the stiff 2e8 default over-packs/interpenetrates". The 2e8 default is acknowledged in-code as not the physiological value. (The 2026-06-23 aggregation re-test argued 2e8 gives cleaner pen_frac, but that is an outcome-tuning argument, not a derivation.) | `dcm_warp_decohesion.py:204` (`rep_strength=2.0e8`); help `:1091` | Physiological-baseline / no-tuning-to-outcome | Medium | Medium-High | Default to the derived ξ≈4e7 and treat 2e8 as an explicit override. If 2e8 is retained for numerical cleanliness, that is a stability-cap argument (like force_cap) and must be documented as such with PI sign-off — it cannot be the silent default while the help text calls it the over-packing value. |
| 6 | **`gamma_surf = 1e-4 N/m` surface-tension default is not anchored to MCF7.** No in-repo derivation hit. Memory record: direct MCF7 cortical tension is ~1e-2 N/m (Moazzeni 2021) and γ=2.7e-4 is explicitly "NOT measured MCF7" (SimuCell3D sub-floor). 1e-4 sits ~100× below the direct MCF7 datum with no citation. Default-off (`surface_tension=False`), which limits blast radius, but if turned on in production it injects an unphysical tension. | `dcm_warp_decohesion.py:218` (`gamma_surf=1.0e-4`); help `:1095` | No-magic-number; physiological-baseline when enabled | Low-Medium | Medium | Anchor γ to the MCF7 cortical-tension datum (~1e-2 N/m direct, or the resolved 3-channel estimator) with a citation, or surface to PI. Until then, document 1e-4 as PROVISIONAL and do not enable `--surface-tension` in a production/data-comparison run. |
| 7 | **`force_cap = 5e-8 N` shapes physics, not just stability.** Applied as a hard clamp in cohesion/cad kernels (`:199-202`) and the wetting cap and lamellipodium tether. A clamp that ever binds changes the force law (caps the repulsion that prevents interpenetration; the cohesion tent peak `ω·A·c_adh/2` can exceed it). Presented as a "BAOAB guard" but with no derivation that it sits above all physiological forces. | `dcm_warp_decohesion.py:207`; cohesion clamp `dcm_neighbor_warp.py:199-202` | No-magic-number (a cap that binds is a tuning constant) | Low | Medium | Verify 5e-8 N is strictly above every physiological per-node force at the operating point (so it never binds in production) and assert/log when it binds. If it binds, it is silently re-shaping the force — surface to PI. |
| 8 | **Cadherin bundle force/k_off uses molecular load, but the kernel applies bundle force — and the host break-test recomputes load on the host without the bundle factor consistently.** Force applied = `k_trans·bundle_n·(L−r0)/L` (`:662`), while the break rate evaluates `k_off` at per-molecule load `F=k_trans·(L−r0)` (`dcm_cadherin_host.py:112`, no `·bundle_n`). That IS the documented intent (force ×N, k_off at molecular F/N). Verified consistent — NOT a bug. Recorded as checked. | `dcm_warp_decohesion.py:662`; `dcm_cadherin_host.py:111-113` | — (refuted) | — | — | No action; the bundle bridge is internally consistent. |

## Summary of real (survived-refutation) findings, ranked

1. **(High) Triple adhesion double-count** under `--cadherin --coupling` (#1) — violates the
   explicit no-double-count contract and the emergent-de-cohesion rule. The node-node tent must
   be off in cadherin mode.
2. **(High) Repulsion double-count** node-node + node-face every step (#2) — same `rep_strength`
   applied twice in overlapping ranges; the node-node center-line model is the superseded one.
3. **(Medium) Deep-penetration zero-force tunnelling** past `c_rep` (#3) — non-monotone penalty;
   mitigations exist but are opt-in/off-by-default on the explicit path.
4. **(Medium) K_vol 7.73e5 vs 1e3 inconsistency** (#4) — derivable but two un-reconciled values
   + a flagged mode mismatch.
5. **(Medium) rep_strength 2e8 vs lit 4e7** (#5) — in-code-acknowledged non-physiological default.
6. **(Low-Med) gamma_surf 1e-4 unanchored** (#6) and **(Low) force_cap 5e-8 may bind** (#7).

Signs (turgor outward, repulsion outward, adhesion inward, gravity −z, surface-tension
area-reducing, Newton-3 reactions) all traced CORRECT. Momentum: internal pairwise forces are
equal-and-opposite; the only intentionally one-sided forces (ECM clutch, lamellipodium tether,
substrate well, gravity) are external by design (rigid dish / pinned actin / body force).
