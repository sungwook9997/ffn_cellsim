# DCM spheroid adhesion diagnosis — why cells stay "bag of marbles" (2026-06-30)

PI (furious, repeatedly-flagged): the simulated spheroid cells do NOT adhere/deform — they stay
rounded with low contact (~0.13–0.20) instead of flattening into a compact tissue; raising cortical
tension γ makes it WORSE; γ=0 looks most natural. PI suspected the face-face adhesion does not scale
with contact AREA. **Diagnosis (12-agent ultracode workflow: code-read + literature + adversarial
verify): PI is right.**

## Verdict

Three cell-cell adhesion paths exist; the **default production path is NOT area-proportional**:

| path | code | area-scaling? | used in production? |
|---|---|---|---|
| node-NODE cohesion tent | `dcm_cohesion_warp.py` / `cohesion_grid_kernel` | **NO** — `A = area_per_node = 4πR0²/npc` frozen at init (`dcm_warp_decohesion.py:341`), never updated as the cell deforms | YES (default) |
| explicit cadherin bonds | `dcm_neighbor_warp.py:438` | **NO** — plain harmonic tether `F=k(L−r0)`, ~6 bonds/junction (point contact) | only with `--cadherin` |
| **node-FACE contact** | `dcm_contact_warp.py:120-133` / `contact_grid_kernel` | **YES** — `amp = ω·(c_adh/min_d−1)·A`, `A = ½‖(b−a)×(c−a)‖` (live face area, every step). = SimuCell3D Model-1 | runs, but **dominated** |

The correct **area-coupled** node-face adhesion (SimuCell3D Model-1, the `W_adh × A_contact`
traction) **already exists**. The reason cells still round:

1. **Balance is wrong (Regime I, tension-dominated):** production runs `rep=2e8` vs `adh=5e7`
   (**repulsion 4× adhesion**) → repulsion wins at contact, cells can't flatten. The node-node
   constant-A cohesion runs alongside and dilutes the area-coupled term.
2. **Cortical tension is area-minimizing and adhesion can't out-scale it:** `surface_tension_kernel`
   applies `F = −γ·∂A/∂r` on ALL faces. As the contact tries to grow, γ's opposition grows ∝ area
   but the (default) adhesion does not → contact pins tiny. Committed γ-sweeps confirm: contact_area
   falls **monotonically to 0** as γ rises (0.078→0.060→0.037→0.002→0); V/V0 crushed 1.09→0.20. This
   is why γ=0 looks "most natural" — it stops killing the contact.

## Real spheroid (the target we're missing)

A real MCF7/T47D spheroid is a **jammed wet active foam** (Steinberg DAH; Manning PNAS 2010;
SimuCell3D Runser 2024; Maître & Heisenberg 2015):
- **Large flat shared faces, ~2-3% core porosity** (T47D EM PMC10212087: interstitial space inner
  12.5%→**2.4%** day5→day20 = essentially space-filling/jammed; confluent → Kelvin ~14 neighbours).
- **Central cells deform MORE than rim** (pressure-driven; T47D inner ~224 µm² vs outer ~164 µm²) —
  a concrete validation target our turgor-pressured DCM should reproduce.
- Adhesion energy ∝ **contact area** (`W_adh × A_contact`) is the load-bearing term that drives
  faceting; cadherin locally **down-regulates actomyosin at the contact** (β_contact < β_free,
  Maître Science 2012) so γ is reduced WHERE cells touch — cadherin binding energy alone (~2-7% of
  cortical tension) is far too small to open a contact.

Our "Ψ≈0.99 bag of marbles / contact 0.13-0.20" is the **opposite** — stuck in Regime I.

## The fix (PI-directed: rebuild the cell-cell interaction to SimuCell3D-faithful)

1. **Area-coupled node-face adhesion as the primary cell-cell path** (the SimuCell3D Model-1 traction
   `F = min{A_a,A_b}·σ(d)`, σ = ξ·d overlap / ω·d softening / ω·(c−d) hardening, barycentric to
   nodes). It exists; make it dominant, demote the node-node constant-A tent + sparse point-bonds.
2. **Regime II balance (a RATIO, derived — NOT tuned to outcome):** `γ̃=γ/(Kℓ)` low, `ω̃=ωℓ/K` high;
   crossover where adhesion flattens = **ω > 2γ at the contact** (Manning γ/β>2). Derive ω from the
   lit adhesion energy density `w_cs = 2.85e-3 J/m²` × area; do NOT fit ω to make a gate pass.
   Repulsion ξ must not dominate ω (ξ≳ω is enough for non-penetration; ξ=4ω is Regime I).
3. **Contact-localized γ reduction** (the γ-floor's other half): make γ a field cadherin engagement
   locally lowers (`polarized_surface_tension_kernel`, Douezan S=w_cs−2γ) — Maître's contact-expansion
   mechanism — so γ stays at its physiological value (physiological-baseline rule) instead of dialed
   to 0. De-cohesion stays emergent (catch-bond rupture, not a binary latch).

## Open / needs the test run (in progress)

- **Metric check:** the "0.13-0.20" is `contact_area_frac` from `gamma_sweep` (a real contact metric),
  NOT `pen_frac` — symptom is real. (pen_frac=8000 in the division run is a separate degenerate-face
  metric artifact; that aggregate was compact/non-ejected.)
- **Decisive test (gbook, N=48):** A current(rep2e8/adh5e7/γ1e-3=marbles) vs B Regime-II(rep5e7=adh5e7,
  γ off) vs C adh-dominant+coupling+cadherin. If B/C contact ≫ A(0.13) → the fix is a balance/wiring
  change (not a from-scratch rebuild); if not → deeper. Results pending.

## Resolution attempts (N=48 tests, gbook)

Tested in order; contact_area_frac (a real metric, not pen_frac) + per-cell asphericity + V/V0:

| lever | contact | cellAsph | V/V0 | verdict |
|---|---|---|---|---|
| baseline (rep2e8 adh5e7 γ1e-3 uniform) | 0.105 | 0.018 | 0.99 | round marbles |
| balance (rep=adh, γ off) | 0.165 | 0.015 | 1.09 | still round |
| adhesion 10× (adh5e8, ω̃≈2.4) | 0.190 | 0.028 | 1.09 | **STILL round — adhesion strength is NOT the lever** |
| differential tension (γ_contact=γ−w_adh/2) | 0.11–0.18 | 0.01–0.05 | 0.4–1.1 | high γ → cells SHRINK apart (S=w_cs−2γ<0) |
| **Maître + stiff turgor** (γ1e-2 + contact cortex dissolved frac0 + K5e4) | 0.168 | **0.105** | **1.03** | **cells become angular POLYHEDRA (asph 6× baseline), volume held** |

**Findings:**
1. **Adhesion strength is NOT the lever** — even 10× lit adhesion (ω̃≈2.4) leaves cells round (contact
   saturates ~0.19). Confirmed visually.
2. **Plain differential tension is insufficient** — at MCF7 γ (~1e-2), the Douezan S = w_cs − 2γ < 0
   (non-wetting), and soft turgor lets γ SHRINK the cell instead of flattening it.
3. **The working mechanism = Maître contact-cortex-dissolution + volume conservation**: high γ_free
   (penalises free area) + γ_contact → 0 (cadherin dissolves the cortex at contacts) + STIFF turgor K
   (holds V/V0 ≈ 1 so the cell cannot shrink → must deform). → cells become **angular polyhedra**
   (cellAsph 0.018 → 0.105) at MCF7 physiological γ. Implemented: `--diff-tension --contact-tension-frac`.
   K must scale with γ to hold volume (K≈5e4 balances the 2γ/R≈2.7kPa Laplace at γ=1e-2 — a derived
   dependency, NOT a free knob; whether MCF7 cytoplasm K is really ~5e4 is a PI/lit anchor question).
4. **Still PARTIAL at N=48** — contact only 0.168 (cells faceted in SHAPE but not yet space-filling;
   the free aggregate has gaps + a free boundary). Full foam (contact → 0.5-0.9) needs **N≥400
   confinement** (interior cells jammed) — validation run launched (`n400_facet/N400_facet_mcf7`).

→ The mechanism (the missing physics) is now in place and visibly faceting cells; the remaining
question is whether N≥400 confinement closes it into a space-filling foam, and the K-anchor.

## VERDICT (N=400 confinement run, the decisive test)

`n400_facet/N400_facet_mcf7` (N=400, MCF7 γ=1e-2, contact-cortex dissolved, K=5e4, 40k steps,
V/V0 held 1.035): **contact = 0.237** (up from 0.168 at N=48 — confinement helped +40%), cellAsph
0.075. Visualized (`_view_N400.png`): a **compact aggregate of somewhat-angular cells**, but each
cell keeps its rounded-polyhedral identity with visible gaps — **NOT a space-filling wet foam**
(which needs contact ~0.5–0.9 with large flat shared faces).

**Conclusion — the patching approach has hit its ceiling.** Confinement + the Maître + stiff-turgor
+ differential-tension stack gets partial deformation (angular cells, contact 0.24) but does NOT
faithfully reproduce SimuCell3D-quality faceting, and getting even this required co-setting several
non-independently-anchored knobs (contact-cortex→0, K 20× to 5e4, MCF7 γ) — edging toward
outcome-tuning, which is not faithful physics.

**RECOMMENDATION (= PI's "전부 다 갈아엎어"):** stop patching the force-superposition model. Do a
**faithful SimuCell3D port** — define the actual ENERGY functional (per-cell volume constraint
`p=−K·ln(V/V0)` + per-face interfacial tension γ_ij that is LOW on cell-cell faces, HIGH on free
faces) and minimize it with their integrator + remesh, so faceting EMERGES rather than being coaxed
by tuned knobs. Salvageable building blocks: this diagnosis, `dcm_interfacial_tension_warp.py` (the
per-face γ_ij idea), the contact-detection + turgor + remesh kernels. PI decision pending on
committing to the energy-based rewrite.

Source: 12-agent workflow `wf_2088b933-005` (code diag + biology research + adversarial verify +
synthesis). Citations: Steinberg DAH; Manning et al. PNAS 2010; Maître & Heisenberg 2015 (Science
2012); SimuCell3D Runser/Vetter/Iber 2024; CellSim3D Madhikar 2018; T47D EM PMC10212087;
arXiv:2408.07551 (pressure-gradient deformation).
