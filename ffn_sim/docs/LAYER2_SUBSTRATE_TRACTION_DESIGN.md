# Layer-2 D-axis — substrate-reacted active traction (CBM extension) — DESIGN

> PI-authorized 2026-06-04 (build the substrate-traction CBM extension at the protrusion 9.4 nN
> anchor; sub-stepped / substrate-react extension OK = freeze-adjacent D3). Design produced by a
> 3-proposal + adversarial-judge workflow; knife-edge resolved by a KB anchor sweep.

## Why (the finding it fixes)

The lamellipodium→CBM motility bridge (`spheroid/motility_bridge.py`, REPORT §C) localised the
~5–9× A/A₀ magnitude gap to the **active** side and showed the protrusion-anchored per-cell
traction `f_active = v0·γ_cell = 9.4 nN` (Bieling barbed-end v0=1.88 µm/min × γ_cell=0.30 N·s/m)
**exceeds the cell-cell cohesion** F_detach = 6.5 nN (Iturri). In the current CBM the active force
is a **cell-cell** force, so at 9.4 nN a rim cell **detaches** rather than spreading — the
center-based 1-particle structural limit. Cohesion is anchored two independent ways (Iturri + D2
γ), so it is NOT the cause; the fix must be on the active side.

## Recommended mechanism — "In-plane substrate FA-clutch load path"

(Synthesis of 3 designs; all three as-written were rejected — P1's viscous damper vanishes at
v→0 so it only *slows* detachment; P2 prescribes position, bypassing cohesion feedback; P3
anchored the clutch along z so it absorbed zero of the in-plane traction. The synthesis grafts
P3's force-partition + Bell-Evans occupancy + freeze-adjacent preprocessing, P1's γ_cell-consistent
anchoring, and fixes the geometry to **in-plane**.)

A rim cell's 9.4 nN protrusion traction has **two in-plane load paths in parallel**:
1. **Cell-cell Morse/catch cohesion** — physically capped at F_detach = 6.5 nN (the well cannot
   pull harder; any *sustained* outward force > 6.5 nN detaches — the structural truth).
2. **An explicit in-plane substrate FA clutch** — an integrin→ligand anchor to a fixed substrate
   point, holding a steady in-plane force as the cell creeps and slips per Bell-Evans (the
   molecular clutch, `bridge/integrin_bonds.py` / `ligand_species.py`). The substrate is the load
   ground: its reaction is on the dish, **not** on neighbours.

Force balance on a creeping rim cell: `f_active = F_cohesion(neighbour) + F_clutch(substrate) +
γ_cell·v_creep`. The clutch carries the **excess above cohesion**, so the cell-cell tension stays
≤ 6.5 nN and the cell **crawls instead of tearing off**; the aggregate spreads collectively. The
crawl speed is the clutch-limited `v0 = f_active/γ_cell` (the motility-bridge identity), and the
clutch turnover (engage→hold→slip→re-engage) is what permits net forward creep.

**Why it bounds tension (where the rejected proposals failed):** it is a HELD elastic/clutch force
(persists at steady state, unlike a viscous damper); a real force in `net_force` (cohesion feeds
back); and **in-plane** (competes for the 9.4 nN).

## The knife-edge — RESOLVED in favour (the cheap decisive check)

Capacity question: can the substrate clutch carry the excess `f_active − F_detach = 2.9 nN`?
- A **single** FA ensemble (N_c=50, F_s(col-I)=18.6 pN) holds only ~0.7–0.9 nN — *below* 2.9 nN.
- But a spreading cell has **~25–60 focal adhesions** (KB-2.2/2.18), each 1–10 nN, so the
  **per-cell substrate capacity is 10–100 nN** (KB-2.12); MCF7-specific collagen-I de-adhesion
  matures to **~20 nN** (Taubenberger 2007, KB-PIV-3, High confidence). **20 nN ≫ 2.9 nN — and ≫
  the full 9.4 nN traction.** So the per-cell substrate clutch EASILY holds the cell ⇒ **no
  detachment; the mechanism is viable.** The anchor is the *per-cell FA capacity* (measured), so
  there is **no `m_rup` magic-number** to pick (the earlier single-FA framing was the wrong scale).

**Ligand axis falls out for free:** col-I (lower k_off0, higher F_s) → larger per-cell capacity →
more protected; laminin (registry 0.61× col-I) → smaller capacity → earlier slip — the measured
single-cell ordering, never engineered.

## The remaining real question (kill_if-2) — only the sim answers it

With the substrate holding the cell (no detachment), does the cell **crawl fast enough** to close
the gap, or does the substrate **pin it (stasis)**? If clutch turnover admits net creep at
v0≈1.88 µm/min, the aggregate spreads at the protrusion scale and A/A₀ may reach ~7–10 (gap
closes). If the substrate over-pins, A/A₀ stays near the in-band 1.6 nN whole-cell value (gap
stays open). **Falsified** if ARM B (clutch on) holds the cells (detached_frac < 0.05) but A/A₀
stays within seed-noise of the clutch-off / whole-cell runs — i.e. substrate traded detachment for
stasis, not spread → the gap is NOT on the substrate-reaction axis (push back to collective
force-amplification / the Lam4 single↔collective split / the density axis).

## Implementation plan (additive, frozen `baoab.py` untouched)

0. Anchor block: per-cell substrate capacity F_cap_cell = n_FA · F_per_FA, anchored to MCF7
   Taubenberger ~20 nN (and the 10–100 nN generic band); k_clutch = clutch_stiffness (KU-2.18,
   reused, not new); k_off0/x_beta from the `ligand_species` registry. NO new fitted constant.
1. `spheroid/substrate_fa_clutch.py`: `resolve_substrate_fa_clutch(resolved, ligand)` (pure
   arithmetic, unit-tested) + `SubstrateFAClutchForce(md.force.Custom)` — per active basal cell an
   in-plane anchor set on substrate contact; held force F=−k_clutch·dx capped by the **emergent**
   ensemble capacity (engaged-fraction evolves per Bell-Evans slip + re-anchor; z-component zero).
2. Wire into `run_growth_pooled` via optional `substrate_fa_clutch=None` (default → bit-identical
   to today); append the force after the z-wall + edge SettableForce; update anchors/held-forces in
   the existing per-epoch recompute block. `baoab_device.py` untouched.
3. Tests: F_cap = anchored capacity (no hidden constant); n_clutch→0 ⇒ bit-identical baseline;
   load saturation (9.4 nN pull → held ≤ capacity then slips); ligand ordering col-I > laminin;
   in-plane only.
4. **Micro 2-cell harness (no GPU):** rim+neighbour on substrate, 9.4 nN outward; clutch ON holds
   T_cohesion < 6.5 nN and the rim creeps; OFF detaches. The cheap decisive mechanism check.
5. **Decisive GPU experiment (gbook A5000, native N):** 3 arms at f_active=9.4 nN (protrusion
   anchor, edge-localised): ARM A clutch OFF (detachment baseline), ARM B clutch ON col-I (the
   test), ARM C clutch ON laminin (ligand ordering). Record A/A0(t), detached_fraction, ejected,
   max|net_force|, held-force distribution; overlay PI A/A0 (overlay-only). Cross-check vs the
   1.6 nN whole-cell run. **Gap closes iff ARM B A/A0 ≫ ARM A and ≫ the 1.6 nN run with
   detached_fraction < 0.05.**

## Open PI flags

- **Density axis (inherited, not new):** the per-cell n_FA scaling with substrate ligand DENSITY
  is the same pV4D4 col-I adsorption-density gap already flagged (`ligand_traction.DENSITY_FACTORS`)
  pending the KAIST collaborator — caps how quantitative Bare-vs-Pre can be. The ligand-IDENTITY
  axis (col-I vs laminin) is fully anchored.
- **Freeze-adjacency confirm:** the new force enters `net_force` and is updated per-epoch in
  `run_growth_pooled`; `baoab.py`/`baoab_device.py` are NOT edited and the overdamped step is
  unchanged — within the D3 authorization.
- **Overlay-only:** F_cap, n_FA, capacities are anchored to single-cell de-adhesion / FA counts,
  never fitted to the PI A/A0.
