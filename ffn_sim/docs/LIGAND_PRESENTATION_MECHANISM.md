# Ligand presentation mechanism — Bare / Pre / Lam4 (the availability axis)

> **STATUS: DESIGN doc — 2026-06-02. Prep, NOT a contract, NO code.** Parallel session
> (`AUTONOMOUS_LOG_2026-06-02.md`). Closes the last open TODO in the PI-exp Layer-1
> mapping table (`PI_EXP_VALIDATION_MAP.md:63` — *"surface-bound (Pre) vs soluble (Lam4)
> … TODO: mechanism choice"*). Uses ONLY anchors already in the PI-exp map (no new web
> sweep — the map says do not burn one). Orthogonal to `LIGAND_IDENTITY_FA_SPEC.md`. No
> edit to `bridge/fa.py`.

## Two orthogonal axes — keep them separate

The PI experiment varies **two independent things**; conflating them is the trap:

| Axis | What changes | Platform knob | Doc |
|---|---|---|---|
| **Identity** | which ligand → which integrin → **per-bond off-rate** | `PereverzevParams` per species | `LIGAND_IDENTITY_FA_SPEC.md` (done) |
| **Presentation** | how the ligand is offered → **availability / on-rate / spatial density** | FA ligand density, `k_on`, FA placement | **this doc** |

Bare / Pre / Lam4 are mostly a **presentation** difference, NOT a per-bond affinity
difference. So model them on the availability axis, leaving the off-rate to the identity
axis. (Bare/Pre = col-I identity; Lam4 = col-I coating + soluble laminin-111 → adds the
α6β1 identity AND a soluble-presentation change.)

## The mechanism choice (affinity vs avidity vs availability)

PI-exp map §Layer-1 affinity-vs-avidity row (med conf, by analogy — Irvine 2002;
Partridge/Luo 2005 PNAS; Steiger 2024; cilengitide agonist@1 nM / antagonist@1 µM):

- **Affinity** (per-bond K_D): TM-disruption raises *monomeric* affinity. → the
  identity axis already carries per-bond kinetics; do NOT re-encode presentation here.
- **Avidity** (clustering / multivalency): clustering raises occupancy **only under
  cooperative binding**. → a clustering/spatial term, not a single-bond change.
- **Availability** (ligand areal density + on-rate / accessibility): the cleanest,
  least-contested platform handle — the FA already has a ligand density + capture radius
  + `k_on`. **RECOMMENDED primary mechanism.**

**Recommendation:** model Bare / Pre / Lam4 as an **availability** difference
(ligand areal density × `k_on` × spatial distribution), with **avidity (clustering) as
a documented optional refinement**. Do NOT model it as a per-bond affinity change.

| Condition | β1 IF (observable) | availability mapping | spatial (→ clutch_spatial) |
|---|---|---|---|
| **Bare** | weak / diffuse | **low** ligand density / `k_on` | diffuse, low `n_engaged` |
| **Pre** | peripheral | **higher** surface-adsorbed density | peripheral (high `f_edge`) |
| **Lam4** | uniform | col-I density + **soluble α6β1 availability** | uniform (low `angular_cv`) |

The β1-IF pattern is exactly the `bridge/clutch_spatial.py` observable — so this axis is
**falsifiable** against the existing β1-distribution metrics.

## Caveats (honest uncertainty — med confidence)

- The presentation mechanism is **by analogy** (αIIbβ3/αvβ3/RGD), extrapolated to
  laminin (PI-exp map). Tag it proxy/med, not high.
- ⚠️ **REFUTED (do NOT use):** "soluble-FN 5–6× adhesion drop" (1-2). Do not encode a
  large soluble-presentation penalty from it.
- **Single-cell vs collective flag** (PI-exp open item): single-cell laminin traction is
  *lower* (lit) while the Lam4 *collective* effect is *higher* (poster). This doc is the
  **single-cell presentation** mapping — reproduce the lower single-cell laminin
  availability/traction; the Lam4 collective elevation is a **Layer-2** (multi-cell)
  question, not a single-cell presentation knob. Keep the layers separate.
- iCVD pV4D4 adsorption density (the physical basis of Pre availability) is **unfilled**
  — route to the Im lab (PI-exp map), not a sweep. Until then `k_on`/density for Pre is
  a *relative* setting (Pre > Bare), not an absolute.

## Sanity Gate (recorded now per CLAUDE.md)

1. **Dimensional.** ligand areal density [1/m²]; `k_on` [s⁻¹]; availability dimensionless.
2. **Boundary.** Bare = the low-availability baseline; setting Pre/Lam4 density → Bare
   density recovers Bare (monotone, no discontinuity).
3. **Orthogonality.** Changing presentation (density/`k_on`) must NOT change the per-bond
   off-rate (identity axis) — the two knobs are independent; assert the off-rate is
   unchanged across Bare/Pre when the ligand identity is the same (col-I).
4. **Falsifiability.** The predicted β1 pattern (Bare diffuse / Pre peripheral / Lam4
   uniform) must register on `clutch_spatial` metrics — the acceptance is the *ordering*,
   not absolutes (qualitative IF).
5. **Sign-sense.** Higher availability ⇒ more engaged clutch ⇒ higher traction (within a
   ligand identity). Wrong sign = bug.

## Open items for PI / Lead

- [ ] Ratify **availability** (density × `k_on` × spatial) as the primary Bare/Pre/Lam4
      mechanism, with avidity/clustering as an optional refinement.
- [ ] Confirm presentation stays **relative** (Pre > Bare) until the Im-lab pV4D4
      adsorption density lands (no absolute magic number).
- [ ] Confirm the single-cell↔Layer-2 split for the laminin traction direction.
- [ ] Sequence with `LIGAND_IDENTITY_FA_SPEC` (same `fa.py` extension; identity = off-
      rate, presentation = density/`k_on`).
