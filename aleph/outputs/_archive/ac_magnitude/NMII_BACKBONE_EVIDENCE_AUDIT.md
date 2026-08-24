# NMII backbone persistence-length evidence audit

Date: 2026-07-21

Verdict: **no value is authorized for the mature MCF7 NMII bipolar-minifilament backbone.** The composed
runtime remains source-blocked.

## Contract-Graph/TAG result

The TAG query over SourceEvidence, KnowledgeClaim, ModelContract, and Parameter records returned no measured
or derivable NMII-minifilament backbone flexural rigidity or persistence length. Its only structural keyword
hit was `SE283 / Tam2021_Jour`, an unclassified simulation-method record whose motor is a zero-rest-length
Hookean spring. That is a lumped proxy and cannot satisfy the project's head-resolved mechanistic contract.

## Primary-source candidates

Two primary papers establish a useful but insufficient proxy chain:

1. Adamovic, Mijailovich & Karplus (2008), *Biophysical Journal*, DOI
   `10.1529/biophysj.107.122028`, estimates the bending elasticity of a **scallop myosin-II S2 coiled-coil
   subdomain** by molecular dynamics and normal-mode analysis. It reports a generic coiled-coil persistence
   range of 130–170 nm and a 60 nm S2 lateral stiffness around 0.008–0.012 pN/nm.
2. Kaufmann & Schwarz (2020), *PLOS Computational Biology*, DOI `10.1371/journal.pcbi.1007801`, uses
   `L_p = 130 nm` for bending of individual NM2 rods in a minifilament-assembly energy model.

Neither paper directly measures the effective bending rigidity of the **mature, multi-tail, bipolar NMIIA/B
minifilament backbone**, and neither is MCF7-specific. The 130 nm value therefore cannot be copied into the
single-chain mature-backbone primitive as a production parameter without an explicit mapping/ratification.

## A5000 falsification probe

The 130 nm single-coiled-coil value was run only as an explicitly non-production proxy through the existing
two-filament NG-1 diagnostic (`ticks=600`, 20 heads, RTX A5000). It **failed 4/6**:

- transmission ratio: 0.9180 (PASS under the existing diagnostic band);
- clamp reactions: -1.7447 and +1.6017 pN versus analytic 4.9627 pN (**both FAIL**);
- mean head load: -0.2490 pN; bound fraction: 0.7000;
- overall verdict: **FAIL 4/6**.

The older broad GAP sweep also flipped multiple rows (0.3 and 3.0 µm and the 10× arm-stiffness row achieved
only 5/6). Consequently, the earlier claim that NG-1 was insensitive to the persistence-length gap is false.
No gate, tolerance, population, or stiffness was changed in response.

## Required resolution

Before a physiological active-cell run, choose one of these evidence-backed routes and obtain PI ratification:

1. direct mature NMIIA/B minifilament bending data;
2. a structurally derived bundle stiffness from explicit tail count, packing, and slip/crosslink mechanics;
3. a more fine-grained explicit-tail backbone whose per-tail coiled-coil elasticity can legitimately use the
   single-rod evidence above.

Until then, `nmii_backbone_lp_um=None` and `BLOCKED_UNSOURCED_PHYSICAL_MAGNITUDE` are the correct runtime state.
