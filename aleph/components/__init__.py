"""Cell parts and the connectors between them — the layer that owns state.

Each subpackage is one biological compartment with a disjoint population: cortex and lamellipodium
and SF-arc are separate owners, not one label-blind network. They couple **only** through explicit
connectors, never through shared nodes and never through a permanent weld. Co-location in an array
is not a connection, and this directory layout is the first place that has to be true.

Layer position, enforced by `tests/architecture/test_layer_directions.py`:

    units -> state -> laws -> engine -> components -> scripts

A component imports `laws/`; a law never imports a component. The same bending law is held by the
cortex, the stress fibres and the lamellipodium, and the moment it knows which one is holding it,
the other two either import that one or grow a copy. `Project_Aleph` lost this rule — its law layer
accumulated 53 imports of its component layer, one good reason at a time — and re-deriving it is
most of why the merge of 2026-08-09 kept this engine rather than that one.

`incumbent/` is the exception worth naming
------------------------------------------
It is the frozen whole-cell assembly that produced every native number this project can quote, and
`engine/` still binds it: 243 imports, of which the canonical layer is one. It is being strangled
out in gated stages and that work is **not finished**, which is why it sits here as a component
rather than in `archive/`. Nothing new goes into it.

Provenance: this package is the `ac/` ("Active Cell") tree of `ffn_cellsim`, reorganised by the
merge. Its build plan — five pillars, fluid-first poroelastic substrate, active-from-root NMII, one
emergent actomyosin network, deformable-mesh compartments, integrated coupling — is
`docs/v2_audit/NEW_ENGINE_BUILD_PLAN_2026-07-16.md`, and is retained because the pillars still
describe what these components are for.
"""
