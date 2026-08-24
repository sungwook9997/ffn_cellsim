"""Where a number came from — the four-class vocabulary, at the bottom of the layer stack.

`PROVENANCE` labels the ORIGIN of a magnitude, not its value: whether it was read from a source,
derived from sourced quantities, chosen for convenience, or is an open gap the PI must close.

Why it lives here and not where it was written. It was defined in
``aleph.engine.forces_manifest`` and imported from there by ``aleph.world.bond``, which made the
arena — CANONICAL since 2026-08-20 — import the port source it replaces. The dependency was only a
drift guard over four strings, but direction is not about weight: an arena that reaches into
``engine/`` has wrapped it rather than replaced it, and ``test_layer_directions`` caught exactly
that on the promotion's first run. ``aleph/units/`` is where it belonged anyway — that package's own
docstring says it "holds the provenance side of numbers ... and where the bound comes from", and its
layer position is the bottom, so everything may read it and it may read nothing.

``forces_manifest`` re-exports the name, so nothing that imported it from there breaks.
"""

from __future__ import annotations

#: Provenance classes, in decreasing strength of claim.
#:
#: * ``SOURCED`` — read from a cited source, and the citation supports THIS quantity. ⚠ The second
#:   half is not automatic: the NMII per-head force audit of 2026-08-20 found a real citation
#:   (Kovács 2003) attached to a number that paper does not contain, because it measures kinetics
#:   and not force.
#: * ``DERIVED`` — computed from sourced quantities by a stated relation.
#: * ``CONVENIENCE`` — chosen to make something run. Not evidence, and must not be quoted as one.
#: * ``PI_GAP`` — no value exists yet; the build refuses rather than defaulting.
PROVENANCE: tuple[str, ...] = ("SOURCED", "DERIVED", "CONVENIENCE", "PI_GAP")

__all__ = ["PROVENANCE"]
