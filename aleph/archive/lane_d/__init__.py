"""Lane D — vertical load-path diagnostic rigs (NON-AUTHORITATIVE DIAGNOSTIC).

This package is *lane-private assembly glue* for the FF-AC Lane D vertical:

  D0  membrane shell  <->  Biot-Darcy cytosol  (pressure traction + moving boundary)
  D1  nucleus  <->  MT / IF / LINC  (internal load-path transmission)

It **reuses the existing membrane / fluid / nucleus / MT / IF / LINC primitives
read-only** (`aleph.components.fluid`, `aleph.components.incumbent.membrane_pressure`,
`aleph.components.incumbent.fsi_coupling`, `aleph.laws.membrane_surface`,
`aleph.components.nucleus.envelope`, `aleph.components.solid.microtubule`,
`aleph.components.solid.intermediate_filament`) and adds only the assembly/adapter code
that wires them into a minimal scene, plus a few lane-local Warp glue kernels
(overdamped relaxation, interior volume-source perturbation).

HARD SCOPE NOTE (PI 2026-07-09/07-16): these rigs are small, coarse, diagnostic
scenes. They are **not** the native physiological population and **cannot** be
used for any production physics/parameter claim. They exist to *visualise and
gate the load-path coupling* before the connected whole-cell native run.

PORTED, NOT MERGED (PI decision D6, 2026-07-28). This code came from
``ac/mt-if-linc-vertical@d2f5ba95``, unmerged since 2026-07-22. The branch was
verified purely additive (24 files, +2422 lines, 0 deletions) so merging was
technically safe, but it self-labels NON-AUTHORITATIVE and merging would have put
its result artifacts in the tree where a later reader would quote them. So the
RIG and its GATES were ported and **its ``outputs/lane_d/`` results deliberately
were not**; they stay on the branch, registered under ``STATE.md`` (c) as not
quotable. The port was not a copy: the branch pinned a CUDA device ORDINAL in two
usage examples, which violates the no-hard-coded-device-ID rule — merging would
have carried that in silently.

Evidence status under the D1-B two-axis vocabulary: rung ``CUDA_UNIT`` at best
(the gates are CUDA-gated and skip on the dev Mac), ``QuantitativeClaim.BLOCKED``
unconditionally — the scenes are coarse by construction, so no magnitude from
them is quotable regardless of how green the gates run.
"""

__all__ = ["membrane_cytosol_rig", "force_path_rig"]
