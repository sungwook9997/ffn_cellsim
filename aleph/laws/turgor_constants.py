"""Resting osmotic turgor constants — extracted so the hot path does not import a retired module.

WHY THIS MODULE EXISTS.  These four constants lived in ``ff/gamma_floor.py``, and ``ff/network_warp.py``
— which owns ``link_spring_kernel``, ~11% of a native inner step — reached through that module to get
them.  ``gamma_floor`` is 684 lines whose own full-native results are RETIRED (``STATE.md`` (c) 2), and
its function-local ``ff.relax`` import chains on into the parked ``dcm/`` engine.  The result was that
every native cortex run had to rsync and build-verify eleven files and ~4,400 lines it never executes,
to reach four numbers.

Nothing about the numbers changed in the move; ``gamma_floor`` re-exports them, so its own uses and
every external caller are untouched.  The only altered field is the ``caller`` string recorded by
``resolve_pi0_pa`` — ``ff.gamma_floor`` becomes ``ff.turgor_constants`` — and that string appears in no
artifact (its only mention in the tree is a comment in ``ff/viz_cell.py``).

Sanity Gate:
    * dimensional: pressures are [pN/µm²]; 1 Pa == 1 pN/µm² in engine units, which is why the 40 Pa
      claim and ``TURGOR_DP0`` are numerically equal and must not be "converted".
    * boundary: ``resolve_pi0_pa`` raises rather than defaulting when the ledger gates Π₀ — the point of
      that gate is that a SILENT turgor is forbidden, and importing this module does not weaken it.
    * conservation/invariant: ``TURGOR_PI_IN0`` is the matched internal osmotic pressure and is ~4
      decades above ``TURGOR_DP0``; only their difference is the turgor, and treating the large value as
      a load is the failure mode the comment below guards.
    * numerical: R·T is evaluated once at import, not per call.
    * sign-sense: both pressures are positive-outward.
    * measurement-protocol: ``TURGOR_DP0`` is a CONVENIENCE claim, not a sourced MCF7 value, and ships
      with ``TURGOR_DP0_PROVENANCE`` so a caller records the classification next to any γ it measures.
"""

from __future__ import annotations

from aleph.laws.turgor_pi0 import PI0_CLAIM_HELA_PROXY_40PA, resolve_pi0_pa

# `TURGOR_DP0` below is NOT a default — it is the EXPLICIT, provenance-classified 40 Pa CONVENIENCE
# claim (HeLa proxy; PARAM_PROVENANCE_AUDIT_2026-07-24.md row 1, "the worst convenience default"),
# kept as a named constant so a caller that wants the historical value must NAME it and gets
# `TURGOR_DP0_PROVENANCE` to record alongside the number. Two other values are live in the tree
# (133 Pa band-implied/circular, ~72 Pa MCF7-geometry estimate) — see the ledger.
_PI0 = resolve_pi0_pa(PI0_CLAIM_HELA_PROXY_40PA, claim="claim_hela_proxy_40pa",
                      caller="ff.turgor_constants")
TURGOR_DP0 = _PI0.value_pa          # pN/µm² (=40 Pa); 1 Pa == 1 pN/µm² — CONVENIENCE, not sourced for MCF7
TURGOR_DP0_PROVENANCE = _PI0.as_artifact_record()   # serialize this next to any γ measured with it
OSMOLYTE_C_MM = 200.0      # mM internal impermeant osmolytes (Guo 2017; ~cytoplasmic salt)
VMIN_FRAC = 0.30           # Vmin/V0 osmotically-inactive volume fraction (Venkova 2022 / Adar 2025 / Guo 2017)
_RT_PN_UM2_PER_MM = 8.314 * 310.0  # R·T at 310 K, per (mol/m³); 1 mM=1 mol/m³, 1 Pa=1 pN/µm² → ×1 in FF units
# Internal osmotic pressure at rest Π_in0 = c·R·T [pN/µm²] (the huge matched value; only the net
# excess vs the medium = TURGOR_DP0 is the turgor). Drives the Guo stiffness, replacing the magic K_vol.
TURGOR_PI_IN0 = OSMOLYTE_C_MM * _RT_PN_UM2_PER_MM   # ≈ 5.15e5 pN/µm²
