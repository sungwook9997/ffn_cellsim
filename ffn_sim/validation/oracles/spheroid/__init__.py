"""Closed-form acceptance oracles for the Layer-2 multicellular spheroid model.

Runtime-import-forbidden (enforced by ``ffn_sim.common.integrity``): these modules
hold only closed-form *reference* math (the analytic answer the simulation is checked
against). They are imported solely from ``tests/`` and ``validation/``.

Contents:
- ``aa0_law``: the MCF7 spreading law ``A/A0 = a + b/R + c/R^2`` + least-squares fit +
  term-dominance decomposition + analytic geometry references for the runtime
  ``ffn_sim.archive.hoomd_legacy.spheroid.observables`` measurement protocol.
"""
