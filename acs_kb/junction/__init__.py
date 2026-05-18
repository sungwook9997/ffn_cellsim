"""Junction Unit 4.x — cell-cell E-cadherin contacts, Young-equation contact angle.

Phase 1 Unit 4.1 deliverables (Worker D — Junction Track first work):
    - EcadherinJunction dataclass (week-9 interface freeze; see ``types``).
    - E-cadherin trans-bond Bell-Evans slip-only off-rate (KU-4.2; see
      ``cadherin``). Catch-bond detail is deferred to Phase 2.
    - Young-equation contact angle (KU-4.4; see ``contact_angle``).

Worker C's :class:`acs_kb.cell.cell.Cell` (week-6 freeze) is the **only**
import boundary for cell-side state: junction code touches cortex, focal
adhesions, or polarity only through Cell's public accessors
(``compute_cortex_boundary_position`` for contact-point lookup). Per the
Unit 4.1 brief, direct imports from ``acs_kb.bridge`` (Worker B) or
``acs_kb.ecm`` (Worker A) are forbidden — those layers are reached
transitively through Cell, never named here.

Source of truth: Notion `Simulation Knowledge Base` Unit 4 page
(KU-4.1 / KU-4.2 / KU-4.4 / KU-4.11 / KU-4.17). If the implementation
diverges from a KU, the KU wins (CLAUDE.md hard rule).
"""
