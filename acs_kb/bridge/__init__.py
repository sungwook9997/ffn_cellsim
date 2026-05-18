"""Bridge Unit 2.x — ECM ↔ cell mechanosensing (focal adhesions, motor-clutch).

Phase 1 Unit 2.1 deliverables:
    - FocalAdhesion dataclass (week-5 interface freeze; see ``types``).
    - Linear-elastic substrate stub (KU-1.21; see ``substrate_stub``).
    - Pereverzev catch-slip off-rate (KU-2.5; see ``catch_bond``).
    - Chan-Odde motor-clutch dynamics (KU-2.4; see ``motor_clutch``).
    - 1-D traction reducer (KU-2.12; see ``traction``).

Worker A's ECM (`acs_kb.ecm`) is not consumed yet; this Unit speaks to a
linear-elastic substrate stub. Replacement by an ECM adapter is deferred
to the next Unit (2.2 brief).

Source of truth: Notion `Simulation Knowledge Base` Unit 2 page
(364120da-ec5d-81ae-a3e5-c747e60524c6). If the implementation diverges
from a KU, the KU wins (CLAUDE.md hard rule).
"""
