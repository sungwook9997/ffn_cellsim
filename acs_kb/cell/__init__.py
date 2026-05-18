"""Cell Unit 3.x — cortex, force balance, Cell dataclass, visualization.

Phase 1 Unit 3.1 ships:
    - discrete bead-spring cortex (circular topology, KU-3.1 / KU-3.17)
    - overdamped quasi-static force balance (KU-3.16)
    - Cell dataclass (week-6 frozen interface for Unit 4 junctions)
    - cortex visualisation + rounding animation

Import submodules directly:
    from acs_kb.cell.cortex import generate_cortex, Cortex
    from acs_kb.cell.force_balance import solve_overdamped_step, relax
    from acs_kb.cell.cell import Cell
"""
