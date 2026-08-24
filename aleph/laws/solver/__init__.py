"""FF single-cell solver core (FF_SINGLE_CELL_OPTIMIZATION_PLAN §4).

Reusable, GPU-resident implicit-solver building blocks shared by the AFM / crawl / resting / ECM workloads
(and, later, by the fleet engine — one implementation, not two). Built SC0-first: the native profile
(`outputs/ff_single_opt/sc0_profiles/REPORT.md`) showed the implicit step is CG-DOMINATED (86-89%), so the
first module here is the CG accelerator, not a matrix-free operator.
"""
