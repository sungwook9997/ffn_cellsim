# ActiveCellSim

**Active poroelastic spheroid spreading on adhesive substrate — a drying-thin-film inspired cellular hydrodynamics framework**

## Quick Summary
3D simulation of MCF7 spheroid spreading on Col1-coated dish, built independently from first principles to provide a new mechanobiological lens (NOT to fit existing data).

## Key Differentiators
1. **Hybrid Eulerian-Lagrangian MPM** with cell-equivalent material points (5-layer integrated model)
2. **Drying-thin-film inspired** — cellular Marangoni, coffee-ring analog, drying-induced concentration via mechano-osmotic coupling
3. **Adhesion network dynamics** (φ ODE) — E-cadherin ↔ Integrin-β1/Laminin transition, dynamically modulating mechanics
4. **Radial vs Full 3D anisotropic** parallel simulations — quantifying the validity domain of radial-symmetry approximations
5. **First-principles validation** from IF≥15 literature, not parameter fitting

## Status
🚧 Planning complete, implementation begins.

## Hardware
Primary: NVIDIA RTX A5000 + Xeon W-11955M workstation (Windows, headless SSH)
Portable: RTX 4090×2 lab workstation, Google Colab fallback

## Workflow Roles
| Tool | Role |
|------|------|
| Claude Desktop App | Brain trust — physics design, model review, paper interpretation, result analysis |
| VSCode + Claude Code (extension) | Code editing, inline diff, file modification |
| Claude Code Terminal | Execution, environment setup, SSH commands, GPU monitoring, batch runs |
| Mac (M1) | Lightweight viz post-processing, interactive review |
| Windows (A5000) | Heavy simulation runs, Blender renders |

## Getting Started (for Claude Code)
1. Read `CLAUDE.md` (this is automatic).
2. Read `STRUCTURE.md` for the file→version map (v1 frozen / v2 active / shared).
3. Read `docs/v2/10_dev_roadmap_v2.md` to understand next milestone (v2 active). The v1 roadmap is at `docs/v1/10_dev_roadmap.md` for historical reference.
4. Read `docs/00_project_vision.md` and `docs/v2/00_project_vision_v2.md` for full framing.
5. Pull other docs only as relevant to the current task.

## Project Lead
PI: Sungwook Yoon (sungwook999@kaist.ac.kr)
Lab: Shin Lab, Dept. of Mechanical Engineering, KAIST
Corresponding: Jennifer H. Shin (j_shin@kaist.ac.kr)
