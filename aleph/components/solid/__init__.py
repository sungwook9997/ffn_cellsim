"""Filament-filament excluded volume (I2b) — the steric pillar of the Active Cell engine.

Cortex / SF / MT filaments occupy real space and cannot interpenetrate. The FF cortex LACKS this: today
only MT-tip<->cortex has a hand-wired steric pair (``ff/network_warp.py::soft_contact_kernel``), so the
cortical mesh has no real steric volume under compression/crowding — one of ENGINE_ARCHITECTURE_PLAN §3's
three defining living-cortex gaps, and the EV-off half of the CLAUDE.md hard worked-example
("Excluded volume: EV-off -> LJ repulsive ON from Phase 1"). I2b closes it.

Increment I2b (this module set) = a soft repulsive **WCA** (purely-repulsive Lennard-Jones) pair potential
over filament NODES via a device ``wp.HashGrid``:

    U_WCA(r) = 4 eps[(sigma/r)^12 - (sigma/r)^6] + eps   for r < r_c = 2^(1/6) sigma   (else 0)
    F(r)     = -dU/dr = (24 eps/r)[2(sigma/r)^12 - (sigma/r)^6] >= 0                    (repulsion, 0 at r_c)

exposed as the ``StericForce.accumulate(state, out_force)`` inner-loop force primitive (AC_PARALLEL_SESSIONS
§1.4). Field-independent — reads only node positions.

Modules:
  * ``wca_analytic``     — pure-NumPy WCA oracle (energy/force/stiffness, compressed-pair root). Ground truth.
  * ``steric_reference`` — pure-NumPy brute-force + cell-list EV force (the algorithm oracle the Warp
                           hash-grid kernel must reproduce). Also the CFL ``k_EV`` stiffness term.
  * ``steric_warp``      — the Warp-CUDA ``StericForce`` + hash-grid WCA kernel (authored source; the lead
                           gates it natively against ``steric_reference`` on the gbook A5000).

Parameters (``params_i0b2b.yaml``): ``sigma`` (steric diameter — physical, filament ~7 nm floor +
coarse-graining note) and ``k_EV`` (contact stiffness — a NUMERICAL repulsion scale, Magic-Number Block,
grid-invariant, NEVER tuned to a crowding outcome; NOT a literature value).

Sanity Gate (before first native execution, per the Sanity Gate Protocol — full detail in wca_analytic /
steric_reference / the module tests):
  * dimensional analysis: sigma [um], epsilon [pN.um], k_EV [pN/um], F [pN]; units close;
  * boundary cases: zero-overlap (r >= r_c) == zero-force == OFF bit-identical; a pair at r_c is C1
    (energy AND force -> 0); the LJ minimum sits exactly at r_c;
  * conservation: net internal force ~ 0 (Newton's third law, own-row + reaction);
  * numerical: F = -dU/dr matches a finite-difference gradient (sign arbiter); k_EV enters the CFL kmax;
    hash-grid == brute-force to summation tolerance (neighbour cell size is a grid-invariant accelerator);
  * sign-sense: F > 0 repulsive for r < r_c; a compressed pair settles at a FINITE r_eq > 0 (never collapses);
  * measurement consistency: WCA vs Weeks-Chandler-Andersen 1971 closed form.

Runtime: NVIDIA Warp on CUDA GPU only (I0-A). No HOOMD, no CPU simulation path.
"""
