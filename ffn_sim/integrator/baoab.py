"""Leimkuhler-Matthews BAOAB-limit overdamped Langevin integrator (D3).

PHASE_0_3_DECISIONS D3 (2026-05-19, PI-ratified): Phase 1 canonical
integrator. Replaces HOOMD's vanilla ``md.methods.Brownian`` (Euler-Maruyama,
O(Δt)) with the Leimkuhler-Matthews BAOAB-limit step (O(Δt²) on harmonic
systems), matching AFINES ``filament::update_positions``
(filament.cpp:200-230, comment references "Leimkuhler 2013").

Update rule per particle::

    r(t+Δt) = r(t) + (F(t)/γ)·Δt + √(kT/(2 γ Δt)) · (W_n + W_{n-1}) · Δt

where ``F(t)`` is the net conservative force at the current configuration,
``γ`` is the per-particle Stokes drag ``γ = 6π η R_bead`` (D4: depends on
bead radius, not bond rest length ℓ₀), and ``W_n``, ``W_{n-1}`` are
independent unit Gaussian draws at the current and previous step. The
"BAOAB-limit" form is the high-friction/large-mass limit of the BAOAB
splitting; averaging two Gaussians is what makes it Leimkuhler-Matthews
rather than Euler-Maruyama, giving "minimal canonical deviations" on
harmonic systems with stiff springs (Leimkuhler & Matthews 2013,
*Appl. Math. Res. eXpress* 2013(1):34-56).

Shared by H.1, H.2, H.3, H.4, H.5 — first physics module of Phase 1, so
the Sanity Gate below is the contract every downstream unit inherits.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule "Sanity Gate Protocol mandatory before first
execution of any physics/numerics module." This block is the validation
contract; checks marked ``RUNTIME`` are asserted in
``LeimkuhlerMatthewsBAOAB.act`` / ``attach``, ``STATIC`` checks live in
``ffn_sim/tests/test_baoab.py``, and ``EMPIRICAL`` checks live in
``ffn_sim/scripts/baoab_polymer_sanity.py``.*

1. **Dimensional analysis**

   - ``F/γ``: [N] / [N·s/m] = [m/s]. Multiplied by Δt [s] → [m]. ✓
   - ``√(kT/(2 γ Δt))``: [√(J / (N·s/m · s))] = [√(J·m/(N·s²))]
     = [√(N·m·m/(N·s²))] = [m/s]. Multiplied by Δt → [m]. ✓
   - Both increment terms are lengths; SI is internally consistent.
     STATIC test asserts the unit balance symbolically using ``sympy``
     where the prefactor ``bd_prefactor = √(kT/(2 γ Δt))`` is recombined
     with Δt and confirmed dimension ``[m]``.

2. **Boundary cases**

   - ``Δt → 0``: position increment → 0 (deterministic part Δt-linear;
     stochastic part ``√(1/Δt)·Δt = √Δt`` — vanishes correctly). ✓
   - ``γ → ∞`` (infinite drag): both terms → 0. Particle freezes. ✓
   - ``kT → 0``: stochastic term → 0; pure gradient descent in (F/γ)·Δt.
     STATIC test: with kT=0 and F=0, positions must not change.
   - ``F = 0, kT > 0``: pure Brownian diffusion, ⟨Δr²⟩ → 2 d D Δt with
     D = kT/γ (Stokes-Einstein) over long time; this is the
     EMPIRICAL check in §6.
   - First step ``prv_rnds = 0``: equivalent to Euler-Maruyama with
     halved noise, then converges to L-M from step 2. Matches AFINES
     post-fracture behaviour and the brief's reference.

3. **Conservation invariants**

   - **Linear momentum is NOT conserved** under stochastic noise; this
     is correct for Langevin (heat bath supplies impulse). EMPIRICAL
     check: in equilibrium without external bias, ``⟨p⟩`` should
     hover near zero with variance ``∝ kT γ Δt`` per particle. We
     remove COM drift only if explicitly enabled (off by default; the
     drift IS a physical signature of finite-sample noise).
   - **Energy is NOT conserved**; equipartition is the relevant
     thermodynamic invariant: ⟨½ k Δx²⟩ = ½ kT per harmonic mode.
     EMPIRICAL gate (§6) on a 100-bead harmonic polymer.
   - **No drift in field-free, force-free system**: ⟨r(t)⟩ should
     remain at the initial mean to within 3σ of √(2 d D t / N).
     EMPIRICAL.

4. **Numerical sanity**

   - **CFL / stability**: explicit overdamped Euler is stable for
     Δt < 2/(k_max/γ) where k_max is the stiffest spring. L-M shares
     this bound; we require ``Δt ≤ cfl_safety · τ_min`` with
     ``τ_min = min(γ_b/k_xl, γ_b ℓ₀/μ, γ_b ℓ₀³/κ)`` and
     ``cfl_safety = 0.1`` (KU-1.26 conventional). RUNTIME check: at
     ``attach`` time, the Updater queries ``sim.operations.integrator.dt``
     and asserts ``dt > 0`` and finite. It does NOT compute τ_min
     itself — that is the caller's contract via
     ``ffn_sim/validation/oracles/common/sanity_gate.py::gate_unit1_2_dynamics``.
   - **Per-particle Stokes drag**: γ must be > 0 and finite for every
     particle type in the filter. RUNTIME check at ``attach``.
   - **NaN/Inf force guard**: AFINES ``bead::update_force``
     (bead.cpp:54-63) aborts on NaN/Inf force. We replicate: RUNTIME
     check inside ``act`` raises ``FloatingPointError`` if any
     ``net_force`` or any newly-written position is non-finite,
     including the timestep so PI can locate the divergence.
   - **Precision**: HOOMD 7 defaults to float64 on CPU; we read
     ``net_force`` and write ``position`` in that dtype. STATIC test
     asserts ``snapshot.particles.position.dtype == np.float64``.

5. **Sign-sense**

   - Gradient descent: ``F = −∇U``. Positive bond stretch → F pulling
     toward rest → r moves toward rest along ``+F/γ·Δt``. STATIC test:
     a single 2-bead harmonic pair stretched by +δ relaxes monotonically
     in mean (deterministic limit kT=0) for at least 100 Δt.
   - Noise sign convention: the same ``W_n`` is used for both deterministic
     and stochastic terms; sign of the diffusion is symmetric so this is
     trivially satisfied. STATIC: ``np.mean(W)`` over N=1e6 draws
     within 5σ of 0.

6. **Measurement-protocol consistency**

   - **Equipartition gate on 100-bead harmonic polymer** (the brief's
     freeze-point validation). Drive the existing
     ``ffn_sim/scripts/hoomd_polymer_sanity.py`` 100-bead chain
     (k=100 ε/σ², kT=1.0, γ=1.0, dt=0.005) with this Updater for
     N_steps ≥ 5e5 after a 5e4-step burn-in. Per-bond
     ``⟨½ k (|r_{i+1}−r_i| − r0)²⟩`` must equal ½ kT to ±5% (block
     average over 10 blocks, last 4.5e5 steps). For the bending
     mode, ``⟨½ k_θ (θ−π)²⟩`` must equal ½ kT to ±10% (looser
     because the harmonic-bend small-angle approximation breaks at
     finite amplitude — same tolerance AFINES quotes in
     §7 of the algorithm notes).
   - **Diffusion gate**: single isolated free bead (k=0, k_θ=0)
     under L-M must give ⟨Δr²(t)⟩/(2 d t) → kT/γ within ±10% for
     t ≥ 1000 Δt (Stokes-Einstein). EMPIRICAL.
   - **L-M order check (D3 deliverable)**: equipartition convergence
     rate vs Δt halving should improve as O(Δt²) on the harmonic
     polymer — i.e., relative error in ⟨U_bond⟩ at Δt/2 should be
     ≤ ¼ of that at Δt. EMPIRICAL. This is the test that distinguishes
     L-M from Euler-Maruyama and the explicit reason D3 picks L-M.
   - **prv_rnds memory**: STATIC test verifies that step n's W_n is
     written into the buffer that step n+1 reads as W_{n-1}, by
     instrumenting the Updater with a deterministic RNG and checking
     the buffer state after two calls.

Failure handling
----------------
If any RUNTIME check fails, ``act`` raises and the Simulation halts.
If a STATIC or EMPIRICAL gate fails, the test prints the offending
quantity and the bound; per CLAUDE.md "no gate-loosening", we surface
to PI rather than relax the tolerance.

Implementation note (HOOMD 7 plumbing)
--------------------------------------
HOOMD 7's ``md.Integrator.methods`` is where position updates live;
``md.methods.Brownian`` is Euler-Maruyama and is exactly what D3
forbids. The cleanest fit for an external position step in HOOMD 7 is
``hoomd.custom.Action`` wrapped in ``hoomd.update.CustomUpdater``,
co-existing with an ``md.Integrator`` whose role is *only* force
evaluation (``forces=[...]``, ``methods=[]``). The Updater runs every
step before the force evaluator's next ``methods.step``, reads
``cpu_local_snapshot.particles.net_force`` (forces from the prior
step's evaluation, i.e. ``F(r(t))``), applies the L-M step, and writes
positions in place. An initial ``sim.run(0)`` (force-eval-only) seeds
``net_force`` before the first integration step. This is the H.1 brief's
"Open implementation question" — benchmarked here on the 100-bead
polymer before scaling to ECM.

References
----------
- Leimkuhler & Matthews 2013, "Rational construction of stochastic
  numerical methods for molecular sampling", *Appl. Math. Res. eXpress*
  2013(1):34-56. BAOAB and its overdamped limit.
- AFINES ``filament.cpp::update_positions`` lines 200-230 (mechanism),
  ``filament.cpp:78`` and ``filament_ensemble.cpp:466`` ("Leimkuhler 2013"
  comments).
- PHASE_0_3_DECISIONS.md §D3.
- ``ffn_sim/docs/briefs/H1_ecm_mikado.md`` §Integration.
- ``ffn_sim/docs/AFINES_ALGORITHM_NOTES.md`` §2.1 (force pipeline) and
  §7 (polymer sanity equipartition expectations).
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

import hoomd
import hoomd.custom

# Optional GPU compute path: ported 2026-05-28 (PI ratified option (i) in
# FIXMAN_LAZY_EVAL_DESIGN.md). On a CPU device this module behaves bit-for-bit
# as before (cupy is not imported). On a GPU device cupy is required and the
# Action runs entirely against HOOMD's gpu_local_snapshot — eliminates the
# per-step GPU↔host snapshot copy that bottlenecked KU-3.5 (3.5 ms/step CPU≈GPU).
try:  # noqa: SIM105
    import cupy as _cp  # type: ignore
    _CUPY_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only on hosts without cupy
    _cp = None
    _CUPY_AVAILABLE = False


class LeimkuhlerMatthewsBAOAB(hoomd.custom.Action):
    """Custom HOOMD ``Action`` implementing the L-M BAOAB-limit step (D3).

    Co-exists with an ``md.Integrator`` carrying ``forces=[...]`` and
    ``methods=[]``: HOOMD evaluates forces every step, this Action reads
    ``cpu_local_snapshot.particles.net_force`` and advances positions per
    the formula in the module docstring. The simulation must not also
    attach an ``md.methods.*`` position updater — that would double-step.

    Parameters
    ----------
    kT : float
        Thermal energy in HOOMD units (matches force units · length).
    gamma : Mapping[str, float]
        Per-particle-type Stokes drag γ_b = 6π η R_bead. Must cover every
        particle type present in the state. γ must be > 0 and finite.
    dt : float
        Integration timestep. Required to equal
        ``sim.operations.integrator.dt``; the Action asserts equality at
        ``attach`` time.
    seed : int
        Seed for the per-Action ``np.random.default_rng``. Independent of
        the HOOMD simulation seed so this Action is RNG-isolated.

    Notes
    -----
    See the module docstring §Sanity Gate for the dimensional / boundary
    / conservation / numerical / sign / measurement checks this
    implementation must satisfy. RUNTIME checks are enforced here; STATIC
    and EMPIRICAL checks live in ``tests/test_baoab.py`` and
    ``scripts/baoab_polymer_sanity.py``.
    """

    def __init__(
        self,
        *,
        kT: float,
        gamma: Mapping[str, float],
        dt: float,
        seed: int = 0,
    ) -> None:
        super().__init__()
        if not (np.isfinite(kT) and kT >= 0.0):
            raise ValueError(f"kT must be finite and ≥ 0, got {kT!r}")
        if not (np.isfinite(dt) and dt > 0.0):
            raise ValueError(f"dt must be finite and > 0, got {dt!r}")
        for typ, g in gamma.items():
            if not (np.isfinite(g) and g > 0.0):
                raise ValueError(
                    f"gamma['{typ}'] must be finite and > 0, got {g!r}"
                )
        self.kT = float(kT)
        self.dt = float(dt)
        self.gamma_map = dict(gamma)
        self._seed = int(seed)
        self._rng = np.random.default_rng(seed)

        # Device-specific dispatch: defaults to numpy/CPU. attach() rewrites
        # _xp/_asarray/_snap_ctx_attr/_rng if the simulation device is a GPU.
        # CPU path is bit-for-bit unchanged from pre-port behaviour.
        self._xp = np
        self._asarray = np.asarray
        self._snap_ctx_attr: str = "cpu_local_snapshot"

        # All per-particle state below is indexed by HOOMD particle TAG
        # (stable identity), not by snapshot row, because HOOMD's
        # ParticleSorter tuner reorders local-snapshot rows between calls.
        # ``_gamma_by_tag`` and ``_prv_rnds`` are allocated at ``attach()``
        # and looked up via ``snap.particles.tag`` each step.
        self._gamma_by_tag = None
        self._bd_prefactor_by_tag = None
        self._prv_rnds = None  # shape (N_tags, 3) on _xp module
        self._sim_ref: hoomd.Simulation | None = None
        self._steps_run: int = 0

    # ------------------------------------------------------------------
    # HOOMD Action lifecycle
    # ------------------------------------------------------------------
    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        """Allocate per-particle buffers; sanity-check integrator wiring."""
        super().attach(simulation)
        self._sim_ref = simulation

        ig = simulation.operations.integrator
        if ig is None:
            raise RuntimeError(
                "LeimkuhlerMatthewsBAOAB requires an md.Integrator with "
                "forces attached (methods=[]). None is set."
            )
        if not np.isclose(float(ig.dt), self.dt, rtol=0, atol=0):
            raise RuntimeError(
                f"Integrator dt={float(ig.dt)!r} must equal Action dt="
                f"{self.dt!r}; the L-M step assumes a fixed dt."
            )
        if len(ig.methods) != 0:
            raise RuntimeError(
                "md.Integrator.methods must be empty; a position-updating "
                "Method would double-step alongside the L-M Action. Got: "
                f"{list(ig.methods)!r}"
            )

        # Device-specific dispatch. The cpu_local_snapshot init code below
        # is intentionally device-agnostic: we read tags/typeid on host
        # (one-shot, cheap) and then transfer the resulting tag-indexed
        # arrays to the target device.
        if isinstance(simulation.device, hoomd.device.GPU):
            if not _CUPY_AVAILABLE:
                raise RuntimeError(
                    "LeimkuhlerMatthewsBAOAB attached to a GPU device but cupy "
                    "is not importable. Install with "
                    "`pip install cupy-cuda12x[ctk]` (CUDA 12 build) and retry."
                )
            self._xp = _cp
            self._asarray = _cp.asarray
            self._snap_ctx_attr = "gpu_local_snapshot"
            self._rng = _cp.random.default_rng(self._seed)

        with simulation.state.cpu_local_snapshot as snap:
            type_ids = np.asarray(snap.particles.typeid)
            tags = np.asarray(snap.particles.tag)
            type_names = list(simulation.state.particle_types)
            N = type_ids.shape[0]

            missing = [t for t in type_names if t not in self.gamma_map]
            if missing:
                raise RuntimeError(
                    f"gamma is missing entries for particle types {missing}; "
                    f"state has types {type_names}."
                )

            # Tag-indexed γ: row r holds the drag for particle whose tag==r.
            # Phase 1 has a single rank, so tags are dense in [0, N).
            if N == 0:
                raise RuntimeError("L-M Action attached to empty state.")
            if int(tags.max()) >= N or int(tags.min()) < 0:
                raise RuntimeError(
                    "L-M Action requires dense particle tags in [0, N); "
                    f"got tag range [{int(tags.min())}, {int(tags.max())}] "
                    f"for N={N}."
                )
            gamma_by_typeid = np.array(
                [self.gamma_map[t] for t in type_names], dtype=np.float64
            )
            gamma_by_tag = np.empty(N, dtype=np.float64)
            gamma_by_tag[tags] = gamma_by_typeid[type_ids]

        # bd_prefactor = √(kT / (2 γ Δt))  (per AFINES filament.cpp:78)
        # Multiplied later by (W_n + W_{n-1}) · Δt to get a displacement.
        bd_prefactor_host = np.sqrt(
            self.kT / (2.0 * gamma_by_tag * self.dt)
        ).reshape(-1, 1)
        # Move tag-indexed buffers to the target device (no-op on CPU).
        self._gamma_by_tag = self._asarray(gamma_by_tag)
        self._bd_prefactor_by_tag = self._asarray(bd_prefactor_host)
        # prv_rnds initialised to zero → first step is effectively
        # Euler-Maruyama with halved noise; converges to L-M from step 2.
        # (Matches AFINES post-fracture initialisation behaviour.)
        self._prv_rnds = self._xp.zeros((N, 3), dtype=np.float64)
        self._steps_run = 0

    def act(self, timestep: int) -> None:
        """Apply one L-M BAOAB-limit step to all filtered particles."""
        if self._prv_rnds is None or self._gamma_by_tag is None:
            raise RuntimeError(
                "LeimkuhlerMatthewsBAOAB.act called before attach(); "
                "wrap this Action in hoomd.update.CustomUpdater and add "
                "it to sim.operations.updaters."
            )

        sim = self._sim_ref
        assert sim is not None  # set in attach()
        xp = self._xp
        snap_ctx = getattr(sim.state, self._snap_ctx_attr)

        with snap_ctx as snap:
            pos = self._asarray(snap.particles.position)      # (N, 3) rw
            F = self._asarray(snap.particles.net_force)       # (N, 3) ro
            image = self._asarray(snap.particles.image)       # (N, 3) rw
            tag = self._asarray(snap.particles.tag)           # (N,) stable id
            N = pos.shape[0]

            if N != self._prv_rnds.shape[0]:
                raise RuntimeError(
                    "Particle count changed between attach() and act() "
                    f"(was {self._prv_rnds.shape[0]}, now {N}); "
                    "topology mutation requires a fresh attach()."
                )

            if not bool(xp.all(xp.isfinite(F))):
                # Sync to host for the error message — only on the failure path.
                F_host = _to_host(F, xp)
                bad = np.argwhere(~np.isfinite(F_host))
                raise FloatingPointError(
                    f"Non-finite net_force at timestep={timestep}; "
                    f"first offending (particle, dim) entries: "
                    f"{bad[:5].tolist()}. AFINES analogue: "
                    "bead::update_force abort (bead.cpp:54-63)."
                )

            # Gather tag-indexed per-particle state into the current row
            # ordering (HOOMD's ParticleSorter reorders rows between calls).
            inv_gamma_row = (1.0 / self._gamma_by_tag[tag]).reshape(-1, 1)
            bd_prefactor_row = self._bd_prefactor_by_tag[tag]
            prv_W_row = self._prv_rnds[tag]

            W_row = self._rng.standard_normal(size=(N, 3))

            # r += (F/γ)·Δt + √(kT/(2 γ Δt))·(W_n + W_{n-1})·Δt
            dr = F * inv_gamma_row * self.dt + (
                bd_prefactor_row * (W_row + prv_W_row) * self.dt
            )
            new_pos = pos + dr

            if not bool(xp.all(xp.isfinite(new_pos))):
                new_pos_host = _to_host(new_pos, xp)
                bad = np.argwhere(~np.isfinite(new_pos_host))
                raise FloatingPointError(
                    f"Non-finite position after L-M step at timestep={timestep}; "
                    f"first offending entries: {bad[:5].tolist()}."
                )

            # Minimum-image wrap into the (possibly Lees-Edwards-sheared)
            # simulation box, updating image flags. HOOMD 7's Python Box
            # exposes no public ``wrap``; we implement the upper-triangular
            # fractional-coord wrap directly so xy-tilt (used by KU-1.30 #2
            # strain stiffening) is handled correctly.
            box = sim.state.box
            wrapped, img_delta = _wrap_into_box(new_pos, box, xp=xp)
            pos[:] = wrapped
            image[:] = image + img_delta

            # Scatter W_row → tag-indexed prv_rnds buffer for next step.
            self._prv_rnds[tag] = W_row

        self._steps_run += 1

    # ------------------------------------------------------------------
    # Introspection (used by tests)
    # ------------------------------------------------------------------
    @property
    def prv_rnds(self) -> np.ndarray | None:
        """Read-only view of the persisted W_{n-1} buffer (None pre-attach).

        On a GPU device the buffer lives on the GPU as a cupy array; this
        property returns a host-side numpy COPY (cupy.ndarray.get()), so
        flipping writeable=False on the returned view does NOT leak back
        into the live integrator buffer.
        """
        if self._prv_rnds is None:
            return None
        if self._xp is np:
            v = self._prv_rnds.view()  # numpy view: independent writeable flag
        else:
            v = self._prv_rnds.get()   # cupy → numpy host copy (always disjoint)
        v.flags.writeable = False
        return v

    @property
    def steps_run(self) -> int:
        return self._steps_run


def _to_host(arr, xp):
    """Return a numpy copy of ``arr`` regardless of whether it is on GPU.

    Used only on error / introspection paths so the per-step hot loop stays
    fully device-resident. ``xp`` is the module the caller's other arrays
    live in (numpy or cupy); cupy arrays expose ``.get()``.
    """
    if xp is np:
        return np.asarray(arr)
    return arr.get()


def _wrap_into_box(pos, box: hoomd.box.Box, *, xp=np):
    """Wrap positions into the box; return (wrapped_pos, image_delta).

    HOOMD's upper-triangular box matrix is

        M = [[Lx, xy·Ly, xz·Lz],
             [0,  Ly,    yz·Lz],
             [0,  0,     Lz   ]]

    Fractional coordinates ``f = M⁻¹ r`` are wrapped to ``[-½, ½)`` and
    rounded; the integer round counts go into ``image_delta``. Returns
    new positions ``M · f_wrapped`` (still in float64) and the per-axis
    image-flag increments as ``int32``.

    ``xp`` is the array module (numpy or cupy) ``pos`` lives in. Box
    geometry (Lx, Ly, Lz, xy, xz, yz) are host scalars; they broadcast
    cleanly under either backend.
    """
    Lx, Ly, Lz = box.Lx, box.Ly, box.Lz
    xy, xz, yz = box.xy, box.xz, box.yz

    rz = pos[:, 2]
    ry = pos[:, 1]
    rx = pos[:, 0]

    fz = rz / Lz
    fy = (ry - yz * Lz * fz) / Ly
    fx = (rx - xy * Ly * fy - xz * Lz * fz) / Lx

    nx = xp.round(fx)
    ny = xp.round(fy)
    nz = xp.round(fz)

    # §4 numerical-sanity guard (PI 2026-05-20): in a well-behaved
    # overdamped step the per-step displacement is ≪ ℓ₀, so the
    # fractional-coordinate magnitude after a single update is O(1).
    # If a force overflow (e.g. inter-fiber LJ overlap before
    # equilibration) drives |f| past ~1e8, the subsequent np.round →
    # int32 cast for ``image_delta`` produces NaN silently and the
    # simulation continues with corrupt image flags. Catch that here
    # — well inside the int32 range (≈ 2.1e9) — and raise so the
    # caller surfaces to PI rather than progressing on garbage state.
    # Batch the three per-axis max() reductions into one host sync so
    # the GPU hot path pays only ~1 sync per step for this guard.
    INT32_GUARD = 1.0e8
    worst = float(
        xp.maximum(
            xp.maximum(xp.abs(nx).max(), xp.abs(ny).max()), xp.abs(nz).max()
        )
    )
    if worst > INT32_GUARD:
        raise FloatingPointError(
            "BAOAB _wrap_into_box: |fractional coord| exceeded the "
            f"int32-image guard (worst |round(f)|={worst:.3e} > "
            f"{INT32_GUARD:.0e}). A per-step displacement many box-lengths "
            "long indicates an unphysical force (e.g. LJ overlap before "
            "equilibration); add a no-shear equilibration prelude or "
            "reduce dt. See ffn_sim/ecm/equilibrate.py."
        )

    fx -= nx
    fy -= ny
    fz -= nz

    out = xp.empty_like(pos)
    out[:, 0] = Lx * fx + xy * Ly * fy + xz * Lz * fz
    out[:, 1] = Ly * fy + yz * Lz * fz
    out[:, 2] = Lz * fz

    img_delta = xp.empty_like(pos, dtype=xp.int32)
    img_delta[:, 0] = nx.astype(xp.int32)
    img_delta[:, 1] = ny.astype(xp.int32)
    img_delta[:, 2] = nz.astype(xp.int32)
    return out, img_delta


def make_baoab_updater(
    *,
    kT: float,
    gamma: Mapping[str, float],
    dt: float,
    seed: int = 0,
) -> tuple[LeimkuhlerMatthewsBAOAB, hoomd.update.CustomUpdater]:
    """Build the L-M Action wrapped in a per-step ``CustomUpdater``.

    Returns the (action, updater) pair so callers can inspect Action
    state (``prv_rnds``, ``steps_run``) after a run while keeping the
    Updater attached to ``sim.operations.updaters``.
    """
    action = LeimkuhlerMatthewsBAOAB(kT=kT, gamma=gamma, dt=dt, seed=seed)
    updater = hoomd.update.CustomUpdater(
        action=action, trigger=hoomd.trigger.Periodic(1)
    )
    return action, updater
