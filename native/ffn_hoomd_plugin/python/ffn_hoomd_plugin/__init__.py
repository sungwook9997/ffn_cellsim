"""Python wrappers for the experimental FFN native HOOMD plugin."""

from __future__ import annotations

import hoomd
from hoomd.operation import Updater
from hoomd.md.force import Force

# The pure-Python host-side fixed-pool (Stage 2 binder) has NO native dependency
# and must import on CPU-only dev boxes for parity testing. Re-export it first,
# before the GPU-only compiled module, so it is available even where
# ``_ffn_native`` is not built.
from .myosin_pool import MyosinAttachmentPool  # noqa: F401

import _ffn_native


class NativeRadialShellForce(Force):
    """Native device-resident radial-shell compartment force (Stage 2 native).

    One ForceCompute for all 3 production compartment forces (the cupy
    ``*_gpu.py`` twins removed the cpu_local_snapshot host-sync but stayed
    launch-bound at a ~233 us/force gpu_local_snapshot floor; this avoids the
    floor entirely with ArrayHandle). law: 0=nucleus, 1=membrane, 2=turgor.
    Use the ``from_*`` constructors. CPU device is unsupported (GPU-only).
    """

    def __init__(self, *, law, tag_start, tag_end, R0=0.0, pa, pb, pc, pd=0.0):
        super().__init__()
        self._law = int(law)
        self._t0 = int(tag_start)
        self._t1 = int(tag_end)
        self._R0 = float(R0)
        self._pa = float(pa)
        self._pb = float(pb)
        self._pc = float(pc)
        self._pd = float(pd)

    @classmethod
    def from_nucleus(cls, p, tag_range):
        return cls(law=0, tag_start=tag_range[0], tag_end=tag_range[1], R0=p.R_nuc,
                   pa=p.k_chrom, pb=p.k_lamin, pc=p.d_knee, pd=p.F_knee)

    @classmethod
    def from_membrane(cls, p, tag_range):
        return cls(law=1, tag_start=tag_range[0], tag_end=tag_range[1],
                   pa=p.gamma_mem, pb=p.K_A, pc=p.A0)

    @classmethod
    def from_turgor(cls, p, tag_range):
        return cls(law=2, tag_start=tag_range[0], tag_end=tag_range[1],
                   pa=p.turgor_dP0, pb=p.K_vol, pc=p.V0)

    def _attach_hook(self) -> None:
        self._cpp_obj = _ffn_native.FFNRadialShellForce(
            self._simulation.state._cpp_sys_def, self._law, self._t0, self._t1,
            self._R0, self._pa, self._pb, self._pc, self._pd)


class NativeAttachmentSpringForce(Force):
    """Fixed-pool binder enabler: K head<->actin springs read from device arrays.

    The binder updater calls ``set_attachments`` on its firing (~1% of steps) to
    toggle the pool (head_tag/actin_tag = bound pair, k=0 = inactive); the force
    sums the springs every step on-device. Replaces the per-firing global
    set_snapshot bond mutation (51.9 ms — h7_native_gate_a_profile). GPU-only.
    """

    def __init__(self, *, pool_size):
        super().__init__()
        self._K = int(pool_size)
        self._pending = None

    def set_attachments(self, head_tag, actin_tag, k, r0):
        import numpy as np
        a = (np.ascontiguousarray(head_tag, np.int32),
             np.ascontiguousarray(actin_tag, np.int32),
             np.ascontiguousarray(k, np.float64),
             np.ascontiguousarray(r0, np.float64))
        if getattr(self, "_cpp_obj", None) is None:
            self._pending = a
        else:
            self._cpp_obj.set_attachments(*a)

    def _attach_hook(self) -> None:
        self._cpp_obj = _ffn_native.FFNAttachmentSpringForce(
            self._simulation.state._cpp_sys_def, self._K)
        if self._pending is not None:
            self._cpp_obj.set_attachments(*self._pending)
            self._pending = None


class NativeNoOpUpdater(Updater):
    """Minimal compiled updater used to verify HOOMD native-operation wiring."""

    def __init__(self, trigger: hoomd.trigger.trigger_like = 1) -> None:
        super().__init__(trigger)

    def _attach_hook(self) -> None:
        self._cpp_obj = _ffn_native.FFNNoOpUpdater(
            self._simulation.state._cpp_sys_def, self.trigger
        )

    @property
    def update_count(self) -> int:
        if self._cpp_obj is None:
            return 0
        return int(self._cpp_obj.get_update_count())

    @property
    def last_timestep(self) -> int:
        if self._cpp_obj is None:
            return 0
        return int(self._cpp_obj.get_last_timestep())


class NativeBaoabUpdater(Updater):
    """Native overdamped Leimkuhler-Matthews BAOAB-limit integrator (Stage 1a).

    Attach to an ``md.Integrator(forces=[...], methods=[])`` (HOOMD computes the
    net force each step; this updater reads it and advances positions in native
    CUDA). Mirrors ``ffn_sim.integrator.baoab_device.OverdampedBAOABDevice`` and
    shares its splitmix64 counter-RNG keyed by ``(seed, timestep, tag)`` — so for
    a fixed dense ``[0, N)`` tag space the two are bit-identical given the same
    force sequence. ``gamma_by_tag`` is the per-tag Stokes drag (length N).
    """

    def __init__(
        self,
        *,
        dt: float,
        kT: float,
        gamma_by_tag,
        seed: int = 0,
        trigger: hoomd.trigger.trigger_like = 1,
    ) -> None:
        super().__init__(trigger)
        self._dt = float(dt)
        self._kT = float(kT)
        self._seed = int(seed)
        self._gamma_by_tag = [float(g) for g in gamma_by_tag]

    def _attach_hook(self) -> None:
        self._cpp_obj = _ffn_native.FFNBaoabUpdater(
            self._simulation.state._cpp_sys_def,
            self.trigger,
            self._dt,
            self._seed,
            self._kT,
            self._gamma_by_tag,
        )


class NativeConstrainedBaoabUpdater(Updater):
    """Native device-resident constrained L-M BAOAB (Fixman + M-SHAKE), Stage 1b.

    Mirrors ``constrained_baoab.ConstrainedLeimkuhlerMatthewsBAOAB`` on a GPU
    device. ``inv_gamma_by_tag`` = per-tag mobility (1/γ, length N); ``chains_tag``
    = flat ``F*(m+1)`` bead TAGS (uniform-length linear chains). Attach to an
    ``md.Integrator(forces=[...], methods=[])`` with the stretch bond OMITTED
    (it is the constraint).
    """

    def __init__(
        self,
        *,
        dt: float,
        kT: float,
        inv_gamma_by_tag,
        chains_tag,
        n_chains: int,
        bonds_per_chain: int,
        rest_length: float,
        seed: int = 0,
        tol: float = 1.0e-10,
        max_iter: int = 100,
        trigger: hoomd.trigger.trigger_like = 1,
    ) -> None:
        super().__init__(trigger)
        self._dt = float(dt)
        self._kT = float(kT)
        self._seed = int(seed)
        self._inv_gamma_by_tag = [float(g) for g in inv_gamma_by_tag]
        self._chains_tag = [int(c) for c in chains_tag]
        self._F = int(n_chains)
        self._m = int(bonds_per_chain)
        self._rest_length = float(rest_length)
        self._tol = float(tol)
        self._max_iter = int(max_iter)

    def _attach_hook(self) -> None:
        self._cpp_obj = _ffn_native.FFNConstrainedBaoabUpdater(
            self._simulation.state._cpp_sys_def,
            self.trigger,
            self._dt,
            self._seed,
            self._kT,
            self._inv_gamma_by_tag,
            self._chains_tag,
            self._F,
            self._m,
            self._rest_length,
            self._tol,
            self._max_iter,
        )

    @property
    def lambda_buf(self):
        """Last step's accumulated per-bond Lagrange multipliers, shape (F, m).

        Consumers convert to scalar bond tension ``T_bond = λ·r0/dt`` for the
        method-of-planes γ (RIGID_LAGRANGE_TENSION_DESIGN.md). ``None`` before
        attach / first step.
        """
        if self._cpp_obj is None:
            return None
        return self._cpp_obj.get_lambda()

    @property
    def nonconverged_count(self) -> int:
        if self._cpp_obj is None:
            return 0
        return int(self._cpp_obj.nonconverged_count())


class NativePositionKickUpdater(Updater):
    """Move all particles by a fixed vector in native code every trigger tick."""

    def __init__(
        self,
        delta: tuple[float, float, float],
        trigger: hoomd.trigger.trigger_like = 1,
    ) -> None:
        super().__init__(trigger)
        self._delta = tuple(float(x) for x in delta)
        if len(self._delta) != 3:
            raise ValueError("delta must contain exactly three values")

    def _attach_hook(self) -> None:
        self._cpp_obj = _ffn_native.FFNPositionKickUpdater(
            self._simulation.state._cpp_sys_def, self.trigger, *self._delta
        )

    @property
    def update_count(self) -> int:
        if self._cpp_obj is None:
            return 0
        return int(self._cpp_obj.get_update_count())

    @property
    def last_timestep(self) -> int:
        if self._cpp_obj is None:
            return 0
        return int(self._cpp_obj.get_last_timestep())
