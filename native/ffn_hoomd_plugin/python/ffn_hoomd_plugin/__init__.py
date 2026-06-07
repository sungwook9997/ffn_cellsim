"""Python wrappers for the experimental FFN native HOOMD plugin."""

from __future__ import annotations

import hoomd
from hoomd.operation import Updater

import _ffn_native


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
