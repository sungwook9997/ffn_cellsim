"""Device-aware simple overdamped Leimkuhler-Matthews BAOAB Action (GPU-main port, B2).

Additive sibling of the FROZEN ``integrator.baoab.LeimkuhlerMatthewsBAOAB`` for the Layer-2
center-based spheroid (CBM) native-N GPU runs. The frozen single-cell BAOAB is shared by
H.1–H.5 and stays byte-for-byte untouched (CLAUDE.md integrator-freeze); this module is a
NEW file the Layer-2 CBM opts into when it runs on a GPU device.

Why a separate Action (not editing baoab.py)
---------------------------------------------
1. **Freeze safety.** The single-cell line's gates depend on ``baoab.py`` bit-for-bit. A new
   file carries zero regression risk to H.1–H.5.
2. **The CBM case is strictly simpler.** The L2.4b growth pool pre-allocates ``max_cells``
   particles at t=0; division only flips a parked ``void`` to a ``cell`` (typeid + position
   via ``set_snapshot``), so the **tag space is fixed for the whole run**. The frozen Action's
   Path-A tag-extension machinery (``_extend_tag_buffers``) is therefore unnecessary here —
   this Action asserts a fixed dense tag space and stays minimal.

What the port buys (the B2 bottleneck)
--------------------------------------
For the CBM the cohesion (``md.pair.Morse`` / ``md.pair.Table``) and the neighbour list run
in HOOMD C++ — already GPU-accelerated for free. The ONLY per-step Python is this Action, and
the frozen one reads ``cpu_local_snapshot`` every step, forcing a device→host sync that
throttles the GPU. On a GPU device this Action uses ``gpu_local_snapshot`` + cupy so
positions / forces / RNG stay device-resident (the sync is removed), exactly as the
constrained Action was ported (``GPU_MAIN_PORT_PHASE1_2026-06-01.md``).

Correctness contract
--------------------
- **CPU path (``xp=np``) is bit-identical to the frozen ``LeimkuhlerMatthewsBAOAB``**: same
  per-Action ``np.random.default_rng(seed)`` stream, same draw order
  (``W_row = rng.standard_normal((N, 3))`` once per step), same tag-indexed gather/scatter,
  same frozen ``_wrap_into_box``. ``tests/test_baoab_device.py`` asserts trajectory parity.
- **GPU path (``xp=cp``)** uses cupy's own RNG stream, so GPU↔CPU agreement is *statistical*
  (within seed noise), not bit-exact — the same policy the constrained port adopted.

Update rule (per particle, unchanged from the frozen Action)::

    r(t+Δt) = r(t) + (F/γ)·Δt + √(kT/(2 γ Δt))·(W_n + W_{n-1})·Δt

Sanity Gate
-----------
- Dimensional / boundary / sign-sense: identical to ``integrator.baoab`` (the formula is the
  same); see that module's Sanity Gate. The only change is the array backend / snapshot
  context, neither of which touches the physics.
- Numerical: NaN/Inf force + position guards retained; the frozen ``_wrap_into_box`` int32
  image-overflow guard is used on CPU, its device sibling ``_wrap_into_box_xp`` on GPU.
- Conservation / measurement: same as the frozen Action — this Action only advances positions.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

import hoomd
import hoomd.custom

from ffn_sim.integrator.baoab import _wrap_into_box
from ffn_sim.integrator.constrained_baoab import _wrap_into_box_xp, array_backend

# --------------------------------------------------------------------------- #
# GPU VRAM budget (derived memory estimate — NOT a physics constant).          #
# --------------------------------------------------------------------------- #
# Per-particle device-resident byte budget at native-N. This counts the arrays
# that scale with N and live on the GPU during a step (the per-step transients
# the snapshot exposes + this Action's persistent per-tag buffers); it is a
# headroom estimate, not an exact allocator trace.
#   HOOMD state (read/written through gpu_local_snapshot, float64 default build):
#     position   3 × f64 = 24 B   net_force  3 × f64 = 24 B
#     image      3 × i32 = 12 B   tag/typeid 2 × u32 =  8 B
#   This Action's persistent per-tag buffers:
#     gamma_by_tag        1 × f64 =  8 B   bd_prefactor_by_tag 1 × f64 =  8 B
#     prv_rnds            3 × f64 = 24 B
#   Per-step transients the kernels materialise (W_row, dr, new_pos, wrapped,
#   img_delta, the (N,1) gather rows) — bounded above by ~8 × (3 × f64) = 192 B.
# Sum ≈ 24+24+12+8+8+8+24+192 = 300 B; round UP to 512 B/particle so the guard
# is conservative (covers cupy pool fragmentation + neighbour-list/Morse C++
# arrays HOOMD also holds, which scale with N).
_BYTES_PER_PARTICLE_GPU: int = 512
# A5000 has 16 GiB; reserve ~2 GiB for the CUDA context, cupy pool slack, and
# the HOOMD C++ pair/nlist working set we do not itemise above. Documented
# numerical-policy margin, not a tuned constant.
_GPU_VRAM_TOTAL_BYTES: int = 16 * 1024**3
_GPU_VRAM_RESERVE_BYTES: int = 2 * 1024**3
_GPU_VRAM_BUDGET_BYTES: int = _GPU_VRAM_TOTAL_BYTES - _GPU_VRAM_RESERVE_BYTES


class OverdampedBAOABDevice(hoomd.custom.Action):
    """Device-aware simple overdamped L-M BAOAB step for a FIXED tag space (CBM).

    Mirrors ``integrator.baoab.LeimkuhlerMatthewsBAOAB`` on CPU (bit-identical) and runs
    device-resident (``gpu_local_snapshot`` + cupy) on a GPU device. Requires a fixed dense
    tag space ``[0, N)`` — the Layer-2 growth pool satisfies this (all particles exist at t=0;
    division flips ``void``→``cell`` in place). Tag-space growth/shrinkage is rejected.

    Parameters
    ----------
    kT, gamma, dt, seed
        As in the frozen Action: ``gamma`` is a per-particle-type Stokes drag map covering
        every type in the state (the CBM passes ``{"cell": γ, "void": γ·1e6}``).
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
            raise ValueError(f"kT must be finite and >= 0, got {kT!r}")
        if not (np.isfinite(dt) and dt > 0.0):
            raise ValueError(f"dt must be finite and > 0, got {dt!r}")
        for typ, g in gamma.items():
            if not (np.isfinite(g) and g > 0.0):
                raise ValueError(f"gamma['{typ}'] must be finite and > 0, got {g!r}")
        self.kT = float(kT)
        self.dt = float(dt)
        self.gamma_map = dict(gamma)
        self._seed = int(seed)

        self._gamma_by_tag: np.ndarray | None = None
        self._bd_prefactor_by_tag: np.ndarray | None = None
        self._prv_rnds: np.ndarray | None = None  # (N, 3)
        self._sim_ref: hoomd.Simulation | None = None
        self._on_gpu: bool = False
        self._xp = np
        self._rng = np.random.default_rng(seed)
        self._steps_run: int = 0

    # ------------------------------------------------------------------
    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        """Allocate tag-indexed buffers; pick the device backend; sanity-check wiring."""
        super().attach(simulation)
        self._sim_ref = simulation

        ig = simulation.operations.integrator
        if ig is None:
            raise RuntimeError(
                "OverdampedBAOABDevice requires an md.Integrator with forces attached "
                "(methods=[]). None is set."
            )
        if not np.isclose(float(ig.dt), self.dt, rtol=0, atol=0):
            raise RuntimeError(
                f"Integrator dt={float(ig.dt)!r} must equal Action dt={self.dt!r}."
            )
        if len(ig.methods) != 0:
            raise RuntimeError(
                "md.Integrator.methods must be empty; a position-updating Method would "
                f"double-step alongside the L-M Action. Got: {list(ig.methods)!r}"
            )

        # Tag-indexed gamma is built from a host snapshot once (cheap, attach-time only).
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
            if N == 0:
                raise RuntimeError("L-M Action attached to empty state.")
            if int(tags.max()) >= N or int(tags.min()) < 0:
                raise RuntimeError(
                    "OverdampedBAOABDevice requires dense particle tags in [0, N); "
                    f"got tag range [{int(tags.min())}, {int(tags.max())}] for N={N}."
                )
            gamma_by_typeid = np.array(
                [self.gamma_map[t] for t in type_names], dtype=np.float64
            )
            gamma_by_tag = np.empty(N, dtype=np.float64)
            gamma_by_tag[tags] = gamma_by_typeid[type_ids]

        self._gamma_by_tag = gamma_by_tag
        self._bd_prefactor_by_tag = np.sqrt(
            self.kT / (2.0 * gamma_by_tag * self.dt)
        ).reshape(-1, 1)
        self._prv_rnds = np.zeros((N, 3), dtype=np.float64)
        self._steps_run = 0

        # Follow the simulation's device (constrained-port policy): a GPU device flips this
        # Action onto gpu_local_snapshot + cupy so the per-step host sync is removed; a CPU
        # device keeps the bit-identical numpy path.
        self._on_gpu = isinstance(simulation.device, hoomd.device.GPU)
        if self._on_gpu:
            # --- GPU VRAM pre-check (graceful fail before any device alloc) ---
            # Estimate the device-resident footprint at this N against the A5000
            # budget and raise BEFORE allocating, so a too-large N is a clear
            # diagnostic instead of a silent CUDA OOM mid-allocation. The byte
            # budget is a derived memory estimate (see _BYTES_PER_PARTICLE_GPU);
            # it is not a physics constant.
            est_bytes = N * _BYTES_PER_PARTICLE_GPU
            if est_bytes > _GPU_VRAM_BUDGET_BYTES:
                raise ValueError(
                    "OverdampedBAOABDevice GPU VRAM pre-check failed: estimated "
                    f"device footprint {est_bytes / 1024**2:.0f} MB for N={N} "
                    f"particles (@ {_BYTES_PER_PARTICLE_GPU} B/particle) exceeds "
                    f"the {_GPU_VRAM_BUDGET_BYTES / 1024**2:.0f} MB budget "
                    f"(16 GB device − {_GPU_VRAM_RESERVE_BYTES / 1024**2:.0f} MB "
                    "reserve). Reduce N, raise the ×40 mesoscale coarse-graining, "
                    "or run on a larger-VRAM device."
                )

            xp = array_backend(True)

            # --- cupy/HOOMD device-match check (multi-GPU safety) ---
            # On a multi-GPU box cupy's current device must be the SAME CUDA
            # device HOOMD runs on, else every gpu_local_snapshot↔cupy op is a
            # silent cross-device copy (wrong results / huge slowdown). The
            # introspection is best-effort: a single-GPU box (the common case)
            # has device 0 == 0 and never trips, and any introspection failure
            # is treated as "cannot disprove a match" so it never breaks the
            # common path.
            try:
                hoomd_gpu_ids = list(simulation.device.gpu_ids)
            except Exception:  # noqa: BLE001 — HOOMD API/version variance is non-fatal here
                hoomd_gpu_ids = []
            try:
                cupy_dev_id = int(xp.cuda.runtime.getDevice())
            except Exception:  # noqa: BLE001 — cupy/driver variance is non-fatal here
                cupy_dev_id = None
            if hoomd_gpu_ids and cupy_dev_id is not None:
                if cupy_dev_id not in hoomd_gpu_ids:
                    raise RuntimeError(
                        "OverdampedBAOABDevice device mismatch: cupy current "
                        f"CUDA device {cupy_dev_id} is not in HOOMD's device "
                        f"set {hoomd_gpu_ids}. Pin cupy to HOOMD's GPU "
                        "(e.g. cupy.cuda.Device(<id>).use()) before attach so "
                        "device-resident positions/forces are not silently "
                        "copied across GPUs."
                    )

            self._xp = xp
            self._rng = xp.random.default_rng(self._seed)  # cupy stream (statistical parity)
            self._gamma_by_tag = xp.asarray(self._gamma_by_tag)
            self._bd_prefactor_by_tag = xp.asarray(self._bd_prefactor_by_tag)
            self._prv_rnds = xp.asarray(self._prv_rnds)

    # ------------------------------------------------------------------
    def act(self, timestep: int) -> None:
        """Apply one overdamped L-M step (device-resident on GPU, numpy on CPU)."""
        if self._prv_rnds is None or self._gamma_by_tag is None:
            raise RuntimeError(
                "OverdampedBAOABDevice.act called before attach(); wrap in "
                "hoomd.update.CustomUpdater and add to sim.operations.updaters."
            )
        sim = self._sim_ref
        assert sim is not None
        xp = self._xp

        snap_ctx = (sim.state.gpu_local_snapshot if self._on_gpu
                    else sim.state.cpu_local_snapshot)
        with snap_ctx as snap:
            pos = xp.asarray(snap.particles.position)      # (N, 3) rw
            F = xp.asarray(snap.particles.net_force)       # (N, 3) ro
            image = xp.asarray(snap.particles.image)       # (N, 3) rw
            tag = xp.asarray(snap.particles.tag)           # (N,) stable id
            N = pos.shape[0]

            if N != self._prv_rnds.shape[0]:
                raise RuntimeError(
                    f"Particle count changed ({self._prv_rnds.shape[0]}->{N}); "
                    "OverdampedBAOABDevice requires a fixed tag space (the CBM pool is "
                    "pre-allocated — division flips void->cell in place, never adds tags)."
                )

            if not bool(xp.all(xp.isfinite(F))):
                bad = xp.argwhere(~xp.isfinite(F))[:5]
                if self._on_gpu:
                    bad = xp.asnumpy(bad)
                raise FloatingPointError(
                    f"Non-finite net_force at timestep={timestep}; first: {bad.tolist()}."
                )

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
                bad = xp.argwhere(~xp.isfinite(new_pos))[:5]
                if self._on_gpu:
                    bad = xp.asnumpy(bad)
                raise FloatingPointError(
                    f"Non-finite position after L-M step at timestep={timestep}; "
                    f"first: {bad.tolist()}."
                )

            box = sim.state.box
            if self._on_gpu:
                wrapped, img_delta = _wrap_into_box_xp(new_pos, box, xp)
            else:
                wrapped, img_delta = _wrap_into_box(new_pos, box)
            pos[:] = wrapped
            image[:] = image + img_delta

            self._prv_rnds[tag] = W_row

        self._steps_run += 1

    # ------------------------------------------------------------------
    @property
    def prv_rnds(self) -> np.ndarray | None:
        """Host-numpy view of the persisted W_{n-1} buffer (None pre-attach)."""
        if self._prv_rnds is None:
            return None
        buf = self._prv_rnds
        if self._on_gpu:
            buf = self._xp.asnumpy(buf)
        v = np.asarray(buf).view()
        v.flags.writeable = False
        return v

    @property
    def steps_run(self) -> int:
        return self._steps_run


def make_baoab_updater_device(
    *,
    kT: float,
    gamma: Mapping[str, float],
    dt: float,
    seed: int = 0,
) -> tuple[OverdampedBAOABDevice, hoomd.update.CustomUpdater]:
    """Build the device-aware overdamped L-M Action wrapped in a per-step ``CustomUpdater``.

    Drop-in for ``integrator.baoab.make_baoab_updater`` for the Layer-2 CBM: identical CPU
    behaviour, plus a GPU-resident path when the sim is built on a GPU device.
    """
    action = OverdampedBAOABDevice(kT=kT, gamma=gamma, dt=dt, seed=seed)
    updater = hoomd.update.CustomUpdater(
        action=action, trigger=hoomd.trigger.Periodic(1)
    )
    return action, updater


def make_baoab_updater_for_device(
    device: "hoomd.device.Device | None",
    *,
    kT: float,
    gamma: Mapping[str, float],
    dt: float,
    seed: int = 0,
):
    """Pick the BAOAB updater by device: frozen numpy Action on CPU, device Action on GPU.

    On a CPU device (or ``None`` → CPU default) this returns the FROZEN
    ``integrator.baoab.make_baoab_updater`` so every existing CPU result stays byte-for-byte
    unchanged (the device Action is bit-identical on CPU, but selecting the frozen path keeps
    the provenance unambiguous). On a GPU device it returns the device-resident Action so the
    per-step host sync is removed — the B2 native-N unlock. Returns ``(action, updater)``.
    """
    import sys
    from ffn_sim.integrator.baoab import make_baoab_updater  # frozen CPU Action
    # Provenance: one line so a production log unambiguously records the path. MUST go to
    # stderr — a data driver (layer2_gpu_scaleup) emits its JSON record on stdout, and a
    # stdout print here pollutes the redirected JSONL (broke the production auto-fit, 2026-06-04).
    if isinstance(device, hoomd.device.GPU):
        print("[baoab_device] BAOAB path: device-aware GPU cupy path", file=sys.stderr)
        return make_baoab_updater_device(kT=kT, gamma=gamma, dt=dt, seed=seed)
    print("[baoab_device] BAOAB path: frozen CPU numpy path", file=sys.stderr)
    return make_baoab_updater(kT=kT, gamma=gamma, dt=dt, seed=seed)
