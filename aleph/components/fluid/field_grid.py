"""Conservative field grid (Warp-CUDA device arrays) + a SEPARATE filament neighbour HashGrid — I1a.

Warp kernel SOURCE for the GPU-resident Eulerian field of the fluid-first substrate. Authored on the
dev Mac (Warp imports in CPU mode) but CONSTRUCTED only on CUDA — the runtime is Warp-CUDA only (I0-A);
the constructor rejects a non-CUDA device. The dev-Mac acceptance path is the pure-NumPy oracle
``fv_reference`` (the identical conservative stencil); this module is what the lead runs on the A5000.

Two spatial accelerators, ONE coordinate frame, NOT one data structure (NEW_ENGINE_BUILD_PLAN §1 P1):
  * the Eulerian Biot/RAD PDE fields live on this structured, conservative finite-volume grid;
  * filament partner search + excluded volume use a Warp ``HashGrid`` (a neighbour structure, NOT the
    PDE field grid). They may share origin/cell size but are distinct allocations.

Device arrays (float64 to match ``ff/biot_fluid_warp.py``; cell-centred on a uniform grid):
  * ``p``          pore pressure field p = p_ext + p_bar + p_excess           [Pa]
  * ``p_new``      double-buffer target (ping-pong, as ``BiotField._p/_p2``)
  * ``p_bar``      mean/turgor channel (resting Pi_0 pre-stress; kept distinct so p_excess is not
                   forced zero-mean during an undrained transient)             [Pa]
  * ``s_water``    per-cell interior fluid source (RAD/reaction/coupling producers)                    [1/s]
  * ``s_membrane`` pressure-dependent outer-membrane hydraulic source, rebuilt each Biot subcycle     [1/s]
  * ``s_total``    device-composed ``s_water + s_membrane`` consumed by the Biot update               [1/s]
  * ``div_vs``     solid-skeleton dilatation rate div(v_s) (frozen outer-step source)   [1/s]
  * ``mask``       cell class OUTSIDE=0 / FLUID=1 / NUCLEUS=2 (domain.Domain owns classification)

Sanity Gate (dimensional + numerical, before first native execution):
  * dx, origin, shape define a uniform conservative grid; face area dx^(d-1), cell volume dx^d;
  * CFL for the explicit consolidation step registered as dt_max = safety * S dx^2 / (2 d mobility);
  * mask codes match ``fv_reference`` so the CPU gate and the device kernel are the same discretisation.
"""

from __future__ import annotations

import numpy as np
import warp as wp  # noqa: E402  (Warp-CUDA runtime; HOOMD is never imported — contract test enforced)

from aleph.components.fluid.fv_reference import (  # shared mask codes (CPU/GPU identical)
    FLUID,
    NUCLEUS,
    OUTSIDE,
)

__all__ = ["FieldGrid", "FLUID", "OUTSIDE", "NUCLEUS"]

_FLUID = wp.constant(1)


@wp.kernel
def _total_content_reduce_kernel(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    storage_S: wp.float64,
    cell_volume: wp.float64,
    out: wp.array(dtype=wp.float64),
) -> None:
    """Reduce ``S * sum_FLUID(p) * dV`` into one device scalar (NG-2/NG-6)."""
    i, j, k = wp.tid()
    if mask[i, j, k] == _FLUID:
        wp.atomic_add(out, 0, storage_S * p[i, j, k] * cell_volume)


@wp.kernel
def _compose_water_sources_kernel(
    interior: wp.array3d(dtype=wp.float64),
    membrane: wp.array3d(dtype=wp.float64),
    total: wp.array3d(dtype=wp.float64),
) -> None:
    """Compose independent source owners without deleting either contribution."""
    i, j, k = wp.tid()
    total[i, j, k] = interior[i, j, k] + membrane[i, j, k]


def _require_cuda(device: str | None) -> str:
    """Resolve and REQUIRE a CUDA device (I0-A: Warp-CUDA is the only runtime; no CPU fallback)."""
    wp.init()
    dev = wp.get_device(device)
    if not dev.is_cuda:
        raise RuntimeError(
            "ac.fluid runtime requires a CUDA GPU (I0-A: Warp-CUDA only, no CPU simulation path). "
            f"Resolved device {dev!r} is not CUDA. Author kernels on the dev Mac; run gates on the A5000. "
            "Host-side acceptance uses aleph.components.fluid.fv_reference (NumPy), never this device class."
        )
    return str(dev)


class FieldGrid:
    """A uniform conservative finite-volume grid backing the Eulerian fields, on CUDA.

    Attributes:
        origin: Cell-centre of index (0,0,0), shape (3,) float [um].
        dx: Uniform spacing [um].
        shape: (nx, ny, nz) cells.
        device: Resolved CUDA device string.
    """

    def __init__(
        self,
        shape: tuple[int, int, int],
        dx: float,
        origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
        *,
        hashgrid_dim: tuple[int, int, int] | None = None,
        device: str | None = None,
    ) -> None:
        self.device = _require_cuda(device)
        self.shape = tuple(int(n) for n in shape)
        self.dx = float(dx)
        self.origin = np.asarray(origin, dtype=np.float64)
        with wp.ScopedDevice(self.device):
            self.p = wp.zeros(self.shape, dtype=wp.float64)
            self.p_new = wp.zeros(self.shape, dtype=wp.float64)
            self.p_bar = wp.zeros(self.shape, dtype=wp.float64)
            self.s_water = wp.zeros(self.shape, dtype=wp.float64)
            self.s_membrane = wp.zeros(self.shape, dtype=wp.float64)
            self.s_total = wp.zeros(self.shape, dtype=wp.float64)
            self.div_vs = wp.zeros(self.shape, dtype=wp.float64)
            self.mask = wp.full(self.shape, FLUID, dtype=wp.int32)
            self._total_content_d = wp.zeros(1, dtype=wp.float64)
            # Separate neighbour-search accelerator for filaments/EV — NOT the PDE field grid.
            hg = hashgrid_dim or self.shape
            self.filament_hashgrid = wp.HashGrid(dim_x=hg[0], dim_y=hg[1], dim_z=hg[2])

    @property
    def dim(self) -> int:
        return len(self.shape)

    def cfl_dt(self, mobility: float, storage_S: float, safety: float = 0.9) -> float:
        """Explicit consolidation CFL ``safety * S dx^2 / (2 d mobility)`` (matches ``fv_reference.cfl_dt``)."""
        return safety * storage_S * self.dx * self.dx / (2.0 * self.dim * mobility)

    def set_pressure(self, p_host: np.ndarray) -> None:
        """Upload a host pressure field (used to seed resting p_bar = Pi_0 or an oracle initial state)."""
        with wp.ScopedDevice(self.device):
            self.p = wp.array(np.ascontiguousarray(p_host, dtype=np.float64), dtype=wp.float64)

    def set_mask(self, mask_host: np.ndarray) -> None:
        """Upload the cell classification (OUTSIDE/FLUID/NUCLEUS) computed by ``domain.Domain``."""
        with wp.ScopedDevice(self.device):
            self.mask = wp.array(np.ascontiguousarray(mask_host, dtype=np.int32), dtype=wp.int32)

    def pressure_to_host(self) -> np.ndarray:
        """Download the pressure field (post-processing / gate readback only — never in the hot loop)."""
        return self.p.numpy()

    def total_content_device(self, storage_S: float) -> wp.array:
        """Return the persistent device scalar ``S sum_FLUID(p) dV`` with no physical-loop readback."""
        with wp.ScopedDevice(self.device):
            self._total_content_d.zero_()
            wp.launch(
                _total_content_reduce_kernel,
                dim=self.shape,
                inputs=[self.p, self.mask, wp.float64(storage_S), wp.float64(self.dx**self.dim)],
                outputs=[self._total_content_d],
                device=self.device,
            )
        return self._total_content_d

    def compose_water_sources(self) -> wp.array:
        """Build the Biot source from independent interior and membrane channels on device."""
        with wp.ScopedDevice(self.device):
            wp.launch(
                _compose_water_sources_kernel,
                dim=self.shape,
                inputs=[self.s_water, self.s_membrane],
                outputs=[self.s_total],
                device=self.device,
            )
        return self.s_total

    def total_content(self, storage_S: float) -> float:
        """Post-loop diagnostic readback of :meth:`total_content_device` (never call from the hot loop)."""
        return float(self.total_content_device(storage_S).numpy()[0])
