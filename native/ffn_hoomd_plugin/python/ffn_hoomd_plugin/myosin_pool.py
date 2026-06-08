"""Host-side fixed-pool data structure for the native myosin grip-walk binder.

Stage 2 of ``NATIVE_HOT_LOOP_MIGRATION_2026-06-07.md``: move the per-step myosin
attach-bond force off the per-step ``get_snapshot()`` / ``set_snapshot()`` host
round-trip that ``cortex.myosin.MyosinStepUpdater`` currently pays.

What the per-step myosin force actually is
------------------------------------------
The production cortical myosin runtime (``ffn_sim/cortex/myosin.py``) splits into
two timescales:

* **Every MD step** HOOMD evaluates the head-actin attach bonds as ordinary
  ``md.bond.Harmonic`` springs.  In the ratified ``grip_walk`` stepping mode the
  attach bond TYPES all carry ``r0 = _GRIP_WALK_R0_EPS = 1e-12 m`` (a force-
  negligible 1 pm; see ``cortex_myosin_attach_bin_rest_lengths``), so the
  delivered force on every engaged head is the harmonic spring

      F_head  = +k_head_actin * (r - r0_eps) * (r_actin - r_head)/r
      F_actin = -F_head                                   (Newton's 3rd law)

  between the head particle tag and the cortical-actin bead it is currently
  gripping (``k = p_myo.k_head_actin``).  This is the entire per-step myosin
  contribution to the net force.

* **Only on the batch tick** (``batch_dt``, ~1% of steps) does the updater run
  Bell-Evans unbinding, KDTree binding, and the Hill grip-walk advance.  Those
  re-point each head onto a (possibly different) actin bead — i.e. they mutate
  *which* ``(head_tag, actin_tag)`` pair carries the spring, not the spring law.
  Today that re-pointing is done by rewriting ``snap.bonds.group`` and calling
  ``set_snapshot`` (the 51.9 ms/firing wall in ``h7_native_gate_a_profile``).

Native fixed-pool model
-----------------------
A pool of ``K = 2 * n_heads_per_side * n_motors_per_cell`` slots, **one slot per
myosin head**, each holding ``(head_tag, actin_tag, k, r0)``:

* engaged head  -> ``actin_tag = bound bead``, ``k = k_head_actin``, ``r0 = eps``
* free head     -> ``actin_tag = -1``, ``k = 0`` (slot contributes no force)

HOOMD's attach-bond *count never changes* during production: the head-actin
attach bonds are removed from the HOOMD topology entirely and replaced by this
device-resident pool, which the native ``FFNAttachmentSpringForce`` ForceCompute
sums every step on-GPU.  Binding state changes become a ``set_attachments`` array
write on the batch tick (the off-hot-path ~1% firing), never a HOOMD bond
mutation or a global snapshot.

This module provides:

* :class:`MyosinAttachmentPool` -- the host-side pool, an
  ``update_from_updater`` adaptor that reads the production
  ``MyosinStepUpdater``'s per-head ``_head_bound_to_actin`` state into the pool
  arrays, and a NumPy ``reference_forces`` that reproduces the exact per-bead
  force the production attach bonds deliver.  ``reference_forces`` is the CPU
  ground-truth the ``.cu`` kernel is parity-tested against and the fallback when
  no GPU is present.

No magic numbers: ``k`` and ``r0`` come straight from the resolved
``ResolvedCortexMyosin`` (``k_head_actin`` and the grip-walk attach-bond rest
length), which are themselves config/KU-derived.  This module introduces no new
physical constant.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _min_image(d: np.ndarray, L: np.ndarray) -> np.ndarray:
    """Minimum-image displacement, matching the ``.cu`` ``as_minimage`` (rint)."""
    return d - L * np.rint(d / L)


@dataclass
class MyosinAttachmentPool:
    """Fixed-size host-side pool of myosin head <-> actin harmonic springs.

    Args:
        pool_size: Number of slots ``K`` (one per myosin head). Constant for the
            life of the run -- HOOMD attach-bond count never changes.
        k_head_actin: Head-actin spring constant [N/m] (``p_myo.k_head_actin``).
        r0: Attach-bond rest length [m] (``_GRIP_WALK_R0_EPS`` in grip_walk mode).

    The slot arrays are public and contiguous so they can be handed straight to
    ``NativeAttachmentSpringForce.set_attachments`` (or the ``.cu`` kernel).
    """

    pool_size: int
    k_head_actin: float
    r0: float

    def __post_init__(self) -> None:
        K = int(self.pool_size)
        self.head_tag = np.full(K, -1, dtype=np.int32)
        self.actin_tag = np.full(K, -1, dtype=np.int32)
        self.k = np.zeros(K, dtype=np.float64)
        self.r0_arr = np.full(K, float(self.r0), dtype=np.float64)
        self._k_head_actin = float(self.k_head_actin)

    # -- construction-time wiring -------------------------------------------
    def set_head_tags(self, head_tags: np.ndarray) -> None:
        """Assign the (fixed) global particle tag of each head slot.

        Slot ``i`` corresponds to head local index ``i`` in the production
        ``MyosinStepUpdater`` numbering (``_head_global_tag(i)``). Head tags are
        fixed for the run; only ``actin_tag``/``k`` toggle per binding event.
        """
        head_tags = np.ascontiguousarray(head_tags, dtype=np.int32)
        if head_tags.shape[0] != self.head_tag.shape[0]:
            raise ValueError(
                f"head_tags length {head_tags.shape[0]} != pool_size "
                f"{self.head_tag.shape[0]}"
            )
        self.head_tag[:] = head_tags

    # -- per-firing toggle ---------------------------------------------------
    def update_from_updater(self, updater) -> None:
        """Read a production ``MyosinStepUpdater``'s per-head bound state.

        Mirrors exactly which ``(head_tag, actin_tag)`` attach bonds the updater
        would have written into ``snap.bonds.group``:

        * ``updater._head_bound_to_actin[i] >= 0`` -> slot ``i`` is engaged with
          ``actin_tag = that bead``, ``k = k_head_actin``.
        * ``< 0`` (free) -> ``actin_tag = -1``, ``k = 0`` (no force).

        Head tags are taken from ``updater._head_global_tag`` so the pool stays
        in lock-step with the updater's numbering even if ``set_head_tags`` was
        not called.
        """
        bound = np.asarray(updater._head_bound_to_actin, dtype=np.int64)
        K = self.head_tag.shape[0]
        if bound.shape[0] != K:
            raise ValueError(
                f"updater has {bound.shape[0]} heads but pool_size is {K}"
            )
        # Head tags (fixed) -- fill if not already wired.
        if np.any(self.head_tag < 0):
            self.head_tag[:] = np.fromiter(
                (updater._head_global_tag(i) for i in range(K)),
                dtype=np.int32,
                count=K,
            )
        engaged = bound >= 0
        self.actin_tag[:] = np.where(engaged, bound, -1).astype(np.int32)
        self.k[:] = np.where(engaged, self._k_head_actin, 0.0)

    def push_to(self, native_force) -> None:
        """Send the current pool arrays to a ``NativeAttachmentSpringForce``.

        Call on the binder firing (off-hot-path) after ``update_from_updater``:
        a single ``set_attachments`` device-array write, NO HOOMD bond mutation
        and NO global ``set_snapshot``. ``native_force`` is a
        :class:`ffn_hoomd_plugin.NativeAttachmentSpringForce` of the same K.
        """
        native_force.set_attachments(
            self.head_tag, self.actin_tag, self.k, self.r0_arr
        )

    # -- exact CPU reference force eval (ground truth + GPU fallback) --------
    def reference_forces(
        self, positions: np.ndarray, box_L, tag_to_row: np.ndarray | None = None
    ) -> np.ndarray:
        """Per-particle force from the pool, identical to the ``.cu`` kernel.

        Args:
            positions: ``(N, 3)`` positions indexed by *row* (HOOMD storage
                order). If ``tag_to_row`` is None, row == tag (dense ``[0, N)``).
            box_L: box edge lengths ``(Lx, Ly, Lz)`` for the minimum image.
            tag_to_row: optional ``(max_tag+1,)`` map; ``tag_to_row[t]`` is the
                row of tag ``t`` (matches the kernel's ``row_of_tag`` scatter).

        Returns:
            ``(N, 3)`` force array. Inactive slots (``k <= 0`` or tag ``< 0``)
            contribute nothing; a bead carrying many heads accumulates (atomic
            add in the kernel; ``np.add.at`` here).
        """
        positions = np.ascontiguousarray(positions, dtype=np.float64)
        N = positions.shape[0]
        L = np.asarray(box_L, dtype=np.float64)
        if tag_to_row is None:
            tag_to_row = np.arange(N, dtype=np.int64)
        force = np.zeros((N, 3), dtype=np.float64)

        active = (self.k > 0.0) & (self.head_tag >= 0) & (self.actin_tag >= 0)
        if not np.any(active):
            return force
        h_rows = tag_to_row[self.head_tag[active]]
        a_rows = tag_to_row[self.actin_tag[active]]
        kk = self.k[active]
        r0 = self.r0_arr[active]

        d = _min_image(positions[a_rows] - positions[h_rows], L)  # actin - head
        r = np.linalg.norm(d, axis=1)
        good = r > 0.0
        if not np.any(good):
            return force
        fmag = np.zeros_like(r)
        fmag[good] = kk[good] * (r[good] - r0[good]) / r[good]
        fvec = fmag[:, None] * d  # on head, toward actin

        np.add.at(force, h_rows[good], fvec[good])
        np.add.at(force, a_rows[good], -fvec[good])
        return force
