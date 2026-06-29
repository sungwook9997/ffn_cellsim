"""H.4 emergent FA growth tracker (KU-2.17 / KU-2.2).

The H.4 brief explicitly forbids a Hill wrapper at runtime: the FA
"size" is the direct count of bonded integrins. The Hill form
``k_g(F) = k_g^0 · F^n / (F^n + F_th^n)`` (KU-2.17, Walcott & Sun 2010)
appears only in ``ffn_sim/tests/validation/test_ku217_fa_growth.py`` as
the oracle we compare the emergent dynamics against.

This module ships the **runtime tracker**:

  - ``FAGrowthMonitor`` : per-FA aggregator that reads the engaged mask
    from ``IntegrinBondUpdater`` and the per-FA total force from the
    HOOMD bond-force snapshot, and surfaces the emergent
    ``N_engaged`` / ``F_per_FA`` time series at the
    integrator-friendly tick rate.

  - Aggregator helpers used by the KU-2.17 validation test:
    ``per_fa_force_sum`` reduces per-bond forces onto a per-FA total
    matching the brief's "F^{per-FA} = 50 pN with N_engaged ≈ 10"
    measurement protocol.

Note that the monitor does NOT mutate the simulation; it is a read-only
observer. The growth dynamics (new FAs nucleating, mature FAs
disassembling) are *emergent* from the integrin-clutch dynamics —
NO Hill ODE is integrated in the runtime path.

Sanity Gate
-----------
1. Dimensional analysis: F [N], A [m²], N_engaged [count]. The
   area-equivalent ``A_emergent = (N_engaged / N_total) · A_baseline``
   has units [m²] ✓; ``A_baseline`` is the nascent FA area at
   construction.
2. Boundary cases: ``N_engaged = 0`` → A_emergent = 0; the FA is
   functionally absent but its particle slots remain (week-1 — FA
   nucleation/disassembly enter at week-3 H.5 integration).
3. Conservation: not a conservation law — emergent population dynamics.
4. Numerical sanity: float64; pure NumPy on the engaged-mask + bond-
   force arrays read from the snapshot.
5. Sign / sense: positive F_per_FA = retrograde-flow loading; positive
   N_engaged = bonded integrins. Both monotonic w.r.t. the engaged mask.
6. Measurement-protocol consistency: the KU-2.17 oracle gate
   (``gate_unit2_2_fa_growth``) measures ``F^{per-FA}`` and
   ``N_engaged`` at the trigger frame (first frame where the FA's
   engaged count crosses ``N_engaged_target``). The monitor's
   ``trigger_frames`` table is built to match this protocol exactly.

References
----------
- Walcott & Sun 2010, *PNAS* 107(17):7757-62. KU-2.17 Hill form.
- Tan et al. 2020 (FA growth threshold reference).
- Stricker et al. 2013 (FA size distribution reference).
- H.4 brief ``ffn_sim/docs/briefs/H4_fa_motor_clutch.md`` §FA growth.
- ``ffn_sim/validation/oracles/common/sanity_gate.py::
  gate_unit2_2_fa_growth``.
- v1 archived oracle: ``ffn_sim/validation/oracles/bridge`` (no
  direct fa_growth oracle there; v1 acs_kb/bridge/fa_growth.py is
  not in the v2 oracles tree, so the KU-2.17 closed form is brought
  in inline by the validation test).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ffn_sim.archive.hoomd_legacy.bridge.fa import FALayout, ResolvedH4
    from ffn_sim.archive.hoomd_legacy.bridge.integrin_bonds import IntegrinBondUpdater


@dataclass(slots=True)
class FAGrowthSample:
    """One sample of (per-FA N_engaged, per-FA force-sum, timestep)."""

    timestep: int
    n_engaged: np.ndarray            # shape (n_FAs,) int64
    F_per_FA: np.ndarray             # shape (n_FAs,) float64 [N]


@dataclass
class FAGrowthMonitor:
    """Per-FA growth tracker (read-only observer of the H.4 sim).

    Usage::

        monitor = FAGrowthMonitor(p=p, layouts=layouts)
        for step_block in range(n_blocks):
            sim.run(n_block_steps)
            monitor.sample(sim, integrin_action, timestep=sim.timestep)

    After the run, ``monitor.trigger_frame_for_fa(fa_id)`` returns the
    first sample where ``n_engaged[fa_id] >= N_engaged_target``; the
    KU-2.17 validation test reads ``F_per_FA[fa_id]`` at that frame and
    checks both within tolerance.
    """

    p: "ResolvedH4"
    layouts: "list[FALayout]"
    samples: list[FAGrowthSample] = field(default_factory=list)

    def __post_init__(self) -> None:
        self._n_FAs = len(self.layouts)
        self._N_engaged_target = float(
            self.p.fa_growth["N_engaged_at_threshold"]
        )
        # Map each integrin tag to its FA index (mirror of
        # IntegrinBondUpdater._fa_for_integrin so we don't reach into
        # private state).
        self._fa_for_integrin = np.full(
            sum(L.n_total for L in self.layouts), -1, dtype=np.int64
        )
        for fa_idx, L in enumerate(self.layouts):
            lo = L.integrin_tag_start
            hi = lo + L.n_total
            self._fa_for_integrin[lo:hi] = fa_idx

    @property
    def n_FAs(self) -> int:
        return self._n_FAs

    def sample(
        self,
        sim,
        integrin_action: "IntegrinBondUpdater",
        *,
        timestep: int | None = None,
    ) -> FAGrowthSample:
        """Read N_engaged and F_per_FA at the current sim state.

        Parameters
        ----------
        sim : hoomd.Simulation
            The active simulation. Used to read particle positions for
            the per-FA force aggregate.
        integrin_action : IntegrinBondUpdater
            The dynamic-bond Action; exposes ``n_engaged_per_fa``.
        timestep : int, optional
            If None, reads ``sim.timestep``.

        Returns
        -------
        FAGrowthSample
        """
        if timestep is None:
            timestep = int(sim.timestep)

        # Per-FA engaged count — provided by the integrin Action.
        n_eng = integrin_action.n_engaged_per_fa.astype(np.int64)

        # Per-FA force sum: read the integrin bonds from the snapshot,
        # compute F_i = k_int · |r_int − r_lig|, reduce by FA id.
        F_per_FA = np.zeros(self._n_FAs, dtype=np.float64)
        snap = sim.state.get_snapshot()
        if snap.communicator.rank == 0:
            pos = np.asarray(snap.particles.position, dtype=np.float64)
            bg = np.asarray(snap.bonds.group, dtype=np.int64)
            bt = np.asarray(snap.bonds.typeid, dtype=np.uint32)
            bond_type_names = list(snap.bonds.types)
            try:
                int_type_id = bond_type_names.index("integrin_ligand")
            except ValueError:
                int_type_id = -1
            if int_type_id >= 0 and bg.shape[0] > 0:
                is_int = bt == int_type_id
                int_bonds = bg[is_int]
                if int_bonds.shape[0] > 0:
                    r_int = pos[int_bonds[:, 0]]
                    r_lig = pos[int_bonds[:, 1]]
                    F_mag = self.p.k_int_bare * np.linalg.norm(
                        r_int - r_lig, axis=1
                    )
                    fa_ids = self._fa_for_integrin[int_bonds[:, 0]]
                    np.add.at(F_per_FA, fa_ids, F_mag)

        sample = FAGrowthSample(
            timestep=int(timestep),
            n_engaged=n_eng,
            F_per_FA=F_per_FA,
        )
        self.samples.append(sample)
        return sample

    def trigger_frame_for_fa(self, fa_id: int) -> FAGrowthSample | None:
        """First sample where ``n_engaged[fa_id] >= N_engaged_target``.

        Returns None if no such sample exists in the recorded history —
        which the KU-2.17 test takes as "FA did not reach the growth
        threshold during the run" and surfaces as a gate FAIL with a
        diagnostic about the engagement saturation level achieved.
        """
        target = self._N_engaged_target
        for s in self.samples:
            if s.n_engaged[fa_id] >= target:
                return s
        return None

    def emergent_area(
        self, sample: FAGrowthSample, fa_id: int
    ) -> float:
        """A_emergent = (N_engaged / N_total) · A_baseline  [m²]."""
        L = self.layouts[fa_id]
        baseline = self.p.A_mature if L.is_mature else self.p.A_nascent
        return float(sample.n_engaged[fa_id]) / L.n_total * baseline


def per_fa_force_sum(
    pos: np.ndarray,
    bg: np.ndarray,
    fa_for_integrin: np.ndarray,
    n_FAs: int,
    k_int: float,
) -> np.ndarray:
    """Reduce per-bond F_i = k_int·|r_int−r_lig| onto per-FA totals.

    Helper used by the KU-2.17 validation test when it constructs a
    synthetic engaged-bond population (no full HOOMD run needed for the
    oracle-comparison part).

    Parameters
    ----------
    pos : (N_particles, 3) float64
    bg : (N_int_bonds, 2) int64 — integrin tag in column 0, ligand tag in column 1
    fa_for_integrin : (N_particles,) int64 — FA id for each integrin tag
    n_FAs : int
    k_int : float [N/m]

    Returns
    -------
    F_per_FA : (n_FAs,) float64 [N]
    """
    if bg.shape[0] == 0:
        return np.zeros(n_FAs, dtype=np.float64)
    r_int = pos[bg[:, 0]]
    r_lig = pos[bg[:, 1]]
    F_mag = k_int * np.linalg.norm(r_int - r_lig, axis=1)
    F_per_FA = np.zeros(n_FAs, dtype=np.float64)
    fa_ids = fa_for_integrin[bg[:, 0]]
    np.add.at(F_per_FA, fa_ids, F_mag)
    return F_per_FA
