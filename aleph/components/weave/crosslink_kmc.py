r"""Topology-reforming crosslink KMC (host NumPy oracle) + the single-channel double-count guard — I4.

A crosslinker (alpha-actinin / filamin / fascin) is a dynamic bond: it UNBINDS at a load-dependent off-rate
and REATTACHES to a new near-partner on a DIFFERENT filament found by a neighbour search (a host KD-tree here;
a device ``wp.HashGrid`` on the GPU). Because the rebind partner can differ from the old one, the network
TOPOLOGY reforms over time (the SF/cap/arc bundles condense because crosslinkers hop toward the
myosin-aligned, lower-energy configuration) while the bound crosslinker COUNT stays statistically steady. This
is the crosslinker-relaxation mechanism that makes the woven network alive at I4.

RATE LAWS (identical closed forms to the I3 hand API, ``ac.motor.hand``):
    p_off(f) = k_off0 * exp(|f|/f0)          (Bell 1978 slip; f0 = kT/x_beta)
    p_detach = 1 - exp(-tau * p_off)         (Poisson per tick)
    p_attach = 1 - exp(-tau * k_on)          (Poisson per tick; a free end within ``reach`` of a cross-fiber
                                              node binds)
At the unloaded steady state the bound fraction -> k_on/(k_on + k_off0) (detailed-balance check).

DOUBLE-COUNT GUARD (build-plan §6 instance "I2-crosslinker-relaxation"; the I4 named guard). There are TWO
candidate crosslinker-relaxation channels and running BOTH counts crosslinker fluidization twice:
  * ``"kmc_reattach"``  — THIS module: unbind + rebind to a NEW partner (topology reforms).
  * ``"r0_creep"``      — ``ff.motility_warp.xl_turnover_kernel``: a rebinding crosslink is force-free at its
                          CURRENT length, so its rest length Maxwell-creeps ``rest += frac*(L - rest)``
                          (topology FIXED, stress relaxes).
``CrosslinkRelaxationChannel`` enforces EXACTLY ONE is active; :func:`assert_single_channel` refuses both.
I4 selects ``kmc_reattach`` (it reforms topology, which the r0-creep cannot — needed for bundle condensation);
the lead must then DISABLE ``xl_turnover_kernel`` in the native driver (INTEGRATION.md patch-note).

engine units: force pN, rate 1/s, length um, time s.

Sanity Gate (self-tested in tests/ac/weave/test_crosslink_kmc_oracle.py):
  * dimensional/bounds: p_detach, p_attach in [0,1]; Bell off-rate monotone increasing in |load|.
  * detailed balance: unloaded steady bound fraction -> k_on/(k_on+k_off0) within stochastic spread.
  * conservation: no crosslinker created/destroyed — a detached crosslinker rebinds (or waits); the total
    crosslinker count is invariant across a run.
  * topology reforms: after many ticks the bound-pair SET differs from the initial set (a non-trivial fraction
    of bonds have hopped) while the bound COUNT is steady -> genuine reattachment, not just on/off flicker.
  * single-channel guard: enabling both ``kmc_reattach`` and ``r0_creep`` raises (double-count refused).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from aleph.components.weave.branch_angle import KT_310K_PN_UM

__all__ = [
    "CrosslinkRelaxationChannel",
    "assert_single_channel",
    "bell_off_rate",
    "detach_prob",
    "attach_prob",
    "steady_bound_fraction",
    "reattach_step",
]

_CHANNELS = ("kmc_reattach", "r0_creep")


@dataclass(frozen=True, slots=True)
class CrosslinkRelaxationChannel:
    """The single-selected crosslinker-relaxation channel (double-count guard, build-plan §6).

    Exactly one of the two mechanisms may be active. I4 selects ``kmc_reattach``; the competing
    ``r0_creep`` (``ff.motility_warp.xl_turnover_kernel``) must then be OFF in the native driver.
    """

    name: str

    def __post_init__(self) -> None:
        if self.name not in _CHANNELS:
            raise ValueError(f"channel must be one of {_CHANNELS}; got {self.name!r}")


def assert_single_channel(kmc_reattach: bool, r0_creep: bool) -> CrosslinkRelaxationChannel:
    """Return the single active channel, or RAISE if both (or neither) are on — the double-count guard.

    Args:
        kmc_reattach: whether the topology-reforming KMC (this module) is active.
        r0_creep: whether the ff ``xl_turnover_kernel`` rest-length creep is active.

    Raises:
        ValueError: if both are active (double count) or neither is (crosslinkers frozen).
    """
    if kmc_reattach and r0_creep:
        raise ValueError(
            "double-count: both crosslinker-relaxation channels active "
            "(kmc_reattach AND r0_creep/xl_turnover_kernel). Build-plan §6 requires EXACTLY ONE."
        )
    if not kmc_reattach and not r0_creep:
        raise ValueError("no crosslinker-relaxation channel active (crosslinkers would be permanent springs)")
    return CrosslinkRelaxationChannel("kmc_reattach" if kmc_reattach else "r0_creep")


def bell_off_rate(load: float | npt.NDArray[np.float64], k_off0: float, f0: float) -> npt.NDArray[np.float64]:
    """Bell slip off-rate ``k_off0 * exp(|load|/f0)`` [1/s] (matches ``ac.motor.hand.bell_off_rate``)."""
    if k_off0 < 0.0 or f0 <= 0.0:
        raise ValueError("k_off0 >= 0 and f0 > 0 required")
    return k_off0 * np.exp(np.abs(np.asarray(load, float)) / f0)


def detach_prob(tau: float, p_off: float | npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Per-tick detach probability ``1 - exp(-tau p_off)``."""
    return 1.0 - np.exp(-tau * np.asarray(p_off, float))


def attach_prob(tau: float, k_on: float) -> float:
    """Per-tick attach probability ``1 - exp(-tau k_on)`` (Poisson)."""
    return float(1.0 - np.exp(-tau * k_on))


def steady_bound_fraction(k_on: float, k_off: float) -> float:
    """Unloaded steady-state bound fraction ``k_on / (k_on + k_off)`` of a two-state crosslinker."""
    if k_on < 0.0 or k_off < 0.0 or (k_on + k_off) <= 0.0:
        raise ValueError("k_on, k_off >= 0 and not both zero")
    return float(k_on / (k_on + k_off))


def reattach_step(
    bound: npt.NDArray[np.bool_],
    partner: npt.NDArray[np.int64],
    pos: npt.NDArray[np.float64],
    fiber_of_node: npt.NDArray[np.int64],
    node_of_end: npt.NDArray[np.int64],
    loads: npt.NDArray[np.float64],
    *,
    reach: float,
    k_on: float,
    k_off0: float,
    f0: float,
    tau: float,
    rng: np.random.Generator,
) -> tuple[npt.NDArray[np.bool_], npt.NDArray[np.int64]]:
    """One topology-reforming KMC tick: bound crosslinkers Bell-detach, free ones rebind to a NEW cross-fiber node.

    A crosslinker has one moving end anchored at node ``node_of_end[e]``; when bound it holds ``partner[e]``
    (a node on a DIFFERENT filament) with load ``loads[e]``. Per tick: (1) each bound end detaches with
    ``detach_prob(tau, bell_off_rate(load))``; (2) each free end rebinds to the nearest cross-fiber node
    within ``reach`` with ``attach_prob(tau, k_on)``. The rebind partner may DIFFER from the old one -> the
    topology reforms. Crosslinker count is conserved (ends are never created/destroyed).

    Args:
        bound: (E,) bound flags (mutated copy returned).
        partner: (E,) current partner node index (or -1 if free) (mutated copy returned).
        pos: (N, 3) node positions [um].
        fiber_of_node: (N,) fiber id of every node (for the cross-fiber constraint).
        node_of_end: (E,) the anchored moving-end node of each crosslinker.
        loads: (E,) current per-crosslinker load [pN] (0 for free ends).
        reach: crosslinker capture radius [um].
        k_on, k_off0, f0: attach rate, zero-load off-rate, Bell force [1/s,1/s,pN].
        tau: KMC tick [s].
        rng: NumPy generator.

    Returns:
        Updated ``(bound, partner)``.
    """
    bound = np.asarray(bound, bool).copy()
    partner = np.asarray(partner, np.int64).copy()
    pos = np.asarray(pos, float)
    node_of_end = np.asarray(node_of_end, np.int64)
    bound_start = bound.copy()   # detach and reattach are SEPARATE KMC events in a tick (no same-tick rebind)

    # (1) Bell-detach ends bound at the start of the tick.
    is_bound = np.where(bound_start)[0]
    if is_bound.size:
        p_off = bell_off_rate(loads[is_bound], k_off0, f0)
        p_det = detach_prob(tau, p_off)
        det = rng.random(is_bound.size) < p_det
        bound[is_bound[det]] = False
        partner[is_bound[det]] = -1

    # (2) rebind ends free at the START of the tick to a NEW cross-fiber node within reach, chosen
    # STOCHASTICALLY among the eligible neighbours (not deterministically the nearest) so the partner set
    # genuinely turns over -> the topology reforms (the device HashGrid pass picks a candidate the same way).
    p_att = attach_prob(tau, k_on)
    free = np.where(~bound_start)[0]
    for e in free:
        if rng.random() >= p_att:
            continue
        a = int(node_of_end[e])
        d = np.linalg.norm(pos - pos[a], axis=1)
        eligible = (fiber_of_node != fiber_of_node[a]) & (d <= reach)
        eligible[a] = False
        cand = np.where(eligible)[0]
        if cand.size == 0:
            continue
        partner[e] = int(rng.choice(cand))
        bound[e] = True
    return bound, partner
