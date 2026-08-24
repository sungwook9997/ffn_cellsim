r"""α-actinin bind/unbind for a transient filament crosslink — the delegate, not a datum.

WHY THIS IS NOT A PI DECISION.  ``sf_cortex_transient`` declares ``kinetics=True``, and the runtime
refuses to bind it without a rate law because *a transient crosslink that can bind and never unbind is
a permanent weld*, which ``CLAUDE.md`` forbids by name.  That refusal fired on 2026-08-11 (job 91) and
looked like another PI-GAP.  **It is not.**  Every constant the law needs is already SOURCED on
:data:`~aleph.laws.hand_kmc.ALPHA_ACTININ`:

  ``k_on``               10 /s        binding rate with an acceptor in range
  ``p0`` (``k_off0``)    0.066 /s     zero-force off-rate prefactor — Ferrer 2008 PNAS AFM
  ``f0``                 Bell force   from ``x_β`` = 0.4 nm via :func:`bell_f0_from_x_beta`
  ``capture_radius_um``  0.06 µm      the same ε the pairing topology already uses
  ``link_k``             4.6e5 pN/µm  PI-approved 2026-06-30

So what was missing is a delegate.  Writing one costs no decision, and surfacing it to PI would have
spent one on a value that exists — which is exactly the sorting
``PI_DECISION_CARD_COMPOSITION_2026-08-11.md`` was written to make possible.

WHAT THE LAW IS.  A slip bond, Bell form, the repository's own:
``k_off(f) = p0·exp(f/f0)`` via :meth:`HandParams.off_rate`, and per-tick probabilities via
:func:`attach_probability` / :func:`detach_probability` — ``1 − exp(−k·dt)``, not ``k·dt``, so a large
``dt`` cannot produce a probability above one.  Nothing is linearised and no rate is capped.

WHAT IT PROPOSES, AND WHAT IT MAY NOT DO.  It writes ``candidate_state_d`` **only**.  ``state_d`` is
untouched: a kinetic connector commits ONLY on an accepted physical step, under the transaction, and
this delegate is not the transaction.  A rejected candidate therefore leaves no trace, which is the
property that makes the proposal reversible at all.

The device RNG is ``wp.rand_init(rng_seed, joint)`` — one stream per joint, seeded by the caller's
accepted-step seed, the same construction ``events.py`` and ``microtubule_rig`` use.  No host draw, no
host state written, nothing read back inside the loop.

⚠ **WHAT THIS DOES NOT DECIDE.**  How MANY crosslinks there should be is not here and is not sourced:
the joint population comes from geometry inside the capture radius, and the SF inventory axis is a
standing PI-GAP.  This law governs the state of the joints that exist; it does not license their count.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — force pN, rates 1/s, ``dt`` s, probabilities dimensionless; ``f0`` is a force and is
    read from the sourced record rather than formed here.
  * boundary — ``dt = 0`` proposes nothing; an inactive joint is skipped; the probability form is
    saturating, so no ``dt`` yields ``p > 1``.
  * conservation/invariant — only ``candidate_state_d`` is written, so an unaccepted step is a no-op
    and the committed state cannot drift from the transaction.
  * CFL/precision — no integration; float64 throughout.
  * sign sense — the off-rate RISES with load (slip bond): a loaded crosslink is more likely to
    release, never less.  A catch bond would need the Pereverzev branch and α-actinin is not one
    (``catch_slip=False`` on the sourced record), so this refuses to run for a catch-slip hand rather
    than silently applying the wrong law.
  * measurement protocol — device-resident; the caller supplies the accepted-step seed.

engine units: force pN, time s, length µm.  Runtime: NVIDIA Warp on CUDA only.
"""

from __future__ import annotations

from dataclasses import dataclass

import warp as wp

from aleph.laws.hand_kmc import ALPHA_ACTININ, HandParams

__all__ = ["AlphaActininKinetics", "build_alpha_actinin_kinetics"]

#: ``JointState.FREE`` / ``JointState.ACTIN_ENGAGED`` as device constants — the two states a transient
#: actin crosslink occupies.  Mirrored here rather than imported into kernel scope, which Warp cannot do.
_FREE = wp.constant(0)
_ENGAGED = wp.constant(2)


@wp.kernel
def _alpha_actinin_propose_kernel(
    active: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    load: wp.array(dtype=wp.float64),
    candidate_state: wp.array(dtype=wp.int32),
    k_on: wp.float64,
    p0: wp.float64,
    f0: wp.float64,
    dt: wp.float64,
    rng_seed: wp.int32,
) -> None:
    """Propose one joint's bind/unbind for this candidate. Writes ``candidate_state`` only."""
    j = wp.tid()
    if active[j] == 0:
        return

    current = state[j]
    candidate_state[j] = current  # default: unchanged, so a skipped branch is not a silent FREE

    rstate = wp.rand_init(rng_seed, j)
    u = wp.float64(wp.randf(rstate))

    if current == _ENGAGED:
        # Bell slip bond: the off-rate rises with the load the joint is carrying.
        k_off = p0 * wp.exp(wp.abs(load[j]) / f0)
        p_detach = wp.float64(1.0) - wp.exp(-k_off * dt)
        if u < p_detach:
            candidate_state[j] = _FREE
    elif current == _FREE:
        p_attach = wp.float64(1.0) - wp.exp(-k_on * dt)
        if u < p_attach:
            candidate_state[j] = _ENGAGED


@dataclass(frozen=True, slots=True)
class AlphaActininKinetics:
    """Bell slip-bond bind/unbind for α-actinin crosslinks, proposing onto the candidate state.

    Attributes:
        params: the SOURCED hand record. Defaults to :data:`ALPHA_ACTININ`; a caller may substitute a
            different sourced crosslinker but may not hand-assemble rates.
        launch: injected launcher, so a CUDA-free structural gate can record launches.
    """

    params: HandParams = ALPHA_ACTININ
    launch: object = wp.launch

    def __post_init__(self) -> None:
        if self.params.catch_slip:
            raise ValueError(
                f"{self.params.name!r} is a catch-slip hand; this delegate implements the Bell slip "
                "law only. Applying it would be the wrong rate law wearing the right name."
            )
        if not (self.params.p0 > 0.0 and self.params.f0 > 0.0 and self.params.k_on > 0.0):
            raise ValueError(
                f"{self.params.name!r} must carry positive k_on, p0 and f0; got "
                f"k_on={self.params.k_on!r}, p0={self.params.p0!r}, f0={self.params.f0!r}"
            )

    @property
    def event_channels(self) -> frozenset[str]:
        """The channels this law can propose on."""
        return frozenset({"crosslink_bind", "crosslink_unbind"})

    def propose_events(
        self,
        joints: object,
        rates: object = None,
        dt_phys: float = 0.0,
        rng_seed: int = 0,
        neighbors: object = None,
    ) -> None:
        """Propose bind/unbind for every active joint. Writes ``candidate_state_d`` and nothing else.

        The POSITIONAL order is the transaction's, not this module's choice:
        ``CellTransaction.step`` calls ``proposer.propose_events(rates, dt_phys, rng_seed, neighbors)``
        and ``FilamentCrosslinkConnector`` prepends its joint runtime.  Getting it wrong is silent —
        ``rates`` arrives where ``dt_phys`` is expected and the first native run dies on a ``None``
        comparison, which is exactly what job 93 did.  A keyword-only test would not have caught it,
        and did not; ``test_the_positional_order_matches_the_transaction`` now does.

        Args:
            joints: the ``LoadPathJointRuntime`` holding ``active_d`` / ``state_d`` / ``load_d`` /
                ``candidate_state_d``.
            rates: the transaction's rate bundle. Unused here — α-actinin's rates are its own sourced
                record, not a per-step input — and accepted so the positional contract holds.
            dt_phys: the physical tick [s]. ``0`` proposes nothing rather than dividing or defaulting.
            rng_seed: the accepted-step seed; one device stream per joint is derived from it.
            neighbors: the transaction's neighbour bundle. Unused: this edge's joints are a fixed
                build-time pairing, so there is no re-query per step.
        """
        if dt_phys < 0.0:
            raise ValueError("dt_phys must be nonnegative")
        if dt_phys == 0.0:
            return
        capacity = int(joints.active_d.shape[0])
        if capacity == 0:
            return
        self.launch(
            _alpha_actinin_propose_kernel,
            dim=capacity,
            inputs=[
                joints.active_d, joints.state_d, joints.load_d, joints.candidate_state_d,
                wp.float64(self.params.k_on), wp.float64(self.params.p0),
                wp.float64(self.params.f0), wp.float64(dt_phys), wp.int32(rng_seed),
            ],
            device=str(joints.active_d.device),
        )


def build_alpha_actinin_kinetics(*, params: HandParams | None = None) -> AlphaActininKinetics:
    """Return the α-actinin rate law, refusing anything whose rates are not sourced."""
    return AlphaActininKinetics(params=params or ALPHA_ACTININ)
