"""Typed global cell-engine ledger — the accepted-step accounting sink.

Every :class:`~aleph.engine.runtime.LedgerContributor` in the composed world (each component
state-owner, each connector runtime, and the far-field boundary) pushes its device-resident force /
work / mass / topology terms into ONE :class:`GlobalCellLedger` before the outer acceptance predicate
is evaluated.  The ledger owns the running accumulators; it never decides acceptance and never reads an
authoritative device value on the host.

This class is the concrete typed home of the ``add_far_field_reaction`` slot the ECM far-field boundary
runtime (:mod:`aleph.engine.ecm_world`) forwards to — previously a documented duck-typed gap.  Two
force channels are kept *independently sourced* so the global force-balance gate is a genuine check, not a
tautology:

* ``reaction_resultant`` — Σ of the far-field Dirichlet reaction the truncated ECM boundary resists
  (contributed by the boundary runtime via :meth:`add_far_field_reaction`).
* ``traction_resultant`` — Σ of the cell→ECM clutch/contact load actually scattered onto the ECM body
  (contributed by the composite α2β1–collagen joint via :meth:`add_cell_traction`).

For an inertialess (force-equilibrium) ECM sub-body cut at the FA clutches, Newton's third law plus
equilibrium require ``reaction_resultant + traction_resultant = 0``.  :meth:`assemble_balance` computes that
residual on the device and writes a device flag against a **caller-supplied** tolerance — so "sum to zero" is
CHECKED, never assumed, and no numerical tolerance is baked in here.

WHAT THE TWO CHANNELS ARE, IN GENERAL.  The ECM cut above is the canonical instance, not the definition.
The invariant the gate encodes is *two independently accumulated force channels of ONE Newton pair sum
to the zero vector*, and any cut of the composed cell into two sub-bodies instantiates it.  A two-body
vertical slice does so with its two never-merged force arrays — for ``nmii_sf_motor``, the ``sf_arc``
node array and the ``nmii`` particle array, whose only coupling is the split crossbridge — and
:meth:`add_body_force` is the accessor for that reading.  What must NOT happen is the two channels being
computed from one another (``traction := −reaction``), which turns a check into a tautology; they have
to come from separate accumulations, which is why the accessor takes a body's own force array.

The generic force / work / mass / topology slots (:meth:`add_force_resultant`, :meth:`add_work`,
:meth:`add_mass`, :meth:`add_topology_delta`) are the accumulation surface a ``LedgerContributor`` implies.

Sanity Gate:
    * ownership: every accumulator is a distinct-storage CUDA array on one device; no host authoritative read.
    * dimensional: reaction/traction/force resultants are vec3d [pN]; work/mass are float64 scalars; topology
      counts are int32 deltas.
    * conservation: :meth:`add_far_field_reaction` reduces the reaction array AND the boundary-work scalar in
      one call; the two force channels stay separately sourced so the balance gate cannot self-satisfy.
    * boundary/sign: the balance residual is ``|reaction + traction|`` — equal-and-opposite resultants cancel
      to the zero vector, a sign error does not.
    * numerical: the "sum to zero" tolerance is a caller-owned device scalar (no magic number); the gate writes
      a device int flag, never a Python branch on device state.  Its VALUE is derived, not chosen (PI D8,
      2026-07-28): :func:`assemble_balance_tolerance_ratio` returns the float64 summation bound ``γ_n`` for the
      contribution COUNT, and :meth:`GlobalCellLedger.derive_balance_tolerance` scales it on the device by the
      ledger's own ``|reaction| + |traction|``.  The gate therefore accepts no supplied force magnitude at all
      — a tolerance that took one would have relocated the magic number instead of removing it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

__all__ = ["GlobalCellLedger", "make_global_cell_ledger"]


def _device_is_cuda(array: object) -> bool:
    """Return whether an array-like object declares a CUDA device."""
    device = getattr(array, "device", None)
    return bool(getattr(device, "is_cuda", False))


def _storage_key(array: object) -> tuple[object, ...]:
    """Return a conservative storage identity without copying device data."""
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


def _validate_device_array(
    array: object,
    *,
    label: str,
    dtype: object,
    ndim: int,
    nonempty: bool = True,
) -> None:
    """Validate CUDA array metadata without reading authoritative state."""
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != dtype:
        raise TypeError(f"{label} must have dtype {dtype}")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != ndim:
        raise ValueError(f"{label} must be a {ndim}-dimensional array")
    if nonempty and any(int(size) <= 0 for size in shape):
        raise ValueError(f"{label} must have no empty dimension")


def _require_same_device(arrays: tuple[object, ...], *, label: str) -> None:
    """Require all arrays to occupy one CUDA device."""
    devices = {str(getattr(array, "device", None)) for array in arrays}
    if len(devices) != 1:
        raise ValueError(f"{label} arrays must share one CUDA device")


def _require_distinct_storage(arrays: tuple[object, ...], *, label: str) -> None:
    """Reject authoritative accumulators that alias device storage."""
    keys = tuple(_storage_key(array) for array in arrays)
    if len(set(keys)) != len(keys):
        raise ValueError(f"{label} arrays must own distinct storage")


@wp.kernel
def _reduce_vec3_kernel(
    src: wp.array(dtype=wp.vec3d),
    dst: wp.array(dtype=wp.vec3d),
) -> None:
    """Atomically add every ``src`` vector into the single ``dst[0]`` resultant accumulator."""
    t = wp.tid()
    wp.atomic_add(dst, 0, src[t])


@wp.kernel
def _reduce_scalar_kernel(
    src: wp.array(dtype=wp.float64),
    dst: wp.array(dtype=wp.float64),
) -> None:
    """Atomically add every ``src`` scalar into the single ``dst[0]`` accumulator."""
    t = wp.tid()
    wp.atomic_add(dst, 0, src[t])


@wp.kernel
def _reduce_int_counts_kernel(
    src: wp.array(dtype=wp.int32),
    dst: wp.array(dtype=wp.int32),
) -> None:
    """Add per-slot topology-count deltas index-wise into the running count accumulator."""
    t = wp.tid()
    wp.atomic_add(dst, t, src[t])


#: Unit round-off for float64: ``eps64 / 2``.  Not a tunable — it is a property of the format.
_U_FLOAT64: float = float(np.finfo(np.float64).eps) / 2.0


def assemble_balance_tolerance_ratio(n_terms: int) -> float:
    """Return the float64 accumulation floor ``γ_n`` for a sum of ``n_terms`` force contributions.

    This is the DERIVATION behind the ``assemble_balance`` tolerance (PI decision **D8**, approved
    2026-07-28; the derivation was pre-approved 2026-07-25, the value needed sign-off).  What is
    signed is this formula, not a number — the same principle as D4: never freeze a constant, let the
    run generate its own scale.

    ``reaction_resultant + traction_resultant`` is exactly zero in exact arithmetic (Newton's third law
    on an inertialess sub-body).  Any nonzero residual is therefore floating-point accumulation error,
    and the tolerance is simply a bound on that error.  For summation of ``n`` terms in a format with
    unit round-off ``u``, the standard forward bound (Higham, *Accuracy and Stability of Numerical
    Algorithms*, §4.2) is ``|Σ̂ − Σ| ≤ γ_n · Σ|x_i|`` with::

        γ_n = (n − 1)·u / (1 − (n − 1)·u)

    The **order-independent** bound is the right one here and not a conservative habit: both resultants
    are accumulated with ``wp.atomic_add``, whose arrival order is nondeterministic, so no summation
    order may be assumed (this is the same property that makes bit-identity unavailable as an
    acceptance test — see the plan's §8 defect 8).

    The ``Σ|x_i|`` factor is deliberately NOT a parameter of this function.  It is supplied on the
    device by :meth:`GlobalCellLedger.derive_balance_tolerance` from ``|reaction| + |traction|``, which
    the ledger already holds — so the gate takes **no supplied magnitude at all**, only this count.
    A tolerance that needed a magnitude would have relocated the magic number rather than removed it.

    Args:
        n_terms: Number of accumulated force contributions across BOTH channels — pinned far-field
            nodes plus clutch/contact endpoints.  A count, never a magnitude, and known exactly at
            build time.  Must be at least 1.

    Returns:
        The dimensionless ratio ``γ_n``.  The gate is then ``|reaction + traction| ≤ γ_n · (|reaction| +
        |traction|)``, i.e. a RELATIVE criterion carrying no force unit.

    Raises:
        ValueError: If ``n_terms`` is below 1, or so large that ``(n−1)·u ≥ 1`` and the bound has no
            meaning — at float64 that needs ~9e15 terms, so it reports a wiring error, not a real run.

    Note:
        This bounds ROUND-OFF only.  It does not bound discretisation error, an unconverged solve, or a
        missing force channel — a run can sit far above this floor for entirely physical reasons, and
        that is the gate doing its job.
    """
    n = int(n_terms)
    if n < 1:
        raise ValueError("n_terms must count at least one accumulated force contribution")
    slack = float(n - 1) * _U_FLOAT64
    if slack >= 1.0:
        raise ValueError(
            f"n_terms={n} makes the float64 summation bound vacuous ((n-1)*u >= 1); "
            "this is a wiring error, not a physical population"
        )
    return slack / (1.0 - slack)


@wp.kernel
def _reduce_vec3_magnitude_kernel(
    source: wp.array(dtype=wp.vec3d), total: wp.array(dtype=wp.float64)
) -> None:
    """Accumulate ``Σ|x_i|`` [pN] — the Higham bound's own scale, not the resultant's magnitude."""
    i = wp.tid()
    wp.atomic_add(total, 0, wp.length(source[i]))


@wp.kernel
def _derive_balance_tolerance_kernel(
    reaction_resultant: wp.array(dtype=wp.vec3d),
    traction_resultant: wp.array(dtype=wp.vec3d),
    accumulated_magnitude: wp.array(dtype=wp.float64),
    gamma_n: wp.array(dtype=wp.float64),
    tol_sq: wp.array(dtype=wp.float64),
) -> None:
    """Write the SQUARED balance tolerance from the ledger's own accumulated magnitudes.

    One thread (``dim=1``).  The Higham bound is ``|Σ̂ − Σ| ≤ γ_n · Σ|x_i|`` and the scale in it is the
    sum of the terms' MAGNITUDES, which is not the magnitude of their sum: a channel whose contributions
    largely cancel — a bipolar minifilament is a force dipole, so its resultant is near zero while its
    traffic is not — has ``Σ|x_i| ≫ |Σx_i|``, and bracketing the one by the other would set a tolerance
    orders of magnitude below the arithmetic error it exists to admit, rejecting every physically correct
    step.  So ``accumulated_magnitude`` carries ``Σ|x_i|`` when a contributor accumulated it, and the
    resultant bracket ``|reaction| + |traction|`` remains the fallback for channels that do not.

    Either way the gate takes no SUPPLIED force scale: both come from the ledger's own contributions.
    Squared because :meth:`GlobalCellLedger.assemble_balance` compares squared norms and never needs a
    square root.
    """
    scale = accumulated_magnitude[0]
    if scale <= wp.float64(0.0):
        scale = wp.length(reaction_resultant[0]) + wp.length(traction_resultant[0])
    tol = gamma_n[0] * scale
    tol_sq[0] = tol * tol


@wp.kernel
def _assemble_force_balance_kernel(
    reaction_resultant: wp.array(dtype=wp.vec3d),
    traction_resultant: wp.array(dtype=wp.vec3d),
    tol_sq: wp.array(dtype=wp.float64),
    balance_residual_sq: wp.array(dtype=wp.float64),
    balance_ok: wp.array(dtype=wp.int32),
) -> None:
    """Compute ``|reaction + traction|²`` and set the device gate flag against a caller tolerance.

    One thread (``dim=1``).  The residual is the resultant of the two independently sourced force channels;
    at equilibrium Newton's third law drives it to the zero vector.  ``balance_ok[0]`` is written on the
    device — the acceptance predicate reads it there, never on the host.
    """
    r = reaction_resultant[0] + traction_resultant[0]
    n2 = wp.dot(r, r)
    balance_residual_sq[0] = n2
    if n2 <= tol_sq[0]:
        balance_ok[0] = wp.int32(1)
    else:
        balance_ok[0] = wp.int32(0)


@dataclass(frozen=True, slots=True)
class GlobalCellLedger:
    """Device-resident accepted-step accounting sink for the composed cell world.

    Args:
        reaction_resultant_d: ``(1,)`` vec3d — Σ far-field boundary reaction [pN].
        traction_resultant_d: ``(1,)`` vec3d — Σ cell→ECM clutch/contact load scattered onto the ECM [pN].
        force_resultant_d: ``(1,)`` vec3d — generic net-force accumulator [pN] (any contributor's resultant).
        work_d: ``(1,)`` float64 — accumulated boundary + internal work [pN·µm].
        mass_d: ``(1,)`` float64 — accumulated fluid/monomer mass term (conservation channel).
        topology_counts_d: ``(K,)`` int32 — running topology-count deltas (bonds formed/broken, remeshes …).
        balance_residual_sq_d: ``(1,)`` float64 — last ``|reaction + traction|²`` written by the gate.
        balance_ok_d: ``(1,)`` int32 — device gate flag (1 = balanced within the caller tolerance).
        balance_scale_d: optional ``(1,)`` float64 — ``Σ|x_i|`` over the contributions that fed the two
            balance channels, i.e. the Higham bound's OWN scale.  Optional so an existing construction is
            unaffected: when absent, :meth:`derive_balance_tolerance` falls back to bracketing the scale
            by ``|reaction| + |traction|``, which is only adequate for a channel whose contributions do
            not largely cancel.  :meth:`add_body_force` fills it, because a force dipole's resultant is
            near zero while its arithmetic traffic is not.

    All accumulators must be distinct-storage CUDA arrays on one device.  No method reads an authoritative
    device value on the host; the reductions and the balance gate run entirely in Warp kernels.
    """

    reaction_resultant_d: wp.array
    traction_resultant_d: wp.array
    force_resultant_d: wp.array
    work_d: wp.array
    mass_d: wp.array
    topology_counts_d: wp.array
    balance_residual_sq_d: wp.array
    balance_ok_d: wp.array
    balance_scale_d: wp.array | None = None

    def __post_init__(self) -> None:
        vec_scalars = (
            ("reaction_resultant_d", self.reaction_resultant_d),
            ("traction_resultant_d", self.traction_resultant_d),
            ("force_resultant_d", self.force_resultant_d),
        )
        for label, array in vec_scalars:
            _validate_device_array(array, label=f"ledger.{label}", dtype=wp.vec3d, ndim=1)
            if array.shape != (1,):
                raise ValueError(f"ledger.{label} must be a single-entry resultant accumulator")
        f64_scalars = (
            ("work_d", self.work_d),
            ("mass_d", self.mass_d),
            ("balance_residual_sq_d", self.balance_residual_sq_d),
        )
        for label, array in f64_scalars:
            _validate_device_array(array, label=f"ledger.{label}", dtype=wp.float64, ndim=1)
            if array.shape != (1,):
                raise ValueError(f"ledger.{label} must hold exactly one accepted scalar")
        _validate_device_array(
            self.topology_counts_d, label="ledger.topology_counts_d", dtype=wp.int32, ndim=1
        )
        _validate_device_array(
            self.balance_ok_d, label="ledger.balance_ok_d", dtype=wp.int32, ndim=1
        )
        if self.balance_ok_d.shape != (1,):
            raise ValueError("ledger.balance_ok_d must hold exactly one device gate flag")
        if self.balance_scale_d is not None:
            _validate_device_array(
                self.balance_scale_d, label="ledger.balance_scale_d", dtype=wp.float64, ndim=1
            )
            if self.balance_scale_d.shape != (1,):
                raise ValueError("ledger.balance_scale_d must hold exactly one accumulated scale")
        arrays = tuple(
            array for array in (
                self.reaction_resultant_d,
                self.traction_resultant_d,
                self.force_resultant_d,
                self.work_d,
                self.mass_d,
                self.topology_counts_d,
                self.balance_residual_sq_d,
                self.balance_ok_d,
                self.balance_scale_d,
            ) if array is not None
        )
        _require_same_device(arrays, label="global cell ledger")
        _require_distinct_storage(arrays, label="global cell ledger")

    @property
    def device(self) -> str:
        """CUDA device shared by every ledger accumulator."""
        return str(self.reaction_resultant_d.device)

    @property
    def n_topology_slots(self) -> int:
        """Number of topology counters from host-visible array metadata."""
        return int(self.topology_counts_d.shape[0])

    def reset(self) -> None:
        """Zero every accumulator before a fresh ledger pass (device memset, no host read)."""
        for array in (
            self.reaction_resultant_d,
            self.traction_resultant_d,
            self.force_resultant_d,
            self.work_d,
            self.mass_d,
            self.topology_counts_d,
            self.balance_residual_sq_d,
            self.balance_ok_d,
            *((self.balance_scale_d,) if self.balance_scale_d is not None else ()),
        ):
            array.zero_()

    def add_far_field_reaction(self, reaction_d: wp.array, work_d: wp.array) -> None:
        """Reduce the far-field boundary reaction resultant and boundary work into the ledger.

        ``reaction_d`` is the boundary runtime's committed per-pinned-node reaction ``R = -f`` (vec3d); its
        resultant accumulates into :attr:`reaction_resultant_d`.  ``work_d`` is the ``(1,)`` boundary-work
        scalar.  This is the concrete slot the ECM far-field runtime forwards to.
        """
        _validate_device_array(reaction_d, label="reaction_d", dtype=wp.vec3d, ndim=1)
        _validate_device_array(work_d, label="work_d", dtype=wp.float64, ndim=1)
        _require_same_device(
            (reaction_d, work_d, self.reaction_resultant_d), label="far-field reaction ledger"
        )
        wp.launch(
            _reduce_vec3_kernel,
            dim=int(reaction_d.shape[0]),
            inputs=[reaction_d, self.reaction_resultant_d],
            device=self.device,
        )
        wp.launch(
            _reduce_scalar_kernel,
            dim=int(work_d.shape[0]),
            inputs=[work_d, self.work_d],
            device=self.device,
        )

    def add_cell_traction(self, traction_d: wp.array) -> None:
        """Reduce the cell→ECM clutch/contact load resultant into the traction channel.

        ``traction_d`` is the per-endpoint force the composite α2β1–collagen joint (or membrane contact)
        scatters onto the ECM body [pN, vec3d].  Kept separate from the reaction channel so the balance gate
        is Newton's-third-law-genuine.
        """
        _validate_device_array(traction_d, label="traction_d", dtype=wp.vec3d, ndim=1)
        _require_same_device(
            (traction_d, self.traction_resultant_d), label="cell traction ledger"
        )
        wp.launch(
            _reduce_vec3_kernel,
            dim=int(traction_d.shape[0]),
            inputs=[traction_d, self.traction_resultant_d],
            device=self.device,
        )

    def add_body_force(self, force_d: wp.array, *, side: str) -> None:
        """Reduce ONE sub-body's whole force-array resultant into one side of the balance gate.

        The gate's invariant is that its two channels are independently accumulated halves of one Newton
        pair (see this module's header).  A vertical slice cut into two bodies whose only coupling is a
        connector satisfies it exactly: every force is either internal to one array — where the pairs
        cancel within that array's own resultant — or an adjoint pair scattered across the two, so the
        two resultants are equal and opposite whatever the configuration and however unconverged the
        solve.  Which is the point: the gate then tests the ADJOINT WIRING (Newton's third law across the
        connector), not the convergence, and it fails exactly when a scatter is one-sided, double
        counted, or sign-flipped.

        No new accumulator: this writes the same ``reaction_resultant_d`` / ``traction_resultant_d`` the
        ECM lane fills through its own named slots, so one lane's balance cannot be assembled from
        channels the other has meanwhile filled with a different meaning.

        Args:
            force_d: The body's per-node force array (``vec3d`` [pN]) — the SAME array its kernels
                accumulate into and its integrator reads, never a copy computed from the other side.
            side: ``"reaction"`` or ``"traction"`` — which channel this body occupies.  Arbitrary but
                fixed per lane; the gate is symmetric in the two.

        Raises:
            ValueError: If the array metadata is wrong, the device differs, or ``side`` is not one of
                the two channel names.
        """
        targets = {"reaction": self.reaction_resultant_d, "traction": self.traction_resultant_d}
        if side not in targets:
            raise ValueError(f"side must be 'reaction' or 'traction'; got {side!r}")
        target = targets[side]
        _validate_device_array(force_d, label=f"body_force_d[{side}]", dtype=wp.vec3d, ndim=1)
        _require_same_device((force_d, target), label=f"{side} body-force ledger")
        wp.launch(
            _reduce_vec3_kernel,
            dim=int(force_d.shape[0]),
            inputs=[force_d, target],
            device=self.device,
        )
        if self.balance_scale_d is not None:
            # Σ|f_i| for the SAME contributions, because that — not |Σf_i| — is the scale of the float64
            # summation bound this channel's residual will be judged against.  A body whose forces largely
            # cancel (any dipole) has a resultant far below its arithmetic traffic, and judging its
            # round-off against the resultant would reject every correct step.
            wp.launch(
                _reduce_vec3_magnitude_kernel,
                dim=int(force_d.shape[0]),
                inputs=[force_d, self.balance_scale_d],
                device=self.device,
            )

    def add_force_resultant(self, force_d: wp.array) -> None:
        """Reduce an arbitrary contributor's net force into the generic force channel [pN, vec3d]."""
        _validate_device_array(force_d, label="force_d", dtype=wp.vec3d, ndim=1)
        _require_same_device((force_d, self.force_resultant_d), label="force ledger")
        wp.launch(
            _reduce_vec3_kernel,
            dim=int(force_d.shape[0]),
            inputs=[force_d, self.force_resultant_d],
            device=self.device,
        )

    def add_work(self, work_d: wp.array) -> None:
        """Reduce a contributor's work term into the work channel [pN·µm, float64]."""
        _validate_device_array(work_d, label="work_d", dtype=wp.float64, ndim=1)
        _require_same_device((work_d, self.work_d), label="work ledger")
        wp.launch(
            _reduce_scalar_kernel,
            dim=int(work_d.shape[0]),
            inputs=[work_d, self.work_d],
            device=self.device,
        )

    def add_mass(self, mass_d: wp.array) -> None:
        """Reduce a contributor's mass term into the mass-conservation channel [float64]."""
        _validate_device_array(mass_d, label="mass_d", dtype=wp.float64, ndim=1)
        _require_same_device((mass_d, self.mass_d), label="mass ledger")
        wp.launch(
            _reduce_scalar_kernel,
            dim=int(mass_d.shape[0]),
            inputs=[mass_d, self.mass_d],
            device=self.device,
        )

    def add_topology_delta(self, counts_d: wp.array) -> None:
        """Add per-slot topology-count deltas index-wise into the topology channel [int32]."""
        _validate_device_array(counts_d, label="counts_d", dtype=wp.int32, ndim=1)
        if counts_d.shape != self.topology_counts_d.shape:
            raise ValueError("topology-count delta must match the ledger's topology-slot count")
        _require_same_device((counts_d, self.topology_counts_d), label="topology ledger")
        wp.launch(
            _reduce_int_counts_kernel,
            dim=self.n_topology_slots,
            inputs=[counts_d, self.topology_counts_d],
            device=self.device,
        )

    def derive_balance_tolerance(self, gamma_n_d: wp.array, tol_sq_d: wp.array) -> None:
        """Write the squared balance tolerance on the device from this ledger's own resultants.

        Call immediately before :meth:`assemble_balance` and pass it the same ``tol_sq_d``.  Together
        they make the force-balance gate self-scaling: the only input is a COUNT (via ``gamma_n_d``,
        from :func:`assemble_balance_tolerance_ratio`), and the force scale is read from what this
        ledger already accumulated.  Nothing numerical is supplied by a caller who could have chosen it
        to make the gate pass.

        The scale is :attr:`balance_scale_d` (``Σ|x_i|``, the Higham bound's own) when a contributor
        accumulated it, else the resultant bracket ``|reaction| + |traction|``.

        Args:
            gamma_n_d: ``(1,)`` float64 device scalar — the dimensionless float64 accumulation floor
                ``γ_n`` for the run's contribution count.
            tol_sq_d: ``(1,)`` float64 device scalar — written with ``(γ_n · scale)²`` [pN²].

        Raises:
            ValueError: If either array is not a single-entry float64 device scalar, or they are not
                on this ledger's device.
        """
        _validate_device_array(gamma_n_d, label="gamma_n_d", dtype=wp.float64, ndim=1)
        _validate_device_array(tol_sq_d, label="tol_sq_d", dtype=wp.float64, ndim=1)
        if gamma_n_d.shape != (1,) or tol_sq_d.shape != (1,):
            raise ValueError("gamma_n_d and tol_sq_d must be single-entry scalars")
        _require_same_device(
            (gamma_n_d, tol_sq_d, self.reaction_resultant_d, self.traction_resultant_d),
            label="balance tolerance derivation",
        )
        scale_d = self.balance_scale_d
        if scale_d is None:
            scale_d = wp.zeros(1, dtype=wp.float64, device=self.device)
        wp.launch(
            _derive_balance_tolerance_kernel,
            dim=1,
            inputs=[self.reaction_resultant_d, self.traction_resultant_d, scale_d,
                    gamma_n_d, tol_sq_d],
            device=self.device,
        )

    def assemble_balance(self, tol_sq_d: wp.array) -> None:
        """Check ``reaction_resultant + traction_resultant → 0`` against a caller tolerance (device-only).

        Args:
            tol_sq_d: ``(1,)`` float64 device scalar — the squared force-residual tolerance [pN²] the global
                acceptance gate contracts for.  Supplied by the caller so no numerical tolerance is baked in.

        Writes ``|reaction + traction|²`` into :attr:`balance_residual_sq_d` and the pass flag into
        :attr:`balance_ok_d`.  Reads no device value on the host; the predicate consumes ``balance_ok_d`` on
        the device.
        """
        _validate_device_array(tol_sq_d, label="tol_sq_d", dtype=wp.float64, ndim=1)
        if tol_sq_d.shape != (1,):
            raise ValueError("tol_sq_d must be a single-entry squared-tolerance scalar")
        _require_same_device(
            (tol_sq_d, self.reaction_resultant_d, self.balance_ok_d), label="force-balance gate"
        )
        wp.launch(
            _assemble_force_balance_kernel,
            dim=1,
            inputs=[
                self.reaction_resultant_d,
                self.traction_resultant_d,
                tol_sq_d,
                self.balance_residual_sq_d,
                self.balance_ok_d,
            ],
            device=self.device,
        )


def make_global_cell_ledger(*, n_topology_slots: int = 1, device: str = "cuda") -> GlobalCellLedger:
    """Allocate a zeroed :class:`GlobalCellLedger` on ``device`` (CUDA lane convenience constructor).

    Args:
        n_topology_slots: number of int32 topology counters to allocate (>= 1).
        device: CUDA device string; kernels and accumulators live here.  Never a hard-coded device id in
            production — the caller passes the resolved device.

    Returns:
        A ledger with every accumulator zeroed.  Not importable-safe to call on a CUDA-free host (it
        allocates device memory); the CPU structural gates construct the dataclass with array doubles instead.
    """
    if isinstance(n_topology_slots, bool) or not isinstance(n_topology_slots, int) or n_topology_slots < 1:
        raise ValueError("n_topology_slots must be a positive integer")
    return GlobalCellLedger(
        reaction_resultant_d=wp.zeros(1, dtype=wp.vec3d, device=device),
        traction_resultant_d=wp.zeros(1, dtype=wp.vec3d, device=device),
        force_resultant_d=wp.zeros(1, dtype=wp.vec3d, device=device),
        work_d=wp.zeros(1, dtype=wp.float64, device=device),
        mass_d=wp.zeros(1, dtype=wp.float64, device=device),
        topology_counts_d=wp.zeros(n_topology_slots, dtype=wp.int32, device=device),
        balance_residual_sq_d=wp.zeros(1, dtype=wp.float64, device=device),
        balance_ok_d=wp.zeros(1, dtype=wp.int32, device=device),
        balance_scale_d=wp.zeros(1, dtype=wp.float64, device=device),
    )
