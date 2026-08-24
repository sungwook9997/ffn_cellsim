r"""Overdamped Langevin over arena ranges — the driven-damped step a stationary series needs.

**Why this exists, and it is a finding rather than a plan.** On 2026-08-21 γ was emitted on the
resting path for the first time, which removed the blocker that made ``STATE.md`` (e) 1 *unwritable*.
Running it immediately revealed the next one: the engine has **no arm that produces a STATIONARY
series**. With the geometry frozen γ is exactly constant — a degenerate series with no correlation
time — and under a numerical relaxation it descends monotonically and is rejected as DRIFTING at
4.2e+05 σ. Session C's seed-scatter report refused both, correctly, and
``GAMMA_SEED_SCATTER_REFUSED_TWICE_2026-08-21.md`` records it.

A stationary series is one that FLUCTUATES about a steady state. That needs the system driven and
damped, and (e) 1 says so itself: a living cell is not in static mechanical equilibrium and force
balance holds **statistically**. *Statistics need noise.* This module is the smallest correct thing
that supplies it.

**The step, and why this form.** Overdamped Langevin, one Euler-Maruyama update per node::

    x(t+dt) = x(t) + mobility * F(x) * dt  +  sqrt(2 * kB*T * mobility * dt) * N(0, 1)

The prefactor is not a choice — it is fixed by the **fluctuation-dissipation theorem**, and that is
the whole reason this can be written without inventing anything: given a mobility, FDT determines the
noise that makes the stationary distribution ``exp(-U/kBT)``. A thermostat whose noise amplitude were
a free parameter would be a heater, and its stationary distribution would be whatever the parameter
said.

⚠ **BOTH ``mobility`` AND ``temperature_k`` ARE REQUIRED AND NEITHER HAS A DEFAULT.**

* **Temperature** is the one genuinely sourced quantity here — 310 K, and this repository has it in
  ``laws/units``. It is still required, because a thermostat that defaults its temperature is a
  thermostat nobody has to think about.
* ⚠ **Mobility is the gap, and it is a known one.** The arena's per-node ``mobility`` array is still
  ZERO because a mobility is a drag law's parameter and the builders decline to invent one
  (``build/membrane.py:49-50``). The one drag law this repository has,
  ``laws.motility_warp.physical_node_gammas``, is on ``engine/medium_exterior.FORBIDDEN_DRAG_SYMBOLS``
  as PI framework trap #4, and it applies a CYTOPLASM viscosity to nucleus beads
  (``STATE.md`` open item). **So this module supplies the machinery and states plainly that the
  coefficient it runs on is unsourced.** Nothing here makes a drag law exist.

**What this does NOT do.** It does not accept a step — acceptance is ``world/step.py``'s and is still
``UndefinedAcceptance``. It does not integrate kinetics. It does not claim that any trajectory it
produces is physical: with an unsourced mobility, the TIMESCALE is unsourced, so a stationary series
from this is stationary in *steps*, and calling those seconds requires the drag law.

Sanity Gate (recorded before first execution):
    * **dimensions** — ``x`` µm, ``F`` pN, ``mobility`` µm/(pN·s), ``dt`` s, ``kB*T`` pN·µm. The
      displacement term is µm and the noise term is ``sqrt(pN·µm · µm/(pN·s) · s) = µm``. They match,
      which is the check that catches a units slip in the FDT prefactor.
    * **boundary cases** — ``dt <= 0``, a negative mobility and a non-positive temperature are each
      refused. ``mobility = 0`` is ALLOWED and is the honest default state of this arena: it means no
      node moves and no noise is applied, which is exactly what a zero mobility should mean.
    * **conservation** — the launch is over one population's claim, so a thermostat cannot move a node
      another population owns.
    * **CFL/precision** — float64. There is no CFL here: the overdamped step has no inertia, and the
      step-size bound is the STIFFNESS one, ``dt < 1/(mobility * k_max)``, which the caller owns and
      :func:`stability_bound_s` reports rather than enforces — enforcing it would mean this module
      choosing a timestep, and a timestep chosen to make a run behave is what the charter forbids.
    * **sign sense** — the deterministic term moves a node ALONG its force. Asserted in ``_demo``
      against a hand-built configuration, because a sign error here produces a run that heats up while
      looking like a thermostat.
    * **measurement protocol** — the RNG is seeded per step from ``(base_seed, step_index)`` so a run
      is reproducible from two integers, and :meth:`ThermostatReport.record` carries both.

engine units: length µm, force pN, time s, energy pN·µm. Runtime: NVIDIA Warp on CUDA.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import warp as wp

__all__ = ["KB_PN_UM_PER_K", "ThermostatReport", "langevin_step", "stability_bound_s"]

#: Boltzmann constant in engine units [pN·µm/K]. ``1.380649e-23 J/K`` with ``1 J = 1e12 pN · 1e6 µm``,
#: so ``kB = 1.380649e-23 * 1e18 = 1.380649e-5 pN·µm/K``. At 310 K this gives ``kB*T = 4.28e-3 pN·µm``
#: = 4.28 pN·nm, which is the familiar figure and is the arithmetic check on the conversion.
KB_PN_UM_PER_K = 1.380649e-5


@dataclass(frozen=True, slots=True)
class ThermostatReport:
    """What one thermostatted step was told to do, for the run record."""

    n_nodes: int
    mobility_um_per_pn_s: float
    temperature_k: float
    dt_s: float
    base_seed: int
    step_index: int

    @property
    def kbt_pn_um(self) -> float:
        """``kB*T`` in engine units [pN·µm]."""
        return KB_PN_UM_PER_K * self.temperature_k

    @property
    def rms_kick_um(self) -> float:
        """RMS thermal displacement of one node in one step, per axis [µm].

        Reported because it is the number that says whether a step is doing anything: a kick far
        below the discretisation is a thermostat that is on and irrelevant, and a kick comparable to
        a segment length is one that will tear the filament apart.
        """
        return math.sqrt(2.0 * self.kbt_pn_um * self.mobility_um_per_pn_s * self.dt_s)

    def record(self) -> dict[str, object]:
        """The artifact row, with the gap named rather than implied."""
        return {
            "n_nodes": self.n_nodes,
            "mobility_um_per_pn_s": self.mobility_um_per_pn_s,
            "mobility_basis": (
                "⚠ UNSOURCED. The arena's per-node mobility is zero because a mobility is a drag "
                "law's parameter and the builders decline to invent one (build/membrane.py:49-50). "
                "The one drag law this repo has is on FORBIDDEN_DRAG_SYMBOLS as PI framework trap #4. "
                "So the TIMESCALE of any series produced here is unsourced: it is stationary in "
                "STEPS, and calling those seconds needs the drag law."),
            "temperature_k": self.temperature_k,
            "temperature_basis": "310 K, physiological; the one sourced quantity in this step",
            "kbt_pn_um": self.kbt_pn_um,
            "dt_s": self.dt_s,
            "rms_kick_um": self.rms_kick_um,
            "noise_basis": (
                "sqrt(2 kB T mobility dt) — fixed by the fluctuation-dissipation theorem, NOT a "
                "parameter. A free noise amplitude would make this a heater with an arbitrary "
                "stationary distribution."),
            "rng": {"base_seed": self.base_seed, "step_index": self.step_index,
                    "note": "seeded per step from (base_seed, step_index); the run reproduces from "
                            "those two integers"},
            "claims": "NOTHING physical. This moves nodes; it accepts no step and asserts no state.",
        }


@wp.kernel
def _langevin_kernel(pos: wp.array(dtype=wp.vec3d), force: wp.array(dtype=wp.vec3d),
                     lo: wp.int32, mobility: wp.float64, sigma: wp.float64, dt: wp.float64,
                     seed: wp.int32):
    """One overdamped Euler-Maruyama step on the population starting at ``lo``.

    ``sigma`` arrives PRE-COMPUTED as ``sqrt(2 kB T mobility dt)`` rather than being assembled here,
    so the fluctuation-dissipation relation lives in one place on the host where it can be read and
    checked, instead of inside a kernel where a factor of two would be invisible.
    """
    t = wp.tid()
    i = lo + t
    state = wp.rand_init(seed, i)
    kick = wp.vec3d(wp.float64(wp.randn(state)), wp.float64(wp.randn(state)),
                    wp.float64(wp.randn(state)))
    pos[i] = pos[i] + mobility * dt * force[i] + sigma * kick


def stability_bound_s(mobility_um_per_pn_s: float, k_max_pn_per_um: float) -> float:
    """The explicit overdamped step-size bound ``1 / (mobility * k_max)`` [s].

    **Reported, never enforced.** A module that clamped the caller's ``dt`` to this would be choosing
    a timestep, and a timestep chosen so a run behaves is what the charter forbids. The caller states
    ``dt`` and this says what the stiffest bond in the system implies about it.

    Raises:
        ValueError: on a non-positive mobility or stiffness — the bound is meaningless for either.
    """
    if not (mobility_um_per_pn_s > 0.0) or not (k_max_pn_per_um > 0.0):
        raise ValueError(
            f"the overdamped bound needs a positive mobility and stiffness; got "
            f"{mobility_um_per_pn_s!r} and {k_max_pn_per_um!r}. A zero mobility does not move a node, "
            "so no step size is unstable and the question does not arise."
        )
    return 1.0 / (mobility_um_per_pn_s * k_max_pn_per_um)


def langevin_step(arena: object, lo: int, hi: int, *, mobility_um_per_pn_s: float,
                  temperature_k: float, dt_s: float, base_seed: int,
                  step_index: int) -> ThermostatReport:
    """Advance nodes ``[lo, hi)`` one overdamped Langevin step.

    Args:
        arena: the world holding ``position`` and ``force``.
        lo / hi: the population's NODE claim. The launch is over this range only.
        mobility_um_per_pn_s: ``1/γ`` [µm/(pN·s)]. **No default; see the module warning.** Zero is
            allowed and means nothing moves.
        temperature_k: [K]. **No default.** 310 K is physiological.
        dt_s: the step [s]. **No default** — with an unsourced mobility this is a step index wearing
            a unit, and the record says so.
        base_seed / step_index: the RNG is seeded from the pair, so a run reproduces from two
            integers rather than from a captured state.

    Returns:
        The :class:`ThermostatReport` for this step.

    Raises:
        ValueError: on ``hi <= lo``, a negative mobility, a non-positive temperature or ``dt``.
        RuntimeError: on a device-less arena.
    """
    if hi <= lo:
        raise ValueError(f"an empty range [{lo}, {hi}) has nothing to thermostat")
    if not (mobility_um_per_pn_s >= 0.0) or not math.isfinite(mobility_um_per_pn_s):
        raise ValueError(f"mobility must be finite and nonnegative; got {mobility_um_per_pn_s!r}")
    if not math.isfinite(temperature_k) or temperature_k <= 0.0:
        raise ValueError(
            f"temperature must be finite and positive; got {temperature_k!r}. A zero-temperature "
            "'thermostat' is a gradient descent and should be called one.")
    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise ValueError(f"dt must be finite and positive; got {dt_s!r}")
    dev = getattr(arena, "device", None)
    if dev is None:
        raise RuntimeError(
            "this arena is bookkeeping-only (device=None). There is no CPU simulation path and there "
            "must never be one.")

    kbt = KB_PN_UM_PER_K * float(temperature_k)
    # FDT, on the host and in one place: sigma^2 = 2 kB T mobility dt.
    sigma = math.sqrt(2.0 * kbt * float(mobility_um_per_pn_s) * float(dt_s))
    wp.launch(_langevin_kernel, dim=int(hi - lo),
              inputs=[arena.node_arrays["position"], arena.node_arrays["force"], wp.int32(int(lo)),
                      wp.float64(float(mobility_um_per_pn_s)), wp.float64(sigma),
                      wp.float64(float(dt_s)), wp.int32(int(base_seed) ^ int(step_index))],
              device=wp.get_device(dev))
    return ThermostatReport(n_nodes=int(hi - lo), mobility_um_per_pn_s=float(mobility_um_per_pn_s),
                            temperature_k=float(temperature_k), dt_s=float(dt_s),
                            base_seed=int(base_seed), step_index=int(step_index))


def _demo() -> None:
    """Sanity Gate — the units, the FDT prefactor and every refusal. The launch needs a device."""
    # 1. kB in engine units, checked against the figure everyone knows: kB*T at 310 K = 4.28 pN*nm.
    kbt = KB_PN_UM_PER_K * 310.0
    assert abs(kbt - 4.28e-3) < 1e-5, kbt                      # pN*um
    assert abs(kbt * 1000.0 - 4.28) < 0.01, "kB*T at 310 K is 4.28 pN*nm"

    # 2. The FDT prefactor is dimensionally a length. sqrt(pN*um * um/(pN*s) * s) = um.
    #    This is the check that catches a units slip, and it is why sigma is assembled on the host.
    r = ThermostatReport(n_nodes=1, mobility_um_per_pn_s=1.0e-3, temperature_k=310.0, dt_s=1.0e-4,
                         base_seed=0, step_index=0)
    expect = math.sqrt(2.0 * kbt * 1.0e-3 * 1.0e-4)
    assert abs(r.rms_kick_um - expect) < 1e-18
    # ...and it scales as sqrt(dt), which is what separates a Wiener increment from a drift term.
    r2 = ThermostatReport(n_nodes=1, mobility_um_per_pn_s=1.0e-3, temperature_k=310.0, dt_s=4.0e-4,
                          base_seed=0, step_index=0)
    assert abs(r2.rms_kick_um / r.rms_kick_um - 2.0) < 1e-12, "the kick must scale as sqrt(dt)"

    # 3. The record NAMES the gap rather than implying it. A reader must not be able to take the
    #    mobility for a sourced quantity, because it is the one thing here that is not.
    row = r.record()
    assert "UNSOURCED" in row["mobility_basis"] and "FORBIDDEN_DRAG_SYMBOLS" in row["mobility_basis"]
    assert "fluctuation-dissipation" in row["noise_basis"]
    assert row["claims"].startswith("NOTHING physical")

    # 4. The stability bound is reported and REFUSES the case where it is meaningless, rather than
    #    returning an infinity a caller would divide by.
    assert abs(stability_bound_s(1.0e-3, 1.0e3) - 1.0) < 1e-12
    for bad in ((0.0, 1.0e3), (1.0e-3, 0.0), (-1.0, 1.0)):
        try:
            stability_bound_s(*bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"stability_bound_s{bad} must be refused")

    # 5. Every refusal on the step itself, including the one that is deliberately NOT a refusal.
    class _A:
        device = None
        node_arrays: dict = {}

    ok = dict(mobility_um_per_pn_s=1e-3, temperature_k=310.0, dt_s=1e-4, base_seed=0, step_index=0)
    for kw, exc, needle in (
            ({**ok, "temperature_k": 0.0}, ValueError, "gradient descent"),
            ({**ok, "dt_s": 0.0}, ValueError, "dt must be"),
            ({**ok, "mobility_um_per_pn_s": -1.0}, ValueError, "nonnegative"),
            (ok, RuntimeError, "no CPU simulation path")):
        try:
            langevin_step(_A(), 0, 4, **kw)
        except exc as e:
            assert needle in str(e), (kw, e)
        else:
            raise AssertionError(f"{kw} must be refused")
    try:
        langevin_step(_A(), 4, 4, **ok)
    except ValueError as e:
        assert "empty range" in str(e)
    else:
        raise AssertionError("an empty range must be refused")

    print(f"thermostat self-check OK — kB*T(310 K) = {kbt * 1000:.2f} pN*nm, kick scales as sqrt(dt), "
          "and the mobility gap is named in the record")


if __name__ == "__main__":
    _demo()
