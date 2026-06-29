"""Generate + commit the HOOMD reference fixture for the B1 BAOAB Warp parity test.

Runs the FROZEN ``integrator.baoab.LeimkuhlerMatthewsBAOAB`` Action inside a real
``hoomd.Simulation`` (CPU) driven by a position-independent ``md.force.Constant``,
for K steps at a fixed seed, on a TRICLINIC box (non-zero xy/xz/yz) so the wrap's
shear handling is exercised. The resulting trajectory endpoint is the committed
bit-parity TARGET for ``tests/warp_port/test_baoab_warp_parity.py`` — the Warp
kernel is graded against THIS, never against a fresh self-authored oracle.

Why ``md.force.Constant``: the force is position-independent, so ``net_force`` is
the same constant every step regardless of force-eval ordering. That isolates the
integrator step (the B1 piece) from any force port, and lets the Warp run consume
the identical constant force.

The per-step noise the Action draws (``np.random.default_rng(seed).standard_normal``)
is reproduced deterministically here and committed too — the ParticleSorter is
disabled so snapshot rows stay tag-ordered (row i == tag i), making the draw order
== tag order and thus exactly reproducible. At kT=0 the noise is multiplied by 0
(bit-for-bit parity); at kT>0 the SAME noise is injected into the Warp kernel so
only float-op ordering can differ.

Run (commits two .npz fixtures next to this file)::

    python ffn_sim/warp_port/fixtures/generate_baoab_fixture.py
"""

from __future__ import annotations

import os

import numpy as np

import hoomd

from ffn_sim.integrator.baoab import make_baoab_updater

HERE = os.path.dirname(os.path.abspath(__file__))

# --- fixed problem definition (committed alongside the fixture) ---------------
SEED = 12345
DT = 0.005
K = 400
TYPES = ["A", "B", "C", "D"]
N_PER_TYPE = 16
N = len(TYPES) * N_PER_TYPE  # 64
# Distinct constant force per type (position-independent) and distinct drag, so
# the per-type γ / bd_prefactor logic is genuinely exercised.
FORCE_BY_TYPE = {
    "A": (0.70, -0.30, 0.50),
    "B": (-0.40, 0.60, -0.20),
    "C": (0.25, 0.15, -0.65),
    "D": (-0.55, -0.45, 0.35),
}
GAMMA_BY_TYPE = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.5}
# Triclinic box small enough that the constant-force drift crosses image
# boundaries during the run (exercises image-flag accumulation + shear wrap).
BOX = (1.5, 1.5, 1.5, 0.10, 0.05, 0.07)  # Lx, Ly, Lz, xy, xz, yz


def _build_initial_snapshot() -> "hoomd.Snapshot":
    snap = hoomd.Snapshot()
    snap.particles.N = N
    snap.particles.types = TYPES
    Lx, Ly, Lz, xy, xz, yz = BOX
    snap.configuration.box = [Lx, Ly, Lz, xy, xz, yz]

    rng = np.random.default_rng(987)  # init-position RNG (independent of the BAOAB seed)
    # Spread positions across the inner box so the drift wraps but nothing starts
    # outside the box. Tag order == construction order (row i -> tag i).
    pos = np.zeros((N, 3), dtype=np.float64)
    typeid = np.zeros(N, dtype=np.uint32)
    row = 0
    for tid, tname in enumerate(TYPES):
        for _ in range(N_PER_TYPE):
            pos[row] = rng.uniform(-0.4, 0.4, size=3)
            typeid[row] = tid
            row += 1
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = typeid
    snap.particles.image[:] = np.zeros((N, 3), dtype=np.int32)
    return snap


def _reproduce_noise() -> np.ndarray:
    """The exact per-step W_n sequence the Action draws (tag-ordered, K×N×3).

    The Action draws ``self._rng.standard_normal((N, 3))`` once per ``act()`` from
    ``np.random.default_rng(SEED)``; with the sorter disabled, snapshot rows stay
    tag-ordered, so this reproduces the draws bit-for-bit.
    """
    rng = np.random.default_rng(SEED)
    return np.stack([rng.standard_normal((N, 3)) for _ in range(K)]).astype(np.float64)


def _run_reference(kT: float) -> dict[str, np.ndarray]:
    device = hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=device, seed=1)
    sim.create_state_from_snapshot(_build_initial_snapshot())

    # Disable the ParticleSorter so local-snapshot rows stay tag-ordered — makes
    # the Action's row-ordered noise draw == tag order and thus reproducible.
    sim.operations.tuners.clear()

    integrator = hoomd.md.Integrator(dt=DT)
    const = hoomd.md.force.Constant(filter=hoomd.filter.All())
    for tname, fvec in FORCE_BY_TYPE.items():
        const.constant_force[tname] = fvec
    integrator.forces.append(const)
    integrator.methods = []  # forces-only; the L-M Action does the position step
    sim.operations.integrator = integrator

    action, updater = make_baoab_updater(
        kT=kT, gamma=GAMMA_BY_TYPE, dt=DT, seed=SEED
    )
    sim.operations.updaters.append(updater)

    sim.run(0)   # seed net_force before the first L-M step
    sim.run(K)   # K L-M steps
    assert action.steps_run == K, (action.steps_run, K)

    # The GLOBAL snapshot from get_snapshot() is canonically tag-ordered
    # (row i == tag i), independent of the (disabled) local-row sorter.
    snap = sim.state.get_snapshot()
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    image = np.asarray(snap.particles.image, dtype=np.int32)
    prv = np.asarray(action.prv_rnds, dtype=np.float64)  # already tag-indexed

    # Verify the constant force HOOMD applied matches our per-type table (so the
    # committed `force` array the Warp run consumes is exactly the reference's).
    typeid_sorted = np.asarray(snap.particles.typeid)
    force = np.array(
        [FORCE_BY_TYPE[TYPES[t]] for t in typeid_sorted], dtype=np.float64
    )
    gamma = np.array(
        [GAMMA_BY_TYPE[TYPES[t]] for t in typeid_sorted], dtype=np.float64
    )
    return {
        "pos": pos,
        "image": image,
        "prv": prv,
        "force": force,
        "gamma": gamma,
        "typeid": typeid_sorted.astype(np.int32),
    }


def main() -> None:
    init_snap = _build_initial_snapshot()
    tag0 = np.arange(N)  # construction order
    pos0 = np.asarray(init_snap.particles.position, dtype=np.float64)
    image0 = np.asarray(init_snap.particles.image, dtype=np.int32)
    noise = _reproduce_noise()

    for label, kT in [("kt0", 0.0), ("ktpos", 1.5)]:
        ref = _run_reference(kT)
        out = os.path.join(HERE, f"baoab_ref_{label}.npz")
        np.savez(
            out,
            # problem definition
            seed=SEED, dt=DT, K=K, kT=kT, box=np.array(BOX, dtype=np.float64),
            pos0=pos0, image0=image0, force=ref["force"], gamma=ref["gamma"],
            noise=noise, tag0=tag0, typeid=ref["typeid"],
            # reference endpoint (the parity TARGET)
            ref_pos=ref["pos"], ref_image=ref["image"], ref_prv=ref["prv"],
        )
        drift = np.linalg.norm(ref["pos"] - pos0, axis=1)
        print(
            f"[{label}] kT={kT} wrote {out}  "
            f"max|drift|={drift.max():.4f}  "
            f"max|image|={np.abs(ref['image']).max()}"
        )


if __name__ == "__main__":
    main()
