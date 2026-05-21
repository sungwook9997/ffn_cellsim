"""Parallel ensemble runner — utility module for running independent
seed / parameter sweep tasks across M1 Max's multiple cores.

Each HOOMD CPU simulation is single-threaded, but **independent**
simulations (different seeds, different integrators, different dt) can
be run in parallel via ``multiprocessing.Pool``.  This module provides
the boilerplate.

Why we need this
----------------
``feedback_acs_no_abstractions`` says "fine-grained mechanistic"
options.  Ensemble averaging is the mechanistic way to estimate
ensemble means for Mikado random-network observables (Storm-MacKintosh
convention).  Sequential ensembles waste 5×+ of available wall time on
M1 Max (10 performance cores).

Constraints
-----------
- Each subprocess imports HOOMD = ~1-2 GB memory + ~5 s startup.
  Don't use Pool for trivially-small tasks; reserve for ≥ 30 s tasks.
- For N=65 982 Mikado runs, each subprocess uses ~3-5 GB.  Limit
  ``n_workers`` to ``min(n_tasks, 8)`` on a 32 GB M1 Max.
- Pool uses ``fork`` start method on macOS; HOOMD's RNG seed handling
  is per-process so seeds must be passed explicitly (not relying on
  process-shared state).

Usage
-----
::

    from ffn_sim.scripts.parallel_ensemble import run_seed_pool

    def task(seed):
        # Build sim with this seed, run, return result dict
        ...
        return {"seed": seed, "L_p": 17e-6, ...}

    results = run_seed_pool(task, seeds=list(range(20)), n_workers=8)
"""

from __future__ import annotations

import multiprocessing as mp
import os
import time
from typing import Any, Callable, Sequence


def run_seed_pool(
    task_fn: Callable[..., Any],
    seeds: Sequence[int],
    n_workers: int | None = None,
    *,
    extra_args: tuple = (),
    extra_kwargs: dict | None = None,
) -> list[Any]:
    """Run ``task_fn(seed, *extra_args, **extra_kwargs)`` for each seed.

    Parameters
    ----------
    task_fn
        Top-level (importable, picklable) function to run per seed.
        Must accept ``seed`` as first positional arg.
    seeds
        Iterable of integer seeds.
    n_workers
        Number of worker processes.  Defaults to
        ``min(len(seeds), os.cpu_count() - 2, 8)`` — leave a couple
        cores free for the parent process + OS.
    extra_args, extra_kwargs
        Additional arguments forwarded to ``task_fn``.

    Returns
    -------
    List of task_fn return values in the same order as ``seeds``.
    """
    if extra_kwargs is None:
        extra_kwargs = {}
    if n_workers is None:
        n_workers = min(len(seeds), max((os.cpu_count() or 4) - 2, 1), 8)
    if n_workers < 1:
        n_workers = 1

    t0 = time.time()

    # Build (seed, *extra_args) tuples for starmap; kwargs handled via
    # functools.partial since starmap doesn't take kwargs natively.
    from functools import partial
    bound = partial(task_fn, **extra_kwargs)
    arg_tuples = [(s,) + tuple(extra_args) for s in seeds]

    # macOS default ``spawn`` is safer than ``fork`` for HOOMD-heavy
    # subprocesses (avoid fork-time fd leakage from HOOMD's CUDA-less
    # device init).
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=n_workers) as pool:
        results = pool.starmap(bound, arg_tuples)

    elapsed = time.time() - t0
    print(
        f"[parallel_ensemble] {len(seeds)} tasks × {n_workers} workers "
        f"completed in {elapsed:.1f} s "
        f"(speedup ≈ {len(seeds) * (elapsed / len(seeds)) / elapsed:.1f}×)"
    )
    return results
