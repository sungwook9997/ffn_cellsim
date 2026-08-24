# ffn_cellsim — Claude Code Project Context

> **This file holds only what does NOT change.** Every count, measurement, inventory, cell line, hardware
> model and status belongs in `STATE.md` and **must be read there, never quoted from here.** A number that
> appears in this file is a bug: it will go stale, it will be injected into every boot, and it will be
> believed. Sections below that unavoidably touch changing things are marked **⚠ VERIFY** — treat those
> as pointers to go look, not as facts.

## What this is

**Project Aleph** is a classical probabilistic mechanistic virtual cell. The product is not a
trajectory and not a table of fitted parameters — it is an **evidence-bounded posterior** over latent
cell state, mechanism and parameters, **with an honest account of what the observations cannot
determine.**

The motivating problem: every method for measuring cell mechanics approximates away what its designer
judged ignorable, and the same cells give order-of-magnitude different answers across methods. This
makes those viewpoints executable in one sandbox and locates what each cannot identify.

The instrument is a fine-grained mechanistic **NVIDIA Warp/CUDA** simulation where every cytoskeletal
filament, motor head, adhesion clutch and ECM cross-link is explicit. **It is not the product** — it
is one of three layers, beside a network learned from its own accepted runs (`inner/`) and one
learned from external research (`outer/`).

This supersedes *"infers molecular parameters per cell type and state"*, which it contains: a
parameter table is a posterior with its uncertainty and unidentifiability deleted. PI-ratified
2026-08-09 by the rename (A1).

## Architectural principle (PI 2026-05-19)

**Mechanistic over lumped, at every design decision.** A closed form (Bell-Evans, Hill, Pereverzev) is an
acceptance ORACLE, never the runtime mechanism. The v1 codebase used published models AS the runtime;
this one inverts that.

**This is not a cost tradeoff.** Measured at native, filament physics is a small minority of a step and
scaling the filament count barely moves it — the cost is the fluid grid and the moving membrane boundary.
So when a budget is tight, fidelity is NOT what you trade away: the step count and the solver are.
Coarsening the cell buys almost nothing and forfeits the population the conclusion needs. (⚠ VERIFY the
current numbers in `STATE.md`. If they ever stop supporting this, the principle needs re-arguing — not
quiet coarsening.)

**The one worked example still open — compartment coupling.** Co-location in a shared array is NEVER a
connection. Components own state and couple only through explicit bidirectional + adjoint connectors
(Newton's 3rd law), and a kinetic connector commits ONLY on an accepted physical step under one
device-resident transaction. Every other example from the original table is IMPLEMENTED and
closed — it lives in the code and in the tests, not in this file.

### Hard runtime contract (PI 2026-07-16)

- **Warp CUDA GPU is the only simulation runtime.** HOOMD is never imported or executed — not for
  production, development, parity, fallback, benchmarks or validation. Python may configure, orchestrate
  and postprocess; every physics kernel, PDE/KMC update, neighbour query and mechanical iteration is
  GPU-resident Warp with no authoritative per-step host state. A static contract test enforces this.
- **Full native population, always** — derived from physiological density × geometry, never a count typed
  into this file. Cortex is only one compartment; lamellipodium, SF/arcs, filopodia, MT, IF, nuclear
  lamina/chromatin, ECM, motor heads, crosslinkers and clutches each add from their own density × geometry.
  **Never lower a biological density to fit memory** — optimise device layout or use larger/multi-GPU
  hardware. Track unique active filament IDs, dormant allocation, node/head/state totals and exact peak
  GPU bytes.
- **One filament, one component.** Cortex / SF-arc / lamellipodium / filopodium are separate state-owning
  components owning DISJOINT populations — not one label-blind network. They couple only through
  connectors, never shared nodes and never a permanent weld.
- **Develop on a slice, conclude at native.** A conclusion drawn below native is not a weak result; it is
  not a result. **Structural tests passing ≠ production.**
- Components climb an 8-state evidence ladder (`CONTRACTED → … → PRODUCTION`, skipping forbidden). Current
  rungs and the component/connector census: **⚠ VERIFY in `STATE.md` (a)** — they change every session.

When in doubt, write down the abstraction the option introduces and ask: *does this replace a mechanistic
process with a lumped one?* If yes, prefer the mechanistic alternative.

## Stack

- **Python 3.13.13**; `conda activate ffn_sim`. `python` is NOT on PATH on either machine — always invoke
  the env interpreter explicitly.
- **Simulation core: NVIDIA Warp on CUDA. There is no CPU simulation path and there must never be one.**
  CPU execution is FORBIDDEN — not discouraged — for production, development, debugging, parity, fallback,
  smoke tests and benchmarks alike. **A number measured on CPU is not a weak result; it is not a result**,
  and it may not be quoted, compared, or used to close a gate. The answer to memory pressure is larger or
  multi-GPU hardware, **NEVER a smaller cell**. No hard-coded device IDs. A driver that resolves a
  non-CUDA device must RAISE, never proceed. The dev machine has no CUDA: code is authored there, and
  every physics claim is measured on the GPU host. (Which device, and what is minimum: **⚠ VERIFY in
  `STATE.md`.**)
- **GPU-only inside the loop.** No per-step host snapshot or host state mutation, no hidden fallback. A
  profiler gate must show zero authoritative GPU→CPU roundtrips inside the physical-time loop; host
  readback belongs BETWEEN accepted steps — and is not free there either.
- **Time architecture**: one outer physical-time loop (Biot/RAD/KMC/water flux/loads) + a Warp-GPU inner
  mechanical solve.
- **Active filament model**: Cytosim/NF2007 mechanics + head-resolved Stam-Hocky NMII with Hill
  force-velocity, as Warp kernels.
- **Engine**: `aleph/engine/` binds `aleph/laws/`; parts live in `aleph/components/`. Which is
  canonical and which is the frozen incumbent: **`STATE.md` (a) — single source.** The HOOMD tree was
  DELETED 2026-07-29; nothing may import it.

## Running on the GPU (read before any native run)

**The device is the shared workstation `sungwook@100.110.26.26`** — WSL2 + Slurm, two RTX 4090 and one
RTX 3090, reached over Tailscale. **`gbook` and its A5000 no longer exist** (PI retired it 2026-08-04);
measurements taken on it stay valid evidence, labelled *(RTX A5000, retired)*, and it is not a target.

**Never run GPU work outside a Slurm allocation.** Submit with
`MEM_GB=<GB> CPUS_PER_TASK=<n> gpu-submit <card> <time> <command>`. Outside an allocation `cuInit()`
returns `CUDA_ERROR_NO_DEVICE`, so this is enforced by the operating system rather than by us — which
is why `ffn_gpu.py`'s lease was archived rather than ported (decision B2): a second lock bolted to
the outside of a door the kernel already holds shut, expiring on wall-clock rather than on work.

**Before quoting a native number, check the host held the code the record stamps** —
`aleph/scripts/run_provenance.py`. A record whose commit does not describe what ran cannot be
corrected afterwards: the run is over and the number is in a file.

**The machine is shared with other people.** Never take a card that already has someone else's job on
it, and **ask the PI before anything heavy, every time**. The PI states a card and a duration in
session; the session transcribes it *with transcript path and timestamp*, so anyone can check it. An
agent **may** write that citation and **may never** write a grant.

## Repo layout

Map: `STRUCTURE.md`, rewritten 2026-08-09. **Do not keep a copy here** — the copy that preceded it
had drifted on all four of its entries.

## Knowledge base

Literature and decision knowledge lives in a **Notion Contract-Graph — the single source of truth.** Two
read layers are materialised FROM it: an Obsidian vault (`outputs/obsidian_rag_full/`) and a DuckDB +
TAG query engine (`outputs/tag_kb/`). **Never hand-edit the vault or the `.duckdb` — fix Notion and refresh.**

- structure, "what links to X" -> Obsidian. relational / aggregate / multi-hop -> `tag_query.py`.
  a passage from a paper -> the BM25 content layer.
- **Before citing a KB source in a deliverable, confirm its `source_audit.verdict` is OK.**
- **New `ValidationGate` / `ModelContract` rows are PI-authored, never auto-created.** The harvester only
  FLAGS unmatched references; creating one is a gate-contract change and needs PI sign-off.
- Code->contract edges come ONLY from an explicit `Implements: KU-x.y` docstring line. A bare KU mention
  creates no edge — before this was tightened it linked two modules to the wrong contract.
- The harvest is out-of-band: skip it and the graph goes stale, your run is unaffected.

**⚠ VERIFY — everything countable about the KB changes and must be read from the database, never from a
document.** Row counts, edge counts, PDF and chunk totals, vault size, audit verdict tallies: run
`bash outputs/tag_kb/refresh.sh` (non-destructive; reports drift) or query `kb.duckdb` directly. The
version of this section that stood until 2026-07-28 hard-coded six such numbers and **every one of them
was wrong** — every one, and `STATE.md` in the same boot context said so. Operational detail for
this subsystem lives in `outputs/tag_kb/README.md`.

> TAG = **Table-Augmented Generation** (Biswal et al. 2024), not "태그/label".

## Code conventions

- Type hints, Google-style docstrings. `pytest` collects `aleph/tests/{ac,coordination}/`.
- **Constants are Python literals, not YAML.** The `params_i0b*.yaml` ledgers and
  `validation/oracles/configs/` are read by tests and KB gates only — **editing one changes nothing at
  runtime.** Change the constant where it is defined.

## Hard rules

- **Nothing is chosen to make a gate pass** — not a tuning constant, not a threshold. A constant must be
  derivable and grid-invariant; a gate is a contract written BEFORE the run and may not be edited inline
  or re-thresholded after seeing the data. If either cannot hold, halt and surface to PI.
- **Start at the physiological operating point.** Every parameter and boundary condition begins at its
  real in-vivo value, never a convenient null — a model started from a non-physical baseline cannot be
  compared to real data, because the baseline is what is wrong. `default-off` is code hygiene, not a
  licence to run physics from an unphysical zero. If a physiological value is unknown, surface to PI.
- **Conclude only at the full native population and compartments.** Coarse configurations are for
  development and may not support a physics conclusion.
- **Sanity Gate before the first execution of any physics module**: dimensions, boundary cases,
  conservation invariants, CFL/precision, sign sense, measurement protocol. Record it in the docstring.
- **Every run writes its own figures** beside its record, committed with the data; a component reaching
  `CUDA_UNIT`/`CONNECTED` also renders itself + its connectors in isolation and refreshes the cumulative
  render — full native resolution, **verified in a real browser**, not by grepping the HTML. No axis
  truncation, no log/linear flip without a note, reference band overlaid, units annotated.
- **Never claim more than the evidence class supports.** Each of these is a way a result gets promoted
  past what produced it, and every one has happened here:
  - a **force-accepted** trajectory is not an accepted physical trajectory — it is what the solver
    reached in the iterations it was given, and must be reported as that;
  - a **surrogate / reduced / tensor-compressed** result may not be recorded as a reference-physics
    result, and its approximation error belongs in the artifact, not in a footnote;
  - a gate scores a **scalar from a projection declared before the run**, taken from the stored field —
    never a scalar the run chose to keep;
  - a parameter invisible in the observation may not be reported as a point estimate;
  - a law established on one cell archetype is not a universal law until it is shown to transfer;
  - a model may not be scored against answers its own machinery produced;
  - nothing is called **quantum** unless quantum dynamics were computed — `quantum-inspired` is the
    ceiling for tensor-network and path-integral methods applied to classical dynamics.
- **Ask the PI before deviating from any rule here.**

## Roadmap

The destination is §What this is. The staged path to it — what each stage must PROVE, and the arithmetic
showing a brute-force parameter sweep is impossible by three to five orders of magnitude — is
**`ROADMAP.md`**, deliberately outside the boot context.

**⚠ VERIFY — everything in `ROADMAP.md` changes, so read it rather than recalling it.** Which stage is
current, which gates have passed, the cell types and states, the stage-0 observable set and every number
in the budget arithmetic are all status; **the gate wording itself is under an amendment that is PROPOSED
and NOT ratified.** `ROADMAP.md` states what a stage must prove and never whether it has been proved —
that is `STATE.md` (b), which wins if the two disagree.

**What is NOT a gate, ever** — each of these has passed while the physics underneath was wrong: a run
completing; a green test suite; a conclusion drawn on a slice; a magnitude quoted without the axis that
says whether it may be quoted at all.

## Working model

One **Lead session** drives the work end to end — writes code, runs production on the GPU host,
maintains Notion at closeout — and spawns subagents on demand (`Explore`/`Plan` for heavy search, a
review pass at closeout). **The PI is an approver, not a relay**, and signs off on gate-contract changes,
magic-number triggers, DONE ratification and `ffn/foundation` pushes.

### Session boot

1. Read this file, then `STATE.md`.
2. `git log --oneline -5`; confirm branch and latest commit.
3. `conda activate ffn_sim`; confirm the device is CUDA Warp
   (`$CONDA_PREFIX/bin/python -c "import warp as wp; wp.init(); print(wp.get_device())"`).
4. Restate the next task in one or two sentences. Ask the PI if it is ambiguous.

### Session closeout — MUST

Steps 1-4 complete before the final user-facing message; step 5 is the PI's receipt.

1. Commit. Do not push to `ffn/foundation` without PI sign-off. **Use a pathspec**
   (`git commit -o <paths>`): the tree is shared and a bare `git commit` sweeps another session's stage.
2. Update `STATE.md` — mandatory whenever a result lands, a claim is retired, or a rung changes.
3. Update the Notion **Development Logs** status board + Open items, and append a milestone Day-log
   (start/end commit hashes, gate PASS/FAIL, KU cross-reference).
4. Run `make kb-check`. On **DRIFT** — a result or constant claimed better than the disk and the
   citation audit support — **halt and surface to PI**. Never close over a red gate. If a store cannot
   be written, halt and surface rather than silently skipping.
5. **Receipt**: the final message ends with the literal line **`Notion 업데이트 완료`**. Then stop; do
   not speculate beyond what was done.

### File ownership

**Parallel sessions share ONE tree; isolation is a declaration, not a second checkout.** Each claims path
globs in `ownership.yaml`, enforced by `.githooks/pre-commit` against `$FFN_SESSION` — a commit outside the
claim is refused. Worktree-per-worker was abandoned: it moved the integration problem instead of removing
it. `aleph/validation/oracles/`, `aleph/docs/briefs/`, `STRUCTURE.md`, `README.md` and `pyproject.toml`
are read-only during normal work.

**Carve-out — `STATE.md` is EXPLICITLY EDITABLE by any session at every milestone** (`free` in
`ownership.yaml`), with no PI pre-approval, and MUST be re-dated in the same commit: a stale `STATE.md` is
a defect, not a formality. Its tier-(a) block is GENERATED — edit `docs/state_rows.yaml` and run
`make state`; the hook refuses a drifted block.
