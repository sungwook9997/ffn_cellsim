# External GPU options for ffn_cellsim production runs (2026-06-01)

Workload constraints that decide the option: **HOOMD-blue 7.0.1 GPU build
(conda-forge only, not pip-friendly) + cupy + the ffn_sim repo + long unattended
runs (hours→days) + a single GPU (A5000-class is enough).** Turnkey setup on any
box: `aleph/scripts/cloud_bootstrap.sh` (+ `environment-gpu.yml`) → ~10-15 min
to "ready", with a GPU smoke that asserts the GPU build + the per-function
GPU-vs-CPU match before you commit to a multi-hour run.

## Ranked options (for THIS sim)

1. **Persistent GPU VM (SSH + conda) — best fit.** Like gbook in the cloud: run
   for days under tmux/nohup.
   - **RunPod** — community 4090/A5000/A6000 ~$0.2–0.8/hr (4090 fast+cheap),
     network volume for checkpoints. **Best price/perf.**
   - **Lambda Labs** — A6000/A100/H100 on-demand ~$0.5–1.3/hr, clean Ubuntu.
   - **Vast.ai** — cheapest (~$0.2/hr) marketplace; multi-day reliability varies.
   - **Paperspace (DigitalOcean)** — A4000–A100 persistent, mid-price.
   - AWS g5(A10G)/GCP g2,a2/Azure NC — enterprise, pricier + more setup; only if
     you have credits.
2. **KAIST / KISTI HPC (if allocation exists) — likely #1 for an academic.** GPU
   nodes + Slurm + conda, free, built for long jobs. Check before paying.
3. **Modal — containerised job (paid, great reproducibility).** Define an image
   with hoomd, run as a Python function, per-second billing, volume for
   checkpoints. More upfront, best automation.
4. **Colab / Kaggle — NOT for multi-day.** Free T4 slower than A5000; session
   caps (~12-24 h); conda non-native (condacolab); notebook-locked. Short
   interactive checks only.

## When cloud is worth it vs gbook

- A single ~6 h run (e.g. the v0×100 native pilot): **gbook A5000 is enough**
  (free, unattended). No cloud needed.
- Cloud wins for: **(a) parallel config/seed sweeps** (N cheap pods = N× throughput,
  same wall-time), **(b) the literal-v0 gold-standard run (~days)** (faster GPU +
  doesn't tie up gbook).
- Speed-up from a faster GPU is **modest** (~1.5–3×) — at n_fil~3000 the sim is
  partly Python/kernel-launch bound, not pure FLOPs. The real cloud value is
  parallelism + always-on, not single-run speed.
- Cost: ~$0.5/hr GPU → ~$3 for 6 h, ~$30–50 for a multi-day run; 8 seeds in
  parallel = 8× cost, same wall-time. Cheap relative to the science.

## Checkpointing (for runs > a session or on preemptible/community hosts)

- `h3_ku35_tension.py` already has fingerprinted checkpoint/resume — prefer it for
  long/preemptible runs.
- `h3_ku35_gripwalk_tier1.py` (the micro-diagnostic) has no checkpoint; its
  per-sample flushed log preserves the trajectory even if the run dies. On a
  non-preemptible on-demand VM this is usually fine.

## 3D structural visualisation

Production drivers can now dump a **GSD trajectory** (`--gsd-period N` on
`h3_ku35_gripwalk_tier1.py`): HOOMD-native, device-resident positions + box +
topology, no extra per-step sync. Load directly in **OVITO**, or import into
**Blender** (the GSD-writer was the "viz unlock" in the optimization/viz
decisions). 2D metric trajectories (s_grip / adv / γ) are plotted at run
closeout per the visualize-at-closeout rule.
