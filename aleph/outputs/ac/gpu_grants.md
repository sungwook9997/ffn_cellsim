# GPU grant citations

`CLAUDE.md` §Running on the GPU: *"The PI states a card and a duration in session; the session
transcribes it **with transcript path and timestamp**, so anyone can check it. An agent **may** write
that citation and **may never** write a grant."*

This file is the transcription, nothing more. Every row points at a transcript a third party can open.
An entry here is **evidence that the PI said it**, never authority to take a card.

---

## 2026-08-15 14:31 KST — RTX 4090-1, 24 h

| Field | Value |
|---|---|
| Card | `4090-1` (Slurm licence token `rtx4090_1`, physical CUDA index 1) |
| Duration | `24:00:00` |
| Stated by | PI (sw1), in session |
| Transcript | `~/.claude/projects/-Users-sw1-ffn-cellsim/b2dddf6e-7820-4ec9-9de4-95de1d5daafa.jsonl` |
| Timestamp | 2026-08-15 14:31 KST (duration), 14:35 KST (card + hold form, answered in-session) |
| Slurm job | `259` — `MEM_GB=8 CPUS_PER_TASK=8 gpu-submit 4090-1 24:00:00 sleep 86400` |
| Form | **Idle hold.** Work runs inside it via `~/bin/ffn-run` (`srun --jobid --overlap`). |

⚠ **Job `258` was the same grant submitted at `MEM_GB=48` and was cancelled the same hour.** 48 was
copied from the driver docstrings, where it had never been measured. It reserved 48 of the box's 53 GB
`Slurm max RAM` for 24 h, which would have blocked every other user for a day to hold 0.70 GB of
actual need. Cancelled and resubmitted as `259` at `MEM_GB=8`. The docstrings are corrected in the same
change (`ac_resting_residual_curve.py`, `ac_connector_devicerun_census.py`).

The PI stated the duration first and named the card and the hold form when asked; both halves are in the
transcript above. Other active users on the box at submit time: `ham`, `jimin` — neither on `4090-1`,
and the remaining two cards (`4090-0`, `3090`) were left free.

### Verified at grant time (job 257, then inside job 258)

```
Warp 1.14.0 | CUDA Toolkit 12.9, Driver 13.3
"cuda:0" : NVIDIA GeForce RTX 4090 (24 GiB, sm_89, mempool enabled)
DEVICE: cuda:0 | is_cuda: True | device alloc OK
```

### ⚠ `srun --overlap` does not inherit the card pinning

The batch wrap exports `CUDA_VISIBLE_DEVICES=<idx>`, but `srun --jobid=N --overlap` starts a **new
step** and does not inherit it — all three GPUs become visible and `cuda:0` silently resolves to
**4090-0, a card we did not reserve**. The Slurm licence token reserves the slot; it does not isolate
the device. Observed directly here: a bare `--overlap` step listed `cuda:0`/`cuda:1`/`cuda:2`.

`~/bin/ffn-run` is the guard. It **derives** the index from the reservation's own job name
(`rtx4090_1` → 1) rather than taking one as an argument, so it cannot drift from the card actually
held, and it refuses when no hold is running. Use it for every run inside a hold:

```bash
~/bin/ffn-run /home/sungwook/miniforge3/envs/ffn_sim/bin/python <driver> ...
```

Its self-check asserts exactly one visible CUDA device and that the index matches the held card.

---

## 2026-08-20 14:45 KST — 4090-1, 12 hours — session `d380041d`

**TRANSCRIBED, NOT GRANTED.** `CLAUDE.md`: *"An agent **may** write that citation and **may never**
write a grant."* This is a citation of what the PI said in session; the authority is theirs.

| field | value |
|---|---|
| card | **4090-1** |
| duration | **12 hours** |
| PI wording, verbatim | `4090_1 12시간 승인` |
| stated at | **2026-08-20 14:45 KST** |
| session | `d380041d-06a6-49ea-9daf-1e4e57e76bf0` |
| transcript | `~/.claude/projects/-Users-sw1-ffn-cellsim/d380041d-06a6-49ea-9daf-1e4e57e76bf0.jsonl` |

**Other people's work on the box at the time of transcription** — `squeue` read at 13:4x KST:
job **59**, user **`ham`**, name `rtx3090`, RUNNING with ~8 h 39 m left of a 12 h limit. **The 3090 is
theirs and is not touched.** No other job was queued or running, so 4090-0 and 4090-1 were unclaimed —
⚠ inferred from the naming convention and an empty queue, **not measured**: `Gres=(null)` means Slurm
does not track the GPUs as a resource, and `nvidia-smi` outside an allocation returns *"GPU access
blocked by the operating system"*.

**What this hold is intended for:** the L0 question — does a force-balanced resting configuration exist
at all — under the resting bound-myosin setpoint declared as a TEST AXIS, with the negative control
declared in advance: with the tension source OFF the residual must NOT descend.

Run inside the hold with `~/bin/ffn-run`, never a bare `srun --overlap` — see the warning above.

---

## 2026-08-22 16:46 KST — 4090-1, 12 hours — session `9f70b354`

**TRANSCRIBED, NOT GRANTED.** `CLAUDE.md`: *"An agent **may** write that citation and **may never**
write a grant."* This is a citation of what the PI said in session; the authority is theirs.

| field | value |
|---|---|
| card | **4090-1** (`rtx4090_1`, `CUDA_VISIBLE_DEVICES=1`) |
| duration | **12 hours** |
| PI wording, verbatim | `일단 gpu 남는거 있으면 4090 하나 잡아서 12시간 승인` |
| stated at | **2026-08-22 16:46 KST** |
| session | `9f70b354-71ef-4f62-b406-c2ae5915408d` |
| transcript | `~/.claude/projects/-Users-sw1-ffn-cellsim/9f70b354-71ef-4f62-b406-c2ae5915408d.jsonl` |
| Slurm | job **74**, `StartTime 2026-08-22T16:47:47`, `EndTime 2026-08-23T04:47:47` |

**⚠ The grant is conditional and the condition was checked before submitting.** The PI said *"if there
is a spare GPU"*. `squeue` at 16:47:11 KST was **EMPTY** — job 73 (yesterday's 4090-1 hold) had ended
at its own `EndTime 2026-08-21T17:00:46`, and `sinfo` showed the node **idle** on both partitions. No
other user's job was displaced, and this time that is measured rather than inferred: `gpu-submit`'s own
pre-flight printed `Selected GPU active compute processes: none` before `[APPROVED]`.

**What this hold is intended for**, from the four rulings the PI issued in the same message:

* **queue 18 → (a)** — widen the candidate range in BOTH `stationarity.py` copies and re-analyse.
  ⚠ **Zero GPU**: the records are on disk. Listed here only so the hold is not read as covering it.
* **queue 14 → (a)** — fix the excursion guard FIRST, then rule the radius. The re-measure after the
  fix is TIER 0, seconds of card time.
* **queue 3 → (a)** — sphere-in-contact, *with the standing condition the PI attached*: a spread
  archetype must take **(b)**, and (b) is a different builder plus all four call sites.
* **queue 11 → (c)** — resolve the `capture_radius` contradiction (0.05 vs 0.210 µm) before anything
  else in that lane. Diagnosis, not a run.

⚠ **So most of tonight's work needs no card.** The hold exists so the re-measures after 14 and 11 do
not wait for one, and so a long native run can follow if the rulings open one. Run inside it with
`~/bin/ffn-run`, never a bare `srun --overlap`.
