# GPU citation — AFM autonomous continuation, 12 hours, 2026-08-10

> **This is a citation, not a grant authored by the agent.** It transcribes the PI's two adjacent
> instructions in the active Codex task and does not broaden them.

| Field | Value |
|---|---|
| Duration | **12 hours**, 2026-08-10 00:59 KST → expires **2026-08-10 12:59 KST** |
| Host | `sungwook@100.110.26.26` — shared WSL2 + Slurm workstation |
| Default route | Use an available approved 4090 through Slurm |
| Conditional route | If Claude is using a 4090, inspect occupancy and use the **idle 3090** |
| PI message 1 | *"어떻게든 완성시키고, 클로드가 4090 쓴다고 하면은 너는 빈 상황 보고 3090 써"* |
| PI message 2 | *"12시간 동안 자율 진행 나 자러 갈테니깐 화이팅"* |
| Timestamp basis | Host clock read immediately after the second message: `2026-08-10T00:59:24+09:00` |
| Transcript | Active Codex task; filesystem transcript path / record UUID is not exposed and is not invented |

## Operational scope

- Every GPU process is submitted with `gpu-submit`; never run CUDA outside a Slurm allocation.
- Inspect `squeue` and the launcher's occupancy report before each submission.
- Do not take a card occupied by another user. The 3090 route is conditional on the PI's stated 4090
  contention case and on the 3090 being idle.
- This authorization covers continued C-2 diagnosis/fix, native validation, and the AFM mechanism-demo
  workflow already in scope. It does not authorize changing PI-gated physical parameters, gate thresholds,
  or the frozen incumbent biological contract.
- At 12:59 KST the citation lapses; stop new submissions rather than extending it.
