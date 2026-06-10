# 다음 세션 부트 프롬프트 — ffn_cellsim / h7/full-cell-integration (2026-06-11 작성)

부트 순서: CLAUDE.md → 이 파일 → Dev Logs 상태판 → `git log -8` → `conda activate ffn_sim`.
env: 비대화형/백그라운드는 `PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python` 직접
([[reference_ffn_run_invocation]]). 현재 HEAD `9cf4091` (Mac), gbook 동기화됨.

## 이 라인의 정체
fine-grained 단일세포 magnitude 라인 (route a). cortical active-γ floor + ventral SF
traction의 절대값을 책임짐. active-γ는 **generation/engagement-limited**으로 수렴
(loop24b: per-SF coherent traction가 seed-불안정 → "decisive +131 pN"은 비-robust;
transmission/aggregation/coherence는 벽 아님). band 도달엔 density datum(MCF7) 또는
측정 프로토콜 수정 필요.

## 2026-06-11 세션이 한 것 (커밋)
- `7a84b60` **loop24d M-band targeting A/B-REJECTED** — 개선 없음(paired Δ +62±156 pN, NS),
  per-SF traction seed-불안정 재확인. 코드 보존·default-off(`--mband-mode off`).
- `82fd7e8` 고아 magnitude figure 12개 커밋(visualize-at-closeout).
- `92b4029` **gbook 동기화 복구** — gbook→Mac SSH 키 셋업(미설정이었음), gbook을 dirty
  h3-cortex → clean h7 `@`로 정리(백업 `~/ffn_cellsim_dirtybackup_2026-06-11.tar`).
  gbook: A5000 16GB, hoomd gpu_enabled=True, cupy 14.1.0. [[project_ffn_production_infra]]
- `c41a65a` H.2 driver `--device`(GPU-main) + P2c sign-off PROCEDURE.
- `070e038`/`e8d6ce9` **P2c GPU gate 재검증 sign-off — evidence 완료**. Tier0 integrator
  parity PASS, Tier1 H.3 L_p GPU=CPU=17.46µm in-band PASS, Tier2 H.2 GPU≡CPU(round-off)
  PASS-by-transfer. `outputs/gpu_signoff/REPORT.md`. **→ PI RATIFY 대기** (GPU를 H.2/H.3
  L_p production device로 비준; gate-contract = PI 고유권한).
- `c852e5b` Taeyoon-Kim 그룹 새 ref 3편 인게스트 + ab_xlink 재설계 브리프 + Slater oracle 초안.
- `9231b17`/`9cf4091` **ab_xlink_v2 예비** — off-rate 축 VALID(g_actin +4.5%/총 γ +0.7%,
  10× 느림), **density 축 미검증**(connected_mesh에서 n_xl=seeding output라 `--xl-n` 무시;
  cell.py:745). magnitude 판정 불변(engagement-limited). `outputs/h7/production/ab_xlink_v2/REPORT.md`.

## 자율 가능 다음 스텝 (PI 결정 불필요)
1. **ab_xlink_v2 density 손잡이 수정 후 재시도**: connected_mesh에서 density는 `--cm-z-struct`
   (또는 connected_mesh=False + dynamic n_xl)로 변해야 함 — 그 축으로 transmission sweep 재실행.
2. **Slater stress(r) 추정기 구현**: 반경 r 면 가로지르는 chain의 axial force/면적 — coherent-FA
   -traction과 독립적 2차 추정기. cortical_tension.py method-of-planes를 원통 shell로 적응
   (brief AB_XLINK §3.2 / cortical_tension.py:256-394).

## HALT → PI (결정 대기)
- **GPU sign-off ratification** (P2c) — 위 070e038 evidence 검토 후 비준.
- **결정2 density datum**: MCF7 cortical myosin 밀도 1차 실험데이터 부재(band엔 ρ~수십/µm² 필요).
- **SF active gate 재정의** (platform과 공유).
- **native-scale GPU 확정** (P2c bench: native full-cell ETA/N-scaling — 별도 evidence run).
- **per-SF traction seed-불안정** 해소(full-GPU seed-ensemble 또는 측정 프로토콜) — magnitude 특성화.
- **Slater oracle gate 구현** (cross-model 정성 검증 — 새 gate라 PI 비준).

## 하드룰 / 인프라 주의
- physiological-baseline init, no magic numbers, no gate-loosening. closeout 시 figure +
  3-store Notion + "Notion 업데이트 완료" 수령선. 결정은 PI 고유권한.
- **outputs/ = git-tracked + Syncthing 동시** → gbook git pull이 untracked-outputs 충돌.
  설계는 코드=git / 산출물=Syncthing. gbook은 code@git + outputs@Syncthing으로 완비.
- Dev Logs 페이지 126K(아카이빙 진행 중), tag_kb 옛 duckdb 백업 정리됨.
