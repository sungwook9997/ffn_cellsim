# 🚀 첫 메시지 — Claude Code Terminal에 붙여넣기

아래 전체를 그대로 복사해서 Claude Code Terminal 첫 입력으로 사용:

---

새 프로젝트 시작이야. **ActiveCellSim** — 3D MCF7 spheroid spreading 시뮬레이션.

## 환경
- **Workstation**: Windows 11 + Xeon W-11955M + RTX A5000 (24GB VRAM), via SSH from Mac
- **Backup hardware**: 연구실 RTX 4090 ×2 워크스테이션, Google Colab fallback
- **Mac (M1 Max)**: 분석/시각화 review용
- **이미 셋팅됨**: Tailscale SSH, conda env (`spheroid` on Mac, `pytorch 2.7.1+cu118` on Win), `train_win`/`trainlog` aliases

## 프로젝트 폴더
프로젝트 풀 폴더는 **별도 zip 파일**로 받아서 압축 풀었어.
구조:
```
ActiveCellSim/
├── CLAUDE.md            # 이 파일을 먼저 읽어줘 (전체 컨텍스트의 onboarding)
├── README.md            # 프로젝트 한눈 요약
├── docs/                # 14개 markdown 파일 — 필요할 때만 읽어 (progressive disclosure)
│   ├── 00_project_vision.md          ⭐ framing & assumptions 마스터
│   ├── 01_model_architecture.md
│   ├── 02_force_models.md            (모든 force term + IF≥15 references)
│   ├── 03_adhesion_dynamics.md       (φ ODE)
│   ├── 04_simulation_setup.md
│   ├── 05_radial_approximation.md    (radial vs full 비교 derivation)
│   ├── 06_radial_full_comparison.md
│   ├── 07_internal_flow_dynamics.md  (Marangoni / nematic / vortex)
│   ├── 08_mechano_osmotic.md         (volume regulation Tier 2)
│   ├── 09_visualization.md           (vispy + Blender + plotly)
│   ├── 10_dev_roadmap.md             ⭐ 다음 무엇을 할지 (canonical)
│   ├── 11_performance_protocol.md    (GPU portability, benchmark)
│   ├── 12_validation.md              ⭐ assumption 마스터 리스트
│   ├── 13_data_schema.md
│   └── references.bib
└── data/
    └── experimental/    # PI 실험 데이터 (READ-ONLY, fitting 금지)
```

## 너의 첫 작업
1. `CLAUDE.md`만 먼저 읽어 (200줄 이하, 전체 onboarding).
2. 그 다음 `docs/10_dev_roadmap.md`만 읽어 — Stage 0 시작이야.
3. 다른 docs는 **필요할 때만** 읽어 (Progressive Disclosure 원칙). 한꺼번에 다 로드하지 마.

## 핵심 제약 (절대 어기지 마)
- 실험 데이터 (`data/experimental/*.csv`) 는 **fitting 금지**, 시각화 오버레이만 허용
- 모든 force term은 IF≥15 학술 reference 필요 (`02_force_models.md` 참고)
- 5,000 cells × 80hr 시뮬레이션, pilot 1,000 cells × 4-8hr → wall-clock 30분 안 끝나면 단순화
- GPU portable 필수 (A5000 / 4090×2 / Colab 호환)
- Headless SSH 환경 — GUI 없음, 시각화는 frame export → MP4

## 워크플로우 분담
- **Claude Desktop App (별도 대화창)**: 물리 모델 review, 결과 해석, 논문 작업 — brain trust
- **VSCode + Claude Code 확장**: 코드 inline 수정
- **너 (Claude Code Terminal)**: 환경 셋업, 시뮬 실행, GPU 모니터링, 배치 운영 — execution

## 시작 명령
Stage 0 (환경 셋업)부터 시작해. `docs/10_dev_roadmap.md`의 Stage 0 체크리스트를 따라.

질문은 작업 진행 중에 막힐 때만. 일반적으로는 자율 판단해서 진행해. 막연한 결정이 필요할 때만 확인 요청.

준비 됐으면:
1. `CLAUDE.md` 읽기 → 한 줄로 "WHY/WHAT/HOW 이해됨" 확인
2. `docs/10_dev_roadmap.md` Stage 0 읽기
3. Stage 0 시작 (requirements.txt, conda env, Taichi GPU 검증, configs/ 폴더, 등)

시작해줘.

---

(끝)
