# Deck assets — how to use with Claude Design

**Claude Design은 로컬 절대경로를 직접 못 읽는다.** 아래 파일들을 Finder에서 끌어다 업로드하면 된다.

## 두 가지 사용법
1. **빠른 길 — 큐레이션 세트:** `deck_assets/curated/` 폴더를 통째로 열어서, 슬라이드 번호 순서대로 드래그해 넣는다. 이미 슬라이드에 매칭된 그림만 25개.
2. **더 고르고 싶으면 — 갤러리:** `deck_assets/GALLERY.html` 를 브라우저로 연다 (더블클릭). 테스트하면서 나온 **418개 전부**(355 png · 55 movie · 8 gif)를 폴더별 썸네일로 훑어볼 수 있다. 마음에 드는 타일 밑 회색 경로를 보고 Finder에서 그 파일을 찾아 드래그.

## 동영상
Claude Design이 mp4를 못 받을 수 있어서 각 핵심 무비를 **포스터 PNG + 압축 GIF** 둘 다 만들어 뒀다 (`curated/movie_*`). 슬라이드엔 GIF를 넣으면 움직이고, 정적이 필요하면 `_poster.png`를 쓰면 된다.

## 슬라이드 ↔ 파일 매핑

| 슬라이드 | 큐레이션 파일 | 원본 |
|---|---|---|
| 1 (타이틀/히어로) | `slide01_hero_spheroid_n1138.png` | warp_decohesion/figs/ipc_sphere_clean_n1000_mesh.png |
| 2 (PI 실험) | `slide02_pi_experiment_overlay.png` | layer2/figs/fig_layer2_pi_overlay.png |
| 7–8 (cortex 검증) | `slide07-08_cortex_validation_montage.png` | h3/figs/fig_h3_evidence_montage.png |
| 8 (Young–Dupré) | `slide08_young_dupre_doublet.png` | warp_decohesion/figs/young_dupre_doublet_gate.png |
| 9 (R-법칙 r²=0.998) | `slide09_rlaw_r2_0998.png` | layer2/figs/fig_layer2_prod_rlaw.png |
| 9 (리간드 순서) | `slide09_ligand_ordering.png` | layer2/figs/fig_layer2_ligand_production.png |
| 10 (크기 갭) | `slide10_magnitude_gap.png` | layer2/figs/fig_layer2_magnitude_gap.png |
| 10 (결정적 견인/이젝션) | `slide10_decisive_traction_ejection.png` | layer2/figs/fig_layer2_decisive_traction.png |
| 12 (two-stage 스토리보드) | `slide12_two_stage_storyboard_n400.png` | h_dcm_two_stage/figs/two_stage_n400.png |
| 13 (retraction: peeling) | `slide13_retraction_peeling_vs_cohesive.png` | warp_decohesion/figs/cohesion_fix_peeling_vs_cohesive.png |
| 13 (7 mobile / 93 frozen) | `slide13_actual_cells_7mobile_93frozen.png` | warp_decohesion/figs/A1_actual_cells_3d.png |
| 14 (faceting) | `slide14_faceting_comparison.png` | warp_decohesion/figs/faceting_comparison.png |
| 14 (단면 surf/core) | `slide14_xsection_surf_core.png` | warp_decohesion/figs/RECONSIDER_xsec_surfcore.png |
| 16 (엔진 스케일) | `slide16_engine_scale_n1128.png` | warp_decohesion/figs/ipc_cleave_opt_n1000_mesh.png |

### 무비 (포스터 + GIF)
| 용도 | 파일 | 원본 |
|---|---|---|
| 오프너/응집 | `movie_aggregation_n400.gif` / `_poster.png` | h_dcm_two_stage/figs/two_stage_n400_agg.mp4 |
| 슬라이드 12 (압축) | `movie_spread_compaction_n800.gif` / `_poster.png` | h_dcm_two_stage/figs/two_stage_n800_spread.mp4 |
| 백업 (de-cohesion) | `movie_decohesion_junction_switch.gif` / `_poster.png` | h_dcm_gpu_lod/figs/junction_switch_spatial.mp4 |
| 백업 (necrosis) | `movie_necrosis_3zone.gif` / `_poster.png` | h_dcm_gpu_lod/figs/necrosis_3zone_spatial.mp4 |

### 백업 스틸
- `backup_cleave_mechanism.png` — 세포 분열 메커니즘
- `backup_decohesion_montage_n100.png` — de-cohesion 몽타주
- `backup_necrosis_3zone.png` — necrosis 3-zone

## 재생성
원본 그림이 갱신되면: `bash <scratchpad>/stage_assets.sh` 다시 실행 (또는 새 그림은 GALLERY.html에서 바로 보임 — 갤러리만 다시 빌드하려면 build_gallery.py).
