---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# FF 단일세포 엔진 상태 감사 (2026-07-07)

> 36-agent 워크플로 감사(컴포넌트 그룹 10개 병렬 → 각 "작동함" 주장 adversarial 재검증 → 종합). 58 컴포넌트 감사, 8건의 "active" 주장이 재검증에서 뒤집힘(REFUTE). 원칙: **"loop에서 launch됨 ≠ 세포에 힘이 도달함"**.
>
> **Lead가 코드로 직접 재확인한 load-bearing 주장 2건 (agent 출력 맹신 금지, MT 교훈):**
> - **Protrusion이 GPU-implicit 경로에서 버려짐** — `leading_edge_push_kernel`/`protrusion_reaction_kernel`은 main loop(driver:380/382)과 host-implicit `full_force`(driver:253/255)에서만 launch되고, GPU-implicit RHS `gpu_force_fn`(driver:287-332)은 line 291에서 f_d를 재-zero한 뒤 protrusion을 재계산하지 않음(spreading만). ⇒ `--implicit --device cuda`(네이티브 crawl에서 내가 쓴 config)에서 protrusion이 solve에 안 들어감. **네이티브 disp=0의 직접 원인.** (내 앞선 "traction 피스" 설명은 불완전했음.)
> - **핵 bead는 오직 핵 인덱스 `[Ne,Ne+n_nuc)`에만 힘을 씀** (`nucleus_shell_kernel`: `i = off + j`, network_warp.py) — cortex 반작용 0 ⇒ MT와 동일한 inert passenger.

---

Driver 기준: `aleph/scripts/ff_crawl_on_substrate.py`. 모든 주장은 audit의 file:line 근거에 고정. "loop에서 launch됨" ≠ "세포에 작동함" 원칙으로 분류.

---

## 1. 엔진 한 줄 요약

FF는 Cytosim 물리(NF2007)를 **NVIDIA Warp**로 재구현한 GPU-native 단일세포 엔진이다. 세포 몸체는 구면 위 great-circle actin 섬유 + α-actinin/filamin 가교 + myosin 링크로 짜인 **cortex fiber-network** (`build_crosslinked_cortex`, driver:73-74) 하나가 전부이고, 이것이 `pos_all[:Nc]` (driver:92)로 매 스텝 적분되는 유일한 탄성체다. 두 개의 적분 경로가 있다: **default = explicit overdamped** (`axpy_physical_kernel`, driver:431, dt≈5.5µs, CFL-bound) 와 **`--implicit --device cuda` = GPU-resident semi-implicit** (`ff_implicit_step_gpu`, driver:403, dt=1e-2s, native-scale용). Bending·crosslink만 analytic K에 들어가고(`assemble_K_current_cupy`) 나머지 힘은 explicit RHS로 남는다. **결정적 사실: driver의 문서화된 production 명령(driver:18)은 `--device cuda:0`는 쓰지만 `--implicit`를 넘기지 않는다** → 기본 실행은 explicit 경로다. 이 비대칭이 아래 (B)의 여러 함정을 만든다.

---

## 2. (A) 잘 구현 + 배선 + 실제 작동 (verify-confirmed `yes-active`)

Cortex node block `[0,Nc)`에 매 스텝 힘을 쓰고, 그 힘이 적분되어 위치가 바뀌는 것이 코드로 추적된(verify=CONFIRMED) 컴포넌트들. 이들이 "진짜 세포"다.

| Component | 무엇 | Force path (file:line) | Gate |
|---|---|---|---|
| **Cortex fiber-network** (`build_crosslinked_cortex`) | 탄성 몸체 그 자체 — actin 섬유+가교+myosin 링크 | driver:73-74 빌드 → `pos_all[:Nc]`(92)가 매 스텝 적분되는 몸체; 관측 disp/v_crawl = `pos[:Nc].mean` (468) | on-shell residual + areal-density (cortex_assembly.py:213-229); test_cortex_assembly.py |
| **Cytosim bending** (`cytosim_bending_kernel`) | NF2007 곡률 복원력 | driver:353(explicit)/292(gpu-impl)/227(host-impl), dim=nT; cortex 삼중항만 인덱싱 → f_d → 적분 | PASS: Cytosim parity + Euler π²κ/L² anchor |
| **Crosslinker springs** (`link_spring_kernel`) | α-actinin/filamin Hookean 가교 | driver:354, dim=n_xl; `xl_i/xl_j`(136-140) 모두 `<Nc`; k_xl.max()가 CFL dt 결정(142) | PASS: numpy bit-parity (test_network_warp.py) |
| **Myosin prestress** (`myosin_kernel`, DEFAULT) | 등력 수축 dipole (NMIIA) | driver:355-359; `myosin_linear=False` default(107); has_myo=True; f_myo=5pN | γ-floor는 REFUTE(~530× under)지만 힘은 실제 도달. **⚠️ 5 links만(n_myo=900//160), floor-magnitude** |
| **Osmotic turgor + 막 안쪽 장력** (`turgor_kernel` ×2) | Guo-2017 삼투압 + Young-Laplace 막장력 | driver:360(삼투)/361(막), dim=Nc; V는 `cortex_volume_kernel`(347) cortex faces만 | PASS: ΔP·R/2 consistency, K_vol=Π_in0 유도 |
| **Inextensibility reshape** (`reshape_kernel`) | NF2007 §5.3 hard 길이 제약 투영 | driver:436, 매 20스텝; cortex 노드 직접 재배치; 유일한 비신축 제약 | PASS: NF2007 project-then-reshape |
| **Template-1 lumped 막장력** (γ_mem) | **production의 유일한 막 힘** — ΔP=−2γ_mem/R 내향 압력 | driver:341(dP_mem_area)→361(turgor_kernel); γ_mem=10 pN/µm band-guarded | PASS param-band (test_compartments_shared.py). ⚠️ 상수 baseline만 |
| **FA integrin clutch spring** (`clutch_spring_kernel`) | 기저 cortex 노드↔고정 기질 앵커 traction | driver:369-371; `clutches=True` 하드코딩(531); actin=`S['basal']`⊂[0,Nc)(188/81); 앵커 고정 → **net external 힘** | `--audit` clutch-OFF drift≈0 oracle(554-560) |
| **Clutch catch-slip KMC (detach)** (`clutch_catchslip_kmc_kernel`) | Pereverzev off-rate로 loaded clutch 파열 → bd_d flip | driver:444-448; rupture=not mature=True default; bd_d가 clutch_spring+K diag(399-401) gate | 같은 --audit oracle |
| **Basal FA clutch treadmill (traction)** | 위 두 clutch를 crawl traction으로 조합 | 모든 경로: explicit 369-371 / gpu 307-309 / host 242-244 | traction 2.5-3.8 nN 기록됨. ⚠️ bound_frac=1.0 → treadmill 순환 미검증 |
| **Crosslink turnover** (`xl_turnover_kernel`) | α-actinin k_off Maxwell 이완 → cortex 유동화 | driver:437-438 무조건; r0_d 변경 → link_spring가 매 스텝 읽음(354/293/228) | frac_xl=1−e^(−koff·Δt) 유도 |
| **Physical η-overdamped 적분기** (`axpy_physical_kernel`) | **default 적분기** — x+=(dt/γ_i)·F_i, per-node cytoplasm drag | driver:431, dim=N; f_d의 모든 힘을 위치로 변환; found_spread.json dt=5.48e-6 실측 | CFL-stable by construction |
| **hand_kmc INTEGRIN_A5B1 preset** (parameter provider) | clutch 상수 공급 (Kong2009, KB-2.5) | fa_clutch_warp.py:26→resolve_clutch:97-103→cp; **데이터는 active, 런타임 KMC 클래스는 not-wired** | KB-2.5 anchoring; F*≈7 vs 30 pN flag |

**작동하는 motility (부분):** `spreading_push_kernel`(driver:372-377 explicit **및** gpu_force_fn 310-315)는 default-off(`--spread`)지만 두 기록 run(found_spread/spread_impl)에서 end-to-end 실행됨 — net-zero라 translocation 아닌 flatten/widen. `crosslink turnover`와 함께 이것이 실제 굴러본 유일한 motility 모드.

**핵심:** 실제 작동하는 것은 **cortex 하나 + 그 위의 bending/crosslink/myosin(floor)/turgor/막장력(lumped)/reshape/turnover + 기저 FA clutch traction** 뿐이다. 이것만이 "세포"다.

---

## 3. (B) 구현됐지만 세포에 안 붙음 — MT-style 함정

verify가 prior audit을 **뒤집은(REFUTE)** 항목이 특히 위험하다. 세 부류: passenger(K/loop 안이지만 decoupled) · default-off(production이 안 켬) · diagnostic/not-wired(힘 경로 자체 없음).

### 3-1. Passenger — loop 안에 있는데 cortex trajectory에 0 영향

**① 핵 bead-cloud** (`nucleus_shell_kernel` + `_seed_nucleus_cloud`) — **verify가 "partial"→`no-inert-passenger`로 REFUTE.**
- 하는 일: 3000 bead를 R_nuc=0.70R 구면에 bilinear(chromatin+lamin) 복원력으로 유지. driver:83-85 무조건 seed, 매 스텝 launch(238/303/365).
- 왜 안 붙나: 커널이 힘을 **핵 인덱스 `[Ne,Ne+n_nuc)`에만** 쓴다(network_warp.py:389,412). cortex `[0,Nc)`에 반작용 0. 핵-cortex bond·excluded-volume·LINC **전무**. `centre`는 cortex-only centroid의 **일방향 read**. turgor는 `ConvexHull(pos[:Nc])`(131)라 핵이 cortex 부피조차 안 밀어냄. Implicit K는 핵 DOF에 bending/xlink/vol_g 항이 없어 block-diagonal → coupling 소멸. **핵 on/off와 무관하게 cortex 궤적 동일 = MT-aster 패턴 그대로.**
- 활성화 배선: (a) 최소 — cortex 부피 결합 `V_cyto=V_hull−V_nuc`. 이미 `simulate_whole_cell_compression_on_device`(network_warp.py:436-441)에 존재하나 crawl driver가 호출 안 함. (b) 정식 — Thread-C off-diagonal K (cortex↔lamina LINC bond) 또는 핵bead↔cortex soft contact/EV 커널. 파일: `ff_crawl_on_substrate.py` + `network_warp.py`.

**② MT aster 전체 (6개 컴포넌트)** — default-off **이면서** 켜도 disconnected passenger.
- geometry(build_microtubule_aster) · merge into K(merge_aster_into_cortex) · MT bending · MTOC hub crosslink · MT drag diag 보정 · 검증 gate.
- 왜 안 붙나: `concat_fiber_networks`가 bend_triple을 nbase만큼 offset(fiber_network.py:114-115) → cortex 삼중항 `[0,Nc)`, MT 삼중항 `[Nc,Ne)`, **cortex↔MT를 잇는 삼중항 0개**. hub crosslink는 `hub_i=MTOC`↔`hub_j=arm-base` **둘 다 MT 도메인**(microtubule.py:95-99). cortex crosslink는 MT 노드를 절대 참조 안 함. 결과: merged K가 {cortex+nucleus} vs {MT+MTOC}로 **block-diagonal** → cortex 궤적이 MT block과 수학적으로 독립. docstring이 약속한 tensegrity(액토미오신 균형·핵 위치잡기·confined migration 하중)는 **하나도 안 지어짐**.
- 활성화 배선: (a) **MT-tip↔cortex contact** — tip 노드가 cortex 표면 침투 시 soft 일방향 반발을 `gpu_force_fn`의 `F_total`에 추가(K 아님, Thread-C 계획 §4). (b) **MTOC↔핵 envelope bond** (LINC/nesprin). 둘 중 하나 없으면 merge는 영구 inert. 파일: `ff_crawl_on_substrate.py` gpu_force_fn + `microtubule.py`.

**③ Leading-edge protrusion** (`leading_edge_push_kernel` + `protrusion_reaction_kernel`) — **verify "partial" 확정, 그러나 production GPU 경로에서 passenger.**
- 왜 위험: explicit·host-implicit 경로에선 f_d에 실려 front cortex 노드를 움직임(driver:378-382→431). **그러나 GPU-implicit 경로**(native-scale에 CLAUDE.md가 요구하는 유일한 config)는 `ff_implicit_step_gpu`(403)가 `gpu_force_fn`(287-329)에서 f_d를 **재-zero(291)하고 protrusion을 재계산 안 함** — spreading_push만 있음(310-315). loop body가 쓴 f_d는 버려진다(다음 read는 pos_d, 384). → **`--implicit --device cuda`에서 protrusion on/off와 cortex 궤적 동일 = 조용한 MT 패턴.** 이것이 crawl driver의 headline 힘인데 production solver에서 죽어 있음.
- 활성화 배선: **`gpu_force_fn`(driver:287-329)에 `leading_edge_push_kernel`+`protrusion_reaction_kernel` launch 추가** — 바로 위 spreading_push 블록(310-315)을 그대로 미러링. 한 곳 수정. 이게 (5)의 최우선 권고다.

### 3-2. Default-off — production이 절대 안 켜는 모듈 (physiological-baseline 규칙 위반)

| Component | 상태 & 왜 죽었나 | 활성화 배선 |
|---|---|---|
| **minifilament_kernel** (Stam-Hocky FV, fine-grained myosin) | `--myosin-linear` default-off(509). **verify가 "partial"→`not-wired` REFUTE.** 켜도 `vslide_d`가 zeros로 alloc(186)되고 **어디서도 갱신 안 됨** → Fmag=f_stall·(1−0)=60pN 상수, FV 법칙 dead code. | flag default 켜기 + **v_slide를 실제 앵커 kinematics로 계산**(dead 로직 살리기). 파일: `ff_crawl_on_substrate.py`:186 |
| **FA maturation** (talin→vinculin→growth) | `--fa-maturation` default-off(511). **verify "partial"→`diagnostic-only` REFUTE.** 켜도 `nv_d` write-only(458, 다시 read 안 함), `k_int_effective`가 driver에 import 안 됨(47-48) → talin/vinculin 힘 0. 유일 실경로 `fa_disassemble_kernel`(461)도 8000스텝에 3틱뿐이라 실질 미발동. | flag 켜기 + `k_int_effective` import해 clutch stiffness에 feedback. 파일: `fa_maturation.py`:127 → clutch_spring |
| **Compliant substrate** (Winkler 앵커) | `substrate_E=0.0` default(108/510) → substrate=None(193) → equilibrium 커널 무발동(432). **verify가 `not-wired` REFUTE.** 배선은 정확(432-434→clutch_spring가 앵커 read)하나 gate 꺼짐. | **physiological-baseline 규칙대로 E=5kPa(PAA, KB-1.5)로 default 켜기** + implicit K diag를 series stiffness `k_int·k_sub/(k_int+k_sub)`로(현재 full k_int, driver:401). Thread-A. |
| **Per-filament barbed-end growth** (`--growth`) | default-off(108), 두 기록 run에서 미발동. 켜도 reshape 통한 간접 경로. | `--growth` 켜기 + seg_max 넘으면 bead insertion(미구현, 216) |
| **Actin assembly area growth** (`--assembly`) | default-off(108). `--assembly`+`--implicit` 조합에서만 volume projection 통해 효과. | flag 켜기 + host-implicit volume 결합 |
| **GPU/Host implicit solver + K assembly** (`ff_implicit_step_gpu`, `assemble_K_current[_cupy]`, `implicit_step_current`) | **verify 4건 모두 "yes-active"→`partial` REFUTE.** `--implicit` unset이라 default는 explicit else-branch(430). production 명령(18)도 `--implicit` 생략 → K 조립 자체가 안 돎. host 버전은 CPU-dev 전용(GPU mandate가 배제). | production 명령에 `--implicit` 추가 (native-scale·MT에 필수, 523-524) |

### 3-3. Diagnostic-only / Not-wired — 힘 경로 자체가 없거나 driver가 import 안 함

- **Piezo1 reporter** (`piezo.py`): `--piezo` default-off, run() **끝난 뒤** 단발 print(543-545), 상수 g_mem에 평가, `CA_FEEDBACK_GAIN=0.0`(piezo.py:30). per-face `piezo_popen_kernel`은 무발동. docstring도 "bears no mechanical load". → 순수 reporter, 세포에 0 영향.
- **γ-floor harness** (`gamma_floor.py`의 measure_gamma/production/equilibrate): driver는 상수+builder만 import(30-31). 별도 offline 실험. γ_myo ~530× under band = REFUTE(정직 보고, tune 안 함).
- **γ 추정기** (`gamma_estimator.py`), **network_contractility patch**, **kim_network oracle**, **MT euler/L_p gate**, **constraints.py**, **relax.py**, **polymerization_warp 커널**, **polarization_activegel**: 전부 driver가 import 안 함 (측정기/oracle/중복 구현). polarization은 특히 — crawl의 극성축 `phat`이 **하드코딩 (1,0,0)**(driver:57), Bois2011 instability의 emergent cap이 phat에 전혀 안 먹임.
- **Template-2 독립 막 sheet** (`membrane_surface.py`: MembraneMesh/area kernel/**erm_tether_kernel**/**containment_kernel**/reservoir_tension): crawl driver가 통째로 import 안 함(막 import는 `resolve_membrane`만, 50). ERM은 **막↔cortex를 잇는 유일한 경로**인데 production에 부재 → docstring의 "막이 cortex를 담는다" 주장이 production에선 **아무 힘도 안 함**. `ff_membrane_cell.py` 데모에만 존재.
- **nucleus_envelope.py**: `_gpu_port_check.py`의 detached sphere smoke test만 호출. rupture flag는 read되나 set되는 곳 없음, bending 커널 없음.
- **fa_anchor.py** (host prototype), **FFImplicitStepper** (prefactored LU), **weave()/architecture_spec**, **Arp2/3 branch·soft_contact·AFM plate·Bell turnover 등 미사용 network_warp 커널**: 전부 dead code w.r.t. crawl.
- **ecm_mikado.py**: crawl이 import 안 함 (별도 driver `ff_protrusion_into_ecm.py`에만). crawl 세포가 닿는 기질은 **매끈한 rigid plane**(`substrate_plane_kernel`)뿐, fiber ECM 없음.

---

## 4. (C) 놓친 것 / 미완성 — master plan 7 compartment 대비

Master plan(COMPARTMENT_BUILD_PLAN)은 20 unit(감사: 6 present/16 partial/10 absent)을 세운다. 7 핵심 compartment 축으로 정리:

| Compartment | 계획된 것 | 현재 상태 | 미완/누락 |
|---|---|---|---|
| **① Substrate (기질)** | compliant Winkler 앵커 + ligand coat(Bare/Pre/Lam4) + Thread-A | 코드 존재하나 **default E=0 rigid** | Thread-A(implicit 앵커 DOF) 미구축; ligand preset 미구축; production이 rigid로 실행 |
| **② Cortex (액토미오신)** | Hill/Stam-Hocky myosin minifilament | **LIVE (유일한 실compartment)**, 단 myosin은 5 links·5pN floor | Hill FV 미배선(minifilament v_slide≡0); myosin 밀도 floor |
| **③ FA (부착반)** | 성숙(talin-vinculin) + ECM 왕복 traction | clutch spring/KMC는 LIVE | maturation default-off+inert; **ECM로의 Newton-3rd 반작용 미배선** |
| **④ ECM (fibrillar)** | Mikado collagen에 clutch 반작용 | 빌드됨(`ecm_mikado.py`), 별도 driver에만 | crawl에 미import; 왕복 remodeling 없음 |
| **⑤ Microtubule** | aster + MTOC + 압축 strut + centering | 빌드됨, **disconnected passenger** | MT-tip↔cortex contact 없음, MTOC↔핵 bond 없음, DI kinetics 미구현 |
| **⑥ Nucleus (핵)** | envelope surface + lamina split + poroelastic 내부 + LINC + nucleolus + IF cage | bead-cloud만, **inert passenger** | **cortex 결합 전무**; envelope surface(Thread-C)·lamina·poroelastic·LINC·nucleolus·IF cage 전부 미배선 |
| **⑦ Membrane (막)** | 독립 sheet + ERM + containment + reservoir + Helfrich κ + spectrin | **Template-1 lumped γ_mem만** | Template-2 sheet·ERM·containment·reservoir·bending·spectrin 전부 not-wired (Thread-D) |
| (추가) **Cytosol 점탄성** | G′-memory(SLS)+Darcy poroelastic | — | **아예 미구축**; 현재는 per-node drag 상수만 |
| (추가) **Piezo feedback** | P_open→Ca→contractility | reporter만, gain=0 stub | 신호 loop 미구축 |
| (추가) **Stress fiber / Glycocalyx / Fascin / Septin** | — | spread morphology·KB에 blocked | 계획상 후순위(§3 BLOCKED) |

**미구축 cross-cutting Threads (각각 여러 compartment를 gate):**
- **Thread A** (implicit 앵커 DOF) — Substrate·ECM. 미구축.
- **Thread B** (nonlinear per-bond tangent) — IF cage·NE lamina. 미구축.
- **Thread C** (multi-node-set off-diagonal K) — **핵↔cortex LINC 결합의 핵심**. 미구축 → 핵이 passenger인 근본 원인.
- **Thread D** (독립 막 sheet를 production에 승격) — ERM·reservoir·spectrin·septin·local Piezo. 미구축.
- **Thread E** (per-node mobility array) — nucleolus·poroelastic 내부. 미구축.

---

## 5. 우선순위 권고 — 최소 배선, 최대 기능 효과

파라미터 튜닝은 제외. 아래는 전부 "죽은 경로 살리기 / 누락된 실제 힘 추가 / 실버그 수정"에 해당.

**P0 — Protrusion을 GPU-implicit 경로에 배선 (한 곳, 최대 ROI).**
현재 crawl의 headline 힘(leading-edge push)이 **production solver에서 버려진다**: loop body가 f_d에 쓰지만 `ff_implicit_step_gpu`가 `gpu_force_fn`(driver:287-329)에서 f_d를 재-zero하고 protrusion을 재계산 안 함(spreading_push만 존재, 310-315). → **`gpu_force_fn`에 `leading_edge_push_kernel`+`protrusion_reaction_kernel` launch 추가**, 바로 위 spreading 블록을 미러링. 파일 `ff_crawl_on_substrate.py`, ~10줄. 효과: `--implicit --device cuda` production run에서 protrusion이 실제로 살아난다 (지금은 explicit dev 경로에서만 작동).

**P1 — 핵↔cortex 결합 (Thread-C 최소판).**
핵이 완전 passenger인 근본 원인. 최소 배선: 이미 존재하는 hydrostatic 결합 `V_cyto=V_hull−V_nuc`(network_warp.py:436-441, `simulate_whole_cell_compression_on_device`)를 crawl loop에 끌어오기 — cortex turgor의 V0를 핵 부피만큼 보정. 정식: 핵 bead↔cortex soft contact/EV 커널 또는 LINC bond(off-diagonal K). 파일 `ff_crawl_on_substrate.py` + `network_warp.py`. 효과: 3000 bead가 처음으로 cortex에 반작용 → 압축/이동 시 핵이 실제 하중을 실음.

**P2 — 독립 막 sheet + ERM를 production에 승격 (Thread-D).**
production 막은 lumped 스칼라(γ_mem) 하나뿐이고, 막↔cortex 결합(`erm_tether_kernel`)은 crawl이 import조차 안 한다. `ff_membrane_cell.py`에 이미 프로토타입된 sheet+ERM+containment를 driver에 배선(THREADC_UNIFIED_PLAN Stage-2가 정확한 배선 순서 제공: `--membrane` off면 bit-identical, on이면 lumped proxy 대체). 파일 `ff_crawl_on_substrate.py`가 `ff/membrane_surface.py` import. 효과: 막이 스칼라 압력에서 실제 compartment(면적탄성+ERM tether+bleb 파열+containment)로 승격.

**P3 — Substrate compliance를 physiological E로 켜기 (Thread-A, 규칙 준수).**
배선은 이미 정확(driver:432-434→clutch_spring가 movable 앵커 read)하나 `substrate_E=0.0` default라 rigid pin으로 죽어 있다. **physiological-baseline HARD 규칙**(additive 모듈은 production에서 생리값으로 ON)에 따라 E=5kPa(PAA)로 default 승격 + implicit K diagonal을 series stiffness `k_int·k_sub/(k_int+k_sub)`로 수정(현재 full k_int, driver:401 — 앵커 갱신과 Jacobian 불일치). 이건 gate 통과용 튜닝이 아니라 규칙 준수 + 실버그 수정. 파일 `ff_crawl_on_substrate.py`.

**P4 (낮음) — MT-tip↔cortex contact + MTOC↔핵 bond.**
tensegrity를 살리려면 필요하나 crawl motility에는 ROI 낮음(계획도 Stage-1 이후로 defer). tip soft 반발을 `gpu_force_fn`의 F_total에 추가(K 아님).

**요약:** 가장 작은 변경으로 가장 큰 차이는 **P0 (protrusion을 GPU 경로에 배선)** — 이건 새 물리가 아니라 production solver에서 이미 계산되는 힘이 버려지는 것을 잇는 한 곳 수정이며, 현재 native-scale crawl이 protrusion 없이 도는 조용한 실패를 고친다. 그 다음 P1·P2가 "passenger 3000 bead"와 "lumped 막"을 실제 compartment로 승격시키는 최소 배선이다.
