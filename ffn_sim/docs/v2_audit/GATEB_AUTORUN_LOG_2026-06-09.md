# Gate-B autorun log — 2026-06-09 (autonomous session)

One line per loop step. Frontier: the **active-myosin-stress vs turgor balance**
(H7_GATE_B_RESULT_2026-06-09.md §5 PRIMARY) — Gate-A REFUTE (γ_soft≈3.06e-3 mN/m,
114× under Hosseini band [0.18,0.40]); Gate-B buckling REFUTE; active floor persists
on connected mesh. Why does myosin contribute ~0 to cortical tension (active stress
~400× < turgor)? Candidates: (1) force/density too low, (2) turgor too high, (3)
local-contraction→spanning-tension conversion failing. Contract: H7_GATE_B_CONTRACT_2026-06-08.md (band LOCKED).

- boot: read AGENTS.md + Gate-A/B result+contract+buckling#1/#2 docs; env ffn_sim hoomd 7.0.1; branch h7/full-cell-integration @ 2ebf437. WIP committed (07a75b8). Key insight: k_on=50/s, k_off0=0.35/s → kinetic equil bound-frac ~99%, yet only ~1-2% engage → under-binding is GEOMETRIC/availability, not kinetics. grip_walk force = k·r delivered straight (r0≈0); MOP faithful. Next: force-budget diagnostic to split engagement vs per-head-force vs transmission.
- loop1 (force-budget audit): built h7_active_force_budget.py — decomposes active-γ floor into engagement × per-head-force × generation, reading live updater + gate's own MOP. Smoke (160 fil, CPU, loading phase): engagement 1.4→8.1% (climbing, geometric not kinetic), per-head T 0.90pN (stall 8.48), γ_soft 1.7e-4≈γ_IK 6.3e-5. **DECISIVE: analytic param envelope γ_active≈½·n2D·f_minifil·ℓ = 0.0101 mN/m = 17.8× UNDER band_lo even at FULL engage+stall.** ⇒ active-γ floor is GENERATION-bound at literature params (Nie-2015 HeLa 0.6/µm², 2pN/head), NOT a transmission/buckling/mesh code bug. Coherent stack: band 0.18 > envelope 0.010 (18×) > Gate-A plateau 0.0031 (60×) > loading 1.7e-4. PI-surface: MCF7 myosin density (HeLa proxy), gel-formula effective-length, band passive/active split. doc H7_ACTIVE_FORCE_BUDGET_2026-06-09.md. NEXT: literature reconciliation of the 3 open items before any param surface.
