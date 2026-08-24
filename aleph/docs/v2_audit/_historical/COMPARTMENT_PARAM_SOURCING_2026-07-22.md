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

# Compartment parameter sourcing — DOI-verified dossier (2026-07-22)

> Overnight parallel literature sourcing (7 research agents) for the missing/partial `ac/` compartments, to
> unblock implementation under the magic-number rule. Every value is DOI-verified against a primary source or
> explicitly flagged `[unverified]`/`[proxy]`/`[MCF7-absent]`. **Before citing any of these in a deliverable,
> register as a Notion SourceEvidence/KnowledgeClaim and confirm verdict=OK (`verify_sources.py`)** — this doc
> is the staging ground, not the SoT. Units: µm-pN-s (1 Pa = 1 pN/µm²; k_BT₃₇ = 4.28e-3 pN·µm).
> Companion to [`COMPARTMENT_COMPLETION_ROADMAP_2026-07-22.md`](COMPARTMENT_COMPLETION_ROADMAP_2026-07-22.md).

## Microtubule (MT) — IMPLEMENTED 2026-07-22 (`ac/solid/microtubule.py`)
- **EI = 22 pN·µm²** (Gittes 1993, `10.1083/jcb.120.4.923`); L_p = 5.2 mm; E ≈ 1200 pN/µm². ff KAPPA_MT=20 in-band.
- **Length-dependent softening** short MTs: L_p 110 µm@2.6µm → 5 mm@47µm (Pampaloni 2006, `10.1073/pnas.0603931103`) — deferred refinement.
- **In-cell buckling reinforcement ~100×, λ≈3 µm** (Brangwynne 2006, `10.1083/jcb.200601060`) — NEVER bare Euler in the full cell; emergent from cortex.
- **Motors**: dynein ~1.1 pN, kinesin 5–7 pN (Leidel 2013, `10.1073/pnas.1219961110`) — deferred.
- Geometry: OD 25 nm / 13-PF / 1625 dimers·µm⁻¹. **Count MCF7-absent → proxy ~250–600/cell (PI GAP).**

## Intermediate filaments (IF) — keratin(MCF7) vs vimentin(MDA-231) EMT split — NEXT
- **Persistence length**: vimentin l_p ≈ 0.5 µm (Lin 2010 `10.1016/j.jmb.2010.04.054`); keratin K8/K18 ≈ 0.30 µm (+Mg 0.48) (Lichtenstern 2012 `10.1016/j.jsb.2011.11.003`). → κ = k_BT·l_p: vimentin 2.0e-3, keratin 1.2e-3 pN·µm².
- **Stretch modulus** ~6–9 MPa (network, Lin 2010) / ≤~1 GPa (stabilized, Guzmán 2006 `10.1016/j.jmb.2006.05.030`) — axial-sliding regime split, implement BOTH not lumped.
- **Extensibility ≥160% (mean 2.6×), up to 260–300%** via α→β transition (Kreplak 2005 `10.1016/j.jmb.2005.09.092`; Qin-Buehler 2009 `10.1371/journal.pone.0007294`; rate-dependent multi-state Block 2017 `10.1103/PhysRevLett.118.048101`). ← the nonlinear force-extension primitive.
- **EMT**: keratin loss softens epithelial ~30–60% + invasive (Ramms/Seltmann 2013 `10.1073/pnas.1310493110`); vimentin doubles cytoplasm G′ ~5→10 Pa (Guo 2013 `10.1016/j.bpj.2013.08.037`), large-strain safety net (Hu 2019 `10.1073/pnas.1903890116`), protects nucleus in confinement (Patteson 2019 `10.1083/jcb.201902046`).
- Perinuclear plectin cage + B-lamin–VIF LINC (Vahabikashi 2022 `10.1073/pnas.2121816119`). **Mesh size = [unverified] PI GAP.**

## Nucleus completion — chromatin WLC + nucleoplasm (Stephens–Marko two-regime)
- **Chromatin governs SMALL strain (<3 µm ≈ <30%)**: k ≈ 0.70 nN/µm = **700 pN/µm** (HeLa; Stephens 2017 `10.1091/mbc.E16-09-0653`). Crosslinked-polymer interior + shell (Banigan 2017 `10.1016/j.bpj.2017.08.034`).
- **Lamin-A/C governs LARGE strain**, strain-stiffens 1.5–2.5× (long-ext 0.85 nN/µm); knee ~30% (matches existing `knee_strain~0.10`… actually onset ~30% global — keep distinct).
- **Nucleoplasm**: bulk viscoelastic η≈52 Pa·s, G′≈18 Pa solid-like (Tseng 2004 `10.1242/jcs.01073`); **pore viscosity ≈ 2–8× water** (nucGEMs preprint `10.1101/2021.11.18.469159` [CHECK]). Wire pore→fluid channel, 18 Pa→bulk.
- **Lamina area modulus 25 mN/m** (Dahl 2004 `10.1242/jcs.01357`). **NE rupture ~1% areal transient tensile; lamin-A/C protects, lamin-B doesn't** (Zhang-Lele 2018 `10.1091/mbc.E18-09-0604`) — distinct from the strain-stiffening knee.
- **MCF7 geometry**: R_nuc ≈ 5.1–6.3 µm (diam 10–13); **6 µm better-centered than provisional 5**. ⚠️N:C=0.68 is DIAMETER ratio → volume N:C≈0.30 (pin the definition!). Oblate aspect [unverified] PI GAP. (Strohm 2016 `10.1007/s10765-016-2129-y`.)

## ECM + integrin clutch + focal adhesion (highest impact; ff/ physics exists)
- **Collagen-I fibril E ≈ 0.5 GPa** hydrated (van der Rijt 2006 `10.1002/mabi.200600063`; sweep 0.1–0.8); d ≈ 150–200 nm; molecule L_p 14.5 nm (Sun 2002); **pore ξ ≈ 2 µm** @2 mg/mL (Mickel 2008 `10.1529/biophysj.108.132738`); **gel G′(2 mg/mL,37°) ≈ 15 Pa**, G′∝c^2.1 (Yang-Kaufman 2009 `10.1016/j.bpj.2009.07.035`; Licup 2015).
- **⚠️α2β1–collagen = SLIP bond only** (no catch data): **k_off0 = 0.44 s⁻¹, x_β = 0.7 nm (f_β = 5.9 pN)** (Attwood 2013 `10.3390/ijms14022832`). Catch shape = α5β1-FN PROXY peak ~30 pN (Kong 2009 `10.1083/jcb.200810002`). **k_on = GAP (PI)**; Chan-Odde default 0.3 s⁻¹ placeholder.
- **Clutch spring κ_c ≈ 0.8 pN/nm** (Chan-Odde 2008 `10.1126/science.1163595` / Bangasser 2013; range 0.5–5). Companion table: F_bond 2 pN, k_off0 0.1/s, motor stall 2 pN, v0 120 nm/s, n_c=n_m=50.
- **Talin unfolds ~5 pN** (R3, del Rio 2009 `10.1126/science.1162912`), step 40–50 nm, vinculin ratchet (Elosegui-Artola 2016 `10.1038/ncb3336`; Yao 2016 `10.1038/ncomms11966`).
- **FA**: traction **5.5 nN/µm² = 5.5 kPa** (Balaban 2001 `10.1038/35074532`, Tan 2003 `10.1073/pnas.0235407100`); mature area ~1 µm²; maturation 2–5 min (Choi 2008 `10.1038/ncb1763`); integrin density ~1000/µm². All fibroblast/CHO **proxy**. MCF7 cell-scale traction ~0.2–0.5 kPa (Gil-Redondo 2023, weak).

## Filopodium (fascin bundle) + lamellipodium (Arp2/3) + polymerization ratchet
- **Fascin bundle**: crossband ~36 nm, axis spacing 8–12 nm (Gong-Alushin 2024 `10.1038/s41594-024-01477-2`); FRAP t½<10 s (Vignjevic 2006 `10.1083/jcb.200603013`); **10–30 filaments/bundle, self-limited ~20, optimum ~30** (Claessens 2008 `10.1073/pnas.0711149105`; Mogilner 2005 `10.1529/biophysj.104.056515`). Bundle EI: fully-coupled ∝N² / decoupled ∝N (Bathe 2008 `10.1529/biophysj.107.119743`) — use coupling model. Single fascin pN/nm **[unverified]** (~1–2 model). Filopodium d 100–300 nm, L up to >35 µm (Mattila 2008 `10.1038/nrm2406`).
- **Filopodial protrusive force ~1–3 pN** (Cojoc 2007 `10.1371/journal.pone.0001072`) — reject "~100 pN" folklore; retraction ceiling ~50 pN (Bornschlögl 2013). Velocity ~17–33 nm/s (Mallavarapu 1999 `10.1083/jcb.146.5.1097`).
- **Arp2/3 branch 70±7°** (Mullins 1998 `10.1073/pnas.95.11.6181`); MD fluctuation ±2–4° (Pfaendtner 2012 `10.1016/j.jmb.2011.12.025`); **debranch ~4e-4 s⁻¹, force-sensitized F½=0.054 pN** (Pandit 2020 `10.1073/pnas.1911183117`; Mahaffy 2006 `10.1529/biophysj.106.080937`); 1 branch/6.2 µm (Vinzenz 2012 `10.1242/jcs.107623`); lamellipodium thickness ~175 nm (Svitkina 1995/1999). Torsional k **[unverified]** → use ±2–4°.
- **Polymerization ratchet**: single-filament stall ~1 pN measured (Footer 2007 `10.1073/pnas.0607052104`) / ~5–9 pN thermodynamic at 10–20 µM G-actin (derived, Mogilner-Oster 2003 `10.1016/S0006-3495(03)74969-8`); barbed-end k_on=11.6 µM⁻¹s⁻¹, k_off=1.4/s, C_c=0.12 µM, δ=2.7 nm, v~310 nm/s @10µM (Pollard 1986 `10.1083/jcb.103.6.2747`); ~100 filaments/µm leading edge (Prass 2006 `10.1083/jcb.200601159`).

## Already-sourced this session (reference)
- **k_erm = 4.6 pN/nm = 4600 pN/µm** (Braunger 2014 `10.1074/jbc.M113.530659`) → KB-3.B1.6 (registered, PI-ratified).

## KB registration plan (PI-gated; proxies/models tagged, not MCF7 measurements)
ECM→KB-1.x, adhesion→KB-2.x, compartments→KB-3.x, cell-type→KB-6.x. ~30 KnowledgeClaims + SourceEvidence.
Hard PI-surface GAPs (never tuned): α2β1–collagen **k_on** + **catch-bond** (α5β1 proxy only); MT **count**;
MCF7 **oblate aspect** + **N:C definition**; IF **cage mesh size**; fascin single-crosslink **stiffness**.
