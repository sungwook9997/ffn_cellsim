# ALEPH-PORT-3635 — importing `ff/ecm_library`, and the one thing an import must not become

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3635` |
| Lane | `46143f30` Lane W2 |
| Status | `PROPOSED` |
| Written | `2026-08-05` — **before the code**, as `CLAUDE.md` §3 requires |
| Port class | **IMPORT, not copy.** No file is transcribed. `aleph/**` calls `ffn_sim.ff.ecm_library` at runtime through a single named seam. |
| Provider | `/Users/sw1/ffn_cellsim/ffn_sim/ff/ecm_library.py` (770 lines), READ ONLY, untouched |
| Aleph target | one new adapter module + `aleph/scenarios/whole_cell.py`'s `_ecm` builder |
| Exists because | The PI recommended importing it rather than re-deriving it (this session's transcript, 2026-08-05 18:54 KST). This lane had recommended a from-scratch port and was asked to reconsider. It reconsidered; the reasons are §2. |

---

## 0. The rule this sits next to, stated rather than skirted

`CLAUDE.md` §3: *"`/Users/sw1/ffn_cellsim` is READ ONLY. **Copying a file is a failure even when
the file is good.** Every port needs a ledger entry in `ports/ledger/` written before the code."*

An import is not a copy, so the first sentence is not what is at issue. **The second sentence is,
and it is honoured here**: this file exists before the code, and it exists so that every ECM number
the engine later reports can be traced to a declaration rather than appearing in the tree unsourced.
That is what the porting rule is protecting, and an import can violate it just as easily as a copy —
by pulling parameters in with no record of where they came from.

`decided_by` is absent, as it must be. The PI's recommendation is cited with its coordinates and can
be checked by opening the transcript; it is not transcribed here as an authorisation.

## 1. What is imported, exactly

Public surface, from the module's own `__all__`:

| symbol | what it is |
|---|---|
| `build_ecm(material, box_lo, box_hi, **kw) -> ECMNetwork` | the dispatcher; fibrillar or continuum by material |
| `REGISTRY` / `get_spec(key)` | six materials: `collagen_I`, `fibrin`, `agarose`, `pa_gel`, `hyaluronic_acid`, `matrigel` |
| `ECMSpec` | the material declaration — literature-sourced microstructure |
| `ECMNetwork` | the built matrix: `.links (M,2)`, `.link_k`, `.link_rest`, `.pinned`, `.net.positions`, and the **emergent** `mesh_size_um`, `connectivity_z`, `S_measured` |

Verified importable from this project's own interpreter with only
`sys.path` extended — `/Users/sw1/miniconda3/envs/aleph/bin/python`, not the `ffn_sim` env. Its
transitive dependencies are `numpy`, `scipy.spatial.cKDTree`, `ffn_sim.ff.units` and
`ffn_sim.ff.fiber_network`, and nothing else.

## 2. Why importing is the better call here, and this lane was wrong to resist it

Three reasons, and the third is the one that changed this lane's position.

1. **The physics is sourced and says so.** The module's own docstring carries the citations
   (Wilhelm & Frey 2003 PRL 91:108103; Head–Levine–MacKintosh 2003 PRE 68:061907; Broedersz &
   MacKintosh 2014 Rev.Mod.Phys. 86:995) and states a hard rule this project would have had to
   invent for itself: *"we never tune the network to hit a target Pa (hard rule:
   no-parameter-tuning-to-outcome)"*. The macroscopic modulus **emerges** from fibre length density,
   connectivity `⟨z⟩`, `κ = k_BT·L_p`, `EA` and `k_xl`. Re-deriving that would produce the same
   equations with fewer citations.
2. **Units already match.** FF is `µm · pN · s` and `1 pN/µm² = 1 Pa` exactly, which is Aleph's
   convention. There is no conversion layer to get wrong.
3. **It takes a box, not a band — so it fixes the defect this lane had just measured.**
   `build_ecm(material, box_lo, box_hi, …)` is constructed in world coordinates and
   `pin_faces=("z_lo",)` already exists. Aleph's ECM today is placed by
   `CYTOPLASMIC_PLACEMENT["ecm"] = (5.2, 7.5)` onto a **radial band midpoint**, which puts it at
   `z = −6.35` while the cell's ventral pole is at `z = −5.00` — a 1.340 µm gap, with adhesions
   spanning 7.251 and 8.611 µm against integrin's 0.020–0.030. A box builder removes the band from
   the load path entirely: the substrate is placed where a substrate goes, because the caller says
   where.

   **This lane argued the opposite an hour earlier** — that a material port could not fix a
   placement defect. That was wrong: it assumed the imported thing would arrive as another
   band-placed compartment. It does not.

## 3. The one thing this import must not become

`aleph/**` gains a runtime dependency on a tree outside the repository. The recorded boundary
principle is *Aleph owns interfaces only; never absorb the provider's dispatch internals, and
project through an allowlist.* Concretely:

- **One seam, named.** Exactly one Aleph module imports `ffn_sim.*`. Everything else in `aleph/**`
  talks to Aleph types. A second import site is the failure mode this section exists to prevent.
- **An allowlist, not a wildcard.** Only the four symbols in §1 cross. `FiberNetwork`,
  `build_fiber_network`, the Watson-distribution helpers and the `_private` functions do not.
- **The dependency is declared, not implicit.** The seam states where the provider is and raises a
  typed refusal naming this ledger when it is absent, rather than failing at import time inside an
  unrelated test.
- **`ECMNetwork` is not an Aleph owner.** It is converted at the seam into whatever Aleph's owner
  contract requires. An imported dataclass appearing in a connector signature is absorption.

## 4. What has to be true before this is called done

| # | control |
|---|---|
| C1 | The seam is the only `ffn_sim` import in `aleph/**`. A grep over the tree is the test. |
| C2 | The cell's ventral pole and the substrate top are in contact to within the adhesion standoff, measured — not asserted from the box arguments. |
| C3 | Adhesion spans fall to the integrin scale. Today: 7.251 / 8.611 µm. The control names the number it expects and fails if the geometry silently reverts. |
| C4 | The emergent modulus the library reports is recorded in `docs/design/AXIS_REGISTRY.md` **with its source**, so this import removes `UNSOURCED` rows rather than adding them. |
| C5 | Every existing whole-cell control still passes, or its change is explained. This is a physics fixture change: committed numbers move. |

## 5. What is NOT claimed

- **Not that this fixes convergence.** `ecm` carries 2,020.78 pN of the 2,622.967 pN residual, but
  the largest free-node residual is `lamellipodium` at 2,456 pN and that is a separate problem.
- **Not that the band guard is obsolete.** `CYTOPLASMIC_PLACEMENT["ecm"] = (5.2, 7.5)` was written
  because four ECM nodes were inside the nucleus. Moving the substrate to a box **replaces** that
  guard and must not merely relax it.
- **Not that the provider's 2D path is exercised.** `dim=2` exists and has not been run here.
- **Not that the import cost is measured.** Build time and memory for a real matrix at cell scale
  are unknown.
