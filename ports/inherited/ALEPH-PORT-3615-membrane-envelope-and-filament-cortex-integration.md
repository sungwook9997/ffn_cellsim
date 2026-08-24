# ALEPH-PORT-3615 — membrane osmotic envelope and explicit filament-cortex integration

| Field | Value |
|---|---|
| Aleph port ID | `ALEPH-PORT-3615` |
| Lane | `L36 integration — vertical ownership and independent topology` |
| Status | `PROPOSED` |
| Written | `2026-08-03` — before the code, per `PLAN.md` §0.2.5 |
| Port class | `RE-DERIVED` |
| Authorises | `aleph/vertical/cortex_surface_coupling.py`, the corresponding builder changes in `aleph/vertical/assembly.py`, and the osmotic-owner declaration in `aleph/state/census_environment_surface.py` |

---

## 1. Aleph API

This entry authorises the following exact surface. The existing constitutive laws in
`cortex_filaments.py`, `membrane.py`, and `pressure.py` are composed, not reimplemented.

```python
from aleph.vertical.cortex_surface_coupling import (
    CortexSurfaceContact,
    ErmSurfaceTether,
    build_radial_cortex_network,
    radial_surface_mapping,
)
```

It also authorises `build_vertical(...)` to return a `CortexFilamentNetwork`, apply
`OsmoticEnvelope` directly to `HelfrichMembrane`, and bind the two independent-topology couplings.

## 2. Source identity

| Field | Value |
|---|---|
| Source repository | `/Users/sw1/ffn_cellsim` (READ ONLY; never modified) |
| Source commit | `be0e58760baaddc04460bb6b34b0476b5d8aa6c5` (`be0e5876`) |
| Source path | `ffn_sim/ac/engine/cortex_population.py`; `ffn_sim/ac/engine/cortex_state.py`; `ffn_sim/ac/solid/gamma_floor.py`; `ffn_sim/ff/forces_warp.py`; `ffn_sim/ff/constraints.py` |
| Source symbols | architectural roles of the cortex population, membrane surface, ERM/contact coupling, and the reference axial-constraint split; no source symbol is ported |
| Read from | working tree, during a bounded read-only audit explicitly requested by the user |
| Working tree == commit? | no — the source tree was dirty; the immutable commit is cited as its nearest identity, and no claim is made that every read byte equals it |
| Lines or constants taken | **0 lines; 0 identifiers; 0 parameters** |

## 3. Why source-derived porting beats clean-room

It does not. This is a clean-room, re-derived integration. The source is useful only as prior-art
evidence for the architectural split: a closed membrane can own pressure while a distinct filament
graph owns cortical mechanics. Its implementation is not suitable for direct porting: the audited
tree splits axial inextensibility out of the cortex accumulator, its legacy floor module labels
under-connected crosslinks and disabled overlap resolution, and its working tree differs from the
cited commit.

The Aleph implementation therefore composes its already-controlled axial/bending/crosslink cortex,
Helfrich membrane, and osmotic volume gradient. The new code is only the geometry map and
equal-and-opposite scatter between disjoint topologies.

## 4. Physical or mathematical law represented

There are three laws and one ownership rule.

1. **Osmotic load:** for membrane vertex `i`, `F_i = delta_p * grad_i(V)`, where the oriented closed
   membrane owns `V`. Pressure is never accumulated on the filament graph.
2. **ERM tether:** a membrane site and its mapped cortical material point interact through a
   tensile-only spring `E = 1/2 k max(d - l0, 0)^2` by default. Forces are the negative gradient and
   scatter to the two owners with exactly opposite resultants.
3. **Containment contact:** if the signed inward clearance of a cortical node relative to its mapped
   membrane site is below the declared clearance, a one-sided quadratic penalty acts along the
   membrane normal. The inactive branch is exactly zero.
4. **Ownership:** `membrane` is the sole osmotic envelope; `cortex` is an explicit filament graph and
   is forbidden from also owning a volume load.

The map is deterministic but not constitutive: membrane vertices map to named cortex material
points by radial direction, and all scatters use that immutable index map. No equal-array-length or
same-index assumption is allowed.

## 5. Units, domains, singular cases, invariants

| Quantity | Aleph unit | SI unit | Domain |
|---|---|---|---|
| position, rest length, clearance | µm | m | finite; rest length and clearance non-negative |
| force | pN | N | finite |
| energy | pN·µm | J | finite and non-negative for the two couplings |
| pressure | pN/µm² | Pa | finite |
| spring stiffness | pN/µm | N/m | finite and non-negative |

Singular cases are refused: a zero-length active pair has no defined direction; an open or
mis-oriented pressure surface has no valid signed volume; an out-of-range map cannot scatter; and a
stale cortex material-point handle cannot be reused after a topology transaction. Contact outside
its active branch returns exact zeros, not a small residual.

Invariants:

- **I1.** Membrane and cortex arrays have different sizes, and coupling still closes force and moment.
- **I2.** Pressure is accumulated on membrane only; cortex receives load only through connectors.
- **I3.** The assembled cortex is `CortexFilamentNetwork`, has a bound axial law, and has no triangles
  or volume method masquerading as a surface.
- **I4.** Every connector evaluation returns full-owner force blocks and equal-and-opposite pair
  resultants.
- **I5.** A uniformly contracted straight filament costs axial energy, preserving the control earned
  by `ALEPH-PORT-2301`.

## 6. Source evidence class and known retractions

Source evidence class is `AUDIT_READ`: code and self-labelled defect records were read, but the
source was not executed and no GPU job ran. The source contains contradictory live descriptions of
where pressure enters and separates axial constraint enforcement from the cortex force accumulator.
Those are known limitations, not Aleph authority. The legacy floor module additionally marks
under-connected crosslinks, disabled overlap resolution, and a mono-length population as live
structural limitations. None of its parameters or prose survives this re-derivation.

The Aleph decision provenance is stronger and local: the PI selected the filament-graph cortex on
2026-07-30 and then directed that osmotic ownership move to the membrane. The user reconfirmed this
implementation after the present read-only audit on 2026-08-03.

## 7. Independent oracle or derivation

The pressure oracle is the discrete volume-gradient identity and the spherical Laplace balance
`delta_p = 2 sigma / R`, measured by recovering `sigma` from the assembled membrane's external
virial. The connector oracle is exact translational and rotational closure: pair scatter must have
zero net force and zero net moment about any origin. Both are independent of the source tree.

The topology oracle is structural: changing the cortex node count while keeping the membrane mesh
fixed must leave the map valid and the pressure field shape equal to the membrane shape. A coupling
that relies on corresponding raw array indices fails this check before dynamics begin.

## 8. Positive control

| Control | Establishes |
|---|---|
| `tests/vertical/test_cortex_surface_integration.py::test_assembled_cortex_is_an_explicit_graph_with_independent_topology` | The assembled cortex is a 3-node-per-site filament graph with crosslinks, axial mechanics, and no triangles or volume API. |
| `tests/vertical/test_cortex_surface_integration.py::test_pressure_is_owned_by_membrane_and_has_membrane_shape_only` | The sole pressure field has membrane shape and targets the membrane. |
| `tests/vertical/test_cortex_surface_integration.py::test_mapped_connectors_close_force_and_moment_across_different_node_counts` | Both mapped couplings close resultant force and moment across unequal owner-array sizes. |
| `tests/vertical/test_cortex_surface_integration.py::test_signed_normal_contact_force_is_the_gradient_of_its_energy` | Signed-normal contact includes the vertex-normal derivative and matches the directional energy gradient. |
| `tests/vertical/test_vertical_controls.py::TestLaplace::test_recovered_tension_matches_the_declared_tension` | The membrane-owned pressure path retains the independent spherical Laplace oracle. |

## 9. Deliberately failing negative control

| Control | Deliberately planted failure |
|---|---|
| `tests/vertical/test_cortex_surface_integration.py::test_pressure_routed_to_filament_cortex_is_refused_structurally` | Routes the osmotic field to the graph and requires an assembly-time refusal. |
| `tests/vertical/test_cortex_surface_integration.py::test_same_index_shortcut_is_rejected_for_the_real_owner_arrays` | Hands unequal full-owner arrays to the same-index primitive and requires a shape refusal. |
| `tests/vertical/test_cortex_surface_integration.py::test_stale_material_map_is_refused_after_topology_commit` | Commits a graph topology change and requires the old map epoch to be rejected. |
| `tests/vertical/test_cortex_surface_integration.py::test_without_connectors_external_membrane_load_is_exactly_pressure` | Removes both membrane-cortex connectors and proves no hidden pressure relay remains. |

## 10. Numerical and precision envelope

All host mechanics use NumPy `float64`. Geometry maps are integer arrays and must match exactly.
Force/moment closure is scaled by constituent magnitudes, not the near-zero resultant, and is
required below `1e-13` relative. The signed-contact directional gradient is required below `2e-7`
relative on its deterministic active fixture. Existing Laplace and mesh-family thresholds remain
unchanged. Domain errors are raised rather than hidden by clamps.

## 11. Production-backend residency and transfer

This entry is host-side CPU integration. It performs no GPU execution and does not claim CUDA
residency. The existing resident-world demonstration remains a separate approximation and is not
evidence for this topology. A future device port must carry the immutable index map and material
point interpolation weights to device without a per-step host round-trip; that work is outside this
entry.

## 12. Comments and docstrings to discard

Discard all source package names, internal gate labels, working-tree assumptions, CUDA kernel
details, and claims that the source's external projector makes its owner independently complete.
Discard the legacy floor's parameter values and defect-status prose. New prose names only Aleph
owners, laws, maps, domains, and controls, so the package remains readable without the source tree.

## 13. Acceptance

| Field | Value |
|---|---|
| Acceptance result | Implemented on 2026-08-03. All controls in §§8–9, the full indentation and spreading scenario files, affected vertical/state/CellWorld controls, and new-file Ruff checks pass. Status remains `PROPOSED` pending PI ratification. |
| Reviewer | Codex implementation proposed under the user's direct instruction; final status remains unratified unless the PI ratifies it |
| Rollback | Revert the new coupling module and builder/census/audit edits together; the prior placeholder shell world returns, while `CortexFilamentNetwork` itself remains available under `ALEPH-PORT-2301` |

## 14. Honest limits

- The radial map is a deterministic coarse-grained attachment model, not measured ERM density or
  cortical remodelling kinetics.
- The builder does not claim a biological filament-length distribution, branching population, or
  crosslink kinetics.
- No GPU path, large-cell scaling result, or remeshing-after-topology-change result is established.
- Moving material-point attachments across a committed cortex topology change require rebuilding the
  map; automatic rebinding is outside this entry.
- Passing a spherical Laplace oracle does not validate non-spherical biological cortex parameters.
