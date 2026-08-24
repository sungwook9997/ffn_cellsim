# ECM GPU topology handoff

## Frozen scope

This slice was built from `6f91aec8` on `codex/ecm-gpu-topology`.  It does not
modify `ac/cell`, the coupled solver, convergence, seams, connector scheduling,
or `ac/engine/ecm_world.py`.  It supplies a concrete ECM-owned CUDA topology
behind that already-ratified engine contract.

The collagen constitutive modulus gate remains **HELD**: the existing FF sources
contain conflicting reference bands of 30–100 Pa and 5–100 Pa.  This package
does not import either value, does not implement constitutive mechanics, and
cannot be used for a material-production closeout.

## Implemented ownership

`ECMTopologyState` owns fixed-capacity Warp-CUDA arrays for:

- current/reference node position and force;
- active node/fiber/segment/bend flags and explicit ownership indices;
- ordered fiber allocation offsets, segment pairs, and bending triples;
- persistent fiber/node/segment identities, slot generations, parent-segment
  lineage plus dyadic root/path lineage, topology epoch/dirty state, and a
  device-side next-ID cursor;
- segment material coordinates/rest lengths, refinement level, damage, and
  graph-endpoint refcounts;
- boundary-face tags only (the boundary connector still owns reaction/work);
- active-population counts, refinement-derived dormant capacities, and device
  free stacks for node/segment/bend slots.

`max_refinement_level` is a required resolution contract.  Initial biological
population is never reduced: reserve capacity is derived as binary segment
splits through that level and is tracked separately as dormant allocation.

## GPU Mikado initialization

`MikadoTopologyBuilder.initialize()` performs the following entirely on CUDA:

1. samples uniform centers and isotropic rod directions from a deterministic
   per-fiber Warp RNG stream;
2. clips each straight rod against the 3-D box without bending it at a face;
3. writes the ordered nodes, consecutive axial segments, and consecutive bend
   triples into per-fiber fixed-capacity blocks;
4. assigns persistent identities and zero generations;
5. tags prescribed far-field face bands; and
6. packs every dormant refinement slot into device free stacks.

It does not create static crosslink bonds.  Crosslink chemistry and accepted-step
association remain graph-owned.

## Accepted-step remodeling transaction

`ECMRemodelRuntime` implements ECM-internal remodeling without taking ownership
of a force law or damage criterion.  An external CUDA owner writes only:

- per-segment `NONE / REFINE / COARSEN` proposals;
- nonnegative accepted-step damage increments; and
- per-fiber `NONE / SLEEP / WAKE` proposals.

The execution order is:

1. `snapshot_candidate()` clears stale proposal/plan state.  It does not copy the
   full ECM because authoritative topology is commit-only and remains its own
   accepted snapshot.
2. The external physics/KMC owner writes proposal arrays on CUDA.
3. `prepare_candidate()` validates inputs, greedily selects edits with a
   one-segment exclusion halo, performs device prefix scans for deterministic
   fresh IDs, assigns dormant resources inside each owning fiber's fixed block,
   and publishes an endpoint-remap table.
4. The graph-owned connector remaps all live endpoints and writes acknowledgement
   plus target-refcount arrays.  `validate_remap_acknowledgements()` publishes the
   device `candidate_valid_d` scalar.
5. The global scheduler calls both transaction hooks.  Rejection only clears the
   proposal.  `commit_irreversible()` changes damage/topology/sleep state only when
   both `accepted_d[0]` and `candidate_valid_d[0]` are nonzero.

No accepted predicate or validity scalar is read on the host.  Refinement splits
the current/reference geometry, rest length, material interval, damage state, and
endpoint population.  Coarsening performs the inverse chain/bend edit and uses a
contour-length-weighted damage average.  Every created node/segment gets a fresh
persistent ID; slot generations change on reuse or retirement.  Root/path lineage
preserves sibling identity through repeated multilevel refine/coarsen cycles.

### Connector endpoint remap

For a source material point `u ∈ [0,1]`, the connector graph applies:

- `SPLIT`: left target with `u' = 2u` for `u ≤ 0.5`; otherwise right target with
  `u' = 2u - 1`;
- `MERGE_LEFT`: merged target with `u' = 0.5u`;
- `MERGE_RIGHT`: merged target with `u' = 0.5 + 0.5u`.

The table includes source and planned target persistent IDs and generations.
Commit is disabled unless every source acknowledgement equals its pre-edit
`endpoint_refcount` and split/merge target counts conserve that population.  This
package deliberately does not mutate connector chemistry or FA/clutch state.

An accepted topology or sleep/wake edit increments `topology_epoch_d` and sets
`topology_dirty_d`.  That scalar is the graph-owned dynamic broadphase handoff:
the connector refreshes its fixed-capacity candidate workspace, then calls
`mark_candidate_query_refreshed()` to clear the dirty bit.  ECM never reaches
into the Lead-owned connector graph.

## Material-point capture candidates

`initial_contact_candidates()` builds a float64 Warp `HashGrid` over packed
active segment midpoints.  The conservative query radius is

`max_segment_length + crosslink_capture_radius`.

Every culled pair is then tested by exact segment–segment closest-point geometry.
Each emitted candidate carries both segment indices, both persistent segment
IDs, both generations, both dimensionless material coordinates `u_a/u_b`, both
fiber indices, and separation.  Same-fiber pairs are excluded and the canonical
pair ordering emits each segment pair once.

The prefix count/scan/fill is CUDA-resident.  One initialization-only scalar
read sizes the exact candidate array before physical time begins.  This method
is not a hot-loop refresh path; dynamic accepted-step crosslink search must use
a graph-owned fixed-capacity workspace.

## Gates run

- Local macOS authoring environment: all `aleph/tests/ac` tests passed; CUDA
  gates skipped under the I0-A rule.
- Isolated gbook directory, RTX A5000 Laptop GPU: initialization/contact and
  remodeling native gates passed (22/22).
- The CUDA candidate set matched a brute-force all-active-segment oracle exactly,
  including pair keys, both material coordinates, and separation.
- Same-seed initialization produced bit-identical geometry, topology, IDs,
  boundary tags, and free stacks.
- Rejected and invalid-ack transactions left topology bit-identical; accepted
  refine/coarsen conserved material interval, rest length, endpoint population,
  active/free counts, and chain/bend connectivity.
- Two refinement levels coarsened back to the original root/path/rest state, with
  monotone non-reused persistent IDs.
- Static gates reject NumPy/SciPy/cKDTree, legacy FF ECM builders, `ac/cell`, and
  `ecm_world` dependencies from the runtime package.
- The full A5000 `tests/ac` suite passed with one explicitly deselected,
  ECM-unrelated load-path force assertion.  That assertion fails identically at
  the untouched B1 commit `0ad3b40d`; no engine/load-path file was changed here.

## Deferred to the Lead-owned coupled core

- collagen stretch/bend/damage constitutive laws and the held modulus decision;
- crosslink selection, one-per-fiber-pair policy if ratified, and Bell/KMC state;
- far-field reaction/work mechanics;
- FA-clutch and membrane-contact connector binding;
- graph consumption of the endpoint remap and global acceptance-ledger wiring;
- graph-owned dynamic broadphase workspace and native physiological-population
  memory/peak ledger closeout.
