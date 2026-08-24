# Active Cell NMII transmission and segment-runtime integration

Status: **point-to-segment CUDA mechanism and transaction gates PASS; physiological NMII production remains
blocked by unsourced magnitudes and the unrun full-native active baseline.**

Integration base: `a599e368` on `codex/ff-ac-codex`
Native device: NVIDIA RTX A5000 Laptop GPU, Warp 1.14.0 CUDA

## Outcome

The inherited NMII handoff added angle-harmonic backbone/head-arm mechanics, separated the Hill tangential
load from the Bell full-force load, and supplied barycentric point-to-segment kernels. This integration pass
closed the remaining runtime gap:

- the assembled cell builds an actin-segment topology with one polarity per segment;
- a GPU HashGrid queries live segment midpoints, then exact point-to-segment distance selects the anchor;
- the crossbridge reaction is split over both actin endpoints by the barycentric weights;
- bound-head walk directions follow the current segment geometry and barbed polarity;
- segment attach, Hill step, and Bell detach run only from the scheduler's final accepted outer-step hook;
- rejected attempts do not change motor state or the device RNG epoch;
- the analytic implicit operator and block-Jacobi diagonal include the three-node segment crossbridge tangent.

The node-anchor path remains available only for isolated diagnostics. The composed cell ledger records
`myosin_attachment_primitive = POINT_TO_SEGMENT_BARYCENTRIC`.

## Physiological firewall

The inherited 1.0 µm NMII-backbone persistence length was explicitly marked `GAP → PI` but entered the default
runtime stiffness. It is now a diagnostic convergence-sweep fixture only. The composed runtime has no default:
`nmii_backbone_lp_um` and its provenance are all-or-none, positive, and source-gated before CUDA initialization.
Without them, backbone angular stiffness is absent and the ledger reports
`BLOCKED_UNSOURCED_PHYSICAL_MAGNITUDE`; no full-cell production claim is authorized.

Other pre-existing NMII magnitude gaps remain open. Passing the two-filament mechanism gate does not establish
the MCF7 minifilament density, engaged population, cortical tension, or whole-cell physiological baseline.

## Backbone evidence audit and correction

TAG found no Contract-Graph claim or parameter for mature NMII-minifilament backbone rigidity. Primary-source
review found a 130–170 nm persistence range for a single myosin-II S2/coiled-coil and a nonmuscle-minifilament
assembly model that assumes 130 nm for individual rods. Those are informative proxies, not a direct mature
multi-tail backbone measurement and not MCF7-specific.

An A5000 falsification run at the 0.13 µm single-rod proxy failed NG-1 at **4/6**: transmission 0.9180 passed,
but the two clamp reactions were only 1.7447 and 1.6017 pN versus the 4.9627 pN analytic target. The earlier
broad GAP sweep also contains multiple 5/6 rows. Therefore the historical GAP-insensitivity claim is false;
the 1.0 µm 6/6 fixture cannot be promoted to a physiological motor verdict. No value or gate was tuned.
See `NMII_BACKBONE_EVIDENCE_AUDIT.md` and `ng1_lp013_proxy_a5000.log`.

## Transaction and RNG correction

The earlier host commit counter advanced even when an outer attempt was rejected. Although bond state stayed
unchanged, the next accepted attempt would receive a different random stream. NMII and ERM now own
device-resident accepted epochs: rejection preserves both biological state and epoch bit-exactly; acceptance
advances each epoch once after its KMC launch. Post-loop ledgers expose the accepted epochs and bound counts.

## Verification

- Local full `ffn_sim/tests/ac/`: **390 PASS + 30 CUDA-only skip** (420 collected).
- RTX A5000 full suite, split by package to avoid one-process CUDA module accumulation: **420/420 PASS**.
- RTX A5000 focused segment/ERM/implicit/wiring suite: **35/35 PASS**.
- RTX A5000 post-ledger cell/root regression: **100/100 PASS**.
- RTX A5000 0.13 µm single-coiled-coil proxy falsification: **NG-1 FAIL 4/6** (expected finding; production
  remains blocked).
- RTX A5000 composed 100-filament wiring probe: 600 actin segments, 20 heads,
  `POINT_TO_SEGMENT_BARYCENTRIC`, and the unsourced-backbone firewall present.
- Changed-file Ruff: **PASS**.

The focused CUDA runtime gate proves: rejected query/KMC is a bit-exact no-op; accepted high-hazard attachment
selects the segment midpoint (`t = 0.5`), assigns both real endpoints, follows the barbed direction, advances
the power stroke, and increments the accepted RNG epoch exactly once. Test hazards are controls, not
physiological parameters.

## Figures

- `figs/round2_verdict.png` — two-filament transmission and force diagnostics against the active cortical band,
  with the unsourced persistence-length fixture and unresolved full-cell magnitude stated explicitly.

## Open items

1. Register and PI-ratify NMII backbone persistence length (or an equivalent bending-rigidity contract).
2. Resolve the remaining NMII kinetic/stall/density GAPs before enabling a physiological active baseline.
3. Run the full 70,686-cortical-filament accepted active trajectory and method-of-planes cortical-tension
   measurement; the current resting candidate remains unconverged and NG-3 remains open.
4. Administrator-enable NVML accounting and rerun NG-9 in a fresh process.
5. Build and ledger the separately counted I5+ physiological populations.
