# Visual-observation lane report

Status: **complete, AGENT-PROPOSED, non-authoritative**  
Run date: 2026-08-05 KST  
Frozen annotation selection: 300 sources, SHA-256
`f4be334ea3fce6659a7395e77a2107a3b6a3dfdb5547f8b163e535946ff83226`

## Outcome

The two OA JATS stores contain exactly 38,316 figures from 5,635 source families and 57,957
`<graphic>` hrefs. The lane emitted one deterministic visual observation candidate per figure,
without caption text, and compressed the full projection to 6.7 MiB. Every candidate retains the
source XML hash, exact JATS figure/graphic locator, caption hash, proposed modality tags and a typed
refusal to derive scientific values from uncalibrated pixels.

The extraction lane's final figure-join index was bound at SHA-256
`906f1c9e10ea1489104fb3600fca4616833a69d8f060010518d61205be16a65f`.
Exact source-family + normalized JATS figure ID + caption SHA linked **14,302 / 14,302** extraction
rows to canonical `aleph:experiment:*` identifiers, with **0 caption-hash mismatches**. The other
24,014 figures remain explicitly typed figure-scoped candidates; they were not silently attached to
nearby text experiments.

## Frozen 300-source acquisition

All 300 selected sources have canonical receipts. The denominator hierarchy matters:

| Denominator | Count |
|---|---:|
| selected source families | 300 |
| selected figures | 2,250 |
| selected graphic identities, including alternate renderings | 3,315 |
| canonical receipt rows / unique graphic IDs | 3,315 / 3,315 |
| duplicate receipt IDs | 0 |

Only official `pmc.ncbi.nlm.nih.gov` article pages were used to discover assets, and only URLs that
the page itself listed under `https://cdn.ncbi.nlm.nih.gov/pmc/blobs/` were fetched. There were zero
official-URL violations and no publisher or non-OA fallback.

| Receipt state | Graphic references |
|---|---:|
| acquired | 627 |
| official filename unresolved at bounded closeout | 368 |
| deterministically deferred by figure budget | 2,320 |

Acquisition produced 480 unique SHA-256 objects from 116 / 300 sources. Acquired coverage is
480 / 2,250 selected figures (21.33%); acquired-reference coverage is 627 / 3,315 (18.91%). The
verified downloaded payload total is 174,796,895 bytes and content-addressed local storage is about
118 MiB. All object hashes revalidated with zero failures.

The bounded continuation selected caption-tagged figures first and allowed at most two newly
requested figures per source. The earliest resumable phase began before that cap and its already
downloaded, hash-verified assets were retained rather than deleted. Therefore the final receipt
distribution is observed provenance, not a claim that every source has exactly two images. The
2,320 deferred identities make the remaining lazy-acquisition scope explicit.

## Modality coverage

Caption rules are proposals, can overlap, and do not constitute assay truth.

| Proposed modality | All figures | Frozen-set figures | Acquired references |
|---|---:|---:|---:|
| immunofluorescence | 4,480 | 404 | 150 |
| fluorescence microscopy | 3,029 | 289 | 98 |
| western blot / gel | 5,343 | 393 | 143 |
| PCR / qPCR curve | 2,864 | 218 | 87 |
| morphology | 4,070 | 334 | 107 |
| time series | 2,705 | 238 | 71 |
| AFM | 1,482 | 96 | 27 |
| TFM map | 250 | 24 | 13 |
| PIV / vector field | 238 | 15 | 6 |

## CPU descriptors and panels

FFmpeg/ffprobe 8.1.1 decoded a fixed, aspect-preserving 128×128 RGB view; NumPy 2.5.1 computed
source dimensions/pixel format/frame count, channel spread, grayscale mean/std/quantiles, entropy,
edge density, and bright-gutter panel candidates. This took about 13 seconds with three CPU workers
for 480 unique objects; duplicate references reuse the same descriptor.

- 627 receipt rows described, representing 480 unique images.
- Source dimensions span 380–1,985 px wide and 182–2,100 px high.
- 480 unique images comprise 442 `yuvj444p`, 27 `yuvj420p`, and 11 grayscale decodes.
- 273 unique images have panel candidates, totaling 769 candidate boxes. Receipt-weighted counts are
  354 and 987 because alternate graphic identities may share an asset.
- OCR is typed `unavailable_no_local_ocr_backend`.
- Pixel-derived physical or biological values: **0**; calibration-dependent extraction is refused.

Descriptor generation was replayed independently. The gzip bytes were identical with SHA-256
`52e25fe8d7ea7e9d7c36b0b5986022e5c16ef153e85742204bcae7fc6de1f6cd`.

## Runtime and controls

Observed CPU wall times were approximately 25 seconds for the full JATS figure index, 3 seconds for
the final experiment join, and 13 seconds per full descriptor pass. The acquisition window spanned
about 27 minutes including bounded retries, checkpoint restarts, and the policy-cap transition; it
is an operational receipt, not a throughput benchmark.

Focused controls cover deterministic/domain-separated IDs, proposed modality classification,
official-host allowlisting, extension fallback without publisher fallback, bright-gutter panels,
uniform/border adversaries, byte-exact hashes, missing graphics, no committed caption text,
deterministic descriptor replay, and exact figure-locator parsing.

Final audit invariants: 0 duplicate index IDs; 0 duplicate receipt IDs; receipt/descriptor ID sets
equal; 0 object digest failures; 0 official URL violations; 0 non-proposed records; 0 scientific
values extracted from pixels.

## Limits

These are observation candidates and layout descriptors, not scientific annotations. Panel boxes
need human review; channels are decoded color channels rather than fluorophore identities; animation
metadata is not a time calibration; captions can mention an assay without depicting it. No asset,
tag, link or descriptor promotes Aleph evidence or runtime authority.
