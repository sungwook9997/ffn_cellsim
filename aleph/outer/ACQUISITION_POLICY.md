# Acquisition and use policy

## Storage classes

1. `open_access`: retain a local object only after the licence and canonical source are recorded. Redistribution follows that licence; Git still stores metadata rather than bulky binaries.
2. `institution_kaist`: access may be used for local research and extraction subject to publisher terms. The payload is `local_only`, excluded from Git, exports, releases, and model bundles.
3. `metadata_only`: retain DOI, bibliographic metadata, abstract where its source permits, and an acquisition request. Do not infer missing results.

Institutional authentication is evidence of access, not permission to redistribute. A URL that downloads successfully does not override this rule.

## Local object store

Payloads reside under ignored `data/external_training/objects/<sha256>`. The manifest records the SHA-256 digest, media type, retrieval timestamp, canonical URL, access route, licence status, and redistribution class. Cookies, tokens, signed URLs, and KAIST credentials are never recorded.

## Extraction contract

An extractable observation must preserve:

- source family and exact figure, table, panel, or supplement locator;
- biological context: species, cell type, state, perturbation, substrate and geometry;
- measurement protocol and calibration;
- observable definition, value, unit, uncertainty, replicate structure, and exclusions;
- whether the value is reported, digitised, derived, or inferred;
- contradiction and quality flags.

Image-derived values retain the transform receipt (crop, segmentation, calibration, registration, and uncertainty). IF/WB/PCR values are not interchangeable with force or geometry measurements; modality-specific observation models must map them into shared latent biological state.

## Dataset integrity

- Split by source family, laboratory, biological system, and experiment series before augmentation.
- Keep the same paper, supplement, replot, and repository deposit in one split.
- Separate parameter priors, training observations, calibration controls, and held-out validation.
- Never promote a record beyond `proposed` merely because it is highly cited or institutionally accessible.
- Preserve failed fits, non-convergence, null effects, and incompatible protocols as negative evidence.
