# Second-pass OA recoverability screen

This directory applies the unchanged proposed policy and signal implementation in
`review_oa/screen_oa.py` to the independent second-pass OA extraction. It does not alter the first
pass or automatically promote any article to Tier A.

Full locator/hash records are generated under ignored
`data/external_training/review_oa_second_pass/`. Committed outputs contain compact hashes, counts,
manual-review routing, and within-pass/combined dataset-leakage groups only. No article text or
excerpt is committed.

The combined accounting joins both passes by source-family ID and public dataset-accession group
hash. A shared accession is a split constraint and a manual-review lead, not proof of identical
processed data.
