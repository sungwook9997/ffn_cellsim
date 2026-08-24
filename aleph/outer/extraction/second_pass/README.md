# Second-pass OA JATS extraction

This directory records the proposed, digest-only audit of the independently
discovered second-pass OA corpus.  It reuses
`corpus/external_training/extraction/extract_jats.py` without modification and
writes text-bearing chunks only beneath ignored
`data/external_training/derived_second_pass/`.

No full text, long excerpt, figure, table, publisher payload, or authoritative
claim is committed here.  Lexical tags remain non-authoritative retrieval aids.

`audit_second_pass.py` compares two extraction summaries while ignoring only
the expected run-history field `resumed_articles`, then independently verifies
every local chunk-file digest, text digest, chunk digest, and authority marker.

