# Visual observation candidates

This proposed lane turns OA JATS `<fig>` elements into deterministic, text-free observation
records. It links source family, a figure-scoped experiment candidate, exact XML locator, caption
digest, graphic href digest, proposed modality tags, official-OA asset receipt, and reproducible CPU
image descriptors.

The complete text-free JSONL index is deterministically gzip-compressed (`mtime=0`) to keep the
38,316-record committed projection compact; its SHA-256 is recorded in the summary.

It does **not** read scientific values from pixels. Panel boxes are layout candidates, modality tags
are caption-rule proposals, and every record remains below Aleph runtime or evidence authority.
Assets are acquired only for the schema lane's 300-source `gold_candidate` annotation selection,
only through the official PMC article page and NCBI CDN URLs, and remain under ignored `data/`.
The committed receipts cover every selected graphic identity once. Network work is bounded: tagged
figures are preferred, newly requested continuation work is capped at two figures per source, and
the rest remains an explicit deferred state for later annotation-time acquisition.

The descriptor is deliberately dependency-light: `ffmpeg` decodes a fixed 128×128 RGB view and
NumPy computes dimensions, intensity quantiles, entropy, edge density, color spread, and bright
gutter panel candidates. OCR is typed unavailable because no local OCR backend was present.

Final measured results and denominators are in [REPORT.md](REPORT.md); machine-readable counts and
hashes are under `results/`.
