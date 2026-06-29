# ffn_cellsim — developer entrypoints.
#
# The KB integrity gates are pure-Python (stdlib + pyyaml) and read only
# COMMITTED artifacts (config YAMLs, *_manifest.yaml, source_audit_report.md) —
# no Notion token, no network — so they run anywhere in well under a second.
# CI runs the blocking gates on every push (.github/workflows/ci.yml); this
# Makefile is the single local entrypoint so CI, the git hook, and the docs
# never drift apart.

TKB := ffn_sim/outputs/tag_kb

.PHONY: help kb-check kb-figs hooks warp-parity

help:
	@echo "make kb-check  - run the KB integrity gates (results + params blocking; citation audit informational)"
	@echo "make kb-figs   - regenerate all presentation artifacts (figures + 5-hop demo + standalone HTML)"
	@echo "make hooks     - enable the tracked git pre-commit hook (runs kb-check before each commit)"
	@echo "make warp-parity - Warp DCM engine parity tests, CPU backend (GPU verdict: parity_report.py --device cuda:0 on the A5000)"

# All KB auditors. results + params are BLOCKING gates (exit 1 on drift = a
# result/constant claimed better than the disk supports). The citation audit is
# read-only here; the production-cited subset of sources is already enforced via
# the params gate (a fabrication-risk citation on a constant declared 'verified'
# drifts -> blocks).
kb-check:
	@python $(TKB)/verify_runs.py --gate
	@python $(TKB)/verify_params.py --gate
	@python $(TKB)/verify_sources.py --check

# Phase-C Warp DCM engine parity (CPU backend). The GPU-backend verdict is produced
# on the gbook A5000 via parity_report.py --device cuda:0 (see warp_port/ENGINE.md).
warp-parity:
	@python -m pytest ffn_sim/tests/dcm -q

kb-figs:
	@cd $(TKB)/presentation && for f in fig*_*.py; do echo "  render $$f"; python "$$f" >/dev/null; done
	@python $(TKB)/presentation/demo_5hop_query.py >/dev/null && echo "  demo  -> DEMO_5hop_query.md"
	@python $(TKB)/presentation/export_html.py
	@echo "presentation artifacts -> $(TKB)/presentation/"

# Opt-in: route git hooks to the tracked .githooks/ dir so the gates also run
# locally before each commit. Reversible: git config --unset core.hooksPath
hooks:
	@git config core.hooksPath .githooks
	@echo "pre-commit hook enabled (core.hooksPath=.githooks)."
