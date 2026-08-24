# ffn_cellsim — developer entrypoints.
#
# The KB integrity gates are pure-Python (stdlib + pyyaml) and read only
# COMMITTED artifacts (config YAMLs, *_manifest.yaml, source_audit_report.md) —
# no Notion token, no network — so they run anywhere in well under a second.
# CI runs the blocking gates on every push (.github/workflows/ci.yml); this
# Makefile is the single local entrypoint so CI, the git hook, and the docs
# never drift apart.

# `python` is NOT on PATH in this project's environments — the conda env exposes its own
# bin dir and there is no bare `python` shim. Every recipe below used it anyway, so
# `make kb-check` — the gate the closeout protocol requires — died with
# "make: python: No such file or directory" before a single check ran.
#
# The predicate is deliberately "can this interpreter RUN the gates", not "does a python
# exist": the first attempt at this fix resolved /usr/bin/python3 and traded a PATH error
# for a ModuleNotFoundError, which is the same failure wearing a different hat.
# CONDA_PREFIX is the real mechanism; the ~/miniconda3 entry is a last-resort convenience.
# Override explicitly with `make PY=/path/to/python <target>`.
PY := $(shell \
  for c in "$(CONDA_PREFIX)/bin/python" "$$HOME/miniconda3/envs/ffn_sim/bin/python" \
           "$$(command -v python 2>/dev/null)" "$$(command -v python3 2>/dev/null)"; do \
    [ -x "$$c" ] && "$$c" -c 'import yaml' >/dev/null 2>&1 && { echo "$$c"; break; }; \
  done)
PY := $(if $(PY),$(PY),python3)

TKB := aleph/outputs/tag_kb

.PHONY: help boot-check kb-check kb-figs hooks state state-check gpu dashboard own cell cell-png

help:
	@echo "make boot-check - budget + no-magnitudes + memory index + retired-claim cross-check"
	@echo "make kb-check  - run the KB integrity gates (results + params blocking; citation audit informational)"
	@echo "make kb-figs   - regenerate all presentation artifacts (figures + 5-hop demo + standalone HTML)"
	@echo "make hooks     - enable the tracked git pre-commit hook (runs kb-check before each commit)"
	@echo "make state      - render STATE.md tier-(a) from the append-only aleph/docs/state_rows.yaml"
	@echo "make state-check- fail if that table is stale (the pre-commit hook runs this)"
	@echo "make gpu        - how to get a card, and why there is no lease to check"
	@echo "make own        - may THIS session (FFN_SESSION) commit what is staged?"
	@echo "make dashboard  - one page: lease, lanes, evidence axes, runs, doc debt"
	@echo "make cell       - open the newest cell export in the native app (no size limit, nothing thinned)"
	@echo "make cell-png   - the same, cut open, to a png -- for a machine with no display"
	@echo "make tau3       - read the tau-III seeds in the declared order, then table, then figures"

# All KB auditors. results + params are BLOCKING gates (exit 1 on drift = a
# result/constant claimed better than the disk supports). The citation audit is
# read-only here; the production-cited subset of sources is already enforced via
# the params gate (a fabrication-risk citation on a constant declared 'verified'
# drifts -> blocks).
boot-check:
	@$(PY) aleph/scripts/check_boot.py

kb-check:
	@$(PY) -c 'import yaml' 2>/dev/null || { echo "  kb-check needs an env with pyyaml — 'conda activate ffn_sim', or 'make PY=/path/to/python kb-check'"; exit 1; }
	@$(PY) $(TKB)/verify_runs.py --gate
	@$(PY) $(TKB)/verify_params.py --gate
	@$(PY) $(TKB)/verify_sources.py --check

# `warp-parity` was REMOVED 2026-07-29. It ran `aleph/tests/dcm`, deleted with the
# HOOMD tree that produced its fixtures. Per STATE.md (c) 13 the 25 `warp-*` parity
# claims are permanently uncheckable — regenerating them needs the forbidden runtime —
# so there is no target to restore, and a target pointing at an absent directory reads
# as a parity leg that merely has not been run.

kb-figs:
	@cd $(TKB)/presentation && for f in fig*_*.py; do echo "  render $$f"; $(PY) "$$f" >/dev/null; done
	@$(PY) $(TKB)/presentation/demo_5hop_query.py >/dev/null && echo "  demo  -> DEMO_5hop_query.md"
	@$(PY) $(TKB)/presentation/export_html.py
	@echo "presentation artifacts -> $(TKB)/presentation/"

# Opt-in: route git hooks to the tracked .githooks/ dir so the gates also run
# locally before each commit. Reversible: git config --unset core.hooksPath
hooks:
	@git config core.hooksPath .githooks
	@echo "pre-commit hook enabled (core.hooksPath=.githooks)."

# --- parallel-session coordination -------------------------------------------
# One tree, several sessions, one GPU. Ownership is declared (.claude/ownership.yaml)
# rather than enforced by a second worktree, and STATE.md's tier-(a) table is
# GENERATED so concurrent appends land in a YAML list instead of on one markdown line.

state:
	@$(PY) aleph/scripts/render_state_rows.py

state-check:
	@$(PY) aleph/scripts/render_state_rows.py --check

# ⚠ There is nothing to query. Until 2026-08-04 this printed who held the A5000's lease; the PI
#   retired that card and decision B2 archived the lease, because it was a second lock bolted to the
#   outside of a door the kernel already holds shut -- outside a Slurm allocation cuInit() returns
#   CUDA_ERROR_NO_DEVICE. The help line kept saying "who holds the one A5000" for seventeen days
#   after the card stopped existing.
gpu:
	@echo "no lease to check: it was archived with the A5000 (decision B2, PI 2026-08-04)."
	@echo "the device is the shared workstation sungwook@100.110.26.26 -- two RTX 4090, one RTX 3090."
	@echo "ASK THE PI for a card and a duration, EVERY TIME, then:"
	@echo "    MEM_GB=<GB> CPUS_PER_TASK=<n> gpu-submit <card> <time> <command>"
	@echo "an agent may TRANSCRIBE the PI's grant with transcript path and timestamp."
	@echo "an agent may NEVER write one."
	@echo
	@echo "what is running right now (needs the host to be reachable):"
	@-ssh -o BatchMode=yes -o ConnectTimeout=8 sungwook@100.110.26.26 \
	    'squeue -o "%.5i %.9u %.8T %.11L %.20R %j" | head -12' 2>/dev/null \
	    || echo "  (host not reachable from here -- that is not evidence the card is free)"

own:
	@$(PY) aleph/scripts/check_ownership.py

dashboard:
	@$(PY) aleph/scripts/ffn_dashboard.py

# ── the native cell viewer ────────────────────────────────────────────────────────────────────────
# `make cell` opens the most recent export. No path, no flags, no size limit — which is the whole
# reason it is not an HTML page: a 4.59 M-node cell is 94 MB static and 258 MB for four frames, and
# nothing is thinned anywhere in the tool. Keys are printed on start; `c` cuts the cell open, which is
# how the interior is seen at native density without misrepresenting it.
cell:
	@$(PY) aleph/viz/cell_app.py $(ARGS)

# The three tau-III seeds, read in the declared order and rendered, as ONE target.
#
# ⚠ Three separate commands is three chances to do two of them. The night this was written, the runs
#   landed with their fields read and NO FIGURE OF ANY KIND -- the charter's "every run writes its own
#   figures" is easy to satisfy and easy to forget, and forgetting it leaves no trace. The read comes
#   FIRST because contract III says nothing later may be quoted if something earlier refuses, and the
#   figures are drawn either way: a refused series is drawn exactly like an accepted one.
TAU3 := $(wildcard aleph/outputs/ac/world_phase4/tau3_seed*.json)
tau3:
	@test -n "$(TAU3)" || { echo "no tau3_seed*.json records found — nothing to read"; exit 2; }
	@echo "── records: $(TAU3)"
	@$(PY) aleph/scripts/world_tau_read.py $(TAU3)
	@echo
	@$(PY) aleph/scripts/world_tau_transient_table.py $(TAU3)
	@echo
	@$(PY) aleph/scripts/world_tau_figures.py $(TAU3)

# The same, cut open on z and written to a file instead of a window — for a machine with no display.
cell-png:
	@$(PY) aleph/viz/cell_app.py --cut z- --view oblique --headless /tmp/cell.png $(ARGS)
