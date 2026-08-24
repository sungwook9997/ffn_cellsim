"""Package marker so these tests get a qualified module name.

Without it pytest imports `test_contracts.py` here as the bare module ``test_contracts``, which collides
with any same-named file elsewhere under `tests/ac/` — and a collision does not fail one file, it
INTERRUPTS COLLECTION, so the entire suite reports zero tests run.  That happened on 2026-07-29 when a
second `test_contracts.py` appeared under `tests/ac/virtual_cell/`; the suite had been exiting 0 an hour
earlier, so the breakage looked like whatever change happened to be in flight.

`tests/ac/`, `tests/ac/cell/` and `tests/ac/motor/` already carry this file; `engine/` and
`virtual_cell/` did not. Adding it here is enough to disambiguate the pair, and it restores the
convention the sibling directories already follow.
"""
