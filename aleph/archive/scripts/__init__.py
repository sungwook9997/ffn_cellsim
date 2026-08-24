"""Drivers and one-off analyses that are neither a gate nor a tier-(a) producer.

Selected by decision A6: kept if the script is a gate, produced a `STATE.md` (b) row, is named or
imported by a test, or is one hop from something that is. **194 of 247 scripts landed here.**

They are moved rather than deleted because "nothing imports it" has been wrong six times in one
session on this tree, and a driver is inert until somebody runs it. If one of these turns out to be
the only way to reproduce a row, move it back and say why.
"""
