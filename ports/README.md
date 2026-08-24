# ports — the ledger, and the direction it now runs in

**Copy nothing by default. Port only after evidence.**

Every port needs a ledger entry in `ports/ledger/`, written **before** the code lands, and
`tests/ports/test_port_discipline.py` fails the suite when one is missing a mandatory heading or
still carries a placeholder. The template without the test is a convention, and a convention is
what this tree already had.

## The direction reversed on 2026-08-09

This mechanism was built in `Project_Aleph` to carry code **out of** `ffn_cellsim` under a firewall,
because the clean-room rule was the whole point of that tree. The PI decided on 2026-08-09 that
`ffn_cellsim` is the base and is renamed Project Aleph, so the arrow now points the other way:

```text
Project_Aleph @ an immutable commit
  -> port proposal (this ledger, written first)
  -> an oracle or control this tree owns
  -> rewrite here
  -> tests and, where the claim is physical, native measurement
  -> ACCEPTED or REJECTED
```

The 136 entries `Project_Aleph/ports/ledger/` already holds **travel with the mechanism** — they are
the record of why each contract exists, and dropping them would make every rule here look arbitrary.
They keep their `ALEPH-PORT-36xx` and lower IDs and are not renumbered.

**Ports in the new direction start at `ALEPH-PORT-4001`.** A gap between `-3665` and `-4001` is
deliberate: the ID says which way the code moved without reading the entry.

## Two rules that decide whether an entry is worth writing

- **The unit of approval is a function, kernel, dataclass, or small cohesive module.** A directory is
  never the unit of approval. If you cannot name the symbol, stop.
- **Passing the source's tests is not evidence, and matching the source numerically is not evidence.**
  The entry must stand on a control this tree owns — one that includes a **vacuity** control, i.e. a
  demonstration that the test can fail.

## Status values

`PROPOSED` → `AUDITED` → `ACCEPTED`, or `REJECTED`. A rejected entry stays; the reason a port was
refused is worth more than the silence of never having written it down.
