# CITED — four PI approvals given in session `97fb3916`, 2026-08-09 22:04

> **This is a CITATION, not a ratification.** No `decided_by` field appears here and an agent may not
> add one. What is recorded is that the PI said a specific thing at a specific coordinate; whether
> each affected `PROPOSAL-*` may be marked ratified is a separate act, and the last section says why
> two of them should not be yet.
>
> Relayed by session `97fb3916`. **Verified by this session by opening that transcript**, not
> accepted on report — the uuids below are not exposed to the session that relayed them, so this
> record is checkable by a third party in a way the relay was not.

| Field | Value |
|---|---|
| Written | 2026-08-09 22:xx KST, session `5ccfd28d`, lane `contracts` |
| Source transcript | `~/.claude/projects/-Users-sw1-ffn-cellsim/97fb3916-7a35-4be7-b029-d6ad2485617b.jsonl` (6,731,950 B) |
| Verified | 2026-08-09 22:xx by reading the file; every quotation below is byte-checked |

---

## 1. What the PI actually said — one sentence, not four approvals

```
line 1397 · 2026-08-09 22:04:56 KST · uuid 1f846ebe-07eb-43f6-8c08-c34972ccb9ff
"(a) 나머지 승인 그리고 계속 진행"
```

That is the whole of it. **The PI did not enumerate four items.** The four came from the question it
answers — `97fb3916` was asked *"pi 대기 내용 설명"* at 22:02:14
(uuid `549040b9-c061-4309-a684-87b05e4fb37a`) and listed four pending items in its reply at line
1391. The approval's **scope is defined by that list**, which is an assistant message written by the
relaying session, not by the PI.

This is not a complaint about the relay — the relay was accurate, and the enumeration is a fair
reading of what "나머지" covers. It is recorded because **the granularity in the report is the
session's, and a reader six weeks from now would otherwise believe the PI itemised four decisions.**

## 2. The list that "나머지" refers to

| | item | what the approval means for it |
|---|---|---|
| ① | `lamellipodium_nascent_fa` / `filopodium_nascent_fa` — what a nascent adhesion hangs from | **explicit: option (a)** — declare an ECM-side partner edge, so a nascent adhesion is a SERIES like the mature one. Not (b) clamp to `world_boundary`, not (c) a mechanics-free state edge |
| ② | **Card-5**, the interior-column strangler transition — moving membrane / nucleus / cytosol off the shared `cell.pos_d` / `f_d` onto private arrays, which requires editing the feature-frozen `components/incumbent/driver.py` | covered by "나머지 승인" |
| ③ | the five unfrozen `INHERITED-` proposals | covered by "나머지 승인" |
| ④ | `docs/decisions/**` write access for `97fb3916` | covered by "나머지 승인" |

**② is the one that matters most and is easiest to under-read.** `STATE.md` (c) 4 says *"Arrays are
still shared"* and has said so through three engine generations; Card-5 is the change that stops it
being true. Three owners — membrane, cortex, nucleus — are currently constructed against the same
`cell.pos_d` and `cell.f_d`. That is co-location standing in for connection, in the tree whose
charter opens by forbidding exactly that.

## 3. Where each lands, and who owns the surface

| | lands in | owner |
|---|---|---|
| ① | `engine/`, `components/` — and moves the connector census | `97fb3916` (lane `engine`) |
| ② | `components/incumbent/driver.py` | `97fb3916` |
| ③ | `docs/decisions/` | `5ccfd28d` (lane `contracts`) — this session |
| ④ | `ownership.yaml` | this session; **declined, see §5** |

## 4. Two of the five should not be marked ratified yet, and one of the two is the relay's own finding

`97fb3916` flagged this itself rather than acting on it, which is the behaviour the citation
discipline is for:

- **`darcy-or-stokes-for-the-cytosol`** — "approve" reads cleanly as *ratify what is implemented*.
  Biot/Darcy is what runs. Safe.
- **`a-composite-group-is-one-connector-object`** — already true in the code. Safe.
- **`the-far-field-anchor-is-a-clamp-and-the-registry-calls-it-a-connector`** — **not safe.** What is
  implemented *is the thing the proposal disputes*. "Approve" cannot mean "ratify what is
  implemented" for a proposal whose content is that what is implemented is mis-declared. Reading a
  blanket approval onto it would ratify the defect.

Adding to that from this session's side:

- **`focal-adhesion-endpoint-shape`** and **`substrate-adhesion-belongs-to-an-owner`** are the two the
  relay says item ① *is*. If ①'s "(a)" decides them, they are ratified by ① and not by "나머지" — and
  they should be marked with ①'s coordinate, not this file's.

**None of the five is marked ratified by this record.** Each needs its own line, written when the
code that ratifies it lands, citing the coordinate above.

## 5. §4 — `docs/decisions/**` access: declined, and the reason is not authority

`97fb3916` offered either: re-cut the glob to it, or keep the surface and write the ratifications
here as it sends material. **Keeping it**, for the reason `ownership.yaml`'s own comment gives one
line above the lane: *a contract two sessions edit is not a contract.*

That is an engineering answer, not a jurisdictional one. The relay has the PI's approval to write
here and this session is not overriding it — it is taking the offered alternative, because a
decision record with two authors and one hook has no owner when the two disagree. `97fb3916` sends
the material; this session writes it and cites the source. If that becomes a bottleneck, re-cut it.

**And a boundary this session will not cross regardless:** a peer reporting a PI approval is not the
PI approving. This record exists because the transcript was opened and the quotation found. Had the
file been missing or the text absent, the correct output would have been *unverifiable here* — which
is a question for the PI and not a verdict, and is precisely the mistake that cost Aleph ninety
minutes on 2026-07-30 when a session called an authentic instruction fabricated.
