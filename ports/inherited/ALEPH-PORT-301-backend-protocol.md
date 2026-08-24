# Port ledger — ALEPH-PORT-301 · runtime backend protocol + GPU authorization guard

| Field | Value |
|---|---|
| Aleph target | `aleph/runtime/backend.py` |
| Lane | L3 runtime |
| Written | 2026-07-30 (before the code landed, per PLAN §0.2.5) |
| Status | LANDED |
| Evidence class | STRUCTURAL (CPU reference implementation + guard tests). Not physics. |

## What was read in `ffn_cellsim` (read-only)

- `ffn_sim/ac/engine/transaction.py`, `world.py`, `ledger.py` — all three `import warp as wp` at module
  scope. The reference project has no backend seam at all: Warp is a hard, unconditional dependency of
  its transaction layer, so nothing in that layer is importable, let alone testable, without CUDA.

## What was re-derived, not copied

Nothing was ported. The observation is *negative* — the reference has no backend abstraction — so
Aleph's `Backend` Protocol is original work, written to satisfy PLAN §2 `ALEPH-DQ-107`, which ratifies
Warp "on availability + auditability, NOT inherited" and explicitly notes the reversal cost is low
"because the runtime layer is written against a backend Protocol". This file is that Protocol.

Design decisions taken here and their reasons:

1. `NumpyBackend` is the *reference* implementation, not a fallback. Every physical law in Aleph must
   have a CPU expression that a human can read. Warp is an accelerator for a law defined elsewhere.
2. `warp` is imported inside `WarpBackend.__post_init__`, never at module scope. A test asserts
   `"warp" not in sys.modules` after importing `aleph.runtime`.
3. `op_count` — the backend counts compute launches. This is not profiling. It is the *witness* the
   dispatch pipeline (ALEPH-PORT-303) uses to prove a participant actually evaluated something rather than
   returning silently. Host reads (`to_host`) and `synchronize()` deliberately do not count.

## GPU authorization guard — why it exists

PLAN §0.1: the shared workstation is **not authorized**, the PI grants access explicitly and for a
stated duration each time, and "any script that could launch GPU work must fail closed with an explicit
'no GPU authorization on record' error unless a signed authorization file is present."

`require_gpu_authorization()` implements exactly that and is called by `WarpBackend` **before** warp is
imported, so an unauthorized process never even loads a CUDA runtime. It fails closed on: missing file,
unparseable file, absent or unparseable expiry, expiry in the past, hostname mismatch, an absent device
list, and a device id not in the granted set — the same seven rules, with the same wording, as
`scripts/gpu_preflight.py::load_authorization`. A timezone-naive `expires_at` is read as UTC rather
than refused, again matching the preflight: two conventions for one permission slip would be worse
than none.

## Controls

- Positive control: `test_authorization_accepts_a_valid_grant` — a well-formed, unexpired, host-matching file
  returns a `GPUAuthorization` with the granted device set and a positive remaining duration.
- Negative (must fail): `test_authorization_missing_file_refuses`, `test_authorization_expired_refuses`,
  `test_authorization_wrong_host_refuses`, `test_authorization_device_not_granted_refuses`,
  `test_authorization_without_expiry_refuses`, `test_warp_backend_checks_authorization_before_importing_warp`.

## Provenance and envelope

- **Source repository / commit read:** `ffn_cellsim` at commit `be0e5876`, read-only.
- **Why not clean-room:** it *was* clean-room. Nothing was ported; the finding is negative (the
  reference has no backend seam at all), so no source-derived code exists to justify.
- **Discarded prose:** none to discard — no comments or docstrings were carried, because no code was.
  All prose here is new and written for Aleph's situation.
- **Units, domains, invariants:** the backend is dimensionless; units are the caller's. Domain
  invariants enforced here: `op_count` is monotonic and never advanced by a host read; `zeros`
  allocates float64 by default; `scatter_add` accumulates rather than overwrites on repeated indices.
- **Numerical and precision envelope:** float64 throughout on the CPU reference. No tolerance is
  chosen at this layer. The accelerator path (float32 compute / float64 accumulation, `ALEPH-DQ-107`)
  is unimplemented and must be qualified against `NumpyBackend` before it is used for any number.
- **Residency:** host only. The device path is guarded and unreached in this session.

## Independence

No `ffn_cellsim` import, no `ffn_cellsim` string, no shared identifier beyond ordinary English words.
