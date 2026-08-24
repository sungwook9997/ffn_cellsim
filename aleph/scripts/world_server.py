"""A warm GPU executor that lives inside one Slurm hold, so parallel sessions stop paying for cold starts.

WHAT THIS IS FOR.  The PI's grant is ONE hold on ONE card, and four sessions need it. Without this each
of them would ``gpu-submit`` — fighting over the licence, or drifting onto a card nobody reserved, which
``outputs/ac/gpu_grants.md`` records as a real trap (``srun --overlap`` does not inherit the wrap's
``CUDA_VISIBLE_DEVICES``). So the hold is held once, here, and sessions post work to it.

AND THE COLD START IS NOT SMALL.  Measured on job 69: ``build_wall_s`` 23.7 s for a 551,434-node cell,
before any step. Warp's kernel cache compiles on first launch on top of that. Paid once per process
here instead of once per submission.

WHAT IT DELIBERATELY IS NOT, YET.  Not a world-resident API. Keeping a built arena alive across requests
is the point of the exercise, but the populations do not exist yet — PHASE 1 is being written by two
sessions right now — and an API designed before its object is the thing ``arena.py`` warns about
("guessing would bake an answer into the layout"). This runs whole drivers in a warm process and
records what they cost. The resident-world endpoints get added when there is a world to keep resident.

SECURITY POSTURE, stated because it is a server.  It binds **127.0.0.1 only** and is reached by an SSH
tunnel; it is not exposed on Tailscale or anywhere else. It executes a driver path under the repo root
by an allow-list of one directory (``aleph/scripts``), with arguments passed through — this is a
single-user box and the caller is another session of the same user, but the allow-list means a typo
cannot run an arbitrary binary.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — wall times [s], bytes [B]. Nothing physical is computed here.
  * boundary — a request naming a driver outside ``aleph/scripts`` is refused, not sanitised; a request
    while another job is running is refused with the running job's name rather than queued silently.
  * conservation/invariant — one job at a time, by a lock. Two concurrent CUDA drivers in one process
    would share a context and their peak-memory readings would be each other's.
  * CFL/precision — not applicable.
  * sign sense — not applicable.
  * measurement protocol — ``/probe`` reads the CUDA driver's own free/total, not a computed figure, and
    reports the delta across a job. It is a host-side read between jobs, never inside one.

Runtime: NVIDIA Warp on CUDA, inside a Slurm allocation. Refuses to start without a CUDA device.
"""

from __future__ import annotations

import json
import subprocess
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DRIVER_DIR = REPO / "aleph" / "scripts"
_LOCK = threading.Lock()
_STATE: dict = {"running": None, "history": []}


def _cuda_or_die() -> str:
    """Refuse to start off a CUDA device. The charter forbids a CPU path; a server is not an exception."""
    import warp as wp

    wp.init()
    devices = [d for d in wp.get_devices() if d.is_cuda]
    if not devices:
        raise RuntimeError(
            "world_server requires a CUDA device. Outside a Slurm allocation cuInit() returns "
            "CUDA_ERROR_NO_DEVICE — start it with gpu-submit, never bare."
        )
    if len(devices) != 1:
        raise RuntimeError(
            f"expected exactly one visible CUDA device inside the hold, saw {len(devices)}: "
            f"{[str(d) for d in devices]}. A second visible card means the pinning did not take."
        )
    return str(devices[0])


def _mem() -> dict:
    """Free/total from the CUDA driver itself, not arithmetic over what we think we allocated."""
    try:
        import warp as wp

        free, total = wp.get_device("cuda").memory_info()
        return {"free_bytes": int(free), "total_bytes": int(total)}
    except Exception as exc:  # noqa: BLE001 — a probe that cannot read says so rather than guessing
        return {"error": f"{type(exc).__name__}: {exc}"}


def _run(driver: str, args: list[str]) -> dict:
    """Run one driver in this warm process's environment. One at a time, by lock."""
    path = (DRIVER_DIR / driver).resolve()
    if path.parent != DRIVER_DIR or not path.is_file():
        return {"error": f"refused: {driver!r} is not a file directly under aleph/scripts/"}
    if not _LOCK.acquire(blocking=False):
        return {"error": f"busy: {_STATE['running']} is running. One job at a time — two CUDA drivers "
                         "in one context would read each other's peak memory."}
    try:
        _STATE["running"] = driver
        before = _mem()
        t0 = time.time()
        done = subprocess.run([sys.executable, str(path), *args], cwd=str(REPO),
                              capture_output=True, text=True)
        wall = time.time() - t0
        rec = {"driver": driver, "args": args, "returncode": done.returncode, "wall_s": round(wall, 2),
               "mem_before": before, "mem_after": _mem(),
               "stdout_tail": done.stdout[-4000:], "stderr_tail": done.stderr[-4000:]}
        _STATE["history"].append({k: rec[k] for k in ("driver", "returncode", "wall_s")})
        return rec
    finally:
        _STATE["running"] = None
        _LOCK.release()


class Handler(BaseHTTPRequestHandler):
    def _send(self, code: int, payload: dict) -> None:
        body = json.dumps(payload, indent=2, default=str).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt, *a):  # noqa: A003 — quieter than the default one-line-per-request
        print(f"[world_server] {self.address_string()} {fmt % a}", flush=True)

    def do_GET(self) -> None:  # noqa: N802 — BaseHTTPRequestHandler's interface
        if self.path == "/health":
            self._send(200, {"ok": True, "device": _STATE["device"], "running": _STATE["running"],
                             "jobs_done": len(_STATE["history"])})
        elif self.path == "/probe":
            self._send(200, {"memory": _mem(), "running": _STATE["running"],
                             "history": _STATE["history"][-20:]})
        elif self.path == "/drivers":
            self._send(200, {"drivers": sorted(p.name for p in DRIVER_DIR.glob("*.py"))})
        else:
            self._send(404, {"error": "GET /health | /probe | /drivers ; POST /run"})

    def do_POST(self) -> None:  # noqa: N802
        if self.path != "/run":
            self._send(404, {"error": "POST /run with {\"driver\": \"x.py\", \"args\": [...]}"})
            return
        n = int(self.headers.get("Content-Length", 0))
        try:
            req = json.loads(self.rfile.read(n) or b"{}")
        except json.JSONDecodeError as exc:
            self._send(400, {"error": f"bad JSON: {exc}"})
            return
        driver = req.get("driver")
        if not isinstance(driver, str):
            self._send(400, {"error": "need {\"driver\": \"<name>.py\"}"})
            return
        rec = _run(driver, [str(a) for a in req.get("args", [])])
        self._send(200 if "error" not in rec else 409, rec)


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--port", type=int, default=8080)
    args = ap.parse_args()
    _STATE["device"] = _cuda_or_die()
    print(f"[world_server] device={_STATE['device']} repo={REPO}", flush=True)
    print(f"[world_server] memory at start: {_mem()}", flush=True)
    srv = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"[world_server] listening on 127.0.0.1:{args.port} — tunnel in with "
          f"ssh -L {args.port}:localhost:{args.port}", flush=True)
    srv.serve_forever()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
