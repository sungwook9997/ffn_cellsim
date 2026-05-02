"""Tk desktop room for the PI / Claude / Codex collaboration channel.

This local app uses the same SQLite DB and markdown ledger as the MCP server
and browser room. It is intentionally message-only: no command execution, no
file reads from artifact paths, and no attempt to control Claude/Codex sessions.
"""

from __future__ import annotations

import json
import tkinter as tk
from tkinter import messagebox, ttk
from typing import Any

try:
    from room import (
        ALLOWED_STATUSES,
        DB_PATH,
        LEDGER_PATH,
        PI_AUTHOR,
        _init_schema,
        _insert_pi_message,
        _load_state,
        _parse_refs,
    )
except ImportError:  # pragma: no cover - supports python -m execution later.
    from .room import (
        ALLOWED_STATUSES,
        DB_PATH,
        LEDGER_PATH,
        PI_AUTHOR,
        _init_schema,
        _insert_pi_message,
        _load_state,
        _parse_refs,
    )


REFRESH_MS = 2000


class DesktopRoom(tk.Tk):
    """Small desktop client over the shared collab SQLite DB."""

    def __init__(self) -> None:
        super().__init__()
        self.title("ACS Workroom")
        self.geometry("1180x760")
        self.minsize(900, 580)
        self._last_rendered_ids: tuple[int, ...] = ()
        self._refresh_job: str | None = None
        self._build_style()
        self._build_layout()
        _init_schema()
        self.refresh()

    def _build_style(self) -> None:
        style = ttk.Style(self)
        if "clam" in style.theme_names():
            style.theme_use("clam")
        style.configure("TFrame", background="#f5f6f8")
        style.configure("Panel.TFrame", background="#ffffff", relief="solid", borderwidth=1)
        style.configure("TLabel", background="#f5f6f8", foreground="#1d232f")
        style.configure("Panel.TLabel", background="#ffffff", foreground="#1d232f")
        style.configure("Muted.TLabel", background="#ffffff", foreground="#697386")
        style.configure("Title.TLabel", background="#f5f6f8", font=("Helvetica", 18, "bold"))
        style.configure("TButton", padding=(10, 6))

    def _build_layout(self) -> None:
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        header = ttk.Frame(self, padding=(14, 10))
        header.grid(row=0, column=0, sticky="ew")
        header.columnconfigure(1, weight=1)
        ttk.Label(header, text="ACS Workroom", style="Title.TLabel").grid(row=0, column=0, sticky="w")
        db_text = f"DB: {DB_PATH}"
        if LEDGER_PATH is not None:
            db_text += f"   Ledger: {LEDGER_PATH}"
        ttk.Label(header, text=db_text).grid(row=1, column=0, columnspan=2, sticky="w")
        ttk.Button(header, text="Refresh", command=self.refresh).grid(row=0, column=2, rowspan=2, padx=(12, 0))

        body = ttk.Frame(self, padding=0)
        body.grid(row=1, column=0, sticky="nsew", padx=12, pady=(0, 12))
        body.columnconfigure(0, weight=1)
        body.columnconfigure(1, minsize=320)
        body.rowconfigure(0, weight=1)

        left = ttk.Frame(body, style="Panel.TFrame", padding=10)
        left.rowconfigure(1, weight=1)
        left.columnconfigure(0, weight=1)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 12))

        composer = ttk.Frame(left, style="Panel.TFrame", padding=10)
        composer.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        for idx in range(6):
            composer.columnconfigure(idx, weight=1)

        ttk.Label(composer, text="Topic", style="Panel.TLabel").grid(row=0, column=0, sticky="w")
        self.topic_var = tk.StringVar(value="v2-layer-1")
        ttk.Entry(composer, textvariable=self.topic_var).grid(row=1, column=0, sticky="ew", padx=(0, 8))

        ttk.Label(composer, text="To", style="Panel.TLabel").grid(row=0, column=1, sticky="w")
        self.to_var = tk.StringVar(value="claude,codex")
        ttk.Entry(composer, textvariable=self.to_var).grid(row=1, column=1, sticky="ew", padx=(0, 8))

        ttk.Label(composer, text="Status", style="Panel.TLabel").grid(row=0, column=2, sticky="w")
        self.status_var = tk.StringVar(value="open-question")
        status_box = ttk.Combobox(
            composer,
            textvariable=self.status_var,
            values=sorted(ALLOWED_STATUSES),
            state="readonly",
        )
        status_box.grid(row=1, column=2, sticky="ew", padx=(0, 8))

        ttk.Label(composer, text="Refs", style="Panel.TLabel").grid(row=0, column=3, columnspan=2, sticky="w")
        self.refs_var = tk.StringVar()
        ttk.Entry(composer, textvariable=self.refs_var).grid(row=1, column=3, columnspan=2, sticky="ew", padx=(0, 8))

        ttk.Button(composer, text="Send as PI", command=self.send_message).grid(row=1, column=5, sticky="ew")

        self.message_text = tk.Text(composer, height=4, wrap="word", undo=True)
        self.message_text.grid(row=2, column=0, columnspan=6, sticky="ew", pady=(8, 0))

        timeline_frame = ttk.Frame(left, style="Panel.TFrame")
        timeline_frame.grid(row=1, column=0, sticky="nsew")
        timeline_frame.rowconfigure(0, weight=1)
        timeline_frame.columnconfigure(0, weight=1)

        self.timeline = tk.Text(
            timeline_frame,
            wrap="word",
            state="disabled",
            padx=10,
            pady=10,
            background="#ffffff",
            foreground="#1d232f",
            relief="flat",
        )
        self.timeline.grid(row=0, column=0, sticky="nsew")
        timeline_scroll = ttk.Scrollbar(timeline_frame, orient="vertical", command=self.timeline.yview)
        timeline_scroll.grid(row=0, column=1, sticky="ns")
        self.timeline.configure(yscrollcommand=timeline_scroll.set)

        self.timeline.tag_configure("meta", foreground="#697386", font=("Helvetica", 11))
        self.timeline.tag_configure("pi", background="#fff6d9", spacing1=6, spacing3=6)
        self.timeline.tag_configure("claude", background="#eef6ff", spacing1=6, spacing3=6)
        self.timeline.tag_configure("codex", background="#edf9f1", spacing1=6, spacing3=6)
        self.timeline.tag_configure("body", font=("Helvetica", 13))
        self.timeline.tag_configure("seen", foreground="#2f7a45", font=("Helvetica", 11))
        self.timeline.tag_configure("pending", foreground="#697386", font=("Helvetica", 11))

        right = ttk.Frame(body, style="Panel.TFrame", padding=10)
        right.rowconfigure(1, weight=1)
        right.rowconfigure(3, weight=1)
        right.rowconfigure(5, weight=1)
        right.columnconfigure(0, weight=1)
        right.grid(row=0, column=1, sticky="nsew")

        ttk.Label(right, text="Active Claims", style="Panel.TLabel").grid(row=0, column=0, sticky="w")
        self.claims = tk.Listbox(right, height=7, activestyle="none")
        self.claims.grid(row=1, column=0, sticky="nsew", pady=(4, 12))

        ttk.Label(right, text="Artifacts", style="Panel.TLabel").grid(row=2, column=0, sticky="w")
        self.artifacts = tk.Listbox(right, height=9, activestyle="none")
        self.artifacts.grid(row=3, column=0, sticky="nsew", pady=(4, 12))

        ttk.Label(right, text="Cursors", style="Panel.TLabel").grid(row=4, column=0, sticky="w")
        self.cursors = tk.Listbox(right, height=5, activestyle="none")
        self.cursors.grid(row=5, column=0, sticky="nsew", pady=(4, 0))

    def send_message(self) -> None:
        body = self.message_text.get("1.0", "end").strip()
        try:
            result = _insert_pi_message(
                topic=self.topic_var.get(),
                body=body,
                addressee=self.to_var.get(),
                status=self.status_var.get(),
                refs=_parse_refs(self.refs_var.get()),
            )
        except ValueError as exc:
            messagebox.showerror("Message not sent", str(exc), parent=self)
            return
        self.message_text.delete("1.0", "end")
        self.refs_var.set("")
        self.refresh(force=True)
        self.status_message(f"Sent message #{result['id']} as {PI_AUTHOR}")

    def refresh(self, force: bool = False) -> None:
        if self._refresh_job is not None:
            try:
                self.after_cancel(self._refresh_job)
            except tk.TclError:
                pass
            self._refresh_job = None
        try:
            state = _load_state()
        except Exception as exc:  # pragma: no cover - UI guard.
            self.status_message(f"Refresh failed: {exc}")
            self._refresh_job = self.after(REFRESH_MS, self.refresh)
            return

        message_ids = tuple(int(row[0]) for row in state["messages"])
        if force or message_ids != self._last_rendered_ids:
            self._last_rendered_ids = message_ids
            self._render_timeline(state)
        self._render_sidebars(state)
        self._refresh_job = self.after(REFRESH_MS, self.refresh)

    def _render_timeline(self, state: dict[str, Any]) -> None:
        cursor_map = {author: int(last_seen) for author, last_seen in state["cursors"]}
        self.timeline.configure(state="normal")
        self.timeline.delete("1.0", "end")
        for row in state["messages"]:
            self._append_message(row, cursor_map)
        self.timeline.configure(state="disabled")
        self.timeline.see("end")

    def _append_message(self, row: tuple[Any, ...], cursor_map: dict[str, int]) -> None:
        msg_id, ts, author, addressee, topic, body, status, refs_json = row
        author_tag = author if author in {"pi", "claude", "codex"} else "body"
        meta = f"#{msg_id}  {author} -> {addressee}  [{status}]  {topic}  {ts}\n"
        self.timeline.insert("end", meta, ("meta", author_tag))
        self.timeline.insert("end", f"{body}\n", ("body", author_tag))
        try:
            refs = json.loads(refs_json)
        except json.JSONDecodeError:
            refs = []
        if refs:
            self.timeline.insert("end", f"Refs: {', '.join(refs)}\n", ("meta", author_tag))
        seen_by = [
            name
            for name in ("claude", "codex")
            if name != author and cursor_map.get(name, 0) >= int(msg_id)
        ]
        if seen_by:
            self.timeline.insert("end", f"Seen by {', '.join(seen_by)}\n\n", ("seen", author_tag))
        else:
            self.timeline.insert("end", "Not seen by Claude/Codex yet\n\n", ("pending", author_tag))

    def _render_sidebars(self, state: dict[str, Any]) -> None:
        self._replace_listbox(
            self.claims,
            [
                f"{topic} | {author} | {summary or ''}"
                for topic, author, _claimed_at, summary in state["claims"]
            ],
            empty="No active claims",
        )
        self._replace_listbox(
            self.artifacts,
            [
                f"#{artifact_id} {kind} | {author} | {path}"
                for artifact_id, _ts, author, kind, path, _note in state["artifacts"]
            ],
            empty="No artifacts",
        )
        self._replace_listbox(
            self.cursors,
            [f"{author}: last seen #{last_seen}" for author, last_seen in state["cursors"]],
            empty="No cursors yet",
        )

    def _replace_listbox(self, box: tk.Listbox, values: list[str], empty: str) -> None:
        box.delete(0, "end")
        for value in values or [empty]:
            box.insert("end", value)

    def status_message(self, text: str) -> None:
        self.title(f"ACS Workroom - {text}")
        self.after(4000, lambda: self.title("ACS Workroom"))


def main() -> None:
    app = DesktopRoom()
    app.mainloop()


if __name__ == "__main__":
    main()
