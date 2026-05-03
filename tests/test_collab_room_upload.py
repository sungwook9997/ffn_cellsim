from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


def _load_room(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("COLLAB_MCP_DB", str(tmp_path / "inbox.db"))
    monkeypatch.setenv("COLLAB_ROOM_UPLOAD_DIR", str(tmp_path / "uploads"))
    monkeypatch.setenv("COLLAB_ROOM_TOKEN", "unit-test-room")
    monkeypatch.delenv("COLLAB_MCP_LEDGER", raising=False)
    module_path = Path(__file__).parents[1] / "tools" / "collab_mcp" / "room.py"
    spec = importlib.util.spec_from_file_location("collab_room_under_test", module_path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules["collab_room_under_test"] = module
    spec.loader.exec_module(module)
    module._init_schema()
    return module


def test_upload_is_saved_and_announced_in_message(monkeypatch, tmp_path: Path):
    room = _load_room(monkeypatch, tmp_path)

    upload = room._store_upload("../PI notes.csv", b"time,value\n0,1\n", "text/csv")
    result = room._insert_pi_message_with_uploads(
        topic="v2-layer-1",
        body="검토용 파일입니다.",
        addressee="claude,codex",
        status="FYI",
        refs=[],
        uploads=[upload],
    )

    assert result["id"] == 1
    saved_path = Path(upload["path"])
    assert saved_path.is_file()
    assert saved_path.parent == tmp_path / "uploads"
    assert ".." not in saved_path.name

    message = room._message_to_dict(room._fetch_recent_messages(limit=1)[0])
    assert str(saved_path) in message["body"]
    assert str(saved_path) in message["refs"]
    assert "artifact #1" in message["body"]

    artifacts = room._sidebar_to_dict(room._fetch_sidebar())["artifacts"]
    assert artifacts[0]["kind"] == "pi_upload"
    assert artifacts[0]["path"] == str(saved_path)
    assert "text/csv" in artifacts[0]["note"]


def test_upload_only_message_gets_default_body(monkeypatch, tmp_path: Path):
    room = _load_room(monkeypatch, tmp_path)

    upload = room._store_upload("image.png", b"\x89PNG\r\n", "image/png")
    room._insert_pi_message_with_uploads(
        topic="v2-layer-1",
        body="",
        addressee="codex",
        status="FYI",
        refs=[],
        uploads=[upload],
    )

    message = room._message_to_dict(room._fetch_recent_messages(limit=1)[0])
    assert message["body"].startswith("첨부 파일을 공유합니다.")
    assert upload["path"] in message["refs"]


def test_pi_messages_are_tagged_and_filtered_by_room(monkeypatch, tmp_path: Path):
    room = _load_room(monkeypatch, tmp_path)

    design = room._insert_pi_message(
        topic="routing",
        body="design note",
        addressee="claude,codex",
        status="FYI",
        refs=[],
        room="design-discussion",
    )
    impl = room._insert_pi_message(
        topic="routing",
        body="implementation note",
        addressee="claude,codex",
        status="FYI",
        refs=[],
        room="implementation-work",
    )

    assert design["room"] == "design-discussion"
    assert impl["room"] == "implementation-work"
    design_messages = [
        room._message_to_dict(row)["body"]
        for row in room._fetch_recent_messages(room="design-discussion")
    ]
    impl_messages = [
        room._message_to_dict(row)["body"]
        for row in room._fetch_recent_messages(room="implementation-work")
    ]
    assert design_messages == ["design note"]
    assert impl_messages == ["implementation note"]
    assert "implementation-work" in room._fetch_rooms()


def test_workroom_prefixed_agent_status_is_accepted(monkeypatch, tmp_path: Path):
    room = _load_room(monkeypatch, tmp_path)

    result = room._upsert_agent_status(
        "design-discussion-claude-chat",
        percent=12,
        activity="waiting",
        topic="room-split",
        heartbeat=True,
    )

    assert result["agent"] == "design-discussion-claude-chat"
    agents = room._sidebar_to_dict(room._fetch_sidebar())["agents"]
    assert agents[0]["agent"] == "design-discussion-claude-chat"
    assert room._agent_base("design-discussion-claude-chat") == "claude"


def test_invalid_workroom_prefixed_agent_status_is_rejected(monkeypatch, tmp_path: Path):
    room = _load_room(monkeypatch, tmp_path)

    with pytest.raises(ValueError, match="<workroom>-<pane>"):
        room._upsert_agent_status("BadRoom-claude-chat", 0, "bad")
