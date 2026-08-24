"""Unit tests for the modular TAG pipeline and prompt-cache request shape."""
from __future__ import annotations

import logging
import pathlib
import subprocess
import sys
from types import SimpleNamespace

import duckdb
import pytest

TAG_KB = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TAG_KB))

import tag_backend  # noqa: E402
import tag_exec  # noqa: E402
import tag_query  # noqa: E402
import tag_syn  # noqa: E402


class _Usage:
    def model_dump(self):
        return {
            "input_tokens": 17,
            "cache_creation_input_tokens": 2000,
            "cache_read_input_tokens": 0,
            "output_tokens": 11,
        }


def test_anthropic_request_marks_only_static_prefix(monkeypatch):
    captured = {}

    class _Messages:
        def create(self, **kwargs):
            captured.update(kwargs)
            return SimpleNamespace(
                content=[SimpleNamespace(type="text", text="done")],
                usage=_Usage(),
            )

    class _Client:
        def __init__(self, api_key):
            assert api_key == "test-key"
            self.messages = _Messages()

    import anthropic

    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key")
    monkeypatch.setattr(anthropic, "Anthropic", _Client)
    result = tag_backend.complete(
        "Q: dynamic",
        system="SYS",
        cache_prefix="SCHEMA_SHA256: abc\n\nstatic schema",
        cache_key="abc",
        cache_ttl="1h",
    )

    assert result.text == "done"
    assert captured["model"] == "claude-opus-5"
    assert "temperature" not in captured
    assert captured["thinking"] == {"type": "disabled"}
    blocks = captured["messages"][0]["content"]
    assert blocks[0]["text"].endswith("static schema")
    assert blocks[0]["cache_control"] == {"type": "ephemeral", "ttl": "1h"}
    assert blocks[1] == {"type": "text", "text": "Q: dynamic"}


def test_cli_fallback_warns_and_isolates_completion(monkeypatch, caplog):
    captured = {}

    def _run(cmd, **kwargs):
        captured["cmd"] = cmd
        captured.update(kwargs)
        return subprocess.CompletedProcess(
            cmd,
            0,
            stdout='{"result":"ok","usage":{"input_tokens":9}}',
            stderr="",
        )

    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    monkeypatch.setattr(tag_backend.shutil, "which", lambda _: "/fake/claude")
    monkeypatch.setattr(tag_backend.subprocess, "run", _run)
    caplog.set_level(logging.WARNING)

    result = tag_backend.complete(
        "dynamic",
        system="system",
        cache_prefix="static",
        cache_key="f" * 64,
    )

    assert result.backend == "claude-cli"
    assert result.usage == {"input_tokens": 9}
    assert "prompt-cache BYPASS" in caplog.text
    assert captured["input"] == "static\n\ndynamic"
    assert captured["cmd"][0] == "/fake/claude"
    assert captured["cmd"][captured["cmd"].index("--tools") + 1] == ""
    assert "--system-prompt" in captured["cmd"]
    assert pathlib.Path(captured["cwd"]).name.startswith("tag-kb-claude-")


def test_schema_hash_changes_cache_prefix_when_schema_changes():
    first_prefix, first_hash = tag_syn.static_cache_prefix("table_a(id)")
    second_prefix, second_hash = tag_syn.static_cache_prefix("table_a(id, status)")

    assert first_hash != second_hash
    assert f"SCHEMA_SHA256: {first_hash}" in first_prefix
    assert f"SCHEMA_SHA256: {second_hash}" in second_prefix


def test_schema_text_is_stable_across_connections(tmp_path):
    db_path = tmp_path / "schema.duckdb"
    con = duckdb.connect(str(db_path))
    con.execute(
        "CREATE TABLE edges "
        "(src_id TEXT, src_type TEXT, rel TEXT, dst_id TEXT, dst_type TEXT)"
    )
    con.execute("CREATE TABLE records (id TEXT, status TEXT)")
    con.executemany(
        "INSERT INTO records VALUES (?, ?)",
        [("a", "zeta"), ("b", "alpha"), ("c", "middle")],
    )
    first = tag_syn.schema_text(con)
    con.close()

    con = duckdb.connect(str(db_path))
    second = tag_syn.schema_text(con)
    con.close()

    assert first == second
    assert "status ∈ {alpha, middle, zeta}" in first


def test_syn_keeps_schema_static_and_question_dynamic():
    con = duckdb.connect(":memory:")
    con.execute(
        "CREATE TABLE edges "
        "(src_id TEXT, src_type TEXT, rel TEXT, dst_id TEXT, dst_type TEXT)"
    )
    con.execute("CREATE TABLE source_evidence (id TEXT, status TEXT)")
    expected_schema = tag_syn.schema_text(con)
    captured = {}

    def _llm(prompt, **kwargs):
        captured["prompt"] = prompt
        captured.update(kwargs)
        return "```sql\nSELECT count(*) FROM source_evidence;\n```"

    sql = tag_syn.syn(
        con,
        "How many sources?",
        None,
        llm_fn=_llm,
    )
    con.close()

    assert sql == "SELECT count(*) FROM source_evidence"
    assert captured["prompt"] == "Q: How many sources?\n```sql\n"
    assert "source_evidence(id, status)" in captured["cache_prefix"]
    assert captured["cache_key"] == tag_syn.schema_hash(expected_schema)


@pytest.mark.parametrize(
    "sql",
    [
        "DELETE FROM source_evidence",
        "SELECT 1; SELECT 2",
        "PRAGMA show_tables",
        "WITH x AS (SELECT 1) INSERT INTO t SELECT * FROM x",
    ],
)
def test_exec_rejects_non_readonly_sql(sql):
    with pytest.raises(ValueError):
        tag_exec.validate_select(sql)


def test_tag_query_stays_a_thin_compatibility_surface():
    assert tag_query.exec_sql is tag_exec.exec_sql
    assert tag_query.DEFAULT_MODEL == "claude-opus-5"
    assert callable(tag_query.syn)
    assert callable(tag_query.gen)
