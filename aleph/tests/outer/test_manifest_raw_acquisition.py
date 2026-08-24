"""The manifest's one job that can be wrong quietly: deciding whether a file was verified.

Size and path come from the filesystem and are hard to get wrong. `verified` is a judgement with
three states, and the expensive error is collapsing them into two — a file nothing ever pinned is
not a file that passed. A first version of this module checked filenames instead of records and
called all 64 AllenCell payloads corrupt, so the honest-name case is pinned here too.
"""

from __future__ import annotations

import hashlib
import json

import pytest

from aleph.outer.acquisition import manifest_raw_acquisition as m

PAYLOAD = b"raw acquisition bytes"
DIGEST = hashlib.sha256(PAYLOAD).hexdigest()


def _tree(tmp_path, monkeypatch, *, pin: dict[str, str], names: list[str]):
    """A one-block root holding ``names``, with ``pin`` written as the acquisition record."""
    monkeypatch.setattr(m, "BLOCKS", [("t", "b", "a test")])
    (tmp_path / "b").mkdir()
    for name in names:
        (tmp_path / "b" / name).write_bytes(PAYLOAD)
    (tmp_path / "results").mkdir()
    (tmp_path / m.PINNED).write_text(
        json.dumps({"files": [{"path": f"crop_raw/{k}", "sha256": v} for k, v in pin.items()]})
    )
    return tmp_path


def test_verified_contradicted_and_never_pinned_are_three_states(tmp_path, monkeypatch):
    root = _tree(
        tmp_path,
        monkeypatch,
        pin={"good.tif": DIGEST, "bad.tif": "0" * 64},
        names=["good.tif", "bad.tif", "GSE164378_RAW.tar"],
    )
    got = {r["path"].split("/")[-1]: r for r in m.scan(root)}

    assert got["good.tif"]["verified"] is True
    assert got["bad.tif"]["verified"] is False
    # Never pinned — None, not False, so 1.9 GB of unpinned GEO files cannot read as corrupt,
    # and not True, so it cannot read as checked.
    assert got["GSE164378_RAW.tar"]["verified"] is None
    assert got["GSE164378_RAW.tar"]["verified_against"] is None

    assert m.report(list(got.values())) == 1
    for r in got.values():
        assert r["sha256"] == DIGEST
        assert r["size_bytes"] == len(PAYLOAD)


def test_a_filename_that_looks_like_a_digest_is_not_treated_as_one(tmp_path, monkeypatch):
    """AllenCell's S3 keys lead with 64 hex characters. That is a name, not a claim about bytes."""
    root = _tree(tmp_path, monkeypatch, pin={}, names=[f"{'a' * 64}_raw.ome.tif"])
    assert m.scan(root)[0]["verified"] is None


def test_duplicate_basenames_raise_rather_than_mis_verify(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "BLOCKS", [("t", "b", "a test")])
    for sub in ("x", "y"):
        (tmp_path / "b" / sub).mkdir(parents=True)
        (tmp_path / "b" / sub / "same.tif").write_bytes(PAYLOAD)
    with pytest.raises(ValueError, match="basenames are not unique"):
        m.scan(tmp_path)


def test_absent_block_is_skipped_not_counted_as_empty(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "BLOCKS", [("gone", "nowhere", "a test")])
    assert m.scan(tmp_path) == []
