from aleph.outer.experiment_factory.verification.verify_factory import (
    duplicate_count,
    walk,
)


def test_duplicate_count_is_identity_exact():
    assert duplicate_count(["a", "b", "a"]) == 1
    assert duplicate_count(["A", "a"]) == 0


def test_walk_reaches_nested_raw_keys():
    pairs = list(walk({"safe": [{"evidence_text": "hidden"}]}))
    assert ("evidence_text", "hidden") in pairs


def test_walk_does_not_treat_string_values_as_keys():
    pairs = list(walk({"limitation": "evidence_text is forbidden"}))
    assert ("evidence_text", "evidence_text is forbidden") not in pairs

