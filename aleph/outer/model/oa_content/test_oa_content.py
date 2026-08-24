import hashlib
import importlib.util
import json
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
TRAIN_SPEC = importlib.util.spec_from_file_location("oa_train_test", HERE / "train.py")
train = importlib.util.module_from_spec(TRAIN_SPEC)
TRAIN_SPEC.loader.exec_module(train)
INFER_SPEC = importlib.util.spec_from_file_location("oa_infer_test", HERE / "infer.py")
infer = importlib.util.module_from_spec(INFER_SPEC)
INFER_SPEC.loader.exec_module(infer)


class OAContentControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metrics = json.loads((HERE / "metrics.json").read_text())
        cls.splits = [json.loads(line) for line in (HERE / "splits.jsonl").read_text().splitlines()]
        cls.index = infer.OAContentIndex()

    def test_audit_refuses_2701_as_unique_success_count(self):
        audit = self.metrics["audit"]
        self.assertEqual(audit["chunk_files"], 2701)
        self.assertEqual(audit["nonempty_unique_source_families"], 2700)
        self.assertEqual(audit["usable_source_families"], 2695)
        self.assertEqual(audit["empty_files"], ["b07a03af93de9dbc52ab1808.chunks.jsonl"])

    def test_split_is_source_family_unique_and_all_domains_supported(self):
        families = [row["source_family_id"] for row in self.splits]
        self.assertEqual(len(families), len(set(families)))
        labels = {row["label"] for row in self.splits}
        for split in ("train", "validation", "test"):
            self.assertEqual({row["label"] for row in self.splits if row["split"] == split}, labels)

    def test_content_model_is_fair_same_split_improvement(self):
        self.assertGreater(self.metrics["oa_content_plus_metadata"]["macro_f1"], self.metrics["title_only"]["macro_f1"])
        self.assertGreater(self.metrics["oa_content_plus_metadata"]["retrieval"]["precision_at_k"], self.metrics["title_only"]["retrieval"]["precision_at_k"])

    def test_receipts_match(self):
        receipts = json.loads((HERE / "SHA256SUMS.json").read_text())
        for name, expected in receipts.items():
            self.assertEqual(hashlib.sha256((HERE / name).read_bytes()).hexdigest(), expected, name)

    def test_inference_is_deterministic_unique_and_self_excluding(self):
        family = "doi:10.1038/s41556-025-01807-6"
        first = self.index.query_source(family, 5)
        second = self.index.query_source(family, 5)
        self.assertEqual(first, second)
        retrieved = [row["source_family_id"] for row in first["retrieval"]]
        self.assertEqual(len(retrieved), len(set(retrieved)))
        self.assertNotIn(family, retrieved)

    def test_unknown_source_refuses(self):
        with self.assertRaises(infer.SourceRefusal) as caught:
            self.index.query_source("doi:10.0000/not-present")
        self.assertEqual(caught.exception.code, "source_family_unavailable")

    def test_example_inference_matches_frozen_index(self):
        expected = json.loads((HERE / "example_inference.json").read_text())
        actual = self.index.query_source(expected["source_family_id"], 5)
        self.assertEqual(actual["bounded_section_token_counts"], expected["bounded_section_token_counts"])
        self.assertEqual(actual["weak_domain_probabilities"][0]["domain"], expected["top_weak_domain"]["domain"])
        self.assertAlmostEqual(actual["weak_domain_probabilities"][0]["probability"], expected["top_weak_domain"]["probability"])
        self.assertEqual([row["source_family_id"] for row in actual["retrieval"]], [row["source_family_id"] for row in expected["retrieval"]])

    def test_no_sealed_package_imports(self):
        for path in (HERE / "train.py", HERE / "infer.py"):
            source = path.read_text()
            self.assertNotIn("aleph.learn", source)
            self.assertNotIn("aleph.represent", source)


if __name__ == "__main__":
    unittest.main()
