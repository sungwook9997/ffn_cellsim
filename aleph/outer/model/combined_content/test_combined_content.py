import hashlib
import importlib.util
import json
import unittest
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("combined_infer_test", HERE / "infer.py")
infer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(infer)


class CombinedContentControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metrics = json.loads((HERE / "metrics.json").read_text())
        cls.splits = [json.loads(line) for line in (HERE / "splits.jsonl").read_text().splitlines()]
        cls.index = infer.CombinedContentIndex()

    def test_combined_family_counts_and_support(self):
        self.assertEqual(self.metrics["input"]["usable_unique_families"], 5669)
        self.assertEqual(self.metrics["input"]["derived"]["first"]["usable"], 2695)
        self.assertEqual(self.metrics["input"]["derived"]["second"]["usable"], 2974)
        families = [row["source_family_id"] for row in self.splits]
        self.assertEqual(len(families), len(set(families)))
        labels = {row["label"] for row in self.splits}
        self.assertEqual(len(labels), 12)
        for split in ("train", "validation", "test"):
            self.assertEqual({row["label"] for row in self.splits if row["split"] == split}, labels)

    def test_dataset_components_never_cross_splits(self):
        seen = defaultdict(set)
        for row in self.splits:
            seen[row["group_key"]].add(row["split"])
        self.assertTrue(all(len(values) == 1 for values in seen.values()))
        self.assertEqual(self.metrics["leakage"]["group_split_violations"], 0)
        self.assertEqual(self.metrics["leakage"]["usable_grouped_families"], 140)

    def test_taxonomy_maps_every_observed_weak_label_once(self):
        taxonomy = json.loads((HERE / "taxonomy.json").read_text())
        for pass_name, mapping_name in (("first", "first_pass_mapping"), ("second", "second_pass_mapping")):
            observed = {row["weak_label"] for row in self.splits if row["pass"] == pass_name}
            self.assertEqual(observed, set(taxonomy[mapping_name]))

    def test_content_improves_same_split_baseline(self):
        content = self.metrics["combined_content_plus_metadata"]
        title = self.metrics["title_only"]
        self.assertGreater(content["macro_f1"], title["macro_f1"])
        self.assertGreater(content["retrieval"]["precision_at_k"], title["retrieval"]["precision_at_k"])

    def test_v2_metrics_are_declared_incompatible(self):
        self.assertFalse(self.metrics["v2_coverage_comparison_only"]["metric_comparison_permitted"])

    def test_receipts_match(self):
        receipts = json.loads((HERE / "SHA256SUMS.json").read_text())
        for name, expected in receipts.items():
            self.assertEqual(hashlib.sha256((HERE / name).read_bytes()).hexdigest(), expected, name)

    def test_inference_deterministic_unique_and_self_excluding(self):
        family = "doi:10.1038/ncb3525"
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
        self.assertEqual(actual["coarse_route_probabilities"][0]["route"], expected["top_coarse_route"]["route"])
        self.assertAlmostEqual(actual["coarse_route_probabilities"][0]["probability"], expected["top_coarse_route"]["probability"])
        self.assertEqual([row["source_family_id"] for row in actual["retrieval"]], [row["source_family_id"] for row in expected["retrieval"]])

    def test_no_sealed_package_imports(self):
        for path in (HERE / "train.py", HERE / "infer.py"):
            source = path.read_text()
            self.assertNotIn("aleph.learn", source)
            self.assertNotIn("aleph.represent", source)


if __name__ == "__main__":
    unittest.main()
