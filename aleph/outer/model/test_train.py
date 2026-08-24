import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("train.py")
SPEC = importlib.util.spec_from_file_location("external_model_train", MODULE_PATH)
train = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(train)

INFER_PATH = Path(__file__).with_name("infer.py")
INFER_SPEC = importlib.util.spec_from_file_location("external_model_infer", INFER_PATH)
infer = importlib.util.module_from_spec(INFER_SPEC)
INFER_SPEC.loader.exec_module(infer)


class ExternalModelControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index = infer.ExternalMetadataIndex()

    def test_split_is_deterministic_and_source_family_exclusive(self):
        values = [train.split_name("doi:example", 7, 0.7, 0.15) for _ in range(3)]
        self.assertEqual(values, [values[0]] * 3)

    def test_historical_year_is_not_used_as_label(self):
        record = {"title": "Atomic force microscopy of living cells"}
        self.assertEqual(train.historical_domain(record), "mechanical_measurement_modalities")
        self.assertNotIn("2010", train.historical_domain(record))

    def test_domain_normalisation_joins_temporal_duplicate_lanes(self):
        self.assertEqual(train.canonical_domain("2010s_epithelial_states"), "epithelial_states")
        self.assertEqual(train.canonical_domain("stem_cell_differentiation_states"), "stem_cell_states")

    def test_doi_source_family_deduplication_is_case_insensitive(self):
        upper = train.canonical_source_family("doi:10.7554/eLife.72381")
        lower = train.canonical_source_family("doi:10.7554/elife.72381")
        self.assertEqual(upper, lower)
        self.assertEqual(lower, "doi:10.7554/elife.72381")

    def test_repository_config_declares_external_authority(self):
        config = json.loads(Path(__file__).with_name("config.json").read_text())
        self.assertEqual(config["authority"], "external_non_authoritative_weak_supervision")
        self.assertEqual(config["split"]["group_key"], "source_family_id")

    def test_committed_split_receipt_is_unique_and_complete(self):
        split_path = Path(__file__).with_name("splits.jsonl")
        rows = [json.loads(line) for line in split_path.read_text().splitlines()]
        families = [row["source_family_id"] for row in rows]
        self.assertEqual(len(rows), 7213)
        self.assertEqual(len(families), len(set(families)))
        self.assertEqual({row["split"] for row in rows}, {"train", "validation", "test"})

    def test_external_trainer_does_not_import_sealed_packages(self):
        for path in (MODULE_PATH, INFER_PATH):
            source = path.read_text()
            self.assertNotIn("aleph.learn", source)
            self.assertNotIn("aleph.represent", source)

    def test_inference_is_deterministic(self):
        first = self.index.query("actin cortical tension membrane", top_k=3)
        second = self.index.query("actin cortical tension membrane", top_k=3)
        self.assertEqual(first, second)

    def test_empty_and_unknown_queries_refuse(self):
        with self.assertRaisesRegex(infer.QueryRefusal, "non-empty") as empty:
            self.index.query("  ")
        self.assertEqual(empty.exception.code, "empty_query")
        with self.assertRaisesRegex(infer.QueryRefusal, "No query token") as unknown:
            self.index.query("zzzxqv qqqwvx")
        self.assertEqual(unknown.exception.code, "unknown_query")

    def test_physical_parameter_and_authority_requests_refuse(self):
        for query in ("estimate physical parameter value for cortical tension", "give authoritative ground truth"):
            with self.assertRaises(infer.QueryRefusal) as caught:
                self.index.query(query)
            self.assertEqual(caught.exception.code, "physical_or_authority_claim_refused")

    def test_retrieval_source_families_are_unique(self):
        result = self.index.query("traction force microscopy focal adhesion", top_k=20)
        families = [row["source_family_id"] for row in result["retrieval"]]
        self.assertEqual(len(families), len(set(families)))

    def test_example_results_match_frozen_api(self):
        queries = {row["id"]: row["query"] for row in json.loads(Path(__file__).with_name("example_queries.json").read_text())}
        expected = json.loads(Path(__file__).with_name("example_results.json").read_text())
        for item in expected["results"]:
            actual = self.index.query(queries[item["id"]], top_k=3)
            self.assertEqual(actual["weak_domain_probabilities"][0]["domain"], item["top_weak_domain"]["domain"])
            self.assertAlmostEqual(actual["weak_domain_probabilities"][0]["probability"], item["top_weak_domain"]["probability"])
            self.assertEqual(
                [row["source_family_id"] for row in actual["retrieval"]],
                [row["source_family_id"] for row in item["retrieval"]],
            )
            for actual_row, expected_row in zip(actual["retrieval"], item["retrieval"]):
                self.assertAlmostEqual(actual_row["cosine_score"], expected_row["score"])


if __name__ == "__main__":
    unittest.main()
