import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
BRIDGE_SPEC = importlib.util.spec_from_file_location("combined_retrieval_bridge_test", HERE / "bridge.py")
bridge = importlib.util.module_from_spec(BRIDGE_SPEC)
BRIDGE_SPEC.loader.exec_module(bridge)
BUILD_SPEC = importlib.util.spec_from_file_location("combined_retrieval_build_test", HERE / "build_index.py")
build_index = importlib.util.module_from_spec(BUILD_SPEC)
BUILD_SPEC.loader.exec_module(build_index)


class CombinedRetrievalControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index = bridge.CombinedExternalCorpusRetriever()

    def test_receipt_covers_both_passes_without_committed_text(self):
        receipt = json.loads((HERE / "receipt.json").read_text())
        self.assertEqual(receipt["authority_status"], "proposed")
        self.assertFalse(receipt["text_committed"])
        self.assertEqual(receipt["local_source_families"], 5675)
        self.assertEqual(receipt["model"]["neural_routable_source_families"], 5669)
        self.assertEqual(receipt["pass_source_families"], {"first": 2701, "second": 2974})
        self.assertEqual(receipt["chunk_count"], 284620)
        self.assertEqual(receipt["ragtagcag"]["combined_source_families"], 10208)
        self.assertEqual(hashlib.sha256((HERE / "index.jsonl").read_bytes()).hexdigest(), receipt["index_sha256"])

    def test_index_is_byte_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            receipt = build_index.build(Path(directory))
            self.assertEqual(receipt, json.loads((HERE / "receipt.json").read_text()))
            self.assertEqual((Path(directory) / "index.jsonl").read_bytes(), (HERE / "index.jsonl").read_bytes())

    def test_query_is_deterministic_and_cross_pass(self):
        query = "cortical tension actomyosin membrane"
        first = self.index.query(query, 3)
        second = self.index.query(query, 3)
        self.assertEqual(first, second)
        self.assertEqual(first["authority_status"], "proposed")
        self.assertEqual(first["proposed_coarse_domain_probabilities"][0]["proposed_coarse_domain"], "cortical_membrane_pressure")
        self.assertEqual({row["pass"] for row in first["sources"]}, {"first", "second"})
        families = [row["source_family_id"] for row in first["sources"]]
        self.assertEqual(len(families), len(set(families)))
        for source in first["sources"]:
            chunk = source["exact_locator_chunk"]
            self.assertTrue(source["pmcid"].startswith("PMC"))
            self.assertTrue(chunk["section_locator"])
            self.assertEqual(len(chunk["chunk_sha256"]), 64)
            self.assertEqual(len(chunk["text_sha256"]), 64)
            self.assertLessEqual(len(chunk["snippet"]), 280)

    def test_requested_modalities_have_digest_verified_local_results(self):
        cases = (
            "traction force microscopy focal adhesion TFM",
            "immunofluorescence western blot PCR cell state",
            "particle image velocimetry PIV cell migration state",
        )
        for query in cases:
            with self.subTest(query=query):
                result = self.index.query(query, 3)
                self.assertEqual(len(result["sources"]), 3)
                self.assertTrue(all(row["exact_locator_chunk"]["snippet"] for row in result["sources"]))

    def test_committed_query_hashes_are_reproducible(self):
        benchmark = json.loads((HERE / "benchmark.json").read_text())
        for name, query in benchmark["queries"].items():
            with self.subTest(name=name):
                result = self.index.query(query, 3)
                payload = json.dumps(result, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
                self.assertEqual(hashlib.sha256(payload).hexdigest(), benchmark["combined"]["queries"][name]["result_sha256"])

    def test_hardened_typed_refusals_are_preserved(self):
        cases = {
            "": "empty_query",
            "estimate cortical tension parameter 0.5 nN per um": "physical_parameter_estimation_refused",
            "fit Young's modulus": "physical_parameter_estimation_refused",
            "infer viscosity": "physical_parameter_estimation_refused",
            "calibrate membrane tension": "physical_parameter_estimation_refused",
            "predict elastic modulus value": "physical_parameter_estimation_refused",
            "train neural network model on these papers": "training_refused",
            "give authoritative ground truth evidence for this physics": "authority_claim_refused",
            "zzzxqv qqqwvx": "unknown_query",
        }
        for query, expected in cases.items():
            with self.subTest(expected=expected), self.assertRaises(bridge.RetrievalRefusal) as caught:
                self.index.query(query)
            self.assertEqual(caught.exception.code, expected)

    def test_measurement_literature_queries_remain_allowed(self):
        for query in ("methods measuring cortical tension", "viscosity measurement methods cell", "traction force microscopy methods"):
            with self.subTest(query=query):
                self.assertEqual(len(self.index.query(query, 1)["sources"]), 1)

    def test_combined_scope_does_not_import_sealed_or_physics_modules(self):
        for path in (HERE / "bridge.py", HERE / "build_index.py"):
            source = path.read_text()
            self.assertNotIn("aleph.learn", source)
            self.assertNotIn("aleph.represent", source)
            self.assertNotIn("aleph.vertical", source)
            self.assertNotIn("aleph.scenarios", source)


if __name__ == "__main__":
    unittest.main()
