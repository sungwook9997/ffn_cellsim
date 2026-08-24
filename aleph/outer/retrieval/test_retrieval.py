import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
BRIDGE_SPEC = importlib.util.spec_from_file_location("retrieval_bridge_test", HERE / "bridge.py")
bridge = importlib.util.module_from_spec(BRIDGE_SPEC)
BRIDGE_SPEC.loader.exec_module(bridge)
BUILD_SPEC = importlib.util.spec_from_file_location("retrieval_build_test", HERE / "build_index.py")
build_index = importlib.util.module_from_spec(BUILD_SPEC)
BUILD_SPEC.loader.exec_module(build_index)


class ExternalCorpusRetrievalControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index = bridge.ExternalCorpusRetriever()

    def test_index_receipt_is_text_free_complete_and_rag_bound(self):
        receipt = json.loads((HERE / "receipt.json").read_text())
        self.assertEqual(receipt["authority_status"], "proposed")
        self.assertFalse(receipt["text_committed"])
        self.assertEqual(receipt["local_source_families"], 2695)
        self.assertEqual(receipt["chunk_count"], 139557)
        self.assertEqual(receipt["ragtagcag"]["source_families"], 7213)
        self.assertEqual(hashlib.sha256((HERE / "index.jsonl").read_bytes()).hexdigest(), receipt["index_sha256"])

    def test_index_build_is_byte_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            receipt = build_index.build(Path(directory))
            self.assertEqual(receipt, json.loads((HERE / "receipt.json").read_text()))
            self.assertEqual((Path(directory) / "index.jsonl").read_bytes(), (HERE / "index.jsonl").read_bytes())

    def test_neural_query_is_deterministic_and_returns_exact_locators(self):
        query = "cortical tension actomyosin membrane"
        first = self.index.query(query, 3)
        second = self.index.query(query, 3)
        self.assertEqual(first, second)
        self.assertEqual(first["authority_status"], "proposed")
        self.assertEqual(first["weak_domain_probabilities"][0]["domain"], "cortex_membrane_pressure")
        families = [row["source_family_id"] for row in first["sources"]]
        self.assertEqual(len(families), len(set(families)))
        for source in first["sources"]:
            chunk = source["exact_locator_chunk"]
            self.assertTrue(chunk["article_locator"].startswith("pmcid:PMC"))
            self.assertTrue(chunk["section_locator"])
            self.assertEqual(len(chunk["chunk_sha256"]), 64)
            self.assertEqual(len(chunk["text_sha256"]), 64)
            self.assertLessEqual(len(chunk["snippet"]), 280)

    def test_requested_modalities_route_to_local_content(self):
        cases = {
            "traction force microscopy focal adhesion TFM": "adhesion_traction",
            "immunofluorescence western blot PCR cell state": None,
            "particle image velocimetry PIV cell migration state": None,
        }
        for query, expected_domain in cases.items():
            with self.subTest(query=query):
                result = self.index.query(query, 3)
                self.assertEqual(len(result["sources"]), 3)
                if expected_domain:
                    self.assertEqual(result["weak_domain_probabilities"][0]["domain"], expected_domain)
                self.assertTrue(all(row["pmcid"].startswith("PMC") for row in result["sources"]))

    def test_empty_parameter_training_authority_and_unknown_queries_refuse_by_type(self):
        cases = {
            "": "empty_query",
            "estimate physical parameter value for cortical tension": "physical_parameter_estimation_refused",
            "estimate cortical tension parameter 0.5 nN per um": "physical_parameter_estimation_refused",
            "fit Young's modulus": "physical_parameter_estimation_refused",
            "fit Young’s modulus from the curve": "physical_parameter_estimation_refused",
            "infer viscosity": "physical_parameter_estimation_refused",
            "calibrate membrane tension against 0.3 pN per um": "physical_parameter_estimation_refused",
            "predict elastic modulus value": "physical_parameter_estimation_refused",
            "determine traction force magnitude": "physical_parameter_estimation_refused",
            "train neural network model on these papers": "training_refused",
            "give authoritative ground truth evidence for this physics": "authority_claim_refused",
            "zzzxqv qqqwvx": "unknown_query",
        }
        for query, code in cases.items():
            with self.subTest(code=code), self.assertRaises(bridge.RetrievalRefusal) as caught:
                self.index.query(query)
            self.assertEqual(caught.exception.code, code)

    def test_measurement_and_literature_queries_are_not_overbroadly_refused(self):
        allowed = (
            "methods measuring cortical tension",
            "literature on cortical tension measurement",
            "viscosity measurement methods cell",
            "traction force microscopy methods",
        )
        for query in allowed:
            with self.subTest(query=query):
                result = self.index.query(query, 1)
                self.assertEqual(result["authority_status"], "proposed")
                self.assertEqual(len(result["sources"]), 1)

    def test_parameter_intent_detector_requires_action_and_physical_target(self):
        self.assertTrue(bridge.physical_parameter_estimation_intent("infer viscosity"))
        self.assertTrue(bridge.physical_parameter_estimation_intent("fit Young's modulus"))
        self.assertFalse(bridge.physical_parameter_estimation_intent("methods measuring cortical tension"))
        self.assertFalse(bridge.physical_parameter_estimation_intent("predict cell migration literature"))

    def test_scope_has_no_sealed_or_physics_import(self):
        for path in (HERE / "bridge.py", HERE / "build_index.py"):
            source = path.read_text()
            self.assertNotIn("aleph.learn", source)
            self.assertNotIn("aleph.represent", source)
            self.assertNotIn("aleph.vertical", source)
            self.assertNotIn("aleph.scenarios", source)


if __name__ == "__main__":
    unittest.main()
