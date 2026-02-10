
import unittest
from pathlib import Path
from yleaf.predict_haplogroup import HaplogroupPredictor, HgMarkersLinker
from yleaf.configuration import Configuration
from yleaf.tree import Tree

class TestPredictHaplogroup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # We need a dummy configuration and tree for testing
        cls.config = Configuration()
        # Mock tree or use None if methods don't strictly require it for the logic being tested
        # But get_most_likely_haplotype uses tree.get(k).depth
        # So we should probably mock the tree or rely on file existence.
        # Since we are refactoring, we assume the environment has the data files (as per original tests).

        # We can construct a minimal tree object without loading file if Tree class supports it,
        # or just load the real tree.
        try:
            tree = Tree(cls.config.data_folder / "hg_prediction_tables" / "tree.json")
        except:
            tree = None # Should handle if file not found in test env

        cls.predictor = HaplogroupPredictor(cls.config, tree)

        # We need to manually populate backbone groups if file IO fails or to ensure consistent state
        cls.predictor.backbone_groups = {"R", "R1", "R1b"}
        cls.predictor.main_haplo_groups = {"R", "R1", "R1b"}

        # Mock expected states for R1b
        cls.predictor.expected_states_cache = {
            "R1b": {
                "R1": {"D"},
                "R": {"D"}
            }
        }
        cls.predictor.qc1_score_cache = {}

    def test_get_qc1_score_logic(self):
        # Setup
        path = ["R1b", "R1", "R"]

        haplotype_dict = {}

        # R1b_int.txt has "R1" -> "D".
        # We add R1 as DERIVED to haplotype_dict.
        linker = HgMarkersLinker()
        linker.add("marker1", HgMarkersLinker.DERIVED)
        haplotype_dict["R1"] = linker

        # Reset cache
        self.predictor.qc1_score_cache = {}

        intermediate_states = {
            value: haplotype_dict[value]
            for value in self.predictor.backbone_groups
            if value in haplotype_dict
        }
        score = self.predictor.get_qc1_score(path, intermediate_states)

        # R1 is in intermediate_states. expected state is {D}. Actual state is D.
        # So score[0] (matching) should be 1.
        # score[1] (total) should be 1.
        # Result should be 1.0.

        self.assertEqual(score, 1.0)

    def test_get_qc1_score_miss(self):
        path = ["R1b", "R1", "R"]
        haplotype_dict = {}

        # R1 -> Ancestral. Expected D.
        linker = HgMarkersLinker()
        linker.add("marker1", HgMarkersLinker.ANCESTRAL)
        haplotype_dict["R1"] = linker

        self.predictor.qc1_score_cache = {}

        intermediate_states = {
            value: haplotype_dict[value]
            for value in self.predictor.backbone_groups
            if value in haplotype_dict
        }
        score = self.predictor.get_qc1_score(path, intermediate_states)

        # Matching: 0. Total: 1. Score: 0.0.
        self.assertEqual(score, 0.0)

if __name__ == '__main__':
    unittest.main()
