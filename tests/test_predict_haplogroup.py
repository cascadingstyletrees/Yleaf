
import unittest
from yleaf import predict_haplogroup
from yleaf.predict_haplogroup import get_qc1_score, HgMarkersLinker, read_backbone_groups

class TestPredictHaplogroup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Load the backbone groups once for the test class
        read_backbone_groups()

    def test_get_qc1_score_logic(self):
        # Setup
        path = ["R1b", "R1", "R"]

        haplotype_dict = {}

        # R1b_int.txt has "R1" -> "D".
        # We add R1 as DERIVED to haplotype_dict.
        linker = HgMarkersLinker()
        linker.add("marker1", HgMarkersLinker.DERIVED)
        haplotype_dict["R1"] = linker

        # We also need to make sure we don't hit the per-sample cache across tests if we had multiple
        predict_haplogroup.QC1_SCORE_CACHE = {}

        score = get_qc1_score(path, haplotype_dict)

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

        predict_haplogroup.QC1_SCORE_CACHE = {}
        score = get_qc1_score(path, haplotype_dict)

        # Matching: 0. Total: 1. Score: 0.0.
        self.assertEqual(score, 0.0)

if __name__ == '__main__':
    unittest.main()
