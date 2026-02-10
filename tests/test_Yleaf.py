import unittest
from yleaf.Yleaf import YleafPipeline
from yleaf.configuration import Configuration

class TestYleaf(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = Configuration()
        cls.pipeline = YleafPipeline(cls.config)

    def test_basic_counts(self):
        # Test basic base counting
        seq = "AACTTGGG"
        counts = self.pipeline.get_frequencies(seq)
        self.assertEqual(counts["A"], 2)
        self.assertEqual(counts["C"], 1)
        self.assertEqual(counts["T"], 2)
        self.assertEqual(counts["G"], 3)
        self.assertEqual(counts["+"], 0)
        self.assertEqual(counts["-"], 0)

    def test_start_end_markers(self):
        # Test removal of start (^.) and end ($) markers
        seq = "^^!A$G^F$"
        counts = self.pipeline.get_frequencies(seq)
        self.assertEqual(counts["A"], 1)
        self.assertEqual(counts["G"], 1)

    def test_indels(self):
        # Test insertions and deletions
        seq = "A+1TGC-2GTA"
        counts = self.pipeline.get_frequencies(seq)
        self.assertEqual(counts["A"], 2)
        self.assertEqual(counts["G"], 1)
        self.assertEqual(counts["C"], 1)
        self.assertEqual(counts["+"], 1)
        self.assertEqual(counts["-"], 1)

    def test_indel_skipping_bug(self):
        seq = "A+2TTG"
        counts = self.pipeline.get_frequencies(seq)
        self.assertEqual(counts["A"], 1)
        self.assertEqual(counts["G"], 1)
        self.assertEqual(counts["T"], 0)  # Should be 0 if skipped correctly
        self.assertEqual(counts["+"], 1)

    def test_case_sensitivity(self):
        seq = "acgt"
        counts = self.pipeline.get_frequencies(seq)
        self.assertEqual(counts["A"], 1)
        self.assertEqual(counts["C"], 1)
        self.assertEqual(counts["G"], 1)
        self.assertEqual(counts["T"], 1)

    def test_mixed_pileup(self):
        # Complex string
        seq = "^^!A$A+1TA-1G"
        counts = self.pipeline.get_frequencies(seq)
        self.assertEqual(counts["A"], 3)
        self.assertEqual(counts["+"], 1)
        self.assertEqual(counts["-"], 1)


if __name__ == "__main__":
    unittest.main()
