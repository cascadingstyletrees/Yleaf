import unittest
from collections import Counter
from yleaf.Yleaf import get_frequencies

class TestYleaf(unittest.TestCase):
    def test_basic_counts(self):
        # Test basic base counting
        seq = "AACTTGGG"
        counts = get_frequencies(seq)
        self.assertEqual(counts["A"], 2)
        self.assertEqual(counts["C"], 1)
        self.assertEqual(counts["T"], 2)
        self.assertEqual(counts["G"], 3)
        self.assertEqual(counts["+"], 0)
        self.assertEqual(counts["-"], 0)

    def test_start_end_markers(self):
        # Test removal of start (^.) and end ($) markers
        seq = "^^!A$G^F$"
        # ^! should be removed
        # ^F should be removed
        # $ should be removed
        # Remaining: A, G
        counts = get_frequencies(seq)
        self.assertEqual(counts["A"], 1)
        self.assertEqual(counts["G"], 1)

    def test_indels(self):
        # Test insertions and deletions
        # +1A means insertion of A
        # -2GT means deletion of GT
        seq = "A+1TGC-2GTA"
        # A -> A
        # +1T -> Indel +, skip T
        # G -> G
        # C -> C
        # -2GT -> Indel -, skip GT
        # A -> A

        # Expected: A: 2, G: 1, C: 1, +: 1, -: 1
        counts = get_frequencies(seq)
        self.assertEqual(counts["A"], 2)
        self.assertEqual(counts["G"], 1)
        self.assertEqual(counts["C"], 1)
        self.assertEqual(counts["+"], 1)
        self.assertEqual(counts["-"], 1)

    def test_indel_skipping_bug(self):
        # This tests the specific bug fix where the indel sequence was not properly skipped.
        # Example: +2TT
        # If logic is correct: counts +: 1, T: 0.

        seq = "A+2TTG"
        # A -> A (1)
        # +2TT -> + (1), skip TT
        # G -> G (1)

        counts = get_frequencies(seq)
        self.assertEqual(counts["A"], 1)
        self.assertEqual(counts["G"], 1)
        self.assertEqual(counts["T"], 0) # Should be 0 if skipped correctly
        self.assertEqual(counts["+"], 1)

    def test_case_sensitivity(self):
        seq = "acgt"
        counts = get_frequencies(seq)
        self.assertEqual(counts["A"], 1)
        self.assertEqual(counts["C"], 1)
        self.assertEqual(counts["G"], 1)
        self.assertEqual(counts["T"], 1)

    def test_mixed_pileup(self):
        # Complex string
        seq = "^^!A$A+1TA-1G"
        # ^^! -> remove ^!
        # A -> A
        # $ -> remove
        # A -> A
        # +1T -> + (skip T)
        # A -> A
        # -1G -> - (skip G)

        # Total: A: 3, +: 1, -: 1
        counts = get_frequencies(seq)
        self.assertEqual(counts["A"], 3)
        self.assertEqual(counts["+"], 1)
        self.assertEqual(counts["-"], 1)

if __name__ == '__main__':
    unittest.main()
