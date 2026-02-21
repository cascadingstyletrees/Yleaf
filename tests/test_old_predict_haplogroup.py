
import unittest
import pandas as pd
import numpy as np
from yleaf.old_predict_haplogroup import calc_score_one

class TestOldPredictHaplogroup(unittest.TestCase):
    def test_calc_score_one(self):
        # Create dummy dataframes
        # df_intermediate: columns [0, 1] -> haplogroup, state
        data_int = {
            0: ["H1", "H2", "H3", "H4"],
            1: ["A", "D", "A/D", "D"]
        }
        df_intermediate = pd.DataFrame(data_int)

        # df_haplogroup: columns ["haplogroup", "state"]
        data_hap = {
            "haplogroup": ["H1", "H2", "H3", "H5"],
            "state": ["A", "A", "D", "D"]
        }
        df_haplogroup = pd.DataFrame(data_hap)

        # Expected calculation:
        # H1: present in haplogroup. state in int is 'A', in hap is 'A'. Match. +1 correct, +1 total.
        # H2: present in haplogroup. state in int is 'D', in hap is 'A'. Mismatch. +0 correct, +1 total.
        # H3: present in haplogroup. state in int is 'A/D'. Match (wildcard). +1 correct, +1 total.
        # H4: not present in haplogroup. Ignored.
        # H5: in haplogroup but not intermediate. Ignored.

        # Total correct: 1 (H1) + 1 (H3) = 2
        # Total total: 1 (H1) + 1 (H2) + 1 (H3) = 3
        # Score: 2/3 = 0.6666... -> rounded to 3 decimals = 0.667

        score = calc_score_one(df_intermediate, df_haplogroup)
        self.assertEqual(score, 0.667)

    def test_calc_score_one_empty(self):
        df_intermediate = pd.DataFrame()
        df_haplogroup = pd.DataFrame({"haplogroup": [], "state": []})
        score = calc_score_one(df_intermediate, df_haplogroup)
        self.assertEqual(score, 0.0)

    def test_calc_score_one_zero_division(self):
        # No matches found -> total = 0
        data_int = {0: ["H1"], 1: ["A"]}
        df_intermediate = pd.DataFrame(data_int)
        df_haplogroup = pd.DataFrame({"haplogroup": ["H2"], "state": ["A"]})

        score = calc_score_one(df_intermediate, df_haplogroup)
        self.assertEqual(score, 0.0)

if __name__ == "__main__":
    unittest.main()
