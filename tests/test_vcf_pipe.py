import unittest
from unittest.mock import MagicMock, patch, ANY
import pandas as pd
import sys
from pathlib import Path
import subprocess

# Add the project root to sys.path
sys.path.append(".")
from yleaf.Yleaf import YleafPipeline
from yleaf.configuration import Configuration

class TestVcfPipe(unittest.TestCase):

    @patch('pandas.DataFrame.to_csv')
    @patch('yleaf.Yleaf.pd.read_csv')
    @patch('yleaf.Yleaf.subprocess.Popen')
    @patch('yleaf.Yleaf.YleafPipeline.safe_create_dir')
    @patch('yleaf.Yleaf.YleafPipeline.write_info_file')
    @patch('yleaf.Yleaf.Tree')
    def test_run_vcf_pipe_logic(self, mock_tree, mock_write_info, mock_create_dir, mock_popen, mock_read_csv, mock_to_csv):
        # Setup mocks
        args = MagicMock()
        args.force = True
        args.reads_treshold = 10
        args.base_majority = 90
        args.use_old = False
        args.force = True # Avoid input prompt

        base_out_folder = Path("tmp_out")
        sample_vcf_file = Path("sample.vcf.gz")
        path_markerfile = Path("markers.txt")

        # Mock marker file dataframe
        marker_df = pd.DataFrame({
            "chr": ["Y"], "marker_name": ["M1"], "haplogroup": ["R"],
            "pos": [100], "mutation": ["A->G"], "anc": ["A"], "der": ["G"]
        })

        # Mock pileup file dataframe (what would come from the pipe)
        # The code expects header=None, so columns are 0,1,2,3,4 initially
        pileup_df = pd.DataFrame({
            0: ["Y"], 1: ["100"], 2: ["A"], 3: ["G"], 4: ["10,5"]
        })

        # Configure read_csv to return marker_df first, then pileup_df
        def read_csv_side_effect(*args, **kwargs):
            src = args[0]
            if src == path_markerfile:
                return marker_df.copy()
            else:
                # This should be the pipe
                return pileup_df.copy()

        mock_read_csv.side_effect = read_csv_side_effect

        # Mock subprocess process context manager
        process_mock = MagicMock()
        process_mock.stdout = "PIPE_HANDLE" # Just a placeholder
        process_mock.returncode = 0

        # Popen context manager
        mock_popen.return_value.__enter__.return_value = process_mock

        # Mock open for stderr file and output file writing
        # Use a context manager mock for open
        mock_file_ctx = MagicMock()
        mock_file_obj = MagicMock()
        mock_file_ctx.__enter__.return_value = mock_file_obj

        with patch("builtins.open", return_value=mock_file_ctx):
             # Run the function
            config = Configuration()
            pipeline = YleafPipeline(config)
            pipeline.run_vcf(path_markerfile, base_out_folder, args, sample_vcf_file)

        # Verify subprocess.Popen was called correctly
        # stdout is PIPE, stderr is the mocked file object (result of __enter__)
        expected_cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n", str(sample_vcf_file)]

        # Check that it was called with our mock file for stderr
        mock_popen.assert_called_with(expected_cmd, stdout=subprocess.PIPE, stderr=mock_file_obj, text=True)

        # Verify pd.read_csv was called with the pipe handle
        found_call = False
        for call in mock_read_csv.call_args_list:
            if call[0][0] == "PIPE_HANDLE":
                found_call = True
                # verify other args if needed
                _, kwargs = call
                self.assertEqual(kwargs.get('dtype'), str)
                self.assertEqual(kwargs.get('header'), None)
                self.assertEqual(kwargs.get('sep'), "\t")
                break
        self.assertTrue(found_call, "pd.read_csv was not called with the process stdout pipe")

        # Verify process.wait() was called (not communicate, because we use the file)
        process_mock.wait.assert_called_once()

if __name__ == '__main__':
    unittest.main()
