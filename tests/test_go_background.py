import importlib.util
from io import StringIO
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location(
    "background", Path(__file__).resolve().parents[1] / "src/new/derive_go_background.py")
background = importlib.util.module_from_spec(spec)
spec.loader.exec_module(background)


class BackgroundTests(unittest.TestCase):
    def test_all_targets_including_zero_and_negative_weight_edges(self):
        network = StringIO("TF1,B,2\nTF2,A,0\nTF3,C,-5\nTF1,A,1\n")
        self.assertEqual(background.target_genes(network), ["A", "B", "C"])

    def test_tab_delimiter_and_blank_lines(self):
        self.assertEqual(background.target_genes(StringIO("\nTF1\tA\t1\n"), "\t"), ["A"])

    def test_invalid_and_empty_networks(self):
        for text in ("", "TF1\n", "TF1,,1\n", "TF1,NA,1\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                background.target_genes(StringIO(text))

    def test_output_matches_enrichment_column_contract(self):
        output = StringIO()
        background.write_background(["A", "B"], output)
        self.assertEqual(output.getvalue(), "1\tA\n2\tB\n")


if __name__ == "__main__":
    unittest.main()
