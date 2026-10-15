"""Run with: python -m unittest discover -s tests -p 'test_selected_qscores.py'."""
import csv
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "src/new/plot_selected_qscores.py"
spec = importlib.util.spec_from_file_location("selected_qscores", SCRIPT)
scores = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scores)


class ExactCommunityMatchingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.gmt = self.root / "communities.gmt"
        self.gmt.write_text("final_A\t2 genes\tg1\tg2\nfinal_B\t1 genes\tg4\n", encoding="utf-8")
        self.csv = self.root / "qscore_tar.csv"

    def write_scores(self, rows):
        with self.csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["tar", "community", "qscore", "resolution"])
            writer.writerows(rows)

    def test_exact_membership_not_ids_and_resolution_is_retained(self):
        self.write_scores([
            ["g2", 7, .2, 1], ["g1", 7, .3, 1],
            ["g1", 9, .4, 2], ["g2", 9, "NaN", 2],
            # Same numeric ID, but a superset: must not match.
            ["g1", 7, .1, 3], ["g2", 7, .1, 3], ["g3", 7, .3, 3],
            # A subset is not a match either.
            ["g1", 8, .5, 4],
        ])
        communities = scores.read_gmt(self.gmt)
        matches = scores.exact_matches(communities, self.csv)
        self.assertEqual(matches, {(1., "7"): ["final_A"], (2., "9"): ["final_A"]})
        mapped, coverage = self.root / "mapped.csv", self.root / "coverage.csv"
        distributions = scores.write_mapped_scores(communities, self.csv, matches, mapped, coverage)
        self.assertTrue(math.isnan(distributions["final_A"][(2., "9")][1]))
        with mapped.open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 4)
        self.assertEqual({r["selected_community"] for r in rows}, {"final_A"})
        self.assertEqual({r["condor_community"] for r in rows}, {"7", "9"})
        with coverage.open() as handle:
            report = {r["selected_community"]: r for r in csv.DictReader(handle)}
        self.assertEqual(report["final_A"]["n_undefined"], "1")
        self.assertEqual(report["final_B"]["status"], "no_exact_match")

    def test_duplicate_score_is_rejected(self):
        self.write_scores([["g1", 0, .2, 1], ["g1", 0, .2, 1]])
        with self.assertRaisesRegex(ValueError, "Duplicate qscore"):
            scores.exact_matches(scores.read_gmt(self.gmt), self.csv)

    def test_identical_final_sets_can_have_multiple_ids(self):
        self.gmt.write_text("A\tdescription\tg1\tg2\nB\tdescription\tg2\tg1\n")
        self.write_scores([["g1", "001", .25, 1], ["g2", "001", .25, 1]])
        self.assertEqual(scores.exact_matches(scores.read_gmt(self.gmt), self.csv),
                         {(1., "001"): ["A", "B"]})

    def test_empty_selection_and_no_matches_still_write_headers(self):
        self.gmt.write_text("")
        self.write_scores([["g1", 0, .5, 1]])
        communities = scores.read_gmt(self.gmt)
        matches = scores.exact_matches(communities, self.csv)
        mapped, coverage = self.root / "mapped.csv", self.root / "coverage.csv"
        scores.write_mapped_scores(communities, self.csv, matches, mapped, coverage)
        with mapped.open() as handle:
            self.assertEqual(list(csv.DictReader(handle)), [])
        self.assertIn("status", coverage.read_text())

    def test_missing_header_and_infinite_scores_fail_clearly(self):
        self.csv.write_text("gene,score\ng1,0.5\n")
        with self.assertRaisesRegex(ValueError, "must contain"):
            list(scores.score_rows(self.csv))
        self.write_scores([["g1", 0, "inf", 1]])
        with self.assertRaisesRegex(ValueError, "Infinite"):
            list(scores.score_rows(self.csv))


if __name__ == "__main__":
    unittest.main()
