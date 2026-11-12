"""Regression checks for signed PANDA input and comparison scores."""
import ast
from io import StringIO
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


# Exercise numerical helpers without requiring the Leiden/igraph runtime.
source = Path(__file__).resolve().parents[1] / "src/new/compare_hedgehog_leiden.py"
tree = ast.parse(source.read_text())
helpers = ast.Module(body=[node for node in tree.body if isinstance(node, ast.FunctionDef)
                          and node.name in {"load_network", "community_scores"}], type_ignores=[])
namespace = {"np": np, "pd": pd}
exec(compile(helpers, str(source), "exec"), namespace)


class Graph:
    def __init__(self, weight):
        self.es = {"weight": [weight]}
        self.vs = {"name": ["A", "B"]}
        self.weight = weight

    def strength(self, weights):
        return [self.weight, self.weight]

    def get_adjacency_sparse(self, attribute):
        return np.array([[0., self.weight], [self.weight, 0.]])


class SignedWeightsTests(unittest.TestCase):
    def test_negative_and_zero_sum_networks_are_accepted(self):
        for weights in ((-2, -1), (-1, 1)):
            frame = namespace["load_network"](
                StringIO(f"TF,A,{weights[0]}\nTF,B,{weights[1]}\n"), ",")
            self.assertEqual(frame.weight.tolist(), list(weights))

    def test_nonfinite_weights_still_fail(self):
        for weight in ("nan", "inf", "-inf"):
            with self.subTest(weight=weight), self.assertRaises(ValueError):
                namespace["load_network"](StringIO(f"TF,A,{weight}\n"), ",")

    def test_signed_scores_use_cpm_and_mark_modularity_undefined(self):
        row, = namespace["community_scores"](Graph(-2), {"all": {"A", "B"}}, "leiden", .5)
        self.assertTrue(np.isnan(row["modularity_contribution"]))
        self.assertEqual(row["cpm_contribution"], -5)

    def test_positive_scores_preserve_standard_modularity(self):
        row, = namespace["community_scores"](Graph(2), {"all": {"A", "B"}}, "leiden", .5)
        self.assertEqual(row["modularity_contribution"], 0)
        self.assertEqual(row["cpm_contribution"], 3)


if __name__ == "__main__":
    unittest.main()
