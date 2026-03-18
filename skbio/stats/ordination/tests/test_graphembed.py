
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import unittest
from unittest import skipIf
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import graphembed
from skbio.util import get_package


ge = get_package("graphembed_rs", raise_error=False)


def _make_mock_ge(ids, dims=2):
    mock_ge = MagicMock()

    def fake_embed_sketching(filepath, decay, dim, nbiter, symetric, output, **kw):
        df = pd.DataFrame(np.ones((len(ids), dim)), index=ids)
        df.to_csv(output, sep="\t", header=False)

    mock_ge.embed_sketching.side_effect = fake_embed_sketching
    mock_ge.embed_hope_rank.return_value = np.ones((len(ids), dims))

    return mock_ge


class TestGraphEmbedMissingDependency(unittest.TestCase):

    def test_missing_import_raises(self):
        adj = np.array([[0, 1], [1, 0]], dtype=float)
        with patch.dict(sys.modules, {"graphembed_rs": None}):
            with self.assertRaises(ImportError):
                graphembed(adj)


class TestGraphEmbedWithMock(unittest.TestCase):

    def setUp(self):
        self.ids = ["A", "B", "C"]
        self.dense = np.array(
            [[0, 1, 0.5], [1, 0, 0.2], [0.5, 0.2, 0]], dtype=float
        )
        self.dims = 2

    def _mock_ge(self):
        return _make_mock_ge(self.ids, self.dims)

    def test_distance_matrix_input(self):
        dm = DistanceMatrix(self.dense, ids=self.ids)
        mock_ge = self._mock_ge()
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            res = graphembed(dm, method="sketching", dimensions=self.dims, nbiter=1)
        self.assertEqual(res.samples.shape, (3, self.dims))
        self.assertListEqual(list(res.samples.index), self.ids)

    def test_dataframe_input(self):
        df = pd.DataFrame(self.dense, index=self.ids, columns=self.ids)
        mock_ge = self._mock_ge()
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            res = graphembed(df, method="sketching", dimensions=self.dims, nbiter=1)
        self.assertEqual(res.samples.shape, (3, self.dims))
        self.assertListEqual(list(res.samples.index), self.ids)

    def test_ndarray_input(self):
        mock_ge = self._mock_ge()
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            res = graphembed(
                self.dense, method="sketching", dimensions=self.dims, nbiter=1
            )
        self.assertEqual(res.samples.shape, (3, self.dims))

    def test_sparse_input(self):
        sparse = csr_matrix(self.dense)
        mock_ge = self._mock_ge()
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            res = graphembed(
                sparse, method="sketching", dimensions=self.dims, nbiter=1
            )
        self.assertEqual(res.samples.shape, (3, self.dims))

    def test_sketching_method(self):
        mock_ge = self._mock_ge()
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            res = graphembed(
                self.dense, method="sketching", dimensions=self.dims, nbiter=1
            )
        self.assertEqual(res.short_method_name, "GraphEmbed")
        self.assertIn("Sketching", res.long_method_name)
        mock_ge.embed_sketching.assert_called_once()

    def test_hope_method_ndarray(self):
        coords = np.ones((3, self.dims))
        mock_ge = MagicMock()
        mock_ge.embed_hope_rank.return_value = coords
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            res = graphembed(
                self.dense, method="hope", dimensions=self.dims, nbiter=1
            )
        self.assertEqual(res.samples.shape, (3, self.dims))
        self.assertIn("Hope", res.long_method_name)

    def test_hope_method_dict(self):
        ids_int = ["0", "1", "2"]
        coords_dict = {i: np.ones(self.dims) for i in ids_int}
        mock_ge = MagicMock()
        mock_ge.embed_hope_rank.return_value = coords_dict
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            res = graphembed(
                self.dense, method="hope", dimensions=self.dims, nbiter=1
            )
        self.assertEqual(res.samples.shape, (3, self.dims))

    def test_hope_method_other_type(self):
        coords_list = [[1.0] * self.dims] * 3
        mock_ge = MagicMock()
        mock_ge.embed_hope_rank.return_value = coords_list
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            res = graphembed(
                self.dense, method="hope", dimensions=self.dims, nbiter=1
            )
        self.assertEqual(res.samples.shape, (3, self.dims))

    def test_invalid_method_raises(self):
        mock_ge = self._mock_ge()
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            with self.assertRaises(ValueError):
                graphembed(self.dense, method="invalid_method")

    def test_ordination_results_structure(self):
        mock_ge = self._mock_ge()
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            res = graphembed(
                self.dense, method="sketching", dimensions=self.dims, nbiter=1
            )
        self.assertEqual(len(res.eigvals), self.dims)
        self.assertEqual(len(res.proportion_explained), self.dims)

    def test_sketching_output_at_base_path(self):
        ids = self.ids
        dims = self.dims

        def fake_embed_sketching(filepath, decay, dim, nbiter, symetric, output, **kw):
            df = pd.DataFrame(np.ones((len(ids), dim)), index=ids)
            df.to_csv(output, sep="\t", header=False)

        mock_ge = MagicMock()
        mock_ge.embed_sketching.side_effect = fake_embed_sketching

        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            with patch("os.path.isfile", return_value=False):
                res = graphembed(
                    self.dense, method="sketching", dimensions=dims, nbiter=1
                )
        self.assertEqual(res.samples.shape, (3, dims))

    def test_sketching_no_output_file_raises(self):
        mock_ge = MagicMock()
        mock_ge.embed_sketching.return_value = None
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            with self.assertRaises(Exception):
                graphembed(
                    self.dense, method="sketching", dimensions=self.dims, nbiter=1
                )

    def test_hope_none_result_raises(self):
        mock_ge = MagicMock()
        mock_ge.embed_hope_rank.return_value = None
        with patch.dict(sys.modules, {"graphembed_rs": mock_ge}):
            with self.assertRaises(RuntimeError):
                graphembed(
                    self.dense, method="hope", dimensions=self.dims, nbiter=1
                )


@skipIf(ge is None, "graphembed_rs is not installed.")
class TestGraphEmbedLive(unittest.TestCase):

    def test_sketching_live(self):
        adj = np.array([[0, 1], [1, 0]], dtype=float)
        res = graphembed(adj, method="sketching", dimensions=2, nbiter=1)
        self.assertEqual(res.samples.shape[1], 2)

    def test_sparse_live(self):
        adj = csr_matrix([[0, 1.5], [1.5, 0]])
        res = graphembed(adj, method="hope", dimensions=2)
        self.assertEqual(res.samples.shape[1], 2)


if __name__ == "__main__":
    unittest.main()
