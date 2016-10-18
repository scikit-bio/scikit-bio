# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio import TreeNode, DistanceMatrix
from skbio.diversity import beta_diversity, block_beta_diversity
from skbio.diversity._block import (_block_party, _generate_id_blocks,
                                    _pairs_to_compute, _block_compute,
                                    _block_kwargs, _map, _reduce)


class ParallelBetaDiversity(TestCase):
    def setUp(self):
        self.table1 = [[1, 5],
                       [2, 3],
                       [0, 1]]
        self.sids1 = list('ABC')
        self.tree1 = TreeNode.read([
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'])
        self.oids1 = ['O1', 'O2']

    def test_block_kwargs(self):
        kws = {'ids': [1, 2, 3, 4, 5], 'foo': 'bar', 'k': 2}
        exp = [{'row_ids': np.array((0, 1)),
                'col_ids': np.array((0, 1)),
                'id_pairs': [(0, 1)],
                'ids': [1, 2, 3, 4, 5]},

               {'row_ids': np.array((0, 1)),
                'col_ids': np.array((2, 3)),
                'id_pairs': [(0, 2), (0, 3), (1, 2), (1, 3)],
                'ids': [1, 2, 3, 4, 5]},

               {'row_ids': np.array((0, 1)),
                'col_ids': np.array((4,)),
                'id_pairs': [(0, 4), (1, 4)],
                'ids': [1, 2, 3, 4, 5]},

               {'row_ids': np.array((2, 3)),
                'col_ids': np.array((2, 3)),
                'id_pairs': [(2, 3), ],
                'ids': [1, 2, 3, 4, 5]},

               {'row_ids': np.array((2, 3)),
                'col_ids': np.array((4,)),
                'id_pairs': [(2, 4), (3, 4)],
                'ids': [1, 2, 3, 4, 5]}]

        obs = list(_block_kwargs(**kws))
        npt.assert_equal(obs, exp)

    def test_block_compute(self):
        def mock_metric(u, v):
            return (u + v).sum()

        counts = np.array([[0, 1, 2, 3, 4, 5],
                           [1, 2, 3, 4, 5, 0],
                           [2, 3, 4, 5, 0, 1],
                           [10, 2, 3, 6, 8, 2],
                           [9, 9, 2, 2, 3, 4]])

        kwargs = {'metric': mock_metric,
                  'counts': counts,
                  'row_ids': np.array((2, 3)),
                  'col_ids': np.array((4, )),
                  'id_pairs': [(2, 4), (3, 4)],
                  'ids': [1, 2, 3, 4, 5]}

        exp = DistanceMatrix(np.array([[0, 0, 44],
                                       [0, 0, 60],
                                       [44, 60, 0]]), (2, 3, 4))

        obs = _block_compute(**kwargs)
        npt.assert_equal(obs.data, exp.data)
        self.assertEqual(obs.ids, exp.ids)

    def test_map(self):
        def func(a, b, c=5):
            return a + b + c

        kwargs = [{'a': 0, 'b': 1, 'c': 0},
                  {'a': 2, 'b': 3}]
        exp = [1, 10]
        obs = list(_map(func, kwargs))
        self.assertEqual(obs, exp)

    def test_reduce(self):
        dm1 = DistanceMatrix(np.array([[0, 0, 44],
                                       [0, 0, 60],
                                       [44, 60, 0]]), (2, 3, 4))
        dm2 = DistanceMatrix(np.array([[0, 123],
                                       [123, 0]]), (1, 5))
        dm3 = DistanceMatrix(np.array([[0, 1, 2, 3],
                                       [1, 0, 4, 5],
                                       [2, 4, 0, 6],
                                       [3, 5, 6, 0]]), (0, 3, 4, 5))
        exp = DistanceMatrix(np.array([[0, 0, 0, 1, 2, 3],
                                       [0, 0, 0, 0, 0, 123],
                                       [0, 0, 0, 0, 44, 0],
                                       [1, 0, 0, 0, 64, 5],
                                       [2, 0, 44, 64, 0, 6],
                                       [3, 123, 0, 5, 6, 0]]), list(range(6)))

        obs = _reduce([dm1, dm2, dm3])
        npt.assert_equal(obs.data, exp.data)
        self.assertEqual(obs.ids, exp.ids)

    def test_block_beta_diversity(self):
        exp = beta_diversity('unweighted_unifrac', self.table1, self.sids1,
                             tree=self.tree1, otu_ids=self.oids1)
        obs = block_beta_diversity('unweighted_unifrac', self.table1,
                                   self.sids1, otu_ids=self.oids1,
                                   tree=self.tree1, k=2)
        npt.assert_equal(obs.data, exp.data)
        self.assertEqual(obs.ids, exp.ids)

    def test_generate_id_blocks(self):
        ids = [1, 2, 3, 4, 5]
        exp = [(np.array((0, 1)), np.array((0, 1))),
               (np.array((0, 1)), np.array((2, 3))),
               (np.array((0, 1)), np.array((4,))),
               (np.array((2, 3)), np.array((2, 3))),
               (np.array((2, 3)), np.array((4,))),
               (np.array((4,)),   np.array((4,)))]

        obs = list(_generate_id_blocks(ids, 2))
        npt.assert_equal(obs, exp)

    def test_block_party_notree(self):
        counts = np.arange(15).reshape(5, 3)
        exp = [{'counts': np.array([[0, 1, 2], [3, 4, 5]]),
                'ids': np.array([0, 1])},
               {'counts': np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8],
                                    [9, 10, 11]]),
                'ids': np.array([0, 1, 2, 3])},
               {'counts': np.array([[0, 1, 2], [3, 4, 5], [12, 13, 14]]),
                'ids': np.array([0, 1, 4])},
               {'counts': np.array([[6, 7, 8], [9, 10, 11]]),
                'ids': np.array([2, 3])},
               {'counts': np.array([[6, 7, 8], [9, 10, 11], [12, 13, 14]]),
                'ids': np.array([2, 3, 4])},
               {'counts': np.array([[12, 13, 14]]), 'ids': np.array([4])}]
        obs = [_block_party(counts, rids, cids) for rids, cids in
               _generate_id_blocks(list(range(5)), 2)]
        npt.assert_equal(obs, exp)

    def test_block_party_tree(self):
        counts = np.array([[1, 1, 1],
                           [1, 0, 1],
                           [1, 0, 1],
                           [0, 0, 1],
                           [0, 1, 1]])
        tree = TreeNode.read(['(a:1,b:2,c:3);'])
        otu_ids = ['a', 'b', 'c']

        kw = {'tree': tree, 'otu_ids': otu_ids}
        kw_no_a = {'tree': tree.shear(['b', 'c']), 'otu_ids': ['b', 'c']}
        kw_no_b = {'tree': tree.shear(['a', 'c']), 'otu_ids': ['a', 'c']}

        # python >= 3.5 supports {foo: bar, **baz}
        exp = [dict(counts=np.array([[1, 1, 1], [1, 0, 1]]), **kw),
               dict(counts=np.array([[1, 1, 1], [1, 0, 1], [1, 0, 1],
                                     [0, 0, 1]]), **kw),
               dict(counts=np.array([[1, 1, 1], [1, 0, 1], [0, 1, 1]]), **kw),
               dict(counts=np.array([[1, 1], [0, 1]]), **kw_no_b),
               dict(counts=np.array([[1, 0, 1], [0, 0, 1], [0, 1, 1]]), **kw),
               dict(counts=np.array([[1, 1]]), **kw_no_a)]

        obs = [_block_party(counts, rids, cids, **kw) for rids, cids in
               _generate_id_blocks(list(range(5)), 2)]

        for okw, ekw in zip(obs, exp):
            npt.assert_equal(okw['counts'], ekw['counts'])
            npt.assert_equal(okw['otu_ids'], ekw['otu_ids'])
            self.assertEqual(str(okw['tree']), str(ekw['tree']))

    def test_pairs_to_compute_rids_are_cids(self):
        rids = np.array([0, 1, 2, 10])
        cids = rids
        exp = [(0, 1), (0, 2), (0, 10), (1, 2), (1, 10), (2, 10)]
        self.assertEqual(_pairs_to_compute(rids, cids), exp)

    def test_pairs_to_compute_rids_are_not_cids(self):
        rids = np.array([0, 1, 2])
        cids = np.array([3, 4, 5])
        exp = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4),
               (2, 5)]
        self.assertEqual(_pairs_to_compute(rids, cids), exp)

    def test_pairs_to_compute_rids_overlap_cids(self):
        rids = np.array([0, 1, 2])
        cids = np.array([0, 10, 20])
        with self.assertRaises(ValueError):
            _pairs_to_compute(rids, cids)


if __name__ == "__main__":
    main()
