# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from io import StringIO

import skbio
import numpy.testing as npt
import numpy as np

from skbio.tree import (to_skbio_treenode, from_skbio_treenode, parse_newick,
                        to_skbio_treearray)


class ConversionTests(TestCase):
    def setUp(self):
        self.tstr = "(((a:1,b:2.5)c:6,d:8,(e),(f,g,(h:1,i:2)j:1)k:1.2)l,m:2)r;"
        self.bp = parse_newick(self.tstr)
        self.sktn = skbio.TreeNode.read(StringIO(self.tstr))

    def test_to_skbio_treenode_with_edge_numbers(self):
        in_ = '((A:.01{0}, B:.01{1})D:.01{3}, C:.01{4}) {5};'
        obs = parse_newick(in_)
        obs_sk = to_skbio_treenode(obs)
        self.assertEqual(obs_sk.find('A').edge_num, 0)
        self.assertEqual(obs_sk.find('B').edge_num, 1)
        self.assertEqual(obs_sk.find('D').edge_num, 3)
        self.assertEqual(obs_sk.find('C').edge_num, 4)
        self.assertEqual(obs_sk.edge_num, 5)

    def test_to_skbio_treenode(self):
        obs = to_skbio_treenode(self.bp)
        for o, e in zip(obs.traverse(), self.sktn.traverse()):
            if e.length is None:
                self.assertEqual(o.length, None if e.is_root() else 0.0)
            else:
                self.assertEqual(o.length, e.length)
            self.assertEqual(o.name, e.name)

        self.assertEqual(obs.ascii_art(), self.sktn.ascii_art())

    def test_from_skbio_treenode(self):
        obs_bp = from_skbio_treenode(self.sktn)
        exp_bp = self.bp

        npt.assert_equal(obs_bp.data, exp_bp.data)
        for i in range(len(self.bp.data)):
            self.assertEqual(exp_bp.name(i), obs_bp.name(i))
            self.assertEqual(exp_bp.length(i), obs_bp.length(i))

    def test_to_array(self):
        t = parse_newick('(((a:1,b:2,c:3)x:4,(d:5)y:6)z:7,(e:8,f:9)z:10);')

        exp_child_index = np.array([[4, 0, 2], [5, 3, 3], [8, 4, 5], [9, 6, 7],
                                    [10, 8, 9]], dtype=np.uint32)
        exp_length = np.array([1, 2, 3, 5, 4, 6, 8, 9, 7, 10, 0.0],
                              dtype=np.double)
        exp_id_index = {0: True, 1: True, 2: True, 3: True, 4: False, 5: False,
                        6: True, 7: True, 8: False, 9: False, 10: False}
        exp_name = np.array(['a', 'b', 'c', 'd', 'x', 'y', 'e', 'f', 'z', 'z',
                            None])
        obs = to_skbio_treearray(t)

        obs_child_index = obs['child_index']
        obs_length = obs['length']
        obs_id_index = obs['id_index']
        obs_name = obs['name']

        npt.assert_equal(obs_child_index, exp_child_index)
        npt.assert_equal(obs_length, exp_length)
        self.assertEqual(obs_id_index.keys(), exp_id_index.keys())
        npt.assert_equal(obs_name, exp_name)

        for k in obs_id_index:
            self.assertEqual(obs_id_index[k].is_tip(), exp_id_index[k])


if __name__ == '__main__':
    main()
