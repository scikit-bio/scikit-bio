from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from future.builtins import zip

from skbio.core.tree import TreeNode
from skbio.core.tree.util import shuffle


class UtilTests(TestCase):
    def setUp(self):
        self.simple_tree = TreeNode.from_newick("((a,b),(c,d))")
        self.rev_f = lambda items: items[::-1]
        self.rotate_f = lambda items: [items[-1]] + items[:-1]
        self.complex_tree = TreeNode.from_newick("(((a,b)int1,(x,y,(w,z)int2,"
                                                 "(c,d)int3)int4),(e,f)int5);")

    def test_shuffle_n(self):
        exp = ["((a,b),(d,c));",
               "((a,b),(c,d));",
               "((a,b),(d,c));",
               "((a,b),(c,d));",
               "((a,b),(d,c));"]

        obs_g = shuffle(self.simple_tree, n=2, shuffle_f=self.rev_f)
        obs = [next(obs_g) for i in range(5)]
        self.assertEqual([o.to_newick() for o in obs], exp)
        self.assertFalse(id(self.simple_tree) in [id(o) for o in obs])

        exp = ["((d,c),(b,a));",
               "((a,b),(c,d));",
               "((d,c),(b,a));",
               "((a,b),(c,d));"]
        obs_g = shuffle(self.simple_tree, shuffle_f=self.rev_f)
        obs = [next(obs_g) for i in range(4)]
        self.assertEqual([o.to_newick() for o in obs], exp)

    def test_shuffle_complex(self):
        exp = ["(((a,b)int1,(x,y,(w,z)int2,(f,e)int3)int4),(d,c)int5);",
               "(((a,b)int1,(x,y,(w,z)int2,(c,d)int3)int4),(e,f)int5);",
               "(((a,b)int1,(x,y,(w,z)int2,(f,e)int3)int4),(d,c)int5);",
               "(((a,b)int1,(x,y,(w,z)int2,(c,d)int3)int4),(e,f)int5);"]

        obs_g = shuffle(self.complex_tree, shuffle_f=self.rev_f,
                        names=['c', 'd', 'e', 'f'], inplace=True)

        for e, o in zip(exp, obs_g):
            self.assertEqual(o.to_newick(), e)
            self.assertEqual(id(o), id(self.complex_tree))

    def test_shuffle_names(self):
        exp = ["((c,a),(b,d));",
               "((b,c),(a,d));",
               "((a,b),(c,d));",
               "((c,a),(b,d));"]

        obs_g = shuffle(self.simple_tree, names=['a', 'b', 'c'],
                        shuffle_f=self.rotate_f)
        obs = [next(obs_g) for i in range(4)]
        self.assertEqual([o.to_newick() for o in obs], exp)

    def test_shuffle_raises(self):
        with self.assertRaises(ValueError):
            next(shuffle(self.simple_tree, n=1))

        with self.assertRaises(ValueError):
            next(shuffle(self.simple_tree, n=5, names=['a', 'b']))

        with self.assertRaises(ValueError):
            next(shuffle(self.simple_tree, names=['x', 'y']))


if __name__ == '__main__':
    main()
