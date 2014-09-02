# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils.six import StringIO

import unittest

from skbio.tree import TreeNode
from skbio.io.newick import (_newick_to_tree_node, _tree_node_to_newick,
                             _newick_sniffer)


class TestNewick(unittest.TestCase):
    def _invert_tree(self, tree):
        tree.children = [self._invert_tree(child) for child
                         in reversed(tree.children)]
        return tree

    def _is_node_equal(self, n1, n2):
        self.assertEqual(n1.name, n2.name)
        self.assertEqual(n1.length, n2.length)
        self.assertEqual(len(n1.children), len(n2.children))
        return True

    def _is_equal(self, n1, n2):
        name = lambda x: (x.name, x.length, len(x.children))
        if self._is_node_equal(n1, n2):
            for c1, c2 in zip(sorted(n1.children, key=name),
                              sorted(n2.children, key=name)):
                self.assertTrue(c1.parent is n1)
                self.assertTrue(c2.parent is n2)
                if not self._is_equal(c1, c2):
                    return False
            return True
        return False

    def _setup_tree(self, kwargs_list):
        trees = []
        for kwargs in kwargs_list:
            trees.append(TreeNode(**kwargs))

        trees[4].extend([trees[2], trees[3]])
        trees[5].extend([trees[0], trees[1], trees[4]])

        return trees[5]

    def _setup_linked_list(self, kwargs_list):
        last_node = None
        for idx, kwargs in enumerate(kwargs_list):
            new_node = TreeNode(**kwargs)

            if last_node is not None:
                new_node.append(last_node)
            last_node = new_node
        return last_node

    def setUp(self):
        tree_blank = (self._setup_tree([
            {}, {}, {}, {}, {}, {}
        ]), [
            "(,,(,));",
            "(,(,),);",
            "((,),,);",
            "   ((,[ this is a comment ])      ,    ,   )    ;  ",
            "((,[ i_can_do_this[0] or escape unmatched '[ ]),[more words],);",
        ])

        tree_leaves_named = (self._setup_tree([
            {'name': 'a_'},
            {'name': 'b'},
            {'name': 'c'},
            {'name': 'd'},
            {},
            {}
        ]), [
            "('a_',b,(c,d));",
            "(b,(c,d),'a_');",
            "(b\n,'a_'\n  ,(d \t,c) )  ;",
        ])

        tree_all_named = (self._setup_tree([
            {'name': 'a'},
            {'name': 'b'},
            {'name': 'c'},
            {'name': '[whaaat!\']'},
            {'name': 'e'},
            {'name': 'f'}
        ]), [
            "(a,b,(c,'[whaaat!'']')e)f;",
            "(b,(c,'[whaaat!'']')e,a)f;",
            "(b,[comment] \na,('[whaaat!'']',c)e)f;",
        ])

        tree_all_but_root_distances = (self._setup_tree([
            {'length': 0.1},
            {'length': 0.2},
            {'length': 0.3},
            {'length': 0.4},
            {'length': 0.5},
            {}
        ]), [
            "(:0.1,:0.2,(:0.3,:0.4):0.5);",
            "(:0.2,(:0.3,:0.4):0.5,:0.1);",
            "(:0.2,:0.1,(:0.4,:0.3):0.5);",
        ])

        tree_all_distances = (self._setup_tree([
            {'length': 0.1},
            {'length': 0.2},
            {'length': 0.3},
            {'length': 0.4},
            {'length': 0.5},
            {'length': 0.0}
        ]), [
            "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;",
            "(:0.2,(:0.3,:0.4):0.5,:0.1):0.0;",
            "(:0.2,\n:0.1,(:0.4,\n:0.3):0.5)\n:0.0;",
        ])

        tree_all_leaves_named_with_distances = (self._setup_tree([
            {'name': 'a', 'length': 0.1},
            {'name': 'b_a\'', 'length': 0.2},
            {'name': 'c', 'length': 0.3},
            {'name': 'de d', 'length': 0.4},
            {'length': 0.5},
            {'length': 0.0}
        ]), [
            "(a:0.1,'b_a''':0.2,(c:0.3,de_d:0.4):0.5):0.0;",
            "('b_a''':0.2,(c:0.3,'de d':0.4):0.5,a:0.1):0.0;",
            "('b_a''':0.2,a:0.1,('de d'[why not]:0.4,c:0.3):0.5):0.0;",
        ])

        tree_all_leaves_named_with_distances_no_root = (self._setup_tree([
            {'name': 'a', 'length': 0.1},
            {'name': 'b_a\'', 'length': 0.2},
            {'name': 'c', 'length': 0.3},
            {'name': 'de d', 'length': 0.4},
            {'length': 0.5},
            {}
        ]), [
            "(a:0.1,'b_a''':0.2,(c:0.3,de_d:0.4):0.5);",
            "('b_a''':0.2\n[comment ahoy]\n,(c:0.3,'de d':0.4):0.5,a:0.1);",
            "('b_a''':0.2,a:0.1,(de_d:0.4,c:0.3):0.5);"
        ])

        tree_all = (self._setup_tree([
            {'name': 'a', 'length': 0.1},
            {'name': 'b_a\'', 'length': 0.2},
            {'name': 'c', 'length': 0.3},
            {'name': 'de d', 'length': 0.4},
            {'name': 'e', 'length': 0.5},
            {'name': 'f', 'length': 0.0}
        ]), [
            "(a:0.1,'b_a''':0.2,(c:0.3,de_d:0.4)e:0.5)f:0.0;",
            "('b_a''':0.2,(c:0.3,de_d:0.4)e:0.5,a:0.1)f:0.0;",
            "((de_d:0.4, c:0.3)e:0.5, 'b_a''':0.2, a:0.1)f:0.0;"
        ])

        linked_list_blank = (self._setup_linked_list([
            {}, {}, {}, {}, {}
        ]), [
            "(((())));",
            "[(((())));](((())));",
            "[[(((())));](((())));](((())));"
        ])

        linked_list_named = (self._setup_linked_list([
            {'name': 'aaa'},
            {'name': 'b_a\''},
            {'name': 'c'},
            {'name': 'de d'},
            {'name': 'e'},
        ]), [
            "((((aaa)'b_a''')c)de_d)e;"
        ])

        inked_list_distances = (self._setup_linked_list([
            {'length': 0.4},
            {'length': 0.3},
            {'length': 0.2},
            {'length': 0.1},
            {'length': 0.0},
        ]), [
            "((((:0.4):0.3):0.2):0.1):0.0;",
            "((((:0.4)[not a label]:0.3):0.2):0.1):0.0;"
        ])

        linked_list_all = (self._setup_linked_list([
            {'name': 'a', 'length': 0.4},
            {'name': 'b_a\'', 'length': 0.3},
            {'name': 'c', 'length': 0.2},
            {'name': 'de d', 'length': 0.1},
            {'name': 'eee', 'length': 0.0},
        ]), [
            "((((a:0.4)'b_a''':0.3)c:0.2)de_d:0.1)eee:0.0;"
        ])

        single_empty = (TreeNode(), [";"])
        single_empty = (TreeNode(), ["[insightful comment about the root and"
                                     " it's properties];"])
        single_named = (TreeNode(name='athing'), ["a thing;"])
        single_distance = (TreeNode(length=200), [":200;"])
        single_all = (TreeNode(name='[a]', length=200), ["'[a]':200;"])
        empty = (TreeNode(), [""])

        self.trees_newick_lists = [
            tree_blank,
            tree_leaves_named,
            tree_all_named,
            tree_all_but_root_distances,
            tree_all_distances,
            tree_all_leaves_named_with_distances,
            tree_all_leaves_named_with_distances_no_root,
            tree_all,
            linked_list_blank,
            linked_list_named,
            inked_list_distances,
            linked_list_all,
            single_empty,
            single_named,
            single_distance,
            single_all,
            empty
        ]

    def test_read_valid_files(self):
        for tree, newicks in self.trees_newick_lists:
            for newick in newicks:
                fh = StringIO(newick)
                read_tree = _newick_to_tree_node(fh)

                self.assertTrue(self._is_equal(tree, read_tree))

                fh.close()

    def test_profile_new(self):
        import cProfile
        cProfile.run("_newick_to_tree_node('/home/evan/Downloads/97_otus.tree')")
        cProfile.run("_newick_to_tree_node('/home/evan/Downloads/97_otus.tree')")
        cProfile.run("_newick_to_tree_node('/home/evan/Downloads/97_otus.tree')")


    def test_profile_from_newick(self):
        import cProfile
        cProfile.run("TreeNode.from_newick(open('/home/evan/Downloads/97_otus.tree'))")
        cProfile.run("TreeNode.from_newick(open('/home/evan/Downloads/97_otus.tree'))")
        cProfile.run("TreeNode.from_newick(open('/home/evan/Downloads/97_otus.tree'))")


    def test_tree_node_to_newick(self):
        for tree, newicks in self.trees_newick_lists:
            newick = newicks[0]
            fh = StringIO()
            _tree_node_to_newick(tree, fh)

            self.assertEqual(newick, fh.getvalue())

            fh.close()

    def test_roundtrip_tree_node_to_newick_to_tree_node(self):
        for tree, newicks in self.trees_newick_lists:
            fh = StringIO()
            _tree_node_to_newick(tree, fh)
            fh.seek(0)
            new_tree = _newick_to_tree_node(fh)

            self.assertTrue(self._is_equal(tree, new_tree))

            fh.close()

    def test_roundtrip_newick_to_tree_node_to_newick(self):
        for tree, newicks in self.trees_newick_lists:
            newick = newicks[0]
            fh = StringIO(newick)
            tree = _newick_to_tree_node(fh)
            fh2 = StringIO()
            _tree_node_to_newick(tree, fh2)

            self.assertTrue(newick, fh2.getvalue())

            fh.close()
            fh2.close()


class TestNewickSniffer(unittest.TestCase):
    pass


if __name__ == '__main__':
    unittest.main()
