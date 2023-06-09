# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from skbio._base import SkbioObject, ElasticLines


class TestSkbioObject(unittest.TestCase):
    def test_no_instantiation(self):
        class Foo(SkbioObject):
            pass

        with self.assertRaises(TypeError):
            Foo()


class TestElasticLines(unittest.TestCase):
    def setUp(self):
        self.el = ElasticLines()

    def test_empty(self):
        self.assertEqual(self.el.to_str(), '')

    def test_add_line(self):
        self.el.add_line('foo')
        self.assertEqual(self.el.to_str(), 'foo')

    def test_add_lines(self):
        self.el = ElasticLines()
        self.el.add_lines(['alice', 'bob', 'carol'])
        self.assertEqual(self.el.to_str(), 'alice\nbob\ncarol')

    def test_add_separator(self):
        self.el.add_separator()
        self.assertEqual(self.el.to_str(), '')

        self.el.add_line('foo')
        self.assertEqual(self.el.to_str(), '---\nfoo')

        self.el.add_separator()
        self.el.add_lines(['bar', 'bazzzz'])
        self.el.add_separator()

        self.assertEqual(self.el.to_str(),
                         '------\nfoo\n------\nbar\nbazzzz\n------')


if __name__ == '__main__':
    unittest.main()
