#!/usr/bin/env python
"""Unit tests for parser support libraries dealing with records.
"""
from cStringIO import StringIO
from unittest import TestCase, main

from bipy.parse.mothur import parse_otu_list


class FunctionTests(TestCase):
    def test_parse_otu_list(self):
        observed = list(parse_otu_list(StringIO(mothur_output)))
        expected = [(0.0, [['cccccc'], ['bbbbbb'], ['aaaaaa']]),
                    (0.62, [['bbbbbb', 'cccccc'], ['aaaaaa']]),
                    (0.67000000000000004, [['aaaaaa', 'bbbbbb', 'cccccc']])]
        self.assertEqual(observed, expected)

mothur_output = """unique\t3\tcccccc\tbbbbbb\taaaaaa
0.62\t2\tbbbbbb,cccccc\taaaaaa
0.67\t1\taaaaaa,bbbbbb,cccccc
"""

if __name__ == '__main__':
    main()
