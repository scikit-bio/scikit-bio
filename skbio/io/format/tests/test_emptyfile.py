# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import io

from skbio.io.format.emptyfile import _empty_file_sniffer


class TestEmptyFile(unittest.TestCase):
    def test_empty(self):
        res, kw = _empty_file_sniffer(io.StringIO())
        self.assertTrue(res)
        self.assertEqual({}, kw)

        res, kw = _empty_file_sniffer(io.StringIO("       \n   \t "))
        self.assertTrue(res)
        self.assertEqual({}, kw)

    def test_not_empty(self):
        res, kw = _empty_file_sniffer(io.StringIO("a"))
        self.assertFalse(res)
        self.assertEqual({}, kw)

        res, kw = _empty_file_sniffer(io.StringIO("                  \n \ta"))
        self.assertFalse(res)
        self.assertEqual({}, kw)


if __name__ == '__main__':
    unittest.main()
