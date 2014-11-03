# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio._base import SkbioObject


class TestSkbioObject(unittest.TestCase):
    def test_no_instantiation(self):
        class Foo(SkbioObject):
            pass

        with self.assertRaises(TypeError):
            Foo()


if __name__ == '__main__':
    unittest.main()
