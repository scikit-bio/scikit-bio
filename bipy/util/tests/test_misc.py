#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from StringIO import StringIO

from bipy.util.misc import safe_md5
from bipy.util.unit_test import TestCase, main

class UtilTests(TestCase):
    """Test object for the miscellaneous utility functions"""
    def test_safe_md5(self):
        """Make sure we have the expected md5"""
        exp = 'ab07acbb1e496801937adfa772424bf7'

        fd = StringIO('foo bar baz')
        obs = safe_md5(fd)
        self.assertEqual(obs.hexdigest(),exp)

        fd.close()


if __name__ == '__main__':
    main()
