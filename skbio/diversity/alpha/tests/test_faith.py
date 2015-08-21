# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import main, TestCase

from skbio.diversity.alpha.tests.test_faith_base import FaithTests
from skbio.diversity.alpha import faith_pd_fast


class FaithTests(FaithTests, TestCase):
    _method = {'faith_pd': faith_pd_fast}


if __name__ == '__main__':
    main()
