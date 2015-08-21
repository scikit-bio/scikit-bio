# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import main, TestCase

from skbio.diversity.beta.tests.test_unifrac_base import StatsTests
from skbio.diversity.beta import unweighted_unifrac_fast, weighted_unifrac_fast


class FastUnifracTests(StatsTests, TestCase):
    _method = {'unweighted_unifrac': unweighted_unifrac_fast,
               'weighted_unifrac': weighted_unifrac_fast}


if __name__ == '__main__':
    main()
