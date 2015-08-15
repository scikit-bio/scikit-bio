
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import main

from skbio.diversity.beta.tests.test_unifrac_base import StatsTests
from skbio.diversity.beta._fast_unifrac import (unweighted_unifrac,
                                                weighted_unifrac)


StatsTests._method = {'unweighted_unifrac': unweighted_unifrac,
                      'weighted_unifrac': weighted_unifrac}


if __name__ == '__main__':
    main()
