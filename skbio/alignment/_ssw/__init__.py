# -----------------------------------------------------------------------------
# Copyright (c) 2014--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from .ssw_wrapper import (
    StripedSmithWaterman, local_pairwise_align_ssw, AlignmentStructure)

from numpy.testing import Tester
test = Tester().test
