# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

import skbio.io
from skbio.util import get_data_path
from skbio.alignment import TabularMSA
from skbio.sequence import Sequence, DNA, Protein

from skbio.alignment._distance import align_dists


class AlignDistsTests(TestCase):
    def setUp(self):
        self.msa1 = TabularMSA([
            DNA("ATC-GTATCGT"),
            DNA("ATGCG--CCGC"),
            DNA("GT-CGTACG-T"),
            DNA("GT-NTTACAGT"),
        ], index=list("abcd"))

    def test_align_dists(self):
        obs = align_dists(self.msa1, "p_dist")


if __name__ == "__main__":
    main()
