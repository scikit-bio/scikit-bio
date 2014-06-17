# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

from skbio import (
    global_pairwise_align_protein, local_pairwise_align_protein,
    global_pairwise_align_nucleotide, local_pairwise_align_nucleotide)


class PairwiseAlignmentTests(TestCase):

    def test_global_pairwise_align_protein(self):
        # expected results derived with assistance from
        # http://www.ebi.ac.uk/Tools/psa/emboss_needle/
        # In some cases, placement of non-gap characters surrounded by gap
        # characters are slighly different between skbio and the emboss server.
        # These differences arise from arbitrary implementation differences,
        # and always result in the same score.
        expected = ("HEAGAWGHEE", "---PAWHEAE", 1.0)
        actual = global_pairwise_align_protein("HEAGAWGHEE", "PAWHEAE",
                                               gap_open_penalty=10.,
                                               gap_extend_penalty=5.)
        self.assertEqual(actual, expected)

        expected = ("HEAGAWGHE-E", "---PAW-HEAE", 24.0)
        # EMBOSS diff: P---AW-HEAE
        actual = global_pairwise_align_protein("HEAGAWGHEE", "PAWHEAE",
                                               gap_open_penalty=5.,
                                               gap_extend_penalty=0.5)
        self.assertEqual(actual, expected)
