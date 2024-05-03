# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import copy
import io
import string
from unittest import TestCase, main
from functools import partial
from pathlib import Path
from skbio.util import get_data_path
from skbio.embedding._protein import ProteinEmbedding
import numpy as np
import numpy.testing as npt

class ProteinEmbeddingtests(TestCase):

    def setUp(self):
        self.emb = np.load(get_data_path('embed1.txt.npy'))
        self.seq = ("IGKEEIQQRLAQFVDHWKELKQLAAARGQRLEESLEYQ"
                    "QFVANVEEEEAWINEKMTLVASED")

    def test_clipping(self):
        emb, s = self.emb, self.seq
        nemb = np.zeros((emb.shape[0] + 2, emb.shape[1]))
        nemb[1:-1] = emb
        p2_emb = ProteinEmbedding(nemb, s, clip_head=True, clip_tail=True)
        npt.assert_array_equal(p2_emb.embedding, emb)
        self.assertEqual(p2_emb.sequence, s)

    def test_str(self):
        emb, s = self.emb, self.seq
        p_emb = ProteinEmbedding(emb, s)
        self.assertEqual(str(p_emb), s)
        self.assertEqual(p_emb.sequence, s)

    def test_str_spaces(self):
        seq = ("I G K E E I Q Q R L A Q F V D H W K E L K Q L A "
               "A A R G Q R L E E S L E Y Q Q F V A N V E E E E "
               "A W I N E K M T L V A S E D")
        p_emb = ProteinEmbedding(self.emb, seq)
        self.assertEqual(str(p_emb), self.seq)
        self.assertEqual(p_emb.sequence, self.seq)

    def test_assert_length(self):
        with self.assertRaises(ValueError):
            ProteinEmbedding(self.emb, self.seq + "A")


if __name__ == '__main__':
    main()
