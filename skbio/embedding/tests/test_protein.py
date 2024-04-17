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
import numpy as np
from unittest import TestCase, main
from functools import partial
from pathlib import Path
from skbio.util import get_data_path
from skbio.embedding._protein import ProteinEmbedding
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

    def test_repr(self):
        emb, s = self.emb, self.seq
        p_emb = ProteinEmbedding(emb, s)
        self.assertEqual(p_emb.embedding.shape, (62, 1024))
        res_rstr = repr(p_emb)
        exp_rstr = (
            ("ProteinEmbedding\n"
             "--------------------------------------------------------------------\n"
             "Stats:\n"
             "    length: 62\n"
             "    embedding dimension: 1024\n"
             "    has gaps: False\n"
             "    has degenerates: False\n"
             "    has definites: True\n"
             "    has stops: False\n"
             "--------------------------------------------------------------------\n"
             "0  IGKEEIQQRL AQFVDHWKEL KQLAAARGQR LEESLEYQQF VANVEEEEAW INEKMTLVAS\n"
             "60 ED"))
        self.assertEqual(res_rstr, exp_rstr)


if __name__ == '__main__':
    main()
