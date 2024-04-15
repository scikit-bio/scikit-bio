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


class ProteinEmbeddingtests(TestCase):

    def test_repr(self):
        self.p_emb = ProteinEmbedding.read(get_data_path('prot.emb'))
        self.assertEqual(self.p_emb.embedding.shape, (62, 1024))
        res_rstr = repr(self.p_emb)
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
