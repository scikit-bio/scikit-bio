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
from skbio.embedding._embedding import SequenceEmbedding, Embedding
import numpy as np
import numpy.testing as npt


class EmbeddingTests(TestCase):
    def setUp(self):
        self.emb = np.random.randn(62, 10)
        self.seq = ("IGKEEIQQRLAQFVDHWKELKQLAAARGQRLEESLEYQ"
                    "QFVANVEEEEAWINEKMTLVASED")

    def test_id(self):
        emb, s = self.emb, self.seq
        p_emb = Embedding(emb, list(s))
        np.testing.assert_array_equal(p_emb.ids, np.array(list(s)))

    def test_embedding(self):
        emb, s = self.emb, self.seq
        p_emb = Embedding(emb, s)
        self.assertEqual(p_emb.embedding.shape, (62, 10))

    def test_assert_length(self):
        with self.assertRaises(ValueError):
            Embedding(self.emb, self.seq + "A")


class SequenceEmbeddingTests(TestCase):
    def setUp(self):
        self.emb = np.random.randn(62, 10)
        self.seq = ("IGKEEIQQRLAQFVDHWKELKQLAAARGQRLEESLEYQ"
                    "QFVANVEEEEAWINEKMTLVASED")

    def test_repr(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)
        self.assertTrue('SequenceEmbedding' in p_emb.__repr__())


    def test_str(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)

        self.assertEqual(p_emb.__str__(), s)
        self.assertEqual(p_emb.sequence, s)
        self.assertEqual(str(p_emb.ids.tobytes().decode('ascii')), s)

    def test_embedding(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)
        self.assertEqual(p_emb.embedding.shape, (62, 10))

    def test_assert_length(self):
        with self.assertRaises(ValueError):
            SequenceEmbedding(self.emb, self.seq + "A")



if __name__ == '__main__':
    main()
