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

import numpy as np

from skbio import Sequence, DNA, RNA, Protein, TabularMSA
from skbio.embedding._protein import ProteinEmbedding
from skbio.io import FASTAFormatError, QUALFormatError
from skbio.io.format.embedding import (
    _embed_sniffer, _embed_to_generator,
    _embed_to_object, _generator_to_embed,
    _object_to_embed
)
from skbio.util import get_data_path


class EmbeddingTests(TestCase):
    def setUp(self):
        # single sequence
        rk = 5  # latent dimension of residues
        self.sequences = (
            [
                (
                    np.load(get_data_path('embed1.txt.npy')),
                    Protein(('IGKEEIQQRLAQFVDHWKELKQLAAARGQRL'
                             'EESLEYQQFVANVEEEEAWINEKMTLVASED'),
                            metadata={"id": "seq1"})
                 ),
                (
                    np.load(get_data_path('embed2.txt.npy')),
                    Protein(('QQNKELNFKLREKQNEIFELKKIAETLRSKL'
                             'EKYVDITKKLEDQNLNLQIKISDLEKKLSDA'),
                            metadata={"id": "seq2"})
                )
            ]
        )

    def test_writer_single(self):
        for emb, seq in self.sequences:
            fh = io.BytesIO()
            obj = ProteinEmbedding(emb, seq)
            print(dir(obj))
            print(repr(obj))
            _object_to_embed(obj, fh)
            fh.seek(0)
            emb2, id_ = _embed_to_object(fh)
            np.testing.assert_array_equal(emb, emb2)
            self.assertEqual(str(seq), id_)


if __name__ == '__main__':
    main()
