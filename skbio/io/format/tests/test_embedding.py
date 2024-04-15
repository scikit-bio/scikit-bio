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

import numpy as np
import tempfile
from skbio import Sequence, DNA, RNA, Protein, TabularMSA
from skbio.embedding._protein import ProteinEmbedding
from skbio.io import FASTAFormatError, QUALFormatError
from skbio.io.format.embedding import (
    _embed_sniffer, _embed_to_generator,
    _embed_to_object, _generator_to_embed,
    _objects_to_embed,
    _embed_to_protein, _protein_to_embed
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
        self.tempdir = tempfile.TemporaryDirectory()
        tempdir = Path(self.tempdir.name)
        self.writable_emb_path = str(tempdir / Path('test.emb'))
        #self.writable_emb_path2 = str(tempdir / Path('test2.emb'))


    def test_read_write_single(self):
        for emb, seq in self.sequences:
            fh = self.writable_emb_path
            obj = ProteinEmbedding(emb, seq)
            _protein_to_embed(obj, fh)
            emb2 = _embed_to_protein(fh)
            np.testing.assert_array_equal(emb, emb2.embedding)
            self.assertEqual(str(seq), str(emb2))

    def test_read_write_generator(self):
        writable_emb_path2 = 'test2.emb'

        objs = [ProteinEmbedding(emb, seq) for emb, seq in self.sequences]
        # for obj in objs:
        #     _object_to_embed(obj, writable_emb_path2, 'a')

        _generator_to_embed(objs, writable_emb_path2)

        objs2 = _embed_to_generator(writable_emb_path2)
        for obj1, obj2 in zip(objs, objs2):
            np.testing.assert_array_equal(obj1.embedding, obj2.embedding)
            self.assertEqual(str(obj1), str(obj2))


if __name__ == '__main__':
    main()
