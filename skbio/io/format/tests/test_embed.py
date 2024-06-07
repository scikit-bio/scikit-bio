# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------
from unittest import TestCase, main
from pathlib import Path
import tempfile

import h5py
import numpy as np

import skbio
from skbio import Protein
from skbio.util import get_data_path
from skbio.embedding._protein import ProteinEmbedding
from skbio.embedding._protein import ProteinVector
from skbio.io.format.embed import (
    _embed_sniffer, _embed_to_generator,
    _embed_to_object, _generator_to_embed,
    _objects_to_embed,
    _embed_to_protein, _protein_to_embed,
    _protein_to_vector, _vector_to_protein
)


class EmbedTests(TestCase):
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
        self.writable_emb_path2 = str(tempdir / Path('test2.emb'))

        self.valid_embed_path = get_data_path('prot.emb')
        self.invalid_embed_path = str(tempdir / Path('invalid'))
        self.nonembed_hdf5_path = str(tempdir / Path('other.hdf5'))

        with open(self.invalid_embed_path, 'wb') as fp:
            fp.write(b'this is not a embed file')

        with h5py.File(self.nonembed_hdf5_path, 'w') as fp:
            fp['stuff'] = [1, 2, 3]

    def test_sniffer(self):
        self.assertEqual(_embed_sniffer(self.valid_embed_path), (True, {}))
        self.assertEqual(_embed_sniffer(self.invalid_embed_path), (False, {}))
        self.assertEqual(_embed_sniffer(self.nonembed_hdf5_path), (False, {}))

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
        objs1 = [ProteinEmbedding(emb, seq) for emb, seq in self.sequences]
        _generator_to_embed(objs1, self.writable_emb_path2)
        objs2 = _embed_to_generator(self.writable_emb_path2)
        for obj1, obj2 in zip(objs1, objs2):
            np.testing.assert_array_equal(obj1.embedding, obj2.embedding)
            self.assertEqual(str(obj1), str(obj2))

    def test_write_generator(self):
        sequences = [
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
        f = lambda x: ProteinEmbedding(*x)
        objs1 = (x for x in map(f, sequences))

        tempdir = Path(tempfile.mkdtemp())
        writable_emb_path = str(tempdir / Path('test.emb'))

        skbio.io.write(objs1, format='embed', into=writable_emb_path)
        objs2 = iter(skbio.io.read(writable_emb_path, format='embed',
                                    constructor=ProteinEmbedding))
        for obj1, obj2 in zip(objs1, objs2):
            np.testing.assert_array_equal(obj1.embedding, obj2.embedding)
            self.assertEqual(str(obj1), str(obj2))


class VectorTests(TestCase):
    def setUp(self):
        # single sequence
        rk = 10  # latent dimension of residues
        self.sequences = (
            [
                (
                    np.random.randn(rk),
                    Protein(('IGKEEIQQRLAQFVDHWKELKQLAAARGQRL'
                             'EESLEYQQFVANVEEEEAWINEKMTLVASED'),
                            metadata={"id": "seq1"})
                 ),
                (
                    np.random.randn(rk),
                    Protein(('QQNKELNFKLREKQNEIFELKKIAETLRSKL'
                             'EKYVDITKKLEDQNLNLQIKISDLEKKLSDA'),
                            metadata={"id": "seq2"})
                )
            ]
        )
        self.tempdir = tempfile.TemporaryDirectory()
        tempdir = Path(self.tempdir.name)
        self.writable_emb_path = str(tempdir / Path('test.emb'))
        self.writable_emb_path2 = str(tempdir / Path('test2.emb'))

        self.valid_embed_path = get_data_path('prot_vec.emb')
        self.invalid_embed_path = str(tempdir / Path('invalid'))
        self.nonembed_hdf5_path = str(tempdir / Path('other.hdf5'))

        with open(self.invalid_embed_path, 'wb') as fp:
            fp.write(b'this is not a embed file')

        with h5py.File(self.nonembed_hdf5_path, 'w') as fp:
            fp['stuff'] = [1, 2, 3]

    def test_sniffer(self):
        # make sure that the sniffer throws errors as expected
        self.assertEqual(_embed_sniffer(self.valid_embed_path), (True, {}))
        self.assertEqual(_embed_sniffer(self.invalid_embed_path), (False, {}))
        self.assertEqual(_embed_sniffer(self.nonembed_hdf5_path), (False, {}))
        emb, seq = self.sequences[0]
        obj = ProteinVector(emb, seq)
        _protein_to_vector(obj, str(Path(self.tempdir.name) / Path("prot_vec.emb")))

    def test_read_write_single(self):
        for emb, seq in self.sequences:
            fh = self.writable_emb_path
            obj = ProteinVector(emb, seq)
            _protein_to_vector(obj, fh)
            emb2 = _vector_to_protein(fh)
            np.testing.assert_array_equal(
                emb, emb2.embedding.ravel())
            self.assertEqual(str(seq), str(emb2))

    def test_read_write_generator(self):
        writable_emb_path2 = 'test2.emb'
        objs1 = [ProteinVector(emb, seq) for emb, seq in self.sequences]
        _generator_to_embed(objs1, self.writable_emb_path2)
        objs2 = _embed_to_generator(self.writable_emb_path2,
                                    constructor=ProteinVector)
        for obj1, obj2 in zip(objs1, objs2):
            np.testing.assert_array_equal(obj1.embedding, obj2.embedding)
            self.assertEqual(str(obj1), str(obj2))

    def test_write_generator(self):
        sequences = self.sequences
        f = lambda x: ProteinVector(*x)
        objs1 = (x for x in map(f, sequences))

        tempdir = Path(tempfile.mkdtemp())
        writable_emb_path = str(tempdir / Path('test.emb'))

        skbio.io.write(objs1, format='embed', into=writable_emb_path)
        objs2 = iter(skbio.io.read(writable_emb_path, format='embed',
                                    constructor=ProteinVector))
        for obj1, obj2 in zip(objs1, objs2):
            np.testing.assert_array_equal(obj1.embedding, obj2.embedding)
            self.assertEqual(str(obj1), str(obj2))


if __name__ == '__main__':
    main()
