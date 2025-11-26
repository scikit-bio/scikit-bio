# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile
import os

import numpy as np
import numpy.testing as npt
import h5py

from skbio.stats.distance import PairwiseMatrix, SymmetricMatrix, DistanceMatrix
from skbio.io import read as io_read
from skbio.io.format.binary_dm import (
    _h5py_mat_to_skbio_mat, _skbio_mat_to_h5py_mat, _get_header, _parse_ids,
    _verify_dimensions, _bytes_decoder, _passthrough_decoder, _set_header,
    _binary_dm_sniffer
)


class BinaryMatrixTests(unittest.TestCase):
    def setUp(self):
        _vlen_dtype = h5py.special_dtype(vlen=str)
        self.mat = np.array([[0, 0.1, 0.2],
                             [0.1, 0, 0.3],
                             [0.2, 0.3, 0]])
        self.ids = ['a', 'b', 'c']

        self.tempdir = tempfile.TemporaryDirectory()

        self.basic_fname = os.path.join(self.tempdir.name, 'basic')
        self.basic = h5py.File(self.basic_fname, 'a')
        ids = self.basic.create_dataset('order', shape=(3, ), dtype=_vlen_dtype)
        ids[:] = self.ids
        self.basic.create_dataset('matrix', data=self.mat)
        _set_header(self.basic)
        self.basic.close()

        self.badids_fname = os.path.join(self.tempdir.name, 'badids')
        self.badids = h5py.File(self.badids_fname, 'a')
        ids = self.badids.create_dataset('order', shape=(2, ), dtype=_vlen_dtype)
        ids[:] = ['a', 'b']
        self.badids.create_dataset('matrix', data=self.mat)
        _set_header(self.badids)
        self.badids.close()

        self.noheader_fname = os.path.join(self.tempdir.name, 'noheader')
        self.noheader = h5py.File(self.noheader_fname, 'a')
        ids = self.noheader.create_dataset('order', shape=(3, ), dtype=_vlen_dtype)
        ids[:] = self.ids
        self.noheader.create_dataset('matrix', data=self.mat)
        self.noheader.close()

        self.rw_fname = os.path.join(self.tempdir.name, 'rw')

    def tearDown(self):
        self.tempdir.cleanup()

    def test_binary_dm_sniffer(self):
        with open(self.basic_fname, 'rb') as fh:
            self.assertEqual((True, {}), _binary_dm_sniffer(fh))
        with open(self.badids_fname, 'rb') as fh:
            self.assertEqual((False, {}), _binary_dm_sniffer(fh))
        with open(self.noheader_fname, 'rb') as fh:
            self.assertEqual((False, {}), _binary_dm_sniffer(fh))

    def test_h5py_mat_to_skbio_mat(self):
        exp = DistanceMatrix(self.mat, self.ids)
        obs = _h5py_mat_to_skbio_mat(DistanceMatrix, h5py.File(self.basic_fname, 'r'))
        self.assertEqual(obs, exp)

    def test_io_h5py_mat_to_skbio_mat(self):
        # test DistanceMatrix
        exp = DistanceMatrix(self.mat, self.ids)
        obs = io_read(self.basic_fname,into=DistanceMatrix)
        self.assertEqual(obs, exp)

    def test_io_h5py_mat_to_skbio_mat_pairwise(self):
        # test PairwiseMatrix
        exp = PairwiseMatrix(self.mat, self.ids)
        obs = io_read(self.basic_fname,into=PairwiseMatrix)
        self.assertEqual(obs, exp)

    def test_io_h5py_mat_to_skbio_mat_symmetric(self):
        # test SymmetricMatrix
        exp = SymmetricMatrix(self.mat, self.ids)
        obs = io_read(self.basic_fname,into=SymmetricMatrix)
        self.assertEqual(obs, exp)

    def test_skbio_mat_to_h5py_mat(self):
        fh1 = h5py.File('f1', 'a', driver='core', backing_store=False)

        mat = DistanceMatrix(self.mat, self.ids)
        _skbio_mat_to_h5py_mat(mat, fh1)
        npt.assert_equal(np.asarray(fh1['order'][:], dtype=str), mat.ids)
        npt.assert_equal(fh1['matrix'], mat.data)

    def test_io_mat_to_h5py_mat(self):
        mat = DistanceMatrix(self.mat, self.ids)
        mat.write(self.rw_fname,format='binary_dm')
        exp = mat
        obs = io_read(self.rw_fname,into=DistanceMatrix)
        self.assertEqual(obs, exp)

    def test_io_mat_to_h5py_mat_symmetric(self):
        # SymmetricMatrix
        mat = SymmetricMatrix(self.mat, self.ids)
        mat.write(self.rw_fname,format='binary_dm')
        exp = mat
        obs = io_read(self.rw_fname,into=SymmetricMatrix)
        self.assertEqual(obs, exp)

    def test_io_mat_to_h5py_mat_pairwise(self):
        # PairwiseMatrix
        mat = PairwiseMatrix(self.mat, self.ids)
        mat.write(self.rw_fname,format='binary_dm')
        exp = mat
        obs = io_read(self.rw_fname,into=PairwiseMatrix)
        self.assertEqual(obs, exp)


    def test_get_header(self):
        self.assertEqual(_get_header(h5py.File(self.basic_fname, 'r')),
                         {'format': b'BDSM', 'version': b'2020.12'})
        self.assertEqual(_get_header(h5py.File(self.noheader_fname, 'r')), None)

    def test_parse_ids(self):
        tests = [(['a', 'b', 'c'], ['a', 'b', 'c']),
                 ([b'a', b'b', b'\xc3\xa9\xc3\xb8asd'],
                  ['a', 'b', 'éøasd'])]

        for test, exp in tests:
            self.assertEqual(_parse_ids(test), exp)

    def test_verify_dimensions(self):
        self.assertTrue(_verify_dimensions(h5py.File(self.basic_fname, 'r')))
        self.assertFalse(_verify_dimensions(h5py.File(self.badids_fname, 'r')))

    def test_bytes_decoder(self):
        test = [b'', b'a', b'\xc3\xa9\xc3\xb8asd']
        exp = ['', 'a', 'éøasd']
        self.assertEqual(_bytes_decoder(test), exp)

    def test_passthrough_decoder(self):
        tests = [('', ''), ('a', 'a'), ('éøasd', 'éøasd')]
        for test, expected in tests:
            self.assertEqual(_passthrough_decoder(test), expected)

    def test_set_header(self):
        def mock():
            obj = h5py.File('bung', 'a', driver='core', backing_store=False)
            return obj

        m = mock()
        _set_header(m)


if __name__ == '__main__':
    unittest.main()
