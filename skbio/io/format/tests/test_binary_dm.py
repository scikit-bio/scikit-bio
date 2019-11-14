# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import unittest

import numpy as np
import numpy.testing as npt
import h5py

from skbio import DissimilarityMatrix, DistanceMatrix
from skbio.io.format.binary_dm import (_h5py_mat_to_skbio_mat,
                                       _skbio_mat_to_h5py_mat, _get_header,
                                       _parse_ids, _verify_dimensions,
                                       _bytes_decoder, _passthrough_decoder,
                                       _set_header, _format_header)


class BinaryMatrixTests(unittest.TestCase):
    def setUp(self):
        self.mat = np.array([[0, 0.1, 0.2],
                             [0.1, 0, 0.3],
                             [0.2, 0.3, 0]])
        self.ids = ['a', 'b', 'c']

        self.basic = h5py.File('foo', driver='core', backing_store=False)
        self.basic.create_dataset('ids', data=self.ids)
        self.basic.create_dataset('matrix', data=self.mat)
        self.basic.attrs['symmetric'] = True
        _set_header(self.basic)

        self.badids = h5py.File('foo', driver='core', backing_store=False)
        self.badids.create_dataset('ids', data=['a', 'b'])
        self.badids.create_dataset('matrix', data=self.mat)
        self.basic.attrs['symmetric'] = True
        _set_header(self.badids)

        self.hoheader = h5py.File('foo', driver='core', backing_store=False)
        self.noheader.create_dataset('ids', data=self.ids)
        self.noheader.create_dataset('matrix', data=self.mat)
        self.basic.attrs['symmetric'] = True

        self.matdis = np.array([[0, 0.1, 0.2],
                                [0.2, 0, 0.3],
                                [0.3, 0.4, 0]])
        self.dis = h5py.File('foo', driver='core', backing_store=False)
        self.dis.create_dataset('ids', data=self.ids)
        self.dis.create_dataset('matrix', data=self.matdis)
        self.basic.attrs['symmetric'] = False
        _set_header(self.dis)

        self.baddis = h5py.File('foo', driver='core', backing_store=False)
        self.baddis.create_dataset('ids', data=self.ids)
        self.baddis.create_dataset('matrix', data=self.matdis)
        self.basic.attrs['symmetric'] = True
        _set_header(self.baddis)

    def test_h5py_mat_to_skbio_mat(self):
        exp = DistanceMatrix(self.mat, self.ids)
        obs = _h5py_mat_to_skbio_mat(self.basic)
        self.assertEqual(obs, exp)

        exp = DissimilarityMatrix(self.dismat, self.ids)
        obs = _h5py_mat_to_skbio_mat(self.dis)
        self.assertEqual(obs, exp)

    def test_skbio_mat_to_h5py_mat(self):
        mat = DistanceMatrix(self.mat, self.ids)
        obs = _skbio_mat_to_h5py_mat(mat)
        self.assertEqual(obs['ids'], mat.ids)
        npt.assert_equal(obs['matrix'], mat.data)

        mat = DissimilarityMatrix(self.dismat, self.ids)
        obs = _skbio_mat_to_h5py_mat(mat)
        self.assertEqual(obs['ids'], mat.ids)
        npt.assert_equal(obs['matrix'], mat.data)

    def test_soft_validate_symmetric(self):
        self.assertTrue(self.basic)
        self.assertTrue(self.dis)
        self.assertFalse(self.baddis)

    def test_get_header(self):
        self.assertEqual(_get_header(self.basic), _format_header)
        self.assertEqual(_get_header(self.noheader), None)

    def test_parse_ids(self):
        def mock(ids):
            obj = h5py.File('foo', driver='core', backing_store=False)
            obj.create_dataset('ids', data=ids)
            return obj

        tests = [(mock(['a', 'b', 'c']), ['a', 'b', 'c']),
                 (mock([b'', b'a', b'\xc3\xa9\xc3\xb8asd']),
                  'a', 'b', 'éøasd')]

        for test, exp in tests:
            self.assertEqual(_parse_ids(test), exp)

    def test_verify_dimensions(self):
        self.assertTrue(_verify_dimensions(self.basic))
        self.assertFalse(_verify_dimensions(self.badids))

    def test_bytes_decoder(self):
        tests = [(b'', ''), (b'a', 'a'), (b'\xc3\xa9\xc3\xb8asd', 'éøasd')]
        for test, expected in tests:
            self.assertEqual(_bytes_decoder(test), expected)

    def test_passthrough_decoder(self):
        tests = [('', ''), ('a', 'a'), ('éøasd', 'éøasd')]
        for test, expected in tests:
            self.assertEqual(_passthrough_decoder(test), expected)

    def test_set_header(self):
        class mock:
            attrs = {}
        m = mock()
        _set_header(m)
        self.assertEqual(m, _format_header)


if __name__ == '__main__':
    unittest.main()
