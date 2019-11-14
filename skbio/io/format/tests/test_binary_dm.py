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

from skbio.stats.distance import DissimilarityMatrix, DistanceMatrix
from skbio.io.format.binary_dm import (_h5py_mat_to_skbio_mat,
                                       _skbio_mat_to_h5py_mat, _get_header,
                                       _parse_ids, _verify_dimensions,
                                       _bytes_decoder, _passthrough_decoder,
                                       _set_header, _format_header,
                                       _soft_validate_symmetric, _vlen_dtype,
                                       _binary_dm_sniffer,
                                       _binary_dm_to_dissimiarlity,
                                       _binary_dm_to_distance,
                                       _distance_to_binary_dm,
                                       _dissimilarity_to_binary_dm)


class BinaryMatrixTests(unittest.TestCase):
    def setUp(self):
        self.mat = np.array([[0, 0.1, 0.2],
                             [0.1, 0, 0.3],
                             [0.2, 0.3, 0]])
        self.ids = ['a', 'b', 'c']

        self.basic = h5py.File('foo', driver='core', backing_store=False)
        ids = self.basic.create_dataset('ids', shape=(3, ), dtype=_vlen_dtype)
        ids[:] = self.ids
        self.basic.create_dataset('matrix', data=self.mat)
        self.basic.attrs['symmetric'] = True
        _set_header(self.basic)
        self.fh = [self.basic]

        self.badids = h5py.File('bar', driver='core', backing_store=False)
        ids = self.badids.create_dataset('ids', shape=(2, ), dtype=_vlen_dtype)
        ids[:] = ['a', 'b']
        self.badids.create_dataset('matrix', data=self.mat)
        self.basic.attrs['symmetric'] = True
        _set_header(self.badids)
        self.fh.append(self.badids)

        self.noheader = h5py.File('baz', driver='core', backing_store=False)
        ids = self.noheader.create_dataset('ids', shape=(3, ),
                                           dtype=_vlen_dtype)
        ids[:] = self.ids
        self.noheader.create_dataset('matrix', data=self.mat)
        self.noheader.attrs['symmetric'] = True
        self.fh.append(self.noheader)

        self.matdis = np.array([[0, 0.1, 0.2],
                                [0.2, 0, 0.3],
                                [0.3, 0.4, 0]])
        self.dis = h5py.File('bing', driver='core', backing_store=False)
        ids = self.dis.create_dataset('ids', shape=(3, ), dtype=_vlen_dtype)
        ids[:] = ['a', 'b', 'c']
        self.dis.create_dataset('matrix', data=self.matdis)
        self.dis.attrs['symmetric'] = False
        _set_header(self.dis)
        self.fh.append(self.dis)

        self.baddis = h5py.File('bang', driver='core', backing_store=False)
        self.baddis.create_dataset('ids', shape=(3, ), dtype=_vlen_dtype)
        ids[:] = self.ids
        self.baddis.create_dataset('matrix', data=self.matdis)
        self.baddis.attrs['symmetric'] = True
        _set_header(self.baddis)
        self.fh.append(self.baddis)

    def tearDown(self):
        for f in self.fh:
            f.close()

    def test_binary_dm_sniffer(fh):
        pass

    def test_binary_dm_to_dissimilarity(fh):
        pass

    def test_binary_dm_to_distance(fh):
        pass

    def test_dissimilarity_to_binary_dm(obj, fh):
        pass

    def test__distance_to_binary_dm(obj, fh):
        pass

    def test_h5py_mat_to_skbio_mat(self):
        exp = DistanceMatrix(self.mat, self.ids)
        obs = _h5py_mat_to_skbio_mat(DistanceMatrix, self.basic)
        self.assertEqual(obs, exp)

        exp = DissimilarityMatrix(self.matdis, self.ids)
        obs = _h5py_mat_to_skbio_mat(DissimilarityMatrix, self.dis)
        self.assertEqual(obs, exp)

    def test_skbio_mat_to_h5py_mat(self):
        fh1 = h5py.File('f1', driver='core', backing_store=False)
        fh2 = h5py.File('f2', driver='core', backing_store=False)
        self.fh.append(fh1)
        self.fh.append(fh2)

        mat = DistanceMatrix(self.mat, self.ids)
        _skbio_mat_to_h5py_mat(mat, fh1)
        npt.assert_equal(fh1['ids'][:], mat.ids)
        npt.assert_equal(fh1['matrix'], mat.data)

        mat = DissimilarityMatrix(self.matdis, self.ids)
        _skbio_mat_to_h5py_mat(mat, fh2)
        npt.assert_equal(fh2['ids'][:], mat.ids)
        npt.assert_equal(fh2['matrix'], mat.data)

    def test_soft_validate_symmetric(self):
        self.assertTrue(_soft_validate_symmetric(self.basic))
        self.assertFalse(_soft_validate_symmetric(self.dis))

    def test_get_header(self):
        self.assertEqual(_get_header(self.basic), _format_header)
        self.assertEqual(_get_header(self.noheader), None)

    def test_parse_ids(self):
        tests = [(['a', 'b', 'c'], ['a', 'b', 'c']),
                 ([b'a', b'b', b'\xc3\xa9\xc3\xb8asd'],
                  ['a', 'b', 'éøasd'])]

        for test, exp in tests:
            self.assertEqual(_parse_ids(test), exp)

    def test_verify_dimensions(self):
        self.assertTrue(_verify_dimensions(self.basic))
        self.assertFalse(_verify_dimensions(self.badids))

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
            obj = h5py.File('bung', driver='core', backing_store=False)
            return obj

        m = mock()
        self.fh.append(m)
        _set_header(m)
        obs = {k: m.attrs[k] for k in _format_header}
        self.assertEqual(obs, _format_header)


if __name__ == '__main__':
    unittest.main()
