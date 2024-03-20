# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import unittest
import tempfile
from pathlib import Path

import h5py
import biom

from skbio import Table
from skbio.io import BIOMFormatError
from skbio.io.format.biom import (
    _biom_to_feature_table, _feature_table_to_biom, _biom_sniffer)


class BIOMFormatTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        tempdir = Path(self.tempdir.name)
        self.valid_biom_path = str(tempdir / Path('valid.biom'))
        self.invalid_biom_path = str(tempdir / Path('invalid'))
        self.writable_biom_path = str(tempdir / Path('write.biom'))

        self.table = biom.example_table.copy()
        with h5py.File(self.valid_biom_path, 'w') as fp:
            self.table.to_hdf5(fp, 'unit-test')

        with open(self.invalid_biom_path, 'wb') as fp:
            fp.write(b'this is not a biom file')

    def tearDown(self):
        self.tempdir.cleanup()

    def test_sniffer(self):
        self.assertEqual(_biom_sniffer(self.valid_biom_path), (True, {}))
        self.assertEqual(_biom_sniffer(self.invalid_biom_path), (False, {}))

    def test_biom_to_feature_table(self):
        tab = _biom_to_feature_table(self.valid_biom_path)
        self.assertEqual(tab, self.table)

    def test_feature_table_to_biom(self):
        _feature_table_to_biom(self.table, self.writable_biom_path)
        roundtrip = _biom_to_feature_table(self.writable_biom_path)
        self.assertEqual(roundtrip, self.table)


if __name__ == '__main__':
    unittest.main()
