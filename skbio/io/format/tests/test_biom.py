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

from skbio.table import Table, example_table
from skbio.io.format.biom import (
    _biom_to_table, _table_to_biom, _biom_sniffer)


class BIOMFormatTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        tempdir = Path(self.tempdir.name)
        self.valid_biom_path = str(tempdir / Path('valid.biom'))
        self.invalid_biom_path = str(tempdir / Path('invalid'))
        self.writable_biom_path = str(tempdir / Path('write.biom'))
        self.nonbiom_hdf5_path = str(tempdir / Path('other.hdf5'))
        self.difbiomver_path = str(tempdir / Path('otherver.biom'))

        self.table = example_table.copy()
        with h5py.File(self.valid_biom_path, 'w') as fp:
            self.table.to_hdf5(fp, 'unit-test')

        with open(self.invalid_biom_path, 'wb') as fp:
            fp.write(b'this is not a biom file')

        with h5py.File(self.nonbiom_hdf5_path, 'w') as fp:
            fp['stuff'] = [1, 2, 3]

        self.table = example_table.copy()
        with h5py.File(self.difbiomver_path, 'w') as fp:
            self.table.to_hdf5(fp, 'unit-test')
            fp.attrs['format-version'] = [3, 0]

    def tearDown(self):
        self.tempdir.cleanup()

    def test_sniffer(self):
        self.assertEqual(_biom_sniffer(self.valid_biom_path), (True, {}))
        self.assertEqual(_biom_sniffer(self.invalid_biom_path), (False, {}))
        self.assertEqual(_biom_sniffer(self.nonbiom_hdf5_path), (False, {}))
        self.assertEqual(_biom_sniffer(self.difbiomver_path), (False, {}))

    def test_biom_to_table(self):
        tab = _biom_to_table(self.valid_biom_path)
        self.assertEqual(tab, self.table)

    def test_table_to_biom(self):
        _table_to_biom(self.table, self.writable_biom_path)
        roundtrip = _biom_to_table(self.writable_biom_path)
        self.assertEqual(roundtrip, self.table)


if __name__ == '__main__':
    unittest.main()
