# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from pathlib import Path

from skbio.io._iosources import IOSource, Compressor


class TestIOSource(unittest.TestCase):

    def setUp(self):
        self.file = 'somepath'
        self.file_path = Path('somepath')
        self.options = {'a': 1, 'b': 2}

        self.source = IOSource(self.file, self.options)
        self.source_path = IOSource(self.file_path, self.options)

    def test_attributes(self):
        self.assertEqual(self.source.file, self.file)
        self.assertEqual(self.source.options, self.options)
        self.assertEqual(self.source_path.file, self.file_path)
        self.assertEqual(self.source_path.options, self.options)

    def test_can_read(self):
        self.assertEqual(self.source.can_read(), False)
        self.assertEqual(self.source_path.can_read(), False)

    def test_can_write(self):
        self.assertEqual(self.source.can_write(), False)
        self.assertEqual(self.source_path.can_write(), False)

    def test_get_reader(self):
        with self.assertRaises(NotImplementedError):
            self.source.get_reader()
            self.source_path.get_reader()

    def test_get_writer(self):
        with self.assertRaises(NotImplementedError):
            self.source.get_writer()
            self.source_path.get_writer()


class TestCompressor(TestIOSource):
    def setUp(self):
        super(TestCompressor, self).setUp()
        self.compressor = Compressor(self.file, self.options)

    def test_can_write(self):
        self.assertEqual(self.compressor.can_write(), True)


if __name__ == "__main__":
    unittest.main()
