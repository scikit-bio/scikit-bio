# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from future.utils.six import StringIO, BytesIO

import unittest
import tempfile

from skbio.util.io import open_filepath_or, _is_string_or_bytes


class TestFilePathOpening(unittest.TestCase):
    def test_is_string_or_bytes(self):
        self.assertTrue(_is_string_or_bytes('foo'))
        self.assertTrue(_is_string_or_bytes(u'foo'))
        self.assertTrue(_is_string_or_bytes(b'foo'))
        self.assertFalse(_is_string_or_bytes(StringIO('bar')))
        self.assertFalse(_is_string_or_bytes([1]))

    def test_file_closed(self):
        """File gets closed in decorator"""
        f = tempfile.NamedTemporaryFile('r')
        filepath = f.name
        with open_filepath_or(filepath) as fh:
            pass
        self.assertTrue(fh.closed)

    def test_file_closed_harder(self):
        """File gets closed in decorator, even if exceptions happen."""
        f = tempfile.NamedTemporaryFile('r')
        filepath = f.name
        try:
            with open_filepath_or(filepath) as fh:
                raise TypeError
        except TypeError:
            self.assertTrue(fh.closed)
        else:
            # If we're here, no exceptions have been raised inside the
            # try clause, so the context manager swallowed them. No
            # good.
            raise Exception("`open_filepath_or` didn't propagate exceptions")

    def test_filehandle(self):
        """Filehandles slip through untouched"""
        with tempfile.TemporaryFile('r') as fh:
            with open_filepath_or(fh) as ffh:
                self.assertTrue(fh is ffh)
            # And it doesn't close the file-handle
            self.assertFalse(fh.closed)

    def test_StringIO(self):
        """StringIO (useful e.g. for testing) slips through."""
        f = StringIO("File contents")
        with open_filepath_or(f) as fh:
            self.assertTrue(fh is f)

    def test_BytesIO(self):
        """BytesIO (useful e.g. for testing) slips through."""
        f = BytesIO(b"File contents")
        with open_filepath_or(f) as fh:
            self.assertTrue(fh is f)

if __name__ == '__main__':
    unittest.main()
