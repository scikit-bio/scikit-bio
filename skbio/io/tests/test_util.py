# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from six import StringIO, BytesIO
from requests import HTTPError

import unittest
import tempfile
import gzip
import shutil
from os import path

from skbio.io.util import (open_file, open_files,
                           _is_string_or_bytes, get_filemode)


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
        with open_file(filepath) as fh:
            pass
        self.assertTrue(fh.closed)

    def test_file_closed_harder(self):
        """File gets closed in decorator, even if exceptions happen."""
        f = tempfile.NamedTemporaryFile('r')
        filepath = f.name
        try:
            with open_file(filepath) as fh:
                raise TypeError
        except TypeError:
            self.assertTrue(fh.closed)
        else:
            # If we're here, no exceptions have been raised inside the
            # try clause, so the context manager swallowed them. No
            # good.
            raise Exception("`open_file` didn't propagate exceptions")

    def test_filehandle(self):
        """Filehandles slip through untouched"""
        with tempfile.TemporaryFile('r') as fh:
            with open_file(fh) as ffh:
                self.assertTrue(fh is ffh)
            # And it doesn't close the file-handle
            self.assertFalse(fh.closed)

    def test_StringIO(self):
        """StringIO (useful e.g. for testing) slips through."""
        f = StringIO("File contents")
        with open_file(f) as fh:
            self.assertTrue(fh is f)

    def test_BytesIO(self):
        """BytesIO (useful e.g. for testing) slips through."""
        f = BytesIO(b"File contents")
        with open_file(f) as fh:
            self.assertTrue(fh is f)


class TestFilePathsOpening(unittest.TestCase):

    def test_files_closed(self):
        """File gets closed in decorator"""
        f = tempfile.NamedTemporaryFile('r')
        f2 = tempfile.NamedTemporaryFile('r')
        filepath = f.name
        filepath2 = f2.name
        with open_files([filepath, filepath2]) as fhs:
            pass
        for fh in fhs:
            self.assertTrue(fh.closed)

    def test_files_closed_harder(self):
        """File gets closed in decorator, even if exceptions happen."""
        f = tempfile.NamedTemporaryFile('r')
        f2 = tempfile.NamedTemporaryFile('r')
        filepath = f.name
        filepath2 = f2.name
        try:
            with open_files([filepath, filepath2]) as fhs:
                raise TypeError
        except TypeError:
            for fh in fhs:
                self.assertTrue(fh.closed)
        else:
            # If we're here, no exceptions have been raised inside the
            # try clause, so the context manager swallowed them. No
            # good.
            raise Exception("`open_file` didn't propagate exceptions")

    def test_filehandle(self):
        """Filehandles slip through untouched"""
        with tempfile.TemporaryFile('r') as fh:
            with tempfile.TemporaryFile('r') as fh2:
                with open_file([fh, fh2]) as fhs:
                    self.assertTrue(fh is fhs[0])
                    self.assertTrue(fh2 is fhs[1])
                # And it doesn't close the file-handle
                for fh in fhs:
                    self.assertFalse(fh.closed)

    def test_StringIO(self):
        """StringIO (useful e.g. for testing) slips through."""
        f = StringIO("File contents")
        with open_files([f]) as fhs:
            self.assertTrue(fhs[0] is f)

    def test_BytesIO(self):
        """BytesIO (useful e.g. for testing) slips through."""
        f = BytesIO(b"File contents")
        with open_files([f]) as fhs:
            self.assertTrue(fhs[0] is f)

    def test_remote_failing_fna(self):
        with self.assertRaises(HTTPError) as e:
            with open_files(['http://google.com/foo-seqs.fna']) as fhs:
                for f in fhs:
                    f.read()
        self.assertEquals(str(e.exception), '404 Client Error: Not Found')

    def test_remote_fna(self):
        url = ('http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmax=1'
               '00&retmode=text&tool=skbio&db=nucleotide&id=459567&rettype=fas'
               'ta&retstart=0&email=foo@bar.com')
        with open_files([url]) as fhs:
            for f in fhs:
                self.assertEqual(f.read(), FASTA)

    def test_remote_fna_kwargs(self):
        url = ('http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmax=1'
               '00&retmode=text&tool=skbio&db=nucleotide&id=459567&rettype=fas'
               'ta&retstart=0&email=foo@bar.com')
        with open_files([url], stream=True) as fhs:
            for f in fhs:
                self.assertEqual(f.read(), FASTA)

    def test_gzip(self):
        try:
            tmp_dir = tempfile.mkdtemp()
            tmp_path = path.join(tmp_dir, 'test.fa.gz')

            content = FASTA.decode('utf-8')
            with gzip.open(tmp_path, 'wt') as fh:
                fh.write(content)

            # Vanilla.
            with open_files([fh.name]) as fhs:
                self.assertEqual(fhs[0].read(), content)

            # Mode in kwargs.
            with open_files([fh.name], mode='U') as fhs:
                self.assertEqual(fhs[0].read(), content)

            # Mode in args.
            with open_files([fh.name], 'rU') as fhs:
                self.assertEqual(fhs[0].read(), content)
        finally:
            shutil.rmtree(tmp_dir)


FASTA = (b'>gi|459567|dbj|D28543.1|HPCNS5PC Hepatitis C virus gene for NS5 pr'
         b'otein, partial cds, isolate: B4/92\nGAGCACGACATCTACCAATGTTGCCAACTG'
         b'AACCCAGAGGCCAAGAAAGCCATAACATCCTTGACAGAGA\nGGCTTTACCTTGGTGGTCCCATGT'
         b'TTAACTCGCGAGGTCAGCTCTGCGGGACACGCAGATGCCGGGCGAG\nCGGGGTTCTTCCAACCAG'
         b'CATGGGCAATACCCTCACATGTTACCTGAAAGCACAGGCAGCTTGCCGTGCA\nGCAGGCCTCACC'
         b'AATTCTGACATGTTGGTTTGCGGAGATGATTTGGTAGTCATCACTGAGAGTGCCGGAG\nTC\n\n')


class TestGetFilemode(unittest.TestCase):

    def test_file(self):
        with tempfile.NamedTemporaryFile('r') as fh:
            self.assertEqual(get_filemode(fh), 'r')

    def test_gzip_file(self):
        try:
            tmp_dir = tempfile.mkdtemp()
            tmp_path = path.join(tmp_dir, 'test.fa.gz')

            with gzip.open(tmp_path, 'wb') as fh:
                fh.write(FASTA)

            with gzip.open(tmp_path, 'rb') as fh:
                self.assertEqual(get_filemode(fh), 'rb')

            with gzip.open(tmp_path, 'rt') as fh:
                self.assertEqual(get_filemode(fh), 'rt')
        finally:
            shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    unittest.main()
