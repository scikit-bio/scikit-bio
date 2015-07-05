# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest
import tempfile
import shutil
import io
import os.path

import skbio.io
from skbio.io.registry import open_file
from skbio.util import get_data_path


class TextSourceTests(object):
    pass


class WritableTextSourceTests(TextSourceTests):
    pass


class ReadableTextSourceTests(TextSourceTests):
    def test_open(self):
        pass

    def test_with_open(self):
        pass


class BinarySourceTests(object):
    pass



class ReadableBinarySourceTests(BinarySourceTests):
    def check_closed(self, file, expected):
        if hasattr(file, 'closed'):
            self.assertEqual(file.closed, expected)

    def check_open_state_contents(self, file, contents, is_binary, **kwargs):
        result = skbio.io.open(file, **kwargs)
        if is_binary:
            self.assertIsInstance(result, (io.BufferedReader,
                                           io.BufferedRandom))
        else:
            self.assertIsInstance(result, io.TextIOBase)
        self.assertTrue(result.readable())
        self.assertEqual(result.read(), contents)
        self.assertFalse(result.closed)

        result.close()
        self.assertTrue(result.closed)
        self.check_closed(file, True)

    def check_open_file_state_contents(self, file, contents, is_binary,
                                       **kwargs):
        with open_file(file, **kwargs) as f:
            if is_binary:
                self.assertIsInstance(f, (io.BufferedReader,
                                          io.BufferedRandom))
            else:
                self.assertIsInstance(f, io.TextIOBase)
            self.assertTrue(f.readable())
            self.assertEqual(f.read(), contents)
        self.assertEqual(f.closed, self.expected_close)
        self.check_closed(file, self.expected_close)

        f.close()
        self.assertTrue(f.closed)
        self.check_closed(file, True)

    def check_open_buffer_close_behaviour(self, file, **kwargs):
        if hasattr(file, 'close'):
            wrapped = skbio.io.open(file, **kwargs)
            file.close()
            self.assertTrue(wrapped.closed)

    def check_open_file_buffer_close_behaviour(self, file, **kwargs):
        if hasattr(file, 'close'):
            with open_file(file, **kwargs) as wrapped:
                file.close()
                self.assertTrue(wrapped.closed)

    def check_open_gc_behaviour(self, file, **kwargs):
        def mangle(file):
            skbio.io.open(file, **kwargs)

        f = skbio.io.open(file, encoding='binary')
        mangle(f)
        self.assertFalse(f.closed)

    def check_open_file_gc_behaviour(self, file, **kwargs):
        def mangle(file):
            with open_file(file, **kwargs):
                pass

        with open_file(file, encoding='binary') as f:
            mangle(f)
            self.assertFalse(f.closed)

    def test_open_gc_binary(self):
        self.check_open_gc_behaviour(self.read_file)

    def test_open_gc_encoding(self):
        self.check_open_gc_behaviour(self.encoded_file)

    def test_open_gc_compression(self):
        self.check_open_gc_behaviour(self.gzip_file)
        self.check_open_gc_behaviour(self.bz2_file)

    def test_open_gc_compression_encoding(self):
        self.check_open_gc_behaviour(self.gzip_encoded_file)
        self.check_open_gc_behaviour(self.bz2_encoded_file)

    def test_open_file_gc_binary(self):
        self.check_open_file_gc_behaviour(self.read_file)

    def test_open_file_gc_encoding(self):
        self.check_open_file_gc_behaviour(self.encoded_file)

    def test_open_file_gc_compression(self):
        self.check_open_file_gc_behaviour(self.gzip_file)
        self.check_open_file_gc_behaviour(self.bz2_file)

    def test_open_file_gc_compression_encoding(self):
        self.check_open_file_gc_behaviour(self.gzip_encoded_file)
        self.check_open_file_gc_behaviour(self.bz2_encoded_file)

    def test_open_underclose_binary(self):
        self.check_open_buffer_close_behaviour(self.read_file)

    def test_open_underclose_encoding(self):
        self.check_open_buffer_close_behaviour(self.encoded_file)

    def test_open_underclose_compression(self):
        self.check_open_buffer_close_behaviour(self.gzip_file)
        self.check_open_buffer_close_behaviour(self.bz2_file)

    def test_open_underclose_compression_encoding(self):
        self.check_open_buffer_close_behaviour(self.gzip_encoded_file)
        self.check_open_buffer_close_behaviour(self.bz2_encoded_file)

    def test_open_file_underclose_binary(self):
        self.check_open_file_buffer_close_behaviour(self.read_file)

    def test_open_file_underclose_encoding(self):
        self.check_open_file_buffer_close_behaviour(self.encoded_file)

    def test_open_file_underclose_compression(self):
        self.check_open_file_buffer_close_behaviour(self.gzip_file)
        self.check_open_file_buffer_close_behaviour(self.bz2_file)

    def test_open_file_underclose_compression_encoding(self):
        self.check_open_file_buffer_close_behaviour(self.gzip_encoded_file)
        self.check_open_file_buffer_close_behaviour(self.bz2_encoded_file)

    def test_open_binary(self):
        self.check_open_state_contents(self.read_file, self.binary_contents,
                                       True, mode='r', encoding='binary')

    def test_open_binary_compression_none(self):
        self.check_open_state_contents(self.read_file, self.binary_contents,
                                       True, mode='r', encoding='binary',
                                       compression=None)

    def test_open_encoding(self):
        self.check_open_state_contents(self.encoded_file,
                                       self.decoded_contents, False,
                                       mode='r', encoding=self.encoding)

    def test_open_auto_compression_binary(self):
        self.check_open_state_contents(self.gzip_file,
                                       self.binary_contents, True,
                                       mode='r', encoding='binary',
                                       compression='auto')

        self.check_open_state_contents(self.bz2_file,
                                       self.binary_contents, True,
                                       mode='r', encoding='binary',
                                       compression='auto')

    def test_open_gzip_compression_binary(self):
        self.check_open_state_contents(self.gzip_file,
                                       self.binary_contents, True,
                                       mode='r', encoding='binary',
                                       compression='gzip')

    def test_open_bz2_compression_binary(self):
        self.check_open_state_contents(self.bz2_file,
                                       self.binary_contents, True,
                                       mode='r', encoding='binary',
                                       compression='bz2')

    def test_open_default_compression_encoding(self):
        self.check_open_state_contents(self.gzip_encoded_file,
                                       self.decoded_contents, False,
                                       mode='r', encoding=self.encoding)

        self.check_open_state_contents(self.bz2_encoded_file,
                                       self.decoded_contents, False,
                                       mode='r', encoding=self.encoding)

    def test_open_file_binary(self):
        self.check_open_file_state_contents(self.read_file,
                                            self.binary_contents,
                                            True, mode='r', encoding='binary')

    def test_open_file_binary_compression_none(self):
        self.check_open_file_state_contents(self.read_file,
                                            self.binary_contents,
                                            True, mode='r', encoding='binary',
                                            compression=None)

    def test_open_file_encoding(self):
        self.check_open_file_state_contents(self.encoded_file,
                                            self.decoded_contents, False,
                                            mode='r', encoding=self.encoding)

    def test_open_file_auto_compression_binary(self):
        self.check_open_file_state_contents(self.gzip_file,
                                            self.binary_contents, True,
                                            mode='r', encoding='binary',
                                            compression='auto')

        self.check_open_file_state_contents(self.bz2_file,
                                            self.binary_contents, True,
                                            mode='r', encoding='binary',
                                            compression='auto')

    def test_open_file_gzip_compression_binary(self):
        self.check_open_file_state_contents(self.gzip_file,
                                            self.binary_contents, True,
                                            mode='r', encoding='binary',
                                            compression='gzip')

    def test_open_file_bz2_compression_binary(self):
        self.check_open_file_state_contents(self.bz2_file,
                                            self.binary_contents, True,
                                            mode='r', encoding='binary',
                                            compression='bz2')

    def test_open_file_default_compression_encoding(self):
        self.check_open_file_state_contents(self.gzip_encoded_file,
                                            self.decoded_contents, False,
                                            mode='r', encoding=self.encoding)

        self.check_open_file_state_contents(self.bz2_encoded_file,
                                            self.decoded_contents, False,
                                            mode='r', encoding=self.encoding)


class ReadableSourceTest(unittest.TestCase):
    def setUp(self):
        self.read_file = self.get_fileobj(get_data_path("example_file"))
        self.gzip_file = \
            self.get_fileobj(get_data_path("example_file.gz"))
        self.bz2_file = \
            self.get_fileobj(get_data_path("example_file.bz2"))
        self.encoded_file = self.get_fileobj(get_data_path("big5_file"))
        self.gzip_encoded_file = \
            self.get_fileobj(get_data_path("big5_file.gz"))
        self.bz2_encoded_file = \
            self.get_fileobj(get_data_path("big5_file.bz2"))

        self.binary_contents = (b"This is some content\n"
                                b"It occurs on more than one line\n")
        self.decoded_contents = u'\u4f60\u597d\n'  # Ni Hau
        self.compression = 'gzip'
        self.encoding = "big5"

    def tearDown(self):
        self.safe_close(self.read_file)
        self.safe_close(self.gzip_file)
        self.safe_close(self.bz2_file)
        self.safe_close(self.encoded_file)
        self.safe_close(self.gzip_encoded_file)
        self.safe_close(self.bz2_encoded_file)

    def safe_close(self, f):
        if hasattr(f, 'close'):
            f.close()



class WritableBinarySourceTests(BinarySourceTests):
    def check_closed(self, file, expected):
        if hasattr(file, 'closed'):
            self.assertEqual(file.closed, expected)

    def check_open_state_contents(self, file, contents, expected, is_binary, **kwargs):
        result = skbio.io.open(file, mode='w', **kwargs)
        if is_binary:
            self.assertIsInstance(result, (io.BufferedWriter,
                                           io.BufferedRandom))
        else:
            self.assertIsInstance(result, io.TextIOBase)
        self.assertTrue(result.writeable())

        result.write(contents)
        self.assertFalse(result.closed)

        self.assertEqual(self.get_contents(file), expected)

        result.close()
        self.assertTrue(result.closed)
        self.check_closed(file, True)

    def test_open_binary(self):
        with skbio.io.open(self.binary_file, mode='w', encoding='binary') as fh:
            fh.write(self.binary_contents)
        result = self.get_contents(self.binary_file)
        self.assertEqual(result, self.binary_contents)

class WritableSourceTest(unittest.TestCase):
    def setUp(self):
        self._dir = tempfile.mkdtemp()
        self.write_file = os.path.join(self._dir, "write_file")

        with io.open(get_data_path('example_file'), mode='rb') as f:
            self.binary_contents = f.read()
        self.binary_file = self._make_file('example_file')

        with io.open(get_data_path('big5_file'), mode='rb') as f:
            self.encoded_contents = f.read()
        self.encoded_file = self._make_file('big5_file')

        with io.open(get_data_path('example_file.gz'), mode='rb') as f:
            self.gzip_contents = f.read()
        self.gzip_file = self._make_file('example_file.gz')

        with io.open(get_data_path('example_file.bz2'), mode='rb') as f:
            self.bz2_contents = f.read()
        self.bz2_file = self._make_file('example_file.bz2')

        with io.open(get_data_path('big5_file.gz'), mode='rb') as f:
            self.gzip_encoded_contents = f.read()
        self.gzip_encoded_file = self._make_file('big5_file.gz')

        with io.open(get_data_path('big5_file.bz2'), mode='rb') as f:
            self.bz2_encoded_contents = f.read()
        self.bz2_encoded_file = self._make_file('big5_file.bz2')

        self.decoded_contents = self.encoded_contents.decode('big5')
        self.text_contents = self.binary_contents.decode('utf8')


    def tearDown(self):
        shutil.rmtree(self._dir)

    def _make_file(self, name):
        return self.get_fileobj(os.path.join(self._dir, name))


class TestReadFilepath(ReadableBinarySourceTests, ReadableSourceTest):
    expected_close = True

    def get_fileobj(self, path):
        return path

class TestWriteFilepath(WritableBinarySourceTests, WritableSourceTest):
    expected_close = True

    def get_fileobj(self, path):
        return path

    def get_contents(self, file):
        with io.open(file, mode='rb') as f:
            return f.read()


class TestReadBytesIO(ReadableBinarySourceTests, ReadableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        with io.open(path, mode='rb') as f:
            return io.BytesIO(f.read())


class TestReadBufferedReader(ReadableBinarySourceTests, ReadableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        return io.open(path, mode='rb')


class TestReadReadable(ReadableBinarySourceTests, ReadableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        class Readable(object):
            def __init__(self, file):
                self._file = file

            def __del__(self):
                self._file.close()

            def read(self, n):
                return self._file.read(n)

        return Readable(io.open(path, mode='rb'))

if __name__ == '__main__':
    unittest.main()
