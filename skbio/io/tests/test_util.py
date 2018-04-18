# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile
import shutil
import io
import os.path

try:
    import httpretty
    has_httpretty = True
except ImportError:
    has_httpretty = False

import skbio.io
from skbio.io.registry import open_file
from skbio.util import get_data_path


class TestOpen(unittest.TestCase):
    def test_open_invalid_mode(self):
        with self.assertRaises(ValueError):
            skbio.io.open([], mode='a')

    def test_open_invalid_source(self):
        with self.assertRaises(skbio.io.IOSourceError):
            skbio.io.open(42)

    def test_open_invalid_source_compression(self):
        with self.assertRaises(ValueError):
            skbio.io.open(['foo'], compression='gzip')

    def test_open_invalid_source_encoding(self):
        with self.assertRaises(ValueError):
            skbio.io.open(['foo'], encoding='binary')

        with self.assertRaises(ValueError):
            skbio.io.open(['foo'], encoding='binary', newline='\r')

    def test_open_invalid_compression(self):
        with self.assertRaises(ValueError):
            skbio.io.open(io.BytesIO(), compression='foo')


class ReadableBinarySourceTests:
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
            result = skbio.io.open(file, **kwargs)
            self.assertIsInstance(result, io.TextIOBase)

        f = skbio.io.open(file, encoding='binary')
        mangle(f)
        self.assertFalse(f.closed)
        f.close()

    def check_open_file_gc_behaviour(self, file, **kwargs):
        def mangle(file):
            with open_file(file, **kwargs) as result:
                self.assertIsInstance(result, io.TextIOBase)

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
        self.decoded_contents = '\u4f60\u597d\n'  # Ni Hau
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


class WritableBinarySourceTests:
    def check_closed(self, file, expected):
        if hasattr(file, 'closed'):
            self.assertEqual(file.closed, expected)

    def check_open_state_contents(self, file, contents, is_binary,
                                  **kwargs):
        result = skbio.io.open(file, mode='w', **kwargs)
        if is_binary:
            self.assertIsInstance(result, (io.BufferedWriter,
                                           io.BufferedRandom))
        else:
            self.assertIsInstance(result, io.TextIOBase)
        self.assertTrue(result.writable())

        result.write(contents)
        self.assertFalse(result.closed)

        if self.expected_close:
            result.close()
            self.assertTrue(result.closed)
            self.check_closed(file, True)

    def compare_gzip_file_contents(self, a, b):
        # The first 10 bytes of a gzip header include a timestamp. The header
        # can be followed by other "volatile" metadata, so only compare gzip
        # footers (last 8 bytes) which contain a CRC-32 checksum and the length
        # of the uncompressed data.
        self.assertEqual(a[-8:], b[-8:])

    def test_open_binary(self):
        self.check_open_state_contents(self.binary_file, self.binary_contents,
                                       True, encoding='binary',
                                       compression=None)

        self.assertEqual(self.get_contents(self.binary_file),
                         self.binary_contents)

    def test_open_gzip(self):
        self.check_open_state_contents(self.gzip_file, self.text_contents,
                                       False, compression='gzip')

        self.compare_gzip_file_contents(self.get_contents(self.gzip_file),
                                        self.gzip_contents)

    def test_open_bz2(self):
        self.check_open_state_contents(self.bz2_file, self.text_contents,
                                       False, compression='bz2')

        self.assertEqual(self.get_contents(self.bz2_file),
                         self.bz2_contents)

    def test_open_encoding(self):
        self.check_open_state_contents(self.big5_file, self.decoded_contents,
                                       False, encoding='big5')

        self.assertEqual(self.get_contents(self.big5_file),
                         self.encoded_contents)

    def test_open_gzip_encoding(self):
        self.check_open_state_contents(self.gzip_encoded_file,
                                       self.decoded_contents, False,
                                       compression='gzip', encoding='big5')

        self.compare_gzip_file_contents(
            self.get_contents(self.gzip_encoded_file),
            self.gzip_encoded_contents)

    def test_open_bz2_encoding(self):
        self.check_open_state_contents(self.bz2_encoded_file,
                                       self.decoded_contents, False,
                                       compression='bz2', encoding='big5')

        self.assertEqual(self.get_contents(self.bz2_encoded_file),
                         self.bz2_encoded_contents)


class WritableSourceTest(unittest.TestCase):
    def setUp(self):
        self._dir = tempfile.mkdtemp()

        with io.open(get_data_path('example_file'), mode='rb') as f:
            self.binary_contents = f.read()
        self.binary_file = self._make_file('example_file')

        with io.open(get_data_path('big5_file'), mode='rb') as f:
            self.encoded_contents = f.read()
        self.big5_file = self._make_file('big5_file')

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
        self.safe_close(self.binary_file)
        self.safe_close(self.gzip_file)
        self.safe_close(self.bz2_file)
        self.safe_close(self.big5_file)
        self.safe_close(self.gzip_encoded_file)
        self.safe_close(self.bz2_encoded_file)

    def safe_close(self, f):
        if hasattr(f, 'close'):
            f.close()

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


@unittest.skipIf(not has_httpretty, "HTTPretty not available to mock tests.")
class TestReadURL(ReadableBinarySourceTests, ReadableSourceTest):
    expected_close = True

    def setUp(self):
        super(TestReadURL, self).setUp()
        httpretty.enable()

        for file in (get_data_path('example_file'),
                     get_data_path('big5_file'),
                     get_data_path('example_file.gz'),
                     get_data_path('example_file.bz2'),
                     get_data_path('big5_file.gz'),
                     get_data_path('big5_file.bz2')):

            with io.open(file, mode='rb') as f:
                httpretty.register_uri(httpretty.GET, self.get_fileobj(file),
                                       body=f.read(),
                                       content_type="application/octet-stream")

    def tearDown(self):
        super(TestReadURL, self).setUp()
        httpretty.disable()

    def get_fileobj(self, path):
        return "http://example.com/" + os.path.split(path)[1]


class TestReadBytesIO(ReadableBinarySourceTests, ReadableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        with io.open(path, mode='rb') as f:
            return io.BytesIO(f.read())


class TestWriteBytesIO(WritableBinarySourceTests, WritableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        return io.BytesIO()

    def get_contents(self, file):
        return file.getvalue()

    def test_open_gzip(self):
        self.check_open_state_contents(self.gzip_file, self.text_contents,
                                       False, compression='gzip')

        self.compare_gzip_file_contents(self.get_contents(self.gzip_file),
                                        self.gzip_contents)

    def test_open_gzip_encoding(self):
        self.check_open_state_contents(self.gzip_encoded_file,
                                       self.decoded_contents, False,
                                       compression='gzip', encoding='big5')

        self.compare_gzip_file_contents(
            self.get_contents(self.gzip_encoded_file),
            self.gzip_encoded_contents)


class TestReadBufferedReader(ReadableBinarySourceTests, ReadableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        return io.open(path, mode='rb')


class TestWriteBufferedReader(WritableBinarySourceTests, WritableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        return io.open(path, mode='w+b')

    def get_contents(self, file):
        file.close()
        with io.open(file.name, mode='rb') as f:
            return f.read()


class TestReadNamedTemporaryFile(ReadableBinarySourceTests,
                                 ReadableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        fileobj = tempfile.NamedTemporaryFile()
        with io.open(path, mode='rb') as fh:
            fileobj.write(fh.read())
            fileobj.flush()
            fileobj.seek(0)
        return fileobj


class TestWriteNamedTemporaryFile(WritableBinarySourceTests,
                                  WritableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        return tempfile.NamedTemporaryFile()

    def get_contents(self, file):
        file.flush()
        file.seek(0)
        contents = file.read()
        file.close()
        return contents


class TestReadTemporaryFile(ReadableBinarySourceTests, ReadableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        fileobj = tempfile.TemporaryFile()
        with io.open(path, mode='rb') as fh:
            fileobj.write(fh.read())
            fileobj.flush()
            fileobj.seek(0)
        return fileobj


class TestWriteTemporaryFile(WritableBinarySourceTests, WritableSourceTest):
    expected_close = False

    def get_fileobj(self, path):
        return tempfile.TemporaryFile()

    def get_contents(self, file):
        file.flush()
        file.seek(0)
        contents = file.read()
        file.close()
        return contents


class TestIterableReaderWriter(unittest.TestCase):
    def test_open(self):
        def gen():
            yield from ('a', 'b', 'c')
        list_ = list(gen())

        for input_ in gen(), list_:
            with skbio.io.open(input_) as result:
                self.assertIsInstance(result, io.TextIOBase)
                self.assertEqual(result.read(), 'abc')

    def test_open_with_newline(self):
        lines = ['a\r', 'b\r', 'c\r']
        with skbio.io.open(lines, newline='\r') as result:
            self.assertIsInstance(result, io.TextIOBase)
            self.assertEqual(result.readlines(), lines)

    def test_open_invalid_iterable(self):
        with self.assertRaises(skbio.io.IOSourceError):
            skbio.io.open([1, 2, 3])

    def test_open_empty_iterable(self):
        with skbio.io.open([]) as result:
            self.assertIsInstance(result, io.TextIOBase)
            self.assertEqual(result.read(), '')

    def test_open_write_mode(self):
        lines = []
        with skbio.io.open(lines, mode='w') as fh:
            fh.write('abc')
        self.assertEqual(lines, ['abc'])

        lines = []
        with skbio.io.open(lines, mode='w', newline='\r') as fh:
            fh.write('ab\nc\n')
        self.assertEqual(lines, ['ab\r', 'c\r'])

        self.assertTrue(fh.closed)
        fh.close()
        self.assertTrue(fh.closed)


if __name__ == '__main__':
    unittest.main()
