# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from six import string_types, text_type

import io
import gzip
import bz2file
from tempfile import gettempdir

import requests
from cachecontrol import CacheControl
from cachecontrol.caches import FileCache

from ._fileobject import (ReadableBufferedIO, ReadableTextIO,
                          IterableStringWriterIO, IterableStringReaderIO,
                          WrappedBufferedRandom)


def get_io_sources():
    return (
        # The order of these source is significant as they will short-circuit
        HTTPSource,
        FilePathSource,
        BytesIOSource,
        BufferedIOSource,
        TextIOSource,
        ReadableSource,
        IterableSource,
    )


def _compressors():
    return (
        GzipCompressor,
        BZ2Compressor
    )


def get_compression_handler(name):
    compressors = {c.name: c for c in _compressors()}
    compressors['auto'] = AutoCompressor
    return compressors.get(name, False)


class IOSource(object):
    closeable = True

    def __init__(self, file, options):
        self.file = file
        self.options = options

    def can_read(self):
        return False

    def can_write(self):
        return False

    def get_reader(self):
        raise NotImplementedError()

    def get_writer(self):
        return NotImplementedError()


class Compressor(IOSource):
    streamable = True
    name = ''

    def can_write(self):
        return True


class FilePathSource(IOSource):
    def can_read(self):
        return isinstance(self.file, string_types)

    def can_write(self):
        return self.can_read()

    def get_reader(self):
        return io.open(self.file, mode='rb')

    def get_writer(self):
        return io.open(self.file, mode='wb')


class HTTPSource(IOSource):
    def can_read(self):
        return (
            isinstance(self.file, string_types) and
            requests.compat.urlparse(self.file).scheme in {'http', 'https'})

    def get_reader(self):
        sess = CacheControl(requests.Session(),
                            cache=FileCache(gettempdir()))
        req = sess.get(self.file)

        # if the response is not 200, an exception will be raised
        req.raise_for_status()

        return io.BufferedReader(io.BytesIO(req.content))


class BytesIOSource(IOSource):
    closeable = False

    def can_read(self):
        return isinstance(self.file, io.BytesIO)

    def can_write(self):
        return self.can_read()

    def get_reader(self):
        return WrappedBufferedRandom(self.file)

    def get_writer(self):
        return self.get_reader()


class BufferedIOSource(IOSource):
    closeable = False

    def can_read(self):
        # `peek` is part of the API we want to guarantee, so we can't just look
        # for io.BufferedIOBase. Despite the fact that the C implementation of
        # io.BufferedRandom inherits io.BufferedReader/Writer it is not
        # reflected in an isinstance check, so we need to check for it manually
        return isinstance(self.file, (io.BufferedReader, io.BufferedRandom))

    def can_write(self):
        return isinstance(self.file, (io.BufferedWriter, io.BufferedRandom))

    def get_reader(self):
        return self.file

    def get_writer(self):
        return self.file


class TextIOSource(IOSource):
    closeable = False

    def can_read(self):
        return isinstance(self.file, io.TextIOBase) and self.file.readable()

    def can_write(self):
        return isinstance(self.file, io.TextIOBase) and self.file.writable()

    def get_reader(self):
        return self.file

    def get_writer(self):
        return self.file


class IterableSource(IOSource):
    def can_read(self):
        if hasattr(self.file, '__iter__'):
            iterator = iter(self.file)
            head = next(iterator)
            if isinstance(head, text_type):
                self.repaired = self._repair_iterable(head, iterator)
                return True
            else:
                # We may have mangled a generator at this point, so just abort
                raise Exception()
        return False

    def can_write(self):
        return hasattr(self.file, 'append') and hasattr(self.file, '__iter__')

    def get_reader(self):
        return IterableStringReaderIO(self.repaired)

    def get_writer(self):
        return IterableStringWriterIO(self.file)

    def _repair_iterable(self, head, tail):
        yield head
        for head in tail:
            yield head


class ReadableSource(IOSource):
    closeable = False

    def can_read(self):
        return hasattr(self.file, 'read')

    def get_reader(self):
        buffer_size = io.DEFAULT_BUFFER_SIZE

        raw = self.file.read(buffer_size)
        file = self._repair_readable(raw, self.file)
        if isinstance(raw, text_type):
            return ReadableTextIO(file, newline=self.options['newline'])
        elif isinstance(raw, bytes):
            return ReadableBufferedIO(file, buffer_size=buffer_size)
        else:
            raise Exception()

    def _repair_readable(self, raw, file):
        class Readable(object):
            def __init__(self):
                self._raw_read = False

            def read(self, b):
                if not self._raw_read:
                    self._raw_read = True
                    return raw
                else:
                    return file.read(b)

            def close(self):
                if hasattr(file, 'close'):
                    file.close()

        return Readable()


class GzipCompressor(Compressor):
    name = 'gzip'
    streamable = True

    def can_read(self):
        return self.file.peek(2)[:2] == b'\x1f\x8b'

    def get_reader(self):
        return gzip.GzipFile(fileobj=self.file)

    def get_writer(self):
        return gzip.GzipFile(fileobj=self.file,
                             compresslevel=self.options['compresslevel'])


class BZ2Compressor(Compressor):
    name = 'bz2'
    streamable = False

    def can_read(self):
        return self.file.peek(3)[:3] == b'BZh'

    def get_reader(self):
        return bz2file.BZ2File(self.file, mode='rb')

    def get_writer(self):
        return bz2file.BZ2File(self.file, mode='wb',
                               compresslevel=self.options['compresslevel'])


class AutoCompressor(Compressor):
    streamable = True  # We can' write so it doesn't matter
    name = 'auto'

    def get_reader(self):
        for compression_handler in _compressors():
            compressor = compression_handler(self.file, self.options)
            if compressor.can_read():
                return compressor.get_reader()

        return self.file

    def get_writer(self):
        return self.file


# class IOMiddleware(object):
#     def __init__(self, file):
#         self.file = file
#
#     def can_resolve(self):
#         return False
#
#     def read(self):
#         raise NotImplementedError("")
#
#
# class IOSource(IOMiddleware):
#     def can_resolve(self, mode):
#         return False
#
#     def should_close(self, mode):
#         return True
#
#     def encoding(self):
#         return None
#
#     def write(self, ext):
#         raise NotImplementedError("")
#
# class Compressor(IOMiddleware):
#     ext = ''
#     name = ''
#
#     def write(self, compresslevel):
#         raise NotImplementedError("")
#
#
# class PassThroughSource(IOSource):
#     def should_close(self, mode):
#         return False
#
#     def can_resolve(self, mode):
#         t = type(self.file)
#         if t == io.BufferedRandom:
#             return True
#         elif mode == 'r'
#             return t == io.BufferedReader
#         else:
#             return t == io.BufferedWriter
#
#     def read(self):
#         return self.file
#
#     def write(self):
#         return self.file
#
#
# class FilePathSource(IOSource):
#     def can_resolve(self, mode):
#         return isinstance(self.file, string_types)
#
#     def read(self):
#         return io.open(self.file, mode='rb')
#
#     def write(self, ext):
#         return io.open(self.file + ext, mode='wb')
#
# class _RequestsFileLike(object):
#     def __init__(self, response):
#         self.iterator = None
#
#     def read(self, n):
#         if self.iterator is None:
#             # It is unfortunate that we cannot change the chunk_size each
#             # time,
#             # but it shouldn't be needed in practice anyways.
#             self.iterator = response.iter_content(chunk_size=n,
#                                                   decode_unicode=False)
#         try:
#             return next(self.iterator)
#         except StopIteration:
#             return b''
#
# class HTTPSource(IOSource):
#     def __init__(self, file):
#         super(URLSource, self).__init__(file)
#         self.resolveable = (isinstance(self.file, string_types) and
#                             requests.compat.urlparse(self.file).scheme in
#                             {'http', 'https'})
#         if self.resolveable:
#             self.response = requests.get(self.file, stream=True)
#
#     def can_resolve(self, mode):
#         return mode == 'r' and self.resolveable
#
#     def read(self):
#         self.response.raise_for_status()
#         return FileObjectSource(_RequestsFileLike(self.response)).read()
#
#     def encoding(self):
#         return self.response.encoding
#
#
# class _TempFile(io.FileIO):
#     def __init__(self, mode='r+'):
#         fd, self.path = tempfile.mkstemp()
#         super(TempFile, self).__init__(fd, mode=mode, closefd=True)
#
#     def close(self):
#         super(TempFile, self).close()
#         os.unlink(self.path)
#
# class _WrappedBufferedWriter(io.BufferedRandom):
#     def __init__(self, file, buffer_size=io.DEFAULT_BUFFER_SIZE):
#         self._buffer_size = buffer_size
#         self._user_file = file
#         super(WrappedBufferedWriter, self).__init__(_TempFile(),
#                                                     buffer_size=buffer_size)
#
#     def close(self):
#         self.seek(0)
#         still_content = True
#         while still_content:
#             raw = self.read(self._buffer_size)
#             self._user_file.write(raw)
#             still_content = len(raw) == self._buffer_size
#
#         super(WrappedBufferedWriter, self).close()
#
#
# class FileObjectSource(IOSource):
#     def __init__(self, file):
#         super(FileObjectSource, self).__init__(file)
#         self._encoding = None
#
#     def can_resolve(self, mode):
#         return hasattr(self.file, 'read' if mode == 'r' else 'write')
#
#     def read(self, buffer_size=io.DEFAULT_BUFFER_SIZE):
#         raw = self.file.read(buffer_size + 1)
#         t = type(raw)
#         if t == text_type:
#             raw = raw.decode('utf-8')
#             self._encoding = 'utf-8'
#         elif t != bytes:
#             raise Exception()
#
#         if len(raw) > buffer_size > 0:
#             newfile = TempFile()
#             still_content = True
#             while still_content:
#                 newfile.write(raw)
#                 raw = self.file.read(buffer_size)
#                 if t == text_type:
#                     raw = raw.decode('utf-8')
#                 still_content = len(raw) > 0
#         else:
#             newfile = io.BytesIO(raw)
#
#         return io.BufferedReader(newfile)
#
#     def encoding(self):
#         return self._encoding
#
#     def write(self):
#         return _WrappedBufferedWriter(self.file)
