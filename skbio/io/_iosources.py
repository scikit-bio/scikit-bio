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
