# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import io
import tempfile
import os

import bz2file


def is_binary_file(file):
    return isinstance(file, (io.BufferedReader, io.BufferedWriter,
                             io.BufferedRandom))

# Everything beyond this point will be some kind of hack needed to make
# everything work. It's not pretty and it doesn't make great sense much
# of the time. I am very sorry to the poor soul who has to read beyond.


class StringIO(io.StringIO):
    """Treat Bytes the same as Unicode by decoding ascii, for testing only."""
    def __init__(self, string=None, **kwargs):
        if isinstance(string, bytes):
            string = string.decode()
        super(StringIO, self).__init__(string, **kwargs)


class SaneTextIOWrapper(io.TextIOWrapper):
    def __init__(self, *args, **kwargs):
        super(SaneTextIOWrapper, self).__init__(*args, **kwargs)
        self._should_close_buffer = True

    def __del__(self):
        # Accept the inevitability of the buffer being closed by the destructor
        # because of this line in Python 2.7:
        # https://github.com/python/cpython/blob/2.7/Modules/_io/iobase.c#L221
        self._should_close_buffer = False
        # Actually close for Python 3 because it is an override.
        # We can't call super because Python 2 doesn't actually
        # have a `__del__` method for IOBase (hence this
        # workaround). Close is idempotent so it won't matter
        # that Python 2 will end up calling this twice
        self.close()

    def close(self):
        # We can't stop Python 2.7 from calling close in the deconstructor
        # so instead we can prevent the buffer from being closed with a flag.

        # Based on:
        # https://github.com/python/cpython/blob/2.7/Lib/_pyio.py#L1586
        if self.buffer is not None and not self.closed:
            try:
                self.flush()
            finally:
                if self._should_close_buffer:
                    self.buffer.close()


class CompressedMixin(object):
    """Act as a bridge between worlds"""
    def __init__(self, before_file, *args, **kwargs):
        super(CompressedMixin, self).__init__(*args, **kwargs)
        self._should_close_raw = True
        self._before_file = before_file

    def __del__(self):
        self._should_close_raw = False
        self.close()

    @property
    def closed(self):
        return self.raw.closed or self._before_file.closed

    # Based on:
    # https://github.com/python/cpython/blob/2.7/Lib/_pyio.py#L732
    def close(self):
        if self.raw is not None and not self.closed:
            try:
                # may raise BlockingIOError or BrokenPipeError etc
                self.flush()
            finally:
                if self._should_close_raw:
                    self.raw.close()
                    # The above will not usually close the before_file
                    # We want the decompression to be transparent, so we don't
                    # want users to deal with this edge case. Instead we can
                    # just close the original now that we are being closed.
                    self._before_file.close()


class CompressedBufferedReader(CompressedMixin, io.BufferedReader):
    pass


class CompressedBufferedWriter(CompressedMixin, io.BufferedWriter):
    def flush(self):
        super(CompressedBufferedWriter, self).flush()
        self.raw.flush()


class BZ2File(bz2file.BZ2File):
    def flush(self):
        # HACK because flush does not actually work.
        # based on
        # https://github.com/nvawda/bz2file/blob/master/bz2file.py#L129
        super(BZ2File, self).flush()
        if self._mode == bz2file._MODE_WRITE:
            with self._lock:
                self._fp.write(self._compressor.flush())


class TemporaryFile(io.FileIO):
    """This exists because tempfile.TemporaryFile is not composable with io"""
    def __init__(self, mode='r+'):
        fd, self._path = tempfile.mkstemp()
        super(TemporaryFile, self).__init__(fd, mode=mode, closefd=True)

    def __del__(self):
        if not self.closed:
            self.close()

    def close(self):
        try:
            super(TemporaryFile, self).close()
        finally:
            os.unlink(self._path)


class IterableStringReaderIO(io.StringIO):
    def __init__(self, iterable, newline=None):
        self._iterable = iterable
        super(IterableStringReaderIO, self).__init__(u''.join(iterable),
                                                     newline=newline)


class IterableStringWriterIO(IterableStringReaderIO):
    def close(self):
        if not self.closed:
            backup = self.tell()
            self.seek(0)
            for line in self:
                self._iterable.append(line)
            self.seek(backup)
        super(IterableStringWriterIO, self).close()


class ReadableBufferedIO(io.BufferedReader):
    def __init__(self, file, buffer_size=io.DEFAULT_BUFFER_SIZE):
        self._file = file
        file = TemporaryFile()

        has_content = True
        while has_content:
            raw = self._file.read(buffer_size)
            file.write(raw)
            has_content = len(raw) > 0
        file.seek(0)
        super(ReadableBufferedIO, self).__init__(file)

    def close(self):
        if hasattr(self._file, 'close'):
            self._file.close()
        self.partial_close()

    def partial_close(self):
        super(ReadableBufferedIO, self).close()


class ReadableTextIO(io.TextIOWrapper):
    def __init__(self, file, newline=None):
        self._file = ReadableBufferedIO(_InlineUTF8Decoder(file))
        super(ReadableTextIO, self).__init__(self._file, encoding='utf-8',
                                             newline=newline)

    def partial_close(self):
        self._file.partial_close()


class _InlineUTF8Decoder(object):
    def __init__(self, file):
        self.file = file

    def read(self, b):
        return self.file.read(b).encode('utf-8')

    def close(self):
        if hasattr(self.file, 'close'):
            self.file.close()
