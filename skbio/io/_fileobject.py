# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import io


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


class WrappedBufferedRandom(io.BufferedRandom):
    def __init__(self, *args, **kwargs):
        super(WrappedBufferedRandom, self).__init__(*args, **kwargs)
        self._should_close_raw = True

    def __del__(self):
        self._should_close_raw = False
        self.close()

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


class CompressedMixin(object):
    """Act as a bridge between worlds"""
    def __init__(self, before_file, *args, **kwargs):
        self.streamable = kwargs.pop('streamable', True)
        self._should_close_raw = True
        self._before_file = before_file
        super(CompressedMixin, self).__init__(*args, **kwargs)

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
    pass


class IterableStringReaderIO(io.StringIO):
    def __init__(self, iterable, newline):
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
