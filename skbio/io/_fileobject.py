import six

import io
import traceback
import tempfile
import os

def is_binary_file(file):
    return isinstance(file, (io.BufferedReader, io.BufferedWriter, io.BufferedRandom))

def is_text_file(file):
    return isinstance(file, io.TextIOBase)

class StringIO(io.StringIO):
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
        self.close()  # Actually close for Python 3 because it is an override.
                      # We can't call super because Python 2 doesn't actually
                      # have a `__del__` method for IOBase (hence this
                      # workaround). Close is idempotent so it won't matter
                      # that Python 2 will end up calling this twice

    def close(self):
        # We can't stop Python 2.7 from calling close in the deconstructor
        # so instead we can prevent the buffer from being closed with a flag.

        # Based on:
        # https://github.com/python/cpython/blob/2.7/Lib/_pyio.py#L1586
        # https://github.com/python/cpython/blob/3.4/Lib/_pyio.py#L1615
        if self.buffer is not None and not self.closed:
            try:
                self.flush()
            finally:
                if self._should_close_buffer:
                    self.buffer.close()


class TemporaryFile(io.FileIO):
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
