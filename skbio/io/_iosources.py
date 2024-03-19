# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import gzip
import bz2
import tempfile
import itertools

import requests
import h5py

from skbio.io import IOSourceError
from ._fileobject import (
    IterableStringWriterIO,
    IterableStringReaderIO,
    WrappedBufferedRandom,
)


# NamedTemporaryFile isn't an actual file class, it is a function which
# returns _TemporaryFileWrapper around a normal file object. Instead of
# relying on this implementation, we take whatever the class of the result of
# NamedTemporaryFile is.
with tempfile.NamedTemporaryFile() as fh:
    _WrappedTemporaryFile = type(fh)


def get_io_sources():
    return (
        # The order of these source is significant as they will short-circuit
        HTTPSource,
        FilePathSource,
        BytesIOSource,
        BufferedIOSource,
        TextIOSource,
        WrappedTemporaryFileSource,
        IterableSource,
    )


def _containers():
    return (HDF5Container, )


def _compressors():
    return (GzipCompressor, BZ2Compressor)


def get_compression_handler(name):
    compressors = {c.name: c for c in _compressors()}
    compressors["auto"] = AutoCompressor
    return compressors.get(name, False)


def get_container_handler(name):
    containers = {c.name: c for c in _containers()}
    containers["auto"] = AutoContainer
    return containers.get(name, False)


class IOSource:
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
        raise NotImplementedError()


class Compressor(IOSource):
    streamable = True
    name = ""

    def can_write(self):
        return True


class Container(IOSource):
    streamable = False
    name = ""

    def can_write(self):
        return True


class FilePathSource(IOSource):
    def can_read(self):
        return isinstance(self.file, str)

    def can_write(self):
        return self.can_read()

    def get_reader(self):
        return io.open(self.file, mode="rb")

    def get_writer(self):
        return io.open(self.file, mode="wb")


class HTTPSource(IOSource):
    def can_read(self):
        return isinstance(self.file, str) and requests.compat.urlparse(
            self.file
        ).scheme in {"http", "https"}

    def get_reader(self):
        req = requests.get(self.file)

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


class WrappedTemporaryFileSource(IOSource):
    closeable = False

    def can_read(self):
        return isinstance(self.file, _WrappedTemporaryFile) and self.file.readable()

    def can_write(self):
        return isinstance(self.file, _WrappedTemporaryFile) and self.file.writable()

    def get_reader(self):
        # _TemporaryFileWrapper has a file attribute which is an actual fileobj
        return self.file.file

    def get_writer(self):
        return self.file.file


class IterableSource(IOSource):
    def can_read(self):
        if hasattr(self.file, "__iter__"):
            iterator = iter(self.file)
            head = next(iterator, None)
            if head is None:
                self.repaired = []
                return True
            if isinstance(head, str):
                self.repaired = itertools.chain([head], iterator)
                return True
            else:
                # We may have mangled a generator at this point, so just abort
                raise IOSourceError(
                    "Could not open source: %r (mode: %r)"
                    % (self.file, self.options["mode"])
                )
        return False

    def can_write(self):
        return hasattr(self.file, "append") and hasattr(self.file, "__iter__")

    def get_reader(self):
        return IterableStringReaderIO(self.repaired, newline=self.options["newline"])

    def get_writer(self):
        return IterableStringWriterIO(self.file, newline=self.options["newline"])


class GzipCompressor(Compressor):
    name = "gzip"
    streamable = True

    def can_read(self):
        return self.file.peek(2)[:2] == b"\x1f\x8b"

    def get_reader(self):
        return gzip.GzipFile(fileobj=self.file)

    def get_writer(self):
        return gzip.GzipFile(
            fileobj=self.file, mode="wb", compresslevel=self.options["compresslevel"]
        )


class BZ2Compressor(Compressor):
    name = "bz2"
    streamable = False

    def can_read(self):
        return self.file.peek(3)[:3] == b"BZh"

    def get_reader(self):
        return bz2.BZ2File(self.file, mode="rb")

    def get_writer(self):
        return bz2.BZ2File(
            self.file, mode="wb", compresslevel=self.options["compresslevel"]
        )


class AutoCompressor(Compressor):
    streamable = True  # We can' write so it doesn't matter
    name = "auto"

    def get_reader(self):
        for compression_handler in _compressors():
            compressor = compression_handler(self.file, self.options)
            if compressor.can_read():
                return compressor.get_reader()

        return self.file

    def get_writer(self):
        return self.file


class HDF5Container(Container):
    name = 'hdf5'
    streamable = False

    def can_read(self):
        # From https://en.wikipedia.org/wiki/Hierarchical_Data_Format
        # Note that Wikipedia specifies: "\211HDF\r\n\032\n" which is an ordinal form:
        # >>> ord('\211')
        # 137
        # >>> ord('\x89')
        # 137
        # >>> ord('\032')
        # 26
        # >>> ord('\x1a')
        # 26
        return self.file.peek(8) == b"\x89HDF\r\n\x1a\n"

    # The following modes are permitted:
    # r, r+, w, w-, x, a
    # https://github.com/h5py/h5py/blob/d051d24a02e739bf51e744360e309974f6593c07/h5py/_hl/files.py#L263  # noqa
    def get_reader(self):
        return H5PYFile(self.file, mode='r')

    def get_writer(self):
        return H5PYFile(self.file, mode='w')


class H5PYMixin:
    @property
    def closed(self):
        return self.id.valid != 1

    def seek(self, *args, **kwargs):
        return 0

    def seekable(self):
        # assume seekable if open
        return not self.closed

    def readable(self):
        if not self.closed:
            if self.mode in ('r', 'r+', 'a'):
                return True
        return False

    def writable(self):
        if not self.closed:
            if self.mode in ('r+', 'w', 'w-', 'x', 'a'):
                return True
        return False


class H5PYFile(h5py.File, H5PYMixin):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class AutoContainer(Container):
    streamable = True
    name = "auto"

    def get_reader(self):
        for container_handler in _containers():
            container = container_handler(self.file, self.options)
            if container.can_read():
                return container.get_reader()

        return self.file

    def get_writer(self):
        return self.file
