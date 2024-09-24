r"""I/O utilities (:mod:`skbio.io.util`)
====================================

.. currentmodule:: skbio.io.util

This module provides utility functions to deal with files and I/O in general.

Functions
---------

.. autosummary::
   :toctree:

   open
   open_file
   open_files

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from contextlib import contextmanager, ExitStack
from types import MethodType

from skbio.io import IOSourceError
from skbio.io._iosources import get_io_sources, get_compression_handler
from skbio.io._fileobject import (
    is_binary_file,
    SaneTextIOWrapper,
    CompressedBufferedReader,
    CompressedBufferedWriter,
)
import skbio

_d = dict(
    mode="r",
    encoding=None,
    errors=None,
    newline=None,
    compression="auto",
    compresslevel=9,
)


def _resolve(
    file,
    mode=_d["mode"],
    encoding=_d["encoding"],
    errors=_d["errors"],
    newline=_d["newline"],
    compression=_d["compression"],
    compresslevel=_d["compresslevel"],
):
    arguments = locals().copy()

    if mode not in {"r", "w"}:
        raise ValueError("Unsupported mode: %r, use 'r' or 'w'" % mode)

    newfile = None
    source = None
    for source_handler in get_io_sources():
        source = source_handler(file, arguments)
        if mode == "r" and source.can_read():
            newfile = source.get_reader()
            break
        elif mode == "w" and source.can_write():
            newfile = source.get_writer()
            break

    if newfile is None:
        raise IOSourceError("Could not open source: %r (mode: %r)" % (file, mode))

    return newfile, source, is_binary_file(newfile)


def open(
    file,
    mode=_d["mode"],
    encoding=_d["encoding"],
    errors=_d["errors"],
    newline=_d["newline"],
    compression=_d["compression"],
    compresslevel=_d["compresslevel"],
):
    r"""Convert input into a filehandle.

    Supported inputs:

    +--------------------------------------+--------+---------+-----------+
    | type                                 | can \  | can \   | source \  |
    |                                      | read   | write   | type      |
    +======================================+========+=========+===========+
    | file path                            | True   | True    | Binary    |
    +--------------------------------------+--------+---------+-----------+
    | URL                                  | True   | False   | Binary    |
    +--------------------------------------+--------+---------+-----------+
    | ``["lines list\n"]``                 | True   | True    | Text      |
    +--------------------------------------+--------+---------+-----------+
    | :class:`io.StringIO`                 | True   | True    | Text      |
    +--------------------------------------+--------+---------+-----------+
    | :class:`io.BytesIO`                  | True   | True    | Binary    |
    +--------------------------------------+--------+---------+-----------+
    | :class:`io.TextIOWrapper`            | True   | True    | Text      |
    +--------------------------------------+--------+---------+-----------+
    | :class:`io.BufferedReader`           | True   | False   | Binary    |
    +--------------------------------------+--------+---------+-----------+
    | :class:`io.BufferedWriter`           | False  | True    | Binary    |
    +--------------------------------------+--------+---------+-----------+
    | :class:`io.BufferedRandom`           | True   | True    | Binary    |
    +--------------------------------------+--------+---------+-----------+
    | :func:`tempfile.TemporaryFile`       | True   | True    | Binary    |
    +--------------------------------------+--------+---------+-----------+
    | :func:`tempfile.NamedTemporaryFile`  | True   | True    | Binary    |
    +--------------------------------------+--------+---------+-----------+

    .. note:: When reading a list of unicode (str) lines, the input for
       `newline` is used to determine the number of lines in the resulting file
       handle, not the number of elements in the list. This is to allow
       composition with ``file.readlines()``.

    Parameters
    ----------
    file : filepath, url, filehandle, list
        The input to convert to a filehandle.
    mode : {'r', 'w'}, optional
        Whether to return a readable or writable file. Conversely, this does
        not imply that the returned file will be unwritable or unreadable.
        To get a binary filehandle set `encoding` to binary.
    encoding : str, optional
        The encoding scheme to use for the file. If set to 'binary', no bytes
        will be translated. Otherwise this matches the behavior of
        :func:`io.open`.
    errors : str, optional
        Specifies how encoding and decoding errors are to be handled. This has
        no effect when `encoding` is binary (as there can be no errors).
        Otherwise this matches the behavior of :func:`io.open`.
    newline : {None, "", '\\n', '\\r\\n', '\\r'}, optional
        Matches the behavior of :func:`io.open`.
    compression : {'auto', 'gzip', 'bz2', None}, optional
        Will compress or decompress `file` depending on `mode`. If 'auto' then
        determining the compression of the file will be attempted and the
        result will be transparently decompressed. 'auto' will do nothing
        when writing. Other legal values will use their respective compression
        schemes. `compression` cannot be used with a text source.
    compresslevel : int (0-9 inclusive), optional
        The level of compression to use, will be passed to the appropriate
        compression handler. This is only used when writing.

    Returns
    -------
    filehandle : io.TextIOBase or io.BufferedReader/Writer
        When `encoding='binary'` an :class:`io.BufferedReader` or
        :class:`io.BufferedWriter` will be returned depending on `mode`.
        Otherwise an implementation of :class:`io.TextIOBase` will be returned.

        .. note:: Any underlying resources needed to create `filehandle` are
           managed transparently. If `file` was closeable, garbage collection
           of `filehandle` will not close `file`. Calling `close` on
           `filehandle` will close `file`. Conversely calling `close` on `file`
           will cause `filehandle` to reflect a closed state. **This does not
           mean that a `flush` has occured for `filehandle`, there may still
           have been data in its buffer! Additionally, resources may not have
           been cleaned up properly, so ALWAYS call `close` on `filehandle` and
           NOT on `file`.**

    """
    arguments = locals().copy()
    del arguments["file"]

    file, _, is_binary_file = _resolve(file, **arguments)
    return _munge_file(file, is_binary_file, arguments)


def _munge_file(file, is_binary_file, arguments):
    mode = arguments.get("mode", _d["mode"])
    encoding = arguments.get("encoding", _d["encoding"])
    errors = arguments.get("errors", _d["errors"])
    newline = arguments.get("newline", _d["newline"])
    compression = arguments.get("compression", _d["compression"])
    is_output_binary = encoding == "binary"
    newfile = file

    compression_handler = get_compression_handler(compression)

    if is_output_binary and newline is not _d["newline"]:
        raise ValueError("Cannot use `newline` with binary encoding.")

    if compression is not None and not compression_handler:
        raise ValueError("Unsupported compression: %r" % compression)

    if is_binary_file:
        if compression:
            c = compression_handler(newfile, arguments)
            if mode == "w":
                newfile = CompressedBufferedWriter(
                    file, c.get_writer(), streamable=c.streamable
                )
            else:
                newfile = CompressedBufferedReader(file, c.get_reader())

        if not is_output_binary:
            newfile = SaneTextIOWrapper(
                newfile, encoding=encoding, errors=errors, newline=newline
            )
    else:
        if compression is not None and compression != "auto":
            raise ValueError("Cannot use compression with that source.")
        if is_output_binary:
            raise ValueError("Source is not a binary source")

    return newfile


@contextmanager
def _resolve_file(file, **kwargs):
    file, source, is_binary_file = _resolve(file, **kwargs)
    try:
        yield file, source, is_binary_file
    finally:
        if source.closeable:
            file.close()


@contextmanager
def open_file(file, **kwargs):
    r"""Context manager for :func:`skbio.io.util.open`.

    The signature matches :func:`open`. This context manager will not close
    filehandles that it did not create itself.

    Examples
    --------
    Here our input isn't a filehandle and so `f` will get closed.

    >>> with open_file(['a\n']) as f:
    ...     f.read()
    ...
    'a\n'
    >>> f.closed
    True

    Here we provide an open file and so `f` will not get closed and neither
    will `file`.

    >>> file = io.BytesIO(b'BZh91AY&SY\x03\x89\x0c\xa6\x00\x00\x01\xc1\x00\x00'
    ...                   b'\x108\x00 \x00!\x9ah3M\x1c\xb7\x8b\xb9"\x9c(H\x01'
    ...                   b'\xc4\x86S\x00')
    >>> with open_file(file) as f:
    ...     f.read()
    ...
    'a\nb\nc\n'
    >>> f.closed
    False
    >>> file.closed
    False

    """
    with _resolve_file(file, **kwargs) as (file, source, is_binary_file):
        newfile = _munge_file(file, is_binary_file, source.options)
        try:
            yield newfile
        finally:
            # As soon as we leave the above context manager file will be closed
            # It is important to realize that because we are closing an inner
            # buffer, the outer buffer will reflect that state, but it won't
            # get flushed as the inner buffer is oblivious to the outer
            # buffer's existence.
            if not newfile.closed:
                newfile.flush()
                _flush_compressor(newfile)


def _flush_compressor(file):
    if isinstance(file, io.TextIOBase) and hasattr(file, "buffer"):
        file = file.buffer
    if isinstance(file, CompressedBufferedWriter) and not file.streamable:
        # Some formats like BZ2 compress the entire file, and so they will
        # only flush once they have been closed. These kinds of files do not
        # close their underlying buffer, but only testing can prove that...
        file.raw.close()


@contextmanager
def open_files(files, **kwargs):
    """Context manager for :func:`skbio.io.util.open`.

    A plural form of :func:`open_file`.

    """
    with ExitStack() as stack:
        yield [stack.enter_context(open_file(f, **kwargs)) for f in files]


class ReadWriteDescriptor:
    """A descriptor class to generate read/write methods and their documentation."""

    def __init__(self, method):
        """Initialize the descriptor with the method (read/write)."""
        # self.cls_name = cls_name
        self.method = method

    def __get__(self, instance, owner):
        """Return the generated method (read/write) with its documentation."""
        if self.method == "read":
            return self._generate_read_method(owner)
        elif self.method == "write":
            if instance is None:
                return self._generate_write_method(owner)
            else:
                return MethodType(self._generate_write_method(owner), instance)

    def _generate_read_method(self, cls):
        def read_method(file, format=None, **kwargs):
            return skbio.io.read(file, into=cls, format=format, **kwargs)

        read_method.__doc__ = self._make_docstring(cls, "read")

        return read_method

    def _generate_write_method(self, cls):
        def write_method(self, file, format=None, **kwargs):
            if format is None:
                if hasattr(cls, "default_write_format"):
                    format = cls.default_write_format
                else:
                    raise ValueError(f"{cls.__name__} has no default write format.")
            return skbio.io.write(self, into=file, format=format, **kwargs)

        write_method.__doc__ = self._make_docstring(cls, "write")

        return write_method

    def _make_docstring(self, cls, method):
        """
        Generate a docstring for the read or write method.

        Parameters:
        cls (type): The class for which the method docstring is being generated.
        method (str): The method for which the docstring is generated
        ('read' or 'write').

        Returns:
        str: The dynamically generated docstring.
        """
        if method == "read":
            doc = f"""Create a new ``{cls.__name__}`` instance from a file.

This is a convenience method for :func:`skbio.io.read`. For more information about the
I/O system in scikit-bio, please see :mod:`skbio.io`.

Supported file formats include:

something like cls.list_read_formats

Parameters
----------
file : openable (filepath, URL, filehandle, etc.)
    The location to read the given `format` into. Something that is understood by
    :func"`skbio.io.util.open`. Filehandles are not automatically closed, it is the
    responsibility of the caller.
format : str, optional
    The format of the file. If None, the format will be inferred.
kwargs : dict, optional
    Additional arguments passed to the reader for {cls.__name__}.

Returns
-------
{cls.__name__}
    A new instance of {cls.__name__}.

See Also
--------
skbio.io.read
"""
        elif method == "write":
            doc = f"""Write the {cls.__name__} instance to a file.

This is a convenience method for :func:`skbio.io.write`.

Parameters
----------
file : openable (filepath, URL, filehandle, etc.)
    The location to write the given `format`.
format : str, optional
    The format to write the {cls.__name__} object as.
kwargs : dict, optional
    Additional arguments passed to the writer for {cls.__name__}.

Returns
-------
None

See Also
--------
skbio.io.write
"""
        return doc
