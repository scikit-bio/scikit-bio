r"""
I/O utils (:mod:`skbio.io.util`)
================================

.. currentmodule:: skbio.io.util

This module provides utility functions to deal with files and I/O in
general.

Functions
---------

.. autosummary::
    :toctree: generated/

    open_file
    open_files

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from future.builtins import bytes, str
from six import BytesIO

import copy
import io
from contextlib import contextmanager
from gzip import GzipFile
from tempfile import gettempdir

import requests
from cachecontrol import CacheControl
from cachecontrol.caches import FileCache


def _is_string_or_bytes(s):
    """Returns True if input argument is string (unicode or not) or bytes.
    """
    return isinstance(s, str) or isinstance(s, bytes)


def _is_gzip(filepath_or):
    """Checks the first two bytes of the file for the gzip magic number
    If the first two bytes of the file are 1f 8b (the "magic number" of a
    gzip file), return True; otherwise, return false.
    """
    magic_number = b'\x1f\x8b'

    if hasattr(filepath_or, 'name'):
        filepath_or = filepath_or.name

    if _is_string_or_bytes(filepath_or):
        with open(filepath_or, 'rb') as file_:
            is_gzip = file_.read(2) == magic_number
    else:
        with reopen_file(filepath_or, binary=True, gzip=False) as file_:
            is_gzip = file_.read(2) == magic_number

    return is_gzip


def _get_filehandle(filepath_or, mode='r', binary=False, gzip=None,
                    compresslevel=None, encoding=None, **kwargs):
    """Open file if `filepath_or` looks like a string/unicode/bytes, else
    pass through.
    """

    assert mode in {'r', 'w'}

    if _is_string_or_bytes(filepath_or):
        if requests.compat.urlparse(filepath_or).scheme in {'http', 'https'}:
            if mode == 'w':
                raise ValueError('Writing not supported for URLs')

            sess = CacheControl(requests.Session(),
                                cache=FileCache(gettempdir()))
            req = sess.get(filepath_or, **kwargs)

            # if the response is not 200, an exception will be raised
            req.raise_for_status()

            fh = BytesIO(req.content)
        else:
            fh = open(filepath_or, mode + 'b', **kwargs)

        own_fh = True
    else:
        fh, own_fh = filepath_or, False

    # If file is binary, we can check to determine if file should
    # be opened as gzip and/or decoded into text.
    if _fileobj_is_binary(fh):
        if gzip or (gzip is None and _is_gzip(fh)):
            fh = GzipFile(fileobj=fh, compresslevel=compresslevel)

        if not binary:
            fh = io.TextIOWrapper(fh, encoding=encoding)

    return fh, own_fh


@contextmanager
def open_file(filepath_or, mode='r', binary=False, gzip=None,
              compresslevel=None, encoding=None, **kwargs):
    """Context manager, like ``open``, but lets file handles and file like
    objects pass untouched.

    It is useful when implementing a function that can accept both
    strings and file-like objects (like numpy.loadtxt, etc), with the
    additional benefit that it can load data from an HTTP/HTTPS URL
    and open gzip-compressed files transparently.

    Parameters
    ----------
    filepath_or : str/bytes/unicode string or file-like
        If ``filepath_or`` is a file path to be opened the ``open`` function
        is used and a filehandle is returned. If ``filepath_or`` is a file path
        that refers to a gzip-compressed file, ``gzip.open`` is used instead.
        If ``filepath_or`` is a string that refers to an HTTP or HTTPS URL, a
        GET request is created and a BytesIO object is returned with the
        contents of the URL. Else, if a file-like object is passed, the object
        is returned untouched.

    Other parameters
    ----------------
    args, kwargs : tuple, dict
        When `filepath_or` is a string, any extra arguments are passed
        on to the ``open`` builtin or to ``gzip.open`` if `file_path_or`
        refers to a gzip-compressed file. If `filepath_or` is a URL, then
        only kwargs are passed into `requests.get`.

    Notes
    -----
    When files are retrieved from a URL, they are cached in disk inside a
    temporary directory as generated by tempfile.gettempdir.

    Examples
    --------
    >>> with open_file('filename') as f:  # doctest: +SKIP
    ...     pass
    >>> fh = open('filename')             # doctest: +SKIP
    >>> with open_file(fh) as f:          # doctest: +SKIP
    ...     pass
    >>> fh.closed                         # doctest: +SKIP
    False
    >>> fh.close()                        # doctest: +SKIP
    >>> with open_file('http://foo.bar.com/file.fasta') as f: # doctest: +SKIP
    ...     pass

    See Also
    --------
    requests.get
    gzip.open

    """
    fh, own_fh = _get_filehandle(
        filepath_or, mode=mode, binary=binary, gzip=gzip,
        compresslevel=compresslevel, encoding=encoding, **kwargs)
    try:
        yield fh
    finally:
        if own_fh:
            fh.close()


@contextmanager
def open_files(fp_list, mode='r', binary=True, gzip=None,
               compresslevel=None, encoding=None, **kwargs):
    fhs, owns = zip(*[_get_filehandle(
            f, mode=mode, binary=binary, gzip=gzip,
            compresslevel=compresslevel, encoding=encoding, **kwargs)
        for f in fp_list])
    try:
        yield fhs
    finally:
        for fh, is_own in zip(fhs, owns):
            if is_own:
                fh.close()


@contextmanager
def reopen_file(filepath_or, *args, **kwargs):
    if _is_string_or_bytes(filepath_or):
        cfh, _ = _get_filehandle(filepath_or, *args, **kwargs)
    else:
        cfh = _fileobj_reopen(filepath_or)

    try:
        yield cfh
    finally:
        cfh.close()


def _fileobj_reopen(fileobj):
    # The reason we do a copy is because we need the sniffer to not
    # mutate the orginal file while guessing the format. The
    # naive solution would be to seek to 0 at the end, but that
    # would break an explicit offset provided by the user. Instead
    # we create a shallow copy which works out of the box for
    # file-like object, but does not work for real files. Instead
    # the name attribute is reused in open for a new filehandle.
    # Using seek and tell is not viable because in real files tell
    # reflects the position of the read-ahead buffer and not the
    # true offset of the iterator.

    filename = _fileobj_filename(fileobj)

    if filename is None:
        cfh = copy.copy(fileobj)
        cfh.seek(0)
        return cfh
    else:
        gzip = _fileobj_is_gzip(fileobj)
        binary = _fileobj_is_binary(fileobj)
        cfh, _ = _get_filehandle(
            filename, mode='r', gzip=gzip, binary=binary)
        return cfh


def _fileobj_is_binary(fileobj):
    if isinstance(fileobj, (io.TextIOWrapper, io.StringIO)):
        return False
    elif isinstance(fileobj, io.BytesIO):
        return True
    elif hasattr(fileobj, 'mode') and isinstance(fileobj.mode, str):
        return fileobj.mode[-1] == 'b'
    elif hasattr(fileobj, 'buffer'):
        return _fileobj_is_binary(fileobj.buffer)
    elif hasattr(fileobj, 'fileobj'):
        return _fileobj_is_binary(fileobj.fileobj)
    else:
        raise ValueError()


def _fileobj_is_gzip(fileobj):
    if isinstance(fileobj, GzipFile):
        return True
    elif hasattr(fileobj, 'buffer'):
        return _fileobj_is_gzip(fileobj.buffer)
    else:
        return False


def _fileobj_filename(fileobj):
    if hasattr(fileobj, 'filename'):
        return fileobj.filename
    elif hasattr(fileobj, 'name'):
        return fileobj.name
    elif hasattr(fileobj, 'buffer'):
        return _fileobj_filename(fileobj.buffer)
    else:
        return None
