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

from contextlib import contextmanager
from io import TextIOWrapper
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
    elif isinstance(filepath_or, BytesIO):
        start_pos = filepath_or.tell()

        filepath_or.seek(0)
        is_gzip = filepath_or.read(2) == magic_number

        filepath_or.seek(start_pos)
    else:
        raise ValueError('Unsupported type {}'.format(type(filepath_or)))

    return is_gzip


def get_filehandle(filepath_or, mode='r', binary=False,
                   gzip=None, *args, **kwargs):
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
            fh = open(filepath_or, mode + 'b', *args, **kwargs)

        if gzip or (gzip is None and _is_gzip(fh)):
            fh = GzipFile(fileobj=fh)

        if not binary:
            fh = TextIOWrapper(fh)

        own_fh = True
    else:
        fh, own_fh = filepath_or, False
    return fh, own_fh


def get_filemode(fh):
    if _is_gzip(fh.name):
        if isinstance(fh, TextIOWrapper):
            # Text gzip case in Python 3.x.
            mode = fh.buffer.fileobj.mode
            mode = mode.replace('b', 't')
        elif isinstance(fh, GzipFile):
            # Binary/text gzip case in Python 2.x.
            mode = fh.fileobj.mode
            if mode == 'rtb':
                mode = 'rt'
        else:
            # Binary gzip case in Python 3.x.
            mode = fh.fileobj.mode
    else:
        mode = fh.mode

    return mode


@contextmanager
def open_file(filepath_or, *args, **kwargs):
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
    fh, own_fh = get_filehandle(filepath_or, *args, **kwargs)
    try:
        yield fh
    finally:
        if own_fh:
            fh.close()


@contextmanager
def open_files(fp_list, *args, **kwargs):
    fhs, owns = zip(*[get_filehandle(f, *args, **kwargs) for f in fp_list])
    try:
        yield fhs
    finally:
        for fh, is_own in zip(fhs, owns):
            if is_own:
                fh.close()
