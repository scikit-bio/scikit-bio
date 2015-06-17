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
from __future__ import absolute_import, division, print_function

import io
from contextlib2 import contextmanager, ExitStack
from functools import wraps

from skbio.io._iosources import get_io_sources, get_compression_handler
from skbio.io._fileobject import (is_text_file, is_binary_file,
                                  SaneTextIOWrapper)

_d = dict(mode='r', buffer_size=io.DEFAULT_BUFFER_SIZE, encoding=None,
          errors=None, newline=None, compression='auto', compresslevel=9,
          is_binary_file=None)

def resolve(file, mode=_d['mode'], buffer_size=_d['buffer_size'],
            encoding=_d['encoding'], errors=_d['errors'],
            newline=_d['newline'], compression=_d['compression'],
            compresslevel=_d['compresslevel']):
    arguments = locals().copy()

    if mode not in {'r', 'w'}:
        raise ValueError("Unsupported mode: %r, use 'r' or 'w'" % mode)

    newfile = None
    source = None
    for source_handler in get_io_sources():
        source = source_handler(file, arguments)
        if mode == 'r' and source.can_read():
            newfile = source.get_reader()
            break
        elif mode == 'w' and source.can_write():
            newfile = source.get_writer()
            break

    if not newfile:
        raise Exception("Could not open source: %r (mode: %r)" % (file, mode))

    return newfile, source, is_binary_file(newfile)


def open(file, mode=_d['mode'], buffer_size=_d['buffer_size'],
         encoding=_d['encoding'], errors=_d['errors'], newline=_d['newline'],
         compression=_d['compression'], compresslevel=_d['compresslevel']):
    arguments = locals().copy()
    del arguments['file']

    file, _, is_binary_file = resolve(**arguments)

    return _munge_file(file, is_binary_file, **arguments)


def _munge_file(file, is_binary_file, **arguments):
    mode = arguments.get('mode', _d['mode'])
    encoding = arguments.get('encoding', _d['encoding'])
    errors = arguments.get('errors', _d['errors'])
    newline = arguments.get('newline', _d['newline'])
    compression = arguments.get('compression', _d['compression'])
    is_output_binary = encoding == 'binary'
    newfile = file

    compression_handler = get_compression_handler(compression)

    if is_output_binary and (errors is not None or newline is not None):
        raise ValueError("Cannot use `force_encoding`, `errors`,"
                         " or `newline` with binary encoding.")

    if compression and not compression_handler:
     raise ValueError("Unsupported compression: %r" % compression)

    if is_binary_file:
        if compression:
            c = compression_handler(file, arguments)
            newfile = c.get_writer() if mode == 'w' else c.get_reader()

        if not is_output_binary:
            newfile = SaneTextIOWrapper(newfile, encoding=encoding,
                                        errors=errors, newline=newline)

    else:
        if compression and compression != 'auto':
            raise ValueError("Cannot use compression with that source.")
        if is_output_binary:
            raise ValueError("Source is not a binary source")

    return newfile

@contextmanager
def resolve_file(file, **kwargs):
    file, source, is_binary_file = resolve(file, **kwargs)
    try:
        yield file, is_binary_file
    finally:
        if source.closeable:
            if hasattr(file, 'partial_close'):
                file.partial_close()
            else:
                file.close()

@contextmanager
def open_file(file, **kwargs):
    with resolve_file(file, **kwargs) as (file, is_binary_file):
        yield _munge_file(file, is_binary_file, **kwargs)

@contextmanager
def open_files(files, **kwargs):
    with ExitStack() as stack:
        yield [stack.enter_context(open_file(f, **kwargs)) for f in files]

#
# def open(file, mode='r', buffer_size=io.DEFAULT_BUFFER_SIZE, encoding='auto',
#          force_encoding=True, errors=None, newline=None, compression='auto',
#          compress_level=9, compress_ext=True, **kwargs):
#     arguments = locals().copy()
#
#     if mode not in {'r', 'w'}:
#         raise ValueError("Unsupported mode: %r" % mode)
#
#     is_binary = mode[-1] == 'b'
#     skip = False
#
#     if is_binary and (encoding is not None or errors is not None
#                       or newline is not None or not force_encoding):
#         raise ValueError("Cannot use `encoding`, `force_encoding`, `errors`,"
#                          " or `newline` when in binary mode (%r)" % mode)
#
#     mode = mode[0]
#     compression_handler = get_compressor_lookup().get(compression, False)
#
#     if compression and not compression_handler:
#         raise ValueError("Unsupported compression: %r" % compression)
#
#     ext = compression_handler.ext if compression else ''
#
#     file, closeable, options = _resolve_source(file, mode, arguments)
#     if options['encoding'] != encoding and (encoding is None or not
#                                             force_encoding):
#         encoding = options['encoding']
#
#     if is_text_file(file):
#         if compression and compression != 'auto':
#             raise ValueError("Cannot use compression with that source.")
#         if is_binary:
#             raise ValueError("Source is not a binary source")
#         skip = True
#     elif not is_binary_file(file):
#         raise ValueError("Unkown source")
#
#     if compression and not skip:
#         c = compression_handler(file, arguments)
#         if mode == 'w':
#             file = c.get_writer()
#         elif c.can_read():
#             file = c.get_reader()
#         else:
#             raise ValueError("")
#
#     if not is_binary and not skip:
#         file = io.TextIOWrapper(file, encoding=encoding, errors=errors,
#                                 newline=newline)
#
#
#     return_extra = kwargs.pop('_ret_extra', False)
#     if return_extra:
#         return file, owned, skip
#
#     for key in kwargs:
#         raise TypeError("open() got an unexpected keyword argument %r" % key)
#
#     return file


# def _get_filehandle(filepath_or, mode, encoding, newline):
#     """Open file if `filepath_or` looks like a string/unicode/bytes, else
#     pass through.
#     """
#     if encoding == "binary":
#         if mode[-1] == 't':
#             raise ValueError("binary encoding cannot be used with text-mode")
#         if newline is not None:
#             raise ValueError("binary encoding cannot be used with newline")
#         mode += 'b'
#
#     if isinstance(filepath_or, string_types):
#         if requests.compat.urlparse(filepath_or).scheme in {'http', 'https'}:
#             sess = CacheControl(requests.Session(),
#                                 cache=FileCache(gettempdir()))
#             req = sess.get(filepath_or, **kwargs)
#
#             # if the response is not 200, an exception will be raised
#             req.raise_for_status()
#             if encoding == "binary":
#                 return io.BytesIO(req.content), True
#             else:
#                 return io.StringIO(req.content.decode(encoding),
#                                    newline=newline), True
#         else:
#             return io.open(filepath_or, mode=mode, encoding=encoding,
#                            newline=newline), True
#     else:
#         return filepath_or, False
#
#
# @contextmanager
# def open_file(filepath_or, mode='r', encoding=None, newline=None):
#     """Context manager, like ``open``, but lets file handles and file like
#     objects pass untouched.
#
#     It is useful when implementing a function that can accept both
#     strings and file-like objects (like numpy.loadtxt, etc), with the
#     additional benefit that it can load data from an HTTP/HTTPS URL.
#
#     Parameters
#     ----------
#     filepath_or : str/bytes/unicode string or file-like
#         If ``filpath_or`` is a file path to be opened the ``open`` function is
#         used and a filehandle is returned. If ``filepath_or`` is a string that
#         refers to an HTTP or HTTPS URL, a GET request is created and a BytesIO
#         object is returned with the contents of the URL. Else, if a file-like
#         object is passed, the object is returned untouched.
#
#     Other parameters
#     ----------------
#     args, kwargs : tuple, dict
#         When `filepath_or` is a string, any extra arguments are passed
#         on to the ``open`` builtin. If `filepath_or` is a URL, then only kwargs
#         are passed into `requests.get`.
#
#     Notes
#     -----
#     When files are retrieved from a URL, they are cached in disk inside a
#     temporary directory as generated by tempfile.gettempdir.
#
#     Examples
#     --------
#     >>> with open_file('filename') as f:                      # doctest: +SKIP
#     ...     pass
#     >>> fh = open('filename')                                 # doctest: +SKIP
#     >>> with open_file(fh) as f:                              # doctest: +SKIP
#     ...     pass
#     >>> fh.closed                                             # doctest: +SKIP
#     False
#     >>> fh.close()                                            # doctest: +SKIP
#     >>> with open_file('http://example.com/file.fasta') as f: # doctest: +SKIP
#     ...     pass
#
#     See Also
#     --------
#     requests.get
#
#     """
#     fh, own_fh = _get_filehandle(filepath_or, mode, encoding, newline)
#     try:
#         yield fh
#     finally:
#         if own_fh:
#             fh.close()
#
#
# @contextmanager
# def open_files(fp_list, mode='r', encoding=None, newline=None):
#     fhs, owns = zip(*[_get_filehandle(f, mode, encoding, newline) for f in
#                       fp_list])
#     try:
#         yield fhs
#     finally:
#         for fh, is_own in zip(fhs, owns):
#             if is_own:
#                 fh.close()
