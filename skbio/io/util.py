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

    open
    open_file
    open_files
    resolve
    resolve_file

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

from skbio.io._iosources import get_io_sources, get_compression_handler
from skbio.io._fileobject import (
    is_binary_file, SaneTextIOWrapper, CompressedBufferedReader,
    CompressedBufferedWriter)

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

    if newfile is None:
        raise Exception("Could not open source: %r (mode: %r)" % (file, mode))

    return newfile, source, is_binary_file(newfile)


def open(file, mode=_d['mode'], buffer_size=_d['buffer_size'],
         encoding=_d['encoding'], errors=_d['errors'], newline=_d['newline'],
         compression=_d['compression'], compresslevel=_d['compresslevel']):
    arguments = locals().copy()
    del arguments['file']

    file, _, is_binary_file = resolve(file, **arguments)
    return _munge_file(file, is_binary_file, arguments)


def _munge_file(file, is_binary_file, arguments):
    mode = arguments.get('mode', _d['mode'])
    encoding = arguments.get('encoding', _d['encoding'])
    errors = arguments.get('errors', _d['errors'])
    newline = arguments.get('newline', _d['newline'])
    compression = arguments.get('compression', _d['compression'])
    is_output_binary = encoding == 'binary'
    newfile = file

    compression_handler = get_compression_handler(compression)

    if is_output_binary and (errors is not _d['errors'] or
                             newline is not _d['errors']):
        raise ValueError("Cannot use `errors` or `newline` with binary"
                         " encoding.")

    if compression is not None and not compression_handler:
        raise ValueError("Unsupported compression: %r" % compression)

    if is_binary_file:
        if compression:
            c = compression_handler(newfile, arguments)
            if mode == 'w':
                newfile = CompressedBufferedWriter(file, c.get_writer())
            else:
                newfile = CompressedBufferedReader(file, c.get_reader())

        if not is_output_binary:
            newfile = SaneTextIOWrapper(newfile, encoding=encoding,
                                        errors=errors, newline=newline)
    else:
        if compression is not None and compression != 'auto':
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
            # TODO: Add comment
            if hasattr(file, 'partial_close'):
                file.partial_close()
            else:
                file.close()


@contextmanager
def open_file(file, **kwargs):
    with resolve_file(file, **kwargs) as (file, is_binary_file):
        file = _munge_file(file, is_binary_file, kwargs)
        try:
            yield file
        finally:
            # As soon as we leave the above context manager file will be closed
            # It is important to realize that because we are closing an inner
            # buffer, the outer buffer will reflect that state, but it won't
            # get flushed as the inner buffer is oblivious to the outer
            # buffer's existence.
            if not file.closed:
                file.flush()


@contextmanager
def open_files(files, **kwargs):
    with ExitStack() as stack:
        yield [stack.enter_context(open_file(f, **kwargs)) for f in files]
