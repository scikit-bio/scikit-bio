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

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from future.builtins import bytes, str

from contextlib import contextmanager


def _is_string_or_bytes(s):
    """Returns True if input argument is string (unicode or not) or bytes.
    """
    return isinstance(s, str) or isinstance(s, bytes)


def _get_filehandle(filepath_or, *args, **kwargs):
    """Open file if `filepath_or` looks like a string/unicode/bytes, else
    pass through.
    """
    if _is_string_or_bytes(filepath_or):
        fh, own_fh = open(filepath_or, *args, **kwargs), True
    else:
        fh, own_fh = filepath_or, False
    return fh, own_fh


@contextmanager
def open_file(filepath_or, *args, **kwargs):
    """Context manager, like ``open``, but lets file handles and file like
    objects pass untouched.

    It is useful when implementing a function that can accept both
    strings and file-like objects (like numpy.loadtxt, etc).

    Parameters
    ----------
    filepath_or : str/bytes/unicode string or file-like
         If string, file to be opened using ``open``. Else, it is returned
         untouched.

    Other parameters
    ----------------
    args, kwargs : tuple, dict
        When `filepath_or` is a string, any extra arguments are passed
        on to the ``open`` builtin.

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

    """
    fh, own_fh = _get_filehandle(filepath_or, *args, **kwargs)
    try:
        yield fh
    finally:
        if own_fh:
            fh.close()
