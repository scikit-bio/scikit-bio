#!/usr/bin/env python
r"""
Miscellaneous Utilities (:mod:`bipy.util.misc`)
===============================================

.. currentmodule:: bipy.util.misc

This module provides miscellaneous useful utility classes and methods that do
not fit in any specific module.

Functions
---------

.. autosummary::
   :toctree: generated/

   safe_md5
   remove_files

"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import hashlib
from os import remove


def safe_md5(open_file, block_size=2**20):
    """Computes an md5 sum without loading the file into memory

    Parameters
    ----------
    open_file : file object
        open file handle to the archive to compute the checksum
    block_size : int, optional
        size of the block taken per iteration

    Returns
    -------

    md5 : md5 object from the hashlib module
        object with the loaded file

    Notes
    -----

    This method is based on the answers given in:
    http://stackoverflow.com/a/1131255/379593

    Examples
    --------

    >>> from StringIO import StringIO
    >>> from bipy.util.misc import safe_md5
    >>> fd = StringIO("foo bar baz") # open file like object
    >>> x = safe_md5(fd)
    >>> x.hexdigest()
    'ab07acbb1e496801937adfa772424bf7'
    >>> fd.close()
    """
    md5 = hashlib.md5()
    data = True
    while data:
        data = open_file.read(block_size)
        if data:
            md5.update(data)
    return md5


def remove_files(list_of_filepaths, error_on_missing=True):
    """Remove list of filepaths, optionally raising an error if any are missing

    Parameters
    ----------
    list_of_filepaths : list of strings
        list with filepaths to remove
    error_on_missing : bool, optional
        whether or not the function should raise an ``OSError`` if a file is
        not found

    Raises
    ------

    OSError
        If a filepath in the list does not exist


    Examples
    --------

    >>> from tempfile import NamedTemporaryFile
    >>> from os.path import exists
    >>> from bipy.util.misc import remove_files
    >>> h = NamedTemporaryFile(delete=False)
    >>> exists(h.name) # it exists
    True
    >>> remove_files([h.name])
    >>> exists(h.name) # and now it's gone
    False

    """
    missing = []
    for fp in list_of_filepaths:
        try:
            remove(fp)
        except OSError:
            missing.append(fp)

    if error_on_missing and missing:
        raise OSError("Some filepaths were not accessible: %s" %
                      '\t'.join(missing))
