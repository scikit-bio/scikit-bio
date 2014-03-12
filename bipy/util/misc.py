#!/usr/bin/env python
r"""
Misc (:mod:`bipy.util.misc`)
============================

.. currentmodule:: bipy.util.misc

This module provides miscellaneous useful utility classes and methods that do
not fit in any specific module.

Functions
---------

.. autosummary::
   :toctree: generated/

   safe_md5

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
    http://stackoverflow.com/questions/1131220/get-md5-hash-of-big-files-in-python/1131255

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
