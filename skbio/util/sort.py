#!/usr/bin/env python
r"""
Sort util functions (:mod:`skbio.util.sort`)
============================================

.. currentmodule:: skbio.util.sort

This module provides useful functions for sorting.

Functions
---------

.. autosummary::
   :toctree: generated/

   natsort
   signed_natsort

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from re import split


def _natsort_key(item, case_sensitivity=False):
    """Provides normalized version of item for sorting with digits.

    Adapted from:
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    """
    item = str(item)

    try:
        chunks = split('(\d+(?:\.\d+)?)', item)
    except TypeError:
        # if item is a tuple or list (i.e., indexable, but not a string)
        # work with the first element
        chunks = split('(\d+(?:\.\d+)?)', item[0])
    for ii in range(len(chunks)):
        if chunks[ii] and chunks[ii][0] in '0123456789':
            if '.' in chunks[ii]:
                numtype = float
            else:
                numtype = int
            # wrap in tuple with '0' to explicitly specify numbers come first
            chunks[ii] = (0, numtype(chunks[ii]))
        else:
            chunks[ii] = (1, chunks[ii])
    return (chunks, item)


def _natsort_key_case_insensitive(item):
    """Provides normalized version of item for sorting with digits and is case-
    insensitive.

    Adapted from:
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    """
    # added the lower() call to allow for case-insensitive sorting
    item = str(item).lower()

    try:
        chunks = split('(\d+(?:\.\d+)?)', item)
    except TypeError:
        # if item is a tuple or list (i.e., indexable, but not a string)
        # work with the first element
        chunks = split('(\d+(?:\.\d+)?)', item[0])
    for ii in range(len(chunks)):
        if chunks[ii] and chunks[ii][0] in '0123456789':
            if '.' in chunks[ii]:
                numtype = float
            else:
                numtype = int
            # wrap in tuple with '0' to explicitly specify numbers come first
            chunks[ii] = (0, numtype(chunks[ii]))
        else:
            chunks[ii] = (1, chunks[ii])
    return (chunks, item)


def natsort(seq, case_sensitive=True):
    """Sort an iterable in a way that humans expect.

    Adapted from:
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html

    Parameters
    ----------
    seq : iterable
        The sequence of objects to sort
    case_sensitive : bool, optional
        If true, performs case-sensitive sorting. Otherwise, it performs case-
        insensitive sorting. Default: True.

    Returns
    -------
    list
        The input sequence sorted
    """
    if case_sensitive:
        natsort_key = _natsort_key
    else:
        natsort_key = _natsort_key_case_insensitive

    alist = list(seq)
    alist.sort(key=natsort_key)

    return alist


def signed_natsort(data):
    """Sort an iterable considering the cases where elements are signed in a
    way that humans expect.

    The elements will be assumed to be real numbers, if that assumption fails,
    then the elements will be sorted using a natural sorting algorithm.

    Parameters
    ----------
    data: iterable
        When a string is provided, the string will try to be type-casted to a
        float type, if a tuple is provided, the first element will be used to
        sort the list. If a dict is provided a sorted version of the keys will
        be returned.

    Returns
    -------
    list
        sorted version of data
    """

    # list is empty, do nothing
    if not data:
        return data

    # deal with non-[tuple, dict, list] types of data
    if not all([isinstance(element, tuple) or isinstance(element, list) or
                isinstance(element, dict) for element in data]):
        try:
            return sorted(data, key=float)
        except ValueError:
            return natsort(data)

    # deal with tuples type of data, the first element can be a real number or
    # a string, the second element is a string that won't be accounted
    try:
        return sorted(data, key=lambda tup: float(tup[0]))
    except ValueError:
        return natsort(data)
