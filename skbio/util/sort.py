#!/usr/bin/env python
r"""
Gradient analyses (:mod:`skbio.maths.gradient`)
===============================================

.. currentmodule:: skbio.maths.gradient

This module provides functionality for performing gradient analyses.

Functions
---------

.. autosummary::
   :toctree: generated/

   gradient_analysis

"""
from __future__ import division

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from re import split


def _natsort_key(item, case_sensitivity=False):
    """Provides normalized version of item for sorting with digits.

    From:
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
    """Provides normalized version of item for sorting with digits.

    From:
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
    """Sort a sequence of text strings in a reasonable order.

    From:
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    """
    if case_sensitive:
        natsort_key = _natsort_key
    else:
        natsort_key = _natsort_key_case_insensitive

    alist = list(seq)
    alist.sort(key=natsort_key)

    return alist


def signed_natsort(data):
    """sort an iterable considering the cases where elements are signed

    data: list of tuples (with two strings as elements) or strings. When a
    string is provided, the string will try to be type-casted to a float type,
    if a tuple is provided, the first element will be used to sort the list. If
    a dict is provided a sorted version of the keys will be returned.

    output: sorted version of data

    The elements will be assumed to be real numbers, if that assumption fails,
    then the elements will be sorted using a natural sorting algorithm.

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
