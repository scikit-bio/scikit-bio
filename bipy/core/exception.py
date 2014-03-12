#!/usr/bin/env python
"""
Exceptions (:mod:`bipy.core.exception`)
=======================================

.. currentmodule:: bipy.core.exception

This module defines custom exception classes used throughout the core bipy
codebase.

Exceptions
----------

.. autosummary::
   :toctree: generated/

   BiologicalSequenceError
   DistanceMatrixError
   DistanceMatrixFormatError
   IDMismatchError
   MissingDataError
   MissingHeaderError
   MissingIDError

"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


class BiologicalSequenceError(Exception):
    pass


class DistanceMatrixError(Exception):
    """General error for distance matrix validation failures."""
    pass


class MissingIDError(Exception):
    """Error for ID lookup that doesn't exist in the distance matrix."""

    def __init__(self, missing_id):
        super(MissingIDError, self).__init__()
        self.args = ("The ID '%s' is not in the distance matrix." %
                     missing_id,)


class DistanceMatrixFormatError(Exception):
    """Error for reporting issues in distance matrix file format.

    Typically used during parsing.

    """
    pass


class IDMismatchError(Exception):
    """Error for reporting a mismatch between IDs.

    Typically used during parsing.

    """

    def __init__(self, actual, expected):
        super(IDMismatchError, self).__init__()
        self.args = ("Encountered mismatched IDs while parsing the distance "
                     "matrix file. Found '%s' but expected '%s'. Please "
                     "ensure that the IDs match between the distance matrix "
                     "header (first row) and the row labels (first column)." %
                     (actual, expected),)


class MissingHeaderError(Exception):
    """Error for reporting a missing ID header line during parsing."""

    def __init__(self):
        super(MissingHeaderError, self).__init__()
        self.args = ("Could not find a header line containing IDs in the "
                     "distance matrix file. Please verify that the file is "
                     "not empty.",)


class MissingDataError(Exception):
    """Error for reporting missing data lines during parsing."""

    def __init__(self, actual, expected):
        super(MissingDataError, self).__init__()
        self.args = ("Expected %d row(s) of data, but found %d." % (expected,
                                                                    actual),)
class FastqParseError(Exception):
    pass
    