#!/usr/bin/env python
"""
Exceptions (:mod:`skbio.core.exception`)
========================================

.. currentmodule:: skbio.core.exception

This module defines custom exception classes used throughout the core
scikit-bio codebase.

Exceptions
----------

.. autosummary::
   :toctree: generated/

   FileFormatError
   RecordError
   FieldError
   BiologicalSequenceError
   SequenceCollectionError
   DissimilarityMatrixError
   DistanceMatrixError
   DissimilarityMatrixFormatError
   IDMismatchError
   MissingDataError
   MissingHeaderError
   MissingIDError
   TreeError
   NoLengthError
   DuplicateNodeError
   MissingNodeError
   NoParentError
   FileFormatError
   RecordError
   FastqParseError

"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


class FileFormatError(Exception):
    """Exception raised when a file can not be parsed."""
    pass


class RecordError(FileFormatError):
    """Exception raised when a record is bad."""
    pass


class FieldError(RecordError):
    """Exception raised when a field within a record is bad."""
    pass


class BiologicalSequenceError(Exception):
    """General error for biological sequence validation failures."""
    pass


class SequenceCollectionError(Exception):
    """General error for sequence collection validation failures."""
    pass


class DissimilarityMatrixError(Exception):
    """General error for dissimilarity matrix validation failures."""
    pass


class DistanceMatrixError(DissimilarityMatrixError):
    """General error for distance matrix validation failures."""
    pass


class MissingIDError(DissimilarityMatrixError):
    """Error for ID lookup that doesn't exist in the dissimilarity matrix."""

    def __init__(self, missing_id):
        super(MissingIDError, self).__init__()
        self.args = ("The ID '%s' is not in the dissimilarity matrix." %
                     missing_id,)


class DissimilarityMatrixFormatError(DissimilarityMatrixError):
    """Error for reporting issues in dissimilarity matrix file format.

    Typically used during parsing.

    """
    pass


class IDMismatchError(DissimilarityMatrixFormatError):
    """Error for reporting mismatch between IDs in a dissimilarity matrix file.

    Typically used during parsing.

    """

    def __init__(self, actual, expected):
        super(IDMismatchError, self).__init__()
        self.args = ("Encountered mismatched IDs while parsing the "
                     "dissimilarity matrix file. Found '%s' but expected "
                     "'%s'. Please ensure that the IDs match between the "
                     "dissimilarity matrix header (first row) and the row "
                     "labels (first column)." % (actual, expected),)


class MissingHeaderError(DissimilarityMatrixFormatError):
    """Error for reporting a missing ID header line during parsing."""

    def __init__(self):
        super(MissingHeaderError, self).__init__()
        self.args = ("Could not find a header line containing IDs in the "
                     "dissimilarity matrix file. Please verify that the file "
                     "is not empty.",)


class MissingDataError(DissimilarityMatrixFormatError):
    """Error for reporting missing data lines during parsing."""

    def __init__(self, actual, expected):
        super(MissingDataError, self).__init__()
        self.args = ("Expected %d row(s) of data, but found %d." % (expected,
                                                                    actual),)


class TreeError(Exception):
    """General tree error"""
    pass


class NoLengthError(TreeError):
    """Missing length when expected"""
    pass


class DuplicateNodeError(TreeError):
    """Duplicate nodes with identical names"""
    pass


class MissingNodeError(TreeError):
    """Expecting a node"""
    pass


class NoParentError(MissingNodeError):
    """Missing a parent"""
    pass


class FileFormatError(Exception):
    """Exception raised when a file can not be parsed."""
    pass


class RecordError(FileFormatError):
    """Exception raised when a record is bad."""
    pass


class FastqParseError(FileFormatError):
    pass
