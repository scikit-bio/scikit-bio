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

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


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
