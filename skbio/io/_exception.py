from __future__ import absolute_import, division, print_function

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


class FormatIdentificationError(FileFormatError):
    """Exception raised when a file's format cannot be identified"""
    pass


class DMFormatError(FileFormatError):
    """Exception raised when a dissimilarity matrix file cannot be parsed."""
    pass


class DuplicateRegistrationError(Exception):
    """Exception raised a function is already registered in skbio.io"""
    pass
