r"""
File IO (:mod:`skbio.io`)
=========================================

.. currentmodule:: skbio.io

This package provides general exception/warning definitions used throughout
scikit-bio, as well as several subpackages containing various utility
functionality, including I/O and unit-testing convenience functions.

Functions
---------

.. autosummary::
   :toctree: generated/

   write,
   read,
   guess_format,
   get_writer,
   get_reader,
   get_identifier,
   list_write_formats,
   list_read_formats,
   register_writer,
   register_reader,
   register_identifier

Exceptions
----------

.. autosummary::
   :toctree: generated/

   FileFormatError
   RecordError
   FieldError
   DuplicateRegistrationError
   FormatIdentificationError

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._exception import (DuplicateRegistrationError, RecordError, FieldError,
                         FormatIdentificationError, FileFormatError)

from ._registry import (write, read, guess_format, get_writer, get_reader,
                        get_identifier, list_write_formats, list_read_formats,
                        register_writer, register_reader, register_identifier)

__all__ = ['write', 'read', 'guess_format', 'get_writer', 'get_reader',
           'get_identifier', 'list_write_formats', 'list_read_formats',
           'register_writer', 'register_reader', 'register_identifier',
           'DuplicateRegistrationError', 'RecordError', 'FieldError',
           'FormatIdentificationError', 'FileFormatError']

from numpy.testing import Tester
test = Tester().test
