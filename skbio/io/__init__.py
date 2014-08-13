r"""
File I/O (:mod:`skbio.io`)
=========================================

.. currentmodule:: skbio.io

This package provides I/O functionality for skbio.

For information about extending the I/O functionality of skbio, see the
associated Developer Documentation

User Functions
--------------

.. autosummary::
   :toctree: generated/

   write
   read
   guess_format

Developer Documentation
-----------------------
To extend I/O in skbio, developers should create a submodule in `skbio/io/`
named after the functionality it provides.

For example, if you were to create readers and writers for a `fasta` file, you
would create a submodule `skbio/io/fasta.py`.
In this submodule you would use the following decorator factories:
``register_writer``, ``register_reader``, and ``register_identifier``.
These associate your functionality to a format string and potentially an skbio
class.

Please see the relavant documenation for more information about these
functions.

Once you are satisfied with the functionality, you will need to ensure that
`skbio/io/__init__.py` contains an import of your new submodule, that way the
decorators are executed on importing the user functions above. We reccommend
importing functions that return a generator as they are generally useful and
make sense in the context of `skbio.io`.

Developer Functions
-------------------

.. autosummary::
    :toctree: generated/

    register_writer
    register_reader
    register_identifier
    list_write_formats
    list_read_formats
    get_writer
    get_reader
    get_identifier

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

from ._warning import UnprovenFormatWarning
from ._exception import (DuplicateRegistrationError, RecordError, FieldError,
                         FormatIdentificationError, FileFormatError)
from ._registry import (write, read, guess_format, get_writer, get_reader,
                        get_identifier, list_write_formats, list_read_formats,
                        register_writer, register_reader, register_identifier)

__all__ = ['write', 'read', 'guess_format', 'get_writer', 'get_reader',
           'get_identifier', 'list_write_formats', 'list_read_formats',
           'register_writer', 'register_reader', 'register_identifier',
           'DuplicateRegistrationError', 'RecordError', 'FieldError',
           'FormatIdentificationError', 'FileFormatError',
           'UnprovenFormatWarning']

from numpy.testing import Tester
test = Tester().test
