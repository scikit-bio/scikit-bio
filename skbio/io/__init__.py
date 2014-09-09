r"""
File I/O (:mod:`skbio.io`)
==========================

.. currentmodule:: skbio.io

This package provides I/O functionality for skbio.

For information about extending the I/O functionality of skbio, see the
associated Developer Documentation.

Supported File Formats
----------------------

.. autosummary::
   :toctree: generated/

   dm
   ordres
   phylip

User Functions
--------------

.. autosummary::
   :toctree: generated/

   write
   read
   sniff

User Exceptions
---------------

.. autosummary::
   :toctree: generated/

   FileFormatError
   RecordError
   FieldError
   UnrecognizedFormatError
   DMFormatError
   OrdResFormatError
   PhylipFormatError

User Warnings
-------------

.. autosummary::
   :toctree: generated/

   FormatIdentificationWarning
   ArgumentOverrideWarning

Developer Documentation
-----------------------
To extend I/O in skbio, developers should create a submodule in `skbio/io/`
named after the file format it implements.

For example, if you were to create readers and writers for a `fasta` file, you
would create a submodule `skbio/io/fasta.py`.
In this submodule you would use the following decorators:
``register_writer``, ``register_reader``, and ``register_sniffer``.
These associate your functionality to a format string and potentially an skbio
class.

Please see the relevant documenation for more information about these
functions.

Once you are satisfied with the functionality, you will need to ensure that
`skbio/io/__init__.py` contains an import of your new submodule, that way the
decorators are executed on importing the user functions above. We recommend
importing functions that return a generator as they are generally useful and
make sense in the context of `skbio.io`.

Developer Functions
-------------------

.. autosummary::
    :toctree: generated/

    register_writer
    register_reader
    register_sniffer
    list_write_formats
    list_read_formats
    get_writer
    get_reader
    get_sniffer

Developer Exceptions
--------------------

.. autosummary::
   :toctree: generated/

   DuplicateRegistrationError
   InvalidRegistrationError

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from importlib import import_module

from ._warning import FormatIdentificationWarning, ArgumentOverrideWarning
from ._exception import (DuplicateRegistrationError, InvalidRegistrationError,
                         RecordError, FieldError, UnrecognizedFormatError,
                         FileFormatError, DMFormatError, OrdResFormatError,
                         PhylipFormatError)
from ._registry import (write, read, sniff, get_writer, get_reader,
                        get_sniffer, list_write_formats, list_read_formats,
                        register_writer, register_reader, register_sniffer)

__all__ = ['write', 'read', 'sniff', 'get_writer', 'get_reader',
           'get_sniffer', 'list_write_formats', 'list_read_formats',
           'register_writer', 'register_reader', 'register_sniffer',
           'DuplicateRegistrationError', 'InvalidRegistrationError',
           'RecordError', 'FieldError', 'UnrecognizedFormatError',
           'FileFormatError', 'DMFormatError', 'OrdResFormatError',
           'PhylipFormatError', 'FormatIdentificationWarning',
           'ArgumentOverrideWarning']

# Necessary to import each file format module to have them added to the I/O
# registry. We use import_module instead of a typical import to avoid flake8
# unused import errors.
import_module('skbio.io.dm')
import_module('skbio.io.ordres')
import_module('skbio.io.phylip')

from numpy.testing import Tester
test = Tester().test
