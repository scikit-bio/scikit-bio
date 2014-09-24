r"""
File I/O (:mod:`skbio.io`)
==========================

.. currentmodule:: skbio.io

This package provides I/O functionality for skbio.

Introduction to I/O
-------------------

Reading Files Into scikit-bio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are two ways to read files. The first way is to use the
imperative interface:

``my_obj = skbio.io.read(<filehandle or filepath>, format='<format here>',
into=<class to construct>)``

The second is to use the object-oriented interface which is dynamically
generated:

``my_obj = <class to construct>.read(<filehandle or filepath>,
format='<format here>')``

In the case of ``skbio.io.read`` if `into` is not provided, then a generator
will be returned. When `into` is provided, format may be ommitted and the
registry will use its knowledge of the available formats for the requested
class to infer the correct format. This format inferrence is also available in
the OOP interface, meaning that `format` may be ommitted there as well.

We call format inferrence: `sniffing`, much like the `csv` module of python's
standard library. The goal of a `sniffer` is twofold: to identify if a file
is a specific format, and if it is, to provide `**kwargs` which can be used
to better parse the file.

Writing Files from scikit-bio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Just as when reading files, there are two ways to write files.

Imperative Interface:

``skbio.io.write(my_obj, format='<format here>',
into=<filehandle or filepath>)``

OOP Interface:

``my_obj.write(<filehandle or filepath>, format='<format here>')``

In both interfaces, `format` is required. Without it, scikit-bio does not know
how you want to serialize an object.

Supported File Formats
^^^^^^^^^^^^^^^^^^^^^^
For details on what objects are supported by each format,
see the associated documentation.

.. autosummary::
   :toctree: generated/

   dm
   ordres
   newick
   phylip

User Functions
^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   write
   read
   sniff

User Exceptions
^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   FileFormatError
   RecordError
   FieldError
   UnrecognizedFormatError
   DMFormatError
   OrdResFormatError
   NewickFormatError
   PhylipFormatError

User Warnings
^^^^^^^^^^^^^

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
decorators are executed on importing the user functions above. Use the function
``import_module('skbio.io.my_new_format')``

The following keyword args may not be used when defining new `readers` or
`writers` as they already have special meaning to the registry system:
* `format`
* `into`
* `mode`
* `verify`

Developer Functions
^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^

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
                         NewickFormatError, PhylipFormatError)
from ._registry import (write, read, sniff, get_writer, get_reader,
                        get_sniffer, list_write_formats, list_read_formats,
                        register_writer, register_reader, register_sniffer,
                        initialize_oop_interface)

__all__ = ['write', 'read', 'sniff',
           'list_write_formats', 'list_read_formats',
           'get_writer', 'get_reader', 'get_sniffer',
           'register_writer', 'register_reader', 'register_sniffer',
           'initialize_oop_interface',

           'FormatIdentificationWarning', 'ArgumentOverrideWarning',

           'DuplicateRegistrationError', 'InvalidRegistrationError',
           'RecordError', 'FieldError', 'UnrecognizedFormatError',

           'FileFormatError',
           'DMFormatError',
           'OrdResFormatError',
           'NewickFormatError',
           'PhylipFormatError']

# Necessary to import each file format module to have them added to the I/O
# registry. We use import_module instead of a typical import to avoid flake8
# unused import errors.
import_module('skbio.io.dm')
import_module('skbio.io.ordres')
import_module('skbio.io.newick')
import_module('skbio.io.phylip')

# Now that all of our I/O has loaded, we can add the object oriented methods
# (read and write) to each class which has registered I/O operations.
initialize_oop_interface()

from numpy.testing import Tester
test = Tester().test
