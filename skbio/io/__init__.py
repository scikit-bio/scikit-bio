r"""
File I/O (:mod:`skbio.io`)
==========================

.. currentmodule:: skbio.io

This package provides I/O functionality for skbio.

Introduction to I/O
-------------------
Reading and writing files (I/O) can be a complicated task:

* A file format can sometimes be read into more than one in-memory
  representation (i.e., object). For example, a FASTQ file can be read into an
  :mod:`skbio.alignment.SequenceCollection` or :mod:`skbio.alignment.Alignment`
  depending on the file's contents and what operations you'd like to perform on
  your data.
* A single object might be writeable to more than one file format. For example,
  an :mod:`skbio.alignment.Alignment` object could be written to FASTA, FASTQ,
  QSEQ,PHYLIP, or Stockholm formats, just to name a few.
* You might not know the exact file format of your file, but you want to read
  it into an appropriate object.
* You want to read multiple files into a single object, or write an object to
  multiple files.
* Instead of reading a file into an object, you might want to stream the file
  using a generator (e.g., if the file cannot be fully loaded into memory).

To address these issues (and others), scikit-bio provides a simple, powerful
interface for dealing with I/O. We accomplish this by using a single I/O
registry. Below is a description of how to use the registry and how to extend
it.

Reading Files Into scikit-bio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are two ways to read files. The first way is to use the
procedural interface:

``my_obj = skbio.io.read(<filehandle or filepath>, format='<format here>',
into=<class to construct>)``

The second is to use the object-oriented (OO) interface which is automatically
constructed from the procedural interface:

``my_obj = <class to construct>.read(<filehandle or filepath>,
format='<format here>')``

For example, to read a `newick` file using both interfaces you would type:

>>> from skbio.io import read
>>> from skbio import TreeNode
>>> from io import StringIO
>>> open_filehandle = StringIO(u'(a, b);')
>>> tree = read(open_filehandle, format='newick', into=TreeNode)
>>> tree
<TreeNode, name: unnamed, internal node count: 0, tips count: 2>

For the OO interface:

>>> _ = open_filehandle.seek(0)
>>> tree2 = TreeNode.read(open_filehandle, format='newick')
>>> tree2
<TreeNode, name: unnamed, internal node count: 0, tips count: 2>

In the case of ``skbio.io.read`` if `into` is not provided, then a generator
will be returned. What the generator yields will depend on what format is being
read.

When `into` is provided, format may be omitted and the registry will use its
knowledge of the available formats for the requested class to infer the correct
format. This format inference is also available in the OO interface, meaning
that `format` may be omitted there as well.

As an example:

>>> _ = open_filehandle.seek(0)
>>> tree3 = TreeNode.read(open_filehandle)
>>> tree3
<TreeNode, name: unnamed, internal node count: 0, tips count: 2>

We call format inference `sniffing`, much like the
`csv <https://docs.python.org/2/library/csv.html#csv.Sniffer>`_ module of
Python's standard library. The goal of a `sniffer` is twofold: to identify if a
file is a specific format, and if it is, to provide `**kwargs` which can be
used to better parse the file.

.. note:: There is a built-in `sniffer` which results in a useful error message
   if an empty file is provided as input and the format was omitted.

Writing Files from scikit-bio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Just as when reading files, there are two ways to write files.

Procedural Interface:

``skbio.io.write(my_obj, format='<format here>',
into=<filehandle or filepath>)``

OO Interface:

``my_obj.write(<filehandle or filepath>, format='<format here>')``

In the procedural interface, `format` is required. Without it, scikit-bio does
not know how you want to serialize an object. OO interfaces define a default
`format`, so it may not be necessary to include it.

Supported File Formats
^^^^^^^^^^^^^^^^^^^^^^
For details on what objects are supported by each format,
see the associated documentation.

.. autosummary::
   :toctree: generated/

   clustal
   dm
   fasta
   newick
   ordres
   phylip

Formats are considered to be names which represent a way of encoding a file.
A simple format is just a single name as a string, such as ``'newick'``.
In some cases, objects can be constructed from more than one format, such as
``'fasta'`` and ``'qual'`` for reading/writing FASTA/qual files into an
:mod:`skbio.alignment.SequenceCollection`. In these cases, we use what is
called a *compound* format. It can be written like this:
``format=['fasta', 'qual']``. We also support the shorthand:
``format='fasta, qual'``. In these cases, where you would put a filehandle or a
filepath string, you will replace it with a list of filehandles and/or
filepaths which correspond to the order of your compound format
(``[<filehandle or filepath of fasta>, <filehandle or filepath of qual>]``).
The order does not matter for the compound format, as long as the files are
provided in that same order, so we could have used ``['qual', 'fasta']``
instead as long as the quality file came first.

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
   ClustalFormatError
   DMFormatError
   FASTAFormatError
   NewickFormatError
   OrdResFormatError
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
class. Please see the relevant documenation for more information about these
functions and for the contract that is expected of `readers`, `writers`, and
`sniffers`.

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

.. note:: Keyword arguments are not permitted in `sniffers`. `Sniffers` may not
   raise exceptions; if an exception is thrown by a `sniffer`, the user will be
   asked to report it on our `issue tracker
   <https://github.com/biocore/scikit-bio/issues/>`_.

When raising errors, the error should be a subclass of ``FileFormatError``
specific to your new format.

Writing Unit Tests
^^^^^^^^^^^^^^^^^^
Because scikit-bio handles all of the I/O boilerplate, you only need to test
the actual business logic of your `readers`, `writers`, and `sniffers`. The
easiest way to accomplish this is to create a list of files and their expected
results when deserialized. Then you can iterate through the list ensuring the
expected results occur and that the expected results can be reserialized into
an equivalent file. This process is called 'roundtripping'.

It is also important to test some invalid inputs and ensure that the correct
error is raised by your `readers`. Consider using `assertRaises` as a context
manager like so:

.. code-block:: python

   with self.assertRaises(SomeFileFormatErrorSubclass) as cm:
       do_something_wrong()
   self.assertIn('action verb or subject of an error', str(cm.exception))

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
                         FileFormatError, ClustalFormatError, DMFormatError,
                         FASTAFormatError, NewickFormatError,
                         OrdResFormatError, PhylipFormatError)
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
           'ClustalFormatError',
           'DMFormatError',
           'FASTAFormatError',
           'NewickFormatError',
           'OrdResFormatError',
           'PhylipFormatError']

# Necessary to import each file format module to have them added to the I/O
# registry. We use import_module instead of a typical import to avoid flake8
# unused import errors.
import_module('skbio.io.clustal')
import_module('skbio.io.dm')
import_module('skbio.io.fasta')
import_module('skbio.io.newick')
import_module('skbio.io.ordres')
import_module('skbio.io.phylip')

# Now that all of our I/O has loaded, we can add the object oriented methods
# (read and write) to each class which has registered I/O operations.
initialize_oop_interface()

from numpy.testing import Tester
test = Tester().test
