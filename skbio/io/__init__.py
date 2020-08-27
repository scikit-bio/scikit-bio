r"""
File I/O (:mod:`skbio.io`)
==========================

.. currentmodule:: skbio.io

This package provides I/O functionality for skbio.

Supported file formats
----------------------
For details on what objects are supported by each format,
see the associated documentation.

.. currentmodule:: skbio.io.format
.. autosummary::
   :toctree: generated/

   blast6
   blast7
   clustal
   embl
   fasta
   fastq
   genbank
   gff3
   lsmat
   newick
   ordination
   phylip
   qseq
   stockholm

.. currentmodule:: skbio.io.registry


User functions
--------------

.. autosummary::
   :toctree: generated/

   write
   read
   sniff

.. currentmodule:: skbio.io

User exceptions and warnings
----------------------------

.. autosummary::
   :toctree: generated/

   FormatIdentificationWarning
   ArgumentOverrideWarning
   UnrecognizedFormatError
   IOSourceError
   FileFormatError
   BLAST7FormatError
   ClustalFormatError
   EMBLFormatError
   FASTAFormatError
   FASTQFormatError
   GenBankFormatError
   GFF3FormatError
   LSMatFormatError
   NewickFormatError
   OrdinationFormatError
   PhylipFormatError
   QSeqFormatError
   QUALFormatError
   StockholmFormatError


Subpackages
-----------

.. autosummary::
   :toctree: generated/

   registry
   util

For developer documentation on extending I/O, see :mod:`skbio.io.registry`.

Introduction to I/O
-------------------
Reading and writing files (I/O) can be a complicated task:

* A file format can sometimes be read into more than one in-memory
  representation (i.e., object). For example, a FASTA file can be read into an
  :mod:`skbio.alignment.TabularMSA` or :mod:`skbio.sequence.DNA` depending on
  what operations you'd like to perform on your data.
* A single object might be writeable to more than one file format. For example,
  an :mod:`skbio.alignment.TabularMSA` object could be written to FASTA, FASTQ,
  CLUSTAL, or PHYLIP formats, just to name a few.
* You might not know the exact file format of your file, but you want to read
  it into an appropriate object.
* You might want to read multiple files into a single object, or write an
  object to multiple files.
* Instead of reading a file into an object, you might want to stream the file
  using a generator (e.g., if the file cannot be fully loaded into memory).

To address these issues (and others), scikit-bio provides a simple, powerful
interface for dealing with I/O. We accomplish this by using a single I/O
registry.

What kinds of files scikit-bio can use
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To see a complete list of file-like inputs that can be used for reading,
writing, and sniffing, see the documentation for :func:`skbio.io.util.open`.

Reading files into scikit-bio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are two ways to read files. The first way is to use the
procedural interface:

.. code-block:: python

   my_obj = skbio.io.read(file, format='someformat', into=SomeSkbioClass)

The second is to use the object-oriented (OO) interface which is automatically
constructed from the procedural interface:

.. code-block:: python

   my_obj = SomeSkbioClass.read(file, format='someformat')

For example, to read a `newick` file using both interfaces you would type:

>>> from skbio import read
>>> from skbio import TreeNode
>>> from io import StringIO
>>> open_filehandle = StringIO('(a, b);')
>>> tree = read(open_filehandle, format='newick', into=TreeNode)
>>> tree
<TreeNode, name: unnamed, internal node count: 0, tips count: 2>

For the OO interface:

>>> open_filehandle = StringIO('(a, b);')
>>> tree = TreeNode.read(open_filehandle, format='newick')
>>> tree
<TreeNode, name: unnamed, internal node count: 0, tips count: 2>

In the case of :func:`skbio.io.registry.read` if `into` is not provided, then a
generator will be returned. What the generator yields will depend on what
format is being read.

When `into` is provided, format may be omitted and the registry will use its
knowledge of the available formats for the requested class to infer the correct
format. This format inference is also available in the OO interface, meaning
that `format` may be omitted there as well.

As an example:

>>> open_filehandle = StringIO('(a, b);')
>>> tree = TreeNode.read(open_filehandle)
>>> tree
<TreeNode, name: unnamed, internal node count: 0, tips count: 2>

We call format inference `sniffing`, much like the :class:`csv.Sniffer`
class of Python's standard library. The goal of a `sniffer` is twofold: to
identify if a file is a specific format, and if it is, to provide `**kwargs`
which can be used to better parse the file.

.. note:: There is a built-in `sniffer` which results in a useful error message
   if an empty file is provided as input and the format was omitted.

Writing files from scikit-bio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Just as when reading files, there are two ways to write files.

Procedural Interface:

.. code-block:: python

   skbio.io.write(my_obj, format='someformat', into=file)

OO Interface:

.. code-block:: python

   my_obj.write(file, format='someformat')

In the procedural interface, `format` is required. Without it, scikit-bio does
not know how you want to serialize an object. OO interfaces define a default
`format`, so it may not be necessary to include it.


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
from ._exception import (UnrecognizedFormatError, FileFormatError,
                         BLAST7FormatError, ClustalFormatError,
                         FASTAFormatError, GenBankFormatError, IOSourceError,
                         FASTQFormatError, LSMatFormatError, NewickFormatError,
                         OrdinationFormatError, PhylipFormatError,
                         QSeqFormatError, QUALFormatError,
                         StockholmFormatError, GFF3FormatError,
                         EMBLFormatError)
from .registry import write, read, sniff, create_format, io_registry
from .util import open

__all__ = ['write', 'read', 'sniff', 'open', 'io_registry', 'create_format',

           'FormatIdentificationWarning', 'ArgumentOverrideWarning',
           'UnrecognizedFormatError', 'IOSourceError',

           'FileFormatError',
           'BLAST7FormatError',
           'ClustalFormatError',
           'EMBLFormatError',
           'FASTAFormatError',
           'FASTQFormatError',
           'GenBankFormatError',
           'GFF3FormatError',
           'LSMatFormatError',
           'NewickFormatError',
           'OrdinationFormatError',
           'PhylipFormatError',
           'QSeqFormatError',
           'QUALFormatError',
           'StockholmFormatError']


# Necessary to import each file format module to have them added to the I/O
# registry. We use import_module instead of a typical import to avoid flake8
# unused import errors.
import_module('skbio.io.format.blast6')
import_module('skbio.io.format.blast7')
import_module('skbio.io.format.clustal')
import_module('skbio.io.format.embl')
import_module('skbio.io.format.fasta')
import_module('skbio.io.format.fastq')
import_module('skbio.io.format.lsmat')
import_module('skbio.io.format.newick')
import_module('skbio.io.format.ordination')
import_module('skbio.io.format.phylip')
import_module('skbio.io.format.qseq')
import_module('skbio.io.format.genbank')
import_module('skbio.io.format.gff3')
import_module('skbio.io.format.stockholm')

# This is meant to be a handy indicator to the user that they have done
# something wrong.
import_module('skbio.io.format.emptyfile')

# Now that all of our I/O has loaded, we can add the object oriented methods
# (read and write) to each class which has registered I/O operations.
io_registry.monkey_patch()
