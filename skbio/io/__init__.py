r"""Input and Output (:mod:`skbio.io`)
==================================

.. currentmodule:: skbio.io

This module provides input/output (I/O) functionality for scikit-bio.

In bioinformatics there are many different file formats, and in scikit-bio there are
many different classes which can read and write these formats. The many-to-many
nature of the relationships between scikit-bio objects and file formats inspired
the creation of the scikit-bio ``io`` module, which manages these relationships
transparently.

For general guidance on reading and writing files and working with scikit-bio objects,
see the :ref:`tutorial` section and the
`Reading and writing files <https://github.com/scikit-bio/scikit-bio-cookbook/blob/
master/Reading%20and%20writing%20files.ipynb>`_
notebook. For guidance on a specific format or scikit-bio object,
see the documentation for that format or object.

See the
`IORegistry docs <https://scikit.bio/docs/latest/generated/skbio.io.registry.html
#creating-a-new-format-for-scikit-bio>`_
for guidance on creating custom formats and registering custom readers, writers, and
sniffers.

Supported file formats
----------------------

scikit-bio provides parsers for the following file formats. For details on what objects
are supported by each format, see the associated documentation.

.. currentmodule:: skbio.io.format

.. autosummary::
   :toctree: generated/

   binary_dm
   biom
   blast6
   blast7
   clustal
   embl
   embed
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
   taxdump
   sample_metadata


Read/write files
----------------

.. rubric:: Generic I/O functions

.. currentmodule:: skbio.io.registry

.. autosummary::
   :toctree: generated/

   write
   read
   sniff

.. rubric:: Additional I/O utilities

.. currentmodule:: skbio.io

.. autosummary::
   :toctree: generated/

   util


Develop custom formats
----------------------

.. rubric:: Developer documentation on extending I/O

.. autosummary::
   :toctree: generated/

   registry


Exceptions and warnings
-----------------------

.. currentmodule:: skbio.io

.. rubric:: General exceptions and warnings

.. autosummary::

   FormatIdentificationWarning
   ArgumentOverrideWarning
   UnrecognizedFormatError
   IOSourceError
   FileFormatError

.. rubric:: Format-specific exceptions and warnings

.. autosummary::

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


.. _tutorial:

Tutorial
--------

Reading and writing files (I/O) can be a complicated task:

* A file format can sometimes be read into more than one in-memory representation
  (i.e., object). For example, a FASTA file can be read into a
  :class:`~skbio.alignment.TabularMSA` or :class:`~skbio.sequence.DNA` depending on
  what operations you'd like to perform on your data.
* A single object might be writeable to more than one file format. For example, an
  :class:`~skbio.alignment.TabularMSA` object could be written to FASTA, FASTQ,
  CLUSTAL, or PHYLIP formats, just to name a few.
* You might not know the exact file format of your file, but you want to read
  it into an appropriate object.
* You might want to read multiple files into a single object, or write an
  object to multiple files.
* Instead of reading a file into an object, you might want to stream the file
  using a generator (e.g., if the file cannot be fully loaded into memory).

To address these issues (and others), scikit-bio provides a simple, powerful
interface for dealing with I/O. We accomplish this by using a single I/O
registry defined in :class:`~skbio.io.registry.IORegistry`.

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

Here, ``file`` can be a path to a file, a file handle, or any of the other
objects with read support listed in the :func:`skbio.io.util.open` documentation.

The second way to read files is to use the object-oriented interface, which is
automatically constructed from the procedural interface:

.. code-block:: python

   my_obj = SomeSkbioClass.read(file, format='someformat')

.. note::
   A very common use case in bioinformatics is to read multi-line FASTA and
   FASTQ files. For examples on how to achieve this with scikit-bio, please see the
   `FASTA documentation <https://scikit.bio/docs/dev/generated/skbio.io.format.fasta.html
   #examples>`_
   or the
   `FASTQ documentation <https://scikit.bio/docs/dev/generated/skbio.io.format.fastq.html
   #examples>`_.

As an example, let's read a :mod:`~skbio.io.format.newick` file into a
:class:`~skbio.tree.TreeNode` object using both interfaces. Here we will use Python's
built-in :class:`~io.StringIO` class to mimick an open file:

>>> from skbio import read as sk_read
>>> from skbio import TreeNode
>>> from io import StringIO
>>> open_filehandle = StringIO('(a, b);')
>>> tree = sk_read(open_filehandle, format='newick', into=TreeNode)
>>> tree
<TreeNode, name: unnamed, internal node count: 0, tips count: 2>

Or, using the object-oriented interface:

>>> open_filehandle = StringIO('(a, b);')
>>> tree = TreeNode.read(open_filehandle, format='newick')
>>> tree
<TreeNode, name: unnamed, internal node count: 0, tips count: 2>

In the case of :func:`skbio.io.registry.read` if ``into`` is not provided, then a
generator will be returned. What the generator yields will depend on what
format is being read.

When ``into`` is provided, format may be omitted and the registry will use its
knowledge of the available formats for the requested class to infer (sniff) the
correct format. This format inference is also available in the object-oriented
interface, meaning that ``format`` may be omitted there as well.

As an example:

>>> open_filehandle = StringIO('(a, b);')
>>> tree = TreeNode.read(open_filehandle)
>>> tree
<TreeNode, name: unnamed, internal node count: 0, tips count: 2>

We call format inference **sniffing**, much like the :class:`csv.Sniffer`
class of Python's standard library. The goal of a ``sniffer`` is two-fold: to
identify if a file is a specific format, and if it is, to provide ``**kwargs``
which can be used to better parse the file.

.. note::
   There is a built-in ``sniffer`` which results in a useful error message
   if an empty file is provided as input and the format was omitted. See the
   `sniff documentation <https://scikit.bio/docs/dev/generated/skbio.io.registry.sniff.
   html>`_ for more information.

Writing files from scikit-bio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Just as when reading files, there are two ways to write files.

Procedural Interface:

.. code-block:: python

   skbio.io.write(my_obj, format='someformat', into=file)

Object-oriented Interface:

.. code-block:: python

   my_obj.write(file, format='someformat')

In the procedural interface, ``format`` is required. Without it, scikit-bio does
not know how you want to serialize an object. Object-oriented interfaces define a
default ``format``, so it may not be necessary to include it.

For more information on writing to a specific file format, please see that format's
documentation page.

Streaming files with read and write
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you are working with particularly large files, streaming them might be preferable.
For instance, if your file is larger than your available memory, you won't be able
to read the entire file into memory at once. One way to get around this is to use
streaming. Scikit-bio's ``io`` module offers the ability to contruct a streaming
interface from the ``read`` and ``write`` functions.

``skbio.io.read`` returns a generator, which can then be passed to ``skbio.io.write``
to write only one chunk from the generator at a time.

.. code-block:: python

   seq_gen = skbio.io.read(big_file, format='someformat')
   skbio.io.write(seq_gen, into=write_file, format='someformat')


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from importlib import import_module

from ._warning import FormatIdentificationWarning, ArgumentOverrideWarning
from ._exception import (
    UnrecognizedFormatError,
    FileFormatError,
    BLAST7FormatError,
    ClustalFormatError,
    FASTAFormatError,
    GenBankFormatError,
    IOSourceError,
    FASTQFormatError,
    LSMatFormatError,
    NewickFormatError,
    OrdinationFormatError,
    PhylipFormatError,
    QSeqFormatError,
    QUALFormatError,
    StockholmFormatError,
    GFF3FormatError,
    EMBLFormatError,
    BIOMFormatError,
    EmbedFormatError,
)
from .registry import write, read, sniff, create_format, io_registry
from .util import open

__all__ = [
    "write",
    "read",
    "sniff",
    "open",
    "io_registry",
    "create_format",
    "FormatIdentificationWarning",
    "ArgumentOverrideWarning",
    "UnrecognizedFormatError",
    "IOSourceError",
    "FileFormatError",
    "BLAST7FormatError",
    "ClustalFormatError",
    "EMBLFormatError",
    "FASTAFormatError",
    "FASTQFormatError",
    "GenBankFormatError",
    "GFF3FormatError",
    "LSMatFormatError",
    "NewickFormatError",
    "OrdinationFormatError",
    "PhylipFormatError",
    "QSeqFormatError",
    "QUALFormatError",
    "StockholmFormatError",
    "BIOMFormatError",
    "EmbedFormatError",
]


# Necessary to import each file format module to have them added to the I/O
# registry. We use import_module instead of a typical import to avoid flake8
# unused import errors.
import_module("skbio.io.format.blast6")
import_module("skbio.io.format.blast7")
import_module("skbio.io.format.clustal")
import_module("skbio.io.format.embl")
import_module("skbio.io.format.fasta")
import_module("skbio.io.format.fastq")
import_module("skbio.io.format.lsmat")
import_module("skbio.io.format.newick")
import_module("skbio.io.format.ordination")
import_module("skbio.io.format.phylip")
import_module("skbio.io.format.qseq")
import_module("skbio.io.format.genbank")
import_module("skbio.io.format.gff3")
import_module("skbio.io.format.stockholm")
import_module("skbio.io.format.binary_dm")
import_module("skbio.io.format.taxdump")
import_module("skbio.io.format.sample_metadata")
import_module("skbio.io.format.biom")
import_module("skbio.io.format.embed")

# This is meant to be a handy indicator to the user that they have done
# something wrong.
import_module("skbio.io.format.emptyfile")
