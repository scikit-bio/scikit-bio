r"""
Parse biological sequences (:mod:`skbio.parse.sequences`)
=========================================================

.. currentmodule:: skbio.parse.sequences

This module provides functions for parsing sequence files in a variety of
different formats. Two interfaces are provided for parsing sequence files:
sequence iterators (high-level, recommended interface) and parsing functions
(lower-level interface).

Sequence iterator interface
---------------------------
The sequence iterator interface is the recommended way to parse sequence files.
The ``load`` function provides a standard, high-level interface to iterate over
sequence files regardless of file type or whether they are compressed. The
method accepts single or multiple file paths and employs the correct file
handlers, iterator objects, and parsers for the user.

The benefit of the sequence iterator interface is that the type of the file and
any file format details are abstracted away from the user. In this manner, the
user does not need to worry about whether they're operating on FASTA or FASTQ
files or any differences in the returns from their respective parsers.

Classes
^^^^^^^

.. autosummary::
   :toctree: generated/

   SequenceIterator
   FastaIterator
   FastqIterator
   QseqIterator

Functions
^^^^^^^^^

.. autosummary::
   :toctree: generated/

    load

Examples
^^^^^^^^
For the first set of sequence iterator examples, we're going to use the
``load`` function. The ``load`` function is intended to operate on file paths,
so let's create two files for it to use. The first one will be a regular FASTA
file, and the second will be a gzip'd FASTQ file:

>>> import os
>>> import gzip
>>> out = open('test_seqs.fna', 'w')
>>> out.write(">s1\nATGC\n>s2\nATGGC\n")
>>> out.close()
>>> outgz = gzip.open('test_seqs.fq.gz', 'w')
>>> _ = outgz.write("@s3\nAATTGG\n+\nghghgh\n@s4\nAAA\n+\nfgh\n")
>>> outgz.close()

Now let's see what ``load`` can do:

>>> it = load(['test_seqs.fna', 'test_seqs.fq.gz'], phred_offset=64)
>>> for rec in it:
...     print rec['SequenceID']
...     print rec['Sequence']
...     print rec['Qual']
s1
ATGC
None
s2
ATGGC
None
s3
AATTGG
[39 40 39 40 39 40]
s4
AAA
[38 39 40]

To be polite, let's remove the files we just created:

>>> os.remove('test_seqs.fna')
>>> os.remove('test_seqs.fq.gz')

In the following examples, we'll see how to use the sequence iterators directly
instead of using ``load``.

>>> from StringIO import StringIO
>>> from skbio.parse.sequences import FastaIterator, FastqIterator

In this first example, we're going to construct a FASTA iterator that is also
paired with quality scores (e.g., as in 454 fasta/qual files).

>>> seqs = StringIO(">seq1\n"
...                 "ATGC\n"
...                 ">seq2\n"
...                 "TTGGCC\n")
>>> qual = StringIO(">seq1\n"
...                 "10 20 30 40\n"
...                 ">seq2\n"
...                 "1 2 3 4 5 6\n")
>>> it = FastaIterator(seq=[seqs], qual=[qual])
>>> for record in it:
...     print record['Sequence']
...     print record['Qual']
ATGC
[10 20 30 40]
TTGGCC
[1 2 3 4 5 6]

In the next example, we're going to iterate over multiple FASTQ files at once.

>>> seqs1 = StringIO("@seq1\n"
...                  "ATGC\n"
...                  "+\n"
...                  "hhhh\n")
>>> seqs2 = StringIO("@seq2\n"
...                 "AATTGGCC\n"
...                 ">seq2\n"
...                 "abcdefgh\n")
>>> it = FastqIterator(seq=[seqs1, seqs2], phred_offset=64)
>>> for record in it:
...     print record['Sequence']
...     print record['Qual']
ATGC
[40 40 40 40]
AATTGGCC
[33 34 35 36 37 38 39 40]

Finally, we can apply arbitrary transforms to the sequences during iteration.

>>> seqs1 = StringIO("@seq1\n"
...                  "ATGC\n"
...                  "+\n"
...                  "hhhh\n")
>>> seqs2 = StringIO("@seq2\n"
...                 "AATTGGCC\n"
...                 ">seq2\n"
...                 "abcdefgh\n")
>>> def rev_f(st):
...     st['Sequence'] = st['Sequence'][::-1]
...     st['Qual'] = st['Qual'][::-1] if st['Qual'] is not None else None
>>> it = FastqIterator(seq=[seqs1, seqs2], transform=rev_f, phred_offset=64)
>>> for record in it:
...     print record['Sequence']
...     print record['Qual']
CGTA
[40 40 40 40]
CCGGTTAA
[40 39 38 37 36 35 34 33]

Low-level parsing functions
---------------------------
Lower-level parsing functions are also made available in addition to the
sequence iterator interface. These functions can be used to directly parse a
single sequence file. They accept file paths, file handles, or file-like
objects.

Functions
^^^^^^^^^

.. autosummary::
   :toctree: generated/

    parse_fasta
    parse_fastq
    parse_qual
    parse_qseq

Exceptions
----------

.. autosummary::
   :toctree: generated/

    FastqParseError

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .fasta import parse_fasta, parse_qual
from .fastq import parse_fastq
from .qseq import parse_qseq
from .iterator import (FastaIterator, FastqIterator, QseqIterator,
                       SequenceIterator)
from .factory import load
from ._exception import FastqParseError, QseqParseError

__all__ = ['parse_fasta', 'parse_fastq', 'parse_qual',
           'parse_qseq', 'FastqIterator', 'FastaIterator', 'QseqIterator',
           'SequenceIterator', 'load', 'FastqParseError', 'QseqParseError']

from numpy.testing import Tester
test = Tester().test
