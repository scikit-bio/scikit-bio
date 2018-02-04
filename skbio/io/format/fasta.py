"""
FASTA/QUAL format (:mod:`skbio.io.format.fasta`)
================================================

.. currentmodule:: skbio.io.format.fasta

The FASTA file format (``fasta``) stores biological (i.e., nucleotide or
protein) sequences in a simple plain text format that is both human-readable
and easy to parse. The file format was first introduced and used in the FASTA
software package [1]_. Additional descriptions of the file format can be found
in [2]_ and [3]_.

An example of a FASTA-formatted file containing two DNA sequences::

    >seq1 db-accession-149855
    CGATGTCGATCGATCGATCGATCAG
    >seq2 db-accession-34989
    CATCGATCGATCGATGCATGCATGCATG

The QUAL file format is an additional format related to FASTA. A FASTA file is
sometimes accompanied by a QUAL file, particuarly when the FASTA file contains
sequences generated on a high-throughput sequencing instrument. QUAL files
store a Phred quality score (nonnegative integer) for each base in a sequence
stored in FASTA format (see [4]_ for more details). scikit-bio supports reading
and writing FASTA (and optionally QUAL) file formats.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |generator of :mod:`skbio.sequence.Sequence` objects            |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.alignment.TabularMSA`                              |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.RNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.Protein`                                  |
+------+------+---------------------------------------------------------------+

.. note:: All readers and writers support an optional QUAL file via the
   ``qual`` parameter. If one is provided, quality scores will be read/written
   in addition to FASTA sequence data.

Format Specification
--------------------
The following sections define the FASTA and QUAL file formats in detail.

FASTA Format
^^^^^^^^^^^^
A FASTA file contains one or more biological sequences. The sequences are
stored sequentially, with a *record* for each sequence (also referred to as a
*FASTA record*). Each *record* consists of a single-line *header* (sometimes
referred to as a *defline*, *label*, *description*, or *comment*) followed by
the sequence data, optionally split over multiple lines.

.. note:: Blank or whitespace-only lines are only allowed at the beginning of
   the file, between FASTA records, or at the end of the file. A blank or
   whitespace-only line after the header line, within the sequence (for FASTA
   files), or within quality scores (for QUAL files) will raise an error.

   scikit-bio will ignore leading and trailing whitespace characters on each
   line while reading.

.. note:: scikit-bio does not currently support legacy FASTA format (i.e.,
   headers/comments denoted with a semicolon). The format supported by
   scikit-bio (described below in detail) most closely resembles the
   description given in NCBI's BLAST documentation [3]_. See [2]_ for more
   details on legacy FASTA format. If you would like legacy FASTA format
   support added to scikit-bio, please consider submitting a feature request on
   the
   `scikit-bio issue tracker <https://github.com/biocore/scikit-bio/issues>`_
   (pull requests are also welcome!).

Sequence Header
~~~~~~~~~~~~~~~
Each sequence header consists of a single line beginning with a greater-than
(``>``) symbol. Immediately following this is a sequence identifier (ID) and
description separated by one or more whitespace characters.

.. note:: When reading a FASTA-formatted file, the sequence ID and description
   are stored in the sequence `metadata` attribute, under the `'id'` and
   `'description'` keys, repectively. Both are optional. Each will be
   represented as the empty string (``''``) in `metadata` if it is not present
   in the header.

   When writing a FASTA-formatted file, sequence `metadata` identified by keys
   `'id'` and `'description'` will be converted to strings and written as the
   sequence identifier and description, respectively. Each will be written as
   the empty string if not present in sequence `metadata`.

A sequence ID consists of a single *word*: all characters after the greater-
than symbol and before the first whitespace character (if any) are taken as the
sequence ID. Unique sequence IDs are not strictly enforced by the FASTA format
itself. A single standardized ID format is similarly not enforced by the FASTA
format, though it is often common to use a unique library accession number for
a sequence ID (e.g., NCBI's FASTA defline format [5]_).

If a description is present, it is taken as the remaining characters that
follow the sequence ID and initial whitespace(s). The description is considered
additional information about the sequence (e.g., comments about the source of
the sequence or the molecule that it encodes).

For example, consider the following header::

    >seq1 db-accession-149855

``seq1`` is the sequence ID and ``db-accession-149855`` is the sequence
description.

.. note:: scikit-bio's readers will remove all leading and trailing whitespace
   from the description. If a header line begins with whitespace following the
   ``>``, the ID is assumed to be missing and the remainder of the line is
   taken as the description.

Sequence Data
~~~~~~~~~~~~~
Biological sequence data follows the header, and can be split over multiple
lines. The sequence data (i.e., nucleotides or amino acids) are stored using
the standard IUPAC lexicon (single-letter codes).

.. note:: scikit-bio supports both upper and lower case characters.
   This functionality depends on the type of object the data is
   being read into. For ``Sequence`` objects, sciki-bio doesn't care about the
   case. Other sequence objects do, but all provide the `lowercase` parameter
   to control case functionality. Refer to each class's respective constructor
   documentation for details.

   Both ``-`` and ``.`` are supported as gap characters when reading into
   ``DNA``, ``RNA``, and ``Protein`` sequence objects.

   Validation is performed when reading into scikit-bio sequence objects that
   enforce an alphabet (e.g., ``DNA``, ``RNA``, ``Protein``). If any invalid
   characters are found while reading from the FASTA file, an exception is
   raised.

QUAL Format
^^^^^^^^^^^
A QUAL file contains quality scores for one or more biological sequences stored
in a corresponding FASTA file. QUAL format is very similar to FASTA format: it
stores records sequentially, with each record beginning with a header line
containing a sequence ID and description. The same rules apply to QUAL headers
as FASTA headers (see the above sections for details). scikit-bio processes
FASTA and QUAL headers in exactly the same way.

Instead of storing biological sequence data in each record, a QUAL file stores
a Phred quality score for each base in the corresponding sequence. Quality
scores are represented as nonnegative integers separated by whitespace
(typically a single space or newline), and can span multiple lines.

.. note:: When reading a QUAL-formatted file, quality scores are stored in the
   sequence's `positional_metadata` attribute under the `'quality'` column.

   When writing a QUAL-formatted file, a sequence's `positional_metadata`
   `'quality'` column will be written as the quality scores.

.. note:: When reading FASTA and QUAL files, scikit-bio requires records to be
   in the same order in both files (i.e., each FASTA and QUAL record must have
   the same ID and description after being parsed). In addition to having the
   same order, the number of FASTA records must match the number of QUAL
   records (i.e., missing or additonal records are not allowed). scikit-bio
   also requires that the number of quality scores match the number of bases in
   the corresponding sequence.

   When writing FASTA and QUAL files, scikit-bio will maintain the same
   ordering of records in both files (i.e., using the same ID and description
   in both records) to support future reading.

Format Parameters
-----------------
The following parameters are available to change how FASTA/QUAL files are read
or written in scikit-bio.

QUAL File Parameter (Readers and Writers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``qual`` parameter is available to all FASTA format readers and writers. It
can be any file-like type supported by scikit-bio's I/O registry (e.g., file
handle, file path, etc.). If ``qual`` is provided when reading, quality scores
will be included in each in-memory ``Sequence`` object, in addition
to sequence data stored in the FASTA file. When writing, quality scores will be
written in QUAL format in addition to the sequence data being written in FASTA
format.

Reader-specific Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^
The available reader parameters differ depending on which reader is used.

Generator and TabularMSA Reader Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``constructor`` parameter can be used with the ``Sequence`` generator and
``TabularMSA`` FASTA readers. ``constructor`` specifies the type of in-memory
sequence object to read each sequence into. For example, if you know that the
FASTA file you're reading contains protein sequences, you would pass
``constructor=Protein`` to the reader call.

When reading into a ``Sequence`` generator, ``constructor`` defaults to
``Sequence`` and must be a subclass of ``Sequence`` if supplied.

When reading into a ``TabularMSA``, ``constructor`` is a required format
parameter and must be a subclass of ``GrammaredSequence`` (e.g., ``DNA``,
``RNA``, ``Protein``).

.. note:: The FASTA sniffer will not attempt to guess the ``constructor``
   parameter.

Sequence Reader Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``seq_num`` parameter can be used with the ``Sequence``,
``DNA``, ``RNA``, and ``Protein`` FASTA readers. ``seq_num`` specifies which
sequence to read from the FASTA file (and optional QUAL file), and defaults to
1 (i.e., such that the first sequence is read). For example, to read the 50th
sequence from a FASTA file, you would pass ``seq_num=50`` to the reader call.

Writer-specific Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^
The following parameters are available to all FASTA format writers:

- ``id_whitespace_replacement``: string to replace **each** whitespace
  character in a sequence ID. This parameter is useful for cases where an
  in-memory sequence ID contains whitespace, which would result in an on-disk
  representation that would not be read back into memory as the same ID (since
  IDs in FASTA format cannot contain whitespace). Defaults to ``_``. If
  ``None``, no whitespace replacement is performed and IDs are written as they
  are stored in memory (this has the potential to create an invalid
  FASTA-formatted file; see note below). This parameter also applies to a QUAL
  file if one is provided.

- ``description_newline_replacement``: string to replace **each** newline
  character in a sequence description. Since a FASTA header must be a single
  line, newlines are not allowed in sequence descriptions and must be replaced
  in order to write a valid FASTA file. Defaults to a single space. If
  ``None``, no newline replacement is performed and descriptions are written as
  they are stored in memory (this has the potential to create an invalid
  FASTA-formatted file; see note below). This parameter also applies to a QUAL
  file if one is provided.

- ``max_width``: integer specifying the maximum line width (i.e., number of
  characters) for sequence data and/or quality scores. If a sequence or its
  quality scores are longer than ``max_width``, it will be split across
  multiple lines, each with a maximum width of ``max_width``. Note that there
  are some caveats when splitting quality scores. A single quality score will
  *never* be split across multiple lines, otherwise it would become two
  different quality scores when read again. Thus, splitting only occurs
  *between* quality scores. This makes it possible to have a single long
  quality score written on its own line that exceeds ``max_width``. For
  example, the quality score ``12345`` would not be split across multiple lines
  even if ``max_width=3``. Thus, a 5-character line would be written. Default
  behavior is to not split sequence data or quality scores across multiple
  lines.

- ``lowercase``: String or boolean array. If a string, it is treated as a key
  into the positional metadata of the object. If a boolean array, it
  indicates characters to write in lowercase. Characters in the sequence
  corresponding to `True` values will be written in lowercase. The boolean
  array must be the same length as the sequence.

.. note:: The FASTA format writers will have noticeably better runtime
   performance if ``id_whitespace_replacement`` and/or
   ``description_newline_replacement`` are set to ``None`` so that whitespace
   replacement is not performed during writing. However, this can potentially
   create invalid FASTA files, especially if there are newline characters in
   the IDs or descriptions. For IDs with whitespace, this can also affect how
   the IDs are read into memory in a subsequent read operation. For example, if
   an in-memory sequence ID is ``'seq 1'`` and
   ``id_whitespace_replacement=None``, reading the FASTA file back into memory
   would result in an ID of ``'seq'``, and ``'1'`` would be part of the
   sequence description.

Examples
--------

Reading and Writing FASTA Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Suppose we have the following FASTA file with five equal-length sequences
(example modified from [6]_)::

    >seq1 Turkey
    AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT
    >seq2 Salmo gair
    AAGCCTTGGCAGTGCAGGGTGAGCCGTGG
    CCGGGCACGGTAT
    >seq3 H. Sapiens
    ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA
    >seq4 Chimp
    AAACCCTTGCCG
    TTACGCTTAAAC
    CGAGGCCGGGAC
    ACTCAT
    >seq5 Gorilla
    AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA

.. note:: Original copyright notice for the above example file:

   *(c) Copyright 1986-2008 by The University of Washington. Written by Joseph
   Felsenstein. Permission is granted to copy this document provided that no
   fee is charged for it and that this copyright notice is not removed.*

Note that the sequences are not required to be of equal length in order for the
file to be a valid FASTA file (this depends on the object that you're reading
the file into). Also note that some of the sequences occur on a single line,
while others are split across multiple lines.

Let's define this file in-memory as a ``StringIO``, though this could be a real
file path, file handle, or anything that's supported by scikit-bio's I/O
registry in practice:

>>> fl = [">seq1 Turkey\\n",
...       "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT\\n",
...       ">seq2 Salmo gair\\n",
...       "AAGCCTTGGCAGTGCAGGGTGAGCCGTGG\\n",
...       "CCGGGCACGGTAT\\n",
...       ">seq3 H. Sapiens\\n",
...       "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA\\n",
...       ">seq4 Chimp\\n",
...       "AAACCCTTGCCG\\n",
...       "TTACGCTTAAAC\\n",
...       "CGAGGCCGGGAC\\n",
...       "ACTCAT\\n",
...       ">seq5 Gorilla\\n",
...       "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA\\n"]

Since these sequences are of equal length (presumably because they've been
aligned), let's read the FASTA file into a ``TabularMSA`` object:

>>> from skbio import TabularMSA, DNA
>>> msa = TabularMSA.read(fl, constructor=DNA)
>>> msa
TabularMSA[DNA]
------------------------------------------
Stats:
    sequence count: 5
    position count: 42
------------------------------------------
AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT
AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT
ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA
AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT
AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA

Note that we didn't specify a file format in the ``read`` call. The FASTA
sniffer detected the correct file format for us!

To write the ``TabularMSA`` in FASTA format:

>>> from io import StringIO
>>> with StringIO() as fh:
...     print(msa.write(fh).getvalue())
>seq1 Turkey
AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT
>seq2 Salmo gair
AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT
>seq3 H. Sapiens
ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA
>seq4 Chimp
AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT
>seq5 Gorilla
AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA
<BLANKLINE>

``TabularMSA`` loads all of the sequences from the FASTA file into memory at
once. If the FASTA file is large (which is often the case), this may be
infeasible if you don't have enough memory. To work around this issue, you can
stream the sequences using scikit-bio's generator-based FASTA reader and
writer. The generator-based reader yields ``Sequence`` objects (or subclasses
if ``constructor`` is supplied) one at a time, instead of loading all sequences
into memory. For example, let's use the generator-based reader to process a
single sequence at a time in a ``for`` loop:

>>> import skbio.io
>>> for seq in skbio.io.read(fl, format='fasta'):
...     seq
...     print('')
Sequence
------------------------------------------------
Metadata:
    'description': 'Turkey'
    'id': 'seq1'
Stats:
    length: 42
------------------------------------------------
0 AAGCTNGGGC ATTTCAGGGT GAGCCCGGGC AATACAGGGT AT
<BLANKLINE>
Sequence
------------------------------------------------
Metadata:
    'description': 'Salmo gair'
    'id': 'seq2'
Stats:
    length: 42
------------------------------------------------
0 AAGCCTTGGC AGTGCAGGGT GAGCCGTGGC CGGGCACGGT AT
<BLANKLINE>
Sequence
------------------------------------------------
Metadata:
    'description': 'H. Sapiens'
    'id': 'seq3'
Stats:
    length: 42
------------------------------------------------
0 ACCGGTTGGC CGTTCAGGGT ACAGGTTGGC CGTTCAGGGT AA
<BLANKLINE>
Sequence
------------------------------------------------
Metadata:
    'description': 'Chimp'
    'id': 'seq4'
Stats:
    length: 42
------------------------------------------------
0 AAACCCTTGC CGTTACGCTT AAACCGAGGC CGGGACACTC AT
<BLANKLINE>
Sequence
------------------------------------------------
Metadata:
    'description': 'Gorilla'
    'id': 'seq5'
Stats:
    length: 42
------------------------------------------------
0 AAACCCTTGC CGGTACGCTT AAACCATTGC CGGTACGCTT AA
<BLANKLINE>

A single sequence can also be read into a ``Sequence`` (or subclass):

>>> from skbio import Sequence
>>> seq = Sequence.read(fl)
>>> seq
Sequence
------------------------------------------------
Metadata:
    'description': 'Turkey'
    'id': 'seq1'
Stats:
    length: 42
------------------------------------------------
0 AAGCTNGGGC ATTTCAGGGT GAGCCCGGGC AATACAGGGT AT

By default, the first sequence in the FASTA file is read. This can be
controlled with ``seq_num``. For example, to read the fifth sequence:

>>> seq = Sequence.read(fl, seq_num=5)
>>> seq
Sequence
------------------------------------------------
Metadata:
    'description': 'Gorilla'
    'id': 'seq5'
Stats:
    length: 42
------------------------------------------------
0 AAACCCTTGC CGGTACGCTT AAACCATTGC CGGTACGCTT AA

We can use the same API to read the fifth sequence into a ``DNA`` sequence:

>>> dna_seq = DNA.read(fl, seq_num=5)
>>> dna_seq
DNA
------------------------------------------------
Metadata:
    'description': 'Gorilla'
    'id': 'seq5'
Stats:
    length: 42
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 50.00%
------------------------------------------------
0 AAACCCTTGC CGGTACGCTT AAACCATTGC CGGTACGCTT AA

Individual sequence objects can also be written in FASTA format:

>>> with StringIO() as fh:
...     print(dna_seq.write(fh).getvalue())
>seq5 Gorilla
AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA
<BLANKLINE>

Reading and Writing FASTA/QUAL Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In addition to reading and writing standalone FASTA files, scikit-bio supports
reading and writing FASTA and QUAL files together. Suppose we have the
following FASTA file::

    >seq1 db-accession-149855
    CGATGTC
    >seq2 db-accession-34989
    CATCGTC

Also suppose we have the following QUAL file::

    >seq1 db-accession-149855
    40 39 39 4
    50 1 100
    >seq2 db-accession-34989
    3 3 10 42 80 80 79

>>> fasta_fl = [
...     ">seq1 db-accession-149855\\n",
...     "CGATGTC\\n",
...     ">seq2 db-accession-34989\\n",
...     "CATCGTC\\n"]
>>> qual_fl = [
...     ">seq1 db-accession-149855\\n",
...     "40 39 39 4\\n",
...     "50 1 100\\n",
...     ">seq2 db-accession-34989\\n",
...     "3 3 10 42 80 80 79\\n"]

To read in a single ``Sequence`` at a time, we can use the
generator-based reader as we did above, providing both FASTA and QUAL files:

>>> for seq in skbio.io.read(fasta_fl, qual=qual_fl, format='fasta'):
...     seq
...     print('')
Sequence
----------------------------------------
Metadata:
    'description': 'db-accession-149855'
    'id': 'seq1'
Positional metadata:
    'quality': <dtype: uint8>
Stats:
    length: 7
----------------------------------------
0 CGATGTC
<BLANKLINE>
Sequence
---------------------------------------
Metadata:
    'description': 'db-accession-34989'
    'id': 'seq2'
Positional metadata:
    'quality': <dtype: uint8>
Stats:
    length: 7
---------------------------------------
0 CATCGTC
<BLANKLINE>

Note that the sequence objects have quality scores stored as positional
metadata since we provided a QUAL file. The other FASTA readers operate in a
similar manner.

Now let's load the sequences and their quality scores into a ``TabularMSA``:

>>> msa = TabularMSA.read(fasta_fl, qual=qual_fl, constructor=DNA)
>>> msa
TabularMSA[DNA]
---------------------
Stats:
    sequence count: 2
    position count: 7
---------------------
CGATGTC
CATCGTC

To write the sequence data and quality scores in the ``TabularMSA`` to FASTA
and QUAL files, respectively:

>>> new_fasta_fh = StringIO()
>>> new_qual_fh = StringIO()
>>> _ = msa.write(new_fasta_fh, qual=new_qual_fh)
>>> print(new_fasta_fh.getvalue())
>seq1 db-accession-149855
CGATGTC
>seq2 db-accession-34989
CATCGTC
<BLANKLINE>
>>> print(new_qual_fh.getvalue())
>seq1 db-accession-149855
40 39 39 4 50 1 100
>seq2 db-accession-34989
3 3 10 42 80 80 79
<BLANKLINE>
>>> new_fasta_fh.close()
>>> new_qual_fh.close()

References
----------
.. [1] Lipman, DJ; Pearson, WR (1985). "Rapid and sensitive protein similarity
   searches". Science 227 (4693): 1435-41.
.. [2] http://en.wikipedia.org/wiki/FASTA_format
.. [3] http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
.. [4] https://www.broadinstitute.org/crd/wiki/index.php/Qual
.. [5] Madden T. The BLAST Sequence Analysis Tool. 2002 Oct 9
   [Updated 2003 Aug 13]. In: McEntyre J, Ostell J, editors. The NCBI Handbook
   [Internet]. Bethesda (MD): National Center for Biotechnology Information
   (US); 2002-. Chapter 16. Available from:
   http://www.ncbi.nlm.nih.gov/books/NBK21097/
.. [6] http://evolution.genetics.washington.edu/phylip/doc/sequence.html

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
import textwrap

import numpy as np

from skbio.io import create_format, FASTAFormatError, QUALFormatError
from skbio.io.registry import FileSentinel
from skbio.io.format._base import (_get_nth_sequence,
                                   _parse_fasta_like_header,
                                   _format_fasta_like_records, _line_generator,
                                   _too_many_blanks)
from skbio.util._misc import chunk_str
from skbio.alignment import TabularMSA
from skbio.sequence import Sequence, DNA, RNA, Protein


fasta = create_format('fasta')


@fasta.sniffer()
def _fasta_sniffer(fh):
    # Strategy:
    #   Ignore up to 5 blank/whitespace-only lines at the beginning of the
    #   file. Read up to 10 records. If at least one record is read (i.e.
    #   the file isn't empty) and no errors are thrown during reading, assume
    #   the file is in FASTA format. If a record appears to be QUAL, do *not*
    #   identify the file as FASTA since we don't want to sniff QUAL files as
    #   FASTA (technically they can be read as FASTA since the sequences may
    #   not be validated but it probably isn't what the user wanted). Also, if
    #   we add QUAL as its own file format in the future, we wouldn't want the
    #   FASTA and QUAL sniffers to both positively identify a QUAL file.
    if _too_many_blanks(fh, 5):
        return False, {}

    num_records = 10
    empty = True
    try:
        parser = _parse_fasta_raw(fh, _sniffer_data_parser, FASTAFormatError)
        for _ in zip(range(num_records), parser):
            empty = False
    except FASTAFormatError:
        return False, {}

    if empty:
        return False, {}
    else:
        return True, {}


def _sniffer_data_parser(chunks):
    data = _parse_sequence_data(chunks)
    try:
        _parse_quality_scores(chunks)
    except QUALFormatError:
        return data
    else:
        # used for flow control within sniffer, user should never see this
        # message
        raise FASTAFormatError('Data appear to be quality scores.')


@fasta.reader(None)
def _fasta_to_generator(fh, qual=FileSentinel, constructor=Sequence, **kwargs):
    if qual is None:
        for seq, id_, desc in _parse_fasta_raw(fh, _parse_sequence_data,
                                               FASTAFormatError):
            yield constructor(seq, metadata={'id': id_, 'description': desc},
                              **kwargs)
    else:
        fasta_gen = _parse_fasta_raw(fh, _parse_sequence_data,
                                     FASTAFormatError)
        qual_gen = _parse_fasta_raw(qual, _parse_quality_scores,
                                    QUALFormatError)

        for fasta_rec, qual_rec in itertools.zip_longest(fasta_gen, qual_gen,
                                                         fillvalue=None):
            if fasta_rec is None:
                raise FASTAFormatError(
                    "QUAL file has more records than FASTA file.")
            if qual_rec is None:
                raise FASTAFormatError(
                    "FASTA file has more records than QUAL file.")

            fasta_seq, fasta_id, fasta_desc = fasta_rec
            qual_scores, qual_id, qual_desc = qual_rec

            if fasta_id != qual_id:
                raise FASTAFormatError(
                    "IDs do not match between FASTA and QUAL records: %r != %r"
                    % (str(fasta_id), str(qual_id)))
            if fasta_desc != qual_desc:
                raise FASTAFormatError(
                    "Descriptions do not match between FASTA and QUAL "
                    "records: %r != %r" % (str(fasta_desc), str(qual_desc)))

            # sequence and quality scores lengths are checked in constructor
            yield constructor(
                fasta_seq,
                metadata={'id': fasta_id, 'description': fasta_desc},
                positional_metadata={'quality': qual_scores}, **kwargs)


@fasta.reader(Sequence)
def _fasta_to_sequence(fh, qual=FileSentinel, seq_num=1, **kwargs):
    return _get_nth_sequence(
        _fasta_to_generator(fh, qual=qual, constructor=Sequence, **kwargs),
        seq_num)


@fasta.reader(DNA)
def _fasta_to_dna(fh, qual=FileSentinel, seq_num=1, **kwargs):
    return _get_nth_sequence(
        _fasta_to_generator(fh, qual=qual,
                            constructor=DNA, **kwargs),
        seq_num)


@fasta.reader(RNA)
def _fasta_to_rna(fh, qual=FileSentinel, seq_num=1, **kwargs):
    return _get_nth_sequence(
        _fasta_to_generator(fh, qual=qual,
                            constructor=RNA, **kwargs),
        seq_num)


@fasta.reader(Protein)
def _fasta_to_protein(fh, qual=FileSentinel, seq_num=1, **kwargs):
    return _get_nth_sequence(
        _fasta_to_generator(fh, qual=qual,
                            constructor=Protein, **kwargs),
        seq_num)


@fasta.reader(TabularMSA)
def _fasta_to_tabular_msa(fh, qual=FileSentinel, constructor=None, **kwargs):
    if constructor is None:
        raise ValueError("Must provide `constructor`.")

    return TabularMSA(
        _fasta_to_generator(fh, qual=qual, constructor=constructor, **kwargs))


@fasta.writer(None)
def _generator_to_fasta(obj, fh, qual=FileSentinel,
                        id_whitespace_replacement='_',
                        description_newline_replacement=' ', max_width=None,
                        lowercase=None):
    if max_width is not None:
        if max_width < 1:
            raise ValueError(
                "Maximum line width must be greater than zero (max_width=%d)."
                % max_width)
        if qual is not None:
            # define text wrapper for splitting quality scores here for
            # efficiency. textwrap docs recommend reusing a TextWrapper
            # instance when it is used many times. configure text wrapper to
            # never break "words" (i.e., integer quality scores) across lines
            qual_wrapper = textwrap.TextWrapper(
                width=max_width, break_long_words=False,
                break_on_hyphens=False)

    formatted_records = _format_fasta_like_records(
        obj, id_whitespace_replacement, description_newline_replacement,
        qual is not None, lowercase)
    for header, seq_str, qual_scores in formatted_records:
        if max_width is not None:
            seq_str = chunk_str(seq_str, max_width, '\n')

        fh.write('>%s\n%s\n' % (header, seq_str))

        if qual is not None:
            qual_str = ' '.join(np.asarray(qual_scores, dtype=np.str))
            if max_width is not None:
                qual_str = qual_wrapper.fill(qual_str)
            qual.write('>%s\n%s\n' % (header, qual_str))


@fasta.writer(Sequence)
def _sequence_to_fasta(obj, fh, qual=FileSentinel,
                       id_whitespace_replacement='_',
                       description_newline_replacement=' ', max_width=None,
                       lowercase=None):
    _sequences_to_fasta([obj], fh, qual, id_whitespace_replacement,
                        description_newline_replacement, max_width, lowercase)


@fasta.writer(DNA)
def _dna_to_fasta(obj, fh, qual=FileSentinel, id_whitespace_replacement='_',
                  description_newline_replacement=' ', max_width=None,
                  lowercase=None):
    _sequences_to_fasta([obj], fh, qual, id_whitespace_replacement,
                        description_newline_replacement, max_width, lowercase)


@fasta.writer(RNA)
def _rna_to_fasta(obj, fh, qual=FileSentinel, id_whitespace_replacement='_',
                  description_newline_replacement=' ', max_width=None,
                  lowercase=None):
    _sequences_to_fasta([obj], fh, qual, id_whitespace_replacement,
                        description_newline_replacement, max_width, lowercase)


@fasta.writer(Protein)
def _protein_to_fasta(obj, fh, qual=FileSentinel,
                      id_whitespace_replacement='_',
                      description_newline_replacement=' ', max_width=None,
                      lowercase=None):
    _sequences_to_fasta([obj], fh, qual, id_whitespace_replacement,
                        description_newline_replacement, max_width, lowercase)


@fasta.writer(TabularMSA)
def _tabular_msa_to_fasta(obj, fh, qual=FileSentinel,
                          id_whitespace_replacement='_',
                          description_newline_replacement=' ', max_width=None,
                          lowercase=None):
    _sequences_to_fasta(obj, fh, qual, id_whitespace_replacement,
                        description_newline_replacement, max_width, lowercase)


def _parse_fasta_raw(fh, data_parser, error_type):
    """Raw parser for FASTA or QUAL files.

    Returns raw values (seq/qual, id, description). It is the responsibility of
    the caller to construct the correct in-memory object to hold the data.

    """
    # Skip any blank or whitespace-only lines at beginning of file
    try:
        seq_header = next(_line_generator(fh, skip_blanks=True))
    except StopIteration:
        return

    # header check inlined here and below for performance
    if seq_header.startswith('>'):
        id_, desc = _parse_fasta_like_header(seq_header)
    else:
        raise error_type(
            "Found non-header line when attempting to read the 1st record:"
            "\n%s" % seq_header)

    data_chunks = []
    prev = seq_header
    for line in _line_generator(fh, skip_blanks=False):
        if line.startswith('>'):
            # new header, so yield current record and reset state
            yield data_parser(data_chunks), id_, desc
            data_chunks = []
            id_, desc = _parse_fasta_like_header(line)
        else:
            if line:
                # ensure no blank lines within a single record
                if not prev:
                    raise error_type(
                        "Found blank or whitespace-only line within record.")
                data_chunks.append(line)
        prev = line
    # yield last record in file
    yield data_parser(data_chunks), id_, desc


def _parse_sequence_data(chunks):
    if not chunks:
        raise FASTAFormatError("Found header without sequence data.")
    return ''.join(chunks)


def _parse_quality_scores(chunks):
    if not chunks:
        raise QUALFormatError("Found header without quality scores.")

    qual_str = ' '.join(chunks)
    try:
        quality = np.asarray(qual_str.split(), dtype=int)
    except ValueError:
        raise QUALFormatError(
            "Could not convert quality scores to integers:\n%s"
            % str(qual_str))

    if (quality < 0).any():
        raise QUALFormatError(
            "Encountered negative quality score(s). Quality scores must be "
            "greater than or equal to zero.")
    if (quality > 255).any():
        raise QUALFormatError(
            "Encountered quality score(s) greater than 255. scikit-bio only "
            "supports quality scores in the range 0-255 (inclusive) when "
            "reading QUAL files.")
    return quality.astype(np.uint8, casting='unsafe', copy=False)


def _sequences_to_fasta(obj, fh, qual, id_whitespace_replacement,
                        description_newline_replacement, max_width,
                        lowercase=None):
    def seq_gen():
        yield from obj

    _generator_to_fasta(
        seq_gen(), fh, qual=qual,
        id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        max_width=max_width, lowercase=lowercase)
