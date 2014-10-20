"""
FASTA format (:mod:`skbio.io.fasta`)
====================================

.. currentmodule:: skbio.io.fasta

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

Format Support
--------------
**Has Sniffer: Yes**

+----------+----------+------------------------------------------------------+
|**Reader**|**Writer**|                   **Object Class**                   |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |generator of :mod:`skbio.sequence.BiologicalSequence` |
|          |          |objects                                               |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.alignment.SequenceCollection`             |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.alignment.Alignment`                      |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.sequence.BiologicalSequence`              |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.sequence.NucleotideSequence`              |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.sequence.DNASequence`                     |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.sequence.RNASequence`                     |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.sequence.ProteinSequence`                 |
+----------+----------+------------------------------------------------------+

Format Specification
--------------------
A FASTA file contains one or more biological sequences. The sequences are
stored sequentially, with a *record* for each sequence (also referred to as a
*FASTA record*). Each *record* consists of a single-line *header* (sometimes
referred to as a *defline*, *label*, *description*, or *comment*) followed by
the sequence data, optionally split over multiple lines. Blank or
whitespace-only lines are not allowed anywhere in the FASTA file.

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
^^^^^^^^^^^^^^^
Each sequence header consists of a single line beginning with a greater-than
(``>``) symbol. Immediately following this is a sequence identifier (ID) and
description separated by one or more whitespace characters. Both sequence ID
and description are optional and are represented as the empty string (``''``)
in scikit-bio's objects if they are not present in the header.

A sequence ID consists of a single word; all characters after the greater-than
symbol and before the first whitespace character (if any) are taken as the
sequence ID. Unique sequence IDs are not strictly enforced by the FASTA format
itself. A single standardized ID format is similarly not enforced by FASTA
format, though it is often common to use a unique library accession number for
a sequence ID (e.g., NCBI's FASTA defline format [4]_).

.. note:: scikit-bio will enforce sequence ID uniqueness depending on the type
   of object that the FASTA file is read into. For example, reading a FASTA
   file as a generator of ``BiologicalSequence`` objects will not enforce
   unique IDs since it simply yields each sequence it finds in the FASTA file.
   However, if the FASTA file is read into a ``SequenceCollection`` object, ID
   uniqueness will be enforced because that is a property of a
   ``SequenceCollection``.

If a description is present, it is taken as the remaining characters that
follow the sequence ID and initial whitespace. The description is considered
additional information about the sequence (e.g., comments about the source of
the sequence or the molecule that it encodes).

For example, consider the following header::

    >seq1 db-accession-149855

``seq1`` is the sequence ID and ``db-accession-149855`` is the sequence
description.

.. note:: scikit-bio's readers will remove all leading and trailing whitespace
   from the description. If a header line begins with whitespace following
   ``>``, the ID is assumed to be missing and the remainder of the line is
   taken as the description.

Sequence Data
^^^^^^^^^^^^^
Biological sequence data follows the header, and can be split over multiple
lines. The sequence data (i.e., nucleotides or amino acids) are stored using
the standard IUPAC lexicon (single-letter codes).

.. note:: scikit-bio supports both upper and lower case characters. Both ``-``
   and ``.`` are supported as gap characters. See :mod:`skbio.sequence` for
   more details on how scikit-bio interprets sequence data in its in-memory
   objects.

   scikit-bio will remove leading and trailing whitespace from each line of
   sequence data before joining the sequence chunks into a single sequence.
   Whitespace characters are **not** removed from the middle of the sequence
   chunks. Likewise, other invalid IUPAC characters are **not** removed from
   the sequence data as it is read. Thus, it is possible to create an invalid
   in-memory sequence object (see warning below).

.. warning:: In an effort to maintain reasonable performance while reading
   FASTA files (which can be quite large), validation of sequence data is
   **not** performed during reading. It is the responsibility of the user to
   validate their in-memory representation of the data if desired (e.g., by
   calling ``is_valid`` on the returned object). Thus, it is possible to read
   invalid characters into objects (e.g. whitespace occurring in the middle of
   a sequence, or invalid IUPAC DNA characters in a DNA sequence).

Format Parameters
-----------------
The following parameters are available to change how FASTA files are read or
written in scikit-bio.

Reader Parameters
^^^^^^^^^^^^^^^^^
The available parameters differ depending on which reader is used.

Generator, SequenceCollection, and Alignment Reader Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``constructor`` parameter can be used with the ``BiologicalSequence``
generator, ``SequenceCollection``, and ``Alignment`` FASTA readers.
``constructor`` specifies the in-memory type of each sequence that is parsed,
and defaults to ``BiologicalSequence``. ``constructor`` should be a subclass of
``BiologicalSequence``. For example, if you know that the FASTA file you're
reading contains protein sequences, you would pass
``constructor=ProteinSequence`` to the reader call.

.. note:: The FASTA sniffer will not attempt to guess the ``constructor``
   parameter, so it will always default to ``BiologicalSequence`` if another
   type is not provided to the reader. The sniffer could attempt to infer the
   type of sequences contained in the file, but this process could be
   error-prone since sequence type is not encoded in the FASTA file format
   itself (and can be ambiguous). This could produce strange or unintended
   behavior in certain cases, so we defer to the user to provide more specific
   sequence type information if it is available.

BiologicalSequence and Subclass Reader Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``seq_num`` parameter can be used with the ``BiologicalSequence``,
``NucleotideSequence``, ``DNASequence``, ``RNASequence``, and
``ProteinSequence`` FASTA readers. ``seq_num`` specifies which sequence to read
from the FASTA file, and defaults to 1 (i.e., such that the first sequence is
read). For example, to read the 50th sequence from a FASTA file, you would pass
``seq_num=50`` to the reader call.

.. note:: The FASTA sniffer will not attempt to guess the ``seq_num``
   parameter, so it will always default to reading the first sequence in the
   file unless overridden by the user. The sniffer cannot provide a reasonable
   guess for this parameter as it is entirely up to the user to specify which
   sequence to read.

Writer Parameters
^^^^^^^^^^^^^^^^^
The following parameters are available to all FASTA format writers:

- ``id_whitespace_replacement``: string to replace **each** whitespace
  character in a sequence ID. This parameter is useful for cases where an
  in-memory sequence ID contains whitespace, which would result in an on-disk
  representation that would not be read back into memory as the same ID (since
  IDs in FASTA format cannot contain whitespace). Defaults to ``_``. If
  ``None``, no whitespace replacement is performed and IDs are written as they
  are stored in memory (this has the potential to create an invalid
  FASTA-formatted file; see note below).
- ``description_newline_replacement``: string to replace **each** newline
  character in a sequence description. Since a FASTA header must be a single
  line, newlines are not allowed in sequence descriptions and must be replaced
  in order to write a valid FASTA file. Defaults to a single space. If
  ``None``, no newline replacement is performed and descriptions are written as
  they are stored in memory (this has the potential to create an invalid
  FASTA-formatted file; see note below).

- ``max_width``: integer specifying the maximum line width (i.e., number of
  characters) for sequence data. If a sequence is longer than ``max_width``, it
  will be split across multiple lines, each with a maximum width of
  ``max_width``. Default is to not split across multiple lines.

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
Suppose we have the following FASTA file with five aligned sequences (example
modified from [5]_)::

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

Note that the sequences are not required to be aligned in order for the file to
be a valid FASTA file (this depends on the object that you're reading the file
into). Also note that some of the sequences occur on a single line, while
others are split across multiple lines.

Let's define this file in-memory as a ``StringIO``, though this could be a real
file path, file handle, etc. in practice (anything that's supported by
scikit-bio's I/O registry):

>>> from StringIO import StringIO
>>> fh = StringIO(
...     ">seq1 Turkey\\n"
...     "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT\\n"
...     ">seq2 Salmo gair\\n"
...     "AAGCCTTGGCAGTGCAGGGTGAGCCGTGG\\n"
...     "CCGGGCACGGTAT\\n"
...     ">seq3 H. Sapiens\\n"
...     "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA\\n"
...     ">seq4 Chimp\\n"
...     "AAACCCTTGCCG\\n"
...     "TTACGCTTAAAC\\n"
...     "CGAGGCCGGGAC\\n"
...     "ACTCAT\\n"
...     ">seq5 Gorilla\\n"
...     "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA\\n")

Let's read the FASTA file into a ``SequenceCollection``:

>>> from skbio import SequenceCollection
>>> sc = SequenceCollection.read(fh)
>>> sc.sequence_lengths()
[42, 42, 42, 42, 42]
>>> sc.ids()
['seq1', 'seq2', 'seq3', 'seq4', 'seq5']

We see that all 5 sequences have 42 characters, and that each of the sequence
IDs were successfully read into memory.

Since these sequences are aligned, let's load the FASTA file into a more
appropriate data structure:

>>> from skbio import Alignment
>>> fh.seek(0) # reset position to beginning of file so we can read again
>>> aln = Alignment.read(fh)
>>> aln.sequence_length()
42

Note that we were able to read the FASTA file into two different data
structures (``SequenceCollection`` and ``Alignment``) using the exact same
``read`` method call (and underlying reading/parsing logic). Also note that we
didn't specify a file format in the ``read`` call. The FASTA sniffer detected
the correct file format for us!

Let's inspect the type of sequences stored in the ``Alignment``:

>>> aln[0]
<BiologicalSequence: AAGCTNGGGC... (length: 42)>

By default, sequences are loaded as ``BiologicalSequence`` objects. We can
change the type of sequence via the ``constructor`` parameter:

>>> from skbio import DNASequence
>>> fh.seek(0) # reset position to beginning of file so we can read again
>>> aln = Alignment.read(fh, constructor=DNASequence)
>>> aln[0]
<DNASequence: AAGCTNGGGC... (length: 42)>

We now have an ``Alignment`` of ``DNASequence`` objects instead of
``BiologicalSequence`` objects. Validation of sequence character data is not
performed during reading (see warning above for details). To verify that each
of the sequences are valid DNA sequences:

>>> aln.is_valid()
True

To write the alignment in FASTA format:

>>> new_fh = StringIO()
>>> aln.write(new_fh)
>>> print(new_fh.getvalue())
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
>>> new_fh.close()

Both ``SequenceCollection`` and ``Alignment`` load all of the sequences from
the FASTA file into memory at once. If the FASTA file is large (which is often
the case), this may be infeasible if you don't have enough memory. To work
around this issue, you can stream the sequences using scikit-bio's
generator-based FASTA reader and writer. The generator-based reader yields
``BiologicalSequence`` objects (or subclasses if ``constructor`` is supplied)
one at a time, instead of loading all sequences into memory. For example, let's
use the generator-based reader to process a single sequence at a time in a
``for`` loop:

>>> from skbio.io import read
>>> fh.seek(0) # reset position to beginning of file so we can read again
>>> for seq in read(fh, format='fasta'):
...     seq
<BiologicalSequence: AAGCTNGGGC... (length: 42)>
<BiologicalSequence: AAGCCTTGGC... (length: 42)>
<BiologicalSequence: ACCGGTTGGC... (length: 42)>
<BiologicalSequence: AAACCCTTGC... (length: 42)>
<BiologicalSequence: AAACCCTTGC... (length: 42)>

A single sequence can also be read into a ``BiologicalSequence`` (or subclass):

>>> from skbio import BiologicalSequence
>>> fh.seek(0) # reset position to beginning of file so we can read again
>>> BiologicalSequence.read(fh)
<BiologicalSequence: AAGCTNGGGC... (length: 42)>

By default, the first sequence in the FASTA file is read. This can be
controlled with ``seq_num``. For example, to read the fifth sequence:

>>> fh.seek(0) # reset position to beginning of file so we can read again
>>> BiologicalSequence.read(fh, seq_num=5)
<BiologicalSequence: AAACCCTTGC... (length: 42)>

We can use the same API to read the fifth sequence into a ``DNASequence``:

>>> fh.seek(0) # reset position to beginning of file so we can read again
>>> dna_seq = DNASequence.read(fh, seq_num=5)
>>> dna_seq
<DNASequence: AAACCCTTGC... (length: 42)>

Individual sequence objects can also be written in FASTA format:

>>> new_fh = StringIO()
>>> dna_seq.write(new_fh)
>>> print(new_fh.getvalue())
>seq5 Gorilla
AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA
<BLANKLINE>
>>> new_fh.close()

References
----------
.. [1] Lipman, DJ; Pearson, WR (1985). "Rapid and sensitive protein similarity
   searches". Science 227 (4693): 1435-41.
.. [2] http://en.wikipedia.org/wiki/FASTA_format
.. [3] http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
.. [4] Madden T. The BLAST Sequence Analysis Tool. 2002 Oct 9
   [Updated 2003 Aug 13]. In: McEntyre J, Ostell J, editors. The NCBI Handbook
   [Internet]. Bethesda (MD): National Center for Biotechnology Information
   (US); 2002-. Chapter 16. Available from:
   http://www.ncbi.nlm.nih.gov/books/NBK21097/
.. [5] http://evolution.genetics.washington.edu/phylip/doc/sequence.html

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import range, zip
from future.standard_library import hooks
with hooks():
    from itertools import zip_longest

import re
import textwrap

import numpy as np

from skbio.io import (register_reader, register_writer, register_sniffer,
                      FASTAFormatError, FASTAFormatError, FASTAFormatError)
from skbio.io._base import _chunk_str
from skbio.alignment import SequenceCollection, Alignment
from skbio.sequence import (BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence)
from skbio.util import cardinal_to_ordinal


@register_sniffer('fasta')
def _fasta_sniffer(fh):
    # Strategy:
    #   If the file appears to be valid FASTA, check if it's also valid QUAL.
    #   If it is, do not identify as FASTA since the QUAL sniffer is more
    #   specific (i.e., it attempts to parse the quality scores as integers).
    valid_fasta = _fasta_or_qual_sniffer(fh, 'fasta')[0]
    if valid_fasta:
        valid_fasta = not _qual_sniffer(fh)[0]
    return valid_fasta, {}


@register_sniffer('qual')
def _qual_sniffer(fh):
    return _fasta_or_qual_sniffer(fh, 'qual')


# TODO is this is correct way to handle the constructor kwarg in a compound
# format? should there be a check that the two values are exactly the same?
@register_reader(['fasta', 'qual'])
def _fasta_qual_to_generator(fasta_fh, qual_fh,
                             constructor=(BiologicalSequence,
                                          BiologicalSequence)):
    # TODO is this try/finally necessary anymore?
    try:
        fasta_gen = _fasta_or_qual_to_generator(fasta_fh, 'fasta')
        qual_gen = _fasta_or_qual_to_generator(qual_fh, 'qual')

        for fasta_rec, qual_rec in zip_longest(fasta_gen, qual_gen,
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
                    % (fasta_id, qual_id))
            if fasta_desc != qual_desc:
                raise FASTAFormatError(
                    "Descriptions do not match between FASTA and QUAL "
                    "records: %r != %r" % (fasta_desc, qual_desc))

            yield constructor[0](fasta_seq, id=fasta_id,
                                 description=fasta_desc, quality=qual_scores)
    finally:
        fasta_gen.close()
        qual_gen.close()


@register_reader(['fasta', 'qual'], BiologicalSequence)
def _fasta_qual_to_biological_sequence(fasta_fh, qual_fh, seq_num=(1, 1)):
    return _fasta_or_fasta_qual_to_sequence(
        (fasta_fh, qual_fh), seq_num[0],
        (BiologicalSequence, BiologicalSequence), 'fasta, qual')


@register_reader(['fasta', 'qual'], NucleotideSequence)
def _fasta_qual_to_nucleotide_sequence(fasta_fh, qual_fh, seq_num=(1, 1)):
    return _fasta_or_fasta_qual_to_sequence(
        (fasta_fh, qual_fh), seq_num[0],
        (NucleotideSequence, NucleotideSequence), 'fasta, qual')


@register_reader(['fasta', 'qual'], DNASequence)
def _fasta_qual_to_dna_sequence(fasta_fh, qual_fh, seq_num=(1, 1)):
    return _fasta_or_fasta_qual_to_sequence(
        (fasta_fh, qual_fh), seq_num[0],
        (DNASequence, DNASequence), 'fasta, qual')


@register_reader(['fasta', 'qual'], RNASequence)
def _fasta_qual_to_rna_sequence(fasta_fh, qual_fh, seq_num=(1, 1)):
    return _fasta_or_fasta_qual_to_sequence(
        (fasta_fh, qual_fh), seq_num[0],
        (RNASequence, RNASequence), 'fasta, qual')


@register_reader(['fasta', 'qual'], ProteinSequence)
def _fasta_qual_to_protein_sequence(fasta_fh, qual_fh, seq_num=(1, 1)):
    return _fasta_or_fasta_qual_to_sequence(
        (fasta_fh, qual_fh), seq_num[0],
        (ProteinSequence, ProteinSequence), 'fasta, qual')


@register_reader(['fasta', 'qual'], SequenceCollection)
def _fasta_qual_to_sequence_collection(fasta_fh, qual_fh,
                                       constructor=(BiologicalSequence,
                                                    BiologicalSequence)):
    return SequenceCollection(
        list(_fasta_qual_to_generator(fasta_fh, qual_fh,
                                      constructor=constructor)))


@register_reader(['fasta', 'qual'], Alignment)
def _fasta_qual_to_alignment(fasta_fh, qual_fh,
                             constructor=(BiologicalSequence,
                                          BiologicalSequence)):
    return Alignment(
        list(_fasta_qual_to_generator(fasta_fh, qual_fh,
                                      constructor=constructor)))


@register_writer(['fasta', 'qual'])
def _generator_to_fasta_qual(obj, fasta_fh, qual_fh,
                             id_whitespace_replacement=('_', '_'),
                             description_newline_replacement=(' ', ' '),
                             max_width=(None, None)):
    _generator_to_fasta_or_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement[0],
        description_newline_replacement[0], max_width[0])


@register_writer(['fasta', 'qual'], BiologicalSequence)
def _biological_sequence_to_fasta_qual(obj, fasta_fh, qual_fh,
                                       id_whitespace_replacement=('_', '_'),
                                       description_newline_replacement=(' ',
                                                                        ' '),
                                       max_width=(None, None)):
    _sequence_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], NucleotideSequence)
def _nucleotide_sequence_to_fasta_qual(obj, fasta_fh, qual_fh,
                                       id_whitespace_replacement=('_', '_'),
                                       description_newline_replacement=(' ',
                                                                        ' '),
                                       max_width=(None, None)):
    _sequence_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], DNASequence)
def _dna_sequence_to_fasta_qual(obj, fasta_fh, qual_fh,
                                id_whitespace_replacement=('_', '_'),
                                description_newline_replacement=(' ', ' '),
                                max_width=(None, None)):
    _sequence_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], RNASequence)
def _rna_sequence_to_fasta_qual(obj, fasta_fh, qual_fh,
                                id_whitespace_replacement=('_', '_'),
                                description_newline_replacement=(' ', ' '),
                                max_width=(None, None)):
    _sequence_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], ProteinSequence)
def _protein_sequence_to_fasta_qual(obj, fasta_fh, qual_fh,
                                    id_whitespace_replacement=('_', '_'),
                                    description_newline_replacement=(' ', ' '),
                                    max_width=(None, None)):
    _sequence_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], SequenceCollection)
def _sequence_collection_to_fasta_qual(obj, fasta_fh, qual_fh,
                                       id_whitespace_replacement=('_', '_'),
                                       description_newline_replacement=(' ',
                                                                        ' '),
                                       max_width=(None, None)):
    _sequences_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], Alignment)
def _alignment_to_fasta_qual(obj, fasta_fh, qual_fh,
                             id_whitespace_replacement=('_', '_'),
                             description_newline_replacement=(' ', ' '),
                             max_width=(None, None)):
    _sequences_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_reader('fasta')
def _fasta_to_generator(fh, constructor=BiologicalSequence):
    for seq, id_, desc in _fasta_or_qual_to_generator(fh, format='fasta'):
        yield constructor(seq, id=id_, description=desc)


@register_reader('fasta', BiologicalSequence)
def _fasta_to_biological_sequence(fh, seq_num=1):
    return _fasta_or_fasta_qual_to_sequence((fh,), seq_num, BiologicalSequence,
                                            format='fasta')


@register_reader('fasta', NucleotideSequence)
def _fasta_to_nucleotide_sequence(fh, seq_num=1):
    return _fasta_or_fasta_qual_to_sequence((fh,), seq_num, NucleotideSequence,
                                            format='fasta')


@register_reader('fasta', DNASequence)
def _fasta_to_dna_sequence(fh, seq_num=1):
    return _fasta_or_fasta_qual_to_sequence((fh,), seq_num, DNASequence,
                                            format='fasta')


@register_reader('fasta', RNASequence)
def _fasta_to_rna_sequence(fh, seq_num=1):
    return _fasta_or_fasta_qual_to_sequence((fh,), seq_num, RNASequence,
                                            format='fasta')


@register_reader('fasta', ProteinSequence)
def _fasta_to_protein_sequence(fh, seq_num=1):
    return _fasta_or_fasta_qual_to_sequence((fh,), seq_num, ProteinSequence,
                                            format='fasta')


@register_reader('fasta', SequenceCollection)
def _fasta_to_sequence_collection(fh, constructor=BiologicalSequence):
    return SequenceCollection(
        list(_fasta_to_generator(fh, constructor=constructor)))


@register_reader('fasta', Alignment)
def _fasta_to_alignment(fh, constructor=BiologicalSequence):
    return Alignment(
        list(_fasta_to_generator(fh, constructor=constructor)))


@register_writer('fasta')
def _generator_to_fasta(obj, fh, id_whitespace_replacement='_',
                        description_newline_replacement=' ', max_width=None):
    _generator_to_fasta_or_fasta_qual(
        obj, fh, None, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer('fasta', BiologicalSequence)
def _biological_sequence_to_fasta(obj, fh, id_whitespace_replacement='_',
                                  description_newline_replacement=' ',
                                  max_width=None):
    _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width)


@register_writer('fasta', NucleotideSequence)
def _nucleotide_sequence_to_fasta(obj, fh, id_whitespace_replacement='_',
                                  description_newline_replacement=' ',
                                  max_width=None):
    _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width)


@register_writer('fasta', DNASequence)
def _dna_sequence_to_fasta(obj, fh, id_whitespace_replacement='_',
                           description_newline_replacement=' ',
                           max_width=None):
    _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width)


@register_writer('fasta', RNASequence)
def _rna_sequence_to_fasta(obj, fh, id_whitespace_replacement='_',
                           description_newline_replacement=' ',
                           max_width=None):
    _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width)


@register_writer('fasta', ProteinSequence)
def _protein_sequence_to_fasta(obj, fh, id_whitespace_replacement='_',
                               description_newline_replacement=' ',
                               max_width=None):
    _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width)


@register_writer('fasta', SequenceCollection)
def _sequence_collection_to_fasta(obj, fh, id_whitespace_replacement='_',
                                  description_newline_replacement=' ',
                                  max_width=None):
    _sequences_to_fasta(obj, fh, id_whitespace_replacement,
                        description_newline_replacement, max_width)


@register_writer('fasta', Alignment)
def _alignment_to_fasta(obj, fh, id_whitespace_replacement='_',
                        description_newline_replacement=' ', max_width=None):
    _sequences_to_fasta(obj, fh, id_whitespace_replacement,
                        description_newline_replacement, max_width)


def _fasta_or_qual_sniffer(fh, format):
    # Strategy:
    #   Read up to 10 records. If at least one record is read (i.e. the file
    #   isn't empty) and no errors are thrown during reading, assume the file
    #   is in FASTA or QUAL format.

    # TODO is the finally block still necessary?
    try:
        not_empty = False
        gen = _fasta_or_qual_to_generator(fh, format=format)
        for _ in zip(range(10), gen):
            not_empty = True
        return not_empty, {}
    except (FASTAFormatError, FASTAFormatError):
        return False, {}
    finally:
        gen.close()


def _fasta_or_qual_to_generator(fh, format):
    if format == 'fasta':
        data_parser = _parse_sequence_data
        error_type = FASTAFormatError
        format_label = 'FASTA'
    else:
        data_parser = _parse_quality_scores
        error_type = FASTAFormatError
        format_label = 'QUAL'

    line = next(fh)
    # header check inlined here and below for performance
    if line.startswith('>'):
        id_, desc = _parse_header(line)
    else:
        raise error_type(
            "Found line without a %s header:\n%s" % (format_label, line))

    data_chunks = []
    for line in fh:
        if line.startswith('>'):
            # new header, so yield current record and reset state
            yield data_parser(data_chunks), id_, desc
            data_chunks = []
            id_, desc = _parse_header(line)
        else:
            line = line.strip()
            if line:
                data_chunks.append(line)
            else:
                raise error_type("Found blank or whitespace-only line in "
                                 "%s-formatted file." % format_label)
    # yield last record in file
    yield data_parser(data_chunks), id_, desc


def _parse_header(line):
    id_ = ''
    desc = ''
    header = line[1:].rstrip()
    if header:
        if header[0].isspace():
            # no id
            desc = header.lstrip()
        else:
            header_tokens = header.split(None, 1)
            if len(header_tokens) == 1:
                # no description
                id_ = header_tokens[0]
            else:
                id_, desc = header_tokens
    return id_, desc


def _parse_sequence_data(chunks):
    if not chunks:
        raise FASTAFormatError("Found FASTA header without sequence data.")
    return ''.join(chunks)


def _parse_quality_scores(chunks):
    if not chunks:
        raise FASTAFormatError("Found QUAL header without quality scores.")

    qual_str = ' '.join(chunks)
    try:
        return np.asarray(qual_str.split(), dtype=int)
    except ValueError:
        raise FASTAFormatError(
            "Could not convert quality scores to integers:\n%s" % qual_str)


def _fasta_or_fasta_qual_to_sequence(fh, seq_num, constructor, format):
    if format == 'fasta':
        seq_generator = _fasta_to_generator
        error_type = FASTAFormatError
        format_label = 'FASTA'
    else:
        seq_generator = _fasta_qual_to_generator
        error_type = FASTAFormatError
        format_label = 'FASTA/QUAL'

    if seq_num < 1:
        raise error_type(
            "Invalid sequence number (seq_num=%d). seq_num must be between 1 "
            "and the number of sequences in the %s-formatted file (inclusive)."
            % (seq_num, format_label))

    seq_idx = seq_num - 1
    seq = None
    try:
        gen = seq_generator(*fh, constructor=constructor)
        for idx, curr_seq in enumerate(gen):
            if idx == seq_idx:
                seq = curr_seq
                break
    finally:
        gen.close()

    if seq is None:
        raise error_type(
            "Reached end of %s-formatted file before finding %s biological "
            "sequence." % (format_label, cardinal_to_ordinal(seq_num)))
    return seq


def _generator_to_fasta_or_fasta_qual(obj, fasta_fh, qual_fh,
                                      id_whitespace_replacement,
                                      description_newline_replacement,
                                      max_width):
    if qual_fh is None:
        error_type = FASTAFormatError
        format_label = 'FASTA'
    else:
        error_type = FASTAFormatError
        format_label = 'FASTA/QUAL'
        if max_width is not None:
            # define text wrapper for quality scores here for efficiency.
            # textwrap docs recommend reusing a TextWrapper instance when it is
            # used many times. configure text wrapper to never break words
            # (i.e., integer quality scores)
            qual_wrapper = textwrap.TextWrapper(width=max_width,
                                                break_long_words=False,
                                                break_on_hyphens=False)

    if ((id_whitespace_replacement is not None and
         '\n' in id_whitespace_replacement) or
        (description_newline_replacement is not None and
         '\n' in description_newline_replacement)):
        raise error_type(
            "Newline character (\\n) cannot be used to replace whitespace in "
            "biological sequence IDs, nor to replace newlines in biological "
            "sequence descriptions. Otherwise, the %s-formatted file will be "
            "invalid." % format_label)
    ws_pattern = re.compile(r'\s')
    nl_pattern = re.compile(r'\n')

    for idx, seq in enumerate(obj):
        if len(seq) < 1:
            raise error_type(
                "Cannot write %s biological sequence in %s format because it "
                "does not contain any characters (i.e., it is an empty/blank "
                "sequence). Empty sequences are not supported in the %s file "
                "format." % (cardinal_to_ordinal(idx + 1), format_label,
                             format_label))

        id_ = seq.id
        if id_whitespace_replacement is not None:
            id_ = re.sub(ws_pattern, id_whitespace_replacement, id_)

        desc = seq.description
        if description_newline_replacement is not None:
            desc = re.sub(nl_pattern, description_newline_replacement, desc)

        if desc:
            header = '%s %s' % (id_, desc)
        else:
            header = id_

        seq_str = str(seq)
        if max_width is not None:
            seq_str = _chunk_str(seq_str, max_width, '\n')

        fasta_fh.write('>%s\n%s\n' % (header, seq_str))

        if qual_fh is not None:
            if not seq.has_quality():
                raise error_type(
                    "Cannot write %s biological sequence in %s format because "
                    "it does not have quality scores associated with it." %
                    (cardinal_to_ordinal(idx + 1), format_label))

            qual_str = ' '.join(np.asarray(seq.quality, dtype=np.str))
            if max_width is not None:
                qual_str = qual_wrapper.fill(qual_str)

            qual_fh.write('>%s\n%s\n' % (header, qual_str))


def _sequence_to_fasta_qual(obj, fasta_fh, qual_fh, id_whitespace_replacement,
                            description_newline_replacement, max_width):
    def seq_gen():
        yield obj

    _generator_to_fasta_qual(
        seq_gen(), fasta_fh, qual_fh,
        id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        max_width=max_width)


def _sequences_to_fasta_qual(obj, fasta_fh, qual_fh, id_whitespace_replacement,
                             description_newline_replacement, max_width):
    def seq_gen():
        for seq in obj:
            yield seq

    _generator_to_fasta_qual(
        seq_gen(), fasta_fh, qual_fh,
        id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        max_width=max_width)


def _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width):
    def seq_gen():
        yield obj

    _generator_to_fasta(
        seq_gen(), fh, id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        max_width=max_width)


def _sequences_to_fasta(obj, fh, id_whitespace_replacement,
                        description_newline_replacement, max_width):
    def seq_gen():
        for seq in obj:
            yield seq

    _generator_to_fasta(
        seq_gen(), fh, id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        max_width=max_width)
