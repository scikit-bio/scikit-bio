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
  IDs in FASTA format cannot contain whitespace). Defaults to ``_``.
- ``description_newline_replacement``: string to replace **each** newline
  character in a sequence description. Since a FASTA header must be a single
  line, newlines are not allowed in sequence descriptions and must be replaced
  in order to write a valid FASTA file. Defaults to a single space.
- ``max_width``: integer specifying the maximum line width (i.e., number of
  characters) for sequence data. If a sequence is longer than ``max_width``, it
  will be split across multiple lines, each with a maximum width of
  ``max_width``. Default is to not split across multiple lines.

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
around this issue, you can stream the sequences using scikit-bio's generator
FASTA reader and writer. The generator reader yields ``BiologicalSequence``
objects (or subclasses if ``constructor`` is supplied) one at a time, instead
of loading all sequences into memory. For example, let's use the generator
reader to process a single sequence at a time in a ``for`` loop:

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

import re

from skbio.io import (register_reader, register_writer, register_sniffer,
                      FASTAFormatError)
from skbio.io._base import _chunk_str
from skbio.alignment import SequenceCollection, Alignment
from skbio.sequence import (BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence)
from skbio.util import cardinal_to_ordinal


@register_sniffer('fasta')
def _fasta_sniffer(fh):
    # Strategy:
    #   Read up to 10 sequences (records). If at least one sequence is read
    #   (i.e. the file isn't empty) and no errors are thrown during reading,
    #   assume the file is in FASTA format.
    try:
        gen = _fasta_to_generator(fh)
        num_seqs = 0
        for _ in range(10):
            try:
                next(gen)
            except FASTAFormatError:
                return False, {}
            except StopIteration:
                break
            num_seqs += 1

        if num_seqs < 1:
            return False, {}
        else:
            return True, {}
    finally:
        gen.close()


@register_reader('fasta')
def _fasta_to_generator(fh, constructor=BiologicalSequence):
    line = next(fh)
    if _is_header(line):
        id_, desc = _parse_header(line)
    else:
        raise FASTAFormatError("Found line without a FASTA header:\n%s" % line)

    seq_chunks = []
    for line in fh:
        if _is_header(line):
            # new header, so yield current sequence and reset state
            yield _construct_sequence(constructor, seq_chunks, id_, desc)
            seq_chunks = []
            id_, desc = _parse_header(line)
        else:
            line = line.strip()
            if line:
                seq_chunks.append(line)
            else:
                raise FASTAFormatError("Found blank or whitespace-only line "
                                       "in FASTA-formatted file.")

    # yield last sequence in file
    yield _construct_sequence(constructor, seq_chunks, id_, desc)


@register_reader('fasta', BiologicalSequence)
def _fasta_to_biological_sequence(fh, seq_num=1):
    return _fasta_to_sequence(fh, seq_num, BiologicalSequence)


@register_reader('fasta', NucleotideSequence)
def _fasta_to_nucleotide_sequence(fh, seq_num=1):
    return _fasta_to_sequence(fh, seq_num, NucleotideSequence)


@register_reader('fasta', DNASequence)
def _fasta_to_dna_sequence(fh, seq_num=1):
    return _fasta_to_sequence(fh, seq_num, DNASequence)


@register_reader('fasta', RNASequence)
def _fasta_to_rna_sequence(fh, seq_num=1):
    return _fasta_to_sequence(fh, seq_num, RNASequence)


@register_reader('fasta', ProteinSequence)
def _fasta_to_protein_sequence(fh, seq_num=1):
    return _fasta_to_sequence(fh, seq_num, ProteinSequence)


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
    if '\n' in id_whitespace_replacement or \
            '\n' in description_newline_replacement:
        raise FASTAFormatError(
            "Newline character (\\n) cannot be used to replace whitespace in "
            "biological sequence IDs, nor to replace newlines in biological "
            "sequence descriptions. Otherwise, the FASTA-formatted file will "
            "be invalid.")
    ws_pattern = re.compile(r'\s')
    nl_pattern = re.compile(r'\n')

    for idx, seq in enumerate(obj):
        if len(seq) < 1:
            raise FASTAFormatError(
                "Cannot write %s biological sequence in FASTA format because "
                "it does not contain any characters (i.e., it is an "
                "empty/blank sequence). Empty sequences are not supported in "
                "the FASTA file format." % cardinal_to_ordinal(idx + 1))

        id_ = re.sub(ws_pattern, id_whitespace_replacement, seq.id)
        desc = re.sub(nl_pattern, description_newline_replacement,
                      seq.description)

        if desc:
            header = '%s %s' % (id_, desc)
        else:
            header = id_

        seq_str = str(seq)

        if max_width is not None:
            seq_str = _chunk_str(seq_str, max_width, '\n')

        fh.write('>%s\n%s\n' % (header, seq_str))


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


def _is_header(line):
    return line.startswith('>')


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


def _construct_sequence(constructor, seq_chunks, id_, description):
    if not seq_chunks:
        raise FASTAFormatError("Found FASTA header without sequence data.")
    return constructor(''.join(seq_chunks), id=id_, description=description)


def _fasta_to_sequence(fh, seq_num, constructor):
    if seq_num < 1:
        raise FASTAFormatError(
            "Invalid sequence number (seq_num=%d). seq_num must be between 1 "
            "and the number of sequences in the FASTA-formatted file "
            "(inclusive)." % seq_num)

    seq_idx = seq_num - 1
    seq = None
    try:
        gen = _fasta_to_generator(fh, constructor=constructor)
        for idx, curr_seq in enumerate(gen):
            if idx == seq_idx:
                seq = curr_seq
                break
    finally:
        gen.close()

    if seq is None:
        raise FASTAFormatError(
            "Reached end of FASTA-formatted file before finding %s biological "
            "sequence." % cardinal_to_ordinal(seq_num))
    return seq


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
