r"""
FASTQ format (:mod:`skbio.io.fastq`)
====================================

.. currentmodule:: skbio.io.fastq

The FASTQ file format (``fastq``) stores biological (i.e., nucleotide or
protein) sequences and their quality scores in a simple plain text format that
is both human-readable and easy to parse. The file format was invented by
Jim Mullikin at the Wellcome Trust Sanger Institute but wasn't given a formal
definition, though it has informally become a standard file format for storing
high-throughput sequence data. More information about this format and its
variants can be found in [1]_ and [2]_ (much of the information here is drawn
from these resources).

Conceptually, a FASTQ file is similar to paired FASTA and QUAL files in that
both biological sequences and their quality scores are included. FASTQ differs
from FASTA/QUAL as the quality scores are stored in the same file as the
biological sequence data.

An example of a FASTQ-formatted file containing one DNA sequence and its
quality scores::

    @GAPC_0015:6:1:1314:13295#0/1
    AATATTGCTTTGTCTGAACGATAGTGCTCTTTGAT
    +GAPC_0015:6:1:1314:13295#0/1
    cLcc\\dddddaaYd`T```bLYT\`a```bZccc

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
A FASTQ file contains one or more biological sequences and their corresponding
quality scores stored sequentially as *records*. Each *record* consists of four
sections:

1. Header line consisting of a sequence identifier (ID) and optional
   description
2. Biological sequence data (typically stored using the standard IUPAC lexicon)
3. Line separating sequence data from quality scores (optionally including the
   ID and description from the header line)
4. Quality scores as printable ASCII characters

Each section is described in more detail below.

Header Line
^^^^^^^^^^^
Each record begins with a header line. The header line starts with an ``@``
character. Immediately following this character is a sequence identifier (ID)
and description separated by one or more whitespace characters. Both sequence
ID and description are optional and are represented as the empty string
(``''``) in scikit-bio's objects if they are not present in the header.

TODO finish after fasta docs are reviewed and merged

For example, consider the following header::

    @GAPC_0015:6:1:1314:13295#0/1  db-accession-149855

``GAPC_0015:6:1:1314:13295#0/1`` is the sequence ID and
``db-accession-149855`` is the sequence description

Quality Header
^^^^^^^^^^^^^^^
Each quality header consists of a single line beginning with a (``+``) symbol.
Immediately after this symbol is a quality identifier (ID) and an description
separated by one or more whitespace characters.

For example, consider the following header::

    +GAPC_0015:6:1:1314:13295#0/1  db-accession-149855

``GAPC_0015:6:1:1314:13295#0/1`` is the quality ID and
``db-accession-149855`` is the quality description

Examples
--------
Suppose we have a FASTQ file with two sequences

>>> from StringIO import StringIO
>>> fh = StringIO(
...        "@GAPC_0015:6:1:1259:10413#0/1\n"
...        "AACACCAAACTTCTCCACCACGTGAGCTACAAAAG\n"
...        "+GAPC_0015:6:1:1259:10413#0/1\n"
...        ";;;;;;;;;;;9;7;;.7;393333;;;;;;;;;;\n"
...        "@GAPC_0015:6:1:1283:11957#0/1\n"
...        "TATGTATATATAACATATACATATATACATACATA\n"
...        "+GAPC_0015:6:1:1283:11957#0/1\n"
...        "!''*((((***+))%%%++)(%%%%).1***-+*'\n")


Now read in the FASTQ file

>>> from skbio.io import read
>>> for seq in read(fh, format='fastq'):
...     seq
<BiologicalSequence: AACACCAAAC... (length: 35)>
<BiologicalSequence: TATGTATATA... (length: 35)>

Assume we have a fastq formatted file with the following contents::

    @seq1
    AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
    +
    ````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
    @seq2
    TATGTATATATAACATATACATATATACATACATA
    +
    ]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb

We can use the following code:

>>> from StringIO import StringIO
>>> from skbio.parse.sequences import parse_fastq
>>> fastq_f = StringIO('@seq1\n'
...                     'AACACCAAACTTCTCCACCACGTGAGCTACAAAAG\n'
...                     '+\n'
...                     '````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF\n'
...                     '@seq2\n'
...                     'TATGTATATATAACATATACATATATACATACATA\n'
...                     '+\n'
...                     ']KZ[PY]_[YY^```ac^\\\`bT``c`\\aT``bbb\n')
>>> for record in _fastq_to_generator(fastq_f, phred_offset=64):
...     print(record.id)
...     print(record.sequence)
...     print(record.quality)
seq1
AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
[32 32 32 32 25 30 20 29 32 29 35 30 35 33 34 35 33 35 35 32 30 12 34 30 35
 35 25 20 28 20 28 25 28 23  6]
seq2
TATGTATATATAACATATACATATATACATACATA
[29 11 26 27 16 25 29 31 27 25 25 30 32 32 32 33 35 30 28 28 32 34 20 32 32
 35 32 28 33 20 32 32 34 34 34]

References
----------
.. [1] Peter J. A. Cock, Christopher J. Fields, Naohisa Goto, Michael L. Heuer,
   and Peter M. Rice. The Sanger FASTQ file format for sequences with quality
   scores, and the Solexa/Illumina FASTQ variants. Nucl. Acids Res. (2010) 38
   (6): 1767-1771. first published online December 16, 2009.
   doi:10.1093/nar/gkp1137
.. [2] http://en.wikipedia.org/wiki/FASTQ_format

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
from future import standard_library
with standard_library.hooks():
    from itertools import zip_longest

import re
import numpy as np
from skbio.io import (register_reader, register_writer,
                      register_sniffer,
                      FASTQFormatError)
from skbio.io._base import (_decode_qual_to_phred, _encode_phred_to_qual,
                            _get_nth_sequence, _parse_fasta_like_header)
from skbio.alignment import SequenceCollection, Alignment
from skbio.sequence import (BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence)

from skbio.util import cardinal_to_ordinal


@register_sniffer('fastq')
def _fastq_sniffer(fh):
    # Strategy:
    #   Read up to 10 sequences (records), and checks if one of the two
    #   phred scoring scheme is used.  If at least one sequence is read
    #   (i.e. the file isn't empty) and no errors are thrown during reading,
    #   assume the file is in FASTQ format.
    try:
        not_empty = False
        gen = _fastq_to_generator(fh)
        for _ in zip(range(10), gen):
            not_empty = True
        if not_empty:
            return not_empty, {'phred_offset': 33}
        else:
            return not_empty, {}
    except FASTQFormatError:
        try:
            fh.seek(0)
            not_empty = False
            gen = _fastq_to_generator(fh, phred_offset=64)
            for _ in zip(range(10), gen):
                not_empty = True
            return not_empty, {'phred_offset': 64}
        except FASTQFormatError:
            return False, {}


@register_reader('fastq')
def _fastq_to_generator(fh, variant=None, phred_offset=None,
                        constructor=BiologicalSequence):
    seq_header = next(_line_generator(fh))
    while seq_header is not None:
        if not seq_header.startswith('@'):
            raise FASTQFormatError(
                "Expected sequence header line to start with '@' character: "
                "%r" % seq_header)

        id_, desc = _parse_fasta_like_header(seq_header)
        seq, qual_header = _parse_sequence_data(fh)

        if qual_header != '+' and qual_header[1:] != seq_header[1:]:
            raise FASTQFormatError(
                "Sequence (@) and quality (+) header lines do not match: "
                "%r != %r" % (seq_header[1:], qual_header[1:]))

        phred_scores, seq_header = _parse_quality_scores(fh, len(seq), variant,
                                                         phred_offset)
        yield constructor(seq, id=id_, description=desc, quality=phred_scores)


def _parse_sequence_data(fh):
    seq_chunks = []
    for chunk in _line_generator(fh):
        if chunk.startswith('+'):
            if not seq_chunks:
                raise FASTQFormatError(
                    "Found FASTQ record without sequence data.")
            return ''.join(seq_chunks), chunk
        elif chunk.startswith('@'):
            raise FASTQFormatError(
                "Found FASTQ record that is missing a quality (+) header line "
                "after sequence data.")
        else:
            for c in chunk:
                if c.isspace():
                    raise FASTQFormatError(
                        "Found whitespace in sequence data: %r" % chunk)
            seq_chunks.append(chunk)

    if not seq_chunks:
        raise FASTQFormatError(
            "Found incomplete/truncated FASTQ record at end of file that is "
            "missing sequence data.")
    else:
        raise FASTQFormatError(
            "Found incomplete/truncated FASTQ record at end of file that is "
            "missing a quality (+) header line after sequence data.")


def _parse_quality_scores(fh, seq_len, variant, phred_offset):
    qual_chunks = []
    qual_len = 0
    for chunk in _line_generator(fh):
        if chunk.startswith('@') and qual_len == seq_len:
            phred_scores = _decode_qual_to_phred(''.join(qual_chunks), variant=variant, phred_offset=phred_offset)
            return phred_scores, chunk
        else:
            qual_chunks.append(chunk)
            qual_len += len(chunk)

            if qual_len > seq_len:
                raise FASTQFormatError(
                    "Found more quality score characters than sequence "
                    "characters. Extra quality score characters: %r" %
                    ''.join(qual_chunks)[seq_len:])

    if not qual_chunks:
        raise FASTQFormatError(
            "Found incomplete/truncated FASTQ record at end of file that is "
            "missing quality scores.")
    if qual_len != seq_len:
        raise FASTQFormatError(
            "Found FASTQ record at end of file with different number of "
            "quality score characters than sequence characters: %d != %d" %
            (qual_len, seq_len))

    phred_scores = _decode_qual_to_phred(''.join(qual_chunks), variant=variant, phred_offset=phred_offset)
    return phred_scores, None


def _line_generator(fh):
    for line in fh:
        line = line.rstrip('\n')
        if not line:
            raise FASTQFormatError("Found blank line in FASTQ-formatted file.")
        yield line


@register_reader('fastq', BiologicalSequence)
def _fastq_to_biological_sequence(fh, variant=None, phred_offset=None,
                                  seq_num=1):
    return _get_nth_sequence(
        _fastq_to_generator(fh, variant=variant, phred_offset=phred_offset,
                            constructor=BiologicalSequence),
        seq_num)


@register_reader('fastq', NucleotideSequence)
def _fastq_to_nucleotide_sequence(fh, variant=None, phred_offset=None,
                                  seq_num=1):
    return _get_nth_sequence(
        _fastq_to_generator(fh, variant=variant, phred_offset=phred_offset,
                            constructor=NucleotideSequence),
        seq_num)


@register_reader('fastq', DNASequence)
def _fastq_to_dna_sequence(fh, variant=None, phred_offset=None, seq_num=1):
    return _get_nth_sequence(
        _fastq_to_generator(fh, variant=variant, phred_offset=phred_offset,
                            constructor=DNASequence),
        seq_num)


@register_reader('fastq', RNASequence)
def _fastq_to_rna_sequence(fh, variant=None, phred_offset=None, seq_num=1):
    return _get_nth_sequence(
        _fastq_to_generator(fh, variant=variant, phred_offset=phred_offset,
                            constructor=RNASequence),
        seq_num)


@register_reader('fastq', ProteinSequence)
def _fastq_to_protein_sequence(fh, variant=None, phred_offset=None, seq_num=1):
    return _get_nth_sequence(
        _fastq_to_generator(fh, variant=variant, phred_offset=phred_offset,
                            constructor=ProteinSequence),
        seq_num)


@register_reader('fastq', SequenceCollection)
def _fastq_to_sequence_collection(fh, variant=None, phred_offset=None,
                                  constructor=BiologicalSequence):
    return SequenceCollection(
        list(_fastq_to_generator(fh, variant=variant,
                                 phred_offset=phred_offset,
                                 constructor=constructor)))


@register_reader('fastq', Alignment)
def _fastq_to_alignment(fh, variant=None, phred_offset=None,
                        constructor=BiologicalSequence):
    return Alignment(
        list(_fastq_to_generator(fh, variant=variant,
                                 phred_offset=phred_offset,
                                 constructor=constructor)))


@register_writer('fastq')
def _generator_to_fastq(obj, fh, variant=None, phred_offset=None,
                        id_whitespace_replacement='_',
                        description_newline_replacement=' '):
    if ((id_whitespace_replacement is not None and
         '\n' in id_whitespace_replacement) or
        (description_newline_replacement is not None and
         '\n' in description_newline_replacement)):
        raise FASTQFormatError(
            "Newline character (\\n) cannot be used to replace whitespace in "
            "biological sequence IDs, nor to replace newlines in biological "
            "sequence descriptions. Otherwise, the FASTQ-formatted file will "
            "be invalid.")
    ws_pattern = re.compile(r'\s')
    nl_pattern = re.compile(r'\n')

    for idx, seq in enumerate(obj):
        if len(seq) < 1:
            raise FASTQFormatError(
                "Cannot write %s biological sequence in FASTQ format because "
                "it does not contain any characters (i.e., it is an "
                "empty/blank sequence). Empty sequences are not supported in "
                "the FASTQ file format." % cardinal_to_ordinal(idx + 1))

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

        if not seq.has_quality():
            raise FASTQFormatError(
                "Cannot write %s biological sequence in FASTQ format because "
                "it does not have quality scores associated with it." %
                cardinal_to_ordinal(idx + 1))
        qual_str = _encode_phred_to_qual(seq.quality, variant=variant,
                                         phred_offset=phred_offset)

        fh.write('@')
        fh.write(header)
        fh.write('\n')
        fh.write(str(seq))
        fh.write('\n+\n')
        fh.write(qual_str)
        fh.write('\n')


@register_writer('fastq', BiologicalSequence)
def _biological_sequence_to_fastq(obj, fh, variant=None, phred_offset=None,
                                  id_whitespace_replacement='_',
                                  description_newline_replacement=' '):
    _sequences_to_fastq([obj], fh, variant, phred_offset,
                        id_whitespace_replacement,
                        description_newline_replacement)


@register_writer('fastq', NucleotideSequence)
def _nucleotide_sequence_to_fastq(obj, fh, variant=None, phred_offset=None,
                                  id_whitespace_replacement='_',
                                  description_newline_replacement=' '):
    _sequences_to_fastq([obj], fh, variant, phred_offset,
                        id_whitespace_replacement,
                        description_newline_replacement)


@register_writer('fastq', DNASequence)
def _dna_sequence_to_fastq(obj, fh, variant=None, phred_offset=None,
                           id_whitespace_replacement='_',
                           description_newline_replacement=' '):
    _sequences_to_fastq([obj], fh, variant, phred_offset,
                        id_whitespace_replacement,
                        description_newline_replacement)


@register_writer('fastq', RNASequence)
def _rna_sequence_to_fastq(obj, fh, variant=None, phred_offset=None,
                           id_whitespace_replacement='_',
                           description_newline_replacement=' '):
    _sequences_to_fastq([obj], fh, variant, phred_offset,
                        id_whitespace_replacement,
                        description_newline_replacement)


@register_writer('fastq', ProteinSequence)
def _protein_sequence_to_fastq(obj, fh, variant=None, phred_offset=None,
                               id_whitespace_replacement='_',
                               description_newline_replacement=' '):
    _sequences_to_fastq([obj], fh, variant, phred_offset,
                        id_whitespace_replacement,
                        description_newline_replacement)


@register_writer('fastq', SequenceCollection)
def _sequence_collection_to_fastq(obj, fh, variant=None, phred_offset=None,
                                  id_whitespace_replacement='_',
                                  description_newline_replacement=' '):
    _sequences_to_fastq(obj, fh, variant, phred_offset,
                        id_whitespace_replacement,
                        description_newline_replacement)


@register_writer('fastq', Alignment)
def _alignment_to_fastq(obj, fh, variant=None, phred_offset=None,
                        id_whitespace_replacement='_',
                        description_newline_replacement=' '):
    _sequences_to_fastq(obj, fh, variant, phred_offset,
                        id_whitespace_replacement,
                        description_newline_replacement)


def _sequences_to_fastq(obj, fh, variant, phred_offset,
                        id_whitespace_replacement,
                        description_newline_replacement):
    def seq_gen():
        for seq in obj:
            yield seq

    _generator_to_fastq(
        seq_gen(), fh, variant=variant, phred_offset=phred_offset,
        id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement)
