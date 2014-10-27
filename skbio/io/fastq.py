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

import re
from skbio.io import (register_reader, register_writer, register_sniffer,
                      FASTQFormatError)
from skbio.io._base import (_decode_qual_to_phred, _encode_phred_to_qual,
                            _get_nth_sequence, _parse_fasta_like_header,
                            _format_fasta_like_records)
from skbio.alignment import SequenceCollection, Alignment
from skbio.sequence import (BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence)

_whitespace_regex = re.compile(r'\s')


@register_sniffer('fastq')
def _fastq_sniffer(fh):
    # Strategy:
    #   Read up to 10 records. If at least one record is read (i.e. the file
    #   isn't empty) and the quality scores are in printable ASCII range,
    #   assume the file is FASTQ.
    try:
        not_empty = False
        for _ in zip(range(10), _fastq_to_generator(fh, phred_offset=33)):
            not_empty = True
        return not_empty, {}
    except (FASTQFormatError, ValueError):
        return False, {}


@register_reader('fastq')
def _fastq_to_generator(fh, variant=None, phred_offset=None,
                        constructor=BiologicalSequence):
    seq_header = next(_line_generator(fh))
    if not seq_header.startswith('@'):
        raise FASTQFormatError(
            "Expected sequence (@) header line at start of file: %r"
            % seq_header)

    while seq_header is not None:
        id_, desc = _parse_fasta_like_header(seq_header)
        seq, qual_header = _parse_sequence_data(fh)

        if qual_header != '+' and qual_header[1:] != seq_header[1:]:
            raise FASTQFormatError(
                "Sequence (@) and quality (+) header lines do not match: "
                "%r != %r" % (seq_header[1:], qual_header[1:]))

        phred_scores, seq_header = _parse_quality_scores(fh, len(seq), variant,
                                                         phred_offset)
        yield constructor(seq, id=id_, description=desc, quality=phred_scores)


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
    formatted_records = _format_fasta_like_records(
        obj, id_whitespace_replacement, description_newline_replacement, True)
    for header, seq_str, qual_scores in formatted_records:
        qual_str = _encode_phred_to_qual(qual_scores, variant=variant,
                                         phred_offset=phred_offset)
        fh.write('@')
        fh.write(header)
        fh.write('\n')
        fh.write(seq_str)
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


def _line_generator(fh):
    for line in fh:
        line = line.rstrip('\n')
        if not line:
            raise FASTQFormatError("Found blank line in FASTQ-formatted file.")
        yield line


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
            if _whitespace_regex.search(chunk):
                raise FASTQFormatError(
                    "Found whitespace in sequence data: %r" % chunk)
            seq_chunks.append(chunk)

    raise FASTQFormatError(
        "Found incomplete/truncated FASTQ record at end of file.")


def _parse_quality_scores(fh, seq_len, variant, phred_offset):
    phred_scores = []
    qual_len = 0
    for chunk in _line_generator(fh):
        if chunk.startswith('@') and qual_len == seq_len:
            return phred_scores, chunk
        else:
            qual_len += len(chunk)

            if qual_len > seq_len:
                raise FASTQFormatError(
                    "Found more quality score characters than sequence "
                    "characters. Extra quality score characters: %r" %
                    chunk[-(qual_len - seq_len):])

            phred_scores.extend(
                _decode_qual_to_phred(chunk, variant=variant,
                                      phred_offset=phred_offset))

    if qual_len != seq_len:
        raise FASTQFormatError(
            "Found incomplete/truncated FASTQ record at end of file.")
    return phred_scores, None


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
