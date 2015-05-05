r"""
FASTQ format (:mod:`skbio.io.fastq`)
====================================

.. currentmodule:: skbio.io.fastq

The FASTQ file format (``fastq``) stores biological (e.g., nucleotide)
sequences and their quality scores in a simple plain text format that is both
human-readable and easy to parse. The file format was invented by Jim Mullikin
at the Wellcome Trust Sanger Institute but wasn't given a formal definition,
though it has informally become a standard file format for storing
high-throughput sequence data. More information about the format and its
variants can be found in [1]_ and [2]_.

Conceptually, a FASTQ file is similar to paired FASTA and QUAL files in that it
stores both biological sequences and their quality scores. FASTQ differs from
FASTA/QUAL because the quality scores are stored in the same file as the
biological sequence data.

An example FASTQ-formatted file containing two DNA sequences and their quality
scores::

    @seq1 description 1
    AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
    +
    ````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
    @seq2 description 2
    TATGTATATATAACATATACATATATACATACATA
    +
    ]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |generator of :mod:`skbio.sequence.BiologicalSequence` objects  |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.alignment.SequenceCollection`                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.alignment.Alignment`                               |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.BiologicalSequence`                       |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.NucleotideSequence`                       |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.DNASequence`                              |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.RNASequence`                              |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.ProteinSequence`                          |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
A FASTQ file contains one or more biological sequences and their corresponding
quality scores stored sequentially as *records*. Each *record* consists of four
sections:

1. Sequence header line consisting of a sequence identifier (ID) and
   description (both optional)
2. Biological sequence data (typically stored using the standard IUPAC
   lexicon), optionally split over multiple lines
3. Quality header line separating sequence data from quality scores (optionally
   repeating the ID and description from the sequence header line)
4. Quality scores as printable ASCII characters, optionally split over multiple
   lines. Decoding of quality scores will depend on the specified FASTQ variant
   (see below for more details)

For the complete FASTQ format specification, see [1]_. scikit-bio's FASTQ
implementation follows the format specification described in this excellent
publication, including validating the implementation against the FASTQ examples
provided in the publication's supplementary data.

.. note:: IDs and descriptions will be parsed from sequence header lines in
   exactly the same way as FASTA headers (:mod:`skbio.io.fasta`).

   Whitespace is not allowed in sequence data or quality scores. Leading and
   trailing whitespace is stripped from the file. Blank or whitespace-only
   lines are only permitted at the beginning of the file, between FASTQ
   records, or at the end of the file. A blank or whitespace-only line after
   the header line, within the sequence, or within quality scores is an
   error. If more than 5 blank or whitespace-only lines are at the beginning
   of the file, the sniffer will issue a warning.

   scikit-bio will write FASTQ files in a normalized format, with each record
   section on a single line. Thus, each record will be composed of *exactly*
   four lines. The quality header line won't have the sequence ID and
   description repeated.

Quality Score Variants
^^^^^^^^^^^^^^^^^^^^^^
FASTQ associates quality scores with sequence data, with each quality score
encoded as a single printable ASCII character. In scikit-bio, all quality
scores are decoded as Phred quality scores. This is the most common quality
score metric, though there are others (e.g., Solexa quality scores).
Unfortunately, different sequencers have different ways of encoding quality
scores as ASCII characters, notably Sanger and Illumina. Below is a table
highlighting the different encoding variants supported by scikit-bio, as well
as listing the equivalent variant names used in the Open Bioinformatics
Foundation (OBF) [3]_ projects (e.g., Biopython, BioPerl, etc.).

+-----------+---------+----+--------+-----------------------------------------+
| Variant   | ASCII   |Off\|Quality | Notes                                   |
|           | Range   |set |Range   |                                         |
+===========+=========+====+========+=========================================+
|sanger     |33 to 126|33  |0 to 93 |Equivalent to OBF's fastq-sanger.        |
+-----------+---------+----+--------+-----------------------------------------+
|illumina1.3|64 to 126|64  |0 to 62 |Equivalent to OBF's fastq-illumina. Use  |
|           |         |    |        |this if your data was generated using    |
|           |         |    |        |Illumina 1.3-1.7 software.               |
+-----------+---------+----+--------+-----------------------------------------+
|illumina1.8|33 to 95 |33  |0 to 62 |Equivalent to sanger but with 0 to 62    |
|           |         |    |        |quality score range check. Use this if   |
|           |         |    |        |your data was generated using Illumina   |
|           |         |    |        |1.8 software or later.                   |
+-----------+---------+----+--------+-----------------------------------------+
|solexa     |59 to 126|64  |-5 to 62|Not currently implemented.               |
+-----------+---------+----+--------+-----------------------------------------+

.. note:: When writing, Phred quality scores will be truncated to the maximum
   value in the variant's range and a warning will be issued. This is
   consistent with the OBF projects.

   When reading, an error will be raised if a decoded quality score is outside
   the variant's range.

Format Parameters
-----------------
The following parameters are available to all FASTQ format readers and writers:

- ``variant``: A string indicating the quality score variant used to
  decode/encode Phred quality scores. Must be one of ``sanger``,
  ``illumina1.3``, ``illumina1.8``, or ``solexa``. This parameter is preferred
  over ``phred_offset`` because additional quality score range checks and
  conversions can be performed. It is also more explicit.

- ``phred_offset``: An integer indicating the ASCII code offset used to
  decode/encode Phred quality scores. Must be in the range ``[33, 126]``. All
  decoded scores will be assumed to be Phred scores (i.e., no additional
  conversions are performed). Prefer using ``variant`` over this parameter
  whenever possible.

.. note:: You must provide ``variant`` or ``phred_offset`` when reading or
   writing a FASTQ file. ``variant`` and ``phred_offset`` cannot both be
   provided at the same time.

The following additional parameters are the same as in FASTA format
(:mod:`skbio.io.fasta`):

- ``constructor``: see ``constructor`` parameter in FASTA format

- ``seq_num``: see ``seq_num`` parameter in FASTA format

- ``id_whitespace_replacement``: see ``id_whitespace_replacement`` parameter in
  FASTA format

- ``description_newline_replacement``: see ``description_newline_replacement``
  parameter in FASTA format

Examples
--------
Suppose we have the following FASTQ file with two DNA sequences::

    @seq1 description 1
    AACACCAAACTTCTCCACC
    ACGTGAGCTACAAAAGGGT
    +seq1 description 1
    ''''Y^T]']C^CABCACC
    `^LB^CCYT\T\Y\WF^^^
    @seq2 description 2
    TATGTATATATAACATATACATATATACATACATA
    +
    ]KZ[PY]_[YY^'''AC^\\'BT''C'\AT''BBB

Note that the first sequence and its quality scores are split across multiple
lines, while the second sequence and its quality scores are each on a single
line. Also note that the first sequence has a duplicate ID and description on
the quality header line, while the second sequence does not.

Let's define this file in-memory as a ``StringIO``, though this could be a real
file path, file handle, or anything that's supported by scikit-bio's I/O
registry in practice:

>>> from StringIO import StringIO
>>> fs = '\n'.join([
...     r"@seq1 description 1",
...     r"AACACCAAACTTCTCCACC",
...     r"ACGTGAGCTACAAAAGGGT",
...     r"+seq1 description 1",
...     r"''''Y^T]']C^CABCACC",
...     r"'^LB^CCYT\T\Y\WF^^^",
...     r"@seq2 description 2",
...     r"TATGTATATATAACATATACATATATACATACATA",
...     r"+",
...     r"]KZ[PY]_[YY^'''AC^\\'BT''C'\AT''BBB"])
>>> fh = StringIO(fs)

To load the sequences into a ``SequenceCollection``, we run:

>>> from skbio import SequenceCollection
>>> sc = SequenceCollection.read(fh, variant='sanger')
>>> sc
<SequenceCollection: n=2; mean +/- std length=36.50 +/- 1.50>

Note that quality scores are decoded from Sanger. To load the second sequence
as a ``DNASequence``:

>>> from skbio import DNASequence
>>> fh = StringIO(fs) # reload the StringIO to read from the beginning again
>>> DNASequence.read(fh, variant='sanger', seq_num=2)
<DNASequence: TATGTATATA... (length: 35)>

To write our ``SequenceCollection`` to a FASTQ file with quality scores encoded
using the ``illumina1.3`` variant:

>>> new_fh = StringIO()
>>> sc.write(new_fh, format='fastq', variant='illumina1.3')
>>> print(new_fh.getvalue())
@seq1 description 1
AACACCAAACTTCTCCACCACGTGAGCTACAAAAGGGT
+
FFFFx}s|F|b}b`ab`bbF}ka}bbxs{s{x{ve}}}
@seq2 description 2
TATGTATATATAACATATACATATATACATACATA
+
|jyzox|~zxx}FFF`b}{{FasFFbF{`sFFaaa
<BLANKLINE>
>>> new_fh.close()

Note that the file has been written in normalized format: sequence and quality
scores each only occur on a single line and the sequence header line is
not repeated in the quality header line. Note also that the quality scores are
different because they have been encoded using a different variant.

References
----------
.. [1] Peter J. A. Cock, Christopher J. Fields, Naohisa Goto, Michael L. Heuer,
   and Peter M. Rice. The Sanger FASTQ file format for sequences with quality
   scores, and the Solexa/Illumina FASTQ variants. Nucl. Acids Res. (2010) 38
   (6): 1767-1771. first published online December 16, 2009.
   doi:10.1093/nar/gkp1137
   http://nar.oxfordjournals.org/content/38/6/1767
.. [2] http://en.wikipedia.org/wiki/FASTQ_format
.. [3] http://www.open-bio.org/

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
                            _format_fasta_like_records, _line_generator,
                            _too_many_blanks)
from skbio.alignment import SequenceCollection, Alignment
from skbio.sequence import (BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence)

_whitespace_regex = re.compile(r'\s')


@register_sniffer('fastq')
def _fastq_sniffer(fh):

    if _too_many_blanks(fh, 5):
        return False, {}

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
    # SKip any blank or whitespace-only lines at beginning of file
    seq_header = next(_line_generator(fh, skip_blanks=True))

    if not seq_header.startswith('@'):
        raise FASTQFormatError(
            "Expected sequence (@) header line at start of file: %r"
            % seq_header)

    while seq_header is not None:
        id_, desc = _parse_fasta_like_header(seq_header)
        seq, qual_header = _parse_sequence_data(fh, seq_header)

        if qual_header != '+' and qual_header[1:] != seq_header[1:]:
            raise FASTQFormatError(
                "Sequence (@) and quality (+) header lines do not match: "
                "%r != %r" % (seq_header[1:], qual_header[1:]))

        phred_scores, seq_header = _parse_quality_scores(fh, len(seq),
                                                         variant,
                                                         phred_offset,
                                                         qual_header)
        yield constructor(seq, id=id_, description=desc,
                          quality=phred_scores)


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


def _blank_error(unique_text):
    error_string = ("Found blank or whitespace-only line {} in "
                    "FASTQ file").format(unique_text)
    raise FASTQFormatError(error_string)


def _parse_sequence_data(fh, prev):
    seq_chunks = []
    for chunk in _line_generator(fh):
        if chunk.startswith('+'):
            if not prev:
                _blank_error("before '+'")
            if not seq_chunks:
                raise FASTQFormatError(
                    "Found FASTQ record without sequence data.")
            return ''.join(seq_chunks), chunk
        elif chunk.startswith('@'):
            raise FASTQFormatError(
                "Found FASTQ record that is missing a quality (+) header line "
                "after sequence data.")
        else:
            if not prev:
                _blank_error("after header or within sequence")
            if _whitespace_regex.search(chunk):
                raise FASTQFormatError(
                    "Found whitespace in sequence data: %r" % chunk)
            seq_chunks.append(chunk)
        prev = chunk

    raise FASTQFormatError(
        "Found incomplete/truncated FASTQ record at end of file.")


def _parse_quality_scores(fh, seq_len, variant, phred_offset, prev):
    phred_scores = []
    qual_len = 0
    for chunk in _line_generator(fh):
        if chunk:
            if chunk.startswith('@') and qual_len == seq_len:
                return phred_scores, chunk
            else:
                if not prev:
                    _blank_error("after '+' or within quality scores")
                qual_len += len(chunk)

                if qual_len > seq_len:
                    raise FASTQFormatError(
                        "Found more quality score characters than sequence "
                        "characters. Extra quality score characters: %r" %
                        chunk[-(qual_len - seq_len):])

                phred_scores.extend(
                    _decode_qual_to_phred(chunk, variant=variant,
                                          phred_offset=phred_offset))
        prev = chunk

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
