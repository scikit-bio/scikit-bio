r"""FASTQ format (:mod:`skbio.io.format.fastq`)
===========================================

.. currentmodule:: skbio.io.format.fastq

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
scores:

.. code-block:: none

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
publication, including validating the implementation against the FASTQ example
files provided in the publication's supplementary data.

.. note:: IDs and descriptions will be parsed from sequence header lines in
   exactly the same way as FASTA headers (:mod:`skbio.io.format.fasta`). IDs,
   descriptions, and quality scores are also stored on, and written from,
   sequence objects in the same way as with FASTA.

.. note:: Blank or whitespace-only lines are only allowed at the beginning of
   the file, between FASTQ records, or at the end of the file. A blank or
   whitespace-only line after the header line, within the sequence, or within
   quality scores will raise an error.

   scikit-bio will ignore leading and trailing whitespace characters on each
   line while reading.

.. note:: Validation may be performed depending on the type of object the data
   is being read into. This behavior matches that of FASTA files.

.. note:: scikit-bio will write FASTQ files in a normalized format, with each
   record section on a single line. Thus, each record will be composed of
   *exactly* four lines. The quality header line won't have the sequence ID and
   description repeated.

.. note:: `lowercase` functionality is supported the same as with FASTA.

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
(:mod:`skbio.io.format.fasta`):

- ``constructor``: see ``constructor`` parameter in FASTA format

- ``seq_num``: see ``seq_num`` parameter in FASTA format

- ``id_whitespace_replacement``: see ``id_whitespace_replacement`` parameter in
  FASTA format

- ``description_newline_replacement``: see ``description_newline_replacement``
  parameter in FASTA format

- ``lowercase``: see ``lowercase`` parameter in FASTA format

Examples
--------
Suppose we have the following FASTQ file with two DNA sequences::

    @seq1 description 1
    AACACCAAACTTCTCCACC
    ACGTGAGCTACAAAAG
    +seq1 description 1
    ''''Y^T]']C^CABCACC
    `^LB^CCYT\T\Y\WF
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

>>> from io import StringIO
>>> fs = '\n'.join([
...     r"@seq1 description 1",
...     r"AACACCAAACTTCTCCACC",
...     r"ACGTGAGCTACAAAAG",
...     r"+seq1 description 1",
...     r"''''Y^T]']C^CABCACC",
...     r"'^LB^CCYT\T\Y\WF",
...     r"@seq2 description 2",
...     r"TATGTATATATAACATATACATATATACATACATA",
...     r"+",
...     r"]KZ[PY]_[YY^'''AC^\\'BT''C'\AT''BBB"])
>>> fh = StringIO(fs)

To load the sequences into a ``TabularMSA``, we run:

>>> from skbio import TabularMSA, DNA
>>> msa = TabularMSA.read(fh, constructor=DNA, variant='sanger')
>>> msa
TabularMSA[DNA]
-----------------------------------
Stats:
    sequence count: 2
    position count: 35
-----------------------------------
AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
TATGTATATATAACATATACATATATACATACATA

Note that quality scores are decoded from Sanger. To load the second sequence
as ``DNA``:

>>> fh = StringIO(fs) # reload the StringIO to read from the beginning again
>>> seq = DNA.read(fh, variant='sanger', seq_num=2)
>>> seq
DNA
----------------------------------------
Metadata:
    'description': 'description 2'
    'id': 'seq2'
Positional metadata:
    'quality': <dtype: uint8>
Stats:
    length: 35
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 14.29%
----------------------------------------
0 TATGTATATA TAACATATAC ATATATACAT ACATA

To write our ``TabularMSA`` to a FASTQ file with quality scores encoded using
the ``illumina1.3`` variant:

>>> new_fh = StringIO()
>>> print(msa.write(new_fh, format='fastq', variant='illumina1.3').getvalue())
@seq1 description 1
AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
+
FFFFx}s|F|b}b`ab`bbF}ka}bbxs{s{x{ve
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


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re

import numpy as np

from skbio.io import create_format, FASTQFormatError
from skbio.io.format._base import (
    _decode_qual_to_phred,
    _encode_phred_to_qual,
    _get_nth_sequence,
    _parse_fasta_like_header,
    _format_fasta_like_records,
    _line_generator,
    _too_many_blanks,
)
from skbio.alignment import TabularMSA
from skbio.sequence import Sequence, DNA, RNA, Protein

_whitespace_regex = re.compile(r"\s")


fastq = create_format("fastq")


@fastq.sniffer()
def _fastq_sniffer(fh):
    # Strategy:
    #   Ignore up to 5 blank/whitespace-only lines at the beginning of the
    #   file. Read up to 10 records. If at least one record is read (i.e. the
    #   file isn't empty) and the quality scores are in printable ASCII range,
    #   assume the file is FASTQ.
    if _too_many_blanks(fh, 5):
        return False, {}

    try:
        not_empty = False
        for _, seq in zip(range(10), _fastq_to_generator(fh, phred_offset=33)):
            split_length = len(
                (seq.metadata["id"] + seq.metadata["description"]).split(":")
            )
            description = seq.metadata["description"].split(":")
            if split_length == 10 and description[1] in "YN":
                return True, {"variant": "illumina1.8"}
            not_empty = True
        return not_empty, {}
    except (FASTQFormatError, ValueError):
        return False, {}


@fastq.reader(None)
def _fastq_to_generator(
    fh, variant=None, phred_offset=None, constructor=Sequence, **kwargs
):
    # Skip any blank or whitespace-only lines at beginning of file
    try:
        seq_header = next(_line_generator(fh, skip_blanks=True))
    except StopIteration:
        return

    if not seq_header.startswith("@"):
        raise FASTQFormatError(
            "Expected sequence (@) header line at start of file: %r" % str(seq_header)
        )

    while seq_header is not None:
        id_, desc = _parse_fasta_like_header(seq_header)
        seq, qual_header = _parse_sequence_data(fh, seq_header)

        if qual_header != "+" and qual_header[1:] != seq_header[1:]:
            raise FASTQFormatError(
                "Sequence (@) and quality (+) header lines do not match: "
                "%r != %r" % (str(seq_header[1:]), str(qual_header[1:]))
            )

        phred_scores, seq_header = _parse_quality_scores(
            fh, len(seq), variant, phred_offset, qual_header
        )
        yield constructor(
            seq,
            metadata={"id": id_, "description": desc},
            positional_metadata={"quality": phred_scores},
            **kwargs,
        )


@fastq.reader(Sequence)
def _fastq_to_sequence(fh, variant=None, phred_offset=None, seq_num=1, **kwargs):
    return _get_nth_sequence(
        _fastq_to_generator(
            fh,
            variant=variant,
            phred_offset=phred_offset,
            constructor=Sequence,
            **kwargs,
        ),
        seq_num,
    )


@fastq.reader(DNA)
def _fastq_to_dna(fh, variant=None, phred_offset=None, seq_num=1, **kwargs):
    return _get_nth_sequence(
        _fastq_to_generator(
            fh, variant=variant, phred_offset=phred_offset, constructor=DNA, **kwargs
        ),
        seq_num,
    )


@fastq.reader(RNA)
def _fastq_to_rna(fh, variant=None, phred_offset=None, seq_num=1, **kwargs):
    return _get_nth_sequence(
        _fastq_to_generator(
            fh, variant=variant, phred_offset=phred_offset, constructor=RNA, **kwargs
        ),
        seq_num,
    )


@fastq.reader(Protein)
def _fastq_to_protein(fh, variant=None, phred_offset=None, seq_num=1, **kwargs):
    return _get_nth_sequence(
        _fastq_to_generator(
            fh,
            variant=variant,
            phred_offset=phred_offset,
            constructor=Protein,
            **kwargs,
        ),
        seq_num,
    )


@fastq.reader(TabularMSA)
def _fastq_to_tabular_msa(
    fh, variant=None, phred_offset=None, constructor=None, **kwargs
):
    if constructor is None:
        raise ValueError("Must provide `constructor`.")

    return TabularMSA(
        _fastq_to_generator(
            fh,
            variant=variant,
            phred_offset=phred_offset,
            constructor=constructor,
            **kwargs,
        )
    )


@fastq.writer(None)
def _generator_to_fastq(
    obj,
    fh,
    variant=None,
    phred_offset=None,
    id_whitespace_replacement="_",
    description_newline_replacement=" ",
    lowercase=None,
):
    formatted_records = _format_fasta_like_records(
        obj,
        id_whitespace_replacement,
        description_newline_replacement,
        True,
        lowercase=lowercase,
    )
    for header, seq_str, qual_scores in formatted_records:
        qual_str = _encode_phred_to_qual(
            qual_scores, variant=variant, phred_offset=phred_offset
        )
        fh.write("@")
        fh.write(header)
        fh.write("\n")
        fh.write(seq_str)
        fh.write("\n+\n")
        fh.write(qual_str)
        fh.write("\n")


@fastq.writer(Sequence)
def _sequence_to_fastq(
    obj,
    fh,
    variant=None,
    phred_offset=None,
    id_whitespace_replacement="_",
    description_newline_replacement=" ",
    lowercase=None,
):
    _sequences_to_fastq(
        [obj],
        fh,
        variant,
        phred_offset,
        id_whitespace_replacement,
        description_newline_replacement,
        lowercase=lowercase,
    )


@fastq.writer(DNA)
def _dna_to_fastq(
    obj,
    fh,
    variant=None,
    phred_offset=None,
    id_whitespace_replacement="_",
    description_newline_replacement=" ",
    lowercase=None,
):
    _sequences_to_fastq(
        [obj],
        fh,
        variant,
        phred_offset,
        id_whitespace_replacement,
        description_newline_replacement,
        lowercase=lowercase,
    )


@fastq.writer(RNA)
def _rna_to_fastq(
    obj,
    fh,
    variant=None,
    phred_offset=None,
    id_whitespace_replacement="_",
    description_newline_replacement=" ",
    lowercase=None,
):
    _sequences_to_fastq(
        [obj],
        fh,
        variant,
        phred_offset,
        id_whitespace_replacement,
        description_newline_replacement,
        lowercase=lowercase,
    )


@fastq.writer(Protein)
def _protein_to_fastq(
    obj,
    fh,
    variant=None,
    phred_offset=None,
    id_whitespace_replacement="_",
    description_newline_replacement=" ",
    lowercase=None,
):
    _sequences_to_fastq(
        [obj],
        fh,
        variant,
        phred_offset,
        id_whitespace_replacement,
        description_newline_replacement,
        lowercase=lowercase,
    )


@fastq.writer(TabularMSA)
def _tabular_msa_to_fastq(
    obj,
    fh,
    variant=None,
    phred_offset=None,
    id_whitespace_replacement="_",
    description_newline_replacement=" ",
    lowercase=None,
):
    _sequences_to_fastq(
        obj,
        fh,
        variant,
        phred_offset,
        id_whitespace_replacement,
        description_newline_replacement,
        lowercase=lowercase,
    )


def _blank_error(unique_text):
    error_string = ("Found blank or whitespace-only line {} in " "FASTQ file").format(
        unique_text
    )
    raise FASTQFormatError(error_string)


def _parse_sequence_data(fh, prev):
    seq_chunks = []
    for chunk in _line_generator(fh, skip_blanks=False):
        if chunk.startswith("+"):
            if not prev:
                _blank_error("before '+'")
            if not seq_chunks:
                raise FASTQFormatError("Found FASTQ record without sequence data.")
            return "".join(seq_chunks), chunk
        elif chunk.startswith("@"):
            raise FASTQFormatError(
                "Found FASTQ record that is missing a quality (+) header line "
                "after sequence data."
            )
        else:
            if not prev:
                _blank_error("after header or within sequence")
            if _whitespace_regex.search(chunk):
                raise FASTQFormatError(
                    "Found whitespace in sequence data: %r" % str(chunk)
                )
            seq_chunks.append(chunk)
        prev = chunk

    raise FASTQFormatError("Found incomplete/truncated FASTQ record at end of file.")


def _parse_quality_scores(fh, seq_len, variant, phred_offset, prev):
    phred_scores = []
    qual_len = 0
    for chunk in _line_generator(fh, skip_blanks=False):
        if chunk:
            if chunk.startswith("@") and qual_len == seq_len:
                return np.hstack(phred_scores), chunk
            else:
                if not prev:
                    _blank_error("after '+' or within quality scores")
                qual_len += len(chunk)

                if qual_len > seq_len:
                    raise FASTQFormatError(
                        "Found more quality score characters than sequence "
                        "characters. Extra quality score characters: %r"
                        % chunk[-(qual_len - seq_len) :]
                    )

                phred_scores.append(
                    _decode_qual_to_phred(
                        chunk, variant=variant, phred_offset=phred_offset
                    )
                )
        prev = chunk

    if qual_len != seq_len:
        raise FASTQFormatError(
            "Found incomplete/truncated FASTQ record at end of file."
        )
    return np.hstack(phred_scores), None


def _sequences_to_fastq(
    obj,
    fh,
    variant,
    phred_offset,
    id_whitespace_replacement,
    description_newline_replacement,
    lowercase=None,
):
    def seq_gen():
        yield from obj

    _generator_to_fastq(
        seq_gen(),
        fh,
        variant=variant,
        phred_offset=phred_offset,
        id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        lowercase=lowercase,
    )
