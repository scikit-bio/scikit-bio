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
from skbio.alignment import SequenceCollection, Alignment
from skbio.sequence import (BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence)

from skbio.util import cardinal_to_ordinal


def _ascii_to_phred(s, offset):
    """Convert ascii to Phred quality score with specified ASCII offset."""
    return np.fromstring(s, dtype='|S1').view(np.int8) - offset


def ascii_to_phred33(s):
    """Convert ascii string to Phred quality score with ASCII offset of 33.

    Standard "Sanger" ASCII offset of 33. This is used by Illumina in CASAVA
    versions after 1.8.0, and most other places. Note that internal Illumina
    files still use offset of 64
    """
    return _ascii_to_phred(s, 33)


def ascii_to_phred64(s):
    """Convert ascii string to Phred quality score with ASCII offset of 64.

    Illumina-specific ASCII offset of 64. This is used by Illumina in CASAVA
    versions prior to 1.8.0, and in Illumina internal formats (e.g.,
    export.txt files).
    """
    return _ascii_to_phred(s, 64)


def _drop_id_marker(s, marker):
    """Drop the first character"""
    if s[0] != marker:
        raise FASTQFormatError(
            "Expected header line to start with %r character: %r" %
            (marker, s))
    return s[1:]


def _split_id(s):
    """Split up line into id and description"""
    toks = re.split("\s+", s)
    seqid, description = '', ''
    if len(toks) > 0:
        seqid = toks[0]
    if len(toks) > 1:
        description = ' '.join(toks[1:])
    return seqid, description


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
def _fastq_to_generator(fh, enforce_qual_range=True,
                        phred_offset=33, constructor=BiologicalSequence):
    r"""yields label, seq, and qual from a fastq file.

    Parameters
    ----------
    fh : open file object or str
        An open fastq file (opened in binary mode) or a path to it.
    enforce_qual_range : bool, optional
        Defaults to ``True``. If ``True``, an exception will be raised if a
        quality score outside the range [0, 62] is detected
    phred_offset : {33, 64}, optional
        What Phred offset to use when converting qual score symbols to integers

    Returns
    -------
    label, seq, qual : (str, bytes, np.array)
        yields the label, sequence and quality for each entry

    Examples
    --------
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
    """
    if phred_offset == 33:
        phred_f = ascii_to_phred33
    elif phred_offset == 64:
        phred_f = ascii_to_phred64
    else:
        raise ValueError("Unknown PHRED offset of %s" % phred_offset)

    iters = [iter(fh)] * 4
    for seqid, seq, qualid, qual in zip_longest(*iters):
        # Error if an incomplete record is found
        # Note: seqid cannot be None, because if all 4 values were None,
        # then the loop condition would be false, and we could not have
        # gotten to this point
        if seq is None or qualid is None or qual is None:
            raise FASTQFormatError(
                "Found incomplete/truncated FASTQ record at end of file.")

        if not seqid.strip() or not seq.strip() or not qualid.strip() or \
                not qual.strip():
            raise FASTQFormatError(
                "Found blank or whitespace-only line in FASTQ-formatted file.")

        header = seqid.strip()
        # If the file simply ended in a blankline, do not error
        if header is '':
            continue

        seq = seq.strip()
        for c in seq:
            if c.isspace():
                raise FASTQFormatError(
                    "Found whitespace in sequence data: %r" % seq)

        qualid = qualid.strip()

        header = _drop_id_marker(header, '@')
        qualid = _drop_id_marker(qualid, '+')

        if qualid and header != qualid:
            raise FASTQFormatError(
                "Sequence (@) and quality (+) header lines do not match: "
                "%r != %r" % (header, qualid))

        seqid, description = _split_id(header)
        qualid, _ = _split_id(qualid)

        qual = qual.strip()
        for c in qual:
            ascii_code = ord(c)
            if ascii_code < 33 or ascii_code > 126:
                raise FASTQFormatError(
                    "Found quality score encoded as a non-printable ASCII "
                    "character (ASCII code %d). Each quality score must be in "
                    "the ASCII code range 33-126 (inclusive), regardless of "
                    "FASTQ variant used." % ascii_code)

        # bounds based on illumina limits, see:
        # http://nar.oxfordjournals.org/content/38/6/1767/T1.expansion.html
        qual = phred_f(qual)

        if enforce_qual_range and ((qual < 0).any() or (qual > 62).any()):
            raise FASTQFormatError("Failed qual conversion for seq id: %s."
                                   "This may be because you passed an "
                                   "incorrect value for phred_offset" %
                                   seqid)

        yield constructor(seq, id=seqid, quality=qual,
                          description=description)


@register_reader('fastq', BiologicalSequence)
def _fastq_to_biological_sequence(fh, seq_num=1, phred_offset=33):
    return _fastq_to_sequence(fh, seq_num, BiologicalSequence,
                              phred_offset=phred_offset)


@register_reader('fastq', NucleotideSequence)
def _fastq_to_nucleotide_sequence(fh, seq_num=1, phred_offset=33):
    return _fastq_to_sequence(fh, seq_num, NucleotideSequence,
                              phred_offset=phred_offset)


@register_reader('fastq', DNASequence)
def _fastq_to_dna_sequence(fh, seq_num=1, phred_offset=33):
    return _fastq_to_sequence(fh, seq_num, DNASequence,
                              phred_offset=phred_offset)


@register_reader('fastq', RNASequence)
def _fastq_to_rna_sequence(fh, seq_num=1, phred_offset=33):
    return _fastq_to_sequence(fh, seq_num, RNASequence,
                              phred_offset=phred_offset)


@register_reader('fastq', ProteinSequence)
def _fastq_to_protein_sequence(fh, seq_num=1,
                               phred_offset=33):
    return _fastq_to_sequence(fh, seq_num, ProteinSequence,
                              phred_offset=phred_offset)


@register_reader('fastq', SequenceCollection)
def _fastq_to_sequence_collection(fh, constructor=BiologicalSequence,
                                  phred_offset=33):
    return SequenceCollection(
        list(_fastq_to_generator(fh, constructor=constructor,
                                 phred_offset=phred_offset)))


@register_reader('fastq', Alignment)
def _fastq_to_alignment(fh, constructor=BiologicalSequence,
                        phred_offset=33):
    return Alignment(
        list(_fastq_to_generator(fh, constructor=constructor,
                                 phred_offset=phred_offset)))


@register_writer('fastq')
def _generator_to_fastq(obj, fh, id_whitespace_replacement='_',
                        description_newline_replacement=' ',
                        enforce_qual_range=True,
                        phred_offset=33):
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

        seq_str = seq.sequence
        if phred_offset not in [33, 64]:
            raise ValueError("Unknown PHRED offset of %s" % phred_offset)

        if enforce_qual_range and ((seq.quality < 0).any() or
                                   (seq.quality > 62).any()):
            raise FASTQFormatError("Failed qual conversion for seq id: %s."
                                   "This may be because you passed an "
                                   "incorrect value for phred_offset" %
                                   seq.id)

        qual = seq.quality+phred_offset
        qual_str = ''.join(map(chr, qual))

        fh.write('@%s\n%s\n+%s\n%s\n' % (header, seq_str, header, qual_str))


@register_writer('fastq', BiologicalSequence)
def _biological_sequence_to_fastq(obj, fh, id_whitespace_replacement='_',
                                  description_newline_replacement=' ',
                                  phred_offset=33):

    _sequence_to_fastq(obj, fh, id_whitespace_replacement,
                       description_newline_replacement,
                       phred_offset=phred_offset)


@register_writer('fastq', NucleotideSequence)
def _nucleotide_sequence_to_fastq(obj, fh, id_whitespace_replacement='_',
                                  description_newline_replacement=' ',
                                  phred_offset=33):
    _sequence_to_fastq(obj, fh, id_whitespace_replacement,
                       description_newline_replacement,
                       phred_offset=phred_offset)


@register_writer('fastq', DNASequence)
def _dna_sequence_to_fastq(obj, fh, id_whitespace_replacement='_',
                           description_newline_replacement=' ',
                           phred_offset=33):
    _sequence_to_fastq(obj, fh, id_whitespace_replacement,
                       description_newline_replacement,
                       phred_offset=phred_offset)


@register_writer('fastq', RNASequence)
def _rna_sequence_to_fastq(obj, fh, id_whitespace_replacement='_',
                           description_newline_replacement=' ',
                           phred_offset=33):
    _sequence_to_fastq(obj, fh, id_whitespace_replacement,
                       description_newline_replacement,
                       phred_offset=phred_offset)


@register_writer('fastq', ProteinSequence)
def _protein_sequence_to_fastq(obj, fh, id_whitespace_replacement='_',
                               description_newline_replacement=' ',
                               phred_offset=33):
    _sequence_to_fastq(obj, fh, id_whitespace_replacement,
                       description_newline_replacement,
                       phred_offset=phred_offset)


@register_writer('fastq', SequenceCollection)
def _sequence_collection_to_fastq(obj, fh, id_whitespace_replacement='_',
                                  description_newline_replacement=' ',
                                  phred_offset=33):
    _sequences_to_fastq(obj, fh, id_whitespace_replacement,
                        description_newline_replacement,
                        phred_offset=phred_offset)


@register_writer('fastq', Alignment)
def _alignment_to_fastq(obj, fh, id_whitespace_replacement='_',
                        description_newline_replacement=' ',
                        phred_offset=33):
    _sequences_to_fastq(obj, fh, id_whitespace_replacement,
                        description_newline_replacement,
                        phred_offset=phred_offset)


def _fastq_to_sequence(fh, seq_num, constructor, phred_offset=33):
    if seq_num < 1:
        raise FASTQFormatError(
            "Invalid sequence number (seq_num=%d). seq_num must be between 1 "
            "and the number of sequences in the FASTQ-formatted file "
            "(inclusive)." % seq_num)

    seq_idx = seq_num - 1
    seq = None
    try:
        gen = _fastq_to_generator(fh, constructor=constructor,
                                  phred_offset=phred_offset)
        for idx, curr_seq in enumerate(gen):
            if idx == seq_idx:
                seq = curr_seq
                break
    finally:
        gen.close()

    if seq is None:
        raise FASTQFormatError(
            "Reached end of FASTQ-formatted file before finding %s biological "
            "sequence." % cardinal_to_ordinal(seq_num))
    return seq


def _sequence_to_fastq(obj, fh, id_whitespace_replacement,
                       description_newline_replacement,
                       phred_offset=33):
    def seq_gen():
        yield obj

    _generator_to_fastq(
        seq_gen(), fh, id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        phred_offset=phred_offset)


def _sequences_to_fastq(obj, fh, id_whitespace_replacement,
                        description_newline_replacement,
                        phred_offset=33):
    def seq_gen():
        for seq in obj:
            yield seq

    _generator_to_fastq(
        seq_gen(), fh, id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        phred_offset=phred_offset)
