r"""
FASTQ format :mod:`skbio.io.fastq`
==================================

.. currentmodule:: skbio.io.fastq

The FASTQ file stores biological sequences,in a plain text format.
Each entry in a FASTQ file contains a biological sequence, typically
DNA or RNA, and a quality score for each nucleotide. This format was
first specified in [1]_.  More information about this format can be
found at [2]_.

An example of a FASTQ-formatted file containing one DNA sequence and
their quality scores is shown below.

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
|Yes       |No        |generator of :mod:`skbio.sequence.BiologicalSequence` |
|          |          |objects                                               |
+----------+----------+------------------------------------------------------+

Format Specification
--------------------
A FASTQ file contains one or more biological sequences with their corresponding
quality score. Each entry in the FASTQ file contains a sequence of nucleotides
a set of quality scores and optionally a description and a id.

Sequence Header
^^^^^^^^^^^^^^^
Each sequence header consists of a single line beginning with a (``@``) symbol.
Immediately after this symbol is a sequence identifier (ID) and an description
separated by one or more whitespace characters.

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
...        "@GAPC_0017:6:1:1259:10413#0/1\n"
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
.. [1] http://nar.oxfordjournals.org/content/38/6/1767
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
from skbio.io import (register_reader, register_sniffer,
                      FASTQFormatError)
from skbio.sequence import (BiologicalSequence)
from skbio.io.util import open_file
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


def _drop_id_marker(s):
    """Drop the first character and decode bytes to text"""
    if s[0] not in ["+", "@"]:
        raise FASTQFormatError("Bad FASTQ format for seq id: %s. " % s)

    id_ = s[1:]
    try:
        return str(id_.decode('utf-8'))
    except AttributeError:
        return id_


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
        return not_empty, {'phred_offset':33}
    except FASTQFormatError:
        try:
            fh.seek(0)
            not_empty = False
            gen = _fastq_to_generator(fh, phred_offset=64)
            for _ in zip(range(10), gen):
                not_empty = True
            return not_empty, {'phred_offset':64}
        except FASTQFormatError:
            return False, {}


@register_reader('fastq')
def _fastq_to_generator(fh, strict=False, enforce_qual_range=True,
                        phred_offset=33):
    r"""yields label, seq, and qual from a fastq file.

    Parameters
    ----------
    fh : open file object or str
        An open fastq file (opened in binary mode) or a path to it.
    strict : bool, optional
        Defaults to ``False``. If strict is true a FastqParse error will be
        raised if the seq and qual labels dont' match.
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
        header = seqid.strip()
        # If the file simply ended in a blankline, do not error
        if header is '':
            continue

        # Error if an incomplete record is found
        # Note: seqid cannot be None, because if all 4 values were None,
        # then the loop condition would be false, and we could not have
        # gotten to this point
        if seq is None or qualid is None or qual is None:
            raise FASTQFormatError("Incomplete FASTQ record found at end "
                                   "of file")

        seq = seq.strip()
        qualid = qualid.strip()
        qual = qual.strip()

        header = _drop_id_marker(header)
        seqid, description = _split_id(header)
        try:
            seq = str(seq.decode("utf-8"))
        except AttributeError:
            pass

        qualid = _drop_id_marker(qualid)
        qualid, _ = _split_id(qualid)
        if strict:
            if seqid != qualid:
                raise FASTQFormatError('ID mismatch: {} != {}'.format(
                    seqid, qualid))

        # bounds based on illumina limits, see:
        # http://nar.oxfordjournals.org/content/38/6/1767/T1.expansion.html
        qual = phred_f(qual)

        if enforce_qual_range and ((qual < 0).any() or (qual > 62).any()):
            raise FASTQFormatError("Failed qual conversion for seq id: %s."
                                   "This may be because you passed an "
                                   "incorrect value for phred_offset" %
                                   seqid)

        yield BiologicalSequence(seq, id=seqid, quality=qual,
                                 description=description)


@register_writer('fastq')
def _generator_to_fastq(obj, fh, id_whitespace_replacement='_',
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

        seq_str = seq.sequence
        qual_str = ''.join(map(chr,seq.quality))
        
        qual_str
        fh.write('@%s\n%s\n+%s\n%s\n' % (header, seq_str, header, qual_str))
