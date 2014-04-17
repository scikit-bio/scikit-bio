# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import division

import numpy as np

from skbio.core.exception import FastqParseError
from skbio.util.misc import is_casava_v180_or_later


def _ascii_to_phred(c, offset):
    """Convert ascii to Phred quality score with specified ASCII offset."""
    return ord(c) - offset


def _ascii_to_phred33(c):
    """Convert ascii character to Phred quality score with ASCII offset of 33.

    Standard "Sanger" ASCII offset of 33. This is used by Illumina in CASAVA
    versions after 1.8.0, and most other places. Note that internal Illumina
    files still use offset of 64
    """
    return _ascii_to_phred(c, 33)


def _ascii_to_phred64(c):
    """Convert ascii character to Phred quality score with ASCII offset of 64.

    Illumina-specific ASCII offset of 64. This is used by Illumina in CASAVA
    versions prior to 1.8.0, and in Illumina internal formats (e.g.,
    export.txt files).
    """
    return _ascii_to_phred(c, 64)


def parse_fastq(data, strict=False):
    r"""yields label, seq, and qual from a fastq file.

    Parameters
    ----------
    data : open file object
        An open fastq file.

    strict : bool
        If strict is true a FastqParse error will be raised if the seq and qual
        labels dont' match.

    Returns
    -------
    label, seq, qual : (str, str, np.array)
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
    >>> for label, seq, qual in parse_fastq(fastq_f):
    ...     print label
    ...     print seq
    ...     print qual
    seq1
    AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
    [32 32 32 32 25 30 20 29 32 29 35 30 35 33 34 35 33 35 35 32 30 12 34 30 35
     35 25 20 28 20 28 25 28 23  6]
    seq2
    TATGTATATATAACATATACATATATACATACATA
    [29 11 26 27 16 25 29 31 27 25 25 30 32 32 32 33 35 30 28 28 32 34 20 32 32
     35 32 28 33 20 32 32 34 34 34]

    """
    # fastq format is very simple, defined by blocks of 4 lines
    line_num = -1
    record = []
    for line in data:
        line_num += 1
        if line_num == 4:
            if strict:  # make sure the seq and qual labels match
                if record[0][1:] != record[2][1:]:
                    raise FastqParseError('Invalid format: %s -- %s'
                                          % (record[0][1:], record[2][1:]))

            if is_casava_v180_or_later('@%s' % record[0]):
                phred_f = _ascii_to_phred33
            else:
                phred_f = _ascii_to_phred64

            qual = np.array([phred_f(q) for q in record[3]], dtype=int)
            yield record[0][1:], record[1], qual
            line_num = 0
            record = []
        record.append(line.strip())

    if record:
        if strict and record[0]:  # make sure the seq and qual labels match
            if record[0][1:] != record[2][1:]:
                raise FastqParseError('Invalid format: %s -- %s'
                                      % (record[0][1:], record[2][1:]))

        if record[0]:  # could be just an empty line at eof
            if is_casava_v180_or_later('@%s' % record[0]):
                phred_f = _ascii_to_phred33
            else:
                phred_f = _ascii_to_phred64

            qual = np.array([phred_f(q) for q in record[3]], dtype=int)
            yield record[0][1:], record[1], qual

    if isinstance(data, file):
        data.close()
