#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import division

from skbio.core.exception import FastqParseError


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
    label, seq, qual : string
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
    ````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
    seq2
    TATGTATATATAACATATACATATATACATACATA
    ]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb

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
            yield record[0][1:], record[1], record[3]
            line_num = 0
            record = []
        record.append(line.strip())

    if record:
        if strict and record[0]:  # make sure the seq and qual labels match
            if record[0][1:] != record[2][1:]:
                raise FastqParseError('Invalid format: %s -- %s'
                                      % (record[0][1:], record[2][1:]))

        if record[0]:  # could be just an empty line at eof
            yield record[0][1:], record[1], record[3]

    if isinstance(data, file):
        data.close()

