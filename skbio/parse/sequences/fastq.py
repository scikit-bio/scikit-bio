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
from skbio.util.io import open_file


def _ascii_to_phred(s, offset):
    """Convert ascii to Phred quality score with specified ASCII offset."""
    return np.fromstring(s, dtype='|S1').view(np.int8) - offset


def _ascii_to_phred33(s):
    """Convert ascii string to Phred quality score with ASCII offset of 33.

    Standard "Sanger" ASCII offset of 33. This is used by Illumina in CASAVA
    versions after 1.8.0, and most other places. Note that internal Illumina
    files still use offset of 64
    """
    return _ascii_to_phred(s, 33)


def _ascii_to_phred64(s):
    """Convert ascii string to Phred quality score with ASCII offset of 64.

    Illumina-specific ASCII offset of 64. This is used by Illumina in CASAVA
    versions prior to 1.8.0, and in Illumina internal formats (e.g.,
    export.txt files).
    """
    return _ascii_to_phred(s, 64)


def _drop_id_marker(s):
    """Drop the first character and decode bytes to text"""
    id_ = s[1:]
    try:
        return str(id_.decode('utf-8'))
    except AttributeError:
        return id_


def parse_fastq(data, strict=False, phred_offset=33):
    r"""yields label, seq, and qual from a fastq file.

    Parameters
    ----------
    data : open file object or str
        An open fastq file (opened in binary mode) or a path to it.

    strict : bool
        If strict is true a FastqParse error will be raised if the seq and qual
        labels dont' match.

    phred_offset : int or None
        Force a Phred offset, currently restricted to either 33 or 64.
        Default behavior is to infer the Phred offset.

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
    >>> for label, seq, qual in parse_fastq(fastq_f, phred_offset=64):
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
    # line number for modulus operation
    SEQUENCEID = 0
    SEQUENCE = 1
    QUALID = 2
    QUAL = 3

    with open_file(data, 'rb') as data:
        data = iter(data)
        first_line = next(data).strip()

        if phred_offset == 33:
            phred_f = _ascii_to_phred33
        elif phred_offset == 64:
            phred_f = _ascii_to_phred64
        else:
            raise ValueError("Unknown PHRED offset of %s" % phred_offset)

        seqid = _drop_id_marker(first_line)
        seq = None
        qualid = None
        qual = None

        for idx, line in enumerate(data):
            # +1 due to fetch of line prior to loop
            lineno = idx + 1
            linetype = lineno % 4
            line = line.strip()

            if linetype == SEQUENCEID:
                yield seqid, seq, qual

                seqid = _drop_id_marker(line)
                seq = None
                qualid = None
                qual = None
            elif linetype == SEQUENCE:
                seq = line
                try:
                    seq = str(seq.decode("utf-8"))
                except:
                    pass
            elif linetype == QUALID:
                qualid = _drop_id_marker(line)
                if strict:
                    if seqid != qualid:
                        raise FastqParseError('ID mismatch: {} != {}'.format(
                            seqid, qualid))
            # bounds based on illumina limits, see:
            # http://nar.oxfordjournals.org/content/38/6/1767/T1.expansion.html
            elif linetype == QUAL:
                qual = phred_f(line)
                if (qual < 0).any() or (qual > 62).any():
                    raise FastqParseError("Failed qual conversion for seq "
                                          "id: %s. This may be because you "
                                          "passed an incorrect value for "
                                          "phred_offset." % seqid)
        if seqid:
            yield (seqid, seq, qual)
