# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function
from future.standard_library import hooks
with hooks():
    from itertools import zip_longest

import warnings

import numpy as np

from skbio.io.util import open_file
from ._exception import FastqParseError


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
    id_ = s[1:]
    try:
        return str(id_.decode('utf-8'))
    except AttributeError:
        return id_


def parse_fastq(data, strict=False, enforce_qual_range=True, phred_offset=33):
    r"""yields label, seq, and qual from a fastq file.

    .. note:: Deprecated in scikit-bio 0.2.0-dev
       ``parse_fastq`` will be removed in scikit-bio 0.3.0. It is replaced by
       ``read``, which is a more general method for deserializing
       FASTQ-formatted files. ``read`` supports multiple file formats,
       automatic file format detection, etc. by taking advantage of
       scikit-bio's I/O registry system. See :mod:`skbio.io` for more details.

    Parameters
    ----------
    data : open file object or str
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
    >>> for label, seq, qual in parse_fastq(fastq_f, phred_offset=64):
    ...     print(label)
    ...     print(seq)
    ...     print(qual)
    seq1
    AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
    [32 32 32 32 25 30 20 29 32 29 35 30 35 33 34 35 33 35 35 32 30 12 34 30 35
     35 25 20 28 20 28 25 28 23  6]
    seq2
    TATGTATATATAACATATACATATATACATACATA
    [29 11 26 27 16 25 29 31 27 25 25 30 32 32 32 33 35 30 28 28 32 34 20 32 32
     35 32 28 33 20 32 32 34 34 34]
    """
    warnings.warn(
        "`parse_fastq` is deprecated and will be removed in scikit-bio 0.3.0. "
        "Please update your code to use `skbio.io.read(fh, format='fastq')` "
        "to obtain a generator of `BiologicalSequence` objects (or "
        "subclasses, see the `constructor` parameter).", UserWarning)

    if phred_offset == 33:
        phred_f = ascii_to_phred33
    elif phred_offset == 64:
        phred_f = ascii_to_phred64
    else:
        raise ValueError("Unknown PHRED offset of %s" % phred_offset)

    with open_file(data, 'rb') as data:
        iters = [iter(data)] * 4
        for seqid, seq, qualid, qual in zip_longest(*iters):
            seqid = seqid.strip()
            # If the file simply ended in a blankline, do not error
            if seqid is '':
                continue
            # Error if an incomplete record is found
            # Note: seqid cannot be None, because if all 4 values were None,
            # then the loop condition would be false, and we could not have
            # gotten to this point
            if seq is None or qualid is None or qual is None:
                raise FastqParseError("Incomplete FASTQ record found at end "
                                      "of file")

            seq = seq.strip()
            qualid = qualid.strip()
            qual = qual.strip()

            seqid = _drop_id_marker(seqid)

            try:
                seq = str(seq.decode("utf-8"))
            except AttributeError:
                pass

            qualid = _drop_id_marker(qualid)
            if strict:
                if seqid != qualid:
                    raise FastqParseError('ID mismatch: {} != {}'.format(
                        seqid, qualid))

            # bounds based on illumina limits, see:
            # http://nar.oxfordjournals.org/content/38/6/1767/T1.expansion.html
            qual = phred_f(qual)
            if enforce_qual_range and ((qual < 0).any() or (qual > 62).any()):
                raise FastqParseError("Failed qual conversion for seq id: %s. "
                                      "This may be because you passed an "
                                      "incorrect value for phred_offset." %
                                      seqid)

            yield (seqid, seq, qual)
