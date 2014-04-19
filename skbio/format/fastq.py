"""Formatters for FASTQ"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.util.misc import is_casava_v180_or_later


def _phred_to_ascii(a, offset):
    """Convert Phred quality score to ASCII character with specified offset"""
    return (a + offset).tostring()


def _phred_to_ascii33(a):
    """Convert Phred quality score to ASCII character with offset of 33"""
    return _phred_to_ascii(a, 33)


def _phred_to_ascii64(a):
    """Convert Phred quality score to ASCII character with offset of 64"""
    return _phred_to_ascii(a, 64)


def format_fastq_record(seqid, seq, qual):
    """Format a FASTQ record

    Parameters
    ----------

    seqid : bytes
        The sequence ID
    seq : bytes or subclass of BiologicalSequence
        The sequence
    qual : np.array of int8
        The quality scores

    Returns
    -------

    bytes : a string representation of a single FASTQ record

    Examples
    --------
    >>> from skbio.format.fastq import format_fastq_record
    >>> from numpy import array, int8
    >>> seqid = 'seq1'
    >>> seq = 'AATTGG'
    >>> qual = array([38, 38, 39, 39, 40, 40], dtype=int8)
    >>> print format_fastq_record(seqid, seq, qual),
    @seq1
    AATTGG
    +
    ffgghh
    """
    if is_casava_v180_or_later(b"@" + seqid):
        phred_f = _phred_to_ascii33
    else:
        phred_f = _phred_to_ascii64

    return b'\n'.join([b"@" + seqid, seq, b'+', phred_f(qual), b''])
