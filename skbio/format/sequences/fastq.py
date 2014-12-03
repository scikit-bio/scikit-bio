# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import warnings


def _phred_to_ascii(a, offset):
    """Convert Phred quality score to ASCII character with specified offset"""
    return (a + offset).tostring()


def _phred_to_ascii33(a):
    """Convert Phred quality score to ASCII character with offset of 33"""
    return _phred_to_ascii(a, 33)


def _phred_to_ascii64(a):
    """Convert Phred quality score to ASCII character with offset of 64"""
    return _phred_to_ascii(a, 64)


def format_fastq_record(seqid, seq, qual, phred_offset=33):
    """Format a FASTQ record

    .. note:: Deprecated in scikit-bio 0.2.0-dev
       ``format_fastq_record`` will be removed in scikit-bio 0.3.0. It is
       replaced by ``write``, which is a more general method for serializing
       FASTQ-formatted files. ``write`` supports multiple file formats by
       taking advantage of scikit-bio's I/O registry system. See
       :mod:`skbio.io` for more details.

    Parameters
    ----------
    seqid : bytes
        The sequence ID
    seq : bytes or subclass of BiologicalSequence
        The sequence
    qual : np.array of int8
        The quality scores
    phred_offset : int, either 33 or 64
        Set a phred offset

    Returns
    -------
    bytes : a string representation of a single FASTQ record

    Examples
    --------
    >>> from skbio.format.sequences import format_fastq_record
    >>> from numpy import array, int8
    >>> seqid = 'seq1'
    >>> seq = 'AATTGG'
    >>> qual = array([38, 38, 39, 39, 40, 40], dtype=int8)
    >>> print format_fastq_record(seqid, seq, qual),
    @seq1
    AATTGG
    +
    GGHHII

    """
    warnings.warn(
        "`format_fastq_record` is deprecated and will be removed in "
        "scikit-bio 0.3.0. Please update your code to use `skbio.io.write`.",
        DeprecationWarning)

    if phred_offset == 33:
        phred_f = _phred_to_ascii33
    elif phred_offset == 64:
        phred_f = _phred_to_ascii64
    else:
        raise ValueError("Unknown phred offset: %d" % phred_offset)

    return b'\n'.join([b"@" + seqid, seq, b'+', phred_f(qual), b''])
