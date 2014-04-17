#!/usr/bin/env python
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

    seqid : str
        The sequence ID
    seq : str or subclass of BiologicalSequence
        The sequence
    qual : np.array of int8
        The quality scores
    """
    if is_casava_v180_or_later("@%s" % seqid):
        phred_f = _phred_to_ascii33
    else:
        phred_f = _phred_to_ascii64

    return "@%s\n%s\n+\n%s\n" % (seqid, seq, phred_f(qual))
