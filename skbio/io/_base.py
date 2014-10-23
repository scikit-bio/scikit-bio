# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import range


def _chunk_str(s, n, char):
    """Insert `char` character every `n` characters in string `s`.

    Canonically pronounced "chunkster".

    """
    # Modified from http://stackoverflow.com/a/312464/3776794
    if n < 1:
        raise ValueError(
            "Cannot split string into chunks with n=%d. n must be >= 1." % n)
    return char.join((s[i:i+n] for i in range(0, len(s), n)))


def _decode_qual_to_phred(qual_str, variant=None, phred_offset=None):
    if variant is None and phred_offset is None:
        raise ValueError(
            "Must provide either `variant` or `phred_offset` in order to "
            "decode quality scores.")
    if variant is not None and phred_offset is not None:
        raise ValueError(
            "Cannot provide both `variant` and `phred_offset`.")

    if variant is not None:
        if variant == 'sanger':
            phred_offset = 33
            phred_range = (0, 93)
        elif variant == 'illumina1.3':
            phred_offset = 64
            phred_range = (0, 62)
        elif variant == 'illumina1.8':
            phred_offset = 33
            phred_range = (0, 62)
        elif variant == 'solexa':
            phred_offset = 64
            phred_range = (-5, 62)
            raise NotImplementedError(
                "Decoding Solexa quality scores is not currently supported, "
                "as quality scores are always stored as Phred scores in "
                "scikit-bio. Please see the following scikit-bio issue to "
                "track progress on this:\n\t"
                "https://github.com/biocore/scikit-bio/issues/719")
        else:
            raise ValueError("Unrecognized variant %r." % variant)
    else:
        # Phred scores cannot be less than zero. There is no theoretical max to
        # a Phred score, but more than a byte of information doesn't make
        # sense.
        phred_range = (0, 255)

    phred = []
    for c in qual_str:
        score = ord(c) - phred_offset
        if phred_range[0] <= score <= phred_range[1]:
            phred.append(score)
        else:
            raise ValueError("Decoded Phred score %d is out of range [%d, %d]."
                             % (score, phred_range[0], phred_range[1]))
    return phred
