"""
PHYLIP multiple sequence alignment format (:mod:`skbio.io.phylip`)
==================================================================

.. currentmodule:: skbio.io.phylip

TODO add description

Format Specification
--------------------

TODO add format spec
TODO add references

http://evolution.genetics.washington.edu/phylip/doc/sequence.html
http://www.bioperl.org/wiki/PHYLIP_multiple_alignment_format

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.alignment import Alignment, SequenceCollectionError
from skbio.io import register_writer


@register_writer('phylip', Alignment)
def _alignment_to_phylip(obj, fh):
    if not obj._validate_lengths():
        raise SequenceCollectionError("PHYLIP-formatted string can only "
                                      "be generated if all sequences are "
                                      "of equal length.")

    if obj.is_empty():
        raise SequenceCollectionError("PHYLIP-formatted string can only "
                                      "be generated if there is at least "
                                      "one sequence in the Alignment.")

    sequence_length = obj.sequence_length()
    if sequence_length == 0:
        raise SequenceCollectionError("PHYLIP-formatted string can only "
                                      "be generated if there is at least "
                                      "one position in the Alignment.")

    sequence_count = obj.sequence_count()
    fh.write("%d %d\n" % (sequence_count, sequence_length))

    for seq in obj:
        fh.write("%s %s\n" % (seq.id, str(seq)))
