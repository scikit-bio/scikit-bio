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
http://www.phylo.org/tools/obsolete/phylip.html

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import range

from skbio.alignment import Alignment
from skbio.io import register_writer, PhylipFormatError


@register_writer('phylip', Alignment)
def _alignment_to_phylip(obj, fh):
    if not obj.is_valid():
        raise PhylipFormatError(
            "Alignment can only be written in PHYLIP format if all sequences "
            "are of equal length and contain only valid characters within "
            "their character sets.")

    if obj.is_empty():
        raise PhylipFormatError(
            "Alignment can only be written in PHYLIP format if there is at "
            "least one sequence in the alignment.")

    sequence_length = obj.sequence_length()
    if sequence_length == 0:
        raise PhylipFormatError(
            "Alignment can only be written in PHYLIP format if there is at "
            "least one position in the alignment.")

    chunk_size = 10
    for id_ in obj.ids():
        if len(id_) > chunk_size:
            raise PhylipFormatError(
                "Alignment can only be written in PHYLIP format if all "
                "sequence IDs have %d or fewer characters. Found sequence "
                "with ID '%s' that exceeds this limit." % (chunk_size, id_))

    sequence_count = obj.sequence_count()
    fh.write('{0:d} {1:d}\n'.format(sequence_count, sequence_length))

    fmt = '{0:%d}{1}\n' % chunk_size
    for seq in obj:
        chunked_seq = _chunk_str(str(seq), chunk_size)
        fh.write(fmt.format(seq.id, chunked_seq))


def _chunk_str(s, n):
    """Insert a space every `n` characters in `s`."""
    # Modified from http://stackoverflow.com/a/312464/3776794
    return ' '.join((s[i:i+n] for i in range(0, len(s), n)))
