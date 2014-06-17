#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

__credits__ = "https://github.com/biocore/scikit-bio/graphs/contributors"
__version__ = "0.1.3-dev"

mottos = [
    # 03/15/2014
    "It's gonna get weird, bro.",
    # 05/14/2014
    "no cog yay"
]
motto = mottos[-1]

title = r"""
               _ _    _ _          _     _
              (_) |  (_) |        | |   (_)
      ___  ___ _| | ___| |_ ______| |__  _  ___
     / __|/ __| | |/ / | __|______| '_ \| |/ _ \
     \__ \ (__| |   <| | |_       | |_) | | (_) |
     |___/\___|_|_|\_\_|\__|      |_.__/|_|\___/

"""

art = r"""

           Opisthokonta
                   \  Amoebozoa
                    \ /
                     *    Euryarchaeota
                      \     |_ Crenarchaeota
                       \   *
                        \ /
                         *
                        /
                       /
                      /
                     *
                    / \
                   /   \
        Proteobacteria  \
                       Cyanobacteria
"""

if __doc__ is None:
    __doc__ = title + art
else:
    __doc__ = title + art + __doc__

# imports included for convenience
from skbio.core.sequence import (
    BiologicalSequence, NucleotideSequence, DNA, DNASequence, RNA, RNASequence,
    Protein, ProteinSequence)
from skbio.core.distance import DistanceMatrix
from skbio.core.alignment import (
    align_striped_smith_waterman, SequenceCollection, Alignment)
from skbio.core.alignment.pairwise import (
    local_pairwise_align_nucleotide, local_pairwise_align_protein)
from skbio.core.tree import (
    TreeNode, nj)
from skbio.parse.sequences import (
    parse_fasta, parse_fastq, parse_qual, FastaIterator, FastqIterator,
    SequenceIterator)

__all__ = ['BiologicalSequence', 'NucleotideSequence', 'DNA', 'DNASequence',
           'RNA', 'RNASequence', 'Protein', 'ProteinSequence',
           'DistanceMatrix', 'align_striped_smith_waterman',
           'SequenceCollection', 'Alignment', 'TreeNode', 'nj', 'parse_fasta',
           'parse_fastq', 'parse_qual', 'FastaIterator', 'FastqIterator',
           'SequenceIterator']

from numpy.testing import Tester
test = Tester().test

if __name__ == '__main__':
    print(title)
    print(art)
