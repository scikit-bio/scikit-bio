# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# ruff: noqa: D104

# Add skbio.io to sys.modules to prevent cycles in our imports
import skbio.io  # noqa

# imports included for convenience
from skbio.sequence import Sequence, DNA, RNA, Protein, GeneticCode, SubstitutionMatrix
from skbio.stats.distance import DistanceMatrix
from skbio.alignment import local_pairwise_align_ssw, TabularMSA
from skbio.tree import TreeNode, nj
from skbio.embedding import ProteinEmbedding
from skbio.io import read, write
from skbio.stats.ordination import OrdinationResults
from skbio.table import Table
from skbio.metadata import SampleMetadata
import skbio.diversity  # noqa
import skbio.stats.evolve  # noqa

__all__ = [
    "Sequence",
    "DNA",
    "RNA",
    "Protein",
    "GeneticCode",
    "SubstitutionMatrix",
    "DistanceMatrix",
    "TabularMSA",
    "local_pairwise_align_ssw",
    "TreeNode",
    "nj",
    "read",
    "write",
    "OrdinationResults",
    "Table",
    "SampleMetadata",
]

__credits__ = "https://github.com/scikit-bio/scikit-bio/graphs/contributors"
__version__ = "0.6.2"


mottos = [
    # 03/15/2014
    "It's gonna get weird, bro.",
    # 05/14/2014
    "no cog yay",
    # 03/18/2015
    "bincount!",
]
motto = mottos[-1]

# Created at patorjk.com

title = r"""
*                                                    *
               _ _    _ _          _     _
              (_) |  (_) |        | |   (_)
      ___  ___ _| | ___| |_ ______| |__  _  ___
     / __|/ __| | |/ / | __|______| '_ \| |/ _ \
     \__ \ (__| |   <| | |_       | |_) | | (_) |
     |___/\___|_|_|\_\_|\__|      |_.__/|_|\___/

*                                                    *
"""

# Created by @gregcaporaso

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

__doc__ = title + art + (__doc__ or "")
