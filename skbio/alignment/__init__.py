r"""Sequence Alignments (:mod:`skbio.alignment`)
============================================

.. currentmodule:: skbio.alignment

This module provides functionality for computing and manipulating sequence
alignments. DNA, RNA, and protein sequences can be aligned, as well as
sequences with custom alphabets.


Alignment structure
-------------------

.. autosummary::
   :toctree: generated/

   TabularMSA
   AlignPath
   PairAlignPath


Alignment algorithms
--------------------

.. rubric:: Optimized (i.e., production-ready) algorithms

.. autosummary::
   :toctree: generated/

   StripedSmithWaterman
   AlignmentStructure
   local_pairwise_align_ssw

.. rubric:: Slow (i.e., educational-purposes only) algorithms

.. autosummary::
   :toctree: generated/

   global_pairwise_align_nucleotide
   global_pairwise_align_protein
   global_pairwise_align
   local_pairwise_align_nucleotide
   local_pairwise_align_protein
   local_pairwise_align

Deprecated functionality
^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

    make_identity_substitution_matrix


Tutorial
--------

Alignment data structure
^^^^^^^^^^^^^^^^^^^^^^^^

Load two DNA sequences that have been previously aligned into a ``TabularMSA``
object, using sequence IDs as the MSA's index:

>>> from skbio import TabularMSA, DNA
>>> seqs = [DNA("ACC--G-GGTA..", metadata={'id': "seq1"}),
...         DNA("TCC--G-GGCA..", metadata={'id': "seq2"})]
>>> msa = TabularMSA(seqs, minter='id')
>>> msa
TabularMSA[DNA]
----------------------
Stats:
    sequence count: 2
    position count: 13
----------------------
ACC--G-GGTA..
TCC--G-GGCA..
>>> msa.index
Index(['seq1', 'seq2'], dtype='object')

Using the optimized alignment algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using the convenient ``local_pairwise_align_ssw`` function:

>>> from skbio.alignment import local_pairwise_align_ssw
>>> alignment, score, start_end_positions = local_pairwise_align_ssw(
...     DNA("ACTAAGGCTCTCTACCCCTCTCAGAGA"),
...     DNA("ACTAAGGCTCCTAACCCCCTTTTCTCAGA")
... )
>>> alignment
TabularMSA[DNA]
------------------------------
Stats:
    sequence count: 2
    position count: 30
------------------------------
ACTAAGGCTCTCT-ACCCC----TCTCAGA
ACTAAGGCTC-CTAACCCCCTTTTCTCAGA
>>> score
27
>>> start_end_positions
[(0, 24), (0, 28)]

Using the ``StripedSmithWaterman`` object:

>>> from skbio.alignment import StripedSmithWaterman
>>> query = StripedSmithWaterman("ACTAAGGCTCTCTACCCCTCTCAGAGA")
>>> alignment = query("AAAAAACTCTCTAAACTCACTAAGGCTCTCTACCCCTCTTCAGAGAAGTCGA")
>>> print(alignment)
ACTAAGGCTC...
ACTAAGGCTC...
Score: 49
Length: 28

Using the ``StripedSmithWaterman`` object for multiple targets in an efficient
way and finding the aligned sequence representations:

>>> from skbio.alignment import StripedSmithWaterman
>>> alignments = []
>>> target_sequences = [
...     "GCTAACTAGGCTCCCTTCTACCCCTCTCAGAGA",
...     "GCCCAGTAGCTTCCCAATATGAGAGCATCAATTGTAGATCGGGCC",
...     "TCTATAAGATTCCGCATGCGTTACTTATAAGATGTCTCAACGG",
...     "TAGAGATTAATTGCCACTGCCAAAATTCTG"
... ]
>>> query_sequence = "ACTAAGGCTCTCTACCCCTCTCAGAGA"
>>> query = StripedSmithWaterman(query_sequence)
>>> for target_sequence in target_sequences:
...     alignment = query(target_sequence)
...     alignments.append(alignment)
...
>>> print(alignments[0])
ACTAAGGCTC...
ACT-AGGCTC...
Score: 38
Length: 30
>>> print(alignments[0].aligned_query_sequence)
ACTAAGGCTC---TCTACCCCTCTCAGAGA
>>> print(alignments[0].aligned_target_sequence)
ACT-AGGCTCCCTTCTACCCCTCTCAGAGA

Using the slow alignment algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

scikit-bio also provides pure-Python implementations of Smith-Waterman and
Needleman-Wunsch alignment. These are much slower than the methods described
above, but serve as useful educational examples as they're simpler to
experiment with. Functions are provided for local and global alignment of
protein and nucleotide sequences. The ``global*`` and ``local*`` functions
differ in the underlying algorithm that is applied (``global*`` uses Needleman-
Wunsch while ``local*`` uses Smith-Waterman), and ``*protein`` and
``*nucleotide`` differ in their default scoring of matches, mismatches, and
gaps.

Here we locally align a pair of protein sequences using gap open penalty
of 11 and a gap extend penalty of 1 (in other words, it is much more
costly to open a new gap than extend an existing one).

>>> from skbio import Protein
>>> from skbio.alignment import local_pairwise_align_protein
>>> s1 = Protein("HEAGAWGHEE")
>>> s2 = Protein("PAWHEAE")
>>> alignment, score, start_end_positions = local_pairwise_align_protein(
...     s1, s2, 11, 1)

This returns an ``skbio.TabularMSA`` object, the alignment score, and start/end
positions of each aligned sequence:

>>> alignment
TabularMSA[Protein]
---------------------
Stats:
    sequence count: 2
    position count: 5
---------------------
AWGHE
AW-HE
>>> score
25.0
>>> start_end_positions
[(4, 8), (1, 4)]

Similarly, we can perform global alignment of nucleotide sequences:

>>> from skbio import DNA
>>> from skbio.alignment import global_pairwise_align_nucleotide
>>> s1 = DNA("GCGTGCCTAAGGTATGCAAG")
>>> s2 = DNA("ACGTGCCTAGGTACGCAAG")
>>> alignment, score, start_end_positions = global_pairwise_align_nucleotide(
...     s1, s2)
>>> alignment
TabularMSA[DNA]
----------------------
Stats:
    sequence count: 2
    position count: 20
----------------------
GCGTGCCTAAGGTATGCAAG
ACGTGCCTA-GGTACGCAAG

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._tabular_msa import TabularMSA
from ._pairwise import (
    local_pairwise_align_nucleotide,
    local_pairwise_align_protein,
    local_pairwise_align,
    global_pairwise_align_nucleotide,
    global_pairwise_align_protein,
    global_pairwise_align,
    make_identity_substitution_matrix,
    local_pairwise_align_ssw,
)
from skbio.alignment._ssw_wrapper import StripedSmithWaterman, AlignmentStructure
from skbio.alignment._path import AlignPath, PairAlignPath

__all__ = [
    "TabularMSA",
    "StripedSmithWaterman",
    "AlignmentStructure",
    "AlignPath",
    "PairAlignPath",
    "local_pairwise_align_ssw",
    "global_pairwise_align",
    "global_pairwise_align_nucleotide",
    "global_pairwise_align_protein",
    "local_pairwise_align",
    "local_pairwise_align_nucleotide",
    "local_pairwise_align_protein",
    "make_identity_substitution_matrix",
]
