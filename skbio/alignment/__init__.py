r"""Sequence Alignments (:mod:`skbio.alignment`)
============================================

.. currentmodule:: skbio.alignment

This module provides functionality for computing and manipulating sequence
alignments. DNA, RNA, and protein sequences can be aligned, as well as
sequences with or without custom alphabets.


Alignment structures
--------------------

Tabular format of aligned sequences.

.. autosummary::
   :toctree: generated/

   TabularMSA

Compact format of alignment paths.

.. autosummary::
   :toctree: generated/

   AlignPath
   PairAlignPath


Pairwise alignment
------------------

A versatile, efficient and generalizable pairwise alignment algorithm.

.. autosummary::
   :toctree: generated/

    pair_align

Convenience wrappers with preset scoring schemes.

.. autosummary::
   :toctree: generated/

    pair_align_nucl
    pair_align_prot


Alignment statistics
--------------------

.. autosummary::
   :toctree: generated/

    align_score


Deprecated functionality
------------------------

Pure Python algorithms (slow; educational-purposes only)

.. autosummary::
   :toctree: generated/

   global_pairwise_align_nucleotide
   global_pairwise_align_protein
   global_pairwise_align
   local_pairwise_align_nucleotide
   local_pairwise_align_protein
   local_pairwise_align


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

Alignment Algorithm Examples
----------------------------
scikit-bio provides pure-Python implementations of Smith-Waterman and
Needleman-Wunsch alignment. These are much slower than implementations
in languages such as C, but serve as useful educational examples as they're
simpler to understand for programmers not experienced with C and are easy to
experiment with. These are _far_ too slow to use in production applications
when run on long sequences or run many times.

Functions are provided for local and global alignment of
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
)
from skbio.alignment._path import AlignPath, PairAlignPath
from skbio.alignment._score import align_score
from skbio.alignment._pair import pair_align, pair_align_nucl, pair_align_prot

__all__ = [
    "TabularMSA",
    "AlignPath",
    "PairAlignPath",
    "global_pairwise_align",
    "global_pairwise_align_nucleotide",
    "global_pairwise_align_protein",
    "local_pairwise_align",
    "local_pairwise_align_nucleotide",
    "local_pairwise_align_protein",
    "align_score",
    "pair_align",
    "pair_align_nucl",
    "pair_align_prot",
]
