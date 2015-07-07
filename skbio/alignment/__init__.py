r"""
Sequence collections and alignments (:mod:`skbio.alignment`)
============================================================

.. currentmodule:: skbio.alignment

This module provides functionality for working with biological sequence
collections and alignments. These can be composed of generic sequences,
nucelotide sequences, DNA sequences, and RNA sequences. By default, input is
not validated, except that sequence ids must be unique, but all
contructor methods take a validate option which checks different features of
the input based on ``SequenceCollection`` type.

Data Structures
---------------

.. autosummary::
   :toctree: generated/

   SequenceCollection
   Alignment

Optimized (i.e., production-ready) Alignment Algorithms
-------------------------------------------------------

.. autosummary::
   :toctree: generated/

   StripedSmithWaterman
   AlignmentStructure
   local_pairwise_align_ssw

Slow (i.e., educational-purposes only) Alignment Algorithms
-----------------------------------------------------------

.. autosummary::
   :toctree: generated/

   global_pairwise_align_nucleotide
   global_pairwise_align_protein
   global_pairwise_align
   local_pairwise_align_nucleotide
   local_pairwise_align_protein
   local_pairwise_align

General functionality
---------------------

.. autosummary::
   :toctree: generated/

    make_identity_substitution_matrix

Exceptions
----------

.. autosummary::
   :toctree: generated/

   SequenceCollectionError
   AlignmentError

Data Structure Examples
-----------------------
>>> from skbio import SequenceCollection, Alignment, DNA
>>> seqs = [DNA("ACC--G-GGTA..", metadata={'id':"seq1"}),
...         DNA("TCC--G-GGCA..", metadata={'id':"seqs2"})]
>>> a1 = Alignment(seqs)
>>> a1
<Alignment: n=2; mean +/- std length=13.00 +/- 0.00>

>>> seqs = [DNA("ACCGGG", metadata={'id':"seq1"}),
...         DNA("TCCGGGCA", metadata={'id':"seq2"})]
>>> s1 = SequenceCollection(seqs)
>>> s1
<SequenceCollection: n=2; mean +/- std length=7.00 +/- 1.00>

>>> fasta_lines = [u'>seq1\n',
...                u'CGATGTCGATCGATCGATCGATCAG\n',
...                u'>seq2\n',
...                u'CATCGATCGATCGATGCATGCATGCATG\n']
>>> s1 = SequenceCollection.read(fasta_lines, constructor=DNA)
>>> s1
<SequenceCollection: n=2; mean +/- std length=26.50 +/- 1.50>

Alignment Algorithm Examples
----------------------------

Optimized Alignment Algorithm Examples
--------------------------------------
Using the convenient ``local_pairwise_align_ssw`` function:

>>> from skbio.alignment import local_pairwise_align_ssw
>>> alignment = local_pairwise_align_ssw(
...                 "ACTAAGGCTCTCTACCCCTCTCAGAGA",
...                 "ACTAAGGCTCCTAACCCCCTTTTCTCAGA"
...             )
>>> print(alignment)
>query
ACTAAGGCTCTC-TACCC----CTCTCAGA
>target
ACTAAGGCTC-CTAACCCCCTTTTCTCAGA
<BLANKLINE>

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
ACTAAGGCT-...
ACT-AGGCTC...
Score: 38
Length: 30
>>> print(alignments[0].aligned_query_sequence)
ACTAAGGCT---CTCTACCCCTCTCAGAGA
>>> print(alignments[0].aligned_target_sequence)
ACT-AGGCTCCCTTCTACCCCTCTCAGAGA

Slow Alignment Algorithm Examples
---------------------------------
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

>>> from skbio.alignment import local_pairwise_align_protein
>>> s1 = "HEAGAWGHEE"
>>> s2 = "PAWHEAE"
>>> r = local_pairwise_align_protein(s1, s2, 11, 1)

This returns an ``skbio.Alignment`` object. We can look at the aligned
sequences:

>>> print(str(r[0]))
AWGHE
>>> print(str(r[1]))
AW-HE

We can identify the start and end positions of each aligned sequence
as follows:

>>> r.start_end_positions()
[(4, 8), (1, 4)]

And we can view the score of the alignment using the ``score`` method:

>>> r.score()
25.0

Similarly, we can perform global alignment of nucleotide sequences, and print
the resulting alignment in FASTA format:

>>> from skbio.alignment import global_pairwise_align_nucleotide
>>> s1 = "GCGTGCCTAAGGTATGCAAG"
>>> s2 = "ACGTGCCTAGGTACGCAAG"
>>> r = global_pairwise_align_nucleotide(s1, s2)
>>> print(r)
>0
GCGTGCCTAAGGTATGCAAG
>1
ACGTGCCTA-GGTACGCAAG
<BLANKLINE>


"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util import TestRunner

from ._alignment import Alignment, SequenceCollection
from ._pairwise import (
    local_pairwise_align_nucleotide, local_pairwise_align_protein,
    local_pairwise_align, global_pairwise_align_nucleotide,
    global_pairwise_align_protein, global_pairwise_align,
    make_identity_substitution_matrix, local_pairwise_align_ssw
)
from skbio.alignment._ssw_wrapper import (
    StripedSmithWaterman, AlignmentStructure)
from ._exception import (SequenceCollectionError, AlignmentError)

__all__ = ['Alignment', 'SequenceCollection',
           'StripedSmithWaterman', 'AlignmentStructure',
           'local_pairwise_align_ssw', 'SequenceCollectionError',
           'AlignmentError', 'global_pairwise_align',
           'global_pairwise_align_nucleotide', 'global_pairwise_align_protein',
           'local_pairwise_align', 'local_pairwise_align_nucleotide',
           'local_pairwise_align_protein', 'make_identity_substitution_matrix']

test = TestRunner(__file__).test
