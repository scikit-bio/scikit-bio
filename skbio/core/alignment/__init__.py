r"""
Sequence collections and alignments (:mod:`skbio.core.alignment`)
=================================================================

.. currentmodule:: skbio.core.alignment

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

Optimized Alignment Algorithms
------------------------------

.. autosummary::
   :toctree: generated/

   StripedSmithWaterman
   AlignmentStructure
   align_striped_smith_waterman

Slow Alignment Algorithms
-------------------------

.. autosummary::
   :toctree: generated/

   global_pairwise_align_nucleotide
   global_pairwise_align_protein
   global_pairwise_align
   local_pairwise_align_nucleotide
   local_pairwise_align_protein
   local_pairwise_align

Data Structure Examples
-----------------------
>>> from StringIO import StringIO
>>> from skbio.core.alignment import SequenceCollection, Alignment
>>> from skbio.core.sequence import DNA
>>> seqs = [DNA("ACC--G-GGTA..", id="seq1"),
...     DNA("TCC--G-GGCA..", id="seqs2")]
>>> a1 = Alignment(seqs)
>>> a1
<Alignment: n=2; mean +/- std length=13.00 +/- 0.00>

>>> seqs = [DNA("ACCGGG", id="seq1"),
...     DNA("TCCGGGCA", id="seq2")]
>>> s1 = SequenceCollection(seqs)
>>> s1
<SequenceCollection: n=2; mean +/- std length=7.00 +/- 1.00>

>>> from skbio.parse.sequences import parse_fasta
>>> fasta_f = StringIO('>seq1\n'
...                    'CGATGTCGATCGATCGATCGATCAG\n'
...                    '>seq2\n'
...                    'CATCGATCGATCGATGCATGCATGCATG\n')
>>> s1 = SequenceCollection.from_fasta_records(parse_fasta(fasta_f), DNA)
>>> s1
<SequenceCollection: n=2; mean +/- std length=26.50 +/- 1.50>

Optimized Alignment Algorithm Examples
--------------------------------------
Using the convenient ``align_striped_smith_waterman`` function:

>>> from skbio.core.alignment import align_striped_smith_waterman
>>> alignment = align_striped_smith_waterman(
...                 "ACTAAGGCTCTCTACCCCTCTCAGAGA",
...                 "ACTAAGGCTCCTAACCCCCTTTTCTCAGA"
...             )
>>> print alignment
{
    'optimal_alignment_score': 27,
    'suboptimal_alignment_score': 21,
    'query_begin': 0,
    'query_end': 24,
    'target_begin': 0,
    'target_end_optimal': 28,
    'target_end_suboptimal': 12,
    'cigar': '10M1I2M1D5M4D7M',
    'query_sequence': 'ACTAAGGCTCTCTACCCCTCTCAGAGA',
    'target_sequence': 'ACTAAGGCTCCTAACCCCCTTTTCTCAGA'
}

Using the ``StripedSmithWaterman`` object:

>>> from skbio.core.alignment import StripedSmithWaterman
>>> query = StripedSmithWaterman("ACTAAGGCTCTCTACCCCTCTCAGAGA")
>>> alignment = query("AAAAAACTCTCTAAACTCACTAAGGCTCTCTACCCCTCTTCAGAGAAGTCGA")
>>> print alignment
{
    'optimal_alignment_score': 49,
    'suboptimal_alignment_score': 24,
    'query_begin': 0,
    'query_end': 26,
    'target_begin': 18,
    'target_end_optimal': 45,
    'target_end_suboptimal': 29,
    'cigar': '20M1D7M',
    'query_sequence': 'ACTAAGGCTCTCTACCCCTCTCAGAGA',
    'target_sequence': 'AAAAAACTCTCTAAACTCACTAAGGCTCTCTACCCCTCTTCAGAGAAGTCGA'
}

Using the ``StripedSmithWaterman`` object for multiple targets in an efficient
way and finding the aligned sequence representations:

>>> from skbio.core.alignment import StripedSmithWaterman
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
>>> print alignments[0]
{
    'optimal_alignment_score': 38,
    'suboptimal_alignment_score': 14,
    'query_begin': 0,
    'query_end': 26,
    'target_begin': 4,
    'target_end_optimal': 32,
    'target_end_suboptimal': 15,
    'cigar': '3M1I6M3D17M',
    'query_sequence': 'ACTAAGGCTCTCTACCCCTCTCAGAGA',
    'target_sequence': 'GCTAACTAGGCTCCCTTCTACCCCTCTCAGAGA'
}
>>> print alignments[0].get_aligned_query_sequence()
ACTAAGGCT---CTCTACCCCTCTCAGAGA
>>> print alignments[0].get_aligned_target_sequence()
ACT-AGGCTCCCTTCTACCCCTCTCAGAGA

Slow Alignment Algorithm Examples
---------------------------------
scikit-bio also implements Smith-Waterman and Needleman-Wunsch alignment
in python. These are much slower than the methods described above, but provide
useful education examples as they're simpler to experiment with. Convenience
wrappers are provided for local and global alignment of protein and nucleotide
sequences. The ``global`` and ``local`` functions differ in the underlying
algorithm that is applied (``global`` uses Needleman-Wunsch while local uses
Smith-Waterman), and ``protein`` and ``nucleotide`` differ in their scoring of
matches and mismatches, and the default gap penalties.

Here we locally align a pair of protein sequences using gap open penalty
of 11 and a gap extend penalty of 1 (in other words, it is much more
costly to open a new gap than extend an existing one).
>>> from skbio import local_pairwise_align_protein
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

And we can view the score of the alignment is accessible using the ``score``
method:
>>> r.score()
25.0

Similarly, we can perform global alignment of nucleotide sequences, and print
the resulting alignment as fasta records:
>>> from skbio import global_pairwise_align_nucleotide
>>> s1 = "GCGTGCCTAAGGTATGCAAG"
>>> s2 = "ACGTGCCTAGGTACGCAAG"
>>> r = global_pairwise_align_nucleotide(s1, s2)
>>> print(r.to_fasta())
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

from .alignment import Alignment, SequenceCollection
from .ssw.ssw_wrapper import (
    StripedSmithWaterman, AlignmentStructure, align_striped_smith_waterman)

__all__ = ['Alignment', 'SequenceCollection', 'StripedSmithWaterman',
           'AlignmentStructure', 'align_striped_smith_waterman']

from numpy.testing import Tester
test = Tester().test
