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
   StockholmAlignment

Alignment Algorithms
--------------------

.. autosummary::
   :toctree: generated/

   StripedSmithWaterman
   AlignmentStructure
   align_striped_smith_waterman

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

>>> from skbio.core.sequence import RNA
>>> from skbio.core.alignment import StockholmAlignment
>>> from StringIO import StringIO
>>> sto_in = StringIO("# STOCKHOLM 1.0\\n"
...                   "seq1     ACC--G-GGGU\\nseq2     TCC--G-GGGA\\n"
...                   "#=GC SS_cons (((.....)))\\n//")
>>> sto_records = StockholmAlignment.from_file(sto_in, RNA)
>>> sto = next(sto_records)
>>> print(sto)
# STOCKHOLM 1.0
seq1          ACC--G-GGGU
seq2          TCC--G-GGGA
#=GC SS_cons  (((.....)))
//
>>> sto.gc
{'SS_cons': '(((.....)))'}

Alignment Algorithm Examples
----------------------------
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

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .alignment import Alignment, SequenceCollection, StockholmAlignment
from .ssw.ssw_wrapper import (
    StripedSmithWaterman, AlignmentStructure, align_striped_smith_waterman)

__all__ = ['Alignment', 'SequenceCollection', 'StockholmAlignment',
           'StripedSmithWaterman', 'AlignmentStructure',
           'align_striped_smith_waterman']

from numpy.testing import Tester
test = Tester().test
