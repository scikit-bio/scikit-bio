r"""
Performing Striped Smith Waterman Alignments (:mod:`skbio.core.ssw`)
====================================================================
.. currentmodule:: skbio.core.ssw

This module is a wrapper for a native c implementation of
Striped Smith Waterman Alignment

Classes
-------

.. autosummary::
    :toctree: generated/

    StripedSmithWaterman
    AlignmentStructure

Functions
---------

.. autosummary::
    :toctree: generated/

    align_striped_smith_waterman

Examples
--------
Using the convenient ``align_striped_smith_waterman`` function:

>>> from skbio.core.ssw import align_striped_smith_waterman
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

>>> from skbio.core.ssw import StripedSmithWaterman
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

>>> from skbio.core.ssw import StripedSmithWaterman
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
# -----------------------------------------------------------------------------
# Copyright (c) 2014--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from skbio.core.ssw.ssw_wrapper import (StripedSmithWaterman,
                                        AlignmentStructure,
                                        align_striped_smith_waterman)

__all__ = ['StripedSmithWaterman', 'AlignmentStructure',
           'align_striped_smith_waterman']

from numpy.testing import Tester
test = Tester().test
