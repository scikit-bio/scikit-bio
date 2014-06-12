r"""
Sequence collections and alignments (:mod:`skbio.core.alignment.alignment`)
===========================================================================

.. currentmodule:: skbio.core.alignment.alignment

This module provides functionality for working with biological sequence
collections and alignments. These can be composed of generic sequences,
nucelotide sequences, DNA sequences, and RNA sequences. By default, input is
not validated, except that sequence ids must be unique, but all
contructor methods take a validate option which checks different features of
the input based on ``SequenceCollection`` type.

Classes
-------

.. autosummary::
   :toctree: generated/

   SequenceCollection
   Alignment

Examples
--------
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

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .alignment import Alignment, SequenceCollection

__all__ = ['Alignment', 'SequenceCollection']

from numpy.testing import Tester
test = Tester().test
