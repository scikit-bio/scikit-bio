r"""
Biological sequences (:mod:`skbio.sequence`)
============================================

.. currentmodule:: skbio.sequence

This module provides functionality for working with biological sequences,
including generic sequences, nucelotide sequences, DNA sequences, and RNA
sequences. Class methods and attributes are also available to obtain valid
character sets, complement maps for different sequence types, and for
obtaining degenerate character definitions. Additionaly this module defines the
``GeneticCode`` class, which represents an immutable object that translates RNA
or DNA strings to amino acid sequences.

Classes
-------

.. autosummary::
   :toctree: generated/

   BiologicalSequence
   NucleotideSequence
   DNASequence
   RNASequence
   ProteinSequence
   GeneticCode

Functions
---------

.. autosummary::
   :toctree: generated/

   genetic_code

Exceptions
----------

.. autosummary::
   :toctree: generated/

   BiologicalSequenceError
   GeneticCodeError
   GeneticCodeInitError
   InvalidCodonError

Examples
--------
>>> from skbio.sequence import DNASequence, RNASequence

New sequences are created with optional id and description fields.

>>> d1 = DNASequence('ACC--G-GGTA..')
>>> d1 = DNASequence('ACC--G-GGTA..',id="seq1")
>>> d1 = DNASequence('ACC--G-GGTA..',id="seq1",description="GFP")

New sequences can also be created from existing sequences, for example as their
reverse complement or degapped (i.e., unaligned) version.

>>> d2 = d1.degap()
>>> d1
<DNASequence: ACC--G-GGT... (length: 13)>
>>> d2
<DNASequence: ACCGGGTA (length: 8)>
>>> d3 = d2.reverse_complement()
>>> d3
<DNASequence: TACCCGGT (length: 8)>

It's also straight-forward to compute distances between sequences (optionally
using user-defined distance metrics, default is Hamming distance) for use in
sequence clustering, phylogenetic reconstruction, etc.

>>> d4 = DNASequence('GACCCGCT')
>>> d5 = DNASequence('GACCCCCT')
>>> d3.distance(d4)
0.25
>>> d3.distance(d5)
0.375

Class-level methods contain information about the molecule types.

>>> DNASequence.iupac_degeneracies()['B']
set(['C', 'T', 'G'])

>>> RNASequence.iupac_degeneracies()['B']
set(['C', 'U', 'G'])

>>> DNASequence.is_gap('-')
True

Creating and using a ``GeneticCode`` object

>>> from skbio.sequence import genetic_code
>>> from pprint import pprint
>>> sgc = genetic_code(1)
>>> sgc
GeneticCode(FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG)
>>> sgc['UUU'] == 'F'
True
>>> sgc['TTT'] == 'F'
True
>>> sgc['F'] == ['TTT', 'TTC']          #in arbitrary order
True
>>> sgc['*'] == ['TAA', 'TAG', 'TGA']   #in arbitrary order
True

Retrieving the anticodons of the object

>>> pprint(sgc.anticodons)
{'*': ['TTA', 'CTA', 'TCA'],
 'A': ['AGC', 'GGC', 'TGC', 'CGC'],
 'C': ['ACA', 'GCA'],
 'D': ['ATC', 'GTC'],
 'E': ['TTC', 'CTC'],
 'F': ['AAA', 'GAA'],
 'G': ['ACC', 'GCC', 'TCC', 'CCC'],
 'H': ['ATG', 'GTG'],
 'I': ['AAT', 'GAT', 'TAT'],
 'K': ['TTT', 'CTT'],
 'L': ['TAA', 'CAA', 'AAG', 'GAG', 'TAG', 'CAG'],
 'M': ['CAT'],
 'N': ['ATT', 'GTT'],
 'P': ['AGG', 'GGG', 'TGG', 'CGG'],
 'Q': ['TTG', 'CTG'],
 'R': ['ACG', 'GCG', 'TCG', 'CCG', 'TCT', 'CCT'],
 'S': ['AGA', 'GGA', 'TGA', 'CGA', 'ACT', 'GCT'],
 'T': ['AGT', 'GGT', 'TGT', 'CGT'],
 'V': ['AAC', 'GAC', 'TAC', 'CAC'],
 'W': ['CCA'],
 'Y': ['ATA', 'GTA']}

NucleotideSequences can be translated using a ``GeneticCode`` object.

>>> d6 = DNASequence('ATGTCTAAATGA')
>>> from skbio.sequence import genetic_code
>>> gc = genetic_code(11)
>>> gc.translate(d6)
<ProteinSequence: MSK* (length: 4)>

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._exception import (BiologicalSequenceError, GeneticCodeError,
                         GeneticCodeInitError, InvalidCodonError)
from ._sequence import (BiologicalSequence, NucleotideSequence, DNASequence,
                        RNASequence, ProteinSequence, DNA, RNA, Protein)
from ._genetic_code import GeneticCode, genetic_code

__all__ = ['BiologicalSequenceError', 'GeneticCodeError',
           'GeneticCodeInitError', 'InvalidCodonError', 'BiologicalSequence',
           'NucleotideSequence', 'DNASequence', 'RNASequence',
           'ProteinSequence', 'DNA', 'RNA', 'Protein', 'GeneticCode',
           'genetic_code']

from numpy.testing import Tester
test = Tester().test
