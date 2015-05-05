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

   Sequence
   DNA
   RNA
   Protein
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

   SequenceError
   GeneticCodeError
   GeneticCodeInitError
   InvalidCodonError

Examples
--------
>>> from skbio.sequence import DNA, RNA

New sequences are created with optional id and description fields.

>>> d1 = DNA('ACC--G-GGTA..')
>>> d1 = DNA('ACC--G-GGTA..',id="seq1")
>>> d1 = DNA('ACC--G-GGTA..',id="seq1",description="GFP")

New sequences can also be created from existing sequences, for example as their
reverse complement or degapped (i.e., unaligned) version.

>>> d2 = d1.degap()
>>> d1
DNA('ACC--G-GGTA..', length=13, id='seq1', description='GFP')
>>> d2
DNA('ACCGGGTA', length=8, id='seq1', description='GFP')
>>> d3 = d2.reverse_complement()
>>> d3
DNA('TACCCGGT', length=8, id='seq1', description='GFP')

It's also straight-forward to compute distances between sequences (optionally
using user-defined distance metrics, default is Hamming distance) for use in
sequence clustering, phylogenetic reconstruction, etc.

>>> d4 = DNA('GACCCGCT')
>>> d5 = DNA('GACCCCCT')
>>> d3.distance(d4)
0.25
>>> d3.distance(d5)
0.375

Class-level methods contain information about the molecule types.

>>> DNA.degenerate_map['B']
set(['C', 'T', 'G'])

>>> RNA.degenerate_map['B']
set(['C', 'U', 'G'])

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

>>> d6 = DNA('ATGTCTAAATGA')
>>> from skbio.sequence import genetic_code
>>> gc = genetic_code(11)
>>> gc.translate(d6)
Protein('MSK*', length=4)


.. shownumpydoc
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.util import TestRunner

from ._exception import (SequenceError, GeneticCodeError,
                         GeneticCodeInitError, InvalidCodonError)
from ._sequence import Sequence
from ._protein import Protein
from ._dna import DNA
from ._rna import RNA
from ._genetic_code import GeneticCode, genetic_code

__all__ = ['SequenceError', 'GeneticCodeError',
           'GeneticCodeInitError', 'InvalidCodonError', 'Sequence',
           'Protein', 'DNA', 'RNA', 'GeneticCode',
           'genetic_code']

test = TestRunner(__file__).test
