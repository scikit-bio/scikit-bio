r"""
Sequences (:mod:`skbio.sequence`)
=================================

.. currentmodule:: skbio.sequence

This module provides classes for storing and working with biological sequences,
including generic sequences which have no restrictions on which characters can
be included, and sequences based on IUPAC-defined sets of allowed characters
(including degenerate characters), including ``DNA``, ``RNA`` and ``Protein``
sequences. Common operations are defined as methods, for example computing the
reverse complement of a DNA sequence, or searching for N-glycosylation motifs
in ``Protein`` sequences. Class attributes are available to obtain valid
character sets, complement maps for different sequence types, and for obtaining
degenerate character definitions. Additionally this module defines the
``GeneticCode`` class, which represents an immutable object that translates DNA
or RNA sequences into protein sequences.

Classes
-------

.. autosummary::
   :toctree: generated/

   Sequence
   IUPACSequence
   NucleotideSequence
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

   GeneticCodeError
   GeneticCodeInitError
   InvalidCodonError

Examples
--------
>>> from skbio import DNA, RNA

New sequences are created with optional metadata and positional metadata
fields. Metadata is stored as a Python dict, while positional metadata
becomes a Pandas DataFrame.

>>> d = DNA('ACCGGGTA')
>>> d = DNA('ACCGGGTA', metadata={'id':"my-sequence", 'description':"GFP"},
...          positional_metadata={'quality':[22, 25, 22, 18, 23, 25, 25, 25]})
>>> d = DNA('ACCGGTA', metadata={'id':"my-sequence"})

New sequences can also be created from existing sequences, for example as their
reverse complement or degapped (i.e., unaligned) version.

>>> d1 = DNA('.ACC--GGG-TA...', metadata={'id':'my-sequence'})
>>> d2 = d1.degap()
>>> d2
DNA
-----------------------------
Metadata:
    'id': 'my-sequence'
Stats:
    length: 8
    has gaps: False
    has degenerates: False
    has non-degenerates: True
-----------------------------
0 ACCGGGTA
>>> d3 = d2.reverse_complement()
>>> d3
DNA
-----------------------------
Metadata:
    'id': 'my-sequence'
Stats:
    length: 8
    has gaps: False
    has degenerates: False
    has non-degenerates: True
-----------------------------
0 TACCCGGT

It's also straightforward to compute distances between sequences (optionally
using user-defined distance metrics, the default is Hamming distance which
requires that the sequences being compared are the same length) for use in
sequence clustering, phylogenetic reconstruction, etc.

>>> r1 = RNA('GACCCGCUUU')
>>> r2 = RNA('GCCCCCCUUU')
>>> r1.distance(r2)
0.2

Similarly, you can calculate the percent (dis)similarity between a pair of
aligned sequences.

>>> r3 = RNA('ACCGUUAGUC')
>>> r4 = RNA('ACGGGU--UC')
>>> r3.match_frequency(r4, relative=True)
0.6
>>> r3.mismatch_frequency(r4, relative=True)
0.4

Sequences can be searched for known motif types. This returns the slices that
describe the matches.

>>> r5 = RNA('AGG-GGACUGAA')
>>> for motif in r5.find_motifs('purine-run', min_length=2):
...     motif
slice(0, 3, None)
slice(4, 7, None)
slice(9, 12, None)

Those slices can be used to extract the relevant subsequences.

>>> for motif in r5.find_motifs('purine-run', min_length=2):
...     r5[motif]
...     print('')
RNA
-----------------------------
Stats:
    length: 3
    has gaps: False
    has degenerates: False
    has non-degenerates: True
-----------------------------
0 AGG
<BLANKLINE>
RNA
-----------------------------
Stats:
    length: 3
    has gaps: False
    has degenerates: False
    has non-degenerates: True
-----------------------------
0 GGA
<BLANKLINE>
RNA
-----------------------------
Stats:
    length: 3
    has gaps: False
    has degenerates: False
    has non-degenerates: True
-----------------------------
0 GAA
<BLANKLINE>

And gaps or other features can be ignored while searching, as these may disrupt
otherwise meaningful motifs.

>>> for motif in r5.find_motifs('purine-run', min_length=2, ignore=r5.gaps()):
...     r5[motif]
...     print('')
RNA
-----------------------------
Stats:
    length: 7
    has gaps: True
    has degenerates: False
    has non-degenerates: True
-----------------------------
0 AGG-GGA
<BLANKLINE>
RNA
-----------------------------
Stats:
    length: 3
    has gaps: False
    has degenerates: False
    has non-degenerates: True
-----------------------------
0 GAA
<BLANKLINE>

In the above example, removing gaps from the resulting motif matches is easily
achieved, as the sliced matches themselves are sequences of the same type as
the input.

>>> for motif in r5.find_motifs('purine-run', min_length=2, ignore=r5.gaps()):
...     r5[motif].degap()
...     print('')
RNA
-----------------------------
Stats:
    length: 6
    has gaps: False
    has degenerates: False
    has non-degenerates: True
-----------------------------
0 AGGGGA
<BLANKLINE>
RNA
-----------------------------
Stats:
    length: 3
    has gaps: False
    has degenerates: False
    has non-degenerates: True
-----------------------------
0 GAA
<BLANKLINE>

Sequences can similarly be searched for arbitrary patterns using regular
expressions.

>>> for match in r5.find_with_regex('(G+AC[UT])'):
...     match
slice(4, 9, None)

Class-level methods contain information about the molecule types.

>>> DNA.degenerate_map['B']
set(['C', 'T', 'G'])

>>> RNA.degenerate_map['B']
set(['C', 'U', 'G'])

Creating and using a ``GeneticCode`` object:

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

Nucleotide sequences can be translated using a ``GeneticCode`` object.

>>> d6 = DNA('ATGTCTAAATGA')
>>> from skbio.sequence import genetic_code
>>> gc = genetic_code(11)
>>> gc.translate(d6)
Protein
-----------------------------
Stats:
    length: 4
    has gaps: False
    has degenerates: False
    has non-degenerates: True
-----------------------------
0 MSK*

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.util import TestRunner

from ._exception import (GeneticCodeError, GeneticCodeInitError,
                         InvalidCodonError)
from ._sequence import Sequence
from ._iupac_sequence import IUPACSequence
from ._nucleotide_sequence import NucleotideSequence
from ._protein import Protein
from ._dna import DNA
from ._rna import RNA
from ._genetic_code import GeneticCode

__all__ = ['GeneticCodeError', 'GeneticCodeInitError', 'InvalidCodonError',
           'Sequence', 'IUPACSequence', 'NucleotideSequence',
           'Protein', 'DNA', 'RNA', 'GeneticCode']

test = TestRunner(__file__).test
