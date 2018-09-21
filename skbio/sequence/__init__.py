r"""
Sequences (:mod:`skbio.sequence`)
=================================

.. currentmodule:: skbio.sequence

This module provides classes for storing and working with sequences, including
generic/nonbiological sequences which have no alphabet restrictions
(``Sequence``) and sequences based on IUPAC-defined alphabets (``DNA``,
``RNA``, ``Protein``). Common operations are defined as methods, for example
computing the reverse complement of a DNA sequence, or searching for
N-glycosylation motifs in protein sequences. Class attributes provide valid
character sets, complement maps for different sequence types, and degenerate
character definitions. Additionally this module defines the ``GeneticCode``
class, which represents an immutable object that translates DNA or RNA
sequences into protein sequences.

The primary information stored for each different type of sequence object is
the underlying sequence data itself. This is stored as an immutable numpy
array. Additionally, each type of sequence may include optional metadata
and positional metadata. Note that metadata and positional metadata are
mutable.

Classes
-------

.. autosummary::
   :toctree: generated/

   Sequence
   GrammaredSequence
   DNA
   RNA
   Protein
   GeneticCode

Subpackages
-----------

.. autosummary::
   :toctree: generated/

   distance

Examples
--------
New sequences are created with optional metadata and positional metadata.
Metadata is stored as a Python ``dict``, while positional metadata is stored as
a pandas ``DataFrame``.

>>> from skbio import DNA, RNA
>>> d = DNA('ACCGGGTA', metadata={'id':"my-sequence", 'description':"GFP"},
...          positional_metadata={'quality':[22, 25, 22, 18, 23, 25, 25, 25]})
>>> d
DNA
-----------------------------
Metadata:
    'description': 'GFP'
    'id': 'my-sequence'
Positional metadata:
    'quality': <dtype: int64>
Stats:
    length: 8
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 62.50%
-----------------------------
0 ACCGGGTA

New sequences can also be created from existing sequences, for example as their
reverse complement or degapped (i.e., unaligned) version.

>>> d1 = DNA('.ACC--GGG-TA...', metadata={'id':'my-sequence'})
>>> d2 = d1.degap()
>>> d2
DNA
--------------------------
Metadata:
    'id': 'my-sequence'
Stats:
    length: 8
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 62.50%
--------------------------
0 ACCGGGTA
>>> d3 = d2.reverse_complement()
>>> d3
DNA
--------------------------
Metadata:
    'id': 'my-sequence'
Stats:
    length: 8
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 62.50%
--------------------------
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
--------------------------
Stats:
    length: 3
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 66.67%
--------------------------
0 AGG
<BLANKLINE>
RNA
--------------------------
Stats:
    length: 3
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 66.67%
--------------------------
0 GGA
<BLANKLINE>
RNA
--------------------------
Stats:
    length: 3
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 33.33%
--------------------------
0 GAA
<BLANKLINE>

And gaps or other features can be ignored while searching, as these may disrupt
otherwise meaningful motifs.

>>> for motif in r5.find_motifs('purine-run', min_length=2, ignore=r5.gaps()):
...     r5[motif]
...     print('')
RNA
--------------------------
Stats:
    length: 7
    has gaps: True
    has degenerates: False
    has definites: True
    GC-content: 66.67%
--------------------------
0 AGG-GGA
<BLANKLINE>
RNA
--------------------------
Stats:
    length: 3
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 33.33%
--------------------------
0 GAA
<BLANKLINE>

In the above example, removing gaps from the resulting motif matches is easily
achieved, as the sliced matches themselves are sequences of the same type as
the input.

>>> for motif in r5.find_motifs('purine-run', min_length=2, ignore=r5.gaps()):
...     r5[motif].degap()
...     print('')
RNA
--------------------------
Stats:
    length: 6
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 66.67%
--------------------------
0 AGGGGA
<BLANKLINE>
RNA
--------------------------
Stats:
    length: 3
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 33.33%
--------------------------
0 GAA
<BLANKLINE>

Sequences can similarly be searched for arbitrary patterns using regular
expressions.

>>> for match in r5.find_with_regex('(G+AC[UT])'):
...     match
slice(4, 9, None)

DNA can be transcribed to RNA:

>>> dna = DNA('ATGTGTATTTGA')
>>> rna = dna.transcribe()
>>> rna
RNA
--------------------------
Stats:
    length: 12
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 25.00%
--------------------------
0 AUGUGUAUUU GA

Both DNA and RNA can be translated into a protein sequence. For example, let's
translate our DNA and RNA sequences using NCBI's standard genetic code (table
ID 1, the default genetic code in scikit-bio):

>>> protein_from_dna = dna.translate()
>>> protein_from_dna
Protein
--------------------------
Stats:
    length: 4
    has gaps: False
    has degenerates: False
    has definites: True
    has stops: True
--------------------------
0 MCI*
>>> protein_from_rna = rna.translate()
>>> protein_from_rna
Protein
--------------------------
Stats:
    length: 4
    has gaps: False
    has degenerates: False
    has definites: True
    has stops: True
--------------------------
0 MCI*

The two translations are equivalent:

>>> protein_from_dna == protein_from_rna
True

Class-level methods contain information about the molecule types.

>>> sorted(DNA.degenerate_map['B'])
['C', 'G', 'T']

>>> sorted(RNA.degenerate_map['B'])
['C', 'G', 'U']

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._sequence import Sequence
from ._protein import Protein
from ._dna import DNA
from ._rna import RNA
from ._genetic_code import GeneticCode
from ._grammared_sequence import GrammaredSequence

__all__ = ['Sequence', 'Protein', 'DNA', 'RNA', 'GeneticCode',
           'GrammaredSequence']
