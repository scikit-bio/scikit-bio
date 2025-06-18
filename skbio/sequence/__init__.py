r"""Biological Sequences (:mod:`skbio.sequence`)
============================================

.. currentmodule:: skbio.sequence

This module provides functionality for storing and working with sequences, including
molecular sequences based on IUPAC-defined alphabets (:class:`DNA`, :class:`RNA`,
:class:`Protein`), sequences based on custom alphabets (:class:`GrammaredSequence`),
and generic/non-biological sequences with no alphabet restrictions (:class:`Sequence`).

Additionally, this module defines the :class:`GeneticCode` class, which represents an
immutable object that translates DNA or RNA sequences into protein sequences, and
the :class:`SubstitutionMatrix` class, which stores scores of substitutions between
sequence characters.

See the :ref:`sequence_tutorial` section for working with biological sequences using
scikit-bio.


Sequence types
--------------

.. autosummary::
   :toctree: generated/

   Sequence
   GrammaredSequence
   DNA
   RNA
   Protein


Sequence utilities
------------------

.. autosummary::
   :toctree: generated/

   GeneticCode
   SubstitutionMatrix


Distance calculation
--------------------

.. autosummary::
   :toctree: generated/

   distance


Abstract classes
----------------

.. autosummary::
   :toctree: generated/

   NucleotideMixin


.. _sequence_tutorial:

Tutorial
--------

The primary information stored for each different type of sequence object is the
underlying sequence data itself. This is stored as an immutable NumPy array.
Additionally, each type of sequence may include optional metadata and positional
metadata. Note that metadata and positional metadata are mutable.

Common operations are defined as methods, for example computing the reverse complement
of a DNA sequence, or searching for N-glycosylation motifs in protein sequences. Class
attributes provide valid character sets, complement maps for different sequence types,
and degenerate character definitions.

Create a DNA sequence from the sequence data (nucleotides):

>>> from skbio.sequence import DNA
>>> seq = DNA('ACCGGGTA')
>>> seq
DNA
--------------------------
Stats:
    length: 8
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 62.50%
--------------------------
0 ACCGGGTA

New sequences may be created with optional metadata and positional metadata.
Metadata is stored as a Python ``dict``, while positional metadata is stored as
a pandas ``DataFrame``.

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

>>> from skbio.sequence import RNA
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


Sequence metadata
^^^^^^^^^^^^^^^^^

Metadata describe the properties of a sequence or some of its regions/sites. scikit-bio
supports three types of sequence metadata: **metadata** (for the entire sequence),
**positional metadata** (for sites), and **interval metadata** (for regions).

``metadata`` stores the properties of the entire sequence, such as ID, description,
source, function, etc., as a Python dictionary. The following example creates a DNA
sequence with "gene" and "product" properties.

>>> seq = DNA('TGGTATATAGTTTAAACAAAACGAATGATTTCGACTCATTAAATTATGATAATCATATTTACCAA',
...           metadata={'gene': 'TRNR', 'product': 'tRNA-Arg'})
>>> seq
DNA
--------------------------------------------------------------------
Metadata:
    'gene': 'TRNR'
    'product': 'tRNA-Arg'
Stats:
    length: 65
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 23.08%
--------------------------------------------------------------------
0  TGGTATATAG TTTAAACAAA ACGAATGATT TCGACTCATT AAATTATGAT AATCATATTT
60 ACCAA

Access the metadata by key:

>>> seq.metadata['product']
'tRNA-Arg'

Modify the metadata of an existing sequence:

>>> seq.metadata['organism'] = 'Homo sapiens'

Many of scikit-bio's file format parsers can automatically populate metadata fields.
For example, the :mod:`FASTA <skbio.io.format.fasta>` parser reads the sequence ID and
description into the ``id`` and ``description`` fields:

>>> fasta = '''>NP_000591.1 interleukin-6 isoform 1 precursor [Homo sapiens]
... MNSFSTSAFGPVAFSLGLLLVLPAAFPAPVPPGEDSKDVAAPHRQPLTSSERIDKQIRYILDGISALRKE
... TCNKSNMCESSKEALAENNLNLPKMAEKDGCFQSGFNEETCLVKIITGLLEFEVYLEYLQNRFESSEEQA
... RAVQMSTKVLIQFLQKKAKNLDAITTPDPTTNASLLTKLQAQNQWLQDMTTHLILRSFKEFLQSSLRALR
... QM'''
>>> from skbio.sequence import Protein
>>> seq = Protein.read([fasta], format='fasta')
>>> seq
Protein
---------------------------------------------------------------------
Metadata:
    'description': 'interleukin-6 isoform 1 precursor [Homo sapiens]'
    'id': 'NP_000591.1'
Stats:
    length: 212
    has gaps: False
    has degenerates: False
    has definites: True
    has stops: False
---------------------------------------------------------------------
0   MNSFSTSAFG PVAFSLGLLL VLPAAFPAPV PPGEDSKDVA APHRQPLTSS ERIDKQIRYI
60  LDGISALRKE TCNKSNMCES SKEALAENNL NLPKMAEKDG CFQSGFNEET CLVKIITGLL
120 EFEVYLEYLQ NRFESSEEQA RAVQMSTKVL IQFLQKKAKN LDAITTPDPT TNASLLTKLQ
180 AQNQWLQDMT THLILRSFKE FLQSSLRALR QM

``positional_metadata`` stores per-character properties. Each property has the same
number of values as the sequence length. The following example annotates three types
of sites on the human TP53 protein sequence.

>>> seq = (
...     'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAA'
...     'PPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKT'
...     'CPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRN'
...     'TFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGR'
...     'DRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALEL'
...     'KDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD'
... )
>>> import numpy as np
>>> sites = np.full(len(seq), False)
>>> acetylation = sites.copy()
>>> acetylation[[119, 304, 320, 372, 380, 381]] = True
>>> methylation = sites.copy()
>>> methylation[[332, 334, 336, 369, 371, 372, 381]] = True
>>> phosphorylation = sites.copy()
>>> phosphorylation[[8, 14, 17, 19, 32, 36, 45, 54, 182, 268, 314, 391]] = True
>>> seq = Protein(seq, metadata={'id': 'NP_000537.3', 'description': 'TP53'},
...               positional_metadata={'acetylation': acetylation,
...                                    'methylation': methylation,
...                                    'phosphorylation': phosphorylation})

Once created, positional metadata are stored as a dataframe.

>>> print(seq.positional_metadata)
     acetylation  methylation  phosphorylation
0          False        False            False
1          False        False            False
2          False        False            False
3          False        False            False
4          False        False            False
..           ...          ...              ...
388        False        False            False
389        False        False            False
390        False        False            False
391        False        False             True
392        False        False            False
<BLANKLINE>
[393 rows x 3 columns]

Likewise, scikit-bio's :mod:`FASTQ <skbio.io.format.fastq>` parser can automatically
populate the "quality" column of the positional metadata:

>>> fastq = '''@S00001/1
... CGTGCTCCACGCTCTCCTGCACGCCCGAGCCGTTGCGGTCGGCGAGAGCCATGTCGACGGCGGGCTTGAG
... +
... DDAEGGGGIIIGIJKJKKKHIKJAKKKHCGIKKJAJGICKKKAJEJA:G1EJJKJHJFKJIDIEII?=DC
... '''
>>> seq = DNA.read([fastq], format='fastq', phred_offset=33)
>>> seq
DNA
--------------------------------------------------------------------
Metadata:
    'description': ''
    'id': 'S00001/1'
Positional metadata:
    'quality': <dtype: uint8>
Stats:
    length: 70
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 71.43%
--------------------------------------------------------------------
0  CGTGCTCCAC GCTCTCCTGC ACGCCCGAGC CGTTGCGGTC GGCGAGAGCC ATGTCGACGG
60 CGGGCTTGAG

>>> print(seq.positional_metadata['quality'])
0     35
1     35
2     32
3     36
4     38
      ..
65    40
66    30
67    28
68    35
69    34
Name: quality, Length: 70, dtype: uint8

``interval_metadata`` stores properties that span a range in the sequence, such as
genes and regulatory elements. The following example annotates three regions of the
human H4C6 gene.

>>> fasta = '''>NM_003540.4 Homo sapiens H4 clustered histone 6 (H4C6), mRNA
GCAAAAGTTAAGAGTTGTTGTTTGTCTTCGATCATGTCTGGTAGAGGCAAAGGTGGTAAAGGTTTAGGAA
AGGGAGGCGCCAAGCGCCATCGCAAAGTGCTGCGTGACAACATACAGGGCATCACGAAGCCCGCCATCCG
TCGCTTGGCCCGACGCGGCGGCGTGAAACGCATTTCGGGCCTCATTTATGAGGAGACCCGCGGTGTTCTT
AAGGTGTTCCTGGAGAATGTGATACGGGACGCCGTAACCTACACGGAGCACGCCAAGCGTAAGACAGTCA
CTGCAATGGATGTTGTCTACGCGCTCAAGCGCCAGGGACGCACTCTGTACGGCTTTGGTGGCTGAGCCTC
ACCCCGGCTTTTTATTTAACAGCTCACCCATAAAAGGCCCTTTTCAGGGCC'''
>>> seq = DNA.read([fasta], format='fasta')
>>> seq.interval_metadata.add(bounds=[(0, 401)], metadata={'name': 'exon'})
>>> seq.interval_metadata.add(bounds=[(33, 345)], metadata={'name': 'cds'})
>>> seq.interval_metadata.add(bounds=[(385, 401)], metadata={'name': 'stem_loop'})

scikit-bio's :mod:`GFF3 <skbio.io.format.gff3>` parser can read sequence annotation
into interval metadata. The following code reads a genome's sequence and annotation
from separate files and combine.

>>> from skbio.metadata import IntervalMetadata
>>> seq = DNA.read('genomic.fna', format='fasta')  # doctest: +SKIP
>>> seq.interval_metadata = IntervalMetadata.read(  # doctest: +SKIP
...     'genomic.gff', format='gff3', seq_id=seq.metadata['id'])

If the GFF3 contains the sequence, one can read the sequence and annotation all at
once:

>>> seq = DNA.read('genomic.gff', format='gff3')  # doctest: +SKIP

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._sequence import Sequence
from ._protein import Protein
from ._dna import DNA
from ._rna import RNA
from ._genetic_code import GeneticCode
from ._grammared_sequence import GrammaredSequence
from ._substitution import SubstitutionMatrix
from ._nucleotide_mixin import NucleotideMixin

__all__ = [
    "Sequence",
    "Protein",
    "DNA",
    "RNA",
    "GeneticCode",
    "GrammaredSequence",
    "SubstitutionMatrix",
    "NucleotideMixin",
]
