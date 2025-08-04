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

See the |sequence_tutorial|_ section for working with biological sequences using
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


.. |sequence_tutorial| replace:: **Tutorial**
.. _sequence_tutorial:

Tutorial
--------

Creating sequences
^^^^^^^^^^^^^^^^^^

Create a DNA sequence from a string of nucleotide characters:

>>> from skbio.sequence import DNA
>>> seq = DNA('GAATTC')
>>> seq
DNA
--------------------------
Stats:
    length: 6
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 33.33%
--------------------------
0 GAATTC

Optionally, one may attach metadata to a sequence object. Details of metadata are
introduced in :ref:`annotate_sequences` below.

>>> seq = DNA('GAATTC', metadata={'id': 'S1', 'description': 'EcoRI recognition site'})
>>> seq
DNA
-------------------------------------------
Metadata:
    'description': 'EcoRI recognition site'
    'id': 'S1'
Stats:
    length: 6
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 33.33%
-------------------------------------------
0 GAATTC

Sequences can be read from :mod:`file formats <skbio.io>` recognizable by scikit-bio.
For example, to read a DNA sequence from a :mod:`FASTA <skbio.io.format.fasta>` file:

>>> seq = DNA.read('input.fa', format='fasta')  # doctest: +SKIP

The filename "input.fa" can also be a file handle or other readable objects.

To read all sequences from a multi-FASTA file into a list of sequence objects, you will
need:

>>> from skbio.io import read as sk_read
>>> seqs = list(sk_read('input.fa', format='fasta', constructor=DNA))  # doctest: +SKIP


Sequence grammar
^^^^^^^^^^^^^^^^

The three common biological sequence types: ``DNA``, ``RNA`` and ``Protein``, are
:class:`grammared sequences <GrammaredSequence>`. That is, each of them has a defined
alphabet (character set), which typically consists of:

- :attr:`definite characters <GrammaredSequence.definite_chars>`, such as ``ACGT``,
  which represent the four canonical nucleotides in a DNA sequence,
- :attr:`degenerate characters <GrammaredSequence.degenerate_chars>`, such as ``R``,
  which represents ``A`` or ``G``,
- :attr:`gap character(s) <GrammaredSequence.gap_chars>`, such as ``-``, and
- a :attr:`wildcard character <GrammaredSequence.wildcard_char>`, such as ``N``, which
  is also a degenerate character.

Use of any of these characters in the sequence data is valid. For example:

>>> seq = DNA('GCCRCCATGG', metadata={'name': 'Kozak consensus sequence'})
>>> seq
DNA
--------------------------------------
Metadata:
    'name': 'Kozak consensus sequence'
Stats:
    length: 10
    has gaps: False
    has degenerates: True
    has definites: True
    GC-content: 70.00%
--------------------------------------
0 GCCRCCATGG

One can customize grammared sequence types, for example to include modified amino
acids, methylation sites, or non-biological molecules. Refer to
:class:`GrammaredSequence` for instructions.

Each class contains a degenerate map for lookup:

>>> sorted(DNA.degenerate_map['R'])
['A', 'G']

Check if a sequence contains degenerate characters:

>>> seq.has_degenerates()
True

Identify the positions of degenerate characters:

>>> seq.degenerates()
array([False, False, False,  True, False, False, False, False, False,
       False], dtype=bool)

Likewise, one can create a sequence with gaps, and identify the gap positions.

>>> gapped = DNA('ATG--CAC')
>>> gapped.gaps()
array([False, False, False,  True,  True, False, False, False], dtype=bool)

One can efficiently remove all gaps from a sequence:

>>> ungapped = gapped.degap()
>>> print(ungapped)
ATGCAC


Exploring sequences
^^^^^^^^^^^^^^^^^^^

The primary information stored in a sequence object is the sequence data itself. It can
be retrieved from the string representation of the sequence:

>>> str(seq)
'GCCRCCATGG'

Or:

>>> print(seq)
GCCRCCATGG

Under the hood, the sequence data is stored as a NumPy array of bytes, which permit
compact storage and efficient, vectorized operations.

>>> seq.values
array([b'G', b'C', b'C', b'R', b'C', b'C', b'A', b'T', b'G', b'G'],
      dtype='|S1')

This array can also be viewed as ASCII code points:

>>> seq.values.view('uint8')
array([71, 67, 67, 82, 67, 67, 65, 84, 71, 71], dtype=uint8)

A sequence object shares typical behaviors of a Python string. This means that you can
apply Python string **indexing** and **slicing** approaches on a sequence. For example:

>>> seq[0]
DNA
--------------------------------------
Metadata:
    'name': 'Kozak consensus sequence'
Stats:
    length: 1
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 100.00%
--------------------------------------
0 G

>>> seq[2:5]
DNA
--------------------------------------
Metadata:
    'name': 'Kozak consensus sequence'
Stats:
    length: 3
    has gaps: False
    has degenerates: True
    has definites: True
    GC-content: 66.67%
--------------------------------------
0 CRC

The output of indexing or slicing is also a sequence object. This process doesn't
copy the original sequence data. Instead, it returns a **view** of the same memory
space, which is efficient. Refer to NumPy's `copies and views
<https://numpy.org/doc/stable/user/basics.copies.html>`_ for more information.

>>> seq[2:5].values.base is seq.values.base
True

Some Python string methods, such as ``count`` and ``index``, also work on scikit-bio
sequences.

>>> seq.count('C')
4

>>> seq.index('ATG')
6

To concatenate multiple sequences:

>>> seq = DNA.concat([DNA('GGATCC'), DNA('AAGCTT'), DNA('GAATTC')])
>>> print(seq)
GGATCCAAGCTTGAATTC


Sequence conversion
^^^^^^^^^^^^^^^^^^^

scikit-bio supports convenient conversion between DNA (two strands), RNA and protein
sequences. For example:

>>> dna = DNA('TGCTACATCCAGAACTGCCCCCTGGGA', metadata={'name': 'oxytocin'})
>>> dna
DNA
-------------------------------
Metadata:
    'name': 'oxytocin'
Stats:
    length: 27
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 59.26%
-------------------------------
0 TGCTACATCC AGAACTGCCC CCTGGGA

Reverse complement a DNA sequence:

>>> print(dna.reverse_complement())
TCCCAGGGGGCAGTTCTGGATGTAGCA

Transcribe a DNA sequence into mRNA:

>>> rna = dna.transcribe()
>>> print(rna)
UGCUACAUCCAGAACUGCCCCCUGGGA

Translate an mRNA sequence into protein:

>>> prot = rna.translate()
>>> print(prot)
CYIQNCPLG

By default, the standard genetic code table (1) is used. To use a specific table (see
the :class:`GeneticCode` class for details), simply specify its index:

>>> prot = rna.translate(11)  # Bacteria  # doctest: +SKIP

A DNA sequence can be directly translated into protein. The output is the same as
transcription followed by translation.

>>> prot_ = dna.translate()
>>> prot_ == prot
True


Searching for patterns
^^^^^^^^^^^^^^^^^^^^^^

Sequences can be searched for patterns, a.k.a., **motifs** in molecular biology. Except
for searching for exact substrings using the ``index`` method (see above), one can also
search for flexibly defined patterns using :wiki:`regular expressions
<Regular_expression>`.

>>> from skbio.sequence import RNA
>>> seq = RNA('AGG-GGACUGAA')
>>> for match in seq.find_with_regex('(G+AC[UT])'):
...     match
slice(4, 9, None)

The returned values are slices, with which one can extract the relevant subsequences.

>>> print(seq[match])
GGACU

scikit-bio has a convenience shortcut ``find_motifs``, which searches for contiguous
runs of purines or pyrimidines.

>>> for motif in seq.find_motifs('purine-run', min_length=2):
...     print(seq[motif])
AGG
GGA
GAA


Comparing sequences
^^^^^^^^^^^^^^^^^^^

It is straightforward to compute the distance between two sequences. By default, the
Hamming distance (:func:`~distance.hamming`, i.e., proportion of different sites) is
calculated, This metric is useful in analyses such as sequence clustering and
phylogenetic reconstruction.

>>> seq1 = RNA('GACCCGCUUU')
>>> seq2 = RNA('GCCCCCCUUU')
>>> seq1.distance(seq2)
0.2

Hamming distance requires that the two sequences have the same length. This often
requires that you **align** the sequences before comparison. See :mod:`skbio.alignment`
for how to align sequences.

Similarly, you can calculate the percent (dis)similarity between a pair of
aligned sequences.

>>> seq1 = RNA('ACCGUUAGUC')
>>> seq2 = RNA('ACGGGU--UC')
>>> seq1.match_frequency(seq2, relative=True)
0.6
>>> seq1.mismatch_frequency(seq2, relative=True)
0.4


.. _annotate_sequences:

Annotating sequences
^^^^^^^^^^^^^^^^^^^^

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
For example, the FASTA parser reads the sequence ID and description into the ``id``
and ``description`` fields:

>>> from skbio.sequence import Protein
>>> fasta = '''>NP_000591.1 interleukin-6 isoform 1 precursor [Homo sapiens]
... MNSFSTSAFGPVAFSLGLLLVLPAAFPAPVPPGEDSKDVAAPHRQPLTSSERIDKQIRYILDGISALRKE
... TCNKSNMCESSKEALAENNLNLPKMAEKDGCFQSGFNEETCLVKIITGLLEFEVYLEYLQNRFESSEEQA
... RAVQMSTKVLIQFLQKKAKNLDAITTPDPTTNASLLTKLQAQNQWLQDMTTHLILRSFKEFLQSSLRALR
... QM'''
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
genes and regulatory elements. See the :class:`~skbio.metadata.IntervalMetadata` class
for details. The following example annotates three regions of the human H4C6 gene.

>>> fasta = '''>NM_003540.4 Homo sapiens H4 clustered histone 6 (H4C6), mRNA
... GCAAAAGTTAAGAGTTGTTGTTTGTCTTCGATCATGTCTGGTAGAGGCAAAGGTGGTAAAGGTTTAGGAA
... AGGGAGGCGCCAAGCGCCATCGCAAAGTGCTGCGTGACAACATACAGGGCATCACGAAGCCCGCCATCCG
... TCGCTTGGCCCGACGCGGCGGCGTGAAACGCATTTCGGGCCTCATTTATGAGGAGACCCGCGGTGTTCTT
... AAGGTGTTCCTGGAGAATGTGATACGGGACGCCGTAACCTACACGGAGCACGCCAAGCGTAAGACAGTCA
... CTGCAATGGATGTTGTCTACGCGCTCAAGCGCCAGGGACGCACTCTGTACGGCTTTGGTGGCTGAGCCTC
... ACCCCGGCTTTTTATTTAACAGCTCACCCATAAAAGGCCCTTTTCAGGGCC'''
>>> seq = DNA.read([fasta], format='fasta')
>>> _ = seq.interval_metadata.add(bounds=[(0, 401)], metadata={'name': 'exon'})
>>> _ = seq.interval_metadata.add(bounds=[(33, 345)], metadata={'name': 'cds'})
>>> _ = seq.interval_metadata.add(bounds=[(385, 401)], metadata={'name': 'stem_loop'})

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
