r"""Sequence Alignments (:mod:`skbio.alignment`)
============================================

.. currentmodule:: skbio.alignment

This module provides functionality for computing and manipulating sequence
alignments. DNA, RNA, and protein sequences can be aligned, as well as
sequences with or without custom alphabets.

See the |alignment_tutorial|_ section for performing and working with sequence
alignments using scikit-bio.


Alignment structures
--------------------

Tabular format of aligned sequences.

.. autosummary::
   :toctree: generated/

   TabularMSA

Compact format of alignment paths.

.. autosummary::
   :toctree: generated/

   AlignPath
   PairAlignPath


Pairwise alignment
------------------

A versatile, efficient and generalizable pairwise alignment algorithm.

.. autosummary::
   :toctree: generated/

    pair_align

Convenience wrappers with preset scoring schemes.

.. autosummary::
   :toctree: generated/

    pair_align_nucl
    pair_align_prot


Alignment statistics
--------------------

.. autosummary::
   :toctree: generated/

    align_score
    align_dists


Deprecated functionality
------------------------

Pure Python algorithms (slow; educational-purposes only)

.. autosummary::
   :toctree: generated/

   global_pairwise_align_nucleotide
   global_pairwise_align_protein
   global_pairwise_align
   local_pairwise_align_nucleotide
   local_pairwise_align_protein
   local_pairwise_align


.. |alignment_tutorial| replace:: **Tutorial**
.. _alignment_tutorial:

Tutorial
--------

Sequence alignment is the process of arranging biological sequences such that they
resemble each other when positioned together. The alignment reveals the degree of
similarity between sequences, which may indicate evolutionary relationships (i.e.,
"homology") and/or structural and functional relevance. The alignment process permits
only one operation: shifting characters by inserting gaps (`-`) between them. Shuffling
characters is not permitted.


Quick start
^^^^^^^^^^^

To demonstrate this process, let's create two DNA sequences.

>>> from skbio.sequence import DNA
>>> seq1 = DNA('ACTACCAGATTACTTACGGATCAGGTACTTGCCAACAA')
>>> seq2 = DNA('CGAAACTACTAGATTACGGATCTTACTTTCCAGCAAGG')

Just by eyeballing, we can observe some local similarities between the sequences,
particularly in the middle sections. However, the beginning and ending portions aren't
so similar between the two.

Apply :func:`pair_align_nucl` to perform pairwise alignment of the two nucleotide
sequences.

>>> from skbio.alignment import pair_align_nucl
>>> aln = pair_align_nucl(seq1, seq2)
>>> aln  # doctest: +ELLIPSIS
PairAlignResult(score=22.0, paths=[<PairAlignPath, positions: 44, segments: 7, ...

The result contains ``score`` and ``paths``. The latter is a list of optimal alignment
paths calculated by the algorithm. By default, only one path will be returned, which is
all we need in most applications.

>>> aln.paths[0]
<PairAlignPath, positions: 44, segments: 7, CIGAR: '4I13M4D6M2D13M2I'>

To obtain aligned sequences, we can call the :meth:`~AlignPath.to_aligned` method of
the path.

>>> aln.paths[0].to_aligned((seq1, seq2))
['----ACTACCAGATTACTTACGGATCAGGTACTTGCCAACAA--',
 'CGAAACTACTAGATTAC----GGATCT--TACTTTCCAGCAAGG']

As we can see, by inserting four consecutive gaps in the appropriate positions, the two
sequences appear more similar, or "aligned", with one another.

The score measures the goodness of the alignment. The value is dependent on the
sequences and the alignment parameters (discussed below).

>>> aln.score
22.0

Likewise, we can align protein sequences with :func:`pair_align_prot`:

>>> from skbio.sequence import Protein
>>> from skbio.alignment import pair_align_prot
>>> pair_align_prot(Protein("HEAGAWGHEE"), Protein("PAWHEAE"))  # doctest: +ELLIPSIS
PairAlignResult(score=15.0, paths=[<PairAlignPath, positions: 13, segments: 3, ...


Alignment algorithm
^^^^^^^^^^^^^^^^^^^

The alignment algorithm aims to find the optimal alignment path(s) that yield the
highest possible alignment score for two sequences. For decades,
:wiki:`dynamic programming <Dynamic_programming>` (DP) has been the gold standard for
pairwise sequence alignment, and is widely covered in many bioinformatics textbooks.
scikit-bio also implements this algorithm. ``pair_align_nucl`` and ``pair_align_prot``
are but convenience wrappers for the function :func:`pair_align`, which offers multiple
customizable parameters and comprehensive documentation explaining the algorithm and
everything you need to know for using it.

>>> from skbio.alignment import pair_align

By referring to the documentation of ``pair_align``, we can customize the algorithm's
behavior to best suit our research objectives. The default mode is *global* alignment,
as in the Needleman-Wunsch algorithm (which you may recall from textbooks), with both
ends free of gap penalty (a.k.a., "semi-global alignment", or "overlap alignment", in
some literature).

Next, we will switch to *local* alignment, as in the Smith-Waterman algorithm.

>>> aln = pair_align(seq1, seq2, mode="local")
>>> aln  # doctest: +ELLIPSIS
PairAlignResult(score=13.0, paths=[<PairAlignPath, positions: 25, segments: 3, ...

>>> aln.paths[0].to_aligned((seq1, seq2))
['TTACGGATCAGGTACTTGCCAACAA', 'TTACGGATCT--TACTTTCCAGCAA']

As shown, the alignment is shorter than that of the global alignment, and truncated:
only the highly similar regions in the middle of the two sequences are preserved, while
the less similar regions at both ends are discarded.

We can further tweak the alignment parameters. The wrapper function ``pair_align_nucl``
uses parameters that are consistent with NCBI BLASTN's defaults. To replicate this
behavior using ``pair_align``, we can proceed as follows:

>>> aln = pair_align(seq1, seq2, mode="global", sub_score=(2, -3), gap_cost=(5, 2))
>>> aln  # doctest: +ELLIPSIS
PairAlignResult(score=22.0, paths=[<PairAlignPath, positions: 44, segments: 7, ...

Here, 2 and -3 in ``sub_score`` represent match and mismatch scores (added to the
alignment score). 5 and 2 in ``gap_cost`` represent gap opening and gap extension
penalties (subtracted from the alignment score). The use of two gap penalty values is
known to as *affine* gap penalty.

Let's try an alternative setting. The follow parameters uses a nucleotide substitution
matrix "NUC.4.4" to replace match/mismatch scores. This provides more flexibility in
defining the score between pairs of characters. For example, when there are degenerate
characters in the DNA sequences, such as "R" and "Y", the substitution matrix usually
models the substitution scores more realistically. In addition, this setting uses a
simple *linear* gap penalty scheme to replace the more sophisticated affine gap penalty
scheme.

>>> aln = pair_align(seq1, seq2, mode="global", sub_score="NUC.4.4", gap_cost=3)
>>> aln  # doctest: +ELLIPSIS
PairAlignResult(score=106.0, paths=[<PairAlignPath, positions: 44, segments: 7, ...

Scikit-bio also provides an :func:`align_score` function to calculate the alignment
score of an existing alignment. You will need to supply the same parameters to get the
same outcome.

>>> from skbio.alignment import align_score
>>> align_score((aln.paths[0], (seq1, seq2)), sub_score="NUC.4.4", gap_cost=3)
106.0

Pairwise alignment is often taught in bioinformatics courses. By adding
``keep_matrices=True`` to the function, the alignment matrix calculated during the DP
process will be retained. This can be particularly useful for educational purposes.

>>> score, (path,), (matrix,) = pair_align(
...     "GATCT", "GTAC", mode="global", free_ends=False, keep_matrices=True)
>>> matrix  # doctest: +ELLIPSIS
array([[  0.,  -2.,  -4.,  -6.,  -8.],
       [ -2.,   1.,  -1.,  -3.,  -5.],
       [ -4.,  -1.,   0.,   0.,  -2.],
       [ -6.,  -3.,   0.,  -1.,  -1.],
       [ -8.,  -5.,  -2.,  -1.,   0.],
       [-10.,  -7.,  -4.,  -3.,  -2.]]...


Alignment path
^^^^^^^^^^^^^^

A pairwise alignment path is an instance of the :class:`PairAlignPath` class, which in
turn is a subclass of the :class:`AlignPath` class. This class implements an efficient
data structure to store the "path", i.e., the events of placing gaps inside sequences.
It uses run-length encoding (RLE) under the hood, and does not store the sequence data,
such that its memory cost is minimum. Refer to the class documentation for details of
this data structure.

>>> path = aln.paths[0]
>>> path.lengths
array([ 4, 13,  4,  6,  2, 13,  2])

>>> path.states
array([[1, 0, 2, 0, 2, 0, 1]], dtype=uint8)

The path can be represented by a CIGAR string, which is a common format used by many
high-throughput pairwise alignment tools.

>>> path.to_cigar()
'4I13M4D6M2D13M2I'

The letters in the CIGAR string represent the gap-placing events: ``I`` (insertion):
gap in sequence 1, ``D`` (deletion): gap in sequence 2, ``M`` (mis/match): gap in
neither sequence. The numbers represent the length of each segment. As you may have
noticed, the ``PairAlignPath`` data structure resembles the CIGAR encoding. In this
example, the alignment path contains four consecutive gaps flanking four chunks of
(mis)matched characters.

If we pass the sequences to the :meth:`~PairAlignPath.to_cigar` method, we will get
``=`` (match, i.e., same character) and ``X`` (mismatch, i.e., different characters)
instead of ``M``.

>>> path.to_cigar((seq1, seq2))
'4I5=1X7=4D5=1X2D5=1X3=1X3=2I'

One can construct a ``PairAlignPath`` from a CIGAR string using the
:meth:`~PairAlignPath.from_cigar` method. Likewise, this process doesn't rely on the
sequence data. Therefore, it is suitable for handling many CIGAR strings from large
data files, such as those in SAM/BAM formatted files.

>>> from skbio.alignment import PairAlignPath
>>> path = PairAlignPath.from_cigar('1I8M2D5M2I')
>>> path
<PairAlignPath, positions: 18, segments: 5, CIGAR: '1I8M2D5M2I'>

This method can also be used to parse CIGAR-formatted alignment results of external
sequence alignment tools, such as **Parasail**.

.. code-block:: python

    import parasail
    res = parasail.nw_trace(...)
    path = PairAlignPath.from_cigar(res.cigar.decode)

``AlignPath`` is a more generalized data structure that supports an arbitrary number
of sequences. For example, we can directly construct a path of three aligned sequences
using the :meth:`~AlignPath.from_aligned` method:

>>> from skbio.alignment import AlignPath
>>> path = AlignPath.from_aligned([
...     'CGTCGTGC',
...     'CA--GT-C',
...     'CGTCGT-T',
... ])
>>> path
<AlignPath, sequences: 3, positions: 8, segments: 5>

>>> path.lengths
array([2, 2, 2, 1, 1])

>>> path.states[0]
array([0, 2, 0, 6, 0], dtype=uint8)

The ``AlignPath`` (as well as ``PairAlignPath``) class serves as a central hub for
interacting with different alignment formats and tools. For example, the
:meth:`~AlignPath.from_coordinates` method can parse **BioPython**'s alignment results
(and vice versa with :meth:`~AlignPath.to_coordinates`).

.. code-block:: python

    from Bio import Align
    res = Align.PairwiseAligner().align(...)
    path = AlignPath.from_coordinates(res[0].coordinates)

The :meth:`~AlignPath.from_indices` method can parse **Biotite**'s alignment results
(and vice versa with :meth:`~AlignPath.to_indices`).

.. code-block:: python

    from biotite.sequence.align import align_optimal
    res = align_optimal(...)
    path = AlignPath.from_indices(res[0].trace.T)


Tabular alignment
^^^^^^^^^^^^^^^^^

In addition to alignment paths, scikit-bio provides :class:`TabularMSA`, a tabular data
structure that hosts aligned sequences ("MSA" stands for multiple sequence alignment).

One can convert an alignment path and the original (unaligned) sequences into a
``TabularMSA`` object using the :meth:`~TabularMSA.from_path_seqs` method (and vice
versa with :meth:`AlignPath.from_tabular`).

>>> from skbio.alignment import TabularMSA
>>> score, (path,), _ = pair_align_nucl(seq1, seq2)
>>> msa = TabularMSA.from_path_seqs(path, (seq1, seq2))
>>> msa
TabularMSA[DNA]
--------------------------------------------
Stats:
    sequence count: 2
    position count: 44
--------------------------------------------
----ACTACCAGATTACTTACGGATCAGGTACTTGCCAACAA--
CGAAACTACTAGATTAC----GGATCT--TACTTTCCAGCAAGG

Alternatively, one can read from a multi-FASTA file of aligned sequences:

.. code-block:: python

    msa = TabularMSA.read('input.fasta', constructor=DNA)

Or directly from aligned sequence objects:

>>> msa = TabularMSA([
...     DNA('CGTCGTGC'),
...     DNA('CA--GT-C'),
...     DNA('CGTCGT-T'),
... ])
>>> msa
TabularMSA[DNA]
---------------------
Stats:
    sequence count: 3
    position count: 8
---------------------
CGTCGTGC
CA--GT-C
CGTCGT-T

A ``TabularMSA`` object is iterally like a table. If you are also familiar with pandas
DataFrame, you can apply similar indexing methods to it.

>>> msa.index
RangeIndex(start=0, stop=3, step=1)

>>> msa.iloc[0]
DNA
--------------------------
Stats:
    length: 8
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 75.00%
--------------------------
0 CGTCGTGC

>>> msa.iloc[:, 1]
Sequence
-------------
Stats:
    length: 3
-------------
0 GAG

Like :class:`~skbio.sequence.Sequence`, ``TabularMSA`` can be annotated with metadata.

>>> msa.metadata = {'name': 'test alignment'}
>>> msa.positional_metadata = {
...     'identity': [80, 85, 88, 90, 90, 95, 92, 86]}

The class provides some utilities for downstream analysis. For example, to obtain a
consensus sequence of the alignment:

>>> msa.consensus()
DNA
------------------------------
Positional metadata:
    'identity': <dtype: int64>
Stats:
    length: 8
    has gaps: True
    has degenerates: False
    has definites: True
    GC-content: 71.43%
------------------------------
0 CGTCGT-C

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._tabular_msa import TabularMSA
from ._pairwise import (
    local_pairwise_align_nucleotide,
    local_pairwise_align_protein,
    local_pairwise_align,
    global_pairwise_align_nucleotide,
    global_pairwise_align_protein,
    global_pairwise_align,
)
from skbio.alignment._path import AlignPath, PairAlignPath
from skbio.alignment._score import align_score
from skbio.alignment._distance import align_dists
from skbio.alignment._pair import pair_align, pair_align_nucl, pair_align_prot

__all__ = [
    "TabularMSA",
    "AlignPath",
    "PairAlignPath",
    "global_pairwise_align",
    "global_pairwise_align_nucleotide",
    "global_pairwise_align_protein",
    "local_pairwise_align",
    "local_pairwise_align_nucleotide",
    "local_pairwise_align_protein",
    "align_score",
    "align_dists",
    "pair_align",
    "pair_align_nucl",
    "pair_align_prot",
]
