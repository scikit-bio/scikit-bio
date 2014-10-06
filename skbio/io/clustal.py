r"""
Clustal format (:mod:`skbio.io.clustal`)
========================================

.. currentmodule:: skbio.io.clustal

Clustal format (``clustal``) stores multiple alignment sequences.
This format was originally introduced in Desmond Higgens in the
Clustal package [1]_.

Format Support
--------------
**Has Sniffer: Yes**

+----------+----------+------------------------------------------------------+
|**Reader**|**Writer**|                   **Object Class**                   |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.alignment.Alignment`                      |
+----------+----------+------------------------------------------------------+

Format Specification
--------------------
A Clustal file is a plain text format.  It can optionally have a header,
which states the Clustal version number. This is followed by the
multiple sequence alignment, and optionally conserved sequence information [2]_

Alignment Section
^^^^^^^^^^^^^^^^^
Each sequence in the alignment is divided into subsequences each at most
60 characters long. The sequence ID of the alignment precedes each subsequence.
Optionally, a line containing conservation information about each position
in the multiple alignment.  This information is specified as follows

         *  -- all residues or nucleotides in that column are identical
         :  -- conserved substitutions have been observed
         .  -- semi-conserved substitutions have been observed
            -- no match.


.. note:: scikit-bio does not support writing conservation information for
          multiple alignments
   :: scikit-bio will only write a Clustal-formatted file if the alignment's
   sequence characters are valid IUPAC characters, as defined in
   :mod:`skbio.sequence`. The specific lexicon that is validated against
   depends on the type of sequences stored in the alignment.


Examples
--------

Assume we have a clustal formatted file with the following contents::

    CLUSTAL W (1.82) multiple sequence alignment

    abc   GCAUGCAUCUGCAUACGUACGUACGCAUGCAUCA 60
    def   ----------------------------------
    xyz   ----------------------------------

    abc   GUCGAUACAUACGUACGUCGUACGUACGU-CGAC 11
    def   ---------------CGCGAUGCAUGCAU-CGAU 18
    xyz   -----------CAUGCAUCGUACGUACGCAUGAC 23

We can use the following code to read a clustal file:

>>> from StringIO import StringIO
>>> from skbio.io import read
>>> clustal_f = StringIO('abc   GCAUGCAUCUGCAUACGUACGUACGCAUGCA 60\n'
...                      'def   -------------------------------\n'
...                      'xyz   -------------------------------\n'
...                      '\n'
...                      'abc   GUCGAUACAUACGUACGUCGGUACGU-CGAC 11\n'
...                      'def   ---------------CGUGCAUGCAU-CGAU 18\n'
...                      'xyz   -----------CAUUCGUACGUACGCAUGAC 23\n')
>>> for dna in read(clustal_f,format="clustal", into=Alignment):
...     print(dna.id)
...     print(dna.sequence)
abc
GCAUGCAUCUGCAUACGUACGUACGCAUGCAGUCGAUACAUACGUACGUCGGUACGU-CGAC
def
----------------------------------------------CGUGCAUGCAU-CGAU
xyz
------------------------------------------CAUUCGUACGUACGCAUGAC

We can use the following code to write to a file

>>> from skbio import Alignment, DNA
>>> from skbio.io import write
>>> seqs = [DNA('ACCGTTGTA-GTAGCT', id='seq1'),
...         DNA('A--GTCGAA-GTACCT', id='sequence-2'),
...         DNA('AGAGTTGAAGGTATCT', id='3')]
>>> aln = Alignment(seqs)
>>> from StringIO import StringIO
>>> fh = StringIO()
>>> aln.write(fh, format='clustal')
>>> print(fh.getvalue()) # doctest: +NORMALIZE_WHITESPACE
CLUSTAL
<BLANKLINE>
<BLANKLINE>
seq1        ACCGTTGTA-GTAGCT
sequence-2  A--GTCGAA-GTACCT
3           AGAGTTGAAGGTATCT
<BLANKLINE>
<BLANKLINE>



[1]_ http://www.ncbi.nlm.nih.gov/pubmed/3243435
[2]_ http://web.mit.edu/meme_v4.9.0/doc/clustalw-format.html
"""


# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function

from skbio.parse.record import DelimitedSplitter
from skbio.io import (register_reader, register_writer, register_sniffer,
                      ClustalFormatError)
from skbio.sequence import BiologicalSequence
from skbio.alignment import Alignment


def _label_line_parser(record, splitter, strict=True):
    """Returns dict mapping list of data to labels, plus list with field order.

    Field order contains labels in order encountered in file.

    NOTE: doesn't care if lines are out of order in different blocks. This
    should never happen anyway, but it's possible that this behavior should
    be changed to tighten up validation.
    """
    labels = []
    result = {}
    for line in record:
        try:
            key, val = splitter(line.rstrip())
        except:

            if strict:
                raise ClustalFormatError(
                    "Failed to extract key and value from line %s" %
                    line)
            else:
                continue  # just skip the line if not strict

        if key in result:
            result[key].append(val)
        else:
            result[key] = [val]
            labels.append(key)
    return result, labels


def _is_clustal_seq_line(line):
    """Returns True if line starts with a non-blank character but not 'CLUSTAL'

    Useful for filtering other lines out of the file.
    """
    return line and (not line[0].isspace()) and\
        (not line.startswith('CLUSTAL')) and (not line.startswith('MUSCLE'))

last_space = DelimitedSplitter(None, -1)


def _delete_trailing_number(line):
    """Deletes trailing number from a line.

    WARNING: does not preserve internal whitespace when a number is removed!
    (converts each whitespace run to a single space). Returns the original
    line if it didn't end in a number.
    """
    pieces = line.split()
    try:
        int(pieces[-1])
        return ' '.join(pieces[:-1])
    except ValueError:  # no trailing numbers
        return line


def _check_length(data, labels, num_seqs_check=None):
    """
    Checks the lengths of the clustal sequences to make
    sure that they are lining up right

    num_seqs_check: The number of sequences to check

    Return True if all of the subsequence lengths are equal
                or if data is empty
    Return False if one of the subsequence lengths differs
    """
    if len(labels) == 0:
        return True
    num_subseqs = len(data[labels[0]])
    if num_seqs_check is None:
        num_seqs_check = num_subseqs
    else:
        if num_seqs_check > num_subseqs:
            num_seqs_check = num_subseqs

    subseq_length = len(data[labels[0]][0])

    end_lengths = set()  # subsequence lengths at end of file
    for i in range(num_seqs_check):
        for label in labels:
            seq = data[label][i]
            if len(seq) > subseq_length:
                return False
            elif i+1 == num_subseqs:  # Last subsequence
                end_lengths.add(len(seq))
            elif len(seq) < subseq_length:
                return False
    # All trailing subsequences must be the same
    if len(end_lengths) > 1:
        return False
    return True


@register_sniffer("clustal")
def _clustal_sniffer(fh):
    # Strategy
    #   The following conditions preclude a file from being clustal
    #       * It is an empty file
    #       * The whole sequences have differing lengths
    #       * The sub-sequences have differing lengths
    #       * One of the sequence ids is not immediately
    #         followed by a subsequence
    empty = True
    try:
        records = map(_delete_trailing_number,
                      filter(_is_clustal_seq_line, fh))
        data, labels = _label_line_parser(records, last_space, strict=True)
        empty = False
        # Only check first 50 sequences
        aligned_correctly = _check_length(data, labels, 50)
        if not aligned_correctly:
            raise ClustalFormatError("Sequences not aligned properly")
    except ClustalFormatError:
        return False, {}
    return not empty, {}


@register_writer('clustal', Alignment)
def _alignment_to_clustal(obj, fh):
    r"""writes aligned sequences to a specified file
    Parameters
    ----------
    obj: Alignment object
        An alignment object containing a set of BiologicalSequence objects
    fh: open file handle object
        An open file handle object containing Clustal sequences.

    """
    clen = 60  # Max length of clustal lines
    names, seqs = zip(*[(s.id, s.sequence) for s in obj])
    nameLen = max(map(len, names))
    seqLen = max(map(len, seqs))
    fh.write('CLUSTAL\n\n\n')
    for i in range(0, seqLen, clen):
        for label, seq in zip(names, seqs):
            name = ('{:<%d}' % (nameLen)).format(label)
            fh.write("%s\t%s\n" % (name, seq[i:i+clen]))
        fh.write("\n")


@register_reader('clustal', Alignment)
def _clustal_to_alignment(fh, strict=True):
    r"""yields labels and sequences from msa (multiple sequence alignment)

    Parameters
    ----------

    fh : open file object
        An open Clustal file.
    strict : boolean
        Whether or not to raise a ``ClustalFormatError``
        when no labels are found.

    Returns
    -------
    skbio.Alignment
        Alignment object containing aligned biogical sequences

    Raises
    ------
        skbio.util.exception.ClustalFormatError
            If the sequences in `fh` don't have the same sequence length
            or if the sequence ids don't properly match with the subsequences
    Notes
    -----

    Skips any line that starts with a blank.

    ``_clustal_to_generator`` preserves the order of the sequences from the
    original file.  However, it does use a dict as an intermediate, so
    two sequences can't have the same label. This is probably OK since
    Clustal will refuse to run on a FASTA file in which two sequences have
    the same label, but could potentially cause trouble with manually
    edited files (all the segments of the conflicting sequences would
    be interleaved, possibly in an unpredictable way).

    If the lines have trailing numbers (i.e. Clustal was run with
    `-LINENOS=ON`), silently deletes them. Does not check that the numbers
    actually correspond to the number of chars in the sequence printed so far.

    References
    ----------
    .. [1] Thompson JD, Higgins DG, Gibson TJ,  "CLUSTAL W: improving the
        sensitivity of progressive multiple sequence alignment through sequence
        weighting, position-specific gap penalties and weight matrix choice.
        Thompson", Nucleic Acids Res. 1994 Nov 11;22(22):4673-80.

    """

    records = map(_delete_trailing_number,
                  filter(_is_clustal_seq_line, fh))
    data, labels = _label_line_parser(records, last_space, strict)

    aligned_correctly = _check_length(data, labels)
    if not aligned_correctly:
        raise ClustalFormatError("Sequences not aligned properly")
    alns = []
    for key in labels:
        alns.append(BiologicalSequence(id=key, sequence=''.join(data[key])))
    return Alignment(alns)
