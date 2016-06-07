r"""
Clustal format (:mod:`skbio.io.format.clustal`)
===============================================

.. currentmodule:: skbio.io.format.clustal

Clustal format (``clustal``) stores multiple sequence alignments. This format
was originally introduced in the Clustal package [1]_.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.alignment.TabularMSA`                              |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
A clustal-formatted file is a plain text format. It can optionally have a
header, which states the clustal version number. This is followed by the
multiple sequence alignment, and optional information about the degree of
conservation at each position in the alignment [2]_.

Alignment Section
^^^^^^^^^^^^^^^^^
Each sequence in the alignment is divided into subsequences each at most
60 characters long. The sequence identifier for each sequence precedes each
subsequence. Each subsequence can optionally be followed by the cumulative
number of non-gap characters up to that point in the full sequence (not
included in the examples below). A line containing conservation information
about each position in the alignment can optionally follow all of the
subsequences (not included in the examples below).

.. note:: scikit-bio ignores conservation information when reading and does not
   support writing conservation information.

.. note:: When reading a clustal-formatted file into an
   ``skbio.alignment.TabularMSA`` object, sequence identifiers/labels are
   stored as ``TabularMSA`` index labels (``index`` property).

   When writing an ``skbio.alignment.TabularMSA`` object as a clustal-formatted
   file, ``TabularMSA`` index labels will be converted to strings and written
   as sequence identifiers/labels.

Format Parameters
-----------------
The only supported format parameter is ``constructor``, which specifies the
type of in-memory sequence object to read each aligned sequence into. This must
be a subclass of ``GrammaredSequence`` (e.g., ``DNA``, ``RNA``, ``Protein``)
and is a required format parameter. For example, if you know that the clustal
file you're reading contains DNA sequences, you would pass ``constructor=DNA``
to the reader call.

Examples
--------
Assume we have a clustal-formatted file of RNA sequences:

.. code-block:: none

    CLUSTAL W (1.82) multiple sequence alignment

    abc   GCAUGCAUCUGCAUACGUACGUACGCAUGCAUCA
    def   ----------------------------------
    xyz   ----------------------------------

    abc   GUCGAUACAUACGUACGUCGUACGUACGU-CGAC
    def   ---------------CGCGAUGCAUGCAU-CGAU
    xyz   -----------CAUGCAUCGUACGUACGCAUGAC

We can use the following code to read the clustal file into a ``TabularMSA``:

>>> from skbio import TabularMSA, RNA
>>> clustal_f = ['CLUSTAL W (1.82) multiple sequence alignment\n',
...              '\n',
...              'abc   GCAUGCAUCUGCAUACGUACGUACGCAUGCA\n',
...              'def   -------------------------------\n',
...              'xyz   -------------------------------\n',
...              '\n',
...              'abc   GUCGAUACAUACGUACGUCGGUACGU-CGAC\n',
...              'def   ---------------CGUGCAUGCAU-CGAU\n',
...              'xyz   -----------CAUUCGUACGUACGCAUGAC\n']
>>> msa = TabularMSA.read(clustal_f, constructor=RNA)
>>> msa
TabularMSA[RNA]
--------------------------------------------------------------
Stats:
    sequence count: 3
    position count: 62
--------------------------------------------------------------
GCAUGCAUCUGCAUACGUACGUACGCAUGCAGUCGAUACAUACGUACGUCGGUACGU-CGAC
----------------------------------------------CGUGCAUGCAU-CGAU
------------------------------------------CAUUCGUACGUACGCAUGAC
>>> msa.index
Index(['abc', 'def', 'xyz'], dtype='object')

We can use the following code to write a ``TabularMSA`` to a clustal-formatted
file:

>>> from io import StringIO
>>> from skbio import DNA
>>> seqs = [DNA('ACCGTTGTA-GTAGCT', metadata={'id': 'seq1'}),
...         DNA('A--GTCGAA-GTACCT', metadata={'id': 'sequence-2'}),
...         DNA('AGAGTTGAAGGTATCT', metadata={'id': '3'})]
>>> msa = TabularMSA(seqs, minter='id')
>>> msa
TabularMSA[DNA]
----------------------
Stats:
    sequence count: 3
    position count: 16
----------------------
ACCGTTGTA-GTAGCT
A--GTCGAA-GTACCT
AGAGTTGAAGGTATCT
>>> msa.index
Index(['seq1', 'sequence-2', '3'], dtype='object')
>>> fh = StringIO()
>>> _ = msa.write(fh, format='clustal')
>>> print(fh.getvalue()) # doctest: +NORMALIZE_WHITESPACE
CLUSTAL
<BLANKLINE>
<BLANKLINE>
seq1        ACCGTTGTA-GTAGCT
sequence-2  A--GTCGAA-GTACCT
3           AGAGTTGAAGGTATCT
<BLANKLINE>
<BLANKLINE>

References
----------
.. [1] http://www.sciencedirect.com/science/article/pii/0378111988903307
.. [2] http://web.mit.edu/meme_v4.9.0/doc/clustalw-format.html

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.io import create_format, ClustalFormatError
from skbio.alignment import TabularMSA


clustal = create_format('clustal')


def _label_line_parser(record):
    """Returns dict mapping list of data to labels, plus list with field order.

    Field order contains labels in order encountered in file.

    NOTE: doesn't care if lines are out of order in different blocks. This
    should never happen anyway, but it's possible that this behavior should
    be changed to tighten up validation.
    """
    labels = []
    result = {}
    for line in record:
        split_line = line.strip().rsplit(None, 1)

        if len(split_line) == 2:
            key, val = split_line
        else:
            raise ClustalFormatError(
                "Failed to parse sequence identifier and subsequence from "
                "the following line: %r" % line)

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


@clustal.sniffer()
def _clustal_sniffer(fh):
    # Strategy
    #   The following conditions preclude a file from being clustal
    #       * It is an empty file
    #       * The whole sequences have differing lengths
    #       * The sub-sequences have differing lengths
    #       * One of the sequence ids is not immediately
    #         followed by a subsequence
    empty = True
    if fh.read(7) != 'CLUSTAL':
        return False, {}
    fh.seek(0)
    try:
        records = map(_delete_trailing_number,
                      filter(_is_clustal_seq_line, fh))
        data, labels = _label_line_parser(records)
        if len(data) > 0:
            empty = False
        # Only check first 50 sequences
        aligned_correctly = _check_length(data, labels, 50)
        if not aligned_correctly:
            raise ClustalFormatError("Sequences not aligned properly")
    except ClustalFormatError:
        return False, {}
    return not empty, {}


@clustal.writer(TabularMSA)
def _tabular_msa_to_clustal(obj, fh):
    if not obj.index.is_unique:
        raise ClustalFormatError(
            "TabularMSA's index labels must be unique.")

    clen = 60  # Max length of clustal lines
    seqs = [str(s) for s in obj]
    names = [str(label) for label in obj.index]
    nameLen = max(map(len, names))
    seqLen = max(map(len, seqs))
    fh.write('CLUSTAL\n\n\n')
    for i in range(0, seqLen, clen):
        for label, seq in zip(names, seqs):
            name = ('{:<%d}' % (nameLen)).format(label)
            fh.write("%s\t%s\n" % (name, seq[i:i+clen]))
        fh.write("\n")


@clustal.reader(TabularMSA)
def _clustal_to_tabular_msa(fh, constructor=None):
    r"""yields labels and sequences from msa (multiple sequence alignment)

    Parameters
    ----------
    fh : open file object
        An open Clustal file.

    Returns
    -------
    skbio.TabularMSA
        MSA containing aligned sequences.

    Raises
    ------
    skbio.util.exception.ClustalFormatError
        If the sequences in `fh` don't have the same sequence length
        or if the sequence ids don't properly match with the subsequences

    Notes
    -----
    Skips any line that starts with a blank.

    ``_clustal_to_tabular_msa`` preserves the order of the sequences from the
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
    if constructor is None:
        raise ValueError("Must provide `constructor`.")

    records = map(_delete_trailing_number,
                  filter(_is_clustal_seq_line, fh))
    data, labels = _label_line_parser(records)

    aligned_correctly = _check_length(data, labels)
    if not aligned_correctly:
        raise ClustalFormatError("Sequences not aligned properly")
    seqs = []
    for label in labels:
        seqs.append(constructor(''.join(data[label])))
    return TabularMSA(seqs, index=labels)
