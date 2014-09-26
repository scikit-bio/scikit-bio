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
        num_subseqs = len(data[labels[0]])
        if num_subseqs > 50:  # only check the first 50 subsequences
            num_subseqs = 50
        subseq_length = len(data[labels[0]][0])
        end_lengths = set()  # subsequence lengths at end of file
        for i in range(num_subseqs):
            for label in labels:
                seq = data[label][i]
                if len(seq) > subseq_length:
                    raise ClustalFormatError
                elif len(seq) < subseq_length:
                    end_lengths.add(len(seq))

        # All trailing subsequences must be the same
        if len(end_lengths) > 1:
            raise ClustalFormatError
    except ClustalFormatError:
        return False, {}
    return not empty, {}


@register_writer('clustal')
def _generator_to_clustal(records, fh):
    r"""writes aligned sequences to a specified file
    Parameters
    ----------
    record: iterator
        A generator of aligned sequences
    fh: open file handle object
        An open file handle object containing Clustal sequences.

    """
    clen = 60  # Max length of clustal lines
    names, seqs = zip(*[(s.id, s.sequence) for s in records])
    nameLen = max(map(len, names))
    seqLen = max(map(len, seqs))
    fh.write('CLUSTAL\n\n')
    for i in range(0, seqLen, clen):
        for label, seq in zip(names, seqs):
            name = ('{:<%d}' % (nameLen)).format(label)
            fh.write("%s\t%s\t\n" % (name, seq[i:i+clen]))
        fh.write("\n")


@register_reader('clustal')
def _clustal_to_generator(fh, strict=True):
    r"""yields labels and sequences from msa (multiple sequence alignment)

    Parameters
    ----------

    data : open file object
        An open Clustal file.
    strict : boolean
        Whether or not to raise a ``ClustalFormatError``
        when no labels are found.

    Returns
    -------

    label : str
        label of the sequence
    seq : str
        sequence for each label

    Notes
    -----

    Currently, does not check whether sequences are the same length and are in
    order. Skips any line that starts with a blank.

    ``_clustal_to_msa`` preserves the order of the sequences from the original
    file.  However, it does use a dict as an intermediate, so two sequences
    can't have the same label. This is probably OK since Clustal will refuse to
    run on a FASTA file in which two sequences have the same label, but could
    potentially cause trouble with manually edited files (all the segments of
    the conflicting sequences would be interleaved, possibly in an
    unpredictable way).

    If the lines have trailing numbers (i.e. Clustal was run with
    `-LINENOS=ON`), silently deletes them. Does not check that the numbers
    actually correspond to the number of chars in the sequence printed so far.


    Examples
    --------
    Assume we have a fasta formatted file with the following contents::

        CLUSTAL W (1.82) multiple sequence alignment

        abc   GCAUGCAUCUGCAUACGUACGUACGCAUGCAUCA 60
        def   ----------------------------------
        xyz   ----------------------------------

        abc   GUCGAUACAUACGUACGUCGUACGUACGU-CGAC 11
        def   ---------------CGCGAUGCAUGCAU-CGAU 18
        xyz   -----------CAUGCAUCGUACGUACGCAUGAC 23

    We can use the following code:

    >>> from StringIO import StringIO
    >>> from skbio.io import _clustal_to_msa
    >>> clustal_f = StringIO('abc   GCAUGCAUCUGCAUACGUACGUACGCAUGCA 60\n'
    ...                      'def   -------------------------------\n'
    ...                      'xyz   -------------------------------\n'
    ...                      '\n'
    ...                      'abc   GUCGAUACAUACGUACGUCGGUACGU-CGAC 11\n'
    ...                      'def   ---------------CGUGCAUGCAU-CGAU 18\n'
    ...                      'xyz   -----------CAUUCGUACGUACGCAUGAC 23\n')
    >>> for label, seq in parse_clustal(clustal_f):
    ...     print(label)
    ...     print(seq)
    abc
    GCAUGCAUCUGCAUACGUACGUACGCAUGCAGUCGAUACAUACGUACGUCGGUACGU-CGAC
    def
    ----------------------------------------------CGUGCAUGCAU-CGAU
    xyz
    ------------------------------------------CAUUCGUACGUACGCAUGAC


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

    for key in labels:
        yield BiologicalSequence(id=key, sequence=''.join(data[key]))
