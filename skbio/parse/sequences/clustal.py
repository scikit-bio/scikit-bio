#!/usr/bin/env python
from __future__ import division

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from skbio.core.exception import RecordError
from skbio.parse.record import DelimitedSplitter

from string import strip


def label_line_parser(record, splitter, strict=True):
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
                raise RecordError(
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


def is_clustal_seq_line(line):
    """Returns True if line starts with a non-blank character but not 'CLUSTAL'

    Useful for filtering other lines out of the file.
    """
    return line and (not line[0].isspace()) and\
        (not line.startswith('CLUSTAL')) and (not line.startswith('MUSCLE'))

last_space = DelimitedSplitter(None, -1)


def delete_trailing_number(line):
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


def parse_clustal(record, strict=True):
    r"""Returns a dictionary of the parsed data

    Parameters
    ----------

    data : open file object
        An open Clustal file.


    strict : boolean
        Whether or not to raise a ``RecordError`` when no labels are found.


    Returns
    -------

    labels: list
        list of strings with the labels

    data : dict
        a dict of label to sequence (pieces not joined)

    Notes
    -----

    Currently, does not check whether sequences are the same length and are in
    order. Skips any line that starts with a blank.

    ``parse_clustal`` preserves the order of the sequences from the original
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
    >>> from skbio.parse.sequences import parse_clustal
    >>> clustal_f = StringIO("abc   GCAUGCAUCUGCAUACGUACGUACGCAUGCA 60\n"
    ...                      'def   -------------------------------\n'
    ...                      'xyz   -------------------------------\n'
    ...                      '\n'
    ...                      'abc   GUCGAUACAUACGUACGUCGGUACGU-CGAC 11\n'
    ...                      'def   ---------------CGUGCAUGCAU-CGAU 18\n'
    ...                      'xyz   -----------CAUUCGUACGUACGCAUGAC 23\n')
    >>> data, labels = parse_clustal(clustal_f)
    >>> print labels
    ['abc', 'def', 'xyz']
    >>> for label, seq in data.iteritems():
    ...     print label
    ...     print seq
    xyz
    ['-------------------------------', '-----------CAUUCGUACGUACGCAUGAC']
    abc
    ['GCAUGCAUCUGCAUACGUACGUACGCAUGCA', 'GUCGAUACAUACGUACGUCGGUACGU-CGAC']
    def
    ['-------------------------------', '---------------CGUGCAUGCAU-CGAU']

    References
    ----------

    .. [1] Clustal: Multiple Sequence Alginment http://www.clustal.org

    """
    records = map(delete_trailing_number, filter(is_clustal_seq_line, record))
    return label_line_parser(records, last_space, strict)
