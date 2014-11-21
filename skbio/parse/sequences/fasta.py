# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function

import warnings

import numpy as np

from skbio.io import RecordError
from skbio.parse.record_finder import LabeledRecordFinder


def is_fasta_label(x):
    """Checks if x looks like a FASTA label line."""
    return x.startswith('>')


def is_blank_or_comment(x):
    """Checks if x is blank or a FASTA comment line."""
    return (not x) or x.startswith('#') or x.isspace()


FastaFinder = LabeledRecordFinder(is_fasta_label, ignore=is_blank_or_comment)


def parse_fasta(infile, strict=True, label_to_name=None, finder=FastaFinder,
                label_characters='>', ignore_comment=False):
    r"""Generator of labels and sequences from a fasta file.

    .. note:: Deprecated in scikit-bio 0.2.0-dev
       ``parse_fasta`` will be removed in scikit-bio 0.3.0. It is replaced by
       ``read``, which is a more general method for deserializing
       FASTA-formatted files. ``read`` supports multiple file formats,
       automatic file format detection, etc. by taking advantage of
       scikit-bio's I/O registry system. See :mod:`skbio.io` for more details.

    Parameters
    ----------
    infile : open file object or str
        An open fasta file or a path to a fasta file.

    strict : bool
        If ``True`` a ``RecordError`` will be raised if there is a fasta label
        line with no associated sequence, or a sequence with no associated
        label line (in other words, if there is a partial record). If
        ``False``, partial records will be skipped.

    label_to_name : function
        A function to apply to the sequence label (i.e., text on the header
        line) before yielding it. By default, the sequence label is returned
        with no processing. This function must take a single string as input
        and return a single string as output.

    finder : function
        The function to apply to find records in the fasta file. In general
        you should not have to change this.

    label_characters : str
        String used to indicate the beginning of a new record. In general you
        should not have to change this.

    ignore_comment : bool
        If `True`, split the sequence label on spaces, and return the label
        only as the first space separated field (i.e., the sequence
        identifier). Note: if both ``ignore_comment`` and ``label_to_name`` are
        passed, ``ignore_comment`` is ignored (both operate on the label, so
        there is potential for things to get messy otherwise).

    Returns
    -------
    two-item tuple of str
        yields the label and sequence for each entry.

    Raises
    ------
    RecordError
        If ``strict == True``, raises a ``RecordError`` if there is a fasta
        label line with no associated sequence, or a sequence with no
        associated label line (in other words, if there is a partial record).

    Examples
    --------
    Assume we have a fasta-formatted file with the following contents::

        >seq1 db-accession-149855
        CGATGTCGATCGATCGATCGATCAG
        >seq2 db-accession-34989
        CATCGATCGATCGATGCATGCATGCATG

    >>> from StringIO import StringIO
    >>> fasta_f = StringIO('>seq1 db-accession-149855\n'
    ...                    'CGATGTCGATCGATCGATCGATCAG\n'
    ...                    '>seq2 db-accession-34989\n'
    ...                    'CATCGATCGATCGATGCATGCATGCATG\n')

    We can parse this as follows:

    >>> from skbio.parse.sequences import parse_fasta
    >>> for label, seq in parse_fasta(fasta_f):
    ...     print(label, seq)
    seq1 db-accession-149855 CGATGTCGATCGATCGATCGATCAG
    seq2 db-accession-34989 CATCGATCGATCGATGCATGCATGCATG

    The sequence label or header line in a fasta file is defined as containing
    two separate pieces of information, delimited by a space. The first space-
    separated entry is the sequence identifier, and everything following the
    first space is considered additional information (e.g., comments about the
    source of the sequence or the molecule that it encodes). Often we don't
    care about that information within our code. If you want to just return the
    sequence identifier from that line, you can pass ``ignore_comment=True``:

    >>> from StringIO import StringIO
    >>> fasta_f = StringIO('>seq1 db-accession-149855\n'
    ...                    'CGATGTCGATCGATCGATCGATCAG\n'
    ...                    '>seq2 db-accession-34989\n'
    ...                    'CATCGATCGATCGATGCATGCATGCATG\n')

    >>> from skbio.parse.sequences import parse_fasta
    >>> for label, seq in parse_fasta(fasta_f, ignore_comment=True):
    ...     print(label, seq)
    seq1 CGATGTCGATCGATCGATCGATCAG
    seq2 CATCGATCGATCGATGCATGCATGCATG

    """
    warnings.warn(
        "`parse_fasta` is deprecated and will be removed in scikit-bio 0.3.0. "
        "Please update your code to use `skbio.io.read(fh, format='fasta')` "
        "to obtain a generator of `BiologicalSequence` objects (or "
        "subclasses, see the `constructor` parameter).", DeprecationWarning)

    for rec in finder(infile):
        # first line must be a label line
        if not rec[0][0] in label_characters:
            if strict:
                raise RecordError(
                    "Found Fasta record without label line: %s" % rec)
            else:
                continue
        # record must have at least one sequence
        if len(rec) < 2:
            if strict:
                raise RecordError(
                    "Found label line without sequences: %s" % rec)
            else:
                continue

        # remove the label character from the beginning of the label
        label = rec[0][1:].strip()
        # if the user passed a label_to_name function, apply that to the label
        if label_to_name is not None:
            label = label_to_name(label)
        # otherwise, if the user passed ignore_comment, split the label on
        # spaces, and return the first space separated field (i.e., the
        # sequence identifier)
        elif ignore_comment:
            label = label.split()[0]
        else:
            pass

        # join the sequence lines into a single string
        seq = ''.join(rec[1:])

        yield label, seq


def parse_qual(infile, full_header=False):
    r"""yields label and qual from a qual file.

    .. note:: Deprecated in scikit-bio 0.2.0-dev
       ``parse_qual`` will be removed in scikit-bio 0.3.0. It is replaced by
       ``read``, which is a more general method for deserializing
       FASTA/QUAL-formatted files. ``read`` supports multiple file formats,
       automatic file format detection, etc. by taking advantage of
       scikit-bio's I/O registry system. See :mod:`skbio.io` for more details.

    Parameters
    ----------
    infile : open file object or str
        An open fasta file or path to it.

    full_header : bool
        Return the full header or just the id

    Returns
    -------
    label : str
        The quality label
    qual : array
        The quality at each position

    Examples
    --------
    Assume we have a qual formatted file with the following contents::

        >seq1
        10 20 30 40
        >seq2
        1 2 3 4

    >>> from StringIO import StringIO
    >>> from skbio.parse.sequences import parse_qual
    >>> qual_f = StringIO('>seq1\n'
    ...                   '10 20 30 40\n'
    ...                   '>seq2\n'
    ...                   '1 2 3 4\n')
    >>> for label, qual in parse_qual(qual_f):
    ...     print(label)
    ...     print(qual)
    seq1
    [10 20 30 40]
    seq2
    [1 2 3 4]

    """
    warnings.warn(
        "`parse_qual` is deprecated and will be removed in scikit-bio 0.3.0. "
        "Please update your code to use "
        "`skbio.io.read(fasta_fh, qual=qual_fh, format='fasta')` to obtain a "
        "generator of `BiologicalSequence` objects (or subclasses, see the "
        "`constructor` parameter) with quality scores.", DeprecationWarning)

    for rec in FastaFinder(infile):
        curr_id = rec[0][1:]
        curr_qual = ' '.join(rec[1:])
        try:
            parts = np.asarray(curr_qual.split(), dtype=int)
        except ValueError:
            raise RecordError(
                "Invalid qual file. Check the format of the qual file: each "
                "quality score must be convertible to an integer.")
        if full_header:
            curr_pid = curr_id
        else:
            curr_pid = curr_id.split()[0]
        yield (curr_pid, parts)
