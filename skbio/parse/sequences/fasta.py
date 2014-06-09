#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function

import numpy as np

from skbio.core.exception import RecordError
from skbio.parse.record_finder import LabeledRecordFinder


def is_fasta_label(x):
    """Checks if x looks like a FASTA label line."""
    return x.startswith('>')


def is_blank_or_comment(x):
    """Checks if x is blank or a FASTA comment line."""
    return (not x) or x.startswith('#') or x.isspace()


def is_blank(x):
    """Checks if x is blank."""
    return (not x) or x.isspace()


FastaFinder = LabeledRecordFinder(is_fasta_label, ignore=is_blank_or_comment)


def parse_fasta(infile, strict=True, label_to_name=str, finder=FastaFinder,
                is_label=None, label_characters='>'):
    r"""yields label and seq from a fasta file.


    Parameters
    ----------
    data : open file object or str
        An open fasta file or a path to it.

    strict : bool
        If strict is true a ``RecordError`` will
        be raised if no header line is found

    Returns
    -------
    label, sequence : string
        yields the label and sequence for each entry.

    Examples
    --------
    Assume we have a fasta formatted file with the following contents::

        >seq1
        CGATGTCGATCGATCGATCGATCAG
        >seq2
        CATCGATCGATCGATGCATGCATGCATG

    >>> from StringIO import StringIO
    >>> from skbio.parse.sequences import parse_fasta
    >>> fasta_f = StringIO('>seq1\n'
    ...                    'CGATGTCGATCGATCGATCGATCAG\n'
    ...                    '>seq2\n'
    ...                    'CATCGATCGATCGATGCATGCATGCATG\n')
    >>> for label, seq in parse_fasta(fasta_f):
    ...     print label
    ...     print seq
    seq1
    CGATGTCGATCGATCGATCGATCAG
    seq2
    CATCGATCGATCGATGCATGCATGCATG

    """

    for rec in finder(infile):
        # first line must be a label line
        if not rec[0][0] in label_characters:
            if strict:
                raise RecordError("Found Fasta record without label line: %s" %
                                  rec)
            else:
                continue
        # record must have at least one sequence
        if len(rec) < 2:
            if strict:
                raise RecordError("Found label line without sequences: %s" %
                                  rec)
            else:
                continue

        label = rec[0][1:].strip()
        label = label_to_name(label)
        seq = ''.join(rec[1:])

        yield label, seq


def parse_qual(infile, full_header=False):
    r"""yields label and qual from a qual file.


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
    ...     print label
    ...     print qual
    seq1
    [10 20 30 40]
    seq2
    [1 2 3 4]

    """
    for rec in FastaFinder(infile):
        curr_id = rec[0][1:]
        curr_qual = ' '.join(rec[1:])
        try:
            parts = np.asarray(curr_qual.split(), dtype=int)
        except ValueError:
            raise RecordError(
                "Invalid qual file. Check the format of the qual files.")
        if full_header:
            curr_pid = curr_id
        else:
            curr_pid = curr_id.split()[0]
        yield (curr_pid, parts)
