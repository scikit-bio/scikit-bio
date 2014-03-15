#!/usr/bin/env python

"""
Sequences (:mod:`skbio.parse.sequences`)
========================================

.. currentmodule:: skbio.parse.sequences

This module provides a functions for parsing sequence files.

Functions
---------

.. autosummary::
   :toctree: generated/

    is_fasta_label
    is_gde_label
    is_blank_or_comment
    is_blank
    fasta_parse
    gde_parse
    xmfa_label_to_name
    xmfa_parse
    fastq_parse

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from skbio.core.exception import FastqParseError
from skbio.parse.record_finder import LabeledRecordFinder
from skbio.parse.record import RecordError


def is_fasta_label(x):
    """Checks if x looks like a FASTA label line."""
    return x.startswith('>')


def is_gde_label(x):
    """Checks if x looks like a GDE label line."""
    return x and x[0] in '%#'


def is_blank_or_comment(x):
    """Checks if x is blank or a FASTA comment line."""
    return (not x) or x.startswith('#') or x.isspace()


def is_blank(x):
    """Checks if x is blank."""
    return (not x) or x.isspace()

FastaFinder = LabeledRecordFinder(is_fasta_label, ignore=is_blank_or_comment)


def fasta_parse(infile,
                strict=True,
                label_to_name=str,
                finder=FastaFinder,
                is_label=None,
                label_characters='>'):
    """yields label and seq from a fasta file.


    Parameters
    ----------
    data : open file object
        An open fasta file.

    strict : bool
        If strict is true a RecordError error will
        be raised if no header line is found

    Yields
    ------
        yields the label and sequence for each entry.

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

GdeFinder = LabeledRecordFinder(is_gde_label, ignore=is_blank)


def gde_parse(infile, strict=True, label_to_name=str):
    """Parses a file with GDE label line"""
    return fasta_parse(infile,
                       strict,
                       label_to_name,
                       finder=GdeFinder,
                       label_characters='%#')


def xmfa_label_to_name(line):
    """returns name from xmfa label."""
    (loc, strand, contig) = line.split()
    (sp, loc) = loc.split(':')
    (lo, hi) = [int(x) for x in loc.split('-')]
    if strand == '-':
        (lo, hi) = (hi, lo)
    else:
        assert strand == '+'
    name = '%s:%s:%s-%s' % (sp, contig, lo, hi)
    return name


def is_xmfa_blank_or_comment(x):
    """Checks if x is blank or an XMFA comment line."""
    return (not x) or x.startswith('=') or x.isspace()

XmfaFinder = LabeledRecordFinder(is_fasta_label,
                                 ignore=is_xmfa_blank_or_comment)


def xmfa_parse(infile, strict=True):
    """Parses a file with xmfa label line"""
    # Fasta-like but with header info like ">1:10-1000 + chr1"
    return fasta_parse(infile,
                       strict,
                       label_to_name=xmfa_label_to_name,
                       finder=XmfaFinder)


def fastq_parse(data, strict=True):

    """yields label, seq, and qual from a fastq file.


    Parameters
    ----------
    data : open file object
        An open fastq file.

    strict : bool
        If strict is true a FastqParse error will be raised if the seq and qual
    labels dont' match.

    Yields
    ------
        yields the label, sequence and quality for each entry

    """
    # fastq format is very simple, defined by blocks of 4 lines
    line_num = -1
    record = []
    for line in data:
        line_num += 1
        if line_num == 4:
            if strict:  # make sure the seq and qual labels match
                if record[0][1:] != record[2][1:]:
                    raise FastqParseError('Invalid format: %s -- %s'
                                          % (record[0][1:], record[2][1:]))
            yield record[0][1:], record[1], record[3]
            line_num = 0
            record = []
        record.append(line.strip())

    if record:
        if strict and record[0]:  # make sure the seq and qual labels match
            if record[0][1:] != record[2][1:]:
                raise FastqParseError('Invalid format: %s -- %s'
                                      % (record[0][1:], record[2][1:]))

        if record[0]:  # could be just an empty line at eof
            yield record[0][1:], record[1], record[3]

    if type(data) == file:
        data.close()
