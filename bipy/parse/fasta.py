#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, The bipy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from bipy.parse.record_finder import LabeledRecordFinder
from bipy.parse.record import RecordError


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


def MinimalFastaParser(infile, strict=True,
                       label_to_name=str, finder=FastaFinder,
                       is_label=None, label_characters='>'):
    """Yields successive sequences from infile as (label, seq) tuples.

    If strict is True (default), raises RecordError when label or seq missing.
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


def MinimalGdeParser(infile, strict=True, label_to_name=str):
    return MinimalFastaParser(infile, strict, label_to_name, finder=GdeFinder,
                              label_characters='%#')


def xmfa_label_to_name(line):
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


def MinimalXmfaParser(infile, strict=True):
    # Fasta-like but with header info like ">1:10-1000 + chr1"
    return MinimalFastaParser(infile, strict, label_to_name=xmfa_label_to_name,
                              finder=XmfaFinder)
