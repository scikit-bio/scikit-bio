#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""Provides some classes for treating files as sequences of records.

Typically more useful as subclasses. Covers the three main types of records:

    DelimitedRecordFinder:  Records demarcated by an end line, e.g. '\\'
    LabeledRecordFinder:    Records demarcated by a start line, e.g. '>label'
    LineGrouper:            Records consisting of a certain number of lines.
    TailedRecordFinder:     Records demarcated by an end mark, e.g. 'blah.'

All the first classes ignore/delete blank lines and strip leading and trailing
whitespace.  The TailedRecodeFinder is Functional similar to
DelimitedRecordFinder except that it accept a is_tail function instead of a
str.  Note that its default constuctor is rstrip instead of strip.
"""

from string import strip, rstrip

from skbio.parse.record import RecordError


def is_empty(line):
    """Returns True empty lines and lines consisting only of whitespace."""
    return (not line) or line.isspace()


def never_ignore(line):
    """Always returns False."""
    return False


def DelimitedRecordFinder(delimiter, constructor=strip, ignore=is_empty,
                          keep_delimiter=True, strict=True):
    """Returns function that returns successive delimited records from file.

    Includes delimiter in return value. Returns list of relevant lines.

    Default constructor is string.strip, but can supply another constructor
    to transform lines and/or coerce into correct type. If constructor is None,
    passes along the lines without alteration.

    Skips any lines for which ignore(line) evaluates True (default is to skip
    whitespace).

    keep_delimiter: keep delimiter line at the end of last block if True
    (default), otherwise discard delimiter line.

    strict: when lines found after the last delimiter -- raise error if True
    (default), otherwise yield the lines silently
    """
    def parser(lines):
        curr = []
        for line in lines:
            if constructor:
                line = constructor(line)
            # else:
            #    line = l
            # ignore blank lines
            if ignore(line):
                continue
            # if we find the delimiter, return the line; otherwise, keep it
            if line == delimiter:
                if keep_delimiter:
                    curr.append(line)
                yield curr
                curr = []
            else:
                curr.append(line)
        if curr:
            if strict:
                raise RecordError("Found additional data after records: %s" %
                                  (curr))
            else:
                yield curr
    return parser

# The following is an example of the sorts of iterators RecordFinder returns.
GbFinder = DelimitedRecordFinder('//')


def TailedRecordFinder(is_tail_line, constructor=rstrip, ignore=is_empty,
                       strict=True):
    """Returns function that returns successive tailed records from lines.

    Includes tail line in return value. Returns list of relevant lines.

    constructor: a modifier for each line, default is string.rstrip: to remove
    \n and trailing spaces.

    Skips over any lines for which ignore(line) evaluates True (default is
    to skip empty lines).  note that the line maybe modified by constructor.

    strict: if True(default), raise error if the last line is not a tail.
    otherwise, yield the last lines.
    """
    def parser(lines):
        curr = []
        for line in lines:
            if constructor:
                line = constructor(line)
            if ignore(line):
                continue

            curr.append(line)
            # if we find the label, return the previous record
            if is_tail_line(line):
                yield curr
                curr = []

        # don't forget to return the last record in the file
        if curr:
            if strict:
                raise RecordError('lines exist after the last tail_line '
                                  'or no tail_line at all')
            else:
                yield curr

    return parser


def LabeledRecordFinder(is_label_line, constructor=strip, ignore=is_empty):
    """Returns function that returns successive labeled records from file.

    Includes label line in return value. Returns list of relevant lines.

    Default constructor is string.strip, but can supply another constructor
    to transform lines and/or coerce into correct type. If constructor is None,
    passes along the lines without alteration.

    Skips over any lines for which ignore(line) evaluates True (default is
    to skip empty lines).

    NOTE: Does _not_ raise an exception if the last line is a label line: for
    some formats, this is acceptable. It is the responsibility of whatever is
    parsing the sets of lines returned into records to complain if a record
    is incomplete.
    """
    def parser(lines):
        curr = []
        for l in lines:
            if constructor:
                line = constructor(l)
            else:
                line = l
            if ignore(line):
                continue
            # if we find the label, return the previous record
            if is_label_line(line):
                if curr:
                    yield curr
                    curr = []
            curr.append(line)
        # don't forget to return the last record in the file
        if curr:
            yield curr
    return parser


def is_fasta_label(x):
    """Checks if x looks like a FASTA label line."""
    return x.startswith('>')
# The following is an example of the sorts of iterators RecordFinder returns.
FastaFinder = LabeledRecordFinder(is_fasta_label)


def LineGrouper(num, constructor=strip, ignore=is_empty):
    """Returns num lines at a time, stripping and ignoring blanks.

    Default constructor is string.strip, but can supply another constructor
    to transform lines and/or coerce into correct type. If constructor is None,
    passes along the lines without alteration.

    Skips over any lines for which ignore(line) evaluates True: default is to
    skip whitespace lines.

    """
    def parser(lines):
        curr = []
        for l in lines:
            if constructor:
                line = constructor(l)
            else:
                line = l
            if ignore(line):
                continue
            curr.append(line)
            if len(curr) == num:
                yield curr
                curr = []
        if curr:
            raise RecordError("Non-blank lines not even multiple of %s" % num)
    return parser
