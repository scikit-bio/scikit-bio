#!/usr/bin/env python
r"""
Parse biological sequences (:mod:`skbio.parse.sequences`)
=========================================================

.. currentmodule:: skbio.parse.sequences

This module provides functions for parsing sequence files.

Functions
---------

.. autosummary::
   :toctree: generated/

    parse_fasta
    parse_fastq


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from skbio.core.exception import FastqParseError, RecordError
from skbio.parse.record_finder import (LabeledRecordFinder,
                                       DelimitedRecordFinder)
from skbio.core.sequence import BiologicalSequence
from skbio.core.alignment import (Alignment, SequenceCollection,
                                  SequenceCollectionError)
from skbio.parse.clustal import ClustalParser
from skbio.core.sequence import DNA

def is_fasta_label(x):
    """Checks if x looks like a FASTA label line."""
    return x.startswith('>')


def is_blank_or_comment(x):
    """Checks if x is blank or a FASTA comment line."""
    return (not x) or x.startswith('#') or x.isspace()


def is_blank(x):
    """Checks if x is blank."""
    return (not x) or x.isspace()

def is_rfam_header_line(line):
    """Returns True if line is a header line"""
    return line.startswith('#=GF')

def is_rfam_seq_line(line):
    """Returns True if line is a sequence line"""
    return bool(line) and (not line[0].isspace()) and \
    (not line.startswith('#')) and (not line.startswith('//'))

def is_rfam_structure_line(line):
    """Returns True if line is a structure line"""
    return line.startswith('#=GC SS_cons')

def load_from_clustal(data, seq_constructor=BiologicalSequence, strict=True):
    recs = [(name, seq_constructor(seq)) for name, seq in\
        ClustalParser(data, strict)]
    lengths = [len(i[1]) for i in recs]

    if lengths and max(lengths) == min(lengths):
        return Alignment.from_fasta_records(recs, DNA)
    else:
        return SequenceCollection.from_fasta_records(recs, DNA)

def is_empty_or_html(line):
    """Return True for HTML line and empty (or whitespace only) line.

    line -- string

    The Rfam adaptor that retrieves records inlcudes two HTML tags in
    the record. These lines need to be ignored in addition to empty lines. 
    """
    if line.startswith('<pre') or line.startswith('</pre'):
        return True
    return (not line) or line.isspace()


FastaFinder = LabeledRecordFinder(is_fasta_label, ignore=is_blank_or_comment)
RfamFinder = DelimitedRecordFinder('//', ignore=is_empty_or_html)


def parse_fasta(infile,
                strict=True,
                label_to_name=str,
                finder=FastaFinder,
                is_label=None,
                label_characters='>'):
    r"""yields label and seq from a fasta file.


    Parameters
    ----------
    data : open file object
        An open fasta file.

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


def parse_fastq(data, strict=False):
    r"""yields label, seq, and qual from a fastq file.

    Parameters
    ----------
    data : open file object
        An open fastq file.

    strict : bool
        If strict is true a FastqParse error will be raised if the seq and qual
        labels dont' match.

    Returns
    -------
    label, seq, qual : string
        yields the label, sequence and quality for each entry

    Examples
    --------
    Assume we have a fastq formatted file with the following contents::

        @seq1
        AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
        +
        ````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
        @seq2
        TATGTATATATAACATATACATATATACATACATA
        +
        ]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb

    We can use the following code:

    >>> from StringIO import StringIO
    >>> from skbio.parse.sequences import parse_fastq
    >>> fastq_f = StringIO('@seq1\n'
    ...                     'AACACCAAACTTCTCCACCACGTGAGCTACAAAAG\n'
    ...                     '+\n'
    ...                     '````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF\n'
    ...                     '@seq2\n'
    ...                     'TATGTATATATAACATATACATATATACATACATA\n'
    ...                     '+\n'
    ...                     ']KZ[PY]_[YY^```ac^\\\`bT``c`\\aT``bbb\n')
    >>> for label, seq, qual in parse_fastq(fastq_f):
    ...     print label
    ...     print seq
    ...     print qual
    seq1
    AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
    ````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
    seq2
    TATGTATATATAACATATACATATATACATACATA
    ]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb

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

def ChangedSequence(data, seq_constructor=BiologicalSequence):
    """Returns new BiologicalSequence object, replaces dots with dashes in sequence.
    """
    return str(seq_constructor(str(data).replace('.','-')))

def MinimalRfamParser(infile, strict=True, seq_constructor=ChangedSequence):
    """Yield successive sequences as (header, sequences, structure) tuples.
    
    header is a list of header lines
    sequences is an Alignment object. Sequences are objects keyed by the
        original labels in the database.
    structure is a WussStructure
    """
    for record in RfamFinder(infile):
        header = []
        sequences = []
        structure = []
        for line in record:
            if is_rfam_header_line(line):
                header.append(line.strip())
            elif is_rfam_seq_line(line):
                sequences.append(line)
            elif is_rfam_structure_line(line):
                structure.append(line)
            else:
                continue
        #sequence and structure are required. 
        #for example when looking at the stockholm format of just one family
        if not sequences or not structure:
            if strict:
                error = 'Found record with missing element(s): '
                if not sequences:
                    error += 'sequences '
                if not structure:
                    error += 'structure '
                raise RecordError, error
            else:
                continue
        #join all sequence parts together, construct label
        try:
            new_seqs = load_from_clustal(sequences, strict=strict,
                                         seq_constructor=seq_constructor)
            sequences = new_seqs
        except SequenceCollectionError as e:
            if strict:
                raise RecordError, str(e)
            else:
                continue

        #construct the structure
        try:
            res = load_from_clustal(structure, strict=strict)
            assert len(res.NamedSeqs) == 1 #otherwise multiple keys
            structure = res.NamedSeqs['#=GC SS_cons']
        except SequenceCollectionError as e:
            if strict:
                raise RecordError("Can't parse structure of family: %s" %
                                  str(header))
            else:
                structure = None
        yield header, sequences, structure

