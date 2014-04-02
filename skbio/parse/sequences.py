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
    parse_stockholm


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from skbio.core.exception import (FastqParseError, RecordError,
                                  StockholmParseError)
from skbio.parse.record_finder import LabeledRecordFinder
from skbio.format.stockholm import Stockholm
from skbio.core.alignment import Alignment


def _parse_gf_info(lines):
    """Takes care of parsing GF lines in stockholm including special cases"""
    parsed = {}
    #needed for making each multi-line RT and NH one string
    rt = []
    nh = []
    lastline = ""
    for line in lines:
        try:
            init, feature, content = line.split(None, 2)
        except ValueError:
            raise StockholmParseError("Malformed GF line encountered!\n%s" %
                                      line.split(None, 2))
        if init != "#=GF":
            raise StockholmParseError("Non-GF line encountered!")

        #take care of adding multiline RT to the parsed information
        if lastline == "RT" and feature != "RT":
            #add rt line to the parsed dictionary
            rtline = " ".join(rt)
            rt = []
            if "RT" in parsed:
                parsed["RT"].append(rtline)
            else:
                parsed["RT"] = [rtline]
        elif feature == "RT":
            rt.append(content)
            lastline = feature
            continue
        #Take care of adding multiline NH to the parsed dictionary
        elif lastline == "NH" and feature != "NH":
            nhline = " ".join(nh)
            nh = []
            if "NH" in parsed:
                parsed["NH"].append(nhline)
            else:
                parsed["NH"] = [nhline]
        elif feature == "NH":
            nh.append(content)
            lastline = feature
            continue

        #add current feature to the parsed information
        if feature in parsed:
            parsed[feature].append(content)
        else:
            parsed[feature] = [content]
        lastline = feature

    #clean up parsed info by removing unneccessary lists
    for feature in parsed:
        #list of multi-line features to join into single string if necessary
        if feature in ["CC"]:
            parsed[feature] = ' '.join(parsed[feature])
        elif len(parsed[feature]) == 1:
            parsed[feature] = parsed[feature][0]
    return parsed


def _parse_gc_info(lines, strict=False, seqlen=-1):
    """Takes care of parsing GC lines in stockholm format"""
    parsed = {}
    for line in lines:
        try:
            init, feature, content = line.split(None, 2)
        except ValueError:
            raise StockholmParseError("Malformed GC line encountered!\n%s" %
                                      line.split(None, 2))
        if init != "#=GC":
            raise StockholmParseError("Non-GC line encountered!")

        #add current feature to the parsed information
        if feature in parsed:
            if strict:
                raise StockholmParseError("Should not have multiple lines "
                                          "with same feature: %s" % feature)
            parsed[feature].append(content)
        else:
            parsed[feature] = [content]

    #clean up parsed info by removing unneccessary lists
    for feature in parsed:
        parsed[feature] = ''.join(parsed[feature])
        if strict:
            if len(parsed[feature]) != seqlen:
                raise StockholmParseError("GR must have exactly one char "
                                          "per position in alignment!")

    return parsed


def _parse_gs_gr_info(lines, strict=False, seqlen=-1):
    """Takes care of parsing GS and GR lines in stockholm format"""
    parsed = {}
    for line in lines:
        try:
            init, label, feature, content = line.split(None, 3)
        except ValueError:
            raise StockholmParseError("Malformed GS/GR line encountered!\n%s" %
                                      line.split(None, 3))
        if init != "#=GS" and init != "#=GR":
                raise StockholmParseError("Non-GS/GR line encountered!")

        #parse each line, taking into account we can have interleaved format
        if label in parsed and feature in parsed[label]:
                #interleaved format, so need list of content
                parsed[label][feature].append(content)
        else:
            parsed[label] = {feature: [content]}

    #join all the crazy lists created during parsing
    for label in parsed:
        for feature, content in parsed[label].items():
            parsed[label][feature] = ''.join(content)
            if strict:
                if len(parsed[label][feature]) != seqlen:
                    raise StockholmParseError("GR must have exactly one char "
                                              "per position in alignment!")
    return parsed

def parse_stockholm(infile, seq_constructor, strict=False):
    r"""yields Stockholm objects from a stockholm file.

    Parameters
    ----------
    infile : open file object
        An open stockholm file.

    seq_constructor : BiologicalSequence object
        The biologicalsequence object that corresponds to what the stockholm
        file holds. See skbio.core.sequence

    Returns
    -------
    sto : named tuple
        yields in Stockholm named tuple format. For more information, see
        skbio.format.stockholm.Stockholm
    """
    #make sure first line is corect
    line = infile.readline()
    if not line.startswith("# STOCKHOLM 1.0"):
        raise StockholmParseError("Incorrect header found")
    gs_lines = []
    gf_lines = []
    gr_lines = []
    gc_lines = []
    seqs = {}
    for line in infile:
        line = line.strip()
        if line == "":
            #skip blank lines
            continue
        elif line == "//":
            #create alignment from stockholm file
            aln = Alignment.from_fasta_records(seqs.items(), seq_constructor)
            seqlen = len(aln)
            #parse information lines
            GF = _parse_gf_info(gf_lines)
            GS = _parse_gs_gr_info(gs_lines)
            GR = _parse_gs_gr_info(gr_lines, strict, seqlen)
            GC = _parse_gc_info(gc_lines, strict, seqlen)
            #yield the actual stockholm object
            yield Stockholm(aln, GF, GS, GR, GC)
        elif line.startswith("#=GF"):
            gf_lines.append(line)
        elif line.startswith("#=GS"):
            gs_lines.append(line)
        elif line.startswith("#=GR"):
            gr_lines.append(line)
        elif line.startswith("#=GC"):
            gc_lines.append(line)
        else:
            lineinfo = line.split()
            #assume sequence since nothing else in format is left
            #in case of interleaved format, need to do check
            if lineinfo[0] in seqs:
                seqs[lineinfo[0]] += lineinfo[1]
            else:
                seqs[lineinfo[0]] = lineinfo[1]


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
