#!/usr/bin/env python
r"""
Parse stockholm files (:mod:`skbio.parse.stockholm`)
=========================================================

.. currentmodule:: skbio.parse.stockholm

This module provides functions for parsing stockholm files.

Functions
---------

.. autosummary::
   :toctree: generated/

    parse_stockholm


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from skbio.core.exception import StockholmParseError
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
                raise StockholmParseError("GC must have exactly one char "
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

    strict : bool (optional)
        Turns on strict parsing of GR and GC lines to ensure one char per pos
        Default: False

    Returns
    -------
    sto : named tuple
        yields in Stockholm named tuple format. For more information, see
        skbio.format.stockholm.Stockholm

    Raises
    ------
    StockholmParseError
        If any lines are found that don't conform to stockholm format

    Examples
    --------
    Assume we have a very basic stockholm file with the following contents::

        # STOCKHOLM 1.0
        seq1         ACC--G-GGGU
        seq2         TCC--G-GGGA
        #=GC SS_cons (((.....)))
        //

    >>> from skbio.core.sequence import RNA
    >>> from skbio.parse.stockholm import parse_stockholm
    >>> from StringIO import StringIO
    >>> sto_in = StringIO("# STOCKHOLM 1.0\n"
    ...                  "seq1         ACC--G-GGGU\nseq2         TCC--G-GGGA\n"
    ...                  "#=GC SS_cons (((.....)))\n//")
    >>> sto_records = parse_stockholm(sto_in, RNA)
    >>> for sto in sto_records:
    >>>     print sto
    Stockholm(aln=<Alignment: n=2; mean +/- std length=11.00 +/- 0.00>, GF={},
              GS={}, GR={}, GC={'SS_cons': '(((.....)))'})

    Notes
    -----
    If a single record stockholm file is being parsed, you can add .next()
    to get the stockholm record directly without the generator obj:

    sto = parse_stockholm(sto_in, RNA).next()
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
