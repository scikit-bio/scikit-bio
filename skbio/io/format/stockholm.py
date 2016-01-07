"""
Stockholm format (:mod:`skbio.io.format.stockholm`)
===================================================

.. currentmodule:: skbio.io.format.stockholm

The Stockholm format is a multiple sequence alignment format written in Markup
as opposed to a simple text-based format. Data is stored on seperate lines,
each with a unique data 'name' preceding it. The alignment can just contain
data, or contain data and related metadata.

An example Stockholm file, taken from [1]_:

.. code-block:: none

    # STOCKHOLM 1.0
    #=GF ID    UPSK
    #=GF SE    Predicted; Infernal
    #=GF SS    Published; PMID 9223489
    #=GF RN    [1]
    #=GF RM    9223489
    #=GF RT    The role of the pseudoknot at the 3' end of turnip yellow mosaic
    #=GF RT    virus RNA in minus-strand synthesis by the viral RNA-dependent \
RNA
    #=GF RT    polymerase.
    #=GF RA    Deiman BA, Kortlever RM, Pleij CW;
    #=GF RL    J Virol 1997;71:5990-5996.
    AF035635.1/619-641             UGAGUUCUCGAUCUCUAAAAUCG
    M24804.1/82-104                UGAGUUCUCUAUCUCUAAAAUCG
    J04373.1/6212-6234             UAAGUUCUCGAUCUUUAAAAUCG
    M24803.1/1-23                  UAAGUUCUCGAUCUCUAAAAUCG
    #=GC SS_cons                   .AAA....<<<<aaa....>>>>
    //

Format Support
==============
**Has Sniffer: Yes**

**State: Experimental as of 0.4.0-dev.**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`skbio.alignment.TabularMSA`                              |
+------+------+---------------------------------------------------------------+

Format Specification
====================
Stockholm format contains two types of data: raw DNA, RNA, or Protein data and
associated metadata. Raw data lines begin with an associated 'name', which
often times is an id comprised of letters and numbers, such as
'M24804.1/82-104' ([1]_). After the id and a few tab characters,
the data is displayed. Metadata lines, however, begin with a '#' and a
two-letter marker which describes the metadata type, for example '#=GF' (see
Metadata Types). Each metadata line also contains an extra two-letter tag,
such as 'AS' or 'CC' which tells what type of data it precedes. All metadata
is optional.

Metadata Types
++++++++++++++
GF
--
Data relating to multiple sequence alignment as a whole, such as authors or
number of sequences in the alignment. Starts with #=GF followed by a tag and
data relating to the tag. Typically comes first in a Stockholm file.
For example:

.. code-block:: none

    #=GF DE CBS domain

Example taken from [2]_.

GF Types
^^^^^^^^
+---+------------------------+------------------------------------------------+
|Tag|Name                    |Description                                     |
+===+====+===================+================================================+
|AC |Accession Number        |Accession number in form PFxxxxx (Pfam) or RFxx\|
|   |                        |xxx (Rfam)                                      |
+---+------------------------+------------------------------------------------+
|ID |Identification          |One word name for family                        |
+---+------------------------+------------------------------------------------+
|DE |Definition              |Short description of family                     |
+---+------------------------+------------------------------------------------+
|AU |Author                  |Authors of the entry                            |
+---+------------------------+------------------------------------------------+
|SE |Source of seed          |The source suggesting the seed members belong t\|
|   |                        |o one family                                    |
+---+------------------------+------------------------------------------------+
|SS |Source of structure     |The source (prediction or publication) of the c\|
|   |                        |onsensus RNA secondary structure used by Rfam   |
+---+------------------------+------------------------------------------------+
|BM |Build method            |Command line used to generate the model         |
+---+------------------------+------------------------------------------------+
|SM |Search method           |Command line used to perform the search         |
+---+------------------------+------------------------------------------------+
|GA |Gathering threshold     |Search threshold to build the full alignment    |
+---+------------------------+------------------------------------------------+
|TC |Trusted cutoff          |Lowest sequence score (and domain score for Pfa\|
|   |                        |m) of match in the full alignment               |
+---+------------------------+------------------------------------------------+
|NC |Noise cutoff            |Highest sequence score (and domain score for Pf\|
|   |                        |am) of match not in full alignment              |
+---+------------------------+------------------------------------------------+
|TP |Type                    |Type of family                                  |
+---+------------------------+------------------------------------------------+
|SQ |Sequence                |Number of sequences in alignment                |
+---+------------------------+------------------------------------------------+
|DC |Database comment        |Comment about database reference                |
+---+------------------------+------------------------------------------------+
|DR |Database reference      |Reference to external database                  |
+---+------------------------+------------------------------------------------+
|RC |Reference comment       |Comment about literature reference              |
+---+------------------------+------------------------------------------------+
|RN |Reference number        |Reference number                                |
+---+------------------------+------------------------------------------------+
|RM |Reference medline       |Eight digit medline UI number                   |
+---+------------------------+------------------------------------------------+
|RT |Reference title         |Reference title                                 |
+---+------------------------+------------------------------------------------+
|RA |Reference author        |Reference author                                |
+---+------------------------+------------------------------------------------+
|RL |Reference location      |Journal location                                |
+---+------------------------+------------------------------------------------+
|PI |Previous identifier     |Record of all previous ID lines                 |
+---+------------------------+------------------------------------------------+
|KW |Keywords                |Keywords                                        |
+---+------------------------+------------------------------------------------+
|CC |Comment                 |Comment                                         |
+---+------------------------+------------------------------------------------+
|NE |Pfam accession          |Indicates a nested domain                       |
+---+------------------------+------------------------------------------------+
|NL |Location                |Location of nested domains - sequence ID, start |
|   |                        |and end of insert                               |
+---+------------------------+------------------------------------------------+

Descriptions taken from [1]_.

GS
--
Data relating to a specific sequence in the multiple sequence alignment.
Starts with #=GS followed by the sequence name followed by a tag and data
relating to the tag. Typically comes second in a Stockholm file.
For example:

.. code-block:: none

    #=GS O83071/259-312 AC O83071

Example taken from [2]_.

GS Types
^^^^^^^^
+---+-------------------------------------------+
|Tag|Description                                |
+===+===========================================+
|AC |Accession number                           |
+---+-------------------------------------------+
|DE |Description                                |
+---+-------------------------------------------+
|DR |Database reference                         |
+---+-------------------------------------------+
|OS |Organism (species)                         |
+---+-------------------------------------------+
|OC |Organism classification (clade, etc.)      |
+---+-------------------------------------------+
|LO |Look (color, etc.)                         |
+---+-------------------------------------------+

Descriptions taken from [1]_.

GR
--
Data relating to the columns of a specific sequence in a multiple sequence
alignment. Starts with #=GR followed by the sequence name followed by a tag
and data relating to the tag. Typically comes after the data line it relates
to.
For example:

.. code-block:: none

    #=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEHHH

Example taken from [2]_.

GR Types
^^^^^^^^
+---+-------------------------------------------+
|Tag|Description                                |
+===+===========================================+
|SS |Secondary structure                        |
+---+-------------------------------------------+
|SA |Surface accessability                      |
+---+-------------------------------------------+
|TM |Transmembrane                              |
+---+-------------------------------------------+
|PP |Posterior probability                      |
+---+-------------------------------------------+
|LI |Ligand binding                             |
+---+-------------------------------------------+
|AS |Active site                                |
+---+-------------------------------------------+
|pAS|AS - Pfam predicted                        |
+---+-------------------------------------------+
|sAS|AS - from SwissProt                        |
+---+-------------------------------------------+
|IN |Intron (in or after)                       |
+---+-------------------------------------------+

Descriptions taken from [2]_.

GC
--
Data relating to the columns of the multiple sequence alignment as a whole.
Starts with #=GC followed by a tag and data relating to the tag. Typically
comes at the end of the multiple sequence alignment.
For example:

.. code-block:: none

    #=GC SS_cons            CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEH

Example taken from [2]_.

GC Types
^^^^^^^^
Same as GR Types, but add "_cons" to the end of tag.

Examples
========
Suppose we have a Stockholm file:

>>> import skbio.io
>>> from io import StringIO
>>> from skbio.sequence import RNA
>>> from skbio.alignment import TabularMSA
>>> fs = '\\n'.join([
...         '# STOCKHOLM 1.0',
...         '#=GF RA    Deiman BA, Kortlever RM, Pleij CW;',
...         '#=GF RL    J Virol 1997;71:5990-5996.',
...         'AF035635.1/619-641             UGAGUUCUCGAUCUCUAAAAUCG',
...         'M24804.1/82-104                UGAGUUCUCUAUCUCUAAAAUCG',
...         'J04373.1/6212-6234             UAAGUUCUCGAUCUUUAAAAUCG',
...         'M24803.1/1-23                  UAAGUUCUCGAUCUCUAAAAUCG',
...         '#=GC SS_cons                   .AAA....<<<<aaa....>>>>',
...         '//'
...     ])
>>> fh = StringIO(fs)
>>> msa = skbio.io.read(fh, into=TabularMSA, constructor=RNA)
>>> msa # doctest: +NORMALIZE_WHITESPACE
TabularMSA[RNA]
-------------------------------------------------
Metadata:
'RA': '   Deiman BA, Kortlever RM, Pleij CW;'
'RL': '   J Virol 1997;71:5990-5996.'
Positional metadata:
'SS_cons': <dtype: object>
Stats:
sequence count: 4
position count: 23
-------------------------------------------------
UGAGUUCUCGAUCUCUAAAAUCG
UGAGUUCUCUAUCUCUAAAAUCG
UAAGUUCUCGAUCUUUAAAAUCG
UAAGUUCUCGAUCUCUAAAAUCG

References
==========
.. [1] https://en.wikipedia.org/wiki/Stockholm_format
.. [2] http://sonnhammer.sbc.su.se/Stockholm.html
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from collections import OrderedDict

from skbio.alignment import TabularMSA
from skbio.sequence import Protein
from skbio.io import create_format, StockholmFormatError

stockholm = create_format('stockholm')


@stockholm.sniffer()
def _stockholm_sniffer(fh):
    # Smells a Stockholm file if the first line contains 'Stockholm' and the
    # version number.
    line = next(fh)
    if line[:15] == "# STOCKHOLM 1.0":
        return True, {}

    return False, {}


@stockholm.reader(TabularMSA)
def _stockholm_to_tabular_msa(fh, constructor=Protein):
    # Uses OrderedDict() to make sure dna_data isn't arranged randomly
    dna_data = OrderedDict()
    metadata = {}
    positional_metadata = {}
    seqs = []

    # Retrieves data from file
    for line in fh:
        if not (line.startswith("#") or line.startswith("//") or
                line.isspace() == True):
            dna_data = _parse_stockholm_line_data(line, dna_data)
    fh.seek(0)

    # Retrieves metadata from file
    for line in fh:
        if line.startswith("#=GF"):
            metadata = _parse_stockholm_line_gf(line, metadata)
        if line.startswith("#=GS"):
            dna_data = _parse_stockholm_line_gs(line, dna_data)
        if line.startswith("#=GR"):
            dna_data = _parse_stockholm_line_gr(line, dna_data)
        if line.startswith('#=GC'):
            positional_metadata = _parse_stockholm_line_gc(line,
                                                           positional_metadata)

    for key in dna_data.keys():
        # Sets blank dictionaries and lists to None instead
        if not dna_data[key][1]:
            dna_data[key][1] = None
        if not dna_data[key][2]:
            dna_data[key][2] = None
        # Adds each sequence to the MSA data
        seqs.append(constructor(dna_data[key][0], metadata=dna_data[key][1],
                                positional_metadata=dna_data[key][2]))

    if not seqs:
        raise StockholmFormatError("No data present in file.")
    if not positional_metadata:
        positional_metadata = None
    # Constructs TabularMSA
    return TabularMSA(seqs, metadata=metadata,
                      positional_metadata=positional_metadata)


def _parse_stockholm_line_gf(line, metadata):
    """Takes ``#=GF`` line and returns parsed data."""
    line = _remove_newline(line.split(' ', 2))
    if line[1] in metadata.keys():
        metadata[line[1]] = metadata[line[1]] + ' ' + line[2]
    else:
        metadata[line[1]] = line[2]
    return metadata


def _parse_stockholm_line_gs(line, dna_data):
    """Takes ``#=GS`` line and returns parsed data."""
    line = _remove_newline(line.split(' ', 3))
    if line[1] in dna_data.keys():
        dna_data[line[1]][1][line[2]] = line[3]
    else:
        raise StockholmFormatError("Markup line references nonexistent "
                                   "data %r." % line[1])
    return dna_data


def _parse_stockholm_line_gr(line, dna_data):
    """Takes ``#=GR`` line and returns parsed data."""
    line = _remove_newline(line.split(' '))
    if len(line) != 4:
        del(line[3:len(line)-1])
    if line[1] in dna_data.keys():
        if line[2] in dna_data[line[1]][2].keys():
            raise StockholmFormatError("Found duplicate GR label %r associated"
                                       " with data label %r" % (line[2],
                                                                line[1]))
        dna_data[line[1]][2][line[2]] = list(line[3])
    else:
        raise StockholmFormatError("Markup line references nonexistent "
                                   "data %r." % line[1])
    return dna_data


def _parse_stockholm_line_gc(line, positional_metadata):
    """Takes ``#=GC`` line and returns parsed data."""
    line = _remove_newline(line.split(' '))
    if len(line) > 3:
        del(line[2:len(line)-1])
    if line[1] in positional_metadata.keys():
        raise StockholmFormatError("Found duplicate GC label %r." % (line[1]))
    positional_metadata[line[1]] = list(line[2])
    return positional_metadata


def _parse_stockholm_line_data(line, dna_data):
    """Takes data line and returns parsed data."""
    line = line.split()
    if line[0] not in dna_data.keys():
        dna_data[line[0]] = [line[1], {}, {}]
    elif line[0] in dna_data.keys():
        raise StockholmFormatError("Found multiple data lines under same "
                                   "name: %r" % line[0])
    return dna_data


def _remove_newline(line):
    """Removes '\n' from line and returns line."""
    if '\n' in line[len(line)-1]:
        line[len(line)-1] = line[len(line)-1].rstrip('\n')
    return line
