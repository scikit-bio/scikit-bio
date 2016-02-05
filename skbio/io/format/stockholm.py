"""
Stockholm format (:mod:`skbio.io.format.stockholm`)
===================================================

.. currentmodule:: skbio.io.format.stockholm

The Stockholm format is a multiple sequence alignment format written in Markup
as opposed to a simple text-based format. Data is stored on separate lines,
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

**State: Experimental as of 0.4.1-dev.**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`skbio.alignment.TabularMSA`                              |
+------+------+---------------------------------------------------------------+

Format Specification
====================
The Stockholm format contains two types of data. The first can contain raw DNA,
RNA, or Protein data and the second is comprised of associated metadata.
Raw data lines begin with an associated 'name', which
often times is an id comprised of letters and numbers, such as
'M24804.1/82-104' ([1]_). After the id and a few tab characters,
the data is displayed. Metadata lines, however, begin with a '#' and a
two-letter marker which describes the metadata type, for example '#=GF' (see
Metadata Types). Each metadata line also contains an extra two-letter feature,
such as 'AS' or 'CC' which tells what type of data it precedes. All metadata
is optional.

Metadata Types
++++++++++++++
GF
--
Data relating to the multiple sequence alignment as a whole, such as authors or
number of sequences in the alignment. Starts with #=GF followed by a feature
and data relating to the feature. Typically comes first in a Stockholm file.
For example:

.. code-block:: none

    #=GF DE CBS domain

Example taken from [2]_.

GS
--
Data relating to a specific sequence in the multiple sequence alignment.
Starts with #=GS followed by the sequence name followed by a feature and data
relating to the feature. Typically comes second in a Stockholm file.
For example:

.. code-block:: none

    #=GS O83071/259-312 AC O83071

Example taken from [2]_.

GC
--
Data relating to the columns of the multiple sequence alignment as a whole.
Starts with #=GC followed by a feature and data relating to the feature.
Typically comes at the end of the multiple sequence alignment.
For example:

.. code-block:: none

    #=GC SS_cons            CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEH

Example taken from [2]_.

GR
--
Data relating to the columns of a specific sequence in a multiple sequence
alignment. Starts with #=GR followed by the sequence name followed by a feature
and data relating to the feature. Typically comes after the data line it
relates to.
For example:

.. code-block:: none

    #=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEHHH

Example taken from [2]_.

Examples
========
Suppose we have a Stockholm file:

>>> import skbio.io
>>> from io import StringIO
>>> from skbio import RNA, TabularMSA
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

from collections import OrderedDict, namedtuple

from skbio.alignment import TabularMSA
from skbio.sequence._iupac_sequence import IUPACSequence
from skbio.io import create_format, StockholmFormatError

stockholm = create_format('stockholm')
# Creates _SeqData namedtuple to increase readability and efficiency
_SeqData = namedtuple("SeqData", ["seq", "metadata", "positional_metadata"])
supported_data_types = ['#=GF', '#=GS', '#=GR', '#=GC']


@stockholm.sniffer()
def _stockholm_sniffer(fh):
    # Smells a Stockholm file if the following conditions are met:
    # - File isn't empty
    # - File contains correct header
    # - FIle contains correct footer
    try:
        line = next(fh)
    except StopIteration:
        return False, {}

    header = line
    for line in fh:
        pass
    if line.startswith('//') and header == "# STOCKHOLM 1.0\n":
        return True, {}

    return False, {}


@stockholm.reader(TabularMSA)
def _stockholm_to_tabular_msa(fh, constructor=None):
    # Checks that user has passed required constructor parameter
    if constructor is None:
        raise ValueError("Must provide `constructor` parameter.")
    # Checks that contructor parameter is supported
    elif not issubclass(constructor, IUPACSequence):
        raise TypeError("`constructor` must be a subclass of `IUPACSequence`.")

    # Checks that the file isn't empty
    try:
        line = next(fh)
    except StopIteration:
        raise StockholmFormatError("File is empty.")
    # Checks that the file follows basic format (includes the required header)
    if not line == "# STOCKHOLM 1.0\n":
        raise StockholmFormatError("File missing required header: "
                                   "`# STOCKHOLM 1.0\\n`.")

    # Uses OrderedDict() to make sure dna_data isn't arranged randomly
    dna_data = OrderedDict()
    metadata = {}
    positional_metadata = {}
    seqs = []

    # Retrieves data from file, reads first so that data order will be kept
    # consistent.
    # print(len(fh.readlines()[-1]))
    for line in fh:
        if is_data_line(line):
            dna_data = _parse_stockholm_line_data(line, dna_data)
        if not line.isspace():
            data_type = line.split()[0]
        if (data_type.startswith('#=') and
            data_type not in supported_data_types):
            raise StockholmFormatError("Unrecognized data "
                                       "type %r." % data_type)
        if line.startswith('//'):
            break
    if not line.startswith('//'):
        raise StockholmFormatError('Final line does not conform to Stockholm '
                                   'format. Must contain only `//`.')
    fh.seek(0)

    # Retrieves metadata from file
    for line in fh:
        if line.startswith("#=GF"):
            metadata = _parse_stockholm_line_gf(line, metadata)
        elif line.startswith("#=GS"):
            dna_data = _parse_stockholm_line_gs(line, dna_data)
        elif line.startswith("#=GR"):
            dna_data = _parse_stockholm_line_gr(line, dna_data)
        elif line.startswith('#=GC'):
            positional_metadata = _parse_stockholm_line_gc(line,
                                                           positional_metadata)
        elif line.startswith('//'):
            break

    for key in dna_data.keys():
        # Sets blank dictionaries and lists to None instead
        # Note: _replace is not a private function, see
        # https://docs.python.org/2/library/collections.html#namedtuple-
        # factory-function-for-tuples-with-named-fields
        if not dna_data[key].metadata:
            dna_data[key] = dna_data[key]._replace(metadata=None)
        if not dna_data[key].positional_metadata:
            dna_data[key] = dna_data[key]._replace(positional_metadata=None)
        # Adds each sequence to the MSA data
        seqs.append(
            constructor(
                dna_data[key].seq,
                metadata=dna_data[key].metadata,
                positional_metadata=(dna_data[key].positional_metadata)))

    if not positional_metadata:
        positional_metadata = None
    # Constructs TabularMSA
    return TabularMSA(seqs, metadata=metadata,
                      positional_metadata=positional_metadata,
                      index=list(dna_data.keys()))


def _parse_stockholm_line_gf(line, metadata):
    """Takes ``#=GF`` line and returns parsed data."""
    line = _remove_newline(line.split(' ', 2))
    _check_for_malformed_line(len(line), 3)
    gf_feature = line[1]
    gf_feature_data = line[2]
    # Handles first instance of labelled tree (New Hampshire eXtended format)
    if gf_feature == 'TN' and 'NH' not in metadata.keys():
        metadata['NH'] = OrderedDict()
        metadata['NH'][gf_feature_data] = ''
        return metadata

    # Handles second instance of labelled tree (NHX)
    if gf_feature == 'TN' and 'NH' in metadata.keys():
        if gf_feature_data in metadata['NH'].keys():
            raise StockholmFormatError("Tree name %r used multiple times in "
                                       "file." % gf_feature_data)
        metadata['NH'][gf_feature_data] = ''
        return metadata

    # Handles extra line(s) of an already created tree (NHX)
    if gf_feature == 'NH' and gf_feature in metadata.keys():
        trees = metadata[gf_feature]
        tree_id =  list(trees.keys())[-1]
        metadata[gf_feature][tree_id] = trees[tree_id] + gf_feature_data
        return metadata

    if gf_feature in metadata.keys():
        metadata[gf_feature] = metadata[gf_feature] + gf_feature_data
    else:
        metadata[gf_feature] = gf_feature_data
    return metadata


def _parse_stockholm_line_gs(line, dna_data):
    """Takes ``#=GS`` line and returns parsed data."""
    line = _remove_newline(line.split(' ', 3))
    _check_for_malformed_line(len(line), 4)
    data_seq_name = line[1]
    gs_feature = line[2]
    if data_seq_name in dna_data.keys():
        existing_metadata = dna_data[data_seq_name].metadata
        if not existing_metadata or gs_feature not in existing_metadata.keys():
            dna_data[data_seq_name].metadata[gs_feature] = line[3]
        else:
            existing_data = dna_data[data_seq_name].metadata[gs_feature]
            dna_data[data_seq_name].metadata[gs_feature] = (existing_data +
                                                            line[3])
    else:
        raise StockholmFormatError("Markup line references nonexistent "
                                   "data %r." % data_seq_name)
    return dna_data


def _parse_stockholm_line_gr(line, dna_data):
    """Takes ``#=GR`` line and returns parsed data."""
    line = _remove_newline(line.split())
    _check_for_malformed_line(len(line), 4)
    data_seq_name = line[1]
    gr_feature = line[2]
    if data_seq_name in dna_data.keys():
        if gr_feature in dna_data[data_seq_name].positional_metadata.keys():
            _raise_duplicate_error("Found duplicate GR label %r associated"
                                   " with data label %r" % (gr_feature,
                                                            data_seq_name))
        dna_data[data_seq_name].positional_metadata[gr_feature] = list(line[3])
    else:
        raise StockholmFormatError("Markup line references nonexistent "
                                   "data %r." % data_seq_name)
    return dna_data


def _parse_stockholm_line_gc(line, positional_metadata):
    """Takes ``#=GC`` line and returns parsed data."""
    line = _remove_newline(line.split())
    _check_for_malformed_line(len(line), 3)
    gc_feature = line[1]
    if gc_feature in positional_metadata.keys():
        _raise_duplicate_error("Found duplicate GC label %r." % (gc_feature))
    positional_metadata[gc_feature] = list(line[2])
    return positional_metadata


def _parse_stockholm_line_data(line, dna_data):
    """Takes data line and returns parsed data."""
    line = line.split()
    _check_for_malformed_line(len(line), 2)
    data_seq_name = line[0]
    if data_seq_name not in dna_data.keys():
        dna_data[data_seq_name] = _SeqData(seq=line[1], metadata={},
                                           positional_metadata={})
    elif data_seq_name in dna_data.keys():
        _raise_duplicate_error("Found multiple data lines under same "
                               "name: %r" % data_seq_name)
    return dna_data


def _remove_newline(line):
    """Removes '\n' from line and returns line."""
    n_line = line[len(line)-1]
    if '\n' in n_line:
        line[len(line)-1] = n_line.rstrip('\n')
    return line


def is_data_line(line):
    return not (line.startswith("#") or line.startswith("//") or
                line.isspace())

def _raise_duplicate_error(message):
    raise StockholmFormatError(message+' Note: If the file being used is in '
                                       'Stockholm interleaved format, this '
                                       'is not supported by the reader.')

def _check_for_malformed_line(num_present, num_needed):
    if num_present != num_needed:
        raise StockholmFormatError('Line only contains %r item(s). It must '
                                   'contain %r.' % (num_present, num_needed))
