"""
Stockholm format (:mod:`skbio.io.format.stockholm`)
===================================================

.. currentmodule:: skbio.io.format.stockholm

The Stockholm format is a multiple sequence alignment format (MSA) that
optionally supports storing arbitrary alignment features (metadata). Features
can be placed into four different categories: GF, GS, GR, and GC (described in
more detail below).

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
--------------
**Has Sniffer: Yes**

**State: Experimental as of 0.4.1-dev.**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`skbio.alignment.TabularMSA`                              |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
The Stockholm format consists of a header, a multiple sequence alignment,
associated metadata (features), and a footer.

Header
^^^^^^
The first line of a Stockholm file must be the following header:

.. code-block:: none

   # STOCKHOLM 1.0

Multiple Sequence Alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Sequence lines consist of a sequence name, followed by whitespace, followed by
the aligned sequence. For example::

    seq1 ACG-T-GGT
    seq2 ACCGTTCG-

.. note:: scikit-bio currently supports reading Stockholm files where each
   sequence is contained on a single line. Interleaved/wrap-around Stockholm
   files are not supported.

.. warning:: Sequence names must be unique.

Metadata
^^^^^^^^
Stockholm files support storing arbitrary metadata (features) about the MSA.
All metadata described in the following sections are optional and may appear in
any order. Metadata "mark-up" lines begin with either ``#=GF``, ``#=GS``,
``#=GR``, or ``#=GC``, and each line describes a single feature of the
alignment.

.. note:: Stockholm format supports generic features. [1]_ and [2]_ provide a
   list of common features output by Pfam/Rfam. scikit-bio does not
   require that these features are present. These features are processed in the
   same way as any arbitrary feature would be, as a simple key-value pair of
   strings.

GF metadata
+++++++++++
Data relating to the multiple sequence alignment as a whole, such as authors or
number of sequences in the alignment. Starts with ``#=GF`` followed by a
feature name and data relating to the feature. Typically comes first in a
Stockholm file.

For example (taken from [2]_):

.. code-block:: none

    #=GF DE CBS domain

Where ``DE`` is the feature name and ``CBS Domain`` is the feature data.

GF metadata is stored in the ``TabularMSA`` ``metadata`` dictionary.

.. note:: Duplicate GF feature names will have their values concatenated in the
   order they appear in the file.

.. note:: Trees labelled with ``NH``/``TN`` are handled differently than other
   GF features. When reading a Stockholm file with these features, the reader
   follows the rules described in [2]_.

   A single tree without an identifier will be stored as::

       metadata = {
           'NH': 'tree in NHX format'
       }

   A single tree with an identifier will be stored as::

       metadata = {
           'NH': {
               'tree-id': 'tree in NHX format'
           }
       }

   Multiple trees (which must have identifiers) will be stored as::

       metadata = {
           'NH': {
               'tree-id-1': 'tree in NHX format',
               'tree-id-2': 'tree in NHX format'
           }
       }

GS metadata
+++++++++++
Data relating to a specific sequence in the multiple sequence alignment.
Starts with ``#=GS`` followed by the sequence name followed by a feature name
and data relating to the feature. Typically comes after GF metadata in a
Stockholm file.

For example (taken from [2]_):

.. code-block:: none

    #=GS O83071/259-312 AC O83071

Where ``O83071/259-312`` is the sequence name, ``AC`` is the feature name, and
``083071`` is the feature data.

GS metadata is stored in the sequence-specific ``metadata`` dictionary.

.. note:: Duplicate GS feature names will have their values concatenated in the
   order they appear in the file.

GR metadata
+++++++++++
Data relating to the columns of a specific sequence in a multiple sequence
alignment. Starts with ``#=GR`` followed by the sequence name followed by a
feature name and data relating to the feature, one character per column.
Typically comes after the sequence line it relates to.

For example (taken from [2]_):

.. code-block:: none

    #=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEHHH

Where ``O31698/18-71`` is the sequence name, ``SS`` is the feature name, and
``CCCHHHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEHHH`` is the feature data.

GR metadata is stored in sequence-specific ``positional_metadata``.

.. note:: Duplicate GR feature names attributed to a single sequence are
   disallowed.

GC metadata
+++++++++++
Data relating to the columns of the multiple sequence alignment as a whole.
Starts with ``#=GC`` followed by a feature name and data relating to the
feature, one character per column. Typically comes at the end of the multiple
sequence alignment.

For example (taken from [2]_):

.. code-block:: none

    #=GC SS_cons            CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEH

Where ``SS_cons`` is the feature name and
``CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEH`` is the feature data.

GC metadata is stored in ``TabularMSA`` ``positional_metadata``.

.. note:: Duplicate GC feature names are disallowed.

Footer
^^^^^^
The final line of a Stockholm file must be the following footer::

    //

.. note:: scikit-bio currently supports reading a Stockholm file containing a
   single MSA. If the file contains more than one MSA, only the first MSA will
   be read into a ``TabularMSA``.

Format Parameters
-----------------
The only supported format parameter is ``constructor``, which specifies the
type of in-memory sequence object to read each aligned sequence into. This must
be a subclass of ``GrammaredSequence`` (e.g., ``DNA``, ``RNA``, ``Protein``)
and is a required format parameter. For example, if you know that the Stockholm
file you're reading contains DNA sequences, you would pass ``constructor=DNA``
to the reader call.

Examples
--------
Suppose we have a Stockholm file containing an MSA of protein sequences
(modified from [2]_):

>>> import skbio.io
>>> from io import StringIO
>>> from skbio import Protein, TabularMSA
>>> fs = '\\n'.join([
...         '# STOCKHOLM 1.0',
...         '#=GF CC CBS domains are small intracellular modules mostly'
...         ' found ',
...         '#=GF CC in 2 or four copies within a protein.',
...         '#=GS O83071/192-246 AC O83071',
...         '#=GS O31698/88-139 OS Bacillus subtilis',
...         'O83071/192-246          MTCRAQLIAVPRASSLAE..AIACAQKM....RVSRV',
...         '#=GR O83071/192-246 SA  999887756453524252..55152525....36463',
...         'O83071/259-312          MQHVSAPVFVFECTRLAY..VQHKLRAH....SRAVA',
...         'O31698/18-71            MIEADKVAHVQVGNNLEH..ALLVLTKT....GYTAI',
...         'O31698/88-139           EVMLTDIPRLHINDPIMK..GFGMVINN......GFV',
...         'O31699/88-139           EVMLTDIPRLHINDPIMK..GFGMVINN......GFV',
...         '#=GR O31699/88-139 AS   ________________*____________________',
...         '#=GR O31699/88-139 IN   ____________1______________2_________',
...         '#=GC SS_cons            CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEE',
...         '//'
...     ])
>>> fh = StringIO(fs)
>>> msa = TabularMSA.read(fh, constructor=Protein)
>>> msa # doctest: +NORMALIZE_WHITESPACE
TabularMSA[Protein]
----------------------------------------------------------------------
Metadata:
    'CC': 'CBS domains are small intracellular modules mostly found in
           2 or four copies within a protein.'
Positional metadata:
    'SS_cons': <dtype: object>
Stats:
    sequence count: 5
    position count: 37
----------------------------------------------------------------------
MTCRAQLIAVPRASSLAE..AIACAQKM....RVSRV
MQHVSAPVFVFECTRLAY..VQHKLRAH....SRAVA
MIEADKVAHVQVGNNLEH..ALLVLTKT....GYTAI
EVMLTDIPRLHINDPIMK..GFGMVINN......GFV
EVMLTDIPRLHINDPIMK..GFGMVINN......GFV

The ``TabularMSA`` has GF metadata stored in its ``metadata`` dictionary:

>>> msa.metadata
{'CC': 'CBS domains are small intracellular modules mostly found in 2 or four \
copies within a protein.'}

GC metadata is stored in the ``TabularMSA`` ``positional_metadata``:

>>> msa.positional_metadata  # doctest: +ELLIPSIS
   SS_cons
0        C
1        C
2        C
3        C
4        C
5        H
6        H
7        H
8        H
9        H
...

GS metadata is stored in the sequence-specific ``metadata`` dictionary:

>>> msa[0].metadata
{'AC': 'O83071'}

GR metadata is stored in sequence-specific ``positional_metadata``:

>>> msa[4].positional_metadata  # doctest: +ELLIPSIS
   AS IN
0   _  _
1   _  _
2   _  _
3   _  _
4   _  _
5   _  _
6   _  _
7   _  _
8   _  _
9   _  _
...

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
from skbio.sequence._grammared_sequence import GrammaredSequence
from skbio.io import create_format, StockholmFormatError

stockholm = create_format('stockholm')


@stockholm.sniffer()
def _stockholm_sniffer(fh):
    # Smells a Stockholm file if the following conditions are met:
    # - File isn't empty
    # - File contains correct header
    try:
        line = next(fh)
    except StopIteration:
        return False, {}

    if _is_header(line):
        return True, {}

    return False, {}


@stockholm.reader(TabularMSA)
def _stockholm_to_tabular_msa(fh, constructor=None):
    # Checks that user has passed required constructor parameter
    if constructor is None:
        raise ValueError("Must provide `constructor` parameter.")
    # Checks that contructor parameter is supported
    elif not issubclass(constructor, GrammaredSequence):
        raise TypeError("`constructor` must be a subclass of"
                        "`GrammaredSequence`.")

    # Checks that the file isn't empty
    try:
        line = next(fh)
    except StopIteration:
        raise StockholmFormatError("File is empty.")
    # Checks that the file follows basic format (includes the required header)
    if not _is_header(line):
        raise StockholmFormatError("File missing required Stockholm header "
                                   "line.")

    msa_data = _MSAData()
    for line in fh:
        if line.isspace():
            continue

        line = line.rstrip('\n')

        if _is_sequence_line(line):
            seq_name, seq_data = _parse_sequence_line(line)
            msa_data.add_sequence(seq_name, seq_data)
        elif line.startswith("#=GF"):
            feature_name, feature_data = _parse_gf_line(line)
            msa_data.add_gf_metadata(feature_name, feature_data)
        elif line.startswith("#=GS"):
            seq_name, feature_name, feature_data = _parse_gs_line(line)
            msa_data.add_gs_metadata(seq_name, feature_name, feature_data)
        elif line.startswith("#=GR"):
            seq_name, feature_name, feature_data = _parse_gr_line(line)
            msa_data.add_gr_metadata(seq_name, feature_name, feature_data)
        elif line.startswith('#=GC'):
            feature_name, feature_data = _parse_gc_line(line)
            msa_data.add_gc_metadata(feature_name, feature_data)
        elif _is_footer(line):
            break
        else:
            raise StockholmFormatError("Unrecognized line: %r" % line)

    if not _is_footer(line):
        raise StockholmFormatError('Final line does not conform to Stockholm '
                                   'format. Must contain only "//".')

    return msa_data.build_tabular_msa(constructor)


# For storing intermediate data used to construct a Sequence object.
class _MSAData(object):
    def __init__(self):
        self._seqs = {}
        self._seq_order = []
        self._metadata = {}
        self._positional_metadata = {}

    def add_sequence(self, seq_name, seq_data):
        if seq_name not in self._seqs:
            self._seqs[seq_name] = _SeqData(seq_name)
        self._seqs[seq_name].seq = seq_data
        self._seq_order.append(seq_name)

    def add_gf_metadata(self, feature_name, feature_data):
        # Handles first instance of labelled tree
        if feature_name == 'TN' and 'NH' not in self._metadata:
            self._metadata['NH'] = OrderedDict()
            self._metadata['NH'][feature_data] = ''
        # Handles second instance of labelled tree
        elif feature_name == 'TN' and 'NH' in self._metadata:
            if feature_data in self._metadata['NH']:
                raise StockholmFormatError("Tree name %r used multiple times "
                                           "in file." % feature_data)
            self._metadata['NH'][feature_data] = ''
        # Handles extra line(s) of an already created tree
        elif feature_name == 'NH' and feature_name in self._metadata:
            trees = self._metadata[feature_name]
            tree_id = list(trees.keys())[-1]
            self._metadata[feature_name][tree_id] = (trees[tree_id] +
                                                     feature_data)
        elif feature_name in self._metadata:
            self._metadata[feature_name] = (self._metadata[feature_name] +
                                            feature_data)
        else:
            self._metadata[feature_name] = feature_data

    def add_gc_metadata(self, feature_name, feature_data):
        if feature_name in self._positional_metadata:
            _raise_duplicate_error("Found duplicate GC label %r."
                                   % feature_name)
        self._positional_metadata[feature_name] = feature_data

    def add_gs_metadata(self, seq_name, feature_name, feature_data):
        if seq_name not in self._seqs:
            self._seqs[seq_name] = _SeqData(seq_name)
        self._seqs[seq_name].add_metadata_feature(feature_name, feature_data)

    def add_gr_metadata(self, seq_name, feature_name, feature_data):
        if seq_name not in self._seqs:
            self._seqs[seq_name] = _SeqData(seq_name)
        self._seqs[seq_name].add_positional_metadata_feature(feature_name,
                                                             feature_data)

    def build_tabular_msa(self, constructor):
        if len(self._seqs) != len(self._seq_order):
            invalid_seq_names = set(self._seqs) - set(self._seq_order)
            raise StockholmFormatError('Found GS or GR metadata for '
                                       'nonexistent sequence(s): %r'
                                       % invalid_seq_names)

        seqs = []
        for seq_name in self._seq_order:
            seqs.append(self._seqs[seq_name].build_sequence(constructor))

        positional_metadata = self._positional_metadata
        if not positional_metadata:
            positional_metadata = None

        metadata = self._metadata
        if not metadata:
            metadata = None

        # Constructs TabularMSA
        return TabularMSA(seqs, metadata=metadata,
                          positional_metadata=positional_metadata,
                          index=self._seq_order)


class _SeqData(object):
    def __init__(self, name):
        self.name = name
        self._seq = None
        self.metadata = None
        self.positional_metadata = None

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, seq):
        if self._seq is None:
            self._seq = seq
        else:
            _raise_duplicate_error("Found duplicate sequence name: %r"
                                   % self.name)

    def add_metadata_feature(self, feature_name, feature_data):
        if self.metadata is None:
            self.metadata = {}
        if feature_name in self.metadata:
            self.metadata[feature_name] += feature_data
        else:
            self.metadata[feature_name] = feature_data

    def add_positional_metadata_feature(self, feature_name, feature_data):
        if self.positional_metadata is None:
            self.positional_metadata = {}
        if feature_name in self.positional_metadata:
            _raise_duplicate_error("Found duplicate GR label %r associated "
                                   "with sequence name %r"
                                   % (feature_name, self.name))
        else:
            self.positional_metadata[feature_name] = feature_data

    def build_sequence(self, constructor):
        return constructor(self.seq, metadata=self.metadata,
                           positional_metadata=(self.positional_metadata))


def _parse_gf_line(line):
    line = line.split(None, 2)
    _check_for_malformed_line(line, 3)
    return line[1:]


def _parse_gs_line(line):
    line = line.split(None, 3)
    _check_for_malformed_line(line, 4)
    return line[1:]


def _parse_gr_line(line):
    line = line.split(None, 3)
    _check_for_malformed_line(line, 4)
    seq_name = line[1]
    feature_name = line[2]
    feature_data = list(line[3])
    return seq_name, feature_name, feature_data


def _parse_gc_line(line):
    line = line.split(None, 2)
    _check_for_malformed_line(line, 3)
    feature_name = line[1]
    feature_data = list(line[2])
    return feature_name, feature_data


def _parse_sequence_line(line):
    line = line.split(None, 1)
    _check_for_malformed_line(line, 2)
    return line


def _is_header(line):
    return line == '# STOCKHOLM 1.0\n'


def _is_footer(line):
    return line.rstrip() == '//'


def _is_sequence_line(line):
    return not (line.startswith("#") or _is_footer(line))


def _raise_duplicate_error(message):
    raise StockholmFormatError(message+' Note: If the file being used is in '
                                       'Stockholm interleaved format, this '
                                       'is not supported by the reader.')


def _check_for_malformed_line(line, expected_len):
    if len(line) != expected_len:
        raise StockholmFormatError('Line contains %d item(s). It must '
                                   'contain exactly %d item(s).'
                                   % (len(line), expected_len))
