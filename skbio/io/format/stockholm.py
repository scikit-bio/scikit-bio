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

**State: Experimental as of 0.4.2.**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.alignment.TabularMSA`                              |
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

Sequence names (``seq1``, ``seq2``) are stored in the ``TabularMSA``
``index``.

.. note:: scikit-bio currently supports reading Stockholm files where each
   sequence is contained on a single line. Interleaved/wrap-around Stockholm
   files are not supported. When writing, each sequence will be placed on its
   own line.

.. warning:: Sequence names must be unique in the Stockholm file. Likewise,
   when writing from a ``TabularMSA``, ``index`` must be unique.

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
   strings. When writing, feature names, feature data, and sequence names are
   converted to type ``str``.

.. note:: When writing a Stockholm file, scikit-bio will place the metadata in
   the format's recommended order:

   - GF: Above the alignment
   - GS: Above the alignment (after GF)
   - GR: Below corresponding sequence
   - GC: Below the alignment

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

.. note:: When reading, duplicate GF feature names will have their values
   concatenated in the order they appear in the file. Concatenation will
   also add a space between lines if one isn't already there in order to avoid
   joining words together. When writing, each GF feature will be placed on its
   own line, regardless of length.

.. note:: Trees labelled with ``NH``/``TN`` are handled differently than other
   GF features. When reading a Stockholm file with these features, the reader
   follows the rules described in [2]_. Trees split over multiple lines will
   have their values concatenated. Unlike other GF features, trees will never
   have a space added when they are concatenated.

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

.. note:: References labelled with ``RN``/``RM``/``RT``/``RA``/``RL``/``RC``
   are handled differently than other GF features. When reading a Stockholm
   file with these features, the reader populates a list of dictionaries,
   where each dictionary represents a single reference. The list contains
   references in the order they appear in the file, regardless of the value
   provided for ``RN``. If a reference does not include all possible reference
   tags (e.g. ``RC`` is missing), the dictionary will only contain the
   reference tags present for that reference. When writing, the writer adds a
   reference number (``RN``) line before writing each reference, for example:

   .. code-block:: none

      #=GF RN [1]
      #=GF RA Kestrel Gorlick
      ...
      #=GF RN [2]
      ...

   References will be stored as::

       metadata = {
           'RN': [{
               'RM': 'reference medline',
               'RT': 'reference title',
               'RA': 'reference author',
               'RL': 'reference location',
               'RC': 'reference comment'
           }, {
               'RM': 'reference medline',
               ...
           }]
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

.. note:: When reading, duplicate GS feature names will have their values
   concatenated in the order they appear in the file. Concatenation will
   also add a space between lines if one isn't already there in order to avoid
   joining words together. When writing, each GS feature will be placed on its
   own line, regardless of length.

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
...         ' found',
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

The sequence names are stored in the ``index``:

>>> msa.index
Index(['O83071/192-246', 'O83071/259-312', 'O31698/18-71', 'O31698/88-139',
       'O31699/88-139'],
      dtype='object')

The ``TabularMSA`` has GF metadata stored in its ``metadata`` dictionary:

>>> msa.metadata
OrderedDict([('CC', 'CBS domains are small intracellular modules mostly found \
in 2 or four copies within a protein.')])

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
OrderedDict([('AC', 'O83071')])

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

Let's write this ``TabularMSA`` in Stockholm format:

>>> fh = StringIO()
>>> _ = msa.write(fh, format='stockholm')
>>> print(fh.getvalue())
# STOCKHOLM 1.0
#=GF CC CBS domains are small intracellular modules mostly found in 2 or four \
copies within a protein.
#=GS O83071/192-246 AC O83071
#=GS O31698/88-139 OS Bacillus subtilis
O83071/192-246         MTCRAQLIAVPRASSLAE..AIACAQKM....RVSRV
#=GR O83071/192-246 SA 999887756453524252..55152525....36463
O83071/259-312         MQHVSAPVFVFECTRLAY..VQHKLRAH....SRAVA
O31698/18-71           MIEADKVAHVQVGNNLEH..ALLVLTKT....GYTAI
O31698/88-139          EVMLTDIPRLHINDPIMK..GFGMVINN......GFV
O31699/88-139          EVMLTDIPRLHINDPIMK..GFGMVINN......GFV
#=GR O31699/88-139 AS  ________________*____________________
#=GR O31699/88-139 IN  ____________1______________2_________
#=GC SS_cons           CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEE
//
<BLANKLINE>
>>> fh.close()

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

from collections import OrderedDict

from skbio.alignment import TabularMSA
from skbio.sequence._grammared_sequence import GrammaredSequence
from skbio.io import create_format, StockholmFormatError

stockholm = create_format('stockholm')
_REFERENCE_TAGS = frozenset({'RM', 'RT', 'RA', 'RL', 'RC'})


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
        raise ValueError("Must provide `constructor` parameter indicating the "
                         "type of sequences in the alignment. `constructor` "
                         "must be a subclass of `GrammaredSequence` "
                         "(e.g., `DNA`, `RNA`, `Protein`).")
    # Checks that contructor parameter is supported
    elif not issubclass(constructor, GrammaredSequence):
        raise TypeError("`constructor` must be a subclass of "
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
class _MSAData:
    def __init__(self):
        self._seqs = {}
        self._seq_order = []
        self._metadata = OrderedDict()
        self._positional_metadata = OrderedDict()

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
            if isinstance(trees, OrderedDict):
                tree_id = next(reversed(trees))
                self._metadata[feature_name][tree_id] = (trees[tree_id] +
                                                         feature_data)
            else:
                self._metadata[feature_name] = (self._metadata[feature_name] +
                                                feature_data)
        elif feature_name == 'RN':
            if feature_name not in self._metadata:
                self._metadata[feature_name] = [OrderedDict()]
            else:
                self._metadata[feature_name].append(OrderedDict())
        elif feature_name in _REFERENCE_TAGS:
            if 'RN' not in self._metadata:
                raise StockholmFormatError("Expected 'RN' tag to precede "
                                           "'%s' tag." % feature_name)
            reference_dict = self._metadata['RN'][-1]
            if feature_name not in reference_dict:
                reference_dict[feature_name] = feature_data
            else:
                padding = _get_padding(reference_dict[feature_name])
                reference_dict[feature_name] += padding + feature_data
        elif feature_name in self._metadata:
            padding = _get_padding(self._metadata[feature_name][-1])
            self._metadata[feature_name] = (self._metadata[feature_name] +
                                            padding + feature_data)
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


class _SeqData:
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
            self.metadata = OrderedDict()
        if feature_name in self.metadata:
            padding = _get_padding(self.metadata[feature_name][-1])
            self.metadata[feature_name] += padding + feature_data
        else:
            self.metadata[feature_name] = feature_data

    def add_positional_metadata_feature(self, feature_name, feature_data):
        if self.positional_metadata is None:
            self.positional_metadata = OrderedDict()
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


@stockholm.writer(TabularMSA)
def _tabular_msa_to_stockholm(obj, fh):
    if not obj.index.is_unique:
        raise StockholmFormatError("The TabularMSA's index labels must be"
                                   " unique.")
    # Writes header
    fh.write("# STOCKHOLM 1.0\n")

    # Writes GF data to file
    if obj.has_metadata():
        for gf_feature, gf_feature_data in obj.metadata.items():
            if gf_feature == 'NH' and isinstance(gf_feature_data, dict):
                for tree_id, tree in gf_feature_data.items():
                    fh.write("#=GF TN %s\n" % tree_id)
                    fh.write("#=GF NH %s\n" % tree)
            elif gf_feature == 'RN':
                if not isinstance(gf_feature_data, list):
                    raise StockholmFormatError(
                        "Expected 'RN' to contain a list of reference "
                        "dictionaries, got %r." % gf_feature_data)

                for ref_num, dictionary in enumerate(gf_feature_data, start=1):
                    if not isinstance(dictionary, dict):
                        raise StockholmFormatError(
                            "Expected reference information to be stored as a "
                            "dictionary, found reference %d stored as %r." %
                            (ref_num, type(dictionary).__name__))

                    fh.write("#=GF RN [%d]\n" % ref_num)
                    for feature in dictionary:
                        if feature not in _REFERENCE_TAGS:
                            formatted_reference_tags = ', '.join(
                                [tag for tag in _REFERENCE_TAGS])
                            raise StockholmFormatError(
                                "Invalid reference tag %r found in reference "
                                "dictionary %d. Valid reference tags are: %s."
                                % (feature, ref_num, formatted_reference_tags))

                        fh.write("#=GF %s %s\n" % (feature,
                                                   dictionary[feature]))
            else:
                fh.write("#=GF %s %s\n" % (gf_feature, gf_feature_data))

    unpadded_data = []
    # Writes GS data to file, retrieves GR data, and retrieves sequence data
    for seq, seq_name in zip(obj, obj.index):
        seq_name = str(seq_name)

        if seq.has_metadata():
            for gs_feature, gs_feature_data in seq.metadata.items():
                fh.write("#=GS %s %s %s\n" % (seq_name, gs_feature,
                                              gs_feature_data))

        unpadded_data.append((seq_name, str(seq)))
        if seq.has_positional_metadata():
            df = _format_positional_metadata(seq.positional_metadata,
                                             'Sequence-specific positional '
                                             'metadata (GR)')
            for gr_feature in df.columns:
                gr_feature_data = ''.join(df[gr_feature])
                gr_string = "#=GR %s %s" % (seq_name, gr_feature)
                unpadded_data.append((gr_string, gr_feature_data))

    # Retrieves GC data
    if obj.has_positional_metadata():
        df = _format_positional_metadata(obj.positional_metadata,
                                         'Multiple sequence alignment '
                                         'positional metadata (GC)')
        for gc_feature in df.columns:
            gc_feature_data = ''.join(df[gc_feature])
            gc_string = "#=GC %s" % gc_feature
            unpadded_data.append((gc_string, gc_feature_data))

    # Writes GR, GC, and raw data to file with padding
    _write_padded_data(unpadded_data, fh)

    # Writes footer
    fh.write("//\n")


def _write_padded_data(data, fh):
    max_data_len = 0
    for label, _ in data:
        if len(label) > max_data_len:
            max_data_len = len(label)
    fmt = '{0:%d} {1}\n' % max_data_len
    for label, value in data:
        fh.write(fmt.format(label, value))


def _format_positional_metadata(df, data_type):
    # Asserts positional metadata feature names are unique
    if not df.columns.is_unique:
        num_repeated_columns = len(df.columns) - len(set(df.columns))
        raise StockholmFormatError('%s feature names must be unique. '
                                   'Found %d duplicate names.'
                                   % (data_type, num_repeated_columns))

    str_df = df.astype(str)

    # Asserts positional metadata dataframe items are one character long
    for column in str_df.columns:
        if (str_df[column].str.len() != 1).any():
            raise StockholmFormatError("%s must contain a single character for"
                                       " each position's value. Found value(s)"
                                       " in column %s of incorrect length."
                                       % (data_type, column))
    return str_df


def _get_padding(item):
    return '' if item[-1].isspace() else ' '
