"""GenBank format (:mod:`skbio.io.format.genbank`)
===============================================

.. currentmodule:: skbio.io.format.genbank

GenBank format (GenBank Flat File Format) stores sequence and its
annotation together. The start of the annotation section is marked by
a line beginning with the word "LOCUS". The start of sequence section
is marked by a line beginning with the word "ORIGIN" and the end of
the section is marked by a line with only "//".

The GenBank file usually ends with .gb or sometimes .gbk. The GenBank
format for protein has been renamed to GenPept. The GenBank (for
nucleotide) and Genpept are essentially the same format. An example of
a GenBank file can be seen here [1]_.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.RNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.Protein`                                  |
+------+------+---------------------------------------------------------------+
|Yes   | Yes  | generator of :mod:`skbio.sequence.Sequence` objects           |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------

Sections before ``FEATURES``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All the sections before ``FEATURES`` will be read into the attribute
of ``metadata``. The header and its content of a section is stored as
a pair of key and value in ``metadata``. For the ``REFERENCE``
section, its value is stored as a list, as there are often multiple
reference sections in one GenBank record.

.. _genbank_feature_section:

``FEATURES`` section
^^^^^^^^^^^^^^^^^^^^
The International Nucleotide Sequence Database Collaboration (INSDC
[2]_) is a joint effort among the DDBJ, EMBL, and GenBank. These
organisations all use the same "Feature Table" layout in their plain
text flat file formats, which are documented in detail [3]_. The
feature keys and their qualifiers are also described in this webpage
[4]_.

The ``FEATURES`` section will be stored in ``interval_metadata`` of
``Sequence`` or its sub-class. Each sub-section is stored as an
``Interval`` object in ``interval_metadata``. Each ``Interval`` object
has ``metadata`` keeping the information of this feature in the
sub-section.

To normalize the vocabulary between multiple formats (currently only
the INSDC Feature Table and GFF3) to store metadata of interval
features, we rename some terms in some formats to the same common name
when parsing them into memory, as described in this table:

+-----------+-----------+-----------+---------+------------------------------+
|INSDC      |GFF3       |Key stored |Value    |Description                   |
|feature    |columns or |           |type     |                              |
|table      |attributes |           |stored   |                              |
+===========+===========+===========+=========+==============================+
|inference  |source     |source     |str      |the algorithm or experiment   |
|           |(column 2) |           |         |used to generate this feature |
+-----------+-----------+-----------+---------+------------------------------+
|feature key|type       |type       |str      |the type of the feature       |
|           |(column 3) |           |         |                              |
+-----------+-----------+-----------+---------+------------------------------+
|N/A        |score      |score      |float    |the score of the feature      |
|           |(column 6) |           |         |                              |
+-----------+-----------+-----------+---------+------------------------------+
|N/A        |strand     |strand     |str      |the strand of the feature. +  |
|           |(column 7) |           |         |for positive strand, - for    |
|           |           |           |         |minus strand, and . for       |
|           |           |           |         |features that are not         |
|           |           |           |         |stranded. In addition, ?  can |
|           |           |           |         |be used for features whose    |
|           |           |           |         |strandedness is relevant, but |
|           |           |           |         |unknown.                      |
+-----------+-----------+-----------+---------+------------------------------+
|codon_start|phase      |phase      |int      |the offset at which the first |
|           |(column 8) |           |         |complete codon of a coding    |
|           |           |           |         |feature can be found, relative|
|           |           |           |         |to the first base of that     |
|           |           |           |         |feature. It is 0, 1, or 2 in  |
|           |           |           |         |GFF3 or 1, 2, or 3 in GenBank.|
|           |           |           |         |The stored value is 0, 1, or  |
|           |           |           |         |2, following in GFF3 format.  |
+-----------+-----------+-----------+---------+------------------------------+
|db_xref    |Dbxref     |db_xref    |list of  |A database cross reference    |
|           |           |           |str      |                              |
+-----------+-----------+-----------+---------+------------------------------+
|N/A        |ID         |ID         |str      |feature ID                    |
+-----------+-----------+-----------+---------+------------------------------+
|note       |Note       |note       |str      |any comment or additional     |
|           |           |           |         |information                   |
+-----------+-----------+-----------+---------+------------------------------+
|translation|N/A        |translation|str      |the protein sequence for CDS  |
|           |           |           |         |features                      |
+-----------+-----------+-----------+---------+------------------------------+

``Location`` string
+++++++++++++++++++
There are 5 types of location descriptors defined in Feature
Table. This explains how they will be parsed into the bounds of
``Interval`` object (note it converts the 1-based coordinate to
0-based):

    1. a single base number. e.g. 67. This is parsed to ``(66, 67)``.

    2. a site between two neighboring bases. e.g. 67^68. This
       is parsed to ``(66, 67)``.

    3. a single base from inside a range. e.g. 67.89. This is parsed to
       ``(66, 89)``.

    4. a pair of base numbers defining a sequence span. e.g. 67..89. This
       is parsed to ``(66, 89)``.

    5. a remote sequence identifier followed by a location descriptor
       defined above. e.g. J00123.1:67..89. This will be discarded
       because it is not on the current sequence. When it is combined
       with local descriptor like J00123.1:67..89,200..209, the
       local part will be kept to be ``(199, 209)``.


.. note:: The Location string is fully stored in ``Interval.metadata``
   with key ``__location``.  The key starting with ``__`` is "private"
   and should be modified with care.


``ORIGIN`` section
^^^^^^^^^^^^^^^^^^
The sequence in the ``ORIGIN`` section is always in lowercase for
the GenBank files downloaded from NCBI. For the RNA molecules, ``t``
(thymine), instead of ``u`` (uracil) is used in the sequence. All
GenBank writers follow these conventions while writing GenBank files.


Format Parameters
-----------------

Reader-specific Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``constructor`` parameter can be used with the ``Sequence`` generator
to specify the in-memory type of each GenBank record that is parsed.
``constructor`` should be ``Sequence`` or a sub-class of ``Sequence``.
It is also detected by the unit label on the LOCUS line. For example, if it
is ``bp``, it will be read into ``DNA``; if it is ``aa``, it will be read
into ``Protein``. Otherwise, it will be read into ``Sequence``. This default
behavior is overridden by setting ``constructor``.

``lowercase`` is another parameter available for all GenBank readers.
By default, it is set to ``True`` to read in the ``ORIGIN`` sequence
as lowercase letters. This parameter is passed to ``Sequence`` or
its sub-class constructor.

``seq_num`` is a parameter used with the ``Sequence``, ``DNA``, ``RNA``, and
``Protein`` GenBank readers. It specifies which GenBank record to read from
a GenBank file with multiple records in it.

Examples
--------
Reading and Writing GenBank Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Suppose we have the following GenBank file example modified from [5]_:

>>> gb_str = '''
... LOCUS       3K1V_A       34 bp    RNA     linear   SYN 10-OCT-2012
... DEFINITION  Chain A, Structure Of A Mutant Class-I Preq1.
... ACCESSION   3K1V_A
... VERSION     3K1V_A  GI:260656459
... KEYWORDS    .
... SOURCE      synthetic construct
...   ORGANISM  synthetic construct
...             other sequences; artificial sequences.
... REFERENCE   1  (bases 1 to 34)
...   AUTHORS   Klein,D.J., Edwards,T.E. and Ferre-D'Amare,A.R.
...   TITLE     Cocrystal structure of a class I preQ1 riboswitch
...   JOURNAL   Nat. Struct. Mol. Biol. 16 (3), 343-344 (2009)
...    PUBMED   19234468
... COMMENT     SEQRES.
... FEATURES             Location/Qualifiers
...      source          1..34
...                      /organism="synthetic construct"
...                      /mol_type="other RNA"
...                      /db_xref="taxon:32630"
...      misc_binding    1..30
...                      /note="Preq1 riboswitch"
...                      /bound_moiety="preQ1"
... ORIGIN
...         1 agaggttcta gcacatccct ctataaaaaa ctaa
... //
... '''


Now we can read it as ``DNA`` object:

>>> import io
>>> from skbio import DNA, RNA, Sequence
>>> gb = io.StringIO(gb_str)
>>> dna_seq = DNA.read(gb)
>>> dna_seq
DNA
-----------------------------------------------------------------
Metadata:
    'ACCESSION': '3K1V_A'
    'COMMENT': 'SEQRES.'
    'DEFINITION': 'Chain A, Structure Of A Mutant Class-I Preq1.'
    'KEYWORDS': '.'
    'LOCUS': <class 'dict'>
    'REFERENCE': <class 'list'>
    'SOURCE': <class 'dict'>
    'VERSION': '3K1V_A  GI:260656459'
Interval metadata:
    2 interval features
Stats:
    length: 34
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 35.29%
-----------------------------------------------------------------
0 AGAGGTTCTA GCACATCCCT CTATAAAAAA CTAA


Since this is a riboswitch molecule, we may want to read it as
``RNA``.  As the GenBank file usually have ``t`` instead of ``u`` in
the sequence, we can read it as ``RNA`` by converting ``t`` to ``u``:

>>> gb = io.StringIO(gb_str)
>>> rna_seq = RNA.read(gb)
>>> rna_seq
RNA
-----------------------------------------------------------------
Metadata:
    'ACCESSION': '3K1V_A'
    'COMMENT': 'SEQRES.'
    'DEFINITION': 'Chain A, Structure Of A Mutant Class-I Preq1.'
    'KEYWORDS': '.'
    'LOCUS': <class 'dict'>
    'REFERENCE': <class 'list'>
    'SOURCE': <class 'dict'>
    'VERSION': '3K1V_A  GI:260656459'
Interval metadata:
    2 interval features
Stats:
    length: 34
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 35.29%
-----------------------------------------------------------------
0 AGAGGUUCUA GCACAUCCCU CUAUAAAAAA CUAA

>>> rna_seq == dna_seq.transcribe()
True

>>> with io.StringIO() as fh:
...     print(dna_seq.write(fh, format='genbank').getvalue())
LOCUS       3K1V_A   34 bp   RNA   linear   SYN   10-OCT-2012
DEFINITION  Chain A, Structure Of A Mutant Class-I Preq1.
ACCESSION   3K1V_A
VERSION     3K1V_A  GI:260656459
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            other sequences; artificial sequences.
REFERENCE   1  (bases 1 to 34)
  AUTHORS   Klein,D.J., Edwards,T.E. and Ferre-D'Amare,A.R.
  TITLE     Cocrystal structure of a class I preQ1 riboswitch
  JOURNAL   Nat. Struct. Mol. Biol. 16 (3), 343-344 (2009)
  PUBMED    19234468
COMMENT     SEQRES.
FEATURES             Location/Qualifiers
     source          1..34
                     /db_xref="taxon:32630"
                     /mol_type="other RNA"
                     /organism="synthetic construct"
     misc_binding    1..30
                     /bound_moiety="preQ1"
                     /note="Preq1 riboswitch"
ORIGIN
        1 agaggttcta gcacatccct ctataaaaaa ctaa
//
<BLANKLINE>

References
----------
.. [1] http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
.. [2] http://www.insdc.org/
.. [3] http://www.insdc.org/files/feature_table.html
.. [4] http://www.ebi.ac.uk/ena/WebFeat/
.. [5] http://www.ncbi.nlm.nih.gov/nuccore/3K1V_A


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
from functools import partial

from skbio.io import create_format, GenBankFormatError
from skbio.io.format._base import _get_nth_sequence, _line_generator, _too_many_blanks
from skbio.util._misc import chunk_str
from skbio.sequence import Sequence, DNA, RNA, Protein
from skbio.io.format._sequence_feature_vocabulary import (
    _yield_section,
    _parse_section_default,
    _serialize_section_default,
    _parse_feature_table,
    _serialize_feature_table,
)


genbank = create_format("genbank")

# This list is ordered
# used to read and write genbank file.
_HEADERS = [
    "LOCUS",
    "DEFINITION",
    "ACCESSION",
    "VERSION",
    "DBSOURCE",
    "DBLINK",
    "KEYWORDS",
    "SOURCE",
    "REFERENCE",
    "COMMENT",
    "FEATURES",
    "ORIGIN",
]


@genbank.sniffer()
def _genbank_sniffer(fh):
    # check the 1st real line is a valid LOCUS line
    if _too_many_blanks(fh, 5):
        return False, {}
    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    try:
        _parse_locus([line])
    except GenBankFormatError:
        return False, {}
    return True, {}


@genbank.reader(None)
def _genbank_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_genbanks(fh):
        yield _construct(record, constructor, **kwargs)


@genbank.reader(Sequence)
def _genbank_to_sequence(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_genbanks(fh), seq_num)
    return _construct(record, Sequence, **kwargs)


@genbank.reader(DNA)
def _genbank_to_dna(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_genbanks(fh), seq_num)
    return _construct(record, DNA, **kwargs)


@genbank.reader(RNA)
def _genbank_to_rna(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_genbanks(fh), seq_num)
    return _construct(record, RNA, **kwargs)


@genbank.reader(Protein)
def _genbank_to_protein(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_genbanks(fh), seq_num)
    return _construct(record, Protein, **kwargs)


@genbank.writer(None)
def _generator_to_genbank(obj, fh):
    for obj_i in obj:
        _serialize_single_genbank(obj_i, fh)


@genbank.writer(Sequence)
def _sequence_to_genbank(obj, fh):
    _serialize_single_genbank(obj, fh)


@genbank.writer(DNA)
def _dna_to_genbank(obj, fh):
    _serialize_single_genbank(obj, fh)


@genbank.writer(RNA)
def _rna_to_genbank(obj, fh):
    _serialize_single_genbank(obj, fh)


@genbank.writer(Protein)
def _protein_to_genbank(obj, fh):
    _serialize_single_genbank(obj, fh)


def _construct(record, constructor=None, **kwargs):
    """Construct the object of Sequence, DNA, RNA, or Protein."""
    seq, md, imd = record
    if "lowercase" not in kwargs:
        kwargs["lowercase"] = True
    if constructor is None:
        unit = md["LOCUS"]["unit"]
        if unit == "bp":
            # RNA mol type has T instead of U for genbank from from NCBI
            constructor = DNA
        elif unit == "aa":
            constructor = Protein

    if constructor == RNA:
        return DNA(seq, metadata=md, interval_metadata=imd, **kwargs).transcribe()
    else:
        return constructor(seq, metadata=md, interval_metadata=imd, **kwargs)


def _parse_genbanks(fh):
    data_chunks = []
    for line in _line_generator(fh, skip_blanks=True, strip=False):
        if line.startswith("//"):
            yield _parse_single_genbank(data_chunks)
            data_chunks = []
        else:
            data_chunks.append(line)


def _parse_single_genbank(chunks):
    metadata = {}
    interval_metadata = None
    sequence = ""
    # each section starts with a HEADER without indent.
    section_splitter = _yield_section(lambda x: not x[0].isspace(), strip=False)
    for section in section_splitter(chunks):
        header = section[0].split(None, 1)[0]
        parser = _PARSER_TABLE.get(header, _parse_section_default)

        if header == "FEATURES":
            # This requires 'LOCUS' line parsed before 'FEATURES', which should
            # be true and is implicitly checked by the sniffer.
            parser = partial(parser, length=metadata["LOCUS"]["size"])

        parsed = parser(section)

        # reference can appear multiple times
        if header == "REFERENCE":
            if header in metadata:
                metadata[header].append(parsed)
            else:
                metadata[header] = [parsed]
        elif header == "ORIGIN":
            sequence = parsed
        elif header == "FEATURES":
            interval_metadata = parsed
        else:
            metadata[header] = parsed
    return sequence, metadata, interval_metadata


def _serialize_single_genbank(obj, fh):
    """Write a GenBank record.

    Always write it in NCBI canonical way:
    1. sequence in lowercase
    2. 'u' as 't' even in RNA molecules.

    Parameters
    ----------
    obj : Sequence or its child class
        A Sequence object or its child class representing the biological sequence to be
        serialized.
    fh : file object
        A file object open for writing.

    """
    # write out the headers
    md = obj.metadata
    for header in _HEADERS:
        serializer = _SERIALIZER_TABLE.get(header, _serialize_section_default)
        if header in md:
            out = serializer(header, md[header])
            # test if 'out' is a iterator.
            # cf. Effective Python Item 17
            if iter(out) is iter(out):
                for s in out:
                    fh.write(s)
            else:
                fh.write(out)
        if header == "FEATURES":
            if obj.has_interval_metadata():
                # magic number 21: the amount of indentation before
                # feature table starts as defined by INSDC
                indent = 21
                fh.write(
                    "{header:<{indent}}Location/Qualifiers\n".format(
                        header=header, indent=indent
                    )
                )
                for s in serializer(obj.interval_metadata._intervals, indent):
                    fh.write(s)
    # write out the sequence
    # always write RNA seq as DNA
    if isinstance(obj, RNA):
        obj = obj.reverse_transcribe()

    # always write in lowercase
    seq_str = str(obj).lower()

    for s in _serialize_origin(seq_str):
        fh.write(s)
    fh.write("//\n")


def _parse_locus(lines):
    """Parse the line LOCUS.

    Format:
    #    Positions  Contents
    #    ---------  --------
    #    00:06      LOCUS
    #    06:12      spaces
    #    12:??      Locus name
    #    ??:??      space
    #    ??:29      Length of sequence, right-justified
    #    29:33      space, bp/aa/rc, space
    #    33:41      molecule type (can be blank): DNA, ssDNA, dsRNA, tRNA, etc.
    #    41:42      space
    #    42:51      Blank (implies linear), linear or circular
    #    51:52      space
    #    52:55      The division code (e.g. BCT, VRL, INV)
    #    55:62      space
    #    62:73      Date, in the form dd-MMM-yyyy (e.g., 15-MAR-1991)
    """
    line = lines[0]
    pattern = (
        r"LOCUS"
        r" +([^\s]+)"
        r" +([0-9]+)"
        r" +(bp|aa|rc)"
        r" +(.*DNA|.*RNA)?"
        r" +(linear|circular)?"
        r" +(?!.*DNA|.*RNA)([A-Z]{3})"
        r" +([0-9]{2}-[A-Z]{3}-[0-9]{4})"
    )
    matches = re.match(pattern, line)

    try:
        res = dict(
            zip(
                ["locus_name", "size", "unit", "mol_type", "shape", "division", "date"],
                matches.groups(),
            )
        )
    except Exception:
        raise GenBankFormatError("Could not parse the LOCUS line:\n%s" % line)

    res["size"] = int(res["size"])
    return res


def _serialize_locus(header, obj, indent=12):
    """Serialize LOCUS line.

    Parameters
    ----------
    header : str
        The header of the LOCUS line.
    obj : dict
        A dictionary containing key-value pairs representing the attributes of the
        sequence entry.
    indent : int, optional
        The number of spaces used to indent the serialized LOCUS line. Defaults to 12.

    """
    # use 'or' to convert None to ''
    kwargs = {k: v or "" for k, v in obj.items()}

    return (
        "{header:<{indent}}{locus_name}   {size} {unit}"
        "   {mol_type}   {shape}   {division}   {date}\n"
    ).format(header=header, indent=indent, **kwargs)


def _parse_reference(lines):
    """Parse single REFERENCE field."""
    res = {}
    # magic number 11: the non keyworded lines in REFERENCE
    # are at least indented with 11 spaces.
    feature_indent = " " * 11
    section_splitter = _yield_section(
        lambda x: not x.startswith(feature_indent), skip_blanks=True, strip=False
    )
    for section in section_splitter(lines):
        label, data = _parse_section_default(
            section, join_delimiter=" ", return_label=True
        )
        res[label] = data
    return res


def _serialize_reference(header, obj, indent=12):
    """Serialize REFERENCE.

    Parameters
    ----------
    header : str
        The header of the REFERENCE section.
    obj : list
        A list of dictionaries representing the individual references cited in the
        sequence entry.
    indent : int, optional
        The number of spaces used to indent each line of the serialized REFERENCE
        section. Defaults to 12.

    """
    padding = "  "
    sort_order = {"REFERENCE": 0, "AUTHORS": 1, "TITLE": 2, "JOURNAL": 3, "PUBMED": 4}
    for obj_i in obj:
        ref_i = []
        for h in sorted(obj_i, key=lambda k: sort_order.get(k, 100)):
            if h == header:
                s = "{h:<{indent}}{ref}".format(h=h, indent=indent, ref=obj_i[h])
            else:
                s = "{h:<{indent}}{value}".format(
                    h=padding + h, indent=indent, value=obj_i[h]
                )
            ref_i.append(s)
        yield "%s\n" % "\n".join(ref_i)


def _parse_source(lines):
    """Parse SOURCE field."""
    res = {}
    # magic number 11: the non keyworded lines in SOURCE
    # are at least indented with 11 spaces.
    feature_indent = " " * 11
    section_splitter = _yield_section(
        lambda x: not x.startswith(feature_indent), skip_blanks=True, strip=False
    )
    # SOURCE line is not informative; skip it
    _, organism = list(section_splitter(lines))

    res["ORGANISM"] = organism[0].split(None, 1)[1].strip()
    res["taxonomy"] = " ".join([i.strip() for i in organism[1:]])
    return res


def _serialize_source(header, obj, indent=12):
    """Serialize SOURCE.

    Parameters
    ----------
    header : str
        The header of the SOURCE section,
    obj : dict
        A dictionary containing key-value pairs representing the attributes of the
        SOURCE section.
    indent : int, optional
        The number of spaces used to indent each line of the serialized SOURCE section.
        Defaults to 12.

    """
    s = (
        "{header:<{indent}}{organism}\n"
        "{h:<{indent}}{organism}\n"
        "{space}{taxonomy}\n"
    ).format(
        header=header,
        indent=indent,
        h="  ORGANISM",
        organism=obj["ORGANISM"],
        space=" " * 12,
        taxonomy=obj["taxonomy"],
    )
    return s


def _parse_origin(lines):
    """Parse the ORIGIN section for sequence."""
    sequence = []
    for line in lines:
        if line.startswith("ORIGIN"):
            continue
        # remove the number at the beg of each line
        items = line.split()
        sequence.append("".join(items[1:]))
    return "".join(sequence)


def _serialize_origin(seq, indent=9):
    """Serialize seq to ORIGIN.

    Parameters
    ----------
    seq : str
        The sequence data to be serialized.
    indent : int, optional
        The number of spaces used to indent each line of the serialized ORIGIN section.

    """
    n = 1
    line_size = 60
    frag_size = 10
    for i in range(0, len(seq), line_size):
        line = seq[i : i + line_size]
        s = "{n:>{indent}} {s}\n".format(
            n=n, indent=indent, s=chunk_str(line, frag_size, " ")
        )
        if n == 1:
            s = "ORIGIN\n" + s
        n += line_size
        yield s


_PARSER_TABLE = {
    "LOCUS": _parse_locus,
    "SOURCE": _parse_source,
    "REFERENCE": _parse_reference,
    "FEATURES": _parse_feature_table,
    "ORIGIN": _parse_origin,
}


_SERIALIZER_TABLE = {
    "LOCUS": _serialize_locus,
    "SOURCE": _serialize_source,
    "REFERENCE": _serialize_reference,
    "FEATURES": _serialize_feature_table,
}
