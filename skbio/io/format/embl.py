"""EMBL format (:mod:`skbio.io.format.embl`)
=========================================

.. currentmodule:: skbio.io.format.embl

EMBL format stores sequence and its annotation together. The start of the
annotation section is marked by a line beginning with the word "ID". The start
of sequence section is marked by a line beginning with the word "SQ".
The "//" (terminator) line also contains no data or comments and designates
the end of an entry. More information on EMBL file format can be found
here [1]_.

The EMBL file may end with .embl or .txt extension. An example of EMBL file
can be seen here [2]_.

Feature Level Products
^^^^^^^^^^^^^^^^^^^^^^

As described in [3]_ *"Feature-level products contain nucleotide sequence
and related annotations derived from submitted ENA assembled and annotated
sequences. Data are distributed in flatfile format, similar to that of parent
ENA records, with each flatfile representing a single feature"*.
While only the sequence of the feature is included in such entries, features
are derived from the parent entry, and can't be applied as interval metadata.
For such reason, interval metatdata are ignored from Feature-level products,
as they will be ignored by subsetting a generic Sequence object.

Format Support
--------------
**Has Sniffer: Yes**

**NOTE: No protein support at the moment**

Current protein support development is tracked in issue-1499 [4]_

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.RNA`                                      |
+------+------+---------------------------------------------------------------+
|No    |No    |:mod:`skbio.sequence.Protein`                                  |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   | generator of :mod:`skbio.sequence.Sequence` objects           |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------

Sections before ``FH (Feature Header)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All the sections before ``FH (Feature Header)`` will be read into the attribute
of ``metadata``. The header and its content of a section are stored as
key-value pairs in ``metadata``. For the ``RN (Reference Number)``
section, its value is stored as a list, as there are often multiple
reference sections in one EMBL record.

``FT`` section
^^^^^^^^^^^^^^

See :ref:`Genbank FEATURES section<genbank_feature_section>`

``SQ`` section
^^^^^^^^^^^^^^
The sequence in the ``SQ`` section is always in lowercase for
the EMBL files downloaded from ENA. For the RNA molecules, ``t``
(thymine), instead of ``u`` (uracil) is used in the sequence. All
EMBL writers follow these conventions while writing EMBL files.

Examples
--------
Reading EMBL Files
^^^^^^^^^^^^^^^^^^
Suppose we have the following EMBL file example:

>>> embl_str = '''
... ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
... XX
... AC   X56734; S46826;
... XX
... DT   12-SEP-1991 (Rel. 29, Created)
... DT   25-NOV-2005 (Rel. 85, Last updated, Version 11)
... XX
... DE   Trifolium repens mRNA for non-cyanogenic beta-glucosidase
... XX
... KW   beta-glucosidase.
... XX
... OS   Trifolium repens (white clover)
... OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
... OC   Spermatophyta; Magnoliophyta; eudicotyledons; Gunneridae;
... OC   Pentapetalae; rosids; fabids; Fabales; Fabaceae; Papilionoideae;
... OC   Trifolieae; Trifolium.
... XX
... RN   [5]
... RP   1-1859
... RX   DOI; 10.1007/BF00039495.
... RX   PUBMED; 1907511.
... RA   Oxtoby E., Dunn M.A., Pancoro A., Hughes M.A.;
... RT   "Nucleotide and derived amino acid sequence of the cyanogenic
... RT   beta-glucosidase (linamarase) from white clover
... RT   (Trifolium repens L.)";
... RL   Plant Mol. Biol. 17(2):209-219(1991).
... XX
... RN   [6]
... RP   1-1859
... RA   Hughes M.A.;
... RT   ;
... RL   Submitted (19-NOV-1990) to the INSDC.
... RL   Hughes M.A., University of Newcastle Upon Tyne, Medical School,
... RL   Newcastle
... RL   Upon Tyne, NE2 4HH, UK
... XX
... DR   MD5; 1e51ca3a5450c43524b9185c236cc5cc.
... XX
... FH   Key             Location/Qualifiers
... FH
... FT   source          1..1859
... FT                   /organism="Trifolium repens"
... FT                   /mol_type="mRNA"
... FT                   /clone_lib="lambda gt10"
... FT                   /clone="TRE361"
... FT                   /tissue_type="leaves"
... FT                   /db_xref="taxon:3899"
... FT   mRNA            1..1859
... FT                   /experiment="experimental evidence, no additional
... FT                   details recorded"
... FT   CDS             14..1495
... FT                   /product="beta-glucosidase"
... FT                   /EC_number="3.2.1.21"
... FT                   /note="non-cyanogenic"
... FT                   /db_xref="GOA:P26204"
... FT                   /db_xref="InterPro:IPR001360"
... FT                   /db_xref="InterPro:IPR013781"
... FT                   /db_xref="InterPro:IPR017853"
... FT                   /db_xref="InterPro:IPR033132"
... FT                   /db_xref="UniProtKB/Swiss-Prot:P26204"
... FT                   /protein_id="CAA40058.1"
... FT                   /translation="MDFIVAIFALFVISSFTITSTNAVEASTLLDIGNLSRS
... FT                   SFPRGFIFGAGSSAYQFEGAVNEGGRGPSIWDTFTHKYPEKIRDGSNADITV
... FT                   DQYHRYKEDVGIMKDQNMDSYRFSISWPRILPKGKLSGGINHEGIKYYNNLI
... FT                   NELLANGIQPFVTLFHWDLPQVLEDEYGGFLNSGVINDFRDYTDLCFKEFGD
... FT                   RVRYWSTLNEPWVFSNSGYALGTNAPGRCSASNVAKPGDSGTGPYIVTHNQI
... FT                   LAHAEAVHVYKTKYQAYQKGKIGITLVSNWLMPLDDNSIPDIKAAERSLDFQ
... FT                   FGLFMEQLTTGDYSKSMRRIVKNRLPKFSKFESSLVNGSFDFIGINYYSSSY
... FT                   ISNAPSHGNAKPSYSTNPMTNISFEKHGIPLGPRAASIWIYVYPYMFIQEDF
... FT                   EIFCYILKINITILQFSITENGMNEFNDATLPVEEALLNTYRIDYYYRHLYY
... FT                   IRSAIRAGSNVKGFYAWSFLDCNEWFAGFTVRFGLNFVD"
... XX
... SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
...      aaacaaacca aatatggatt ttattgtagc catatttgct ctgtttgtta ttagctcatt
...      cacaattact tccacaaatg cagttgaagc ttctactctt cttgacatag gtaacctgag
...      tcggagcagt tttcctcgtg gcttcatctt tggtgctgga tcttcagcat accaatttga
...      aggtgcagta aacgaaggcg gtagaggacc aagtatttgg gataccttca cccataaata
...      tccagaaaaa ataagggatg gaagcaatgc agacatcacg gttgaccaat atcaccgcta
...      caaggaagat gttgggatta tgaaggatca aaatatggat tcgtatagat tctcaatctc
...      ttggccaaga atactcccaa agggaaagtt gagcggaggc ataaatcacg aaggaatcaa
...      atattacaac aaccttatca acgaactatt ggctaacggt atacaaccat ttgtaactct
...      ttttcattgg gatcttcccc aagtcttaga agatgagtat ggtggtttct taaactccgg
...      tgtaataaat gattttcgag actatacgga tctttgcttc aaggaatttg gagatagagt
...      gaggtattgg agtactctaa atgagccatg ggtgtttagc aattctggat atgcactagg
...      aacaaatgca ccaggtcgat gttcggcctc caacgtggcc aagcctggtg attctggaac
...      aggaccttat atagttacac acaatcaaat tcttgctcat gcagaagctg tacatgtgta
...      taagactaaa taccaggcat atcaaaaggg aaagataggc ataacgttgg tatctaactg
...      gttaatgcca cttgatgata atagcatacc agatataaag gctgccgaga gatcacttga
...      cttccaattt ggattgttta tggaacaatt aacaacagga gattattcta agagcatgcg
...      gcgtatagtt aaaaaccgat tacctaagtt ctcaaaattc gaatcaagcc tagtgaatgg
...      ttcatttgat tttattggta taaactatta ctcttctagt tatattagca atgccccttc
...      acatggcaat gccaaaccca gttactcaac aaatcctatg accaatattt catttgaaaa
...      acatgggata cccttaggtc caagggctgc ttcaatttgg atatatgttt atccatatat
...      gtttatccaa gaggacttcg agatcttttg ttacatatta aaaataaata taacaatcct
...      gcaattttca atcactgaaa atggtatgaa tgaattcaac gatgcaacac ttccagtaga
...      agaagctctt ttgaatactt acagaattga ttactattac cgtcacttat actacattcg
...      ttctgcaatc agggctggct caaatgtgaa gggtttttac gcatggtcat ttttggactg
...      taatgaatgg tttgcaggct ttactgttcg ttttggatta aactttgtag attagaaaga
...      tggattaaaa aggtacccta agctttctgc ccaatggtac aagaactttc tcaaaagaaa
...      ctagctagta ttattaaaag aactttgtag tagattacag tacatcgttt gaagttgagt
...      tggtgcacct aattaaataa aagaggttac tcttaacata tttttaggcc attcgttgtg
...      aagttgttag gctgttattt ctattatact atgttgtagt aataagtgca ttgttgtacc
...      agaagctatg atcataacta taggttgatc cttcatgtat cagtttgatg ttgagaatac
...      tttgaattaa aagtcttttt ttattttttt aaaaaaaaaa aaaaaaaaaa aaaaaaaaa
... //
... '''


Now we can read it as ``DNA`` object:

>>> import io
>>> from skbio import DNA, RNA, Sequence
>>> embl = io.StringIO(embl_str)
>>> dna_seq = DNA.read(embl)
>>> dna_seq
DNA
----------------------------------------------------------------------
Metadata:
    'ACCESSION': 'X56734; S46826;'
    'CROSS_REFERENCE': <class 'list'>
    'DATE': <class 'list'>
    'DBSOURCE': 'MD5; 1e51ca3a5450c43524b9185c236cc5cc.'
    'DEFINITION': 'Trifolium repens mRNA for non-cyanogenic beta-
                   glucosidase'
    'KEYWORDS': 'beta-glucosidase.'
    'LOCUS': <class 'dict'>
    'REFERENCE': <class 'list'>
    'SOURCE': <class 'dict'>
    'VERSION': 'X56734.1'
Interval metadata:
    3 interval features
Stats:
    length: 1859
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 35.99%
----------------------------------------------------------------------
0    AAACAAACCA AATATGGATT TTATTGTAGC CATATTTGCT CTGTTTGTTA TTAGCTCATT
60   CACAATTACT TCCACAAATG CAGTTGAAGC TTCTACTCTT CTTGACATAG GTAACCTGAG
...
1740 AGAAGCTATG ATCATAACTA TAGGTTGATC CTTCATGTAT CAGTTTGATG TTGAGAATAC
1800 TTTGAATTAA AAGTCTTTTT TTATTTTTTT AAAAAAAAAA AAAAAAAAAA AAAAAAAAA

Since this is a mRNA molecule, we may want to read it as ``RNA``.
As the EMBL file usually have ``t`` instead of ``u`` in
the sequence, we can read it as ``RNA`` by converting ``t`` to ``u``:

>>> embl = io.StringIO(embl_str)
>>> rna_seq = RNA.read(embl)
>>> rna_seq
RNA
----------------------------------------------------------------------
Metadata:
    'ACCESSION': 'X56734; S46826;'
    'CROSS_REFERENCE': <class 'list'>
    'DATE': <class 'list'>
    'DBSOURCE': 'MD5; 1e51ca3a5450c43524b9185c236cc5cc.'
    'DEFINITION': 'Trifolium repens mRNA for non-cyanogenic beta-
                   glucosidase'
    'KEYWORDS': 'beta-glucosidase.'
    'LOCUS': <class 'dict'>
    'REFERENCE': <class 'list'>
    'SOURCE': <class 'dict'>
    'VERSION': 'X56734.1'
Interval metadata:
    3 interval features
Stats:
    length: 1859
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 35.99%
----------------------------------------------------------------------
0    AAACAAACCA AAUAUGGAUU UUAUUGUAGC CAUAUUUGCU CUGUUUGUUA UUAGCUCAUU
60   CACAAUUACU UCCACAAAUG CAGUUGAAGC UUCUACUCUU CUUGACAUAG GUAACCUGAG
...
1740 AGAAGCUAUG AUCAUAACUA UAGGUUGAUC CUUCAUGUAU CAGUUUGAUG UUGAGAAUAC
1800 UUUGAAUUAA AAGUCUUUUU UUAUUUUUUU AAAAAAAAAA AAAAAAAAAA AAAAAAAAA

We can also ``trascribe`` a sequence and verify that it will be a ``RNA``
sequence

>>> rna_seq == dna_seq.transcribe()
True

Reading EMBL Files using generators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Soppose we have an EMBL file with multiple records: we can instantiate a
generator object to deal with multiple records

>>> import skbio
>>> embl = io.StringIO(embl_str)
>>> embl_gen = skbio.io.read(embl, format="embl")
>>> dna_seq = next(embl_gen)

For more informations, see :mod:`skbio.io`

References
----------
.. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/release/doc/usrman.txt
.. [2] http://www.ebi.ac.uk/ena/data/view/X56734&display=text
.. [3] http://www.ebi.ac.uk/ena/browse/feature-level-products
.. [4] https://github.com/scikit-bio/scikit-bio/issues/1499


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# std modules
import re
import copy
import textwrap

from functools import partial

# skbio modules
from skbio.io import create_format, EMBLFormatError
from skbio.io.format._base import _line_generator, _get_nth_sequence
from skbio.io.format._sequence_feature_vocabulary import (
    _yield_section,
    _parse_single_feature,
    _serialize_section_default,
    _serialize_single_feature,
)
from skbio.metadata import IntervalMetadata
from skbio.sequence import Sequence, DNA, RNA, Protein
from skbio.util._misc import chunk_str


# look at skbio.io.registry to have an idea on how to define this class
embl = create_format("embl")

# This list is ordered used to read and write embl file. By processing those
# values one by one, I will write embl sections with the same order
_HEADERS = [
    "LOCUS",
    "ACCESSION",
    "PARENT_ACCESSION",
    "PROJECT_IDENTIFIER",
    "DATE",
    "DEFINITION",
    "GENE_NAME",
    "KEYWORDS",
    "SOURCE",
    "REFERENCE",
    "DBSOURCE",
    "COMMENT",
    "FEATURES",
]


# embl has a series of keys different from genbank; moreover keys are not so
# easy to understand (eg. RA for AUTHORS). I want to use the same keys used by
# genbank both to convert between formats and to use the same methods to get
# info from Sequence and its derived objects Here is a dictionary of keys
# conversion (EMBL->GB). All the unspecified keys will remain in embl format
KEYS_TRANSLATOR = {
    # identification
    "ID": "LOCUS",
    "AC": "ACCESSION",
    # PA means PARENT ACCESSION (?) and applies to
    # feature-level-products entries
    "PA": "PARENT_ACCESSION",
    "PR": "PROJECT_IDENTIFIER",
    "DT": "DATE",
    "DE": "DEFINITION",
    "GN": "GENE_NAME",  # uniprot specific
    "KW": "KEYWORDS",
    # Source (taxonomy and classification)
    "OS": "ORGANISM",
    "OC": "taxonomy",
    "OG": "organelle",
    # reference keys
    "RA": "AUTHORS",
    "RP": "REFERENCE",
    "RC": "REFERENCE_COMMENT",
    "RX": "CROSS_REFERENCE",
    "RG": "GROUP",
    "RT": "TITLE",
    "RL": "JOURNAL",
    # Cross references
    "DR": "DBSOURCE",
    "CC": "COMMENT",
    # features
    "FH": "FEATURES",
    "FT": "FEATURES",
    "SQ": "ORIGIN",
}


# the inverse of KEYS_TRANSLATOR, for semplicity
REV_KEYS_TRANSLATOR = {v: k for k, v in KEYS_TRANSLATOR.items()}


# the original genbank _yield_section divides entries in sections relying on
# spaces (the same section has the same level of indentation). EMBL entries
# have a key for each line, so to divide record in sections I need to define a
# correspondance for each key to section, then I will divide a record in
# sections using these section name.
KEYS_2_SECTIONS = {
    # identification
    "ID": "LOCUS",
    "AC": "ACCESSION",
    # PA means PARENT ACCESSION (?) and applies to
    # feature-level-products entries
    "PA": "PARENT_ACCESSION",
    "PR": "PROJECT_IDENTIFIER",
    "DT": "DATE",
    "DE": "DEFINITION",
    "GN": "GENE_NAME",  # uniprot specific
    "KW": "KEYWORDS",
    # Source (taxonomy and classification)
    "OS": "SOURCE",
    "OC": "SOURCE",
    "OG": "SOURCE",
    # reference keys
    "RA": "REFERENCE",
    "RP": "REFERENCE",
    "RC": "REFERENCE",
    "RX": "REFERENCE",
    "RG": "REFERENCE",
    "RT": "REFERENCE",
    "RL": "REFERENCE",
    # This shuold be Reference Number. However, to split
    # between references with _embl_yield_section I need to
    # change section after reading one reference. So a single
    # reference is completed when I found a new RN. The
    # reference number information will be the reference
    # position in the final REFERENCE list metadata
    "RN": "SPACER",
    # Cross references
    "DR": "DBSOURCE",
    "CC": "COMMENT",
    "AH": "ASSEMBLY",
    "AS": "ASSEMBLY",
    "FH": "FEATURES",
    "FT": "FEATURES",
    # sequence
    "SQ": "ORIGIN",
    "  ": "ORIGIN",
    "CO": "CONSTRUCTED",
    # spacer (discarded)
    "XX": "SPACER",
}


# for convenience: I think such functions are more readadble while accessing
# values in lambda functions
def _get_embl_key(line):
    """Return first part of a string as a embl key (ie 'AC M14399;' -> 'AC')."""
    # embl keys have a fixed size of 2 chars
    return line[:2]


def _get_embl_section(line):
    """Return the embl section from uniprot key(ie 'RA' -> 'REFERENCE')."""
    # get embl key
    key = _get_embl_key(line)

    # get embl section from key
    section = KEYS_2_SECTIONS[key]

    return section


def _translate_key(key):
    """Translate a single key from EMBL to genbank.

    Returns key itself if no traslation is defined.
    """
    return KEYS_TRANSLATOR.get(key, key)


# a method to translate keys from embl to genbank for a dict object. All keys
# not defined in the original dict will remain the same
def _translate_keys(data):
    """Translate keys from EMBL to genbank for a dict object.

    Translate a dictionary of uniprot key->value in a genbank like dictionary
    of key values. Keep old keys if no translation is defined.
    """
    # traslate keys and get a new_data object
    new_data = {_translate_key(k): v for k, v in data.items()}

    return new_data


# define a default textwrap.Wrapper for embl
def _get_embl_wrapper(embl_key, indent=5, subsequent_indent=None, width=80):
    """Return a textwrap.TextWrapper for embl records.

    For example, write <key>   <string> by providing embl key and a string. Wrap text to
    80 column.
    """
    # define the string to prepen (eg "OC   ")
    prepend = "{key:<{indent}}".format(key=embl_key, indent=indent)

    # deal with 2° strings and more
    if subsequent_indent is None:
        subsequent_prepend = prepend

    else:
        subsequent_prepend = "{key:<{indent}}".format(
            key=embl_key, indent=subsequent_indent
        )

    # define a text wrapper object
    wrapper = textwrap.TextWrapper(
        initial_indent=prepend, subsequent_indent=subsequent_prepend, width=width
    )

    return wrapper


def _serialize_list(embl_wrapper, data, sep="\n"):
    """Serialize a list of obj using a textwrap.TextWrapper instance.

    Returns one string of wrapped embl objects.
    """
    # the output array
    output = []

    for line in data:
        output += embl_wrapper.wrap(line)

    # merge dates in one string. Add final newline
    output = sep.join(output) + "\n"

    # return comupted string
    return output


# Method to determine if file is in EMBL format or not. A uniprot embl format
# can't be parsed by this module (at the moment)
@embl.sniffer()
def _embl_sniffer(fh):
    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    try:
        _parse_id([line])

    except EMBLFormatError:
        return False, {}

    return True, {}


@embl.reader(None)
def _embl_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_embls(fh):
        yield _construct(record, constructor, **kwargs)


# Method to read EMBL data as skbio.sequence.DNA
@embl.reader(Sequence)
def _embl_to_sequence(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_embls(fh), seq_num)
    return _construct(record, Sequence, **kwargs)


# Method to read EMBL data as skbio.sequence.DNA
@embl.reader(DNA)
def _embl_to_dna(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_embls(fh), seq_num)
    return _construct(record, DNA, **kwargs)


# Method to read EMBL data as skbio.sequence.DNA
@embl.reader(RNA)
def _embl_to_rna(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_embls(fh), seq_num)
    return _construct(record, RNA, **kwargs)


# No protein support at the moment
@embl.reader(Protein)
def _embl_to_protein(fh, seq_num=1, **kwargs):
    # no protein support, at the moment
    raise EMBLFormatError(
        "There's no protein support for EMBL record. "
        "Current status of EMBL protein support is "
        "described in issue-1499 (https://github.com/"
        "biocore/scikit-bio/issues/1499)"
    )


# Writer methods
@embl.writer(None)
def _generator_to_embl(obj, fh):
    for obj_i in obj:
        _serialize_single_embl(obj_i, fh)


@embl.writer(Sequence)
def _sequence_to_embl(obj, fh):
    _serialize_single_embl(obj, fh)


@embl.writer(DNA)
def _dna_to_embl(obj, fh):
    _serialize_single_embl(obj, fh)


@embl.writer(RNA)
def _rna_to_embl(obj, fh):
    _serialize_single_embl(obj, fh)


@embl.writer(Protein)
def _protein_to_embl(obj, fh):
    # no protein support, at the moment
    raise EMBLFormatError(
        "There's no protein support for EMBL record. "
        "Current status of EMBL protein support is "
        "described in issue-1499 (https://github.com/"
        "biocore/scikit-bio/issues/1499)"
    )


def _construct(record, constructor=None, **kwargs):
    """Construct the object of Sequence, DNA, RNA, or Protein."""
    # sequence, metadata and interval metadata
    seq, md, imd = record

    if "lowercase" not in kwargs:
        kwargs["lowercase"] = True

    if constructor is None:
        unit = md["LOCUS"]["unit"]
        if unit == "bp":
            # RNA mol type has T instead of U for genbank from from NCBI
            constructor = DNA

        elif unit == "aa":
            # no protein support, at the moment
            # constructor = Protein
            raise EMBLFormatError("There's no protein support for EMBL record")

    if constructor == RNA:
        return DNA(seq, metadata=md, interval_metadata=imd, **kwargs).transcribe()
    else:
        return constructor(seq, metadata=md, interval_metadata=imd, **kwargs)


# looks like the genbank _parse_genbank
def _parse_embls(fh):
    """Chunk multiple EMBL records by '//', and returns a generator."""
    data_chunks = []
    for line in _line_generator(fh, skip_blanks=True, strip=False):
        if line.startswith("//"):
            yield _parse_single_embl(data_chunks)
            data_chunks = []
        else:
            data_chunks.append(line)


def _parse_single_embl(chunks):
    metadata = {}
    interval_metadata = None
    sequence = ""

    # define a section splitter with _embl_yield_section function defined in
    # this module (return the embl section by embl key). returns generator for
    # each block with different line type
    section_splitter = _embl_yield_section(
        lambda line: _get_embl_section(line), skip_blanks=True, strip=False
    )

    # process each section, like genbank does.
    for section, section_name in section_splitter(chunks):
        # section is a list of records with the same session (eg RA, RP for
        # for a single reference). section_name is the name of the section
        # (eg REFERENCE for the section of the previous example)

        # search for a specific method in PARSER_TABLE using section_name or
        # set _embl_parse_section_default
        parser = _PARSER_TABLE.get(section_name, _embl_parse_section_default)

        if section_name == "FEATURES":
            # This requires 'ID' line parsed before 'FEATURES', which should
            # be true and is implicitly checked by the sniffer. This is true
            # since the first section is parsed by the last else condition

            if "PARENT_ACCESSION" in metadata:
                # this is a feature-level-products entry and features are
                # relative to parent accession; in the same way a subset of a
                # Sequence object has no interval metadata, I will refuse to
                # process interval metadata here
                continue

            # partials add arguments to previous defined functions, in this
            # case length of Sequence object
            parser = partial(parser, length=metadata["LOCUS"]["size"])

        elif section_name == "COMMENT":
            # mantain newlines in comments
            # partials add arguments to previous defined functions
            parser = partial(parser, join_delimiter="\n")

        # call function on section
        parsed = parser(section)

        # reference can appear multiple times
        if section_name == "REFERENCE":
            # genbank data hasn't CROSS_REFERENCE section, To have a similar
            # metatadata object, I chose to remove CROSS_REFERENCE from
            # each single reference and put them in metadata. Since I could
            # have more references, I need a list of CROSS_REFERENCE, with
            # None values when CROSS_REFERENCE are not defined: there are cases
            # in which some references have a CROSS_REFERENCE and others not.
            # So each reference will have it's cross reference in the same
            # index position, defined or not
            cross_reference = parsed.pop("CROSS_REFERENCE", None)

            # fix REFERENCE metadata. Ask if is the first reference or not
            # I need a reference number as genbank, this could be reference
            # size
            if section_name in metadata:
                RN = len(metadata[section_name]) + 1

            else:
                RN = 1

            # fix reference fields. Get RN->REFERENCE value from dict
            positions = parsed.pop("REFERENCE", None)
            parsed["REFERENCE"] = str(RN)

            # append position to RN (eg "1  (bases 1 to 63)")
            if positions:
                parsed["REFERENCE"] += "  %s" % (positions)

            # cross_reference will be a list of cross reference; Also
            # metadata[REFERENCE] is a list of references
            if section_name in metadata:
                # I've already seen a reference, append new one
                metadata[section_name].append(parsed)
                metadata["CROSS_REFERENCE"].append(cross_reference)

            else:
                # define a list for this first reference and its RX
                metadata[section_name] = [parsed]
                metadata["CROSS_REFERENCE"] = [cross_reference]

        elif section_name == "ORIGIN":
            sequence = parsed

        elif section_name == "FEATURES":
            interval_metadata = parsed

        elif section_name == "DATE":
            # read data (list)
            metadata[section_name] = parsed

            # fix locus metadata using last date. Take only last date
            date = metadata[section_name][-1].split()[0]
            metadata["LOCUS"]["date"] = date

        # parse all the others sections (SOURCE, ...)
        else:
            metadata[section_name] = parsed

    # after metadata were read, add a VERSION section like genbank
    # eval if entry is a feature level product or not
    if "ACCESSION" in metadata:
        metadata["VERSION"] = "{accession}.{version}".format(
            accession=metadata["ACCESSION"].split(";")[0],
            version=metadata["LOCUS"]["version"],
        )

    elif "PARENT_ACCESSION" in metadata:
        # locus name is in the format
        # <accession>.<version>:<feature location>:<feature name>[:ordinal]
        # and ordinal could be present or not, depends on how many features
        # are found in such location. Such entry couldn't be found in others
        # database like NCBI (at the moment) so we will take the version
        # relying on parent accession (hoping that an update in the parent
        # accession will generate an update in all feature level products)
        metadata["VERSION"] = metadata["PARENT_ACCESSION"]

    # return a string, metatdata as a dictionary and IntervalMetadata object
    return sequence, metadata, interval_metadata


def _write_serializer(fh, serializer, embl_key, data):
    """Write serializer to a file. Append 'XX'."""
    # call the serializer function
    out = serializer(embl_key, data)
    # test if 'out' is a iterator.
    # cf. Effective Python Item 17
    if iter(out) is iter(out):
        for s in out:
            fh.write(s)

    else:
        fh.write(out)

    # add spacer between sections
    fh.write("XX\n")


# main function for writer methods
def _serialize_single_embl(obj, fh):
    """Write a EMBL record.

    Always write it in ENA canonical way:
    1. sequence in lowercase (uniprot are uppercase)
    2. 'u' as 't' even in RNA molecules.

    Parameters
    ----------
    obj : Sequence or its child class
        A Sequence object or its child class representing the biological sequence to be
        serialized.
    fh : file object
        A file object open for writing.

    """
    # shortcut to deal with metadata
    md = obj.metadata

    # embl has a different magick number than embl
    serialize_default = partial(_serialize_section_default, indent=5)

    # Now cicle for GB like headers (sections) in _HEADERS.
    for header in _HEADERS:
        # Get appropriate serializer method or default one
        serializer = _SERIALIZER_TABLE.get(header, serialize_default)

        # headers needs to be converted into embl, or matained as they are
        # if no conversion could be defined.
        embl_key = REV_KEYS_TRANSLATOR.get(header, header)

        # this is true also for locus line
        if header in md:
            # deal with special source case, add cross references if needed
            if header == "REFERENCE":
                serializer = partial(
                    serializer, cross_references=md.get("CROSS_REFERENCE")
                )

            elif header == "LOCUS":
                # pass also metadata (in case of entries from genbank)
                serializer = partial(serializer, metadata=md)

            # call the serializer function
            _write_serializer(fh, serializer, embl_key, md[header])

        else:
            # header not in metadata. Could be date read from GB?
            if header == "DATE":
                # Have I date in locus metadata?
                if md["LOCUS"]["date"]:
                    # call serializer on date. Date is a list of values
                    _write_serializer(fh, serializer, embl_key, [md["LOCUS"]["date"]])

        if header == "FEATURES":
            if obj.has_interval_metadata():
                # magic number 21: the amount of indentation before
                # feature table starts as defined by INSDC
                indent = 21
                feature_key = "FH   Key"
                fh.write(
                    "{header:<{indent}}Location/Qualifiers\n".format(
                        header=feature_key, indent=indent
                    )
                )

                # add FH spacer
                fh.write("FH\n")

                for s in serializer(obj.interval_metadata._intervals, indent):
                    fh.write(s)

                # add spacer between sections
                fh.write("XX\n")

    # write out the sequence
    # always write RNA seq as DNA
    if isinstance(obj, RNA):
        obj = obj.reverse_transcribe()

    # serialize sequence from a Sequence object
    for s in _serialize_sequence(obj):
        fh.write(s)

    # terminate a embl record with
    fh.write("//\n")


def _parse_id(lines):
    """Parse the identification line of an EMBL record.

    From EMBL user manual (Release 130, November 2016).
    (ftp://ftp.ebi.ac.uk/pub/databases/embl/release/doc/usrman.txt)

    The ID (IDentification) line is always the first line of an entry. The
    format of the ID line is:
    ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
    The tokens represent:
       1. Primary accession number
       2. Sequence version number
       3. Topology: 'circular' or 'linear'
       4. Molecule type (see note 1 below)
       5. Data class (see section 3.1 of EMBL user manual)
       6. Taxonomic division (see section 3.2 of EMBL user manual)
       7. Sequence length (see note 2 below)

    Note 1 - Molecule type: this represents the type of molecule as stored and
    can be any value from the list of current values for the mandatory mol_type
    source qualifier. This item should be the same as the value in the mol_type
    qualifier(s) in a given entry.
    Note 2 - Sequence length: The last item on the ID line is the length of the
    sequence (the total number of bases in the sequence). This number includes
    base positions reported as present but undetermined (coded as "N").
    An example of a complete identification line is shown below:
    ID   CD789012; SV 4; linear; genomic DNA; HTG; MAM; 500 BP.
    """
    # get only the first line of EMBL record
    line = lines[0]

    # define a specific patter for EMBL
    pattern = re.compile(
        r"ID"
        r" +([^\s]+);"  # ie: CD789012
        r" +SV ([0-9]*);"  # 4
        r" +(\w+);"  # linear
        r" +([^;]+);"  # genomic DNA
        r" +(\w*);"  # HTG
        r" +(\w+);"  # MAM
        r" +(\d+)"  # 500
        r" +(\w+)\.$"
    )  # BP

    # search it
    matches = re.match(pattern, line)

    try:
        res = dict(
            zip(
                [
                    "locus_name",
                    "version",
                    "shape",
                    "mol_type",
                    "class",
                    "division",
                    "size",
                    "unit",
                ],
                matches.groups(),
            )
        )
    except AttributeError:
        raise EMBLFormatError("Could not parse the ID line:\n%s" % line)

    # check for CON entries:
    if res["class"] == "CON":
        # entries like http://www.ebi.ac.uk/ena/data/view/LT357133
        # doesn't have sequence, so can't be read by skbio.sequence
        raise EMBLFormatError(
            "There's no support for embl CON record: for more information "
            "see issue-1506 (https://github.com/scikit-bio/scikit-bio/issues/"
            "1506)"
        )

    # those values are integer
    res["size"] = int(res["size"])

    # version could be integer
    if res["version"]:
        res["version"] = int(res["version"])

    # unit are in lower cases in others modules
    res["unit"] = res["unit"].lower()

    # initialize a date record (for gb compatibility)
    res["date"] = None

    # returning parsed attributes
    return res


def _serialize_id(header, obj, metadata={}, indent=5):
    """Serialize ID line.

    Parameters
    ----------
    header : str
        The header of the ID line. Usually 'ID' for EMBL or 'LOCUS' for GenBank.
    obj : dict
        A dictionary containing key-value pairs representing the attributes of the
        sequence entry.
    metadata : dict, optional
        Additional metadata information, typically extracted from a GenBank entry.
    indent : int, optional
        The number of spaces used to indent the serialized ID line. Defaults to 5.

    """
    # get key->value pairs, or key->'' if values is None
    kwargs = {k: "" if v is None else v for k, v in obj.items()}

    # then unit is in upper cases
    kwargs["unit"] = kwargs["unit"].upper()

    # check for missing keys (eg from gb data). Keys in md are in uppercase
    for key in ["version", "class"]:
        if key not in kwargs:
            if key.upper() in metadata:
                kwargs[key] = metadata[key.upper()]

            else:
                kwargs[key] = ""

    # version from genbank could be "M14399.1  GI:145229". I need an integer
    version = kwargs["version"]

    # version could by empty, integer or text
    if version != "":
        try:
            int(kwargs["version"])

        # could be a text like M14399.1
        except ValueError:
            match = re.search(r"^\w+\.([0-9]+)", version)

            if match:
                kwargs["version"] = match.groups()[0]

    # return first line
    return (
        "{header:<{indent}}{locus_name}; SV {version}; {shape}; "
        "{mol_type}; {class}; {division}; {size} {unit}.\n"
    ).format(header=header, indent=indent, **kwargs)


# similar to skbio.io.format._sequence_feature_vocabulary.__yield_section
# but applies to embl file format
def _embl_yield_section(get_line_key, **kwargs):
    """Return function that returns successive sections from file.

    Parameters
    ----------
    get_line_key : callable
        It takes a string as input and a key indicating the section
        (could be the embl key or embl KEYS_2_SECTIONS)

    kwargs : dict, optional
        Keyword arguments will be passed to `_line_generator`.

    Returns
    -------
    function
        A function accept a list of lines as input and return
        a generator to yield section one by one.

    """

    def parser(lines):
        curr = []
        curr_type = None
        for line in _line_generator(lines, **kwargs):
            # if we find another line, return the previous section
            line_type = get_line_key(line)

            # changed line type
            if line_type != curr_type:
                if curr:
                    # returning block
                    yield curr, curr_type

                    # reset curr after yield
                    curr = []

                # reset curr_type in any cases
                curr_type = line_type

            # don't append record if line type is a spacer
            if "SPACER" not in line_type:
                curr.append(line)

        # don't forget to return the last section in the file
        if curr:
            yield curr, curr_type

    return parser


# replace skbio.io.format._sequence_feature_vocabulary._parse_section_default
def _embl_parse_section_default(
    lines, label_delimiter=None, join_delimiter=" ", return_label=False
):
    """Parse sections in default way.

    Do 2 things:
        1. split first line with label_delimiter for label
        2. join all the lines into one str with join_delimiter.
    """
    data = []
    label = None
    line = lines[0]

    # take the first line, divide the key from the text
    items = line.split(label_delimiter, 1)

    if len(items) == 2:
        label, section = items
    else:
        label = items[0]
        section = ""

    # append the text of the first element in a empty array
    data.append(section)

    # Then process all the elements with the same embl key. remove the key
    # and append all the text in the data array
    data.extend(line.split(label_delimiter, 1)[-1] for line in lines[1:])

    # Now concatenate the text using join_delimiter. All content with the same
    # key will be placed in the same string. Strip final "\n
    data = join_delimiter.join(i.strip() for i in data)

    # finally return the merged text content, and the key if needed
    if return_label:
        return label, data
    else:
        return data


# parse an embl reference record.
def _parse_reference(lines):
    """Parse single REFERENCE field."""
    # parsed reference will be placed here
    res = {}

    # define a section splitter with _embl_yield_section function defined in
    # this module
    section_splitter = _embl_yield_section(
        lambda line: _get_embl_key(line), skip_blanks=True, strip=False
    )

    # now itereta along sections (lines of the same type)
    for section, section_name in section_splitter(lines):
        # this function append all data in the same keywords. A list of lines
        # as input (see skbio.io.format._sequence_feature_vocabulary)
        label, data = _embl_parse_section_default(
            section, join_delimiter=" ", return_label=True
        )

        res[label] = data

    # now RX (CROSS_REFERENCE) is a joined string of multiple values. To get
    # back to a list of values you can use: re.compile("([^;\s]*); ([^\s]*)")

    # search for pubmed record, and add the PUBMED key
    if "RX" in res:
        match = re.search(r"PUBMED; (\d+)\.", res["RX"])

        if match:
            # add pubmed notation
            res["PUBMED"] = match.groups()[0]

    # fix RP field like genbank (if exists), Ie: (bases 1 to 63)
    if "RP" in res:
        match = re.search(r"(\d+)-(\d+)", res["RP"])

        if match:
            # fix rp fields
            res["RP"] = "(bases {start} to {stop})".format(
                start=match.groups()[0], stop=match.groups()[1]
            )

    # return translated keys (EMBL->GB)
    return _translate_keys(res)


def _serialize_reference(header, obj, cross_references, indent=5):
    """Serialize a list of references."""
    reference = []
    sort_order = ["RC", "RP", "RX", "RG", "RA", "RT", "RL"]

    # deal with RX pattern and RP pattern
    RX = re.compile(r"([^;\s]*); ([^\s]*)")
    RP = re.compile(r"bases (\d+) to (\d+)")

    # create a copy of obj, that can be changed. I need to delete values or
    # adding new ones
    obj = copy.deepcopy(obj)

    # obj is a list of references. Now is a copy of metadata[SOURCE]
    for i, data in enumerate(obj):
        # get the reference number (as the iteration number)
        embl_key = "RN"

        # get cross_references
        if cross_references:
            cross_reference = cross_references[i]

            # append cross reference [i] to data (obj[i]) (if they exists)
            if cross_reference:
                data["CROSS_REFERENCE"] = cross_reference

            # delete PUBMED key (already present in CROSS_REFERENCE)
            if "PUBMED" in data:
                del data["PUBMED"]

        else:
            # no cross reference, do I have PUBMED in data?
            if "PUBMED" in data:
                # add a fake CROSS_REFERENCE
                data["CROSS_REFERENCE"] = "PUBMED; %s." % data["PUBMED"]

        # get an embl wrapper
        wrapper = _get_embl_wrapper(embl_key, indent)

        # define wrapped string and add RN to embl data
        reference += wrapper.wrap("[{RN}]".format(RN=i + 1))

        # now process each record for references
        for embl_key in sort_order:
            # get internal key (genbank like key)
            key = _translate_key(embl_key)

            # have I this reference in my reference data?
            if key not in data:
                continue

            # if yes, define wrapper
            wrapper = _get_embl_wrapper(embl_key, indent)

            # data could have newlines
            records = data[key].split("\n")

            for record in records:
                # strip after newlines
                record = record.strip()

                # define wrapped string. beware RX
                if embl_key == "RX":
                    for match in re.finditer(RX, record):
                        source, link = match.groups()
                        # join text
                        cross_reference = "; ".join([source, link])
                        reference += wrapper.wrap(cross_reference)

                # RP case
                elif embl_key == "RP":
                    match = re.search(RP, record)

                    # if I have position, re-define RP key
                    if match:
                        record = "%s-%s" % match.groups()
                        reference += wrapper.wrap(record)

                    # if not, ignore RP key
                    else:
                        continue

                # all the other cases, go in wrapper as they are
                else:
                    reference += wrapper.wrap(record)

        # add a spacer between references (but no at the final reference)
        # cause the caller will add spacer
        if (i + 1) < len(obj):
            reference += ["XX"]

    # now define a string and add a final "\n"
    s = "\n".join(reference) + "\n"

    # and return it
    return s


# parse an embl reference record.
def _parse_source(lines):
    """Parse single SOURCE field."""
    # parsed reference will be placed here
    res = {}

    # define a section splitter with _embl_yield_section function defined in
    # this module
    section_splitter = _embl_yield_section(
        lambda line: _get_embl_key(line), skip_blanks=True, strip=False
    )

    # now itereta along sections (lines of the same type)
    for section, section_name in section_splitter(lines):
        # this function append all data in the same keywords. A list of lines
        # as input (see skbio.io.format._sequence_feature_vocabulary)
        label, data = _embl_parse_section_default(
            section, join_delimiter=" ", return_label=True
        )

        res[label] = data

    # return translated keys
    return _translate_keys(res)


def _serialize_source(header, obj, indent=5):
    """Serialize SOURCE.

    Parameters
    ----------
    header: str
        The section header.
    obj : dict
        A dictionary containing key-value pairs representing the attributes
        of the SOURCE section.
    indent : int, optional
        The number of spaces used to indent the serialized SOURCE section.
        Defaults to 5.

    """
    source = []

    # treat taxonomy and all others keys
    for key in ["ORGANISM", "taxonomy", "organelle"]:
        # get data to serielize
        data = obj.get(key)

        # if key is not defined (eg. organelle, continue)
        if data is None:
            continue

        # get embl key for my key (eg, taxonomy -> OC)
        embl_key = REV_KEYS_TRANSLATOR.get(key, key)

        # get an embl wrapper
        wrapper = _get_embl_wrapper(embl_key, indent)

        # define wrapped string
        source += wrapper.wrap(data)

    # now define a string and add a final "\n"
    s = "\n".join(source) + "\n"

    # and return it
    return s


def _parse_sequence(lines):
    """Parse the sequence section for sequence."""
    # result array
    sequence = []

    for line in lines:
        # ignore record like:
        # SQ   Sequence 275 BP; 64 A; 73 C; 88 G; 50 T; 0 other;
        if line.startswith("SQ"):
            continue

        # remove the numbers inside strings. revome spaces around string
        items = [i for i in line.split() if not i.isdigit()]

        # append each sequence items to sequence list
        sequence += items

    return "".join(sequence)


def _serialize_sequence(obj, indent=5):
    """Serialize seq to SQ.

    Parameters
    ----------
    obj : DNA, RNA, Sequence Obj
        A DNA, RNA, or Sequence object representing the biological sequence to be
        serialized.
    indent : int, optional
        The number of spaces used to indent the serialized sequence. Defaults to 5.

    """
    # a flag to determine if I wrote header or not
    flag_header = False

    # magic numbers: there will be 60 letters (AA, bp) on each line
    chunk_size = 60

    # letters (AA, bp) will be grouped by 10: each group is divided by
    # one space from each other
    frag_size = 10

    # fasta sequence will have indent spaces on the left, chunk_size/frag_size
    # groups of frag_size letters separated by n-1 groups of single spaces,
    # then the sequence length aligned on the right to get a string of
    # line_size. Setting left and right padding for semplicity
    pad_right = 65  # there are also 5 columns for indentation
    pad_left = 10  # sequence number will be in the last 10 columns

    # get sequence as a string with lower letters (uniprot will be upper!)
    seq = str(obj).lower()

    # count bases in sequence. Frequencies returns a dictionary of occurences
    # of A,C,G,T. Sequences are stored always in capital letters
    freq = obj.frequencies()

    # get values instead of popping them: I can't assure that the letter T,
    # for example, is always present
    n_a = freq.get("A", 0)
    n_c = freq.get("C", 0)
    n_g = freq.get("G", 0)
    n_t = freq.get("T", 0)

    # this will be the count of all others letters (more than ACGT)
    n_others = len(obj) - (n_a + n_c + n_g + n_t)

    # define SQ like this:
    # SQ   Sequence 275 BP; 63 A; 72 C; 88 G; 52 T; 0 other;
    SQ = (
        "SQ   Sequence {size} {unit}; {n_a} A; {n_c} C; {n_g} G; "
        + "{n_t} T; {n_others} other;\n"
    )

    # TODO: deal with protein SQ: they have a sequence header like:
    # SQ   SEQUENCE   256 AA;  29735 MW;  B4840739BF7D4121 CRC64;

    # apply format
    SQ = SQ.format(
        size=len(obj),
        unit=obj.metadata["LOCUS"]["unit"].upper(),
        n_a=n_a,
        n_c=n_c,
        n_g=n_g,
        n_t=n_t,
        n_others=n_others,
    )

    for i in range(0, len(seq), chunk_size):
        line = seq[i : i + chunk_size]
        # pad string left and right
        s = "{indent}{s:<{pad_right}}{pos:>{pad_left}}\n".format(
            indent=" " * indent,
            s=chunk_str(line, frag_size, " "),
            pad_left=pad_left,
            pos=i + len(line),
            pad_right=pad_right,
        )

        if not flag_header:
            # First time here. Add SQ header to sequence
            s = SQ + s

            # When I added header, I need to turn off this flag
            flag_header = True

        yield s


def _embl_parse_feature_table(lines, length):
    """Parse embl feature tables."""
    # define interval metadata
    imd = IntervalMetadata(length)

    # get only FT records, and remove key from line
    lines = [line[2:] for line in lines if line.startswith("FT")]

    # magic number 19: after key removal, the lines of each feature
    # are indented with 19 spaces.
    feature_indent = " " * 19

    section_splitter = _yield_section(
        lambda x: not x.startswith(feature_indent), skip_blanks=True, strip=False
    )

    for section in section_splitter(lines):
        _parse_single_feature(section, imd)
    return imd


def _serialize_feature_table(intervals, indent=21):
    """Serialize a list of ``Interval`` objects into EMBL format.

    Parameters
    ----------
    intervals : list of ``Interval``
        A list of Interval objects representing the intervals to be serialized.
    indent : int, optional
        The number of spaces to indent each serialized feature. Defaults to 21.

    """
    # define a embl wrapper object. I need to replace only the first two
    # characters from _serialize_single_feature output
    wrapper = _get_embl_wrapper("FT", indent=2, subsequent_indent=21)

    for intvl in intervals:
        tmp = _serialize_single_feature(intvl, indent)
        output = []

        # I need to remove two spaces, cause I will add a FT key
        for line in tmp.split("\n"):
            output += wrapper.wrap(line[2:])

        # re add newline between elements, and a final "\n"
        yield "\n".join(output) + "\n"


def _parse_date(lines, label_delimiter=None, return_label=False):
    """Parse embl date records."""
    # take the first line, and derive a label
    label = lines[0].split(label_delimiter, 1)[0]

    # read all the others dates and append to data array
    data = [line.split(label_delimiter, 1)[-1] for line in lines]

    # strip returned data
    data = [i.strip() for i in data]

    # finally return data array, and the key if needed
    if return_label:
        return label, data
    else:
        return data


def _serialize_date(embl_key, date_list, indent=5):
    """Serialize date line.

    Parameters
    ----------
    embl_key : str
        The EMBL key ID corresponding to the date line.
    date_list : list
        A list of dates associated with the sequence entry.
    indent : int, optional
        The number of spaces used to indent the serialized date line. Defaults to 5.

    """
    # get an embl wrapper
    wrapper = _get_embl_wrapper(embl_key, indent)

    # # serialize date and return them as a string
    return _serialize_list(wrapper, date_list)


def _serialize_comment(embl_key, obj, indent=5):
    """Serialize comment (like Assembly)."""
    # obj is a string, Split it by newlines
    data = obj.split("\n")

    # get an embl wrapper
    wrapper = _get_embl_wrapper(embl_key, indent)

    # serialize data and return it
    return _serialize_list(wrapper, data)


def _serialize_dbsource(embl_key, obj, indent=5):
    """Serialize DBSOURCE."""
    # data are stored like 'SILVA-LSU; LK021130. SILVA-SSU; LK021130. ...
    # I need to split string after final period (not AAT09660.1)

    # deal with re pattern. A pattern to find a period as end of sentence
    DR = re.compile(r"\.\s")

    # splitting by this pattern, I will have
    # ["SILVA-LSU; LK021130", "SILVA-SSU; LK021130", ...]
    # I need that each of them will be in a DR record.

    # get an embl wrapper
    wrapper = _get_embl_wrapper(embl_key, indent)

    # serialize data and return it. Split dbsource using re. Add a
    # final period between elements since I removed it by splitting
    return _serialize_list(wrapper, re.split(DR, obj), sep=".\n")


def _parse_assembly(lines):
    """Parse embl assembly records."""
    output = []

    # first line is header, skip it
    for line in lines[1:]:
        data = line.split()

        # data could have comp feature or not. First element in data is 'AS'
        if len(data) == 5:
            res = dict(
                zip(
                    ["local_span", "primary_identifier", "primary_span", "comp"],
                    data[1:],
                )
            )

        elif len(data) == 4:
            res = dict(
                zip(
                    ["local_span", "primary_identifier", "primary_span", "comp"],
                    data[1:] + [""],
                )
            )

        else:
            raise EMBLFormatError("Can't parse assembly line %s" % line)

        # append res to output
        output += [res]

    return output


# Map a function to each section of the entry
_PARSER_TABLE = {
    "LOCUS": _parse_id,
    "SOURCE": _parse_source,
    "DATE": _parse_date,
    "REFERENCE": _parse_reference,
    "FEATURES": _embl_parse_feature_table,
    "ORIGIN": _parse_sequence,
    "ASSEMBLY": _parse_assembly,
}

# for writer functions
_SERIALIZER_TABLE = {
    "LOCUS": _serialize_id,
    "SOURCE": _serialize_source,
    "DATE": _serialize_date,
    "REFERENCE": _serialize_reference,
    "FEATURES": _serialize_feature_table,
    "COMMENT": _serialize_comment,
    "DBSOURCE": _serialize_dbsource,
}
