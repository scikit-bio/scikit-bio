"""
EMBL format (:mod:`skbio.io.format.embl`)
=========================================

.. currentmodule:: skbio.io.format.embl

EMBL format stores sequence and its annotation together. The start of the
annotation section is marked by a line beginning with the word "ID". The start
of sequence section is marked by a line beginning with the word "SQ" and the
end of the section is marked by a line with only "//". More information on
EMBL file format can be found here [1]_.

The EMBL file may ends with .embl or .txt extension. An example of EMBL file
can be seen here [2]_.

Feature level products
^^^^^^^^^^^^^^^^^^^^^^

Feature-level products contain nucleotide sequence and related annotation
derived from submitted ENA assembled and annotated sequences. Data are
distributed in flatfile format, similar to that of parent ENA records,
with each flatfile representing a single feature. While only the sequence
of the feature is included in such entries, features are derived from the
parent entry, and can't be applied as interval metadata. For such reason,
interval metatdata are ignored from Feature-level products, as they will be
ignored by subsetting a generic Sequence object. More information on
Feature-level product can be found here [3]_.

Format Support
--------------
**Has Sniffer: Yes**
**NOTE: No protein support at the moment**

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
**State: Experimental as of 0.5.1.**

Sections before ``FH (Feature Header)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All the sections before ``FH (Feature Header)`` will be read into the attribute
of ``metadata``. The header and its content of a section is stored as
a pair of key and value in ``metadata``. For the ``RN (Reference Number)``
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

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# std modules
import re
import textwrap

from functools import partial

# skbio modules
from skbio.io import create_format, EMBLFormatError
from skbio.io.format._base import (_line_generator, _get_nth_sequence)
from skbio.io.format._sequence_feature_vocabulary import (
    _yield_section, _parse_single_feature, _serialize_section_default,
    _serialize_single_feature)
from skbio.metadata import IntervalMetadata
from skbio.sequence import Sequence, DNA, RNA, Protein
from skbio.util._misc import chunk_str


# look at skbio.io.registry to have an idea on how to define this class
embl = create_format('embl')

# This list is ordered used to read and write embl file.
_HEADERS = [
    'LOCUS',
    'ACCESSION',
    'PARENT_ACCESSION',
    'PROJECT_IDENTIFIER',
    'DATE',
    'DEFINITION',
    'GENE_NAME',
    'KEYWORDS',
    'SOURCE',
    'REFERENCE',
    'DBSOURCE',
    'COMMENT',
    'FEATURES'
    ]


# embl has a series of keys different from genbank; moreover keys are not so
# easy to understand (eg. RA for AUTHORS). Here is a dictionary of keys
# conversion (EMBL->GB). All the unspecified keys will remain in embl format
KEYS_TRANSLATOR = {
                   # identification
                   'ID': 'LOCUS',
                   'AC': 'ACCESSION',
                   # PA means PARENT ACCESSION (?) and applies to
                   # feature-level-products entries
                   'PA': 'PARENT_ACCESSION',
                   'PR': 'PROJECT_IDENTIFIER',
                   'DT': 'DATE',
                   'DE': 'DEFINITION',
                   'GN': 'GENE_NAME',  # uniprot specific
                   'KW': 'KEYWORDS',
                   # Source (taxonomy and classification)
                   'OS': 'ORGANISM',
                   'OC': 'taxonomy',
                   'OG': 'organelle',
                   # reference keys
                   'RA': 'AUTHORS',
                   'RP': 'REFERENCE',
                   'RC': 'REFERENCE_COMMENT',
                   'RX': 'CROSS_REFERENCE',
                   'RG': 'GROUP',
                   'RT': 'TITLE',
                   'RL': 'JOURNAL',
                   # Cross references
                   'DR': 'DBSOURCE',
                   'CC': 'COMMENT',
                   # features
                   'FH': 'FEATURES',
                   'FT': 'FEATURES',
                   'SQ': 'ORIGIN',
                   }

# the original _yield_section divides entries in sections relying on spaces
# (the same section has the same level of indentation). Uniprot entries have
# a key for each line, so to divide record in sections I need to define a
# corresponance for each key to section, then I will divide a record in
# sections using these section name. Some keys are commented out since I don't
# find them in example. What I have to do if I find them?
KEYS_2_SECTIONS = {
                   # identification
                   'ID': 'LOCUS',
                   'AC': 'ACCESSION',
                   # PA means PARENT ACCESSION (?) and applies to
                   # feature-level-products entries
                   'PA': 'PARENT_ACCESSION',
                   'PR': 'PROJECT_IDENTIFIER',
                   'DT': 'DATE',
                   'DE': 'DEFINITION',
                   'GN': 'GENE_NAME',  # uniprot specific
                   'KW': 'KEYWORDS',
                   # Source (taxonomy and classification)
                   'OS': 'SOURCE',
                   'OC': 'SOURCE',
                   'OG': 'SOURCE',
                   # reference keys
                   'RA': 'REFERENCE',
                   'RP': 'REFERENCE',
                   'RC': 'REFERENCE',
                   'RX': 'REFERENCE',
                   'RG': 'REFERENCE',
                   'RT': 'REFERENCE',
                   'RL': 'REFERENCE',
                   'RN': 'SPACER',
                   # Cross references
                   'DR': 'DBSOURCE',
                   'CC': 'COMMENT',
                   # 'AH': 'ASSEMBLY,
                   # 'AS': 'ASSEMBLY'
                   'FH': 'FEATURES',
                   'FT': 'FEATURES',
                   # sequence
                   'SQ': 'ORIGIN',
                   '  ': 'ORIGIN',
                   # 'CO': 'ORIGIN,'
                   # spacer (discarded)
                   'XX': 'SPACER'
                  }


# for convenience: I think such functions are more readadble in lambda
# functions
def _get_embl_key(line):
    """Return first part of a string as a embl key (ie 'AC M14399;' -> 'AC')"""

    # embl keys have a fixed size of 2 chars
    return line[:2]


# get embl key from value
# http://stackoverflow.com/questions/8023306/get-key-by-value-in-dictionary
def _get_embl_key_by_value(value, mydict=KEYS_TRANSLATOR):
    """Return the key(s) associated to a value from a dictionary, or the
    value if no key is defined"""

    # try to get a key (keys) from a value in dictionary
    try:
        return list(mydict.keys())[list(mydict.values()).index(value)]

    # returns value if key is not present
    except ValueError:
        return value


def _get_embl_section(line):
    """Return the embl section from uniprot key(ie 'RA' -> 'REFERENCE')"""

    # get embl key
    key = _get_embl_key(line)

    # get embl section from key
    section = KEYS_2_SECTIONS[key]

    # debug
    # print("_get_embl_section line: >%s<"%(line))
    # print("_get_embl_section key: >%s<"%(key))
    # print("_get_embl_section section: >%s<"%(section))

    return section


def _translate_key(key):
    """A method to translate a single key. Return key itself if no traslation
    is defined"""

    return KEYS_TRANSLATOR.get(key, key)


# a method to translate keys for a dict object. All keys not defined will
# remain the same
def _translate_keys(data):
    """Translate a dictionary of uniprot key->value in a genbank like
    dictionary of key values. Keep old keys if no translation is defined"""

    # get all keys to validate
    old_keys = data.keys()

    # a new dictionary of results
    new_data = {}

    # I can't replace keys in original values, sometimes test will fails. So, I
    # create a new copy. This is a strange behaviour, I don't understand
    for old_key in old_keys:
        new_key = _translate_key(old_key)
        new_data[new_key] = data[old_key]

    # returning translated keys
    return new_data


# define a default textwrap.Wrapper for embl
def _get_embl_wrapper(embl_key, indent=5, subsequent_indent=None):
    """Returns a textwrap.TextWrapper for embl records"""

    # define the string to prepen (eg "OC   ")
    prepend = '{key:<{indent}}'.format(key=embl_key, indent=indent)

    # deal with 2° strings and more
    if subsequent_indent is None:
        subsequent_prepend = prepend

    else:
        subsequent_prepend = '{key:<{indent}}'.format(
            key=embl_key, indent=subsequent_indent)

    # define a text wrapper object
    wrapper = textwrap.TextWrapper(
        initial_indent=prepend,
        subsequent_indent=subsequent_prepend,
        width=80
        )

    return wrapper


def _serialize_list(embl_wrapper, data):
    """Serialize a list of obj using a textwrap.TextWrapper instance"""

    # the output array
    output = []

    for line in data:
        output += embl_wrapper.wrap(line)

    # merge dates in one string
    output = "\n".join(output) + "\n"

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
    # record = _get_nth_sequence(_parse_embls(fh), seq_num)
    # return _construct(record, Protein, **kwargs)
    raise EMBLFormatError("There's no protein support for EMBL record")


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
    # _serialize_single_embl(obj, fh)
    raise EMBLFormatError("There's no protein support for EMBL record")


def _construct(record, constructor=None, **kwargs):
    '''Construct the object of Sequence, DNA, RNA, or Protein.'''

    # sequence, metadata and interval metadata
    seq, md, imd = record

    if 'lowercase' not in kwargs:
        kwargs['lowercase'] = True

    if constructor is None:
        unit = md['LOCUS']['unit']
        if unit == 'bp':
            # RNA mol type has T instead of U for genbank from from NCBI
            constructor = DNA

        elif unit == 'aa':
            # no protein support, at the moment
            # constructor = Protein
            raise EMBLFormatError("There's no protein support for EMBL record")

    if constructor == RNA:
        return DNA(
            seq, metadata=md, interval_metadata=imd, **kwargs).transcribe()
    else:
        return constructor(
            seq, metadata=md, interval_metadata=imd, **kwargs)


# looks like the genbank _parse_genbank
def _parse_embls(fh):
    """Chunck multiple EMBL records by '//', and returns a generator"""

    data_chunks = []
    for line in _line_generator(fh, skip_blanks=True, strip=False):
        if line.startswith('//'):
            yield _parse_single_embl(data_chunks)
            data_chunks = []
        else:
            data_chunks.append(line)


def _parse_single_embl(chunks):
    metadata = {}
    interval_metadata = None
    sequence = ''

    # define a section splitter with _embl_yield_section function defined in
    # this module (return the embl section by embl key). returns generator for
    # each block with different line type
    section_splitter = _embl_yield_section(
        lambda line: _get_embl_section(line),
        skip_blanks=True,
        strip=False)

    # process each section, like genbank does
    for section in section_splitter(chunks):
        # key is line type (eg ID, AC, ...)
        embl_key = _get_embl_key(section[0])

        # even section have different keys, (RA, RP, ...), I need to get the
        # correct section to call the appropriate method (eg REFERENCE)
        embl_header = _get_embl_section(embl_key)

        # search for a specific method in PARSER_TABLE or set
        # _embl_parse_section_default
        parser = _PARSER_TABLE.get(
            embl_header, _embl_parse_section_default)

        if embl_header == 'FEATURES':
            # This requires 'ID' line parsed before 'FEATURES', which should
            # be true and is implicitly checked by the sniffer. This is true
            # since the first section is parsed by the last else condition

            # partials add arguments to previous defined functions
            parser = partial(
                parser, metadata=metadata)

        elif embl_header == "COMMENT":
            # mantain newlines in comments
            # partials add arguments to previous defined functions
            parser = partial(
                parser, join_delimiter="\n")

        # call function on section
        parsed = parser(section)

        # reference can appear multiple times
        if embl_header == 'REFERENCE':
            if embl_header in metadata:
                metadata[embl_header].append(parsed)

            else:
                metadata[embl_header] = [parsed]

        elif embl_header == 'ORIGIN':
            sequence = parsed

        elif embl_header == 'FEATURES':
            interval_metadata = parsed

        # parse all the others sections (DATE, SOURCE, ...)
        else:
            metadata[embl_header] = parsed

    return sequence, metadata, interval_metadata


# main function for writer methods
def _serialize_single_embl(obj, fh):
    '''Write a EMBL record.

    Always write it in ENA canonical way:
    1. sequence in lowercase (uniprot are uppercase)
    2. 'u' as 't' even in RNA molecules.

    Parameters
    ----------
    obj : Sequence or its child class

    '''
    # write out the headers
    md = obj.metadata

    # embl has a different magick number than embl
    serialize_default = partial(
        _serialize_section_default, indent=5)

    for header in _HEADERS:
        serializer = _SERIALIZER_TABLE.get(
            header, serialize_default)

        # header needs to be convert into embl, or matained as it is
        # if no conversion could be defined
        embl_header = _get_embl_key_by_value(header)

        # this is true also for locus line
        if header in md:
            # call the serializer function
            out = serializer(embl_header, md[header])
            # test if 'out' is a iterator.
            # cf. Effective Python Item 17
            if iter(out) is iter(out):
                for s in out:
                    fh.write(s)

            else:
                fh.write(out)

            # add spacer between sections
            fh.write("XX\n")

        if header == 'FEATURES':
            if obj.has_interval_metadata():
                # magic number 21: the amount of indentation before
                # feature table starts as defined by INSDC
                indent = 21
                feature_key = "FH   Key"
                fh.write('{header:<{indent}}Location/Qualifiers\n'.format(
                    header=feature_key, indent=indent))

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
    fh.write('//\n')


def _parse_id(lines):
    """
    The ID (IDentification) line is always the first line of an entry. The
    format of the ID line is:
    ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
    The tokens represent:
       1. Primary accession number
       2. Sequence version number
       3. Topology: 'circular' or 'linear'
       4. Molecule type (see note 1 below)
       5. Data class (see section 3.1)
       6. Taxonomic division (see section 3.2)
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
    pattern = re.compile(r'ID'
                         ' +([^\s]+);'     # ie: CD789012
                         ' +SV ([0-9]+);'  # 4
                         ' +(\w+);'        # linear
                         ' +([^;]+);'      # genomic DNA
                         ' +(\w+);'        # HTG
                         ' +(\w+);'        # MAM
                         ' +(\d+)'         # 500
                         ' +(\w+)\.$')     # BP

    # search it
    matches = re.match(pattern, line)

    try:
        res = dict(zip(
            ['accession', 'version', 'shape', 'mol_type',
             'class', 'division', 'size', 'unit'],
            matches.groups()))
    except:
        raise EMBLFormatError(
            "Could not parse the ID line:\n%s" % line)

    # those values are integer
    res['size'] = int(res['size'])
    res['version'] = int(res['version'])

    # unit are in lower cases in others modules
    res['unit'] = res['unit'].lower()

    # returning parsed attributes
    return res


def _serialize_id(header, obj, indent=5):
    '''Serialize ID line.

    Parameters
    ----------
    obj : dict
    '''

    # use 'or' to convert None to ''
    kwargs = {k: v or '' for k, v in obj.items()}

    # then unit is in upper cases
    kwargs["unit"] = kwargs["unit"].upper()

    # return first line
    return ('{header:<{indent}}{accession}; SV {version}; {shape}; '
            '{mol_type}; {class}; {division}; {size} {unit}.\n').format(
                header=header, indent=indent, **kwargs)


# For non-coding, rRNA and spacer records, where a feature-level ID has not
# previously existed, the ID, e.g. AB012758.1:1..40:tRNA, has a complex
# format to ensure that it is unique and unambiguous. The structure of
# the ID may be represented as:
# <accession>.<version>:<feature location>:<feature name>[:ordinal]
def _parse_accession(locus_dict):
    """Parse accession string like :
        <accession>.<version>:<feature location>:<feature name>[:ordinal]"""

    # locus_dict is the dictionary read by _parse_id
    accession = locus_dict.get("accession")

    # define a regular expression to read accession
    pattern = re.compile("(\w+)\.([0-9]+)\:([^:]+)\:(\w+)")
    matches = re.match(pattern, accession)

    # read data
    res = dict(zip(["parent_accession", "version", "feature_location",
                    "feature_name"], matches.groups()))

    # read locations
    start, stop = res.get("feature_location").split("..")

    # fix values. Convert in O base coordinates
    res["start"] = int(start) - 1
    res["stop"] = int(stop)
    res["size"] = abs(res["stop"] - res["start"])
    res["version"] = int(res["version"])

    # return parsed accession
    return res


# replace skbio.io.format._sequence_feature_vocabulary.__yield_section
def _embl_yield_section(get_line_key, **kwargs):
    '''Returns function that returns successive sections from file.

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
    '''
    def parser(lines):
        curr = []
        curr_type = None
        for line in _line_generator(lines, **kwargs):
            # if we find another line, return the previous section
            line_type = get_line_key(line)

            # debug
            # print("_embl_yield_section line: >%s<" %(line))
            # print("_embl_yield_section line_type: >%s<" %(line_type))

            # changed line type
            if line_type != curr_type:
                if curr:
                    # debug
                    # print("_embl_yield_section curr_type: >%s<" %(curr_type))
                    # print("_embl_yield_section curr: >%s<" %(curr))

                    # returning block
                    yield curr

                    # reset curr after yield
                    curr = []

                # reset curr_type in any cases
                curr_type = line_type

            # don't append record if line type is a spacer
            if 'SPACER' not in line_type:
                curr.append(line)

            # debug
            # else:
            #    print("Ignoring %s" %(line))

        # don't forget to return the last section in the file
        if curr:
            yield curr

    return parser


# replace skbio.io.format._sequence_feature_vocabulary._parse_section_default
def _embl_parse_section_default(
        lines, label_delimiter=None, join_delimiter=' ', return_label=False):
    '''Parse sections in default way.

    Do 2 things:
        1. split first line with label_delimiter for label
        2. join all the lines into one str with join_delimiter.
    '''

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


def _embl_parse_section_newlines(
        lines, label_delimiter=None, join_delimiter=' ', return_label=False):

    '''Parse sections with newlines. Keep "\n" when sentences end

    Do 3 things:
        1. split first line with label_delimiter for label
        3. Search for end of sentences, keep their "\n"
        2. join all the lines into one str with join_delimiter.
    '''

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

    # deal with data like this
    # RL   Submitted (27-JUN-2016) to the INSDC.
    # RL   Key Laboratory of Coastal Zone Environment Processes and Ecological
    # RL   Remediation, Yantai Institute of Coastal Zone Research (YIC),
    # RL   Chinese Academy of Sciences (CAS), 17 Chunhui Road, Laishan
    # RL   District, Yantai, Shandong 264003, China

    # in order to put "\n" after the first line
    # define end of sentence pattern
    pattern = re.compile("[\.]\n$")

    # find end of sentence in data
    idx = [True if re.search(pattern, i) else False for i in data]

    # now strip only when sentence continues
    data = [el.strip() if not idx[i] else el for i, el in enumerate(data)]

    # Now concatenate the text using join_delimiter. All content with the same
    # key will be placed in the same string. Strip final "\n"
    data = join_delimiter.join(data).strip()

    # finally return the merged text content, and the key if needed
    if return_label:
        return label, data
    else:
        return data


# parse an embl reference record. Is also applied on source records
def _parse_reference(lines):
    '''Parse single REFERENCE field.
    '''

    # parsed reference will be placed here
    res = {}

    # define a section splitter with _embl_yield_section function defined in
    # this module
    section_splitter = _embl_yield_section(lambda line: _get_embl_key(line),
                                           skip_blanks=True, strip=False)

    # now itereta along sections (lines of the same type)
    for section in section_splitter(lines):
        # this function append all data in the same keywords. A list of lines
        # as input (see skbio.io.format._sequence_feature_vocabulary)
        label, data = _embl_parse_section_newlines(
            section, join_delimiter=' ', return_label=True)

        res[label] = data

    # now RX (CROSS_REFERENCE) is a joined string of multiple values. To get
    # back to a list of values you can use: re.compile("([^;\s]*); ([^\s]*)")

    # return translated keys
    return _translate_keys(res)


def _serialize_reference(header, obj, indent=5):
    """Serialize a list of references"""

    reference = []
    sort_order = ["RC", "RP", "RX", "RG", "RA", "RT", "RL"]

    # deal with rx pattern
    RX = re.compile("([^;\s]*); ([^\s]*)")

    # obj is a list of references
    for i, data in enumerate(obj):
        # get the reference number
        embl_key = "RN"

        # get an embl wrapper
        wrapper = _get_embl_wrapper(embl_key, indent)

        # define wrapped string
        reference += wrapper.wrap("[{RN}]".format(RN=i+1))

        # now process each record for references
        for embl_key in sort_order:
            # get internal key
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

                else:
                    reference += wrapper.wrap(record)

        # add a spacer between references (but no at the final reference)
        # cause the caller will add spacer
        if (i+1) < len(obj):
            reference += ["XX"]

    # now define a string and add a final "\n"
    s = "\n".join(reference) + "\n"

    # and return it
    return s


def _serialize_source(header, obj, indent=5):
    '''Serialize SOURCE.

    Parameters
    ----------
    header: section header
    obj : dict
    indent : indent length
    '''

    source = []

    # treat taxonomy and all others keys
    for key in ["ORGANISM", "taxonomy", "organelle"]:
        # get data to serielize
        data = obj.get(key)

        # if key is not defined (eg. organelle, continue)
        if data is None:
            continue

        # get embl key for my key (eg, taxonomy -> OC)
        embl_key = _get_embl_key_by_value(key)

        # get an embl wrapper
        wrapper = _get_embl_wrapper(embl_key, indent)

        # define wrapped string
        source += wrapper.wrap(data)

    # now define a string and add a final "\n"
    s = "\n".join(source) + "\n"

    # and return it
    return s


def _parse_sequence(lines):
    '''Parse the sequence section for sequence.'''

    # when reading a feature-level-products accession, features are relative
    # to parent accession, so we need to parse accessions like
    # LK021130.1:74067..75610:rRNA to model a table feature like
    # FT   rRNA            LK021130.1:74067..75610
    # FT                   /gene="16S"
    # FT                   /product="16S rRNA"
    # FT                   /note="16S rRNA subunit (checked and believed
    # FT                   to be right,based on 454- and PacBio-sequencing)"
    # I could express features by resizing sequence, however the resulting
    # sequence will be different than data read. I decide to mantain
    # sequences indentical to read values, and to discard features
    # for such entries

    # result array
    sequence = []

    for line in lines:
        # debug
        # print(line)

        # ignore record like:
        # SQ   Sequence 275 BP; 64 A; 73 C; 88 G; 50 T; 0 other;
        if line.startswith('SQ'):
            continue

        # remove the numbers inside strings. revome spaces around string
        line = ''.join([i for i in line.strip() if not i.isdigit()])

        # remove space from sequence
        line = line.replace(" ", "")

        # append each line sequence to sequence list
        sequence.append(line)

    return ''.join(sequence)


def _serialize_sequence(obj, indent=5):
    '''Serialize seq to SQ.

    Parameters
    ----------
    obj : DNA, RNA, Sequence Obj
    '''

    # magic numbers
    n = 1
    line_size = 60
    frag_size = 10

    # count bases in sequence
    seq = str(obj).lower()
    n_a = seq.count("a")
    n_c = seq.count("c")
    n_g = seq.count("g")
    n_t = seq.count("t")
    n_others = len(obj) - (n_a + n_c + n_g + n_t)

    # define SQ like this:
    # SQ   Sequence 275 BP; 63 A; 72 C; 88 G; 52 T; 0 other;
    SQ = "SQ   Sequence {size} {unit}; {n_a} A; {n_c} C; {n_g} G; " +\
         "{n_t} T; {n_others} other;\n"

    # apply format
    SQ = SQ.format(size=len(obj), unit=obj.metadata["LOCUS"]["unit"].upper(),
                   n_a=n_a, n_c=n_c, n_g=n_g, n_t=n_t, n_others=n_others)

    for i in range(0, len(seq), line_size):
        line = seq[i:i+line_size]
        s = '{key:<{indent}}{s}'.format(
            key='', indent=indent, s=chunk_str(line, frag_size, ' '))
        # pad string left and right
        s = '{:<70}'.format(s) + '{:>10}\n'.format(i+len(line))
        if n == 1:
            # Add SQ header to sequence
            s = SQ + s
        n = n + line_size
        yield s


def _embl_parse_feature_table(lines, metadata):
    """Parse embl feature tables"""

    # get size of the feature
    length = metadata["LOCUS"]["size"]

    if "PARENT_ACCESSION" in metadata:
        # this is a feature-level-products entry and features are relative
        # to parent accession; in the same way a subset of a Sequence objcet
        # has no interval metadata, I will refuse to process interval
        # metadata here
        return

    # define interval metadata
    imd = IntervalMetadata(length)

    # remove feature header table
    idxs = [i for i, line in enumerate(lines) if line.startswith("FH")]
    lines = [line for i, line in enumerate(lines) if i not in idxs]

    # remove FT from the header
    lines = [line.replace("FT", "  ", 1) for line in lines]

    # magic number 21: the lines following header of each feature
    # are indented with 21 spaces.
    feature_indent = ' ' * 21

    section_splitter = _yield_section(
        lambda x: not x.startswith(feature_indent),
        skip_blanks=True, strip=False)

    for section in section_splitter(lines):
        _parse_single_feature(section, imd)
    return imd


def _serialize_feature_table(intervals, indent=21):
    '''
    Parameters
    ----------
    intervals : list of ``Interval``
    '''

    # define a embl wrrapper object. I need to replace only the first two
    # characters from _serialize_single_feature output
    wrapper = _get_embl_wrapper("FT", indent=2, subsequent_indent=21)

    for intvl in intervals:
        tmp = _serialize_single_feature(intvl, indent)
        output = []

        # I need to remove two spaces, cause I will add a FT key
        for line in tmp.split("\n"):
            output += wrapper.wrap(re.sub("^  ", "", line))

        # re add newline between elements, and a final "\n"
        yield "\n".join(output) + "\n"


def _parse_date(lines, label_delimiter=None, return_label=False):
    """Parse embl data records"""

    data = []
    line = lines[0]

    # take the first line, divide the key from the text
    label, section = line.split(label_delimiter, 1)

    # add section to data to return
    data += [section]

    # read all the others dates and append to data array
    data.extend(line.split(label_delimiter, 1)[-1] for line in lines[1:])

    # strip returned data
    data = [i.strip() for i in data]

    # finally return data array, and the key if needed
    if return_label:
        return label, data
    else:
        return data


def _serialize_date(embl_key, date_list, indent=5):
    '''Serialize date line.

    Parameters
    ----------
    header : embl key id
    date_list : a list of dates
    '''

    # get an embl wrapper
    wrapper = _get_embl_wrapper(embl_key, indent)

    # # serialize date and return them as a string
    return _serialize_list(wrapper, date_list)


def _serialize_comment(embl_key, obj, indent=5):
    """Serialize comment (like Assembly)"""

    # obj is a string, Split it by newlines
    data = obj.split("\n")

    # get an embl wrapper
    wrapper = _get_embl_wrapper(embl_key, indent)

    # serialize data and return it
    return _serialize_list(wrapper, data)


def _serialize_dbsource(embl_key, obj, indent=5):
    """Serialize DBSOURCE"""

    # deal with rx pattern
    DR = re.compile("([^;\s]*); ([^\s]*)")

    # get an embl wrapper
    wrapper = _get_embl_wrapper(embl_key, indent)

    # output list
    dbsource = []

    # obj is a string
    for match in re.finditer(DR, obj):
        source, link = match.groups()
        # join text
        cross_reference = "; ".join([source, link])
        dbsource += [cross_reference]

    # serialize data and return it
    return _serialize_list(wrapper, dbsource)


# Map a function to each section of the entry
_PARSER_TABLE = {
    'LOCUS': _parse_id,
    'SOURCE': _parse_reference,
    'DATE': _parse_date,
    'REFERENCE': _parse_reference,
    'FEATURES': _embl_parse_feature_table,
    'ORIGIN': _parse_sequence,
    }

# for writer functions
_SERIALIZER_TABLE = {
    'LOCUS': _serialize_id,
    'SOURCE': _serialize_source,
    'DATE': _serialize_date,
    'REFERENCE': _serialize_reference,
    'FEATURES': _serialize_feature_table,
    'COMMENT': _serialize_comment,
    'DBSOURCE': _serialize_dbsource,
    }
