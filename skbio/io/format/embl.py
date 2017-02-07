"""
EMBL format (:mod:`skbio.io.format.embl`)
=========================================

.. currentmodule:: skbio.io.format.embl

EMBL format stores sequence and its annotation together. The start of the
annotation section is marked by a line beginning with the word "ID". The start
of sequence section is marked by a line beginning with the word "SQ" and the
end of the section is marked by a line with only "//". More information on
EMBL file format can be found here [1]_.

The EMBL file may ens with .embl or .txt extension. An example of EMBL file can
be seen here [2]_.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.RNA`                                      |
+------+------+---------------------------------------------------------------+
|No    |No    |:mod:`skbio.sequence.Protein`                                  |
+------+------+---------------------------------------------------------------+
|Yes   |No    | generator of :mod:`skbio.sequence.Sequence` objects           |
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

#TODO: Add a Feature section with table, like genbank

Examples
--------

Reading and Writing EMBL Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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


References
----------
.. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/release/doc/usrman.txt
.. [2] http://www.ebi.ac.uk/ena/data/view/X56734&display=text

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
from functools import partial

# skbio modules
from skbio.io import create_format, EMBLFormatError
from skbio.io.format._base import (_line_generator, _get_nth_sequence)
from skbio.sequence import Sequence, DNA, RNA, Protein

# look at skbio.io.registry to have an idea on how to define this class
embl = create_format('embl')

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
                   # 'PR': 'PROJECT_IDENTIFIER',
                   'DT': 'DATE',
                   'DE': 'DEFINITION',
                   # 'GN': 'GENE_NAME',  # uniprot specific
                   'KW': 'KEYWORDS',
                   # Source (taxonomy and classification)
                   'OS': 'ORGANIMS',
                   'OC': 'taxonomy',
                   # 'OG': 'organelle'
                   # reference keys
                   'RA': 'AUTHORS',
                   'RP': 'REFERENCE',
                   # 'RC': 'COMMENT',
                   'RX': 'CROSS_REFERENCE',
                   # 'RG': 'GROUP',
                   'RT': 'TITLE',
                   'RL': 'JOURNAL',
                   # Cross references
                   'DR': 'DBSOURCE',
                   'CC': 'COMMENT',
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
                   'PA': 'ACCESSION',
                   # 'PR': 'PROJECT_IDENTIFIER',
                   'DT': 'DATE',
                   'DE': 'DESCRIPTION',
                   # 'GN': 'GENE_NAME',  # uniprot specific
                   'KW': 'KEYWORDS',
                   # Source (taxonomy and classification)
                   'OS': 'SOURCE',
                   'OC': 'SOURCE',
                   # 'OG': 'SOURCE',
                   # reference keys
                   'RA': 'REFERENCE',
                   'RP': 'REFERENCE',
                   # 'RC': 'REFERENCE',
                   'RX': 'REFERENCE',
                   # 'RG': 'REFERENCE',
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
        key = _get_embl_key(section[0])

        # even section have different keys, (RA, RP, ...), I need to get the
        # correct section to call the appropriate method
        header = _get_embl_section(key)

        # debug
        # print("header >%s<" %(header))
        # print("section >%s<" %(section))
        # print("metadata >%s<" %(metadata))

        # search for a specific method in PARSER_TABLE or set
        # _embl_parse_section_default
        parser = _PARSER_TABLE.get(
            header, _embl_parse_section_default)

        if header == 'FEATURES':
            # This requires 'ID' line parsed before 'FEATURES', which should
            # be true and is implicitly checked by the sniffer. This is true
            # since the first section is parsed by the last else condition

            # partials add arguments to previous defined functions
            parser = partial(
                parser, length=metadata['LOCUS']['size'])

        # call function on section
        parsed = parser(section)

        # reference can appear multiple times
        if header == 'REFERENCE':
            if header in metadata:
                metadata[header].append(parsed)

            else:
                metadata[header] = [parsed]

        elif header == 'ORIGIN':
            sequence = parsed

#        elif header == 'FEATURES':
#            interval_metadata = parsed

        # parse all the others sections
        else:
            metadata[header] = parsed

    return sequence, metadata, interval_metadata


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
    # key will be placed in the same string
    data = join_delimiter.join(i.strip() for i in data)

    # finally return the merged text content, and the key if needed
    if return_label:
        return label, data
    else:
        return data


# parse an embl reference record
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
        label, data = _embl_parse_section_default(
            section, join_delimiter=' ', return_label=True)

        res[label] = data

    # return translated keys
    return _translate_keys(res)


def _parse_sequence(lines):
    '''Parse the sequence section for sequence.'''

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


# boilerplate for parse feature functions
def _parse_feature_table(lines, length):
    """Parse feature tables"""

    return


# Map a function to each section of the entry
_PARSER_TABLE = {
    'LOCUS': _parse_id,
    # 'SOURCE': _parse_source,
    'REFERENCE': _parse_reference,
    'FEATURES': _parse_feature_table,
    'ORIGIN': _parse_sequence
    }
