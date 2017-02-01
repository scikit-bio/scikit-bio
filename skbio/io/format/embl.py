"""
EMBL format (:mod:`skbio.io.format.embl`)
=========================================

.. currentmodule:: skbio.io.format.embl

EMBL format stores sequence and its annotation together. The start of the
annotation section is marked by a line beginning with the word "ID". The start
of sequence section is marked by a line beginning with the word "SQ" and the
end of the section is marked by a line with only "//". More information on EMBL
file format can be found here [1]_.

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
|Yes   |No    |:mod:`skbio.sequence.Protein`                                  |
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
... OC   Spermatophyta; Magnoliophyta; eudicotyledons; Gunneridae; Pentapetalae;
... OC   rosids; fabids; Fabales; Fabaceae; Papilionoideae; Trifolieae; Trifolium.
... XX
... RN   [5]
... RP   1-1859
... RX   DOI; 10.1007/BF00039495.
... RX   PUBMED; 1907511.
... RA   Oxtoby E., Dunn M.A., Pancoro A., Hughes M.A.;
... RT   "Nucleotide and derived amino acid sequence of the cyanogenic
... RT   beta-glucosidase (linamarase) from white clover (Trifolium repens L.)";
... RL   Plant Mol. Biol. 17(2):209-219(1991).
... XX
... RN   [6]
... RP   1-1859
... RA   Hughes M.A.;
... RT   ;
... RL   Submitted (19-NOV-1990) to the INSDC.
... RL   Hughes M.A., University of Newcastle Upon Tyne, Medical School, Newcastle
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
... FT                   /experiment="experimental evidence, no additional details
... FT                   recorded"
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
... FT                   /translation="MDFIVAIFALFVISSFTITSTNAVEASTLLDIGNLSRSSFPRGFI
... FT                   FGAGSSAYQFEGAVNEGGRGPSIWDTFTHKYPEKIRDGSNADITVDQYHRYKEDVGIMK
... FT                   DQNMDSYRFSISWPRILPKGKLSGGINHEGIKYYNNLINELLANGIQPFVTLFHWDLPQ
... FT                   VLEDEYGGFLNSGVINDFRDYTDLCFKEFGDRVRYWSTLNEPWVFSNSGYALGTNAPGR
... FT                   CSASNVAKPGDSGTGPYIVTHNQILAHAEAVHVYKTKYQAYQKGKIGITLVSNWLMPLD
... FT                   DNSIPDIKAAERSLDFQFGLFMEQLTTGDYSKSMRRIVKNRLPKFSKFESSLVNGSFDF
... FT                   IGINYYSSSYISNAPSHGNAKPSYSTNPMTNISFEKHGIPLGPRAASIWIYVYPYMFIQ
... FT                   EDFEIFCYILKINITILQFSITENGMNEFNDATLPVEEALLNTYRIDYYYRHLYYIRSA
... FT                   IRAGSNVKGFYAWSFLDCNEWFAGFTVRFGLNFVD"
... XX
... SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
...      aaacaaacca aatatggatt ttattgtagc catatttgct ctgtttgtta ttagctcatt        60
...      cacaattact tccacaaatg cagttgaagc ttctactctt cttgacatag gtaacctgag       120
...      tcggagcagt tttcctcgtg gcttcatctt tggtgctgga tcttcagcat accaatttga       180
...      aggtgcagta aacgaaggcg gtagaggacc aagtatttgg gataccttca cccataaata       240
...      tccagaaaaa ataagggatg gaagcaatgc agacatcacg gttgaccaat atcaccgcta       300
...      caaggaagat gttgggatta tgaaggatca aaatatggat tcgtatagat tctcaatctc       360
...      ttggccaaga atactcccaa agggaaagtt gagcggaggc ataaatcacg aaggaatcaa       420
...      atattacaac aaccttatca acgaactatt ggctaacggt atacaaccat ttgtaactct       480
...      ttttcattgg gatcttcccc aagtcttaga agatgagtat ggtggtttct taaactccgg       540
...      tgtaataaat gattttcgag actatacgga tctttgcttc aaggaatttg gagatagagt       600
...      gaggtattgg agtactctaa atgagccatg ggtgtttagc aattctggat atgcactagg       660
...      aacaaatgca ccaggtcgat gttcggcctc caacgtggcc aagcctggtg attctggaac       720
...      aggaccttat atagttacac acaatcaaat tcttgctcat gcagaagctg tacatgtgta       780
...      taagactaaa taccaggcat atcaaaaggg aaagataggc ataacgttgg tatctaactg       840
...      gttaatgcca cttgatgata atagcatacc agatataaag gctgccgaga gatcacttga       900
...      cttccaattt ggattgttta tggaacaatt aacaacagga gattattcta agagcatgcg       960
...      gcgtatagtt aaaaaccgat tacctaagtt ctcaaaattc gaatcaagcc tagtgaatgg      1020
...      ttcatttgat tttattggta taaactatta ctcttctagt tatattagca atgccccttc      1080
...      acatggcaat gccaaaccca gttactcaac aaatcctatg accaatattt catttgaaaa      1140
...      acatgggata cccttaggtc caagggctgc ttcaatttgg atatatgttt atccatatat      1200
...      gtttatccaa gaggacttcg agatcttttg ttacatatta aaaataaata taacaatcct      1260
...      gcaattttca atcactgaaa atggtatgaa tgaattcaac gatgcaacac ttccagtaga      1320
...      agaagctctt ttgaatactt acagaattga ttactattac cgtcacttat actacattcg      1380
...      ttctgcaatc agggctggct caaatgtgaa gggtttttac gcatggtcat ttttggactg      1440
...      taatgaatgg tttgcaggct ttactgttcg ttttggatta aactttgtag attagaaaga      1500
...      tggattaaaa aggtacccta agctttctgc ccaatggtac aagaactttc tcaaaagaaa      1560
...      ctagctagta ttattaaaag aactttgtag tagattacag tacatcgttt gaagttgagt      1620
...      tggtgcacct aattaaataa aagaggttac tcttaacata tttttaggcc attcgttgtg      1680
...      aagttgttag gctgttattt ctattatact atgttgtagt aataagtgca ttgttgtacc      1740
...      agaagctatg atcataacta taggttgatc cttcatgtat cagtttgatg ttgagaatac      1800
...      tttgaattaa aagtcttttt ttattttttt aaaaaaaaaa aaaaaaaaaa aaaaaaaaa       1859
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

# skbio modules
from skbio.io import create_format, EMBLFormatError
from skbio.io.format._base import (_line_generator, _get_nth_sequence)
from skbio.sequence import Sequence, DNA, RNA, Protein

# look at skbio.io.registry to have an idea on how to define this class
embl = create_format('embl')

# Method to determine if file is in EMBL format or not
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
@embl.reader(DNA)
def _embl_to_dna(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_embls(fh), seq_num)
    return _construct(record, DNA, **kwargs)

# Method to read EMBL data as skbio.sequence.DNA
@embl.reader(RNA)
def _embl_to_rna(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_embls(fh), seq_num)
    return _construct(record, RNA, **kwargs)

def _construct(record, constructor=None, **kwargs):
    '''Construct the object of Sequence, DNA, RNA, or Protein.
    '''
    pass
#    seq, md, imd = record
#    if 'lowercase' not in kwargs:
#        kwargs['lowercase'] = True
#    if constructor is None:
#        unit = md['LOCUS']['unit']
#        if unit == 'bp':
#            # RNA mol type has T instead of U for genbank from from NCBI
#            constructor = DNA
#        elif unit == 'aa':
#            constructor = Protein
#
#    if constructor == RNA:
#        return DNA(
#            seq, metadata=md, interval_metadata=imd, **kwargs).transcribe()
#    else:
#        return constructor(
#            seq, metadata=md, interval_metadata=imd, **kwargs)


def _parse_embls(fh):
    """Chunck multiple EMBL records by '//', and returns a generator"""
    pass
#    data_chunks = []
#    for line in _line_generator(fh, skip_blanks=True, strip=False):
#        if line.startswith('//'):
#            yield _parse_single_genbank(data_chunks)
#            data_chunks = []
#        else:
#            data_chunks.append(line)

            
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
    
    Note 1 - Molecule type: this represents the type of molecule as stored and can
    be any value from the list of current values for the mandatory mol_type source
    qualifier. This item should be the same as the value in the mol_type
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
                         ' +([^\s]+);'    # ie: CD789012
                         ' +SV ([0-9]+);' # 4
                         ' +(\w+);'       # linear
                         ' +([^;]+);'     # genomic DNA
                         ' +(\w+);'       # HTG
                         ' +(\w+);'       # MAM
                         ' +(\d+)'        # 500
                         ' +(\w+)\.$')    # BP
    
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

# TODO: define this method
def _parse_reference(lines):
    '''Parse single REFERENCE field.
    '''
    
    pass
    
    

