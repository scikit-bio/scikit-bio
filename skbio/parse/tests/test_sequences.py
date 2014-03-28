#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from __future__ import division

from skbio.parse.sequences import (MinimalRfamParser, RfamFinder,
                                   ChangedSequence, is_empty_or_html,
                                   is_rfam_header_line, is_rfam_seq_line,
                                   is_rfam_structure_line)
from skbio.core.alignment import Alignment


#from cogent.struct.rna2d import WussStructure
#from cogent.core.moltype import BYTES
#Sequence = BYTES.Sequence

from skbio.core.exception import FastqParseError, RecordError
from skbio.parse.sequences import parse_fastq
from skbio.parse.sequences import parse_fasta
from unittest import TestCase, main


class GenericFastaTest(TestCase):
    """Setup data for all the various FASTA parsers."""

    def setUp(self):
        """standard files"""
        self.labels = '>abc\n>def\n>ghi\n'.split('\n')
        self.oneseq = '>abc\nUCAG\n'.split('\n')
        self.multiline = '>xyz\nUUUU\nCC\nAAAAA\nG'.split('\n')
        self.threeseq = '>123\na\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split(
            '\n')
        self.twogood = '>123\n\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split(
            '\n')
        self.oneX = '>123\nX\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split('\n')
        self.nolabels = 'GJ>DSJGSJDF\nSFHKLDFS>jkfs\n'.split('\n')
        self.empty = []


class ParseFastaTests(GenericFastaTest):
    """Tests of parse_fasta: returns (label, seq) tuples."""

    def test_empty(self):
        """parse_fasta should return empty list from 'file' w/o labels
        """
        self.assertEqual(list(parse_fasta(self.empty)), [])
        self.assertEqual(list(parse_fasta(self.nolabels, strict=False)),
                         [])
        self.assertRaises(RecordError, list, parse_fasta(self.nolabels))

    def test_no_labels(self):
        """parse_fasta should return empty list from file w/o seqs"""
        # should fail if strict (the default)
        self.assertRaises(RecordError, list,
                          parse_fasta(self.labels, strict=True))
        # if not strict, should skip the records
        self.assertEqual(list(parse_fasta(self.labels, strict=False)),
                         [])

    def test_single(self):
        """parse_fasta should read single record as (label, seq) tuple
        """
        f = list(parse_fasta(self.oneseq))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('abc', 'UCAG'))

        f = list(parse_fasta(self.multiline))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('xyz', 'UUUUCCAAAAAG'))

    def test_gt_bracket_in_seq(self):
        """parse_fasta handles alternate finder function

            this test also illustrates how to use the parse_fasta
            to handle "sequences" that start with a > symbol, which can
            happen when we abuse the parse_fasta to parse
            fasta-like sequence quality files.
        """
        oneseq_w_gt = '>abc\n>CAG\n'.split('\n')

        def get_two_line_records(infile):
            line1 = None
            for line in infile:
                if line1 is None:
                    line1 = line
                else:
                    yield (line1, line)
                    line1 = None
        f = list(parse_fasta(oneseq_w_gt, finder=get_two_line_records))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('abc', '>CAG'))

    def test_multiple(self):
        """parse_fasta should read multiline records correctly"""
        f = list(parse_fasta(self.threeseq))
        self.assertEqual(len(f), 3)
        a, b, c = f
        self.assertEqual(a, ('123', 'a'))
        self.assertEqual(b, ('abc', 'caggac'))
        self.assertEqual(c, ('456', 'cg'))

    def test_multiple_bad(self):
        """parse_fasta should complain or skip bad records"""
        self.assertRaises(RecordError, list, parse_fasta(self.twogood))
        f = list(parse_fasta(self.twogood, strict=False))
        self.assertEqual(len(f), 2)
        a, b = f
        self.assertEqual(a, ('abc', 'caggac'))


class ParseFastqTests(TestCase):
    def setUp(self):
        """ Initialize variables to be used by the tests """
        self.FASTQ_EXAMPLE = FASTQ_EXAMPLE.split('\n')
        self.FASTQ_EXAMPLE_2 = FASTQ_EXAMPLE_2.split('\n')

    def test_parse(self):
        """sequence and info objects should correctly match"""
        for label, seq, qual in parse_fastq(self.FASTQ_EXAMPLE):
            self.assertTrue(label in DATA)
            self.assertEqual(seq, DATA[label]["seq"])
            self.assertEqual(qual, DATA[label]["qual"])

    def test_parse_error(self):
        """Does this raise a FastqParseError with incorrect input?"""
        with self.assertRaises(FastqParseError):
            list(parse_fastq(self.FASTQ_EXAMPLE_2, strict=True))

class RfamParserTests(TestCase):
    """ Tests componenets of the rfam parser, in the rfam.py file """

    def setUp(self):
        """ Construct some fake data for testing purposes """

        self._fake_headers = []
        temp = list(fake_headers.split('\n'))
        for line in temp:
            self._fake_headers.append(line.strip())
        del temp
        
        self._fake_record_no_headers =\
            list(fake_record_no_headers.split('\n'))

        self._fake_record_no_sequences =\
            list(fake_record_no_sequences.split('\n'))

        self._fake_record_no_structure =\
            list(fake_record_no_structure.split('\n'))

        self._fake_two_records =\
            list(fake_two_records.split('\n'))
            
        self._fake_record =\
            list(fake_record.split('\n'))

        self._fake_record_bad_header_1 =\
            list(fake_record_bad_header_1.split('\n'))
            
        self._fake_record_bad_header_2 =\
            list(fake_record_bad_header_2.split('\n'))

        self._fake_record_bad_sequence_1 =\
            list(fake_record_bad_sequence_1.split('\n'))

        self._fake_record_bad_structure_1 =\
            list(fake_record_bad_structure_1.split('\n'))                                                    
        self._fake_record_bad_structure_2 =\
            list(fake_record_bad_structure_2.split('\n'))

        self.single_family = single_family.split('\n')
            
    def test_is_empty_or_html(self):
        """is_empty_or_html: should ignore empty and HTML line"""
        line = '        '
        self.assertEqual(is_empty_or_html(line), True)
        line = '\n\n'
        self.assertEqual(is_empty_or_html(line), True)
        line = '<pre>'
        self.assertEqual(is_empty_or_html(line), True)
        line = '</pre>\n\n'
        self.assertEqual(is_empty_or_html(line), True)
        line = '\t<//\n'
        self.assertEqual(is_empty_or_html(line), False)

    def test_is_rfam_header_line(self):
        """is_rfam_header_line: functions correctly w/ various lines """
        self.assertEqual(is_rfam_header_line('#=GF'), True)
        self.assertEqual(is_rfam_header_line('#=GF AC   RF00001'), True)
        self.assertEqual(is_rfam_header_line('#=GF CC   until it is\
            required for transcription. '), True)

        self.assertEqual(is_rfam_header_line(''), False)
        self.assertEqual(is_rfam_header_line('X07545.1/505-619 '), False)
        self.assertEqual(is_rfam_header_line('#=G'), False)
        self.assertEqual(is_rfam_header_line('=GF'), False)
        self.assertEqual(is_rfam_header_line('#=GC SS_cons'), False)

    def test_is_rfam_seq_line(self):
        """is_rfam_seq_line: functions correctly w/ various lines """
        s = 'X07545.1/505-619                     .\
            .ACCCGGC.CAUA...GUGGCCG.GGCAA.CAC.CCGG.U.C..UCGUU'
        self.assertTrue(is_rfam_seq_line('s'))
        self.assertTrue(is_rfam_seq_line('X07545.1/505-619'))
        self.assertTrue(is_rfam_seq_line('M21086.1/8-123'))

        self.assertFalse(is_rfam_seq_line(''))
        self.assertFalse(is_rfam_seq_line('#GF='))
        self.assertFalse(is_rfam_seq_line('//blah'))

    def test_is_rfam_structure_line(self):
        """is_rfam_structure_line: functions correctly w/ various lines """
        s = '#=GC SS_cons\
            <<<<<<<<<........<<.<<<<.<...<.<...<<<<.<.<.......'
        self.assertEqual(is_rfam_structure_line(s), True)
        self.assertEqual(is_rfam_structure_line('#=GC SS_cons'), True)
        self.assertEqual(is_rfam_structure_line('#=GC SS_cons '), True)

        self.assertEqual(is_rfam_structure_line(''), False)
        self.assertEqual(is_rfam_structure_line(' '), False)
        self.assertEqual(is_rfam_structure_line('#=GF AC   RF00001'), False)
        self.assertEqual(is_rfam_structure_line('X07545.1/505-619'), False)
        self.assertEqual(is_rfam_structure_line('=GC SS_cons'), False)
        self.assertEqual(is_rfam_structure_line('#=GC'), False)
        self.assertEqual(is_rfam_structure_line('#=GC RF'), False)

    def test_MinimalRfamParser_strict_missing_fields(self):
        """MinimalRfamParser: toggle strict functions w/ missing fields"""
        # strict = True
        
        self.assertRaises(RecordError,list,\
            MinimalRfamParser(self._fake_record_no_sequences))
        
        self.assertRaises(RecordError,list,\
            MinimalRfamParser(self._fake_record_no_structure))

        # strict = False
        # no header shouldn't be a problem
        self.assertEqual(list(MinimalRfamParser(self._fake_record_no_headers,\
            strict=False)), [([],{'Z11765.1/1-89':'GGUC'},'............>>>')])
        # should get empty on missing sequence or missing structure
        self.assertEqual(list(MinimalRfamParser(self._fake_record_no_sequences,\
            strict=False)), [])
        self.assertEqual(list(MinimalRfamParser(self._fake_record_no_structure,\
            strict=False)), [])

    def test_MinimalRfamParser_strict_invalid_sequence(self):
        """MinimalRfamParser: toggle strict functions w/ invalid seq
        """
        #strict = True
        self.assertRaises(RecordError,list,\
            MinimalRfamParser(self._fake_record_bad_sequence_1))

        # strict = False
        # you expect to get back as much information as possible, also
        # half records or sequences
        result = MinimalRfamParser(self._fake_record_bad_sequence_1,strict=False)
        self.assertEqual(len(list(MinimalRfamParser(\
            self._fake_record_bad_sequence_1,strict=False))[0][1].NamedSeqs), 3)            

    def test_MinimalRfamParser_strict_invalid_structure(self):
        """MinimalRfamParser: toggle strict functions w/ invalid structure
        """
        #strict = True
        self.assertRaises(RecordError,list,\
            MinimalRfamParser(self._fake_record_bad_structure_1))

        # strict = False
        self.assertEqual(list(MinimalRfamParser(\
            self._fake_record_bad_structure_1,strict=False))[0][2],None)                                

    def test_MinimalRfamParser_w_valid_data(self):
        """MinimalRfamParser: integrity of output """

        # Some ugly constructions here, but this is what the output of
        # parsing fake_two_records should be
        headers = ['#=GF AC   RF00014','#=GF AU   Mifsud W']
        sequences =\
        {'U17136.1/898-984':\
        ''.join(['AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA',\
            'AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU']),\
        'M15749.1/155-239':\
        ''.join(['AACGCAUCGGAUUUCCCGGUGUAACGAA-UUUUCAAGUGCUUCUUGCAUU',\
            'AGCAAGUUUGAUCCCGACUCCUG-CGAGUCGGGAUUU']),\
        'AF090431.1/222-139':\
        ''.join(['CUCACAUCAGAUUUCCUGGUGUAACGAA-UUUUCAAGUGCUUCUUGCAUA',\
            'AGCAAGUUUGAUCCCGACCCGU--AGGGCCGGGAUUU'])}

        structure = WussStructure(''.join(\
        ['...<<<<<<<.....>>>>>>>....................<<<<<...',\
        '.>>>>>....<<<<<<<<<<.....>>>>>>>>>>..']))
        
        data = []
        for r in MinimalRfamParser(self._fake_two_records, strict=False):
            data.append(r)
        self.assertEqual(data[0],(headers,sequences,structure))
        assert isinstance(data[0][1],Alignment)

        # This line tests that invalid entries are ignored when strict=False
        # Note, there are two records in self._fake_two_records, but 2nd is
        # invalid
        self.assertEqual(len(data),1)            
            
    def test_RfamFinder(self):
        """RfamFinder: integrity of output """
        fake_record = ['a','//','b','b','//']
        num_records = 0
        data = []
        for r in RfamFinder(fake_record):
            data.append(r)
            num_records += 1
        self.assertEqual(num_records, 2)
        self.assertEqual(data[0], ['a','//'])
        self.assertEqual(data[1], ['b','b','//'])

    def test_ChangedSequence(self):
        """ChangedSequence: integrity of output"""
        # Made up input, based on a line that would look like:
        # U17136.1/898-984  AACA..CAU..CAGAUUUCCU..GGUGUAA.CGAA
        
        s_in = 'AACA..CAU..CAGAUUUCCU..GGUGUAA.CGAA'
        s_out = 'AACA--CAU--CAGAUUUCCU--GGUGUAA-CGAA'
        sequence = ChangedSequence(s_in)
        
        self.assertEqual(sequence, s_out)

        # test some extremes on the seq
        # sequence of all blanks
        s_in = '.' * 5
        s_out = '-' * 5
        sequence = ChangedSequence(s_in)

        self.assertEqual(sequence, s_out)

        # sequence of no blanks
        s_in = 'U' * 5
        s_out = 'U' * 5
        sequence = ChangedSequence(s_in)

        self.assertEqual(sequence, s_out)
       

# This is an altered version of some header info from Rfam.seed modified to
# incorporate different cases for testing
fake_headers = """#=GF AC   RF00001
#=GF AU   Griffiths-Jones SR
#=GF ID   5S_rRNA
#=GF RT   5S Ribosomal RNA Database.
#=GF DR   URL; http://oberon.fvms.ugent.be:8080/rRNA/ssu/index.html;
#=GF DR   URL; http://rdp.cme.msu.edu/html/;
#=GF CC   This is a short
#=GF CC   comment
#=GF SQ   606
#=GF PK   not real"""

fake_record_no_headers ="""Z11765.1/1-89                        GGUC
#=GC SS_cons                         ............>>>
//"""

fake_record_no_sequences ="""#=GF AC   RF00006
#=GC SS_cons                         ............>
//"""

fake_record_no_structure ="""#=GF AC   RF00006

Z11765.1/1-89                        GGUCAGC
//"""

fake_two_records ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AU   Mifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//
#=GF AC   RF00015
//"""

fake_record ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AU   Mifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

fake_record_bad_header_1 ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AUMifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

fake_record_bad_header_2 ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GFAUMifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

fake_record_bad_sequence_1 ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AU   Mifsud W

U17136.1/898-984AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

fake_record_bad_structure_1 ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AU   Mifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons...<<<<<<<.....>>>>>>>....................<<<<<...
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

fake_record_bad_structure_2 ="""# STOCKHOLM 1.0

#=GF AC   RF00014
#=GF AU   Mifsud W

U17136.1/898-984               AACACAUCAGAUUUCCUGGUGUAACGAAUUUUUUAAGUGCUUCUUGCUUA
M15749.1/155-239               AACGCAUCGGAUUUCCCGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUU
AF090431.1/222-139             CUCACAUCAGAUUUCCUGGUGUAACGAA.UUUUCAAGUGCUUCUUGCAUA
#=GC SS_cons                   ...<<<<<<<.....>>>>>>>....................<<<<<!!!
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

U17136.1/898-984               AGCAAGUUUCAUCCCGACCCCCUCAGGGUCGGGAUUU
M15749.1/155-239               AGCAAGUUUGAUCCCGACUCCUG.CGAGUCGGGAUUU
AF090431.1/222-139             AGCAAGUUUGAUCCCGACCCGU..AGGGCCGGGAUUU
#=GC SS_cons                   .>>>>>....<<<<<<<<<<.....>>>>>>>>>>..
#=GC RF                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//"""

single_family=\
"""K02120.1/628-682      AUGGGAAAUUCCCCCUCCUAUAACCCCCCCGCUGGUAUCUCCCCCUCAGA
D00647.1/629-683      AUGGGAAACUCCCCCUCCUAUAACCCCCCCGCUGGCAUCUCCCCCUCAGA
#=GC SS_cons          <<<<<<.........>>>>>>.........<<<<<<.............>

K02120.1/628-682      CUGGC
D00647.1/629-683      CUGGC
#=GC SS_cons          >>>>>
//"""

DATA = {
    "GAPC_0015:6:1:1259:10413#0/1":
    dict(seq='AACACCAAACTTCTCCACCACGTGAGCTACAAAAG',
         qual=r'````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF'),
    "GAPC_0015:6:1:1283:11957#0/1":
    dict(seq='TATGTATATATAACATATACATATATACATACATA',
         qual=r']KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb'),
    "GAPC_0015:6:1:1284:10484#0/1":
    dict(seq='TCAGTTTTCCTCGCCATATTTCACGTCCTAAAGCG',
         qual=r'UM_]]U_]Z_Y^\^^``Y]`^SZ]\Ybb`^_LbL_'),
    "GAPC_0015:6:1:1287:17135#0/1":
    dict(seq='TGTGCCTATGGAAGCAGTTCTAGGATCCCCTAGAA',
         qual=r'^aacccL\ccc\c\cTKTS]KZ\]]I\[Wa^T`^K'),
    "GAPC_0015:6:1:1293:3171#0/1":
    dict(seq="AAAGAAAGGAAGAAAAGAAAAAGAAACCCGAGTTA",
             qual=r"b`bbbU_[YYcadcda_LbaaabWbaacYcc`a^c"),
    "GAPC_0015:6:1:1297:10729#0/1":
    dict(seq="TAATGCCAAAGAAATATTTCCAAACTACATGCTTA",
             qual=r"T\ccLbb``bacc]_cacccccLccc\ccTccYL^"),
    "GAPC_0015:6:1:1299:5940#0/1":
    dict(seq="AATCAAGAAATGAAGATTTATGTATGTGAAGAATA",
             qual=r"dcddbcfffdfffd`dd`^`c`Oc`Ybb`^eecde"),
    "GAPC_0015:6:1:1308:6996#0/1":
    dict(seq="TGGGACACATGTCCATGCTGTGGTTTTAACCGGCA",
             qual=r"a]`aLY`Y^^ccYa`^^TccK_X]\c\c`caTTTc"),
    "GAPC_0015:6:1:1314:13295#0/1":
    dict(seq="AATATTGCTTTGTCTGAACGATAGTGCTCTTTGAT",
         qual=r"cLcc\\dddddaaYd`T```bLYT\`a```bZccc"),
    "GAPC_0015:6:1:1317:3403#0/1":
    dict(seq="TTGTTTCCACTTGGTTGATTTCACCCCTGAGTTTG",
         # had to add space in qual line
         qual=r"\\\ZTYTSaLbb``\_UZ_bbcc`cc^[ac\a\Tc ".strip())
}

FASTQ_EXAMPLE = r"""@GAPC_0015:6:1:1259:10413#0/1
AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
+GAPC_0015:6:1:1259:10413#0/1
````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
@GAPC_0015:6:1:1283:11957#0/1
TATGTATATATAACATATACATATATACATACATA
+GAPC_0015:6:1:1283:11957#0/1
]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb
@GAPC_0015:6:1:1284:10484#0/1
TCAGTTTTCCTCGCCATATTTCACGTCCTAAAGCG
+GAPC_0015:6:1:1284:10484#0/1
UM_]]U_]Z_Y^\^^``Y]`^SZ]\Ybb`^_LbL_
@GAPC_0015:6:1:1287:17135#0/1
TGTGCCTATGGAAGCAGTTCTAGGATCCCCTAGAA
+GAPC_0015:6:1:1287:17135#0/1
^aacccL\ccc\c\cTKTS]KZ\]]I\[Wa^T`^K
@GAPC_0015:6:1:1293:3171#0/1
AAAGAAAGGAAGAAAAGAAAAAGAAACCCGAGTTA
+GAPC_0015:6:1:1293:3171#0/1
b`bbbU_[YYcadcda_LbaaabWbaacYcc`a^c
@GAPC_0015:6:1:1297:10729#0/1
TAATGCCAAAGAAATATTTCCAAACTACATGCTTA
+GAPC_0015:6:1:1297:10729#0/1
T\ccLbb``bacc]_cacccccLccc\ccTccYL^
@GAPC_0015:6:1:1299:5940#0/1
AATCAAGAAATGAAGATTTATGTATGTGAAGAATA
+GAPC_0015:6:1:1299:5940#0/1
dcddbcfffdfffd`dd`^`c`Oc`Ybb`^eecde
@GAPC_0015:6:1:1308:6996#0/1
TGGGACACATGTCCATGCTGTGGTTTTAACCGGCA
+GAPC_0015:6:1:1308:6996#0/1
a]`aLY`Y^^ccYa`^^TccK_X]\c\c`caTTTc
@GAPC_0015:6:1:1314:13295#0/1
AATATTGCTTTGTCTGAACGATAGTGCTCTTTGAT
+GAPC_0015:6:1:1314:13295#0/1
cLcc\\dddddaaYd`T```bLYT\`a```bZccc
@GAPC_0015:6:1:1317:3403#0/1
TTGTTTCCACTTGGTTGATTTCACCCCTGAGTTTG
+GAPC_0015:6:1:1317:3403#0/1
\\\ZTYTSaLbb``\_UZ_bbcc`cc^[ac\a\Tc"""

FASTQ_EXAMPLE_2 = r"""@GAPC_0017:6:1:1259:10413#0/1
AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
+GAPC_0015:6:1:1259:10413#0/1
````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
@GAPC_0015:6:1:1283:11957#0/1
TATGTATATATAACATATACATATATACATACATA
+GAPC_0015:6:1:1283:11957#0/1
]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb
@GAPC_0015:6:1:1284:10484#0/1
"""

if __name__ == "__main__":
    main()
