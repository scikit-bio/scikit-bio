#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio import DNA, RNA, Protein
from skbio.sequence import (GeneticCode, genetic_code,
                            GeneticCodeInitError, InvalidCodonError)


class GeneticCodeTests(TestCase):

    """Tests of the GeneticCode class."""

    def setUp(self):
        """Set up some standard genetic code representations."""
        self.sgc = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAAD"
                    "DEEGGGG")
        self.mt = ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADD"
                   "EEGGGG")
        self.allg = ("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                     "GGGGGGGG")

        self.wrong_length = [
            "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
            "",
            "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
            "G",
        ]
        self.ncbi_standard = [
            'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
            1,
            'Standard Nuclear',
            '---M---------------M---------------M----------------------------',
        ]

    def test_init(self):
        """GeneticCode init should work with correct-length sequences"""
        sgc = GeneticCode(self.sgc)
        self.assertEqual(sgc['UUU'], 'F')
        mt = GeneticCode(self.mt)
        self.assertEqual(mt['UUU'], 'F')
        allg = GeneticCode(self.allg)
        self.assertEqual(allg['UUU'], 'G')
        for i in self.wrong_length:
            self.assertRaises(GeneticCodeInitError, GeneticCode, i)

    def test_eq(self):
        gc_1 = GeneticCode(self.sgc)
        gc_2 = GeneticCode(self.sgc)
        self.assertEqual(gc_1, gc_2)

    def test_eq_type_mismatch(self):
        self.assertFalse(GeneticCode(self.sgc) == 'i cracked the code!')

    def test_ne(self):
        gc_1 = GeneticCode(self.sgc)
        gc_2 = GeneticCode(self.sgc)
        # Explicitly using !=
        self.assertFalse(gc_1 != gc_2)

    def test_standard_code(self):
        """Standard genetic code from NCBI should have correct properties"""
        sgc = GeneticCode(*self.ncbi_standard)
        self.assertEqual(sgc.code_sequence, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRI'
                         'IIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
        self.assertEqual(sgc.start_codon_sequence, '---M---------------M------'
                         '---------M----------------------------')
        self.assertEqual(sgc.start_codons, {'TTG': 'M', 'CTG': 'M',
                                            'ATG': 'M'})
        self.assertEqual(sgc.id, 1)
        self.assertEqual(sgc.name, 'Standard Nuclear')
        self.assertEqual(sgc['UUU'], 'F')
        self.assertEqual(sgc.is_start('ATG'), True)
        self.assertEqual(sgc.is_start('AAA'), False)
        self.assertEqual(sgc.is_stop('UAA'), True)
        self.assertEqual(sgc.is_stop('AAA'), False)
        self.assertEqual(len(sgc.sense_codons), 61)
        self.assertTrue('AAA' in sgc.sense_codons)
        self.assertFalse('TGA' in sgc.sense_codons)

    def test_standard_code_lookup(self):
        """genetic_code should hold codes keyed by id as string and number"""
        sgc_new = GeneticCode(*self.ncbi_standard)
        sgc_number = genetic_code(1)
        sgc_string = genetic_code('1')
        sgc_empty = genetic_code()
        for sgc in sgc_new, sgc_number, sgc_string, sgc_empty:
            self.assertEqual(sgc.code_sequence, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQR'
                             'RRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
            self.assertEqual(sgc.start_codon_sequence, '---M---------------M--'
                             '-------------M----------------------------')
            self.assertEqual(
                sgc.start_codons, {'TTG': 'M', 'CTG': 'M', 'ATG': 'M'})
            self.assertEqual(sgc.id, 1)
            self.assertEqual(sgc.name, 'Standard Nuclear')
            self.assertEqual(sgc['TTT'], 'F')
            self.assertEqual(sgc.is_start('ATG'), True)
            self.assertEqual(sgc.is_start('AAA'), False)
            self.assertEqual(sgc.is_stop('TAA'), True)
            self.assertEqual(sgc.is_stop('AAA'), False)

        mtgc = genetic_code(2)
        self.assertEqual(mtgc.name, 'Vertebrate Mitochondrial')
        self.assertEqual(mtgc.is_start('AUU'), True)
        self.assertEqual(mtgc.is_stop('UGA'), False)

        self.assertEqual(sgc_new.changes(mtgc), {'AGA': 'R*', 'AGG': 'R*',
                                                 'ATA': 'IM', 'TGA': '*W'})
        self.assertEqual(mtgc.changes(sgc_new), {'AGA': '*R', 'AGG': '*R',
                                                 'ATA': 'MI', 'TGA': 'W*'})
        self.assertEqual(mtgc.changes(mtgc), {})
        self.assertEqual(mtgc.changes('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTT'
                         'TNNKKSSRRVVVVAAAADDEEGGGG'), {'AGA': '*R',
                         'AGG': '*R', 'ATA': 'MI', 'TGA': 'W*'})

    def test_str(self):
        """GeneticCode str() should return its code string"""
        code_strings = self.sgc, self.mt, self.allg
        codes = map(GeneticCode, code_strings)
        for code, string in zip(codes, code_strings):
            self.assertEqual(str(code), string)
        # check an example directly in case strings are bad
        self.assertEqual(str(self.sgc), "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMT"
                         "TTTNNKKSSRRVVVVAAAADDEEGGGG")

    def test_cmp(self):
        """GeneticCode cmp() should act on code strings"""
        sgc_1 = GeneticCode(self.sgc)
        sgc_2 = GeneticCode(self.sgc)
        self.assertEqual(sgc_1 is sgc_2, False)  # ensure different objects
        # self.assertNotEqual(sgc_1, sgc_2) # GREG
        self.assertEqual(sgc_1, sgc_2)
        mtgc = GeneticCode(self.mt)
        self.assertNotEqual(sgc_1, mtgc)

    def test_getitem_codon(self):
        """GeneticCode getitem should return amino acid for codon"""
        # specific checks of a particular codon in the standard code
        variant_codons = ['AUU', 'AUU', 'AUU', 'ATT', 'ATU', 'ATU']
        sgc = GeneticCode(self.sgc)
        for i in variant_codons:
            self.assertEqual(sgc[i], 'I')
        # full check for the standard code
        codons = [a + b + c for a in 'UCAG' for b in 'TCAG' for c in 'UCAG']
        for codon, aa in zip(codons, self.sgc):
            self.assertEqual(sgc[codon], aa)
        # full check for another code
        allg = GeneticCode(self.allg)
        for codon, aa in zip(codons, self.allg):
            self.assertEqual(allg[codon], aa)
        # check that degenerate codon returns X
        self.assertEqual(sgc['NNN'], 'X')

    def test_getitem_aa(self):
        """GeneticCode getitem should return codon set for aa"""
        # for all G, should return all the codons (in some order)
        allg = GeneticCode(self.allg)
        codons = [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
        g_codons = allg['G']
        codons_copy = codons[:]
        self.assertEqual(g_codons, codons_copy)

        # check some known cases in the standard genetic code
        sgc = GeneticCode(self.sgc)
        exp_ile = ['ATT', 'ATC', 'ATA']
        obs_ile = sgc['I']
        self.assertEqual(obs_ile, exp_ile)

        exp_arg = ['AGA', 'AGG', 'CGT', 'CGC', 'CGA', 'CGG']
        obs_arg = sgc['R']
        if hasattr(self, 'assertItemsEqual'):
            self.assertItemsEqual(obs_arg, exp_arg)
        else:
            self.assertCountEqual(obs_arg, exp_arg)

        exp_leu = ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']
        obs_leu = sgc['L']
        self.assertEqual(obs_leu, exp_leu)

        exp_met = ['ATG']
        obs_met = sgc['M']
        self.assertEqual(obs_met, exp_met)

        # unknown aa should return []
        self.assertEqual(sgc['U'], [])

    def test_getitem_invalid_length(self):
        """GeneticCode getitem raises InvalidCodonError on wrong length"""
        sgc = GeneticCode(self.sgc)
        self.assertRaises(InvalidCodonError, sgc.__getitem__, 'AAAA')
        self.assertRaises(InvalidCodonError, sgc.__getitem__, 'AA')

    def test_blocks(self):
        """GeneticCode blocks should return correct list"""
        sgc = GeneticCode(self.sgc)
        exp_blocks = [
            ['TTT', 'TTC', ],
            ['TTA', 'TTG', ],
            ['TCT', 'TCC', 'TCA', 'TCG'],
            ['TAT', 'TAC'],
            ['TAA', 'TAG'],
            ['TGT', 'TGC'],
            ['TGA'],
            ['TGG'],
            ['CTT', 'CTC', 'CTA', 'CTG'],
            ['CCT', 'CCC', 'CCA', 'CCG'],
            ['CAT', 'CAC'],
            ['CAA', 'CAG'],
            ['CGT', 'CGC', 'CGA', 'CGG'],
            ['ATT', 'ATC'],
            ['ATA', ],
            ['ATG', ],
            ['ACT', 'ACC', 'ACA', 'ACG'],
            ['AAT', 'AAC'],
            ['AAA', 'AAG'],
            ['AGT', 'AGC'],
            ['AGA', 'AGG'],
            ['GTT', 'GTC', 'GTA', 'GTG'],
            ['GCT', 'GCC', 'GCA', 'GCG'],
            ['GAT', 'GAC'],
            ['GAA', 'GAG'],
            ['GGT', 'GGC', 'GGA', 'GGG'],
        ]
        self.assertEqual(sgc.blocks, exp_blocks)

    def test_anticodons(self):
        """GeneticCode anticodons should return correct list"""
        sgc = GeneticCode(self.sgc)
        exp_anticodons = {
            'F': ['AAA', 'GAA', ],
            'L': ['TAA', 'CAA', 'AAG', 'GAG', 'TAG', 'CAG'],
            'Y': ['ATA', 'GTA'],
            '*': ['TTA', 'CTA', 'TCA'],
            'C': ['ACA', 'GCA'],
            'W': ['CCA'],
            'S': ['AGA', 'GGA', 'TGA', 'CGA', 'ACT', 'GCT'],
            'P': ['AGG', 'GGG', 'TGG', 'CGG'],
            'H': ['ATG', 'GTG'],
            'Q': ['TTG', 'CTG'],
            'R': ['ACG', 'GCG', 'TCG', 'CCG', 'TCT', 'CCT'],
            'I': ['AAT', 'GAT', 'TAT'],
            'M': ['CAT', ],
            'T': ['AGT', 'GGT', 'TGT', 'CGT'],
            'N': ['ATT', 'GTT'],
            'K': ['TTT', 'CTT'],
            'V': ['AAC', 'GAC', 'TAC', 'CAC'],
            'A': ['AGC', 'GGC', 'TGC', 'CGC'],
            'D': ['ATC', 'GTC'],
            'E': ['TTC', 'CTC'],
            'G': ['ACC', 'GCC', 'TCC', 'CCC'],
        }
        self.assertEqual(sgc.anticodons, exp_anticodons)

    def test_translate(self):
        """GeneticCode translate should return correct amino acid string"""
        allg = GeneticCode(self.allg)
        sgc = GeneticCode(self.sgc)
        mt = GeneticCode(self.mt)

        seq = 'AUGCAUGACUUUUGA'
        #      .  .  .  .  .        markers for codon start
        self.assertEqual(allg.translate(seq), Protein('GGGGG'))
        self.assertEqual(allg.translate(seq, 1), Protein('GGGG'))
        self.assertEqual(allg.translate(seq, 2), Protein('GGGG'))
        self.assertEqual(allg.translate(seq, 3), Protein('GGGG'))
        self.assertEqual(allg.translate(seq, 4), Protein('GGG'))
        self.assertEqual(allg.translate(seq, 12), Protein('G'))
        self.assertEqual(allg.translate(seq, 14), Protein(''))
        self.assertRaises(ValueError, allg.translate, seq, 15)
        self.assertRaises(ValueError, allg.translate, seq, 20)

        self.assertEqual(sgc.translate(seq), Protein('MHDF*'))
        self.assertEqual(sgc.translate(seq, 3), Protein('HDF*'))
        self.assertEqual(sgc.translate(seq, 6), Protein('DF*'))
        self.assertEqual(sgc.translate(seq, 9), Protein('F*'))
        self.assertEqual(sgc.translate(seq, 12), Protein('*'))
        self.assertEqual(sgc.translate(seq, 14), Protein(''))
        # check shortest translatable sequences
        self.assertEqual(sgc.translate('AAA'), Protein('K'))
        self.assertEqual(sgc.translate(''), Protein(''))

        # check that different code gives different results
        self.assertEqual(mt.translate(seq), Protein('MHDFW'))

        # check translation with invalid codon(s)
        self.assertEqual(sgc.translate('AAANNNCNC123UUU'), Protein('KXXXF'))

    def test_translate_six_frames(self):
        """GeneticCode translate_six_frames provides six-frame translation"""

        class fake_rna(str):

            """Fake RNA class with reverse-complement"""
            def __new__(cls, seq, rev):
                return str.__new__(cls, seq)

            def __init__(self, seq, rev):
                self.seq = seq
                self.rev = rev

            def rc(self):
                return self.rev

        test_rna = fake_rna('AUGCUAACAUAAA', 'UUUAUGUUAGCAU')
        #                    .  .  .  .  .    .  .  .  .  .
        sgc = GeneticCode(self.sgc)
        self.assertEqual(sgc.translate_six_frames(test_rna), [
            Protein('MLT*'), Protein('C*HK'), Protein('ANI'), Protein('FMLA'),
            Protein('LC*H'), Protein('YVS')])

        # should also actually work with an RNA or DNA sequence!!!
        test_rna = RNA('AUGCUAACAUAAA')
        self.assertEqual(sgc.translate_six_frames(test_rna), [
            Protein('MLT*'), Protein('C*HK'), Protein('ANI'), Protein('FMLA'),
            Protein('LC*H'), Protein('YVS')])

    def test_stop_indexes(self):
        """should return stop codon indexes for a specified frame"""
        sgc = GeneticCode(self.sgc)
        seq = DNA('ATGCTAACATAAA')
        expected = [[9], [4], []]
        for frame, expect in enumerate(expected):
            got = sgc.get_stop_indices(seq, start=frame)
            self.assertEqual(got, expect)

    def test_synonyms(self):
        """GeneticCode synonyms should return aa -> codon set mapping."""
        expected_synonyms = {
            'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'C': ['TGT', 'TGC'],
            'D': ['GAT', 'GAC'],
            'E': ['GAA', 'GAG'],
            'F': ['TTT', 'TTC'],
            'G': ['GGT', 'GGC', 'GGA', 'GGG'],
            'H': ['CAT', 'CAC'],
            'I': ['ATT', 'ATC', 'ATA'],
            'K': ['AAA', 'AAG'],
            'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
            'M': ['ATG'],
            'N': ['AAT', 'AAC'],
            'P': ['CCT', 'CCC', 'CCA', 'CCG'],
            'Q': ['CAA', 'CAG'],
            'R': ['AGA', 'AGG', 'CGT', 'CGC', 'CGA', 'CGG'],
            'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            'T': ['ACT', 'ACC', 'ACA', 'ACG'],
            'V': ['GTT', 'GTC', 'GTA', 'GTG'],
            'W': ['TGG'],
            'Y': ['TAT', 'TAC'],
            '*': ['TAA', 'TAG', 'TGA'],
        }
        obs_synonyms = GeneticCode(self.sgc).synonyms
        # note that the lists will be arbitrary-order
        for i in expected_synonyms:
            if hasattr(self, 'assertItemsEqual'):
                self.assertItemsEqual(obs_synonyms[i], expected_synonyms[i])
            else:
                self.assertCountEqual(obs_synonyms[i], expected_synonyms[i])

    def test_genetic_code_with_too_many_args(self):
        with self.assertRaises(TypeError):
            genetic_code(1, 2)

    def test_genetic_code_with_invalid_id(self):
        with self.assertRaises(ValueError):
            genetic_code(30)


if __name__ == '__main__':
    main()
