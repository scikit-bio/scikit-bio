#!/usr/bin/env python
"""Tests for stockholm sequence format writer.
"""
from unittest import TestCase, main
from collections import OrderedDict

from skbio.core.alignment import Alignment
from skbio.core.sequence import DNA
from skbio.format.stockholm import Stockholm, stockholm_from_alignment


class StockholmTests(TestCase):

    """Tests for stockholm writer.
    """

    def setUp(self):
        """Setup for stockholm tests."""
        seqs = [DNA("ACC-G-GGTA", identifier="seq1"),
                DNA("TCC-G-GGCA", identifier="seq2")]
        self.aln = Alignment(seqs)
        self.GF = OrderedDict({"ID": "TEST", "SQ": 2})
        self.GS = OrderedDict({"seq1": {"AC": 111}, "seq2": {"AC": 222}})
        self.GR = OrderedDict({"seq1": {"SS": "1110101111"},
                               "seq2": {"SS": "0110101110"}})
        self.GC = OrderedDict({"SS_cons": "(((....)))"})

    def test_stockholm_from_alignment(self):
        """ Make sure stockholm with all information contained is formatted
        correctly """
        st = Stockholm(aln=self.aln, GC=self.GC, GF=self.GF, GS=self.GS,
                       GR=self.GR)
        obs = stockholm_from_alignment(st)
        exp = ("# STOCKHOLM 1.0\n#=GF SQ 2\n#=GF ID TEST\n#=GS seq2 AC 222\n"
               "#=GS seq1 AC 111\nseq1          ACC-G-GGTA\n"
               "#=GR seq1 SS  1110101111\nseq2          TCC-G-GGCA\n"
               "#=GR seq2 SS  0110101110\n#=GC SS_cons  (((....)))\n//")
        self.assertEqual(obs, exp)

    def test_stockholm_from_alignment_GC(self):
        """ Make sure stockholm with only GC information contained is formatted
        correctly """
        st = Stockholm(aln=self.aln, GC=self.GC, GF=None, GS=None,
                       GR=None)
        obs = stockholm_from_alignment(st)
        exp = ("# STOCKHOLM 1.0\nseq1          ACC-G-GGTA\n"
               "seq2          TCC-G-GGCA\n"
               "#=GC SS_cons  (((....)))\n//")
        self.assertEqual(obs, exp)

    def test_stockholm_from_alignment_GF(self):
        """ Make sure stockholm with only GF information contained is formatted
        correctly """
        st = Stockholm(aln=self.aln, GC=None, GF=self.GF, GS=None,
                       GR=None)
        obs = stockholm_from_alignment(st)
        exp = ("# STOCKHOLM 1.0\n#=GF SQ 2\n#=GF ID TEST\n"
               "seq1          ACC-G-GGTA\n"
               "seq2          TCC-G-GGCA\n//")
        self.assertEqual(obs, exp)

    def test_stockholm_from_alignment_GS(self):
        """ Make sure stockholm with only GS information contained is formatted
        correctly """
        st = Stockholm(aln=self.aln, GC=None, GF=None, GS=self.GS,
                       GR=None)
        obs = stockholm_from_alignment(st)
        exp = ("# STOCKHOLM 1.0\n#=GS seq2 AC 222\n"
               "#=GS seq1 AC 111\nseq1          ACC-G-GGTA\n"
               "seq2          TCC-G-GGCA\n//")
        self.assertEqual(obs, exp)

    def test_stockholm_from_alignment_GR(self):
        """ Make sure stockholm with only GR information contained is formatted
        correctly """
        st = Stockholm(aln=self.aln, GC=None, GF=None, GS=None,
                       GR=self.GR)
        obs = stockholm_from_alignment(st)
        exp = ("# STOCKHOLM 1.0\nseq1          ACC-G-GGTA\n"
               "#=GR seq1 SS  1110101111\nseq2          TCC-G-GGCA\n"
               "#=GR seq2 SS  0110101110\n//")
        self.assertEqual(obs, exp)

    def test_stockholm_from_alignment_refs(self):
        """ Make sure stockholm with references printed correctly"""
        GF = OrderedDict({"RT": ["TITLE1",  "TITLE2"],
                          "AC": "RF00360",
                          "BM": ["cmbuild  -F CM SEED",
                                 "cmsearch  -Z 274931 -E 1000000"],
                          "SQ": 9, "RA": ["Auth1;", "Auth2;"],
                          "RL": ["J Mol Biol", "Cell"],
                          "RM": ["11469857", "12007400"]})
        st = Stockholm(aln=self.aln, GC=None, GF=GF, GS=None,
                       GR=None)
        obs = stockholm_from_alignment(st)
        exp = ("# STOCKHOLM 1.0\n#=GF RN [1]\n#=GF RM 11469857\n#=GF RT TITLE1"
               "\n#=GF RA Auth1;\n#=GF RL J Mol Biol\n#=GF RN [2]\n"
               "#=GF RM 12007400\n#=GF RT TITLE2\n#=GF RA Auth2;\n#=GF RL Cell"
               "\n#=GF AC RF00360\n#=GF BM cmbuild  -F CM SEED\n"
               "#=GF BM cmsearch  -Z 274931 -E 1000000\n#=GF SQ 9\n"
               "seq1          ACC-G-GGTA\nseq2          TCC-G-GGCA\n//")
        self.assertEqual(obs, exp)

    def test_stockholm_from_alignment_trees(self):
        """ Make sure stockholm with trees printed correctly"""
        GF = OrderedDict({"NH": ["IMATREE", "IMATREETOO"],
                          "TN": ["Tree2", "Tree1"]})
        st = Stockholm(aln=self.aln, GC=None, GF=GF, GS=None,
                       GR=None)
        obs = stockholm_from_alignment(st)
        exp = ("# STOCKHOLM 1.0\n#=GF TN Tree2\n#=GF NH IMATREE\n#=GF TN Tree1"
               "\n#=GF NH IMATREETOO\nseq1          ACC-G-GGTA\n"
               "seq2          TCC-G-GGCA\n//")

        self.assertEqual(obs, exp)

if __name__ == "__main__":
    main()
