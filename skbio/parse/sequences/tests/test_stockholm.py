#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import division
from StringIO import StringIO

from skbio.parse.sequences.stockholm import (parse_stockholm, _parse_gf_info,
                                             _parse_gc_info, _parse_gs_gr_info)
from skbio.format.stockholm import Stockholm
from skbio.core.alignment import Alignment
from skbio.core.sequence import DNA
from skbio.core.exception import StockholmParseError

from unittest import TestCase, main


class ParseStockholmTests(TestCase):
    def setUp(self):
        """Setup for stockholm tests."""
        seqs = [DNA("ACC-G-GGTA", identifier="seq1"),
                DNA("TCC-G-GGCA", identifier="seq2")]
        self.aln = Alignment(seqs)
        self.GF = {"RT": ["TITLE1",  "TITLE2"],
                   "AC": "RF00360",
                   "BM": ["cmbuild  -F CM SEED",
                          "cmsearch  -Z 274931 -E 1000000"],
                   "SQ": "9", "RA": ["Auth1;", "Auth2;"],
                   "RL": ["J Mol Biol", "Cell"],
                   "RM": ["11469857", "12007400"],
                   'RN': ['[1]', '[2]']}
        self.GS = {"seq1": {"AC": "111"}, "seq2": {"AC": "222"}}
        self.GR = {"seq1": {"SS": "1110101111"},
                   "seq2": {"SS": "0110101110"}}
        self.GC = {"SS_cons": "(((....)))"}

    """Tests of parse_stockholm: returns Stockholm named tuples."""

    def test_parse_stockholm_alignment(self):
        """make sure can parse basic sto file with interleaved alignment"""
        sto = StringIO("# STOCKHOLM 1.0\n"
                       "seq1      ACC-G\n"
                       "seq2      TCC-G\n\n"
                       "seq1      -GGTA\n"
                       "seq2      -GGCA\n//")
        obs_sto = parse_stockholm(sto, DNA).next()
        exp_sto = Stockholm(self.aln, {}, {}, {}, {})
        self.assertEqual(obs_sto, exp_sto)

    def test_parse_stockholm_GF(self):
        """Make sure GF lines are parsed correctly"""
        sto = StringIO("# STOCKHOLM 1.0\n#=GF RN [1]\n#=GF RM 11469857\n"
                       "#=GF RT TITLE1\n#=GF RA Auth1;\n#=GF RL J Mol Biol\n"
                       "#=GF RN [2]\n#=GF RM 12007400\n#=GF RT TITLE2\n"
                       "#=GF RA Auth2;\n#=GF RL Cell\n#=GF AC RF00360\n"
                       "#=GF BM cmbuild  -F CM SEED\n"
                       "#=GF BM cmsearch  -Z 274931 -E 1000000\n#=GF SQ 9\n"
                       "seq1         ACC-G-GGTA\nseq2         TCC-G-GGCA\n//")
        obs_sto = parse_stockholm(sto, DNA).next()
        exp_sto = Stockholm(self.aln, self.GF, {}, {}, {})
        self.assertEqual(obs_sto, exp_sto)

    def test_parse_stockholm_GC(self):
        """Make sure GC lines are parsed correctly"""
        sto = StringIO("# STOCKHOLM 1.0\n"
                       "seq1         ACC-G-GGTA\nseq2         TCC-G-GGCA\n"
                       "#=GC SS_cons (((....)))\n//")
        obs_sto = parse_stockholm(sto, DNA).next()
        exp_sto = Stockholm(self.aln, {}, {}, {}, self.GC)
        self.assertEqual(obs_sto, exp_sto)

    def test_parse_stockholm_GS(self):
        """Make sure GS lines are parsed correctly"""
        sto = StringIO("# STOCKHOLM 1.0\n#=GS seq2 AC 222\n#=GS seq1 AC 111\n"
                       "seq1          ACC-G-GGTA\n"
                       "seq2          TCC-G-GGCA\n//")
        obs_sto = parse_stockholm(sto, DNA).next()
        exp_sto = Stockholm(self.aln, {}, self.GS, {}, {})
        self.assertEqual(obs_sto, exp_sto)

    def test_parse_stockholm_GR(self):
        """Make sure GR lines are parsed correctly"""
        sto = StringIO("# STOCKHOLM 1.0\nseq1          ACC-G\n"
                       "#=GR seq1 SS  11101\nseq2          TCC-G\n"
                       "#=GR seq2 SS  01101\n\nseq1          -GGTA\n"
                       "#=GR seq1 SS  01111\nseq2          -GGCA\n"
                       "#=GR seq2 SS  01110\n//")
        obs_sto = parse_stockholm(sto, DNA).next()
        exp_sto = Stockholm(self.aln, {}, {}, self.GR, {})
        self.assertEqual(obs_sto, exp_sto)

    def test_parse_stockholm_multi(self):
        """Make sure yield works correctly with multi-alignment sto files"""
        sto = StringIO("# STOCKHOLM 1.0\n#=GS seq2 AC 222\n#=GS seq1 AC 111\n"
                       "seq1          ACC-G-GGTA\n"
                       "seq2          TCC-G-GGCA\n//\n"
                       "# STOCKHOLM 1.0\nseq1          ACC-G-GGTA\n"
                       "#=GR seq1 SS  1110101111\nseq2          TCC-G-GGCA\n"
                       "#=GR seq2 SS  0110101110\n//")
        obs_sto = parse_stockholm(sto, DNA)
        count = 0
        for obs in obs_sto:
            if count == 0:
                exp_sto = Stockholm(self.aln, {}, self.GS, {}, {})
                self.assertEqual(obs, exp_sto)
            elif count == 1:
                exp_sto = Stockholm(self.aln, {}, {}, self.GR, {})
                self.assertEqual(obs, exp_sto)
            else:
                raise AssertionError("More than 2 sto alignments parsed!")
            count += 1

    def test_parse_gf_info_nongf(self):
        """Makes sure error raised if non-GF line passed"""
        sto = ["#=GF AC BLAAAAAAAHHH", "#=GC HUH THIS SHOULD NOT BE HERE"]
        self.assertRaises(StockholmParseError, _parse_gf_info, sto)

    def test_parse_gf_info_malformed(self):
        """Makes sure error raised if too short a line passed"""
        sto = ["#=GF AC BLAAAAAAAHHH", "#=GF SMALL"]
        self.assertRaises(StockholmParseError, _parse_gc_info, sto)

    def test_parse_gc_info_nongf(self):
        """Makes sure error raised if non-GC line passed"""
        sto = ["#=GC AC BLAAAAAAAHHH", "#=GF HUH THIS SHOULD NOT BE HERE"]
        self.assertRaises(StockholmParseError, _parse_gf_info, sto)

    def test_parse_gc_info_strict(self):
        """Make sure error raised if GC lines bad length and strict parsing"""
        sto = ["#=GC SS_cons (((..)))"]
        self.assertRaises(StockholmParseError, _parse_gc_info, sto, seqlen=20,
                          strict=True)

    def test_parse_gc_info_malformed(self):
        """Makes sure error raised if too short a line passed"""
        sto = ["#=GC AC BLAAAAAAAHHH", "#=GC SMALL"]
        self.assertRaises(StockholmParseError, _parse_gc_info, sto)

    def test_parse_gs_gr_info_mixed(self):
        """Makes sure error raised if mixed GS and GR lines passed"""
        sto = ["#=GS seq1 AC BLAAA", "#=GR seq2 HUH THIS SHOULD NOT BE HERE"]
        self.assertRaises(StockholmParseError, _parse_gs_gr_info, sto)

    def test_parse_gs_gr_info_malformed(self):
        """Makes sure error raised if too short a line passed"""
        sto = ["#=GS AC BLAAAAAAAHHH", "#=GS SMALL"]
        self.assertRaises(StockholmParseError, _parse_gs_gr_info, sto)

    def test_parse_gs_gr_info_strict(self):
        """Make sure error raised if GR lines bad length and strict parsing"""
        sto = ["#=GR seq1 SS  10101111", "#=GR seq2 SS  01101"]
        self.assertRaises(StockholmParseError, _parse_gs_gr_info, sto,
                          seqlen=20, strict=True)

if __name__ == "__main__":
    main()
