#!/usr/bin/env python
"""Tests for stockholm sequence format writer.
"""
from unittest import TestCase, main
from skbio.core.alignment import Alignment
from skbio.core.sequence import DNA
from skbio.format.stockholm import Stockholm, stockholm_from_alignment


class FastaTests(TestCase):

    """Tests for stockholm writer.
    """

    def setUp(self):
        """Setup for stockholm tests."""
        seqs = [DNA("ACC-G-GGTA", identifier="seq1"),
                DNA("TCC-G-GGCA", identifier="seq2")]
        self.aln = Alignment(seqs)
        self.GF = {"ID": "TEST", "SQ": 2}
        self.GS = {"seq1": {"AC": 111}, "seq2": {"AC": 222}}
        self.GR = {"seq1": {"SS": "1110101111"}, "seq2": {"SS": "0110101110"}}
        self.GC = {"SS_cons": "(((....)))"}

    def test_stockholm_from_alignment(self):
        st = Stockholm(aln=self.aln, GC=self.GC, GF=self.GF, GS=self.GS,
                       GR=self.GR)
        print stockholm_from_alignment(st)


if __name__ == "__main__":
    main()
