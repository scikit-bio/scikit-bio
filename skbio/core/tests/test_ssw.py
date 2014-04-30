#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2014--, biocore team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from unittest import TestCase, main
from skbio.core.ssw import StripedSmithWaterman, striped_smith_waterman_alignment

class TestStripedSmithWaterman(TestCase):
    align_attributes = [
        "optimal_alignment_score","suboptimal_alignment_score",
        "target_begin","target_end_optimal","target_end_suboptimal",
        "query_begin","query_end","cigar","query_sequence",
        "target_sequence"
    ]
    def _check_alignment(self, alignment, expected):
        
        for attribute in self.align_attributes:
            self.assertEqual((alignment['target_sequence'], alignment[attribute]), 
                             (expected['target_sequence'], expected[attribute]))
                              # The first element of this tuple is to identify the
                              # broken sequence if one should fail

    def test_access(self):
        q_seq = "ACTAGACTCTCTCGAGCGCGCTATATATCG"
        query = StripedSmithWaterman(q_seq)
        alignment = query("ACCTACGCGAGCATCACTAGACTCTCTCGAGCGCGCTATAAAACTCAGCTCAC")
        for accessible in align_attributes:
            self.assertEqual(getattr(alignment, accessible), alignment[accessible])

    def test_many_targets(self):
        q_seq = ""
        expected_alignments = [
            {"optimal_alignment_score":1,
            "suboptimal_alignment_score":1,
            "target_begin":1,
            "target_end_optimal":1,
            "target_end_suboptimal":1,
            "query_begin":1,
            "query_end":1,
            "cigar":1,
            "query_sequence":q_seq,
            "target_sequence":1},
            {"optimal_alignment_score":1,
            "suboptimal_alignment_score":1,
            "target_begin":1,
            "target_end_optimal":1,
            "target_end_suboptimal":1,
            "query_begin":1,
            "query_end":1,
            "cigar":1,
            "query_sequence":q_seq,
            "target_sequence":1},
            {"optimal_alignment_score":1,
            "suboptimal_alignment_score":1,
            "target_begin":1,
            "target_end_optimal":1,
            "target_end_suboptimal":1,
            "query_begin":1,
            "query_end":1,
            "cigar":1,
            "query_sequence":q_seq,
            "target_sequence":1},
            {"optimal_alignment_score":1,
            "suboptimal_alignment_score":1,
            "target_begin":1,
            "target_end_optimal":1,
            "target_end_suboptimal":1,
            "query_begin":1,
            "query_end":1,
            "cigar":1,
            "query_sequence":q_seq,
            "target_sequence":1},
            {"optimal_alignment_score":1,
            "suboptimal_alignment_score":1,
            "target_begin":1,
            "target_end_optimal":1,
            "target_end_suboptimal":1,
            "query_begin":1,
            "query_end":1,
            "cigar":1,
            "query_sequence":q_seq,
            "target_sequence":1}
        ]      
        query = StripedSmithWaterman(q_seq)
        results = []
        for expected in expected_alignments:
            alignment = query(expected['target_sequence'])
            results.append(alignment)

        for result, expected in zip(results, expected_alignments):
            self._check_alignment(result, expected)

    def test_kwargs_on_init_all(self):
        q_seq = "acgatcgcta"
        expected = {
            "optimal_alignment_score":1,
            "suboptimal_alignment_score":1,
            "target_begin":1,
            "target_end_optimal":1,
            "target_end_suboptimal":1,
            "query_begin":1,
            "query_end":1,
            "cigar":1,
            "query_sequence":q_seq,
            "target_sequence":1
        }
        query = StripedSmithWaterman(q_seq)
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

    def test_kwargs_on_init_score_size(self):
        pass

    def test_kwargs_on_init_protein(self):
        pass

    def test_kwargs_on_init_substitution_matrix(self):
        pass

    def test_kwargs_on_init_match(self):
        pass

    def test_kwargs_on_init_mismatch(self):
        pass

    def test_kwargs_on_init_weight_gap_open(self):
        pass

    def test_kwargs_on_init_weight_gap_extension(self):
        pass

    def test_kwargs_on_init_bit_flag(self):
        pass

    def test_kwargs_on_init_score_filter(self):
        pass

    def test_kwargs_on_init_distance_filter(self):
        pass

    def test_kwargs_on_init_mask_length(self):
        pass

    def test_kwargs_on_init_zero_based(self):
        pass

    def test_kwargs_on_call_weight_gap_open(self):
        pass

    def test_kwargs_on_call_weight_gap_extension(self):
        pass

    def test_kwargs_on_call_bit_flag(self):
        pass

    def test_kwargs_on_call_score_filter(self):
        pass

    def test_kwargs_on_call_distance_filter(self):
        pass

    def test_kwargs_on_call_mask_length(self):
        pass

    def test_kwargs_on_call_zero_based(self):
        pass

    def test_kwargs_on_call_all(self):
        q_seq = "acgatcgcta"
        expected = {
            "optimal_alignment_score":1,
            "suboptimal_alignment_score":1,
            "target_begin":1,
            "target_end_optimal":1,
            "target_end_suboptimal":1,
            "query_begin":1,
            "query_end":1,
            "cigar":1,
            "query_sequence":q_seq,
            "target_sequence":1
        }
        query = StripedSmithWaterman(q_seq)
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

    def test_kwargs_override_all(self):  
        q_seq = "acgatcgcta"
        expected = {
            "optimal_alignment_score":1,
            "suboptimal_alignment_score":1,
            "target_begin":1,
            "target_end_optimal":1,
            "target_end_suboptimal":1,
            "query_begin":1,
            "query_end":1,
            "cigar":1,
            "query_sequence":q_seq,
            "target_sequence":1
        }
        query = StripedSmithWaterman(q_seq)
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

    def test_uppercase(self):
        q_seq = "acgatcgcta"
        expected = {
            "optimal_alignment_score":1,
            "suboptimal_alignment_score":1,
            "target_begin":1,
            "target_end_optimal":1,
            "target_end_suboptimal":1,
            "query_begin":1,
            "query_end":1,
            "cigar":1,
            "query_sequence":q_seq,
            "target_sequence":1
        }
        query = StripedSmithWaterman(q_seq)
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

    def test_lowercase(self):
        q_seq = "acgatcgcta"
        expected = {
            "optimal_alignment_score":1,
            "suboptimal_alignment_score":1,
            "target_begin":1,
            "target_end_optimal":1,
            "target_end_suboptimal":1,
            "query_begin":1,
            "query_end":1,
            "cigar":1,
            "query_sequence":q_seq,
            "target_sequence":1
        }
        query = StripedSmithWaterman(q_seq)
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

class TestStripedSmithWatermanAlignment(TestCase):
    pass


if __name__ == '__main__':
    main()




