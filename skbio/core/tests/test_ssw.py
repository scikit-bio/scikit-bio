#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2014--, biocore team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

#Special thanks to http://www.faculty.ucr.edu/~mmaduro/random.htm for the random DNA generator

from unittest import TestCase, main
from skbio.core.ssw import StripedSmithWaterman, striped_smith_waterman_alignment

class TestSSW(TestCase):
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

class TestStripedSmithWaterman(TestSSW):
    def test_object_is_reusable(self):
        """skbio.core.ssw.StripedSmithWaterman: The StripedSmithWaterman object should be reusable
        """
        q_seq = "AGGGTAATTAGGCGTGTTCACCTA"
        expected_alignments = [
            {
                'optimal_alignment_score': 10,
                'suboptimal_alignment_score': 10,
                'query_begin': 4,
                'query_end': 8,
                'target_begin': 3,
                'target_end_optimal': 7,
                'target_end_suboptimal': 34,
                'cigar': '5M',
                'query_sequence': q_seq,
                'target_sequence': 'TTATAATTTTCTTATTATTATCAATATTTATAATTTGATTTTGTTGTAAT'
            },
            {
                'optimal_alignment_score': 36,
                'suboptimal_alignment_score': 16,
                'query_begin': 0,
                'query_end': 23,
                'target_begin': 6,
                'target_end_optimal': 29,
                'target_end_suboptimal': 13,
                'cigar': '8M1D8M1I7M',
                'query_sequence': q_seq,
                'target_sequence': 'AGTCGAAGGGTAATATAGGCGTGTCACCTA'
            },
            {
                'optimal_alignment_score': 16,
                'suboptimal_alignment_score': 0,
                'query_begin': 0,
                'query_end': 7,
                'target_begin': 6,
                'target_end_optimal': 13,
                'target_end_suboptimal': 0,
                'cigar': '8M',
                'query_sequence': q_seq,
                'target_sequence': 'AGTCGAAGGGTAATA'
            },
            {
                'optimal_alignment_score': 8,
                'suboptimal_alignment_score': 8,
                'query_begin': 0,
                'query_end': 3,
                'target_begin': 7,
                'target_end_optimal': 10,
                'target_end_suboptimal': 42,
                'cigar': '4M',
                'query_sequence': q_seq,
                'target_sequence': 'CTGCCTCAGGGGGAGGAAAGCGTCAGCGCGGCTGCCGTCGGCGCAGGGGC'
            },
            {
                'optimal_alignment_score': 48,
                'suboptimal_alignment_score': 16,
                'query_begin': 0,
                'query_end': 23,
                'target_begin': 0,
                'target_end_optimal': 23,
                'target_end_suboptimal': 7,
                'cigar': '24M',
                'query_sequence': q_seq,
                'target_sequence': q_seq
            }
        ]      
        query = StripedSmithWaterman(q_seq)
        results = []
        for expected in expected_alignments:
            alignment = query(expected['target_sequence'])
            results.append(alignment)

        for result, expected in zip(results, expected_alignments):
            self._check_alignment(result, expected)

    def test_many_arguments_working_together(self):
        """skbio.core.ssw.StripedSmithWaterman: Many arguments should work together
        """
        pass

    def test_score_size_is_settable(self):
        """skbio.core.ssw.StripedSmithWaterman: score_size should work
        """
        pass

    def test_protein_sequence_is_usable(self):
        """skbio.core.ssw.StripedSmithWaterman: Protein sequences should be valid input
        """
        pass

    def test_substitution_matrix_is_usable(self):
        """skbio.core.ssw.StripedSmithWaterman: Should be able to provide a substitution_matrix instead of a match/mismatch
        """
        pass

    def test_match_score_is_settable(self):
        """skbio.core.ssw.StripedSmithWaterman: Should be able to set the match score
        """
        pass

    def test_mismatch_score_is_settable(self):
        """skbio.core.ssw.StripedSmithWaterman: Should be able to set the mismatch score
        """
        pass

    def test_weight_gap_open_is_settable(self):
        """skbio.core.ssw.StripedSmithWaterman: Should be able to change the weight gap open
        """
        pass

    def test_weight_gap_extension_is_settable(self):
        """skbio.core.ssw.StripedSmithWaterman: Should be able to change the weight gap extension
        """
        pass

    def test_score_filter_is_settable_and_used(self):
        """skbio.core.ssw.StripedSmithWaterman: Should filter alignments that are less than the score_filter
        """
        pass

    def test_distance_filter_is_settable_and_used(self):
        """skbio.core.ssw.StripedSmithWaterman: Should filter all alignments that are less than the distance_filter
        """
        pass
        
    def test_mask_length_works(self):
        """skbio.core.ssw.StripedSmithWaterman: Mask length should do whatever the hell it does
        """
        pass

    def test_matrix_overrides_match_and_mismatch(self):
        """skbio.core.ssw.StripedSmithWaterman: Providing a substitution_matrix should override match and mismatch
        """
        pass

    def test_zero_index_changes_base_of_index_to_0_or_1(self):
        """skbio.core.ssw.StripedSmithWaterman: Setting zero_index (True/False) should change the base of the index (0/1)
        """
        expected_alignments = [
            ({
                'optimal_alignment_score': 100,
                'suboptimal_alignment_score': 44,
                'query_begin': 5,
                'query_end': 54,
                'target_begin': 0,
                'target_end_optimal': 49,
                'target_end_suboptimal': 21,
                'cigar': '50M',
                'query_sequence': 'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence': 'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, True),
            ({
                'optimal_alignment_score': 100,
                'suboptimal_alignment_score': 44,
                'query_begin':6,
                'query_end': 55,
                'target_begin': 1,
                'target_end_optimal': 50,
                'target_end_suboptimal': 22,
                'cigar': '50M',
                'query_sequence': 'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence': 'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, False)
        ]
        for expected, z in expected_alignments:
            query = StripedSmithWaterman(expected['query_sequence'], zero_index=z)
            alignment = query(expected['target_sequence'])
            self._check_alignment(alignment, expected)

    def test_lowercase_is_valid_sequence(self):
        """skbio.core.ssw.StripedSmithWaterman: Lowercase should be a valid input
        """
        expected = {
            'optimal_alignment_score': 23,
            'suboptimal_alignment_score': 10,
            'query_begin': 0,
            'query_end': 16,
            'target_begin': 0,
            'target_end_optimal': 20,
            'target_end_suboptimal': 4,
            'cigar': '6M4D11M',
            'query_sequence': 'aaacgataaatccgcgta',
            'target_sequence': 'aaacgactactaaatccgcgtgatagggga'
        }
        query = StripedSmithWaterman(expected['query_sequence'])
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

class TestStripedSmithWatermanAlignment(TestSSW):
    def test_same_as_using_StripedSmithWaterman_object(self):
        """skbio.core.ssw.striped_smith_waterman_alignment: Should be identical to using StripedSmithWaterman object
        """
        query_sequence = 'ATGGAAGCTATAAGCGCGGGTGAG'
        target_sequence = 'AACTTATATAATAAAAATTATATATTCGTTGGGTTCTTTTGATATAAATC'
        query = StripedSmithWaterman(query_sequence)
        align1 = query(target_sequence)
        align2 = striped_smith_waterman_alignment(query_sequence, target_sequence)
        self._check_alignment(align2, align1)

    def test_kwargs_are_usable(self):
        """skbio.core.ssw.striped_smith_waterman_alignment: Should accept kwargs
        """
        kwargs = {}
        kwargs['zero_index'] = False
        kwargs['match'] = 5
        query_sequence = 'AGGGTAATTAGGCGTGTTCACCTA'
        target_sequence = 'TACTTATAAGATGTCTCAACGGCATGCGCAACTTGTGAAGTG'
        query = StripedSmithWaterman(query_sequence, **kwargs)
        align1 = query(target_sequence)
        align2 = striped_smith_waterman_alignment(query_sequence, target_sequence, **kwargs)
        self._check_alignment(align2, align1)

class TestAlignmentStructure(TestSSW):
    def test_works_for_dot_and_square_bracket_access(self):
        """skbio.core.ssw.AlignmentStructure: The alignment should be accessible by . and []
        """
        q_seq = "ACTAGACTCTCTCGAGCGCGCTATATATCG"
        query = StripedSmithWaterman(q_seq)
        alignment = query("ACCTACGCGAGCATCACTAGACTCTCTCGAGCGCGCTATAAAACTCAGCTCAC")
        for accessible in self.align_attributes:
            self.assertEqual(getattr(alignment, accessible), alignment[accessible])

    def test_is_zero_based_returns_true_if_index_base_is_zero(self):
        """skbio.core.ssw.AlignmentStructure: is_zero_based should return True if index base is 0 else False
        """
        expected_alignments = [
            ({
                'query_sequence': 'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence': 'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, True),
            ({
                'query_sequence': 'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence': 'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, False)
        ]
        for expected, z in expected_alignments:
            query = StripedSmithWaterman(expected['query_sequence'], zero_index=z)
            alignment = query(expected['target_sequence'])
            self.assertEqual(z, alignment.is_zero_based())

    def test_set_zero_based_changes_the_index_base(self):
        """skbio.core.ssw.AlignmentStructure: set_zero_based should set the index base to 0 if True, 1 if False
        """
        expected_alignments = [
            ({
                'query_sequence': 'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence': 'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, True),
            ({
                'query_sequence': 'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence': 'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, False)
        ]
        for expected, z in expected_alignments:
            query = StripedSmithWaterman(expected['query_sequence'], zero_index=z)
            alignment = query(expected['target_sequence'])
            alignment.set_zero_based(not z)
            self.assertEqual(not z, alignment.is_zero_based())

if __name__ == '__main__':
    main()




