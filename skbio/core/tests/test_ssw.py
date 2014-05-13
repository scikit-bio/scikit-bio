#!/usr/bin/env python
# -----------------------------------------------------------------------------
#  Copyright (c) 2013--, scikit-bio development team.
#
#  Distributed under the terms of the Modified BSD License.
#
#  The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

# Special thanks to http://www.faculty.ucr.edu/~mmaduro/random.htm for the
# random DNA generator.

# These really only test against regression, not correctness.
# It is assumed that ssw.c and ssw.h are correct.
# Furthermore all expected results are created by running StripedSmithWaterma
# the resulting alignments are verified by hand. Creating tests from the base
# C API is impractical at this time.

from unittest import TestCase, main
from skbio.core.ssw import (StripedSmithWaterman,
                            striped_smith_waterman_alignment)


class TestSSW(TestCase):

    align_attributes = [
        "optimal_alignment_score", "suboptimal_alignment_score",
        "target_begin", "target_end_optimal", "target_end_suboptimal",
        "query_begin", "query_end", "cigar", "query_sequence",
        "target_sequence"
    ]

    def _check_alignment(self, alignment, expected):
        for attribute in self.align_attributes:
            # The first element of this tuple is to identify
            # the broken sequence if one should fail
            self.assertEqual((alignment['target_sequence'],
                              alignment[attribute]),
                             (expected['target_sequence'],
                              expected[attribute]))


class TestStripedSmithWaterman(TestSSW):

    # This is temporary: blosum50 does not exist in skbio yet as per
    # issue 161. When the issue is resolved, this should be removed in favor
    # of an import.
    blosum50 = \
        {
            '*': {'*': 1, 'A': -5, 'C': -5, 'B': -5, 'E': -5, 'D': -5, 'G': -5,
                  'F': -5, 'I': -5, 'H': -5, 'K': -5, 'M': -5, 'L': -5,
                  'N': -5, 'Q': -5, 'P': -5, 'S': -5, 'R': -5, 'T': -5,
                  'W': -5, 'V': -5, 'Y': -5, 'X': -5, 'Z': -5},
            'A': {'*': -5, 'A': 5, 'C': -1, 'B': -2, 'E': -1, 'D': -2, 'G': 0,
                  'F': -3, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -2,
                  'N': -1, 'Q': -1, 'P': -1, 'S': 1, 'R': -2, 'T': 0, 'W': -3,
                  'V': 0, 'Y': -2, 'X': -1, 'Z': -1},
            'C': {'*': -5, 'A': -1, 'C': 13, 'B': -3, 'E': -3, 'D': -4,
                  'G': -3, 'F': -2, 'I': -2, 'H': -3, 'K': -3, 'M': -2,
                  'L': -2, 'N': -2, 'Q': -3, 'P': -4, 'S': -1, 'R': -4,
                  'T': -1, 'W': -5, 'V': -1, 'Y': -3, 'X': -1, 'Z': -3},
            'B': {'*': -5, 'A': -2, 'C': -3, 'B': 6, 'E': 1, 'D': 6, 'G': -1,
                  'F': -4, 'I': -4, 'H': 0, 'K': 0, 'M': -3, 'L': -4, 'N': 5,
                  'Q': 0, 'P': -2, 'S': 0, 'R': -1, 'T': 0, 'W': -5, 'V': -3,
                  'Y': -3, 'X': -1, 'Z': 1},
            'E': {'*': -5, 'A': -1, 'C': -3, 'B': 1, 'E': 6, 'D': 2, 'G': -3,
                  'F': -3, 'I': -4, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0,
                  'Q': 2, 'P': -1, 'S': -1, 'R': 0, 'T': -1, 'W': -3, 'V': -3,
                  'Y': -2, 'X': -1, 'Z': 5},
            'D': {'*': -5, 'A': -2, 'C': -4, 'B': 6, 'E': 2, 'D': 8, 'G': -1,
                  'F': -5, 'I': -4, 'H': -1, 'K': -1, 'M': -4, 'L': -4, 'N': 2,
                  'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -5, 'V': -4,
                  'Y': -3, 'X': -1, 'Z': 1},
            'G': {'*': -5, 'A': 0, 'C': -3, 'B': -1, 'E': -3, 'D': -1, 'G': 8,
                  'F': -4, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0,
                  'Q': -2, 'P': -2, 'S': 0, 'R': -3, 'T': -2, 'W': -3, 'V': -4,
                  'Y': -3, 'X': -1, 'Z': -2},
            'F': {'*': -5, 'A': -3, 'C': -2, 'B': -4, 'E': -3, 'D': -5,
                  'G': -4, 'F': 8, 'I': 0, 'H': -1, 'K': -4, 'M': 0, 'L': 1,
                  'N': -4, 'Q': -4, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 1,
                  'V': -1, 'Y': 4, 'X': -1, 'Z': -4},
            'I': {'*': -5, 'A': -1, 'C': -2, 'B': -4, 'E': -4, 'D': -4,
                  'G': -4, 'F': 0, 'I': 5, 'H': -4, 'K': -3, 'M': 2, 'L': 2,
                  'N': -3, 'Q': -3, 'P': -3, 'S': -3, 'R': -4, 'T': -1,
                  'W': -3, 'V': 4, 'Y': -1, 'X': -1, 'Z': -3},
            'H': {'*': -5, 'A': -2, 'C': -3, 'B': 0, 'E': 0, 'D': -1, 'G': -2,
                  'F': -1, 'I': -4, 'H': 10, 'K': 0, 'M': -1, 'L': -3, 'N': 1,
                  'Q': 1, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -3, 'V': -4,
                  'Y': 2, 'X': -1, 'Z': 0},
            'K': {'*': -5, 'A': -1, 'C': -3, 'B': 0, 'E': 1, 'D': -1, 'G': -2,
                  'F': -4, 'I': -3, 'H': 0, 'K': 6, 'M': -2, 'L': -3, 'N': 0,
                  'Q': 2, 'P': -1, 'S': 0, 'R': 3, 'T': -1, 'W': -3, 'V': -3,
                  'Y': -2, 'X': -1, 'Z': 1},
            'M': {'*': -5, 'A': -1, 'C': -2, 'B': -3, 'E': -2, 'D': -4,
                  'G': -3, 'F': 0, 'I': 2, 'H': -1, 'K': -2, 'M': 7, 'L': 3,
                  'N': -2, 'Q': 0, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -1,
                  'V': 1, 'Y': 0, 'X': -1, 'Z': -1},
            'L': {'*': -5, 'A': -2, 'C': -2, 'B': -4, 'E': -3, 'D': -4,
                  'G': -4, 'F': 1, 'I': 2, 'H': -3, 'K': -3, 'M': 3, 'L': 5,
                  'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -1,
                  'W': -2, 'V': 1, 'Y': -1, 'X': -1, 'Z': -3},
            'N': {'*': -5, 'A': -1, 'C': -2, 'B': 5, 'E': 0, 'D': 2, 'G': 0,
                  'F': -4, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -4, 'N': 7,
                  'Q': 0, 'P': -2, 'S': 1, 'R': -1, 'T': 0, 'W': -4, 'V': -3,
                  'Y': -2, 'X': -1, 'Z': 0},
            'Q': {'*': -5, 'A': -1, 'C': -3, 'B': 0, 'E': 2, 'D': 0, 'G': -2,
                  'F': -4, 'I': -3, 'H': 1, 'K': 2, 'M': 0, 'L': -2, 'N': 0,
                  'Q': 7, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -1, 'V': -3,
                  'Y': -1, 'X': -1, 'Z': 4},
            'P': {'*': -5, 'A': -1, 'C': -4, 'B': -2, 'E': -1, 'D': -1,
                  'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -3,
                  'L': -4, 'N': -2, 'Q': -1, 'P': 10, 'S': -1, 'R': -3,
                  'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': -1},
            'S': {'*': -5, 'A': 1, 'C': -1, 'B': 0, 'E': -1, 'D': 0, 'G': 0,
                  'F': -3, 'I': -3, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1,
                  'Q': 0, 'P': -1, 'S': 5, 'R': -1, 'T': 2, 'W': -4, 'V': -2,
                  'Y': -2, 'X': -1, 'Z': 0},
            'R': {'*': -5, 'A': -2, 'C': -4, 'B': -1, 'E': 0, 'D': -2, 'G': -3,
                  'F': -3, 'I': -4, 'H': 0, 'K': 3, 'M': -2, 'L': -3, 'N': -1,
                  'Q': 1, 'P': -3, 'S': -1, 'R': 7, 'T': -1, 'W': -3, 'V': -3,
                  'Y': -1, 'X': -1, 'Z': 0},
            'T': {'*': -5, 'A': 0, 'C': -1, 'B': 0, 'E': -1, 'D': -1, 'G': -2,
                  'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0,
                  'Q': -1, 'P': -1, 'S': 2, 'R': -1, 'T': 5, 'W': -3, 'V': 0,
                  'Y': -2, 'X': -1, 'Z': -1},
            'W': {'*': -5, 'A': -3, 'C': -5, 'B': -5, 'E': -3, 'D': -5,
                  'G': -3, 'F': 1, 'I': -3, 'H': -3, 'K': -3, 'M': -1, 'L': -2,
                  'N': -4, 'Q': -1, 'P': -4, 'S': -4, 'R': -3, 'T': -3,
                  'W': 15, 'V': -3, 'Y': 2, 'X': -1, 'Z': -2},
            'V': {'*': -5, 'A': 0, 'C': -1, 'B': -3, 'E': -3, 'D': -4, 'G': -4,
                  'F': -1, 'I': 4, 'H': -4, 'K': -3, 'M': 1, 'L': 1, 'N': -3,
                  'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 5,
                  'Y': -1, 'X': -1, 'Z': -3},
            'Y': {'*': -5, 'A': -2, 'C': -3, 'B': -3, 'E': -2, 'D': -3,
                  'G': -3, 'F': 4, 'I': -1, 'H': 2, 'K': -2, 'M': 0, 'L': -1,
                  'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -1, 'T': -2, 'W': 2,
                  'V': -1, 'Y': 8, 'X': -1, 'Z': -2},
            'X': {'*': -5, 'A': -1, 'C': -1, 'B': -1, 'E': -1, 'D': -1,
                  'G': -1, 'F': -1, 'I': -1, 'H': -1, 'K': -1, 'M': -1,
                  'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1,
                  'T': -1, 'W': -1, 'V': -1, 'Y': -1, 'X': -1, 'Z': -1},
            'Z': {'*': -5, 'A': -1, 'C': -3, 'B': 1, 'E': 5, 'D': 1, 'G': -2,
                  'F': -4, 'I': -3, 'H': 0, 'K': 1, 'M': -1, 'L': -3, 'N': 0,
                  'Q': 4, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -2, 'V': -3,
                  'Y': -2, 'X': -1, 'Z': 5}
        }

    def test_object_is_reusable(self):
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
                'target_sequence':
                    'TTATAATTTTCTTATTATTATCAATATTTATAATTTGATTTTGTTGTAAT'
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
                'target_sequence':
                    'CTGCCTCAGGGGGAGGAAAGCGTCAGCGCGGCTGCCGTCGGCGCAGGGGC'
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

    def test_regression_on_instantiation_arguments(self):
        expected = {
            'optimal_alignment_score': 23,
            'suboptimal_alignment_score': 10,
            'query_begin': 0,
            'query_end': 16,
            'target_begin': 0,
            'target_end_optimal': 20,
            'target_end_suboptimal': 4,
            'cigar': '6M4D11M',
            'query_sequence': 'AAACGATAAATCCGCGTA',
            'target_sequence': 'AAACGACTACTAAATCCGCGTGATAGGGGA'
        }
        query = StripedSmithWaterman(expected['query_sequence'],
                                     weight_gap_open=5,
                                     weight_gap_extension=2,
                                     score_size=2,
                                     mask_length=15,
                                     mask_auto=True,
                                     score_only=False,
                                     score_filter=None,
                                     distance_filter=None,
                                     override_skip_babp=False,
                                     protein=False,
                                     match=2,
                                     mismatch=3,
                                     substitution_matrix=None,
                                     suppress_sequences=False,
                                     zero_index=True)
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

    def test_protein_sequence_is_usable(self):
        pass

    def test_substitution_matrix(self):
        pass

    def test_match_score(self):
        pass

    def test_mismatch_score(self):
        pass

    def test_weight_gap_open(self):
        pass

    def test_weight_gap_extension(self):
        pass

    def test_score_filter_is_used(self):
        pass

    def test_distance_filter_is_used(self):
        pass

    def test_mask_length_works(self):
        pass

    def test_matrix_overrides_match_and_mismatch(self):
        pass

    def test_zero_index_changes_base_of_index_to_0_or_1(self):
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
                'query_sequence':
                    'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence':
                    'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, True),
            ({
                'optimal_alignment_score': 100,
                'suboptimal_alignment_score': 44,
                'query_begin': 6,
                'query_end': 55,
                'target_begin': 1,
                'target_end_optimal': 50,
                'target_end_suboptimal': 22,
                'cigar': '50M',
                'query_sequence':
                    'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence':
                    'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, False)
        ]
        for expected, z in expected_alignments:
            query = StripedSmithWaterman(expected['query_sequence'],
                                         zero_index=z)
            alignment = query(expected['target_sequence'])
            self._check_alignment(alignment, expected)

    def test_lowercase_is_valid_sequence(self):
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
        query_sequence = 'ATGGAAGCTATAAGCGCGGGTGAG'
        target_sequence = 'AACTTATATAATAAAAATTATATATTCGTTGGGTTCTTTTGATATAAATC'
        query = StripedSmithWaterman(query_sequence)
        align1 = query(target_sequence)
        align2 = striped_smith_waterman_alignment(query_sequence,
                                                  target_sequence)
        self._check_alignment(align2, align1)

    def test_kwargs_are_usable(self):
        kwargs = {}
        kwargs['zero_index'] = False
        kwargs['match'] = 5
        query_sequence = 'AGGGTAATTAGGCGTGTTCACCTA'
        target_sequence = 'TACTTATAAGATGTCTCAACGGCATGCGCAACTTGTGAAGTG'
        query = StripedSmithWaterman(query_sequence, **kwargs)
        align1 = query(target_sequence)
        align2 = striped_smith_waterman_alignment(query_sequence,
                                                  target_sequence, **kwargs)
        self._check_alignment(align2, align1)


class TestAlignmentStructure(TestSSW):

    def test_works_for_dot_and_square_bracket_access(self):
        q_seq = "AGGGTAATTAGGCGTGTTCACCTA"
        query = StripedSmithWaterman(q_seq)
        alignment = query("TACTTATAAGATGTCTCAACGGCATGCGCAACTTGTGAAGTG")
        for accessible in self.align_attributes:
            self.assertEqual(getattr(alignment, accessible),
                             alignment[accessible])

    def test_is_zero_based_returns_true_if_index_base_is_zero(self):
        expected_alignments = [
            ({
                'query_sequence':
                    'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence':
                    'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, True),
            ({
                'query_sequence':
                    'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence':
                    'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, False)
        ]
        for expected, z in expected_alignments:
            query = StripedSmithWaterman(expected['query_sequence'],
                                         zero_index=z)
            alignment = query(expected['target_sequence'])
            self.assertEqual(z, alignment.is_zero_based())

    def test_set_zero_based_changes_the_index_base(self):
        expected_alignments = [
            ({
                'query_sequence':
                    'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence':
                    'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, True),
            ({
                'query_sequence':
                    'AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC',
                'target_sequence':
                    'CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC'
            }, False)
        ]
        for expected, z in expected_alignments:
            query = StripedSmithWaterman(expected['query_sequence'],
                                         zero_index=z)
            alignment = query(expected['target_sequence'])
            alignment.set_zero_based(not z)
            self.assertEqual(not z, alignment.is_zero_based())

if __name__ == '__main__':
    main()
