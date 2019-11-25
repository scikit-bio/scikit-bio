# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# Special thanks to http://www.faculty.ucr.edu/~mmaduro/random.htm for the
# random DNA generator.

# These tests confirm that StripedSmithWaterman returns the same results as
# SSW. We don't test for correctness of those results (i.e., we assume that
# ssw.c and ssw.h are correct) as that testing is beyond the scope of skbio.
# Furthermore all expected results are created by running StripedSmithWaterman
# the resulting alignments are verified by hand. Creating tests from the base
# C API is impractical at this time.

from unittest import TestCase, main

from skbio import (local_pairwise_align_ssw, Sequence, DNA, RNA, Protein,
                   TabularMSA)
from skbio.alignment import StripedSmithWaterman, AlignmentStructure
from skbio.alignment._pairwise import blosum50


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
            self.assertEqual((expected['target_sequence'],
                              expected[attribute]),
                             (alignment['target_sequence'],
                              alignment[attribute]))

    def _check_argument_with_inequality_on_optimal_align_score(
            self,
            query_sequences=None,
            target_sequences=None,
            arg=None,
            default=None,
            i_range=None,
            compare_lt=None,
            compare_gt=None):
        iterable_kwarg = {}
        default_kwarg = {}
        default_kwarg[arg] = default
        for query_sequence in query_sequences:
            for target_sequence in target_sequences:
                for i in i_range:
                    iterable_kwarg[arg] = i
                    query1 = StripedSmithWaterman(query_sequence,
                                                  **iterable_kwarg)
                    align1 = query1(target_sequence)

                    query2 = StripedSmithWaterman(query_sequence,
                                                  **default_kwarg)
                    align2 = query2(target_sequence)

                    if i == default:
                        self.assertEqual(align1.optimal_alignment_score,
                                         align2.optimal_alignment_score)
                    if i < default:
                        compare_lt(align1.optimal_alignment_score,
                                   align2.optimal_alignment_score)
                    if i > default:
                        compare_gt(align1.optimal_alignment_score,
                                   align2.optimal_alignment_score)

    def _check_bit_flag_sets_properties_falsy_or_negative(
            self,
            query_sequences=None,
            target_sequences=None,
            arg_settings=[],
            properties_to_null=[]):
        kwarg = {}

        def falsy_or_negative(alignment, prop):
            if type(alignment[prop]) is int:
                return alignment[prop] < 0
            else:
                return not alignment[prop]

        for query_sequence in query_sequences:
            for target_sequence in target_sequences:
                for arg, setting in arg_settings:
                    kwarg[arg] = setting
                query = StripedSmithWaterman(query_sequence, **kwarg)
                alignment = query(target_sequence)
                for prop in properties_to_null:
                    self.assertTrue(falsy_or_negative(alignment, prop))
                # Every property not in our null list
                for prop in [p for p in self.align_attributes
                             if p not in properties_to_null]:
                    self.assertFalse(falsy_or_negative(alignment, prop))


class TestStripedSmithWaterman(TestSSW):

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
                'target_sequence': ('TTATAATTTTCTTATTATTATCAATATTTATAATTTGATTT'
                                    'TGTTGTAAT')
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
                'target_sequence': ('CTGCCTCAGGGGGAGGAAAGCGTCAGCGCGGCTGCCGTCGG'
                                    'CGCAGGGGC')
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
                                     gap_open_penalty=5,
                                     gap_extend_penalty=2,
                                     score_size=2,
                                     mask_length=15,
                                     mask_auto=True,
                                     score_only=False,
                                     score_filter=None,
                                     distance_filter=None,
                                     override_skip_babp=False,
                                     protein=False,
                                     match_score=2,
                                     mismatch_score=-3,
                                     substitution_matrix=None,
                                     suppress_sequences=False,
                                     zero_index=True)
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

    def test_protein_sequence_is_usable(self):
        expected = {
            'optimal_alignment_score': 316,
            'suboptimal_alignment_score': 95,
            'query_begin': 0,
            'query_end': 52,
            'target_begin': 0,
            'target_end_optimal': 52,
            'target_end_suboptimal': 18,
            'cigar': '15M1D15M1I22M',
            'query_sequence': ('VHLTGEEKSAVAALWGKVNVDEVGGEALGRXLLVVYPWTQRFFESF'
                               'SDLSTPDABVMSNPKVKAHGK'),
            'target_sequence': ('VHLTPEEKSAVTALWBGKVNVDEVGGEALGRLLVVYPWTQRFFES'
                                'FGDLSTPD*')
        }
        query = StripedSmithWaterman(expected['query_sequence'],
                                     protein=True,
                                     substitution_matrix=blosum50)
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

    def test_align_with_N_in_nucleotide_sequence(self):
        expected = {
            'optimal_alignment_score': 9,
            'suboptimal_alignment_score': 0,
            'query_begin': 0,
            'query_end': 8,
            'target_begin': 0,
            'target_end_optimal': 9,
            'target_end_suboptimal': 0,
            'cigar': '4M1D5M',
            'query_sequence': 'ACTCANNATCGANCTAGC',
            'target_sequence': 'ACTCGAAAATGTNNGCA'
        }
        query = StripedSmithWaterman(expected['query_sequence'])
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

    def test_arg_match_score(self):
        query_sequences = [
            "TTTTTTCTTATTATTATCAATATTTATAATTTGATTTTGTTGTAAT",
            "AGTCGAAGGGTCAATATAGGCGTGTCACCTA",
            "AGTCGAAGGGTAATA",
            "CTGCCTCAAGGGGGAGGAAAGCGTCAGCGCGGCTGCCGTCGGCGCAGGGGC",
            "AGGGTAATTTTAGGCGTGTTCACCTA"
        ]
        target_sequences = query_sequences
        self._check_argument_with_inequality_on_optimal_align_score(
            query_sequences=query_sequences,
            target_sequences=target_sequences,
            arg='match_score',
            default=2,
            i_range=range(0, 5),
            compare_lt=self.assertLess,
            compare_gt=self.assertGreater
        )
        # The above is a strict bound, so we don't need a expected align

    def test_arg_mismatch_score(self):
        query_sequences = [
            "TTATAATTAATTCTTATTATTATCAATATTTATAATTTGATTTTGTTGTAAT",
            "AGTCGAAGGGTAAGGGGTATAGGCGTGTCACCTA",
            "AGTCGAAGGGTAATA",
            "CTGCCTCAGGGGCGAGGAAAGCGTCAGCGCGGCTGCCGTCGGCGCAGGGGC",
            "AGGGTAATTAGCGCGTGTTCACCTA"
        ]
        target_sequences = query_sequences
        self._check_argument_with_inequality_on_optimal_align_score(
            query_sequences=query_sequences,
            target_sequences=target_sequences,
            arg='mismatch_score',
            default=-3,
            i_range=range(-6, 1),
            # These are intentionally inverted
            compare_lt=self.assertLessEqual,
            compare_gt=self.assertGreaterEqual
        )
        # The above is not a strict bound, so lets use an expected align
        # to plug the hole where every align is exactly equal to default
        expected = {
            'optimal_alignment_score': 8,
            'suboptimal_alignment_score': 0,
            'query_begin': 5,
            'query_end': 8,
            'target_begin': 10,
            'target_end_optimal': 13,
            'target_end_suboptimal': 0,
            'cigar': '4M',
            'query_sequence': 'AGAGGGTAATCAGCCGTGTCCACCGGAACACAACGCTATCGGGCGA',
            'target_sequence': 'GTTCGCCCCAGTAAAGTTGCTACCAAATCCGCATG'
        }
        query = StripedSmithWaterman(expected['query_sequence'],
                                     mismatch_score=-8)
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

    def test_arg_matrix_overrides_match_and_mismatch(self):
        query_sequences = [
            "TTATAATTAATTCTTATTATTATCAATATTTATAATTTGATTTTGTTGTAAT",
            "AGTCGAAGGGTAAGGGGTATAGGCGTGTCACCTA",
            "AGTCGAAGGGTAATA",
            "CTGCCTCAGGGGCGAGGAAAGCGTCAGCGCGGCTGCCGTCGGCGCAGGGGC",
            "AGGGTAATTAGCGCGTGTTCACCTA"
        ]
        target_sequences = query_sequences
        matrix = {  # This is a biologically meaningless matrix
            "A": {"A": 4,  "T": -1, "C": -2, "G": -3, "N": 4},
            "T": {"A": -1, "T": 1,  "C": -1, "G": -4, "N": 1},
            "C": {"A": -2, "T": -1, "C": 10, "G": 1,  "N": 1},
            "G": {"A": -3, "T": -4, "C": 1,  "G": 3,  "N": 1},
            "N": {"A": 4,  "T": 1,  "C": 1,  "G": 1,  "N": 0}
        }
        for query_sequence in query_sequences:
            for target_sequence in target_sequences:
                query1 = StripedSmithWaterman(query_sequence)
                align1 = query1(target_sequence)

                query2 = StripedSmithWaterman(query_sequence,
                                              substitution_matrix=matrix)
                align2 = query2(target_sequence)

                self.assertNotEqual(align1.optimal_alignment_score,
                                    align2.optimal_alignment_score)

    def test_arg_gap_open_penalty(self):
        query_sequences = [
            "TTATAATTTTCTTAGTTATTATCAATATTTATAATTTGATTTTGTTGTAAT",
            "AGTCCGAAGGGTAATATAGGCGTGTCACCTA",
            "AGTCGAAGGCGGTAATA",
            "CTGCCTCGGCAGGGGGAGGAAAGCGTCAGCGCGGCTGCCGTCGGCGCAGGGGC",
            "AGGGTAATTAAAGGCGTGTTCACCTA"
        ]
        target_sequences = query_sequences
        self._check_argument_with_inequality_on_optimal_align_score(
            query_sequences=query_sequences,
            target_sequences=target_sequences,
            arg='gap_open_penalty',
            default=5,
            i_range=range(1, 12),
            # These are intentionally inverted
            compare_lt=self.assertGreaterEqual,
            compare_gt=self.assertLessEqual
        )
        # The above is not a strict bound, so lets use an expected align
        # to plug the hole where every align is exactly equal to default
        expected = {
            'optimal_alignment_score': 51,
            'suboptimal_alignment_score': 20,
            'query_begin': 0,
            'query_end': 37,
            'target_begin': 0,
            'target_end_optimal': 29,
            'target_end_suboptimal': 9,
            'cigar': '5M4I3M3I1M1I21M',
            'query_sequence': 'TAGAGATTAATTGCCACATTGCCACTGCCAAAATTCTG',
            'target_sequence': 'TAGAGATTAATTGCCACTGCCAAAATTCTG'
        }
        query = StripedSmithWaterman(expected['query_sequence'],
                                     gap_open_penalty=1)
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

    def test_arg_gap_extend_penalty(self):
        query_sequences = [
            "TTATAATTTTCTTATTATTATCAATATTTATAATTTGATTTTGTTGTAAT",
            "AGTCGAAGGGTAATACTAGGCGTGTCACCTA",
            "AGTCGAAGGGTAATA",
            "CTGCCTCAGGGGGAGGCAAAGCGTCAGCGCGGCTGCCGTCGGCGCAGGGGC",
            "AGGGTAATTAGGCGTGTTCACCTA"
        ]
        target_sequences = query_sequences
        self._check_argument_with_inequality_on_optimal_align_score(
            query_sequences=query_sequences,
            target_sequences=target_sequences,
            arg='gap_extend_penalty',
            default=2,
            i_range=range(1, 10),
            # These are intentionally inverted
            compare_lt=self.assertGreaterEqual,
            compare_gt=self.assertLessEqual
        )
        # The above is not a strict bound, so lets use an expected align
        # to plug the hole where every align is exactly equal to default
        expected = {
            'optimal_alignment_score': 9,
            'suboptimal_alignment_score': 8,
            'query_begin': 6,
            'query_end': 12,
            'target_begin': 7,
            'target_end_optimal': 13,
            'target_end_suboptimal': 38,
            'cigar': '7M',
            'query_sequence': 'TCTATAAGATTCCGCATGCGTTACTTATAAGATGTCTCAACGG',
            'target_sequence': 'GCCCAGTAGCTTCCCAATATGAGAGCATCAATTGTAGATCGGGCC'
        }
        query = StripedSmithWaterman(expected['query_sequence'],
                                     gap_extend_penalty=10)
        alignment = query(expected['target_sequence'])
        self._check_alignment(alignment, expected)

    def test_arg_score_only(self):
        query_sequences = [
            "TTATCGTGATTATTATCAATATTTATAATTTGATTTTGTTGTAAT",
            "AGTCGAAGGGTAATACTATAAGGCGTGTCACCTA",
            "AGTCGAAGGGTAATA",
            "AGGGTAATTAGGCGTGCGTGCGTGTTCACCTA",
            "AGGGTATTAGGCGTGTTCACCTA"
        ]
        target_sequences = query_sequences
        self._check_bit_flag_sets_properties_falsy_or_negative(
            query_sequences=query_sequences,
            target_sequences=target_sequences,
            arg_settings=[('score_only', True)],
            properties_to_null=['query_begin', 'target_begin', 'cigar']
        )

    def test_arg_score_filter_is_used(self):
        query_sequences = [
            "TTATCGTGATTATTATCAATATTTATAATTTGATTTTGTTGTAAT",
            "AGTCGAAGGGTAATACTATAAGGCGTGTCACCTA",
            "AGTCGAAGGGTAATA",
            "AGGGTAATTAGGCGTGCGTGCGTGTTCACCTA",
            "AGGGTATTAGGCGTGTTCACCTA"
        ]
        target_sequences = query_sequences
        self._check_bit_flag_sets_properties_falsy_or_negative(
            query_sequences=query_sequences,
            target_sequences=target_sequences,
            # score_filter will force a BABP and cigar to be falsy
            arg_settings=[('score_filter', 9001)],
            properties_to_null=['query_begin', 'target_begin', 'cigar']
        )

    def test_arg_distance_filter_is_used(self):
        query_sequences = [
            "TTATCGTGATTATTATCAATATTTATAATTTGATTTTGTTGTAAT",
            "AGTCGAAGGGTAATACTATAAGGCGTGTCACCTA",
            "AGTCGAAGGGTAATA",
            "AGGGTAATTAGGCGTGCGTGCGTGTTCACCTA",
            "AGGGTATTAGGCGTGTTCACCTA"
        ]
        target_sequences = query_sequences
        self._check_bit_flag_sets_properties_falsy_or_negative(
            query_sequences=query_sequences,
            target_sequences=target_sequences,
            # distance_filter will force cigar to be falsy only
            arg_settings=[('distance_filter', 1)],
            properties_to_null=['cigar']
        )

    def test_arg_override_skip_babp(self):
        query_sequences = [
            "TTATCGTGATTATTATCAATATTTATAATTTGATTTTGTTGTAAT",
            "AGTCGAAGGGTAATACTATAAGGCGTGTCACCTA",
            "AGTCGAAGGGTAATA",
            "AGGGTAATTAGGCGTGCGTGCGTGTTCACCTA",
            "AGGGTATTAGGCGTGTTCACCTA"
        ]
        target_sequences = query_sequences
        self._check_bit_flag_sets_properties_falsy_or_negative(
            query_sequences=query_sequences,
            target_sequences=target_sequences,
            # score_filter will force a BABP and cigar to be falsy if not for
            # override_skip_babp preventing this for all but the cigar
            arg_settings=[('override_skip_babp', True),
                          ('score_filter', 9001)],
            properties_to_null=['cigar']
        )

    def test_arg_zero_index_changes_base_of_index_to_0_or_1(self):
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
                'query_sequence': ('AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCG'
                                   'CCCCGGGCGGGGC'),
                'target_sequence': ('CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCC'
                                    'GGGCGGGGC')
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
                'query_sequence': ('AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCG'
                                   'CCCCGGGCGGGGC'),
                'target_sequence': ('CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCC'
                                    'GGGCGGGGC')
            }, False)
        ]
        for expected, z in expected_alignments:
            query = StripedSmithWaterman(expected['query_sequence'],
                                         zero_index=z)
            alignment = query(expected['target_sequence'])
            self._check_alignment(alignment, expected)

    def test_arg_suppress_sequences(self):
        expected = {
            'optimal_alignment_score': 100,
            'suboptimal_alignment_score': 44,
            'query_begin': 5,
            'query_end': 54,
            'target_begin': 0,
            'target_end_optimal': 49,
            'target_end_suboptimal': 21,
            'cigar': '50M',
            'query_sequence': '',
            'target_sequence': ''
        }
        query = StripedSmithWaterman(
            "AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC",
            suppress_sequences=True)
        alignment = query("CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCCGGGCGGGGC")
        self._check_alignment(alignment, expected)


class TestAlignStripedSmithWaterman(TestSSW):
    def _check_TabularMSA_to_AlignmentStructure(self, alignment, structure,
                                                expected_dtype):
        msa, score, start_end = alignment

        self.assertEqual(score, structure.optimal_alignment_score)
        self.assertEqual(
            msa,
            TabularMSA([expected_dtype(structure.aligned_query_sequence),
                        expected_dtype(structure.aligned_target_sequence)]))
        if structure.query_begin == -1:
            self.assertEqual(start_end, None)
        else:
            for (start, end), (expected_start, expected_end) in \
                zip(start_end,
                    [(structure.query_begin,
                      structure.query_end),
                     (structure.target_begin,
                      structure.target_end_optimal)]):
                self.assertEqual(start, expected_start)
                self.assertEqual(end, expected_end)

    def test_same_as_using_StripedSmithWaterman_object_DNA(self):
        query_sequence = 'ATGGAAGCTATAAGCGCGGGTGAG'
        target_sequence = 'AACTTATATAATAAAAATTATATATTCGTTGGGTTCTTTTGATATAAATC'
        query = StripedSmithWaterman(query_sequence)
        align1 = query(target_sequence)
        align2 = local_pairwise_align_ssw(DNA(query_sequence),
                                          DNA(target_sequence))
        self._check_TabularMSA_to_AlignmentStructure(align2, align1, DNA)

    def test_same_as_using_StripedSmithWaterman_object_Protein(self):
        query_sequence = 'HEAGAWGHEE'
        target_sequence = 'PAWHEAE'
        query = StripedSmithWaterman(query_sequence,
                                     protein=True,
                                     substitution_matrix=blosum50)
        align1 = query(target_sequence)
        align2 = local_pairwise_align_ssw(Protein(query_sequence),
                                          Protein(target_sequence),
                                          substitution_matrix=blosum50)
        self._check_TabularMSA_to_AlignmentStructure(align2, align1, Protein)

    def test_kwargs_are_usable(self):
        kwargs = {}
        kwargs['mismatch_score'] = -2
        kwargs['match_score'] = 5
        query_sequence = 'AGGGTAATTAGGCGTGTTCACCTA'
        target_sequence = 'TACTTATAAGATGTCTCAACGGCATGCGCAACTTGTGAAGTG'
        query = StripedSmithWaterman(query_sequence, **kwargs)
        align1 = query(target_sequence)
        align2 = local_pairwise_align_ssw(DNA(query_sequence),
                                          DNA(target_sequence), **kwargs)
        self._check_TabularMSA_to_AlignmentStructure(align2, align1, DNA)

    def test_invalid_type(self):
        with self.assertRaisesRegex(TypeError, r"not type 'Sequence'"):
            local_pairwise_align_ssw(DNA('ACGT'), Sequence('ACGT'))

        with self.assertRaisesRegex(TypeError, r"not type 'str'"):
            local_pairwise_align_ssw('ACGU', RNA('ACGU'))

    def test_type_mismatch(self):
        with self.assertRaisesRegex(TypeError, r"same type: 'DNA' != 'RNA'"):
            local_pairwise_align_ssw(DNA('ACGT'), RNA('ACGU'))


class TestAlignmentStructure(TestSSW):

    def mock_object_factory(self, dictionary):
        class MockAlignmentStructure(AlignmentStructure):
            def __init__(self, _a, _b, _c):
                for key in dictionary:
                    setattr(self, key, dictionary[key])
        return MockAlignmentStructure(None, None, 0)

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
                'query_sequence': ('AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCG'
                                   'CCCCGGGCGGGGC'),
                'target_sequence': ('CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCC'
                                    'GGGCGGGGC')
            }, True),
            ({
                'query_sequence': ('AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCG'
                                   'CCCCGGGCGGGGC'),
                'target_sequence': ('CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCC'
                                    'GGGCGGGGC')
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
                'query_sequence': ('AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCG'
                                   'CCCCGGGCGGGGC'),
                'target_sequence': ('CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCC'
                                    'GGGCGGGGC')
            }, True),
            ({
                'query_sequence': ('AGTCACGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCG'
                                   'CCCCGGGCGGGGC'),
                'target_sequence': ('CGCGCGCCGCCGGGGGGCCGGCCGGCGCCGGGGGGCGCCCC'
                                    'GGGCGGGGC')
            }, False)
        ]
        for expected, z in expected_alignments:
            query = StripedSmithWaterman(expected['query_sequence'],
                                         zero_index=z)
            alignment = query(expected['target_sequence'])
            alignment.set_zero_based(not z)
            self.assertEqual(not z, alignment.is_zero_based())

    def test__get_aligned_sequences(self):
        generic_sequence = "123456789abcdefghijklmnopqrstuvwxyz"
        tests = [  # `end_after_cigar` is how far end extends beyond the cigar.
                   #  Negative values on this should not be possible with SSW
            {
                'cigar_tuples': [
                    (4, 'M'), (3, 'I'), (1, 'D'), (15, 'M')
                ],
                'begin': 4,
                'end_after_cigar': 2,
                'gap_type': 'I',
                'expected': "5678---9abcdefghijklmnopq"
            },
            {
                'cigar_tuples': [
                    (12, 'M')
                ],
                'begin': 10,
                'end_after_cigar': 0,
                'gap_type': 'D',
                'expected': "bcdefghijklm"
            },
            {
                'cigar_tuples': [
                    (10, 'D'), (1, 'M'), (3, 'I'), (2, 'M')
                ],
                'begin': 0,
                'end_after_cigar': 5,
                'gap_type': 'I',
                'expected': "123456789ab---cdefghi"
            },
            {
                'cigar_tuples': [
                    (10, 'D'), (1, 'M'), (3, 'I'), (2, 'M')
                ],
                'begin': 3,
                'end_after_cigar': 0,
                'gap_type': 'D',
                'expected': "----------456789"
            },
            {
                'cigar_tuples': [
                    (1, 'I'), (4, 'M'), (3, 'I'), (1, 'D'), (8, 'M'), (8, 'D'),
                    (2, 'I'), (6, 'M'), (1, 'I')
                ],
                'begin': 4,
                'end_after_cigar': 3,
                'gap_type': 'I',
                'expected': "-5678---9abcdefghijklmnop--qrstuv-wxy"
            }
        ]
        for test in tests:
            mock_object = self.mock_object_factory({})
            # Because SSW's output is [a, b] and Python's list ranges use
            # [a, b) a 1 is added in the calculation of aligned sequences.
            # We just have to subtract 1 while we are testing with the easy to
            # verify interface of `end_after_cigar` to cancel this range effect
            # out.
            end = test['end_after_cigar'] - 1 + test['begin'] + \
                sum(le if t != test['gap_type'] else 0
                    for le, t in test['cigar_tuples'])
            self.assertEqual(test['expected'],
                             AlignmentStructure._get_aligned_sequence(
                                 mock_object, generic_sequence,
                                 test['cigar_tuples'], test['begin'],
                                 end, test['gap_type']))

    def test_aligned_query_target_sequence(self):
        query = StripedSmithWaterman("AGGGTAATTAGGCGTGTTCACCTA")
        alignment = query("AGTCGAAGGGTAATATAGGCGTGTCACCTA")
        self.assertEqual("AGGGTAATATAGGCGTG-TCACCTA",
                         alignment.aligned_target_sequence)
        self.assertEqual("AGGGTAAT-TAGGCGTGTTCACCTA",
                         alignment.aligned_query_sequence)

    def test_aligned_query_target_sequence_with_suppressed_sequences(self):
        query = StripedSmithWaterman("AGGGTAATTAGGCGTGTTCACCTA",
                                     suppress_sequences=True)
        alignment = query("AGTCGAAGGGTAATATAGGCGTGTCACCTA")
        self.assertEqual(None, alignment.aligned_target_sequence)
        self.assertEqual(None, alignment.aligned_query_sequence)


if __name__ == '__main__':
    main()
