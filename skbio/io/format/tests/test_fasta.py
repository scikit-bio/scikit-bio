# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import copy
import io
import string
from unittest import TestCase, main
from functools import partial

import numpy as np

from skbio import Sequence, DNA, RNA, Protein, TabularMSA
from skbio.io import FASTAFormatError, QUALFormatError
from skbio.io.format.fasta import (
    _fasta_sniffer, _fasta_to_generator, _fasta_to_sequence,
    _fasta_to_dna, _fasta_to_rna, _fasta_to_protein,
    _fasta_to_tabular_msa, _generator_to_fasta,
    _sequence_to_fasta, _dna_to_fasta, _rna_to_fasta, _protein_to_fasta,
    _tabular_msa_to_fasta)
from skbio.sequence import GrammaredSequence
from skbio.util import get_data_path
from skbio.util import classproperty
from skbio.util._decorator import overrides


class CustomSequence(GrammaredSequence):
    @classproperty
    @overrides(GrammaredSequence)
    def gap_chars(cls):
        return set('-.')

    @classproperty
    @overrides(GrammaredSequence)
    def default_gap_char(cls):
        return '-'

    @classproperty
    @overrides(GrammaredSequence)
    def definite_chars(cls):
        return set(string.ascii_letters)

    @classproperty
    @overrides(GrammaredSequence)
    def degenerate_map(cls):
        return {}


class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
            'fasta_5_blanks_start_of_file',
            'fasta_5_ws_lines_start_of_file',
            'fasta_blanks_end_of_file',
            'fasta_ws_lines_end_of_file',
            'fasta_blank_lines_between_records',
            'fasta_3_seqs_defaults',
            'fasta_max_width_1',
            'fasta_single_bio_seq_non_defaults',
            'fasta_single_prot_seq_non_defaults',
            'fasta_3_seqs_non_defaults',
            'fasta_max_width_5',
            'fasta_single_dna_seq_defaults',
            'fasta_single_rna_seq_defaults',
            'fasta_description_newline_replacement_empty_str',
            'fasta_multi_seq',
            'fasta_single_dna_seq_non_defaults',
            'fasta_single_rna_seq_non_defaults',
            'fasta_description_newline_replacement_multi_char',
            'fasta_prot_seqs_odd_labels',
            'fasta_single_seq',
            'fasta_id_whitespace_replacement_empty_str',
            'fasta_tabular_msa_different_type',
            'fasta_id_whitespace_replacement_multi_char',
            'fasta_single_bio_seq_defaults',
            'fasta_single_prot_seq_defaults',
            'fasta_10_seqs',
            'fasta_invalid_after_10_seqs',
            'fasta_mixed_qual_scores',
            'qual_3_seqs_non_defaults'
        ]))

        self.negative_fps = list(map(get_data_path, [
            'empty',
            'whitespace_only',
            'fasta_invalid_missing_header',
            'fasta_invalid_blank_line_after_header',
            'fasta_invalid_blank_sequence',
            'fasta_invalid_blank_line_within_sequence',
            'fasta_invalid_whitespace_only_line_within_sequence',
            'fasta_invalid_whitespace_line_after_header',
            'fasta_invalid_missing_seq_data_first',
            'fasta_invalid_missing_seq_data_middle',
            'fasta_invalid_missing_seq_data_last',
            'fasta_invalid_legacy_format',
            'fasta_invalid_whitespace_only_sequence',
            'fasta_id_whitespace_replacement_none',
            'fasta_description_newline_replacement_none',
            'fasta_6_blanks_start_of_file',
            'fasta_6_ws_lines_start_of_file',
            'qual_2_seqs_defaults',
            'qual_3_seqs_defaults',
            'qual_3_seqs_defaults_desc_mismatch',
            'qual_3_seqs_defaults_extra',
            'qual_3_seqs_defaults_id_mismatch',
            'qual_3_seqs_defaults_length_mismatch',
            'qual_description_newline_replacement_empty_str',
            'qual_description_newline_replacement_multi_char',
            'qual_description_newline_replacement_none',
            'qual_id_whitespace_replacement_empty_str',
            'qual_id_whitespace_replacement_multi_char',
            'qual_id_whitespace_replacement_none',
            'qual_invalid_blank_line_within_seq',
            'qual_invalid_legacy_format',
            'qual_invalid_missing_header',
            'qual_invalid_missing_qual_scores_first',
            'qual_invalid_missing_qual_scores_last',
            'qual_invalid_missing_qual_scores_middle',
            'qual_invalid_whitespace_line_in_seq',
            'qual_invalid_blank_line_after_header',
            'qual_invalid_blank_sequence',
            'qual_invalid_whitespace_only_sequence',
            'qual_invalid_ws_line_after_header',
            'qual_invalid_qual_scores_float',
            'qual_invalid_qual_scores_string',
            'qual_max_width_1',
            'qual_max_width_5',
            'qual_multi_seq',
            'qual_multi_seq_roundtrip',
            'qual_prot_seqs_odd_labels',
            'qual_tabular_msa_different_type',
            'qual_single_bio_seq_non_defaults',
            'qual_single_dna_seq_non_defaults',
            'qual_single_prot_seq_non_defaults',
            'qual_single_rna_seq_non_defaults',
            'qual_single_seq',
            'qual_ws_lines_between_records',
            'qual_blank_lines_between_records',
            'qual_5_blanks_start_of_file',
            'qual_5_ws_lines_start_of_file',
            'qual_6_blanks_start_of_file',
            'qual_6_ws_lines_start_of_file',
            'qual_blanks_end_of_file',
            'qual_ws_lines_end_of_file'
        ]))

    def test_positives(self):
        for fp in self.positive_fps:
            self.assertEqual(_fasta_sniffer(fp), (True, {}))

    def test_negatives(self):
        for fp in self.negative_fps:
            self.assertEqual(_fasta_sniffer(fp), (False, {}))


class ReaderTests(TestCase):
    def setUp(self):
        # each structure stores the sequence generator results (expanded into a
        # list) that we expect to obtain from reading, matched with kwargs to
        # pass to the reader, and fasta and qual filepaths that should
        # deserialize into the expected generator results

        # empty file shouldn't yield sequences
        self.empty = ([], {}, list(map(get_data_path, ['empty',
                                                       'whitespace_only'])),
                      list(map(get_data_path, ['empty', 'whitespace_only'])))

        # single sequence
        self.single = (
            [Sequence(
                'ACGT-acgt.', metadata={'id': 'seq1', 'description': 'desc1'},
                positional_metadata={'quality':
                                     np.asarray([10, 20, 30, 10, 0, 0, 0, 255,
                                                 1, 255], dtype=np.uint8)})],
            {},
            list(map(get_data_path, ['fasta_single_seq',
                                     'fasta_max_width_1'])),
            list(map(get_data_path, ['qual_single_seq', 'qual_max_width_1']))
        )

        # multiple sequences
        self.multi = (
            [Sequence(
                'ACGT-acgt.', metadata={'id': 'seq1', 'description': 'desc1'},
                positional_metadata={'quality':
                                     np.asarray([10, 20, 30, 10, 0, 0, 0, 255,
                                                 1, 255], dtype=np.uint8)}),
             Sequence('A', metadata={'id': '_____seq__2_', 'description': ''},
                      positional_metadata={'quality':
                                           np.asarray([42], dtype=np.uint8)}),
             Sequence(
                'AACGGuA', metadata={'id': '', 'description': 'desc3'},
                positional_metadata={'quality':
                                     np.asarray([0, 0, 0, 0, 0, 0, 0],
                                                dtype=np.uint8)}),
             Sequence(
                'ACGTTGCAccGG',
                metadata={'id': '', 'description': ''},
                positional_metadata={'quality':
                                     np.asarray([55, 10, 0, 99, 1, 1, 8, 77,
                                                 40, 10, 10, 0],
                                                dtype=np.uint8)}),
             Sequence('ACGUU',
                      metadata={'id': '', 'description': ''},
                      positional_metadata={'quality':
                                           np.asarray([10, 9, 8, 7, 6],
                                                      dtype=np.uint8)}),
             Sequence(
                 'pQqqqPPQQQ',
                 metadata={'id': 'proteinseq',
                           'description':
                               'detailed description \t\twith  new  lines'},
                 positional_metadata={'quality':
                                      np.asarray([42, 42, 255, 255, 42, 42, 42,
                                                  42, 42, 43],
                                                 dtype=np.uint8)})],
            {},
            list(map(get_data_path, ['fasta_multi_seq', 'fasta_max_width_5',
                                     'fasta_blank_lines_between_records',
                                     'fasta_ws_lines_between_records',
                                     'fasta_5_blanks_start_of_file',
                                     'fasta_5_ws_lines_start_of_file',
                                     'fasta_6_blanks_start_of_file',
                                     'fasta_6_ws_lines_start_of_file',
                                     'fasta_blanks_end_of_file',
                                     'fasta_ws_lines_end_of_file'])),
            list(map(get_data_path, ['qual_multi_seq', 'qual_max_width_5',
                                     'qual_blank_lines_between_records',
                                     'qual_ws_lines_between_records',
                                     'qual_5_blanks_start_of_file',
                                     'qual_5_ws_lines_start_of_file',
                                     'qual_6_blanks_start_of_file',
                                     'qual_6_ws_lines_start_of_file',
                                     'qual_blanks_end_of_file',
                                     'qual_ws_lines_end_of_file']))

        )

        # test constructor parameter, as well as odd labels (label only
        # containing whitespace, label description preceded by multiple spaces,
        # no id) and leading/trailing whitespace on sequence data. for qual
        # files, in addition to the odd labels, test leading/trailing
        # whitespace on qual scores, as well as strange number formatting.
        # also test that fasta and qual headers do not need to match
        # exactly, only that they need to match exactly after parsing (e.g.,
        # after stripping leading/trailing whitespace from descriptions)
        self.odd_labels_different_type = (
            [Protein('DEFQfp',
                     metadata={'id': '', 'description': ''},
                     positional_metadata={'quality':
                                          np.asarray([0, 0, 1, 5, 44, 0],
                                                     dtype=np.uint8)},
                     validate=False),
             Protein(
                 'SKBI', metadata={'id': '', 'description': 'skbio'},
                 positional_metadata={'quality':
                                      np.asarray([1, 2, 33, 123],
                                                 dtype=np.uint8)})],
            {'constructor': partial(Protein, validate=False)},
            list(map(get_data_path, ['fasta_prot_seqs_odd_labels'])),
            list(map(get_data_path, ['qual_prot_seqs_odd_labels']))
        )

        # sequences that can be loaded into a TabularMSA
        self.tabular_msa_different_type = (
            [RNA('aUG',
                 metadata={'id': '', 'description': ''},
                 positional_metadata={'quality':
                                      np.asarray([20, 20, 21],
                                                 dtype=np.uint8)},
                 lowercase='introns'),
             RNA('AuC',
                 metadata={'id': 'rnaseq-1', 'description': 'rnaseq desc 1'},
                 positional_metadata={'quality':
                                      np.asarray([10, 9, 10], dtype=np.uint8)},
                 lowercase='introns'),
             RNA('AUg',
                 metadata={'id': 'rnaseq-2', 'description': 'rnaseq desc 2'},
                 positional_metadata={'quality':
                                      np.asarray([9, 99, 99], dtype=np.uint8)},
                 lowercase='introns')],
            {'constructor': partial(RNA, lowercase='introns')},
            list(map(get_data_path,
                     ['fasta_tabular_msa_different_type'])),
            list(map(get_data_path,
                     ['qual_tabular_msa_different_type']))
        )

        self.lowercase_seqs = (
            [DNA('TAcg',
                 metadata={'id': 'f-o-o', 'description': 'b_a_r'},
                 positional_metadata={'quality':
                                      np.asarray([0, 1, 2, 3],
                                                 dtype=np.uint8)},
                 lowercase='introns')],
            {'constructor': DNA, 'lowercase': 'introns'},
            list(map(get_data_path,
                     ['fasta_single_dna_seq_non_defaults'])),
            list(map(get_data_path,
                     ['qual_single_dna_seq_non_defaults']))
        )

        # store fasta filepath, kwargs, error type, and expected error message
        # for invalid input.
        #
        # note: there is some duplication in testing that fasta and qual
        # parsers raise expected errors. even though the parsers share the same
        # underlying logic, these tests are here as a safeguard in case the
        # code is refactored in the future such that fasta and qual have
        # different implementations (e.g., if qual is written in cython while
        # fasta remains in python)
        self.invalid_fps = list(map(lambda e: (get_data_path(e[0]),
                                               e[1], e[2], e[3]), [
            # fasta and qual missing header
            ('fasta_invalid_missing_header', {}, FASTAFormatError,
             r'non-header.*1st'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_missing_header')},
             QUALFormatError, r'non-header.*1st'),

            # fasta and qual with blank line within sequence
            ('fasta_invalid_blank_line_within_sequence', {}, FASTAFormatError,
             r'whitespace-only'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_blank_line_within_seq')},
             QUALFormatError, r'whitespace-only'),

            # fasta and qual with blank after header
            ('fasta_invalid_blank_sequence', {}, FASTAFormatError,
             r'without sequence data'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_blank_sequence')},
             QUALFormatError, r'without quality scores'),

            # fasta and qual with whitespace only sequence
            ('fasta_invalid_whitespace_only_sequence', {}, FASTAFormatError,
             r'without sequence data'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_whitespace_only_sequence')},
             QUALFormatError, r'without quality scores'),

            # fasta and qual with blank line within sequence
            ('fasta_invalid_blank_line_after_header', {}, FASTAFormatError,
             r'whitespace-only'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_blank_line_after_header')},
             QUALFormatError, r'whitespace-only'),

            # fasta and qual with whitespace-only line within sequence
            ('fasta_invalid_whitespace_only_line_within_sequence',
             {}, FASTAFormatError, r'whitespace-only'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_whitespace_line_in_seq')},
             QUALFormatError, r'whitespace-only'),

            # fasta and qual with whitespace-only line after header
            ('fasta_invalid_whitespace_line_after_header',
             {}, FASTAFormatError, r'whitespace-only'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_ws_line_after_header')},
             QUALFormatError, r'whitespace-only'),

            # fasta and qual missing record data (first record)
            ('fasta_invalid_missing_seq_data_first', {}, FASTAFormatError,
             r'without sequence data'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_missing_qual_scores_first')},
             QUALFormatError, r'without quality scores'),

            # fasta and qual missing record data (middle record)
            ('fasta_invalid_missing_seq_data_middle', {}, FASTAFormatError,
             r'without sequence data'),
            ('fasta_3_seqs_defaults',
             {'qual':
              get_data_path('qual_invalid_missing_qual_scores_middle')},
             QUALFormatError, r'without quality scores'),

            # fasta and qual missing record data (last record)
            ('fasta_invalid_missing_seq_data_last', {}, FASTAFormatError,
             r'without sequence data'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_missing_qual_scores_last')},
             QUALFormatError, r'without quality scores'),

            # fasta and qual in legacy format (;)
            ('fasta_invalid_legacy_format', {}, FASTAFormatError,
             r'non-header.*1st'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_legacy_format')},
             QUALFormatError, r'non-header.*1st'),

            # qual file with an extra record
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_3_seqs_defaults_extra')},
             FASTAFormatError, r'QUAL file has more'),

            # fasta file with an extra record
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_2_seqs_defaults')},
             FASTAFormatError, r'FASTA file has more'),

            # id mismatch between fasta and qual
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_3_seqs_defaults_id_mismatch')},
             FASTAFormatError,
             r'IDs do not match.*\'s_e_q_2\' != \'s_e_q_42\''),

            # description mismatch between fasta and qual
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_3_seqs_defaults_desc_mismatch')},
             FASTAFormatError,
             r'Descriptions do not match.*\'desc 2\' != \'desc 42\''),

            # sequence and quality score length mismatch between fasta and qual
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_3_seqs_defaults_length_mismatch')},
             ValueError,
             r'Number of positional metadata values \(3\) must match the '
             r'positional metadata axis length \(4\)\.'),

            # invalid qual scores (string value can't be converted to integer)
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_qual_scores_string')},
             QUALFormatError,
             r'quality scores to integers:\n100 0 1a -42'),

            # invalid qual scores (float value can't be converted to integer)
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_qual_scores_float')},
             QUALFormatError,
             r'quality scores to integers:\n42    41.0 39 40'),

            # invalid qual scores (negative integer)
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_qual_scores_negative')},
             QUALFormatError,
             r'Quality scores must be greater than or equal to zero\.'),

            # invalid qual scores (over 255)
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_qual_scores_over_255')},
             QUALFormatError,
             r'quality score\(s\) greater than 255'),

            # misc. invalid files used elsewhere in the tests
            ('fasta_invalid_after_10_seqs', {}, FASTAFormatError,
             r'without sequence data'),
            ('fasta_id_whitespace_replacement_none', {}, FASTAFormatError,
             r'whitespace-only'),
            ('fasta_description_newline_replacement_none', {},
             FASTAFormatError, r'whitespace-only')
        ]))

    # extensive tests for fasta -> generator reader since it is used by all
    # other fasta -> object readers

    def test_fasta_to_generator_valid_files(self):
        test_cases = (self.empty, self.single, self.multi,
                      self.odd_labels_different_type,
                      self.tabular_msa_different_type,
                      self.lowercase_seqs)

        # Strategy:
        #   for each fasta file, read it without its corresponding qual file,
        #   and ensure observed vs. expected match, ignoring quality scores in
        #   expected. next, parse the current fasta file with each
        #   corresponding quality file and ensure that observed vs. expected
        #   match, this time taking quality scores into account. this
        #   sufficiently exercises parsing a standalone fasta file and paired
        #   fasta/qual files
        for exp, kwargs, fasta_fps, qual_fps in test_cases:
            for fasta_fp in fasta_fps:
                obs = list(_fasta_to_generator(fasta_fp, **kwargs))
                self.assertEqual(len(obs), len(exp))
                for o, e in zip(obs, exp):
                    e = copy.copy(e)
                    del e.positional_metadata['quality']
                    self.assertEqual(o, e)

                for qual_fp in qual_fps:
                    obs = list(_fasta_to_generator(fasta_fp, qual=qual_fp,
                                                   **kwargs))

                    self.assertEqual(len(obs), len(exp))
                    for o, e in zip(obs, exp):
                        self.assertEqual(o, e)

    def test_fasta_to_generator_invalid_files(self):
        for fp, kwargs, error_type, error_msg_regex in self.invalid_fps:
            with self.assertRaisesRegex(error_type, error_msg_regex):
                list(_fasta_to_generator(fp, **kwargs))

    # light testing of fasta -> object readers to ensure interface is present
    # and kwargs are passed through. extensive testing of underlying reader is
    # performed above

    def test_fasta_to_any_sequence(self):
        for constructor, reader_fn in ((Sequence,
                                        _fasta_to_sequence),
                                       (partial(DNA, validate=False,
                                                lowercase='introns'),
                                        partial(_fasta_to_dna,
                                                validate=False,
                                                lowercase='introns')),
                                       (partial(RNA, validate=False,
                                                lowercase='introns'),
                                        partial(_fasta_to_rna,
                                                validate=False,
                                                lowercase='introns')),
                                       (partial(Protein, lowercase='introns'),
                                        partial(_fasta_to_protein,
                                                validate=False,
                                                lowercase='introns'))):

            # empty file
            empty_fp = get_data_path('empty')
            with self.assertRaisesRegex(ValueError, r'1st sequence'):
                reader_fn(empty_fp)
            with self.assertRaisesRegex(ValueError, r'1st sequence'):
                reader_fn(empty_fp, qual=empty_fp)

            # the sequences in the following files don't necessarily make sense
            # for each of the sequence object types that they're read into
            # (e.g., reading a protein sequence into a dna sequence object).
            # however, for the purposes of testing the various
            # fasta -> sequence readers, this works out okay as it is valid to
            # construct a sequence object with invalid characters. we're
            # interested in testing the reading logic here, and don't care so
            # much about constructing semantically-meaningful/valid sequence
            # objects

            # file with only 1 seq, get first
            fasta_fps = list(map(get_data_path,
                                 ['fasta_single_seq', 'fasta_max_width_1']))
            for fasta_fp in fasta_fps:
                exp = constructor(
                    'ACGT-acgt.',
                    metadata={'id': 'seq1', 'description': 'desc1'})

                obs = reader_fn(fasta_fp)
                self.assertEqual(obs, exp)

                exp.positional_metadata.insert(
                    0, 'quality',
                    np.asarray([10, 20, 30, 10, 0, 0, 0, 255, 1, 255],
                               dtype=np.uint8))
                qual_fps = list(map(get_data_path,
                                    ['qual_single_seq', 'qual_max_width_1']))
                for qual_fp in qual_fps:
                    obs = reader_fn(fasta_fp, qual=qual_fp)
                    self.assertEqual(obs, exp)

            # file with multiple seqs
            fasta_fps = list(map(get_data_path,
                                 ['fasta_multi_seq', 'fasta_max_width_5']))
            qual_fps = list(map(get_data_path,
                                ['qual_multi_seq', 'qual_max_width_5']))
            for fasta_fp in fasta_fps:
                # get first
                exp = constructor(
                    'ACGT-acgt.',
                    metadata={'id': 'seq1', 'description': 'desc1'})

                obs = reader_fn(fasta_fp)
                self.assertEqual(obs, exp)

                exp.positional_metadata.insert(
                    0, 'quality',
                    np.asarray([10, 20, 30, 10, 0, 0, 0, 255, 1, 255],
                               dtype=np.uint8))
                for qual_fp in qual_fps:
                    obs = reader_fn(fasta_fp, qual=qual_fp)
                    self.assertEqual(obs, exp)

                # get middle
                exp = constructor('ACGTTGCAccGG',
                                  metadata={'id': '', 'description': ''})

                obs = reader_fn(fasta_fp, seq_num=4)
                self.assertEqual(obs, exp)

                exp.positional_metadata.insert(
                    0, 'quality',
                    np.asarray([55, 10, 0, 99, 1, 1, 8, 77, 40, 10, 10, 0],
                               dtype=np.uint8))
                for qual_fp in qual_fps:
                    obs = reader_fn(fasta_fp, seq_num=4, qual=qual_fp)
                    self.assertEqual(obs, exp)

                # get last
                exp = constructor(
                    'pQqqqPPQQQ',
                    metadata={'id': 'proteinseq',
                              'description':
                                  'detailed description \t\twith  new  lines'})

                obs = reader_fn(fasta_fp, seq_num=6)
                self.assertEqual(obs, exp)

                exp.positional_metadata.insert(
                    0, 'quality',
                    np.asarray([42, 42, 255, 255, 42, 42, 42, 42, 42, 43],
                               dtype=np.uint8))
                for qual_fp in qual_fps:
                    obs = reader_fn(fasta_fp, seq_num=6, qual=qual_fp)
                    self.assertEqual(obs, exp)

                # seq_num too large
                with self.assertRaisesRegex(ValueError, r'8th sequence'):
                    reader_fn(fasta_fp, seq_num=8)
                for qual_fp in qual_fps:
                    with self.assertRaisesRegex(ValueError, r'8th sequence'):
                        reader_fn(fasta_fp, seq_num=8, qual=qual_fp)

                # seq_num too small
                with self.assertRaisesRegex(ValueError, r'`seq_num`=0'):
                    reader_fn(fasta_fp, seq_num=0)
                for qual_fp in qual_fps:
                    with self.assertRaisesRegex(ValueError, r'`seq_num`=0'):
                        reader_fn(fasta_fp, seq_num=0, qual=qual_fp)

    def test_fasta_to_tabular_msa(self):
        test_cases = (self.empty, self.single,
                      self.tabular_msa_different_type,
                      self.lowercase_seqs)

        # see comment in test_fasta_to_generator_valid_files (above) for
        # testing strategy
        for exp_list, kwargs, fasta_fps, qual_fps in test_cases:
            if 'constructor' not in kwargs:
                kwargs['constructor'] = CustomSequence
                exp_list = [CustomSequence(seq) for seq in exp_list]

            exp = TabularMSA(exp_list)

            for fasta_fp in fasta_fps:
                obs = _fasta_to_tabular_msa(fasta_fp, **kwargs)

                self.assertEqual(len(obs), len(exp))
                for o, e in zip(obs, exp):
                    e = copy.copy(e)
                    del e.positional_metadata['quality']
                    self.assertEqual(o, e)

                for qual_fp in qual_fps:
                    obs = _fasta_to_tabular_msa(fasta_fp, qual=qual_fp,
                                                **kwargs)
                    self.assertEqual(obs, exp)

    def test_fasta_to_tabular_msa_no_constructor(self):
        with self.assertRaisesRegex(ValueError, r'`constructor`'):
            _fasta_to_tabular_msa(get_data_path('fasta_single_seq'))


class WriterTests(TestCase):
    def setUp(self):
        self.bio_seq1 = DNA(
            'ACGT-acgt.',
            metadata={'id': 'seq1', 'description': 'desc1'},
            positional_metadata={'quality': [10, 20, 30, 10, 0, 0, 0, 255,
                                             1, 255]},
            lowercase='introns')
        self.bio_seq2 = DNA(
            'A',
            metadata={'id': ' \n  \nseq \t2 '},
            positional_metadata={'quality': [42]},
            lowercase='introns')
        self.bio_seq3 = RNA(
            'AACGGuA',
            metadata={'description': 'desc3'},
            positional_metadata={'quality': [0, 0, 0, 0, 0, 0, 0]},
            lowercase='introns')
        self.dna_seq = DNA(
            'ACGTTGCAccGG',
            positional_metadata={'quality': [55, 10, 0, 99, 1, 1, 8, 77, 40,
                                             10, 10, 0]},
            lowercase='introns')
        self.rna_seq = RNA('ACGUU',
                           positional_metadata={'quality': [10, 9, 8, 7, 6]},
                           lowercase='introns')
        self.prot_seq = Protein(
            'pQqqqPPQQQ',
            metadata={'id': 'proteinseq',
                      'description': "\ndetailed\ndescription \t\twith "
                                     " new\n\nlines\n\n\n"},
            positional_metadata={'quality': [42, 42, 255, 255, 42, 42, 42, 42,
                                             42, 43]},
            lowercase='introns')

        seqs = [
            CustomSequence(
                'UUUU',
                metadata={'id': 's\te\tq\t1', 'description': 'desc\n1'},
                positional_metadata={'quality': [1234, 0, 0, 2]},
                lowercase='introns'),
            CustomSequence(
                'CATC',
                metadata={'id': 's\te\tq\t2', 'description': 'desc\n2'},
                positional_metadata={'quality': [1, 11, 111, 11112]}),
            CustomSequence(
                'sits',
                metadata={'id': 's\te\tq\t3', 'description': 'desc\n3'},
                positional_metadata={'quality': [12345, 678909, 999999,
                                                 4242424242]})
        ]
        self.msa = TabularMSA(seqs)

        def empty_gen():
            yield from ()

        def single_seq_gen():
            yield self.bio_seq1

        # generate sequences with descriptions containing newlines (to test
        # description_newline_replacement)
        def newline_description_gen():
            yield self.prot_seq
            yield DNA('AGGAGAATA',
                      metadata={'id': 'foo', 'description': '\n\n\n\n'},
                      positional_metadata={'quality': range(9)},
                      lowercase='introns')

        # generate sequences with ids containing whitespace (to test
        # id_whitespace_replacement)
        def whitespace_id_gen():
            yield self.bio_seq2
            yield RNA('UA', metadata={'id': '\n\t \t', 'description': 'a\nb'},
                      positional_metadata={'quality': [1000, 1]})

        # multiple sequences of mixed types, lengths, and metadata. lengths are
        # chosen to exercise various splitting cases when testing max_width,
        # including exercising the different splitting algorithms used for
        # sequence data vs. quality scores
        def multi_seq_gen():
            yield from (self.bio_seq1, self.bio_seq2, self.bio_seq3,
                        self.dna_seq, self.rna_seq, self.prot_seq)

        # can be serialized if no qual file is provided, else it should raise
        # an error because one seq has qual scores and the other doesn't
        def mixed_qual_score_gen():
            yield self.bio_seq1
            yield DNA('AAAAT',
                      metadata={'id': 'da,dadadada',
                                'description': '10 hours'},
                      lowercase='introns')

        self.mixed_qual_score_gen = mixed_qual_score_gen()

        # store sequence generator to serialize, writer kwargs (if any), and
        # fasta and qual filepaths of expected results
        self.objs_fps = list(map(lambda e: (e[0], e[1], get_data_path(e[2]),
                                            get_data_path(e[3])), [
            (empty_gen(), {}, 'empty', 'empty'),
            (single_seq_gen(), {'lowercase': 'introns'}, 'fasta_single_seq',
             'qual_single_seq'),

            # no splitting of sequence or qual data across lines b/c max_width
            # is sufficiently large
            (single_seq_gen(), {'max_width': 32, 'lowercase': 'introns'},
             'fasta_single_seq',
             'qual_single_seq'),

            # splitting algorithm for sequence and qual scores is different;
            # make sure individual qual scores aren't split across lines even
            # if they exceed max_width
            (single_seq_gen(), {'max_width': 1, 'lowercase': 'introns'},
             'fasta_max_width_1',
             'qual_max_width_1'),
            (multi_seq_gen(),
             {'lowercase': 'introns'}, 'fasta_multi_seq', 'qual_multi_seq'),
            (multi_seq_gen(),
             {'max_width': 5, 'lowercase': 'introns'}, 'fasta_max_width_5',
             'qual_max_width_5'),
            (newline_description_gen(),
             {'description_newline_replacement': ':-)',
              'lowercase': 'introns'},
             'fasta_description_newline_replacement_multi_char',
             'qual_description_newline_replacement_multi_char'),
            (newline_description_gen(),
             {'description_newline_replacement': '',
              'lowercase': 'introns'},
             'fasta_description_newline_replacement_empty_str',
             'qual_description_newline_replacement_empty_str',),
            (newline_description_gen(),
             {'description_newline_replacement': None,
              'lowercase': 'introns'},
             'fasta_description_newline_replacement_none',
             'qual_description_newline_replacement_none'),
            (whitespace_id_gen(),
             {'id_whitespace_replacement': '>:o'},
             'fasta_id_whitespace_replacement_multi_char',
             'qual_id_whitespace_replacement_multi_char'),
            (whitespace_id_gen(),
             {'id_whitespace_replacement': ''},
             'fasta_id_whitespace_replacement_empty_str',
             'qual_id_whitespace_replacement_empty_str'),
            (whitespace_id_gen(),
             {'id_whitespace_replacement': None},
             'fasta_id_whitespace_replacement_none',
             'qual_id_whitespace_replacement_none'),
        ]))

        def blank_seq_gen():
            yield from (self.bio_seq1, Sequence(''))

        # generators or parameter combos that cannot be written in fasta
        # format, paired with kwargs (if any), error type, and expected error
        # message regexp
        self.invalid_objs = [
            (blank_seq_gen(), {}, ValueError, r'2nd.*empty'),
            (single_seq_gen(),
             {'max_width': 0}, ValueError, r'max_width=0'),
            (multi_seq_gen(), {'id_whitespace_replacement': '-\n_'},
             ValueError, r'Newline character'),
            (multi_seq_gen(), {'description_newline_replacement': '-.-\n'},
             ValueError, r'Newline character'),
            (mixed_qual_score_gen(), {'qual': io.StringIO()}, ValueError,
             r'2nd sequence.*does not have quality scores')
        ]

    # extensive tests for generator -> fasta writer since it is used by all
    # other object -> fasta writers

    def test_generator_to_fasta_no_qual(self):
        # test writing standalone fasta (i.e., without a qual file)
        for obj, kwargs, fp, _ in self.objs_fps:
            fh = io.StringIO()
            _generator_to_fasta(obj, fh, **kwargs)
            obs = fh.getvalue()
            fh.close()

            with io.open(fp) as fh:
                exp = fh.read()
            self.assertEqual(obs, exp)

    def test_generator_to_fasta_mixed_qual_scores(self):
        # test writing some sequences with qual scores and some without is
        # possible if no qual output file is specified
        fh = io.StringIO()
        _generator_to_fasta(self.mixed_qual_score_gen, fh, lowercase='introns')
        obs = fh.getvalue()
        fh.close()

        with io.open(get_data_path('fasta_mixed_qual_scores')) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_generator_to_fasta_with_qual(self):
        # test writing fasta and qual files
        for obj, kwargs, fasta_fp, qual_fp in self.objs_fps:
            if qual_fp is not None:
                fasta_fh = io.StringIO()
                qual_fh = io.StringIO()
                _generator_to_fasta(obj, fasta_fh, qual=qual_fh, **kwargs)
                obs_fasta = fasta_fh.getvalue()
                obs_qual = qual_fh.getvalue()
                fasta_fh.close()
                qual_fh.close()

                with io.open(fasta_fp) as fh:
                    exp_fasta = fh.read()
                with io.open(qual_fp) as fh:
                    exp_qual = fh.read()

                self.assertEqual(obs_fasta, exp_fasta)
                self.assertEqual(obs_qual, exp_qual)

    def test_generator_to_fasta_invalid_input(self):
        for obj, kwargs, error_type, error_msg_regexp in self.invalid_objs:
            fh = io.StringIO()
            with self.assertRaisesRegex(error_type, error_msg_regexp):
                _generator_to_fasta(obj, fh, **kwargs)
            fh.close()

    # light testing of object -> fasta writers to ensure interface is present
    # and kwargs are passed through. extensive testing of underlying writer is
    # performed above
    def test_any_sequence_to_fasta(self):
        # store writer function, sequence object to write, expected
        # fasta filepath for default parameters, expected fasta filepath for
        # non-defaults, and expected qual filepath for non-defaults
        id_ = 'f o o'
        desc = 'b\na\nr'
        test_data = (
            (partial(_sequence_to_fasta, lowercase='introns'),
             Sequence('ACgt', metadata={'id': id_, 'description': desc},
                      positional_metadata={'quality': range(1, 5)},
                      lowercase='introns'),
             ('fasta_single_bio_seq_defaults',
              'fasta_single_bio_seq_non_defaults',
              'qual_single_bio_seq_non_defaults')),
            (partial(_dna_to_fasta, lowercase='introns'),
             DNA('TAcg', metadata={'id': id_, 'description': desc},
                 positional_metadata={'quality': range(4)},
                 lowercase='introns'),
             ('fasta_single_dna_seq_defaults',
              'fasta_single_dna_seq_non_defaults',
              'qual_single_dna_seq_non_defaults')),
            (partial(_rna_to_fasta, lowercase='introns'),
             RNA('uaCG', metadata={'id': id_, 'description': desc},
                 positional_metadata={'quality': range(2, 6)},
                 lowercase='introns'),
             ('fasta_single_rna_seq_defaults',
              'fasta_single_rna_seq_non_defaults',
              'qual_single_rna_seq_non_defaults')),
            (partial(_protein_to_fasta, lowercase='introns'),
             Protein('PqQ', metadata={'id': id_, 'description': desc},
                     positional_metadata={'quality': [42, 41, 40]},
                     lowercase='introns'),
             ('fasta_single_prot_seq_defaults',
              'fasta_single_prot_seq_non_defaults',
              'qual_single_prot_seq_non_defaults')))

        for fn, obj, fps in test_data:
            defaults_fp, non_defaults_fasta_fp, non_defaults_qual_fp = fps

            # test writing with default parameters
            fh = io.StringIO()
            fn(obj, fh)
            obs = fh.getvalue()
            fh.close()

            with io.open(get_data_path(defaults_fp)) as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

            # test writing with non-defaults
            fasta_fh = io.StringIO()
            qual_fh = io.StringIO()
            fn(obj, fasta_fh, id_whitespace_replacement='-',
               description_newline_replacement='_', max_width=1, qual=qual_fh)
            obs_fasta = fasta_fh.getvalue()
            obs_qual = qual_fh.getvalue()
            fasta_fh.close()
            qual_fh.close()

            with io.open(get_data_path(non_defaults_fasta_fp)) as fh:
                exp_fasta = fh.read()
            with io.open(get_data_path(non_defaults_qual_fp)) as fh:
                exp_qual = fh.read()

            self.assertEqual(obs_fasta, exp_fasta)
            self.assertEqual(obs_qual, exp_qual)

    def test_any_sequences_to_fasta(self):
        # test writing with default parameters
        fh = io.StringIO()
        _tabular_msa_to_fasta(self.msa, fh)
        obs = fh.getvalue()
        fh.close()

        with io.open(get_data_path('fasta_3_seqs_defaults')) as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

        # test writing with non-defaults
        fasta_fh = io.StringIO()
        qual_fh = io.StringIO()
        _tabular_msa_to_fasta(self.msa, fasta_fh,
                              id_whitespace_replacement='*',
                              description_newline_replacement='+', max_width=3,
                              qual=qual_fh)
        obs_fasta = fasta_fh.getvalue()
        obs_qual = qual_fh.getvalue()
        fasta_fh.close()
        qual_fh.close()

        with io.open(get_data_path('fasta_3_seqs_non_defaults')) as fh:
            exp_fasta = fh.read()
        with io.open(get_data_path('qual_3_seqs_non_defaults')) as fh:
            exp_qual = fh.read()

        self.assertEqual(obs_fasta, exp_fasta)
        self.assertEqual(obs_qual, exp_qual)


class RoundtripTests(TestCase):
    def test_roundtrip_generators(self):
        # test that fasta and qual files can be streamed into memory and back
        # out to disk using generator reader and writer
        fps = list(map(lambda e: list(map(get_data_path, e)),
                       [('empty', 'empty'),
                        ('fasta_multi_seq_roundtrip',
                         'qual_multi_seq_roundtrip')]))

        for fasta_fp, qual_fp in fps:
            with io.open(fasta_fp) as fh:
                exp_fasta = fh.read()
            with io.open(qual_fp) as fh:
                exp_qual = fh.read()

            fasta_fh = io.StringIO()
            qual_fh = io.StringIO()
            _generator_to_fasta(_fasta_to_generator(fasta_fp, qual=qual_fp),
                                fasta_fh, qual=qual_fh)
            obs_fasta = fasta_fh.getvalue()
            obs_qual = qual_fh.getvalue()
            fasta_fh.close()
            qual_fh.close()

            self.assertEqual(obs_fasta, exp_fasta)
            self.assertEqual(obs_qual, exp_qual)

    def test_roundtrip_tabular_msa(self):
        fps = list(map(lambda e: list(map(get_data_path, e)),
                       [('empty', 'empty'),
                        ('fasta_tabular_msa_different_type',
                         'qual_tabular_msa_different_type')]))

        reader = partial(_fasta_to_tabular_msa, constructor=CustomSequence)
        writer = _tabular_msa_to_fasta
        for fasta_fp, qual_fp in fps:
            # read
            obj1 = reader(fasta_fp, qual=qual_fp)

            # write
            fasta_fh = io.StringIO()
            qual_fh = io.StringIO()
            writer(obj1, fasta_fh, qual=qual_fh)
            fasta_fh.seek(0)
            qual_fh.seek(0)

            # read
            obj2 = reader(fasta_fh, qual=qual_fh)
            fasta_fh.close()
            qual_fh.close()

            self.assertEqual(obj1, obj2)

    def test_roundtrip_biological_sequences(self):
        fps = list(map(lambda e: list(map(get_data_path, e)),
                       [('fasta_multi_seq_roundtrip',
                         'qual_multi_seq_roundtrip'),
                        ('fasta_tabular_msa_different_type',
                         'qual_tabular_msa_different_type')]))

        for reader, writer in ((_fasta_to_sequence,
                                _sequence_to_fasta),
                               (partial(_fasta_to_dna,
                                        validate=False),
                                _dna_to_fasta),
                               (partial(_fasta_to_rna,
                                        validate=False),
                                _rna_to_fasta),
                               (partial(_fasta_to_protein,
                                        validate=False),
                                _protein_to_fasta)):
            for fasta_fp, qual_fp in fps:
                # read
                obj1 = reader(fasta_fp, qual=qual_fp)

                # write
                fasta_fh = io.StringIO()
                qual_fh = io.StringIO()
                writer(obj1, fasta_fh, qual=qual_fh)
                fasta_fh.seek(0)
                qual_fh.seek(0)

                # read
                obj2 = reader(fasta_fh, qual=qual_fh)
                fasta_fh.close()
                qual_fh.close()

                self.assertEqual(obj1, obj2)


if __name__ == '__main__':
    main()
