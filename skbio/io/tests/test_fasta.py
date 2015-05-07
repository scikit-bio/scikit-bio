# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import map, range, zip
from six import StringIO

from unittest import TestCase, main
from functools import partial

from skbio import (Sequence, DNA, RNA, Protein, SequenceCollection, Alignment)
from skbio.io import FASTAFormatError
from skbio.io.fasta import (
    _fasta_sniffer, _fasta_to_generator, _fasta_to_biological_sequence,
    _fasta_to_dna_sequence,
    _fasta_to_rna_sequence, _fasta_to_protein_sequence,
    _fasta_to_sequence_collection, _fasta_to_alignment, _generator_to_fasta,
    _biological_sequence_to_fasta,
    _dna_sequence_to_fasta, _rna_sequence_to_fasta, _protein_sequence_to_fasta,
    _sequence_collection_to_fasta, _alignment_to_fasta)
from skbio.util import get_data_path


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
            'fasta_sequence_collection_different_type',
            'fasta_id_whitespace_replacement_multi_char',
            'fasta_single_bio_seq_defaults',
            'fasta_single_prot_seq_defaults',
            'fasta_10_seqs',
            'fasta_invalid_after_10_seqs',
            'fasta_mixed_qual_scores',
            'qual_invalid_qual_scores_float',
            'qual_invalid_qual_scores_string'
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
            'qual_3_seqs_non_defaults',
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
            'qual_max_width_1',
            'qual_max_width_5',
            'qual_multi_seq',
            'qual_multi_seq_roundtrip',
            'qual_prot_seqs_odd_labels',
            'qual_sequence_collection_different_type',
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
                'ACGT-acgt.', id='seq1', description='desc1',
                quality=[10, 20, 30, 10, 0, 0, 0, 88888, 1, 3456])],
            {},
            list(map(get_data_path, ['fasta_single_seq',
                                     'fasta_max_width_1'])),
            list(map(get_data_path, ['qual_single_seq', 'qual_max_width_1']))
        )

        # multiple sequences
        self.multi = (
            [Sequence(
                'ACGT-acgt.', id='seq1', description='desc1',
                quality=[10, 20, 30, 10, 0, 0, 0, 88888, 1, 3456]),
             Sequence('A', id='_____seq__2_', quality=[42]),
             Sequence(
                'AACGGuA', description='desc3', quality=[0, 0, 0, 0, 0, 0, 0]),
             Sequence(
                'ACGTTGCAccGG',
                quality=[55, 10, 0, 999, 1, 1, 8, 775, 40, 10, 10, 0]),
             Sequence('ACGUU', quality=[10, 9, 8, 7, 6]),
             Sequence(
                 'pQqqqPPQQQ', id='proteinseq',
                 description='detailed description \t\twith  new  lines',
                 quality=[42, 42, 442, 442, 42, 42, 42, 42, 42, 43])],
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
            [Protein('DEFQfp', quality=[0, 0, 1, 5, 44, 0], validate=False),
             Protein(
                 'SKBI', description='skbio', quality=[1, 2, 33, 123456789])],
            {'constructor': partial(Protein, validate=False)},
            list(map(get_data_path, ['fasta_prot_seqs_odd_labels'])),
            list(map(get_data_path, ['qual_prot_seqs_odd_labels']))
        )

        # sequences that can be loaded into a SequenceCollection or Alignment.
        # they are also a different type than Sequence in order to
        # exercise the constructor parameter
        self.sequence_collection_different_type = (
            [RNA('AUG', quality=[20, 20, 21]),
             RNA('AUC', id='rnaseq-1', description='rnaseq desc 1',
                 quality=[10, 9, 10]),
             RNA('AUG', id='rnaseq-2', description='rnaseq desc 2',
                 quality=[9, 99, 999])],
            {'constructor': partial(RNA, validate=False)},
            list(map(get_data_path,
                     ['fasta_sequence_collection_different_type'])),
            list(map(get_data_path,
                     ['qual_sequence_collection_different_type']))
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
             'non-header.*1st FASTA'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_missing_header')},
             FASTAFormatError, 'non-header.*1st QUAL'),

            # fasta and qual with blank line within sequence
            ('fasta_invalid_blank_line_within_sequence', {}, FASTAFormatError,
             'whitespace-only.*FASTA'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_blank_line_within_seq')},
             FASTAFormatError, 'whitespace-only.*QUAL'),

            # fasta and qual with blank after header
            ('fasta_invalid_blank_sequence', {}, FASTAFormatError,
             'without sequence data'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_blank_sequence')},
             FASTAFormatError, 'without quality scores'),

            # fasta and qual with whitespace only sequence
            ('fasta_invalid_whitespace_only_sequence', {}, FASTAFormatError,
             'without sequence data'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_whitespace_only_sequence')},
             FASTAFormatError, 'without quality scores'),

            # fasta and qual with blank line within sequence
            ('fasta_invalid_blank_line_after_header', {}, FASTAFormatError,
             'whitespace-only.*FASTA'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_blank_line_after_header')},
             FASTAFormatError, 'whitespace-only.*QUAL'),

            # fasta and qual with whitespace-only line within sequence
            ('fasta_invalid_whitespace_only_line_within_sequence',
             {}, FASTAFormatError, 'whitespace-only.*FASTA'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_whitespace_line_in_seq')},
             FASTAFormatError, 'whitespace-only.*QUAL'),

            # fasta and qual with whitespace-only line after header
            ('fasta_invalid_whitespace_line_after_header',
             {}, FASTAFormatError, 'whitespace-only.*FASTA'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_ws_line_after_header')},
             FASTAFormatError, 'whitespace-only.*QUAL'),

            # fasta and qual missing record data (first record)
            ('fasta_invalid_missing_seq_data_first', {}, FASTAFormatError,
             'without sequence data'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_missing_qual_scores_first')},
             FASTAFormatError, 'without quality scores'),

            # fasta and qual missing record data (middle record)
            ('fasta_invalid_missing_seq_data_middle', {}, FASTAFormatError,
             'without sequence data'),
            ('fasta_3_seqs_defaults',
             {'qual':
              get_data_path('qual_invalid_missing_qual_scores_middle')},
             FASTAFormatError, 'without quality scores'),

            # fasta and qual missing record data (last record)
            ('fasta_invalid_missing_seq_data_last', {}, FASTAFormatError,
             'without sequence data'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_missing_qual_scores_last')},
             FASTAFormatError, 'without quality scores'),

            # fasta and qual in legacy format (;)
            ('fasta_invalid_legacy_format', {}, FASTAFormatError,
             'non-header.*1st FASTA'),
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_legacy_format')},
             FASTAFormatError, 'non-header.*1st QUAL'),

            # qual file with an extra record
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_3_seqs_defaults_extra')},
             FASTAFormatError, 'QUAL file has more'),

            # fasta file with an extra record
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_2_seqs_defaults')},
             FASTAFormatError, 'FASTA file has more'),

            # id mismatch between fasta and qual
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_3_seqs_defaults_id_mismatch')},
             FASTAFormatError,
             'IDs do not match.*\'s_e_q_2\' != \'s_e_q_42\''),

            # description mismatch between fasta and qual
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_3_seqs_defaults_desc_mismatch')},
             FASTAFormatError,
             'Descriptions do not match.*\'desc 2\' != \'desc 42\''),

            # sequence and quality score length mismatch between fasta and qual
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_3_seqs_defaults_length_mismatch')},
             ValueError,
             'Number of quality scores \(3\) must match the number of characte'
             'rs in the sequence \(4\)\.'),

            # invalid qual scores (string value can't be converted to integer)
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_qual_scores_string')},
             FASTAFormatError,
             'quality scores to integers:\n100 0 1a -42'),

            # invalid qual scores (float value can't be converted to integer)
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_qual_scores_float')},
             FASTAFormatError,
             'quality scores to integers:\n42    41.0 39 40'),

            # invalid qual scores (negative integer)
            ('fasta_3_seqs_defaults',
             {'qual': get_data_path('qual_invalid_qual_scores_negative')},
             ValueError,
             'Quality scores must be greater than or equal to zero\.'),

            # misc. invalid files used elsewhere in the tests
            ('fasta_invalid_after_10_seqs', {}, FASTAFormatError,
             'without sequence data'),
            ('fasta_id_whitespace_replacement_none', {}, FASTAFormatError,
             'whitespace-only.*FASTA'),
            ('fasta_description_newline_replacement_none', {},
             FASTAFormatError, 'whitespace-only.*FASTA')
        ]))

    # extensive tests for fasta -> generator reader since it is used by all
    # other fasta -> object readers

    def test_fasta_to_generator_valid_files(self):
        test_cases = (self.empty, self.single, self.multi,
                      self.odd_labels_different_type,
                      self.sequence_collection_different_type)

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
                print(fasta_fp)
                self.assertEqual(len(obs), len(exp))
                for o, e in zip(obs, exp):
                    self.assertTrue(o.equals(e, ignore=['quality']))

                for qual_fp in qual_fps:
                    obs = list(_fasta_to_generator(fasta_fp, qual=qual_fp,
                                                   **kwargs))

                    self.assertEqual(len(obs), len(exp))
                    for o, e in zip(obs, exp):
                        self.assertTrue(o.equals(e))

    def test_fasta_to_generator_invalid_files(self):
        for fp, kwargs, error_type, error_msg_regex in self.invalid_fps:
            with self.assertRaisesRegexp(error_type, error_msg_regex):
                list(_fasta_to_generator(fp, **kwargs))

    # light testing of fasta -> object readers to ensure interface is present
    # and kwargs are passed through. extensive testing of underlying reader is
    # performed above

    def test_fasta_to_any_sequence(self):
        for constructor, reader_fn in ((Sequence,
                                        _fasta_to_biological_sequence),
                                       (partial(DNA, validate=False),
                                        _fasta_to_dna_sequence),
                                       (partial(RNA, validate=False),
                                        _fasta_to_rna_sequence),
                                       (partial(Protein, validate=False),
                                        _fasta_to_protein_sequence)):

            # empty file
            empty_fp = get_data_path('empty')
            with self.assertRaisesRegexp(ValueError, '1st sequence'):
                reader_fn(empty_fp)
            with self.assertRaisesRegexp(ValueError, '1st sequence'):
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
                    'ACGT-acgt.', id='seq1', description='desc1',
                    quality=[10, 20, 30, 10, 0, 0, 0, 88888, 1, 3456])

                obs = reader_fn(fasta_fp)
                self.assertTrue(obs.equals(exp, ignore=['quality']))

                qual_fps = list(map(get_data_path,
                                    ['qual_single_seq', 'qual_max_width_1']))
                for qual_fp in qual_fps:
                    obs = reader_fn(fasta_fp, qual=qual_fp)
                    self.assertTrue(obs.equals(exp))

            # file with multiple seqs
            fasta_fps = list(map(get_data_path,
                                 ['fasta_multi_seq', 'fasta_max_width_5']))
            qual_fps = list(map(get_data_path,
                                ['qual_multi_seq', 'qual_max_width_5']))
            for fasta_fp in fasta_fps:
                # get first
                exp = constructor(
                    'ACGT-acgt.', id='seq1', description='desc1',
                    quality=[10, 20, 30, 10, 0, 0, 0, 88888, 1, 3456])

                obs = reader_fn(fasta_fp)
                self.assertTrue(obs.equals(exp, ignore=['quality']))

                for qual_fp in qual_fps:
                    obs = reader_fn(fasta_fp, qual=qual_fp)
                    self.assertTrue(obs.equals(exp))

                # get middle
                exp = constructor('ACGTTGCAccGG',
                                  quality=[55, 10, 0, 999, 1, 1, 8, 775, 40,
                                           10, 10, 0])

                obs = reader_fn(fasta_fp, seq_num=4)
                self.assertTrue(obs.equals(exp, ignore=['quality']))

                for qual_fp in qual_fps:
                    obs = reader_fn(fasta_fp, seq_num=4, qual=qual_fp)
                    self.assertTrue(obs.equals(exp))

                # get last
                exp = constructor(
                    'pQqqqPPQQQ', id='proteinseq',
                    description='detailed description \t\twith  new  lines',
                    quality=[42, 42, 442, 442, 42, 42, 42, 42, 42, 43])

                obs = reader_fn(fasta_fp, seq_num=6)
                self.assertTrue(obs.equals(exp, ignore=['quality']))

                for qual_fp in qual_fps:
                    obs = reader_fn(fasta_fp, seq_num=6, qual=qual_fp)
                    self.assertTrue(obs.equals(exp))

                # seq_num too large
                with self.assertRaisesRegexp(ValueError, '8th sequence'):
                    reader_fn(fasta_fp, seq_num=8)
                for qual_fp in qual_fps:
                    with self.assertRaisesRegexp(ValueError, '8th sequence'):
                        reader_fn(fasta_fp, seq_num=8, qual=qual_fp)

                # seq_num too small
                with self.assertRaisesRegexp(ValueError, '`seq_num`=0'):
                    reader_fn(fasta_fp, seq_num=0)
                for qual_fp in qual_fps:
                    with self.assertRaisesRegexp(ValueError, '`seq_num`=0'):
                        reader_fn(fasta_fp, seq_num=0, qual=qual_fp)

    def test_fasta_to_sequence_collection_and_alignment(self):
        test_cases = (self.empty, self.single,
                      self.sequence_collection_different_type)

        for constructor, reader_fn in ((SequenceCollection,
                                        _fasta_to_sequence_collection),
                                       (Alignment,
                                        _fasta_to_alignment)):
            # see comment in test_fasta_to_generator_valid_files (above) for
            # testing strategy
            for exp_list, kwargs, fasta_fps, qual_fps in test_cases:
                exp = constructor(exp_list)

                for fasta_fp in fasta_fps:
                    obs = reader_fn(fasta_fp, **kwargs)

                    # TODO remove this custom equality testing code when
                    # SequenceCollection has an equals method (part of #656).
                    # We need this method to include IDs and description in the
                    # comparison (not part of SequenceCollection.__eq__).
                    self.assertEqual(len(obs), len(exp))
                    for o, e in zip(obs, exp):
                        self.assertTrue(o.equals(e, ignore=['quality']))

                    for qual_fp in qual_fps:
                        obs = reader_fn(fasta_fp, qual=qual_fp, **kwargs)

                        # TODO remove this custom equality testing code when
                        # SequenceCollection has an equals method (part of
                        # #656). We need this method to include IDs and
                        # description in the comparison (not part of
                        # SequenceCollection.__eq__).
                        self.assertEqual(obs, exp)
                        for o, e in zip(obs, exp):
                            self.assertTrue(o.equals(e))


class WriterTests(TestCase):
    def setUp(self):
        self.bio_seq1 = Sequence(
            'ACGT-acgt.', id='seq1', description='desc1',
            quality=[10, 20, 30, 10, 0, 0, 0, 88888, 1, 3456])
        self.bio_seq2 = Sequence(
            'A', id=' \n  \nseq \t2 ', quality=[42])
        self.bio_seq3 = Sequence(
            'AACGGuA', description='desc3', quality=[0, 0, 0, 0, 0, 0, 0])
        self.dna_seq = DNA(
            'ACGTTGCAccGG',
            quality=[55, 10, 0, 999, 1, 1, 8, 775, 40, 10, 10, 0],
            validate=False)
        self.rna_seq = RNA('ACGUU', quality=[10, 9, 8, 7, 6])
        self.prot_seq = Protein(
            'pQqqqPPQQQ', id='proteinseq',
            description='\ndetailed\ndescription \t\twith  new\n\nlines\n\n\n',
            quality=[42, 42, 442, 442, 42, 42, 42, 42, 42, 43], validate=False)

        seqs = [
            RNA('UUUU', id='s\te\tq\t1', description='desc\n1',
                quality=[1234, 0, 0, 2]),
            Sequence(
                'CATC', id='s\te\tq\t2', description='desc\n2',
                quality=[1, 11, 111, 11112]),
            Protein('sits', id='s\te\tq\t3', description='desc\n3',
                    quality=[12345, 678909, 999999, 4242424242],
                    validate=False)
        ]
        self.seq_coll = SequenceCollection(seqs)
        self.align = Alignment(seqs)

        def empty_gen():
            raise StopIteration()
            yield

        def single_seq_gen():
            yield self.bio_seq1

        # generate sequences with descriptions containing newlines (to test
        # description_newline_replacement)
        def newline_description_gen():
            yield self.prot_seq
            yield DNA('AGGAGAATA', id='foo', description='\n\n\n\n',
                      quality=range(9))

        # generate sequences with ids containing whitespace (to test
        # id_whitespace_replacement)
        def whitespace_id_gen():
            yield self.bio_seq2
            yield RNA('UA', id='\n\t \t', description='a\nb',
                      quality=[1000, 1])

        # multiple sequences of mixed types, lengths, and metadata. lengths are
        # chosen to exercise various splitting cases when testing max_width,
        # including exercising the different splitting algorithms used for
        # sequence data vs. quality scores
        def multi_seq_gen():
            for seq in (self.bio_seq1, self.bio_seq2, self.bio_seq3,
                        self.dna_seq, self.rna_seq, self.prot_seq):
                yield seq

        # can be serialized if no qual file is provided, else it should raise
        # an error because one seq has qual scores and the other doesn't
        def mixed_qual_score_gen():
            missing_qual_seq = Sequence(
                'AAAAT', id='da,dadadada', description='10 hours')
            for seq in self.bio_seq1, missing_qual_seq:
                yield seq

        self.mixed_qual_score_gen = mixed_qual_score_gen()

        # store sequence generator to serialize, writer kwargs (if any), and
        # fasta and qual filepaths of expected results
        self.objs_fps = list(map(lambda e: (e[0], e[1], get_data_path(e[2]),
                                            get_data_path(e[3])), [
            (empty_gen(), {}, 'empty', 'empty'),
            (single_seq_gen(), {}, 'fasta_single_seq', 'qual_single_seq'),

            # no splitting of sequence or qual data across lines b/c max_width
            # is sufficiently large
            (single_seq_gen(), {'max_width': 32}, 'fasta_single_seq',
             'qual_single_seq'),

            # splitting algorithm for sequence and qual scores is different;
            # make sure individual qual scores aren't split across lines even
            # if they exceed max_width
            (single_seq_gen(), {'max_width': 1}, 'fasta_max_width_1',
             'qual_max_width_1'),

            (multi_seq_gen(), {}, 'fasta_multi_seq', 'qual_multi_seq'),
            (multi_seq_gen(), {'max_width': 5}, 'fasta_max_width_5',
             'qual_max_width_5'),
            (newline_description_gen(),
             {'description_newline_replacement': ':-)'},
             'fasta_description_newline_replacement_multi_char',
             'qual_description_newline_replacement_multi_char'),
            (newline_description_gen(),
             {'description_newline_replacement': ''},
             'fasta_description_newline_replacement_empty_str',
             'qual_description_newline_replacement_empty_str',),
            (newline_description_gen(),
             {'description_newline_replacement': None},
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
            for seq in self.bio_seq1, Sequence(''):
                yield seq

        # generators or parameter combos that cannot be written in fasta
        # format, paired with kwargs (if any), error type, and expected error
        # message regexp
        self.invalid_objs = [
            (blank_seq_gen(), {}, ValueError, '2nd.*empty'),
            (single_seq_gen(),
             {'max_width': 0}, ValueError, 'max_width=0'),
            (multi_seq_gen(), {'id_whitespace_replacement': '-\n_'},
             ValueError, 'Newline character'),
            (multi_seq_gen(), {'description_newline_replacement': '-.-\n'},
             ValueError, 'Newline character'),
            (mixed_qual_score_gen(), {'qual': StringIO()}, ValueError,
             '2nd sequence.*does not have quality scores')
        ]

    # extensive tests for generator -> fasta writer since it is used by all
    # other object -> fasta writers

    def test_generator_to_fasta_no_qual(self):
        # test writing standalone fasta (i.e., without a qual file)
        for obj, kwargs, fp, _ in self.objs_fps:
            fh = StringIO()
            _generator_to_fasta(obj, fh, **kwargs)
            obs = fh.getvalue()
            fh.close()

            with open(fp, 'U') as fh:
                exp = fh.read()
            self.assertEqual(obs, exp)

    def test_generator_to_fasta_mixed_qual_scores(self):
        # test writing some sequences with qual scores and some without is
        # possible if no qual output file is specified
        fh = StringIO()
        _generator_to_fasta(self.mixed_qual_score_gen, fh)
        obs = fh.getvalue()
        fh.close()

        with open(get_data_path('fasta_mixed_qual_scores'), 'U') as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_generator_to_fasta_with_qual(self):
        # test writing fasta and qual files
        for obj, kwargs, fasta_fp, qual_fp in self.objs_fps:
            if qual_fp is not None:
                fasta_fh = StringIO()
                qual_fh = StringIO()
                _generator_to_fasta(obj, fasta_fh, qual=qual_fh, **kwargs)
                obs_fasta = fasta_fh.getvalue()
                obs_qual = qual_fh.getvalue()
                fasta_fh.close()
                qual_fh.close()

                with open(fasta_fp, 'U') as fh:
                    exp_fasta = fh.read()
                with open(qual_fp, 'U') as fh:
                    exp_qual = fh.read()

                self.assertEqual(obs_fasta, exp_fasta)
                self.assertEqual(obs_qual, exp_qual)

    def test_generator_to_fasta_invalid_input(self):
        for obj, kwargs, error_type, error_msg_regexp in self.invalid_objs:
            fh = StringIO()
            with self.assertRaisesRegexp(error_type, error_msg_regexp):
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
            (_biological_sequence_to_fasta,
             Sequence('ACGT', id=id_, description=desc,
                      quality=range(1, 5)),
             ('fasta_single_bio_seq_defaults',
              'fasta_single_bio_seq_non_defaults',
              'qual_single_bio_seq_non_defaults')),
            (_dna_sequence_to_fasta,
             DNA('TACG', id=id_, description=desc, quality=range(4)),
             ('fasta_single_dna_seq_defaults',
              'fasta_single_dna_seq_non_defaults',
              'qual_single_dna_seq_non_defaults')),
            (_rna_sequence_to_fasta,
             RNA('UACG', id=id_, description=desc, quality=range(2, 6)),
             ('fasta_single_rna_seq_defaults',
              'fasta_single_rna_seq_non_defaults',
              'qual_single_rna_seq_non_defaults')),
            (_protein_sequence_to_fasta,
             Protein('PQQ', id=id_, description=desc, quality=[42, 41, 40]),
             ('fasta_single_prot_seq_defaults',
              'fasta_single_prot_seq_non_defaults',
              'qual_single_prot_seq_non_defaults')))

        for fn, obj, fps in test_data:
            defaults_fp, non_defaults_fasta_fp, non_defaults_qual_fp = fps

            # test writing with default parameters
            fh = StringIO()
            fn(obj, fh)
            obs = fh.getvalue()
            fh.close()

            with open(get_data_path(defaults_fp), 'U') as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

            # test writing with non-defaults
            fasta_fh = StringIO()
            qual_fh = StringIO()
            fn(obj, fasta_fh, id_whitespace_replacement='-',
               description_newline_replacement='_', max_width=1, qual=qual_fh)
            obs_fasta = fasta_fh.getvalue()
            obs_qual = qual_fh.getvalue()
            fasta_fh.close()
            qual_fh.close()

            with open(get_data_path(non_defaults_fasta_fp), 'U') as fh:
                exp_fasta = fh.read()
            with open(get_data_path(non_defaults_qual_fp), 'U') as fh:
                exp_qual = fh.read()

            self.assertEqual(obs_fasta, exp_fasta)
            self.assertEqual(obs_qual, exp_qual)

    def test_any_sequences_to_fasta(self):
        for fn, obj in ((_sequence_collection_to_fasta, self.seq_coll),
                        (_alignment_to_fasta, self.align)):
            # test writing with default parameters
            fh = StringIO()
            fn(obj, fh)
            obs = fh.getvalue()
            fh.close()

            with open(get_data_path('fasta_3_seqs_defaults'), 'U') as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

            # test writing with non-defaults
            fasta_fh = StringIO()
            qual_fh = StringIO()
            fn(obj, fasta_fh, id_whitespace_replacement='*',
               description_newline_replacement='+', max_width=3, qual=qual_fh)
            obs_fasta = fasta_fh.getvalue()
            obs_qual = qual_fh.getvalue()
            fasta_fh.close()
            qual_fh.close()

            with open(get_data_path('fasta_3_seqs_non_defaults'), 'U') as fh:
                exp_fasta = fh.read()
            with open(get_data_path('qual_3_seqs_non_defaults'), 'U') as fh:
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
            with open(fasta_fp, 'U') as fh:
                exp_fasta = fh.read()
            with open(qual_fp, 'U') as fh:
                exp_qual = fh.read()

            fasta_fh = StringIO()
            qual_fh = StringIO()
            _generator_to_fasta(_fasta_to_generator(fasta_fp, qual=qual_fp),
                                fasta_fh, qual=qual_fh)
            obs_fasta = fasta_fh.getvalue()
            obs_qual = qual_fh.getvalue()
            fasta_fh.close()
            qual_fh.close()

            self.assertEqual(obs_fasta, exp_fasta)
            self.assertEqual(obs_qual, exp_qual)

    def test_roundtrip_sequence_collections_and_alignments(self):
        fps = list(map(lambda e: list(map(get_data_path, e)),
                       [('empty', 'empty'),
                        ('fasta_sequence_collection_different_type',
                         'qual_sequence_collection_different_type')]))

        for reader, writer in ((_fasta_to_sequence_collection,
                                _sequence_collection_to_fasta),
                               (_fasta_to_alignment,
                                _alignment_to_fasta)):
            for fasta_fp, qual_fp in fps:
                # read
                obj1 = reader(fasta_fp, qual=qual_fp)

                # write
                fasta_fh = StringIO()
                qual_fh = StringIO()
                writer(obj1, fasta_fh, qual=qual_fh)
                fasta_fh.seek(0)
                qual_fh.seek(0)

                # read
                obj2 = reader(fasta_fh, qual=qual_fh)
                fasta_fh.close()
                qual_fh.close()

                # TODO remove this custom equality testing code when
                # SequenceCollection has an equals method (part of #656).
                # We need this method to include IDs and description in the
                # comparison (not part of SequenceCollection.__eq__).
                self.assertEqual(obj1, obj2)
                for s1, s2 in zip(obj1, obj2):
                    self.assertTrue(s1.equals(s2))

    def test_roundtrip_biological_sequences(self):
        fps = list(map(lambda e: list(map(get_data_path, e)),
                       [('fasta_multi_seq_roundtrip',
                         'qual_multi_seq_roundtrip'),
                        ('fasta_sequence_collection_different_type',
                         'qual_sequence_collection_different_type')]))

        for reader, writer in ((_fasta_to_biological_sequence,
                                _biological_sequence_to_fasta),
                               (_fasta_to_dna_sequence,
                                _dna_sequence_to_fasta),
                               (_fasta_to_rna_sequence,
                                _rna_sequence_to_fasta),
                               (_fasta_to_protein_sequence,
                                _protein_sequence_to_fasta)):
            for fasta_fp, qual_fp in fps:
                # read
                obj1 = reader(fasta_fp, qual=qual_fp)

                # write
                fasta_fh = StringIO()
                qual_fh = StringIO()
                writer(obj1, fasta_fh, qual=qual_fh)
                fasta_fh.seek(0)
                qual_fh.seek(0)

                # read
                obj2 = reader(fasta_fh, qual=qual_fh)
                fasta_fh.close()
                qual_fh.close()

                self.assertTrue(obj1.equals(obj2))


if __name__ == '__main__':
    main()
