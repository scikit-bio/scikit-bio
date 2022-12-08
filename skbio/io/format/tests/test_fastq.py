# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import string
import unittest
import warnings
from functools import partial

from skbio import read, write, Sequence, DNA, RNA, Protein, TabularMSA
from skbio.io import FASTQFormatError
from skbio.io.format.fastq import (
    _fastq_sniffer, _fastq_to_generator, _fastq_to_tabular_msa,
    _generator_to_fastq, _tabular_msa_to_fastq)
from skbio.sequence import GrammaredSequence
from skbio.util import get_data_path
from skbio.util import classproperty
from skbio.util._decorator import overrides

import numpy as np

# Note: the example FASTQ files with file extension .fastq are taken from the
# following open-access publication's supplementary data:
#
# P.J.A. Cock, C.J. Fields, N. Goto, M.L. Heuer and P.M. Rice (2009). The
# Sanger FASTQ file format for sequences with quality scores, and the
# Solexa/Illumina FASTQ variants.
#
# See licenses/fastq-example-files-readme.txt for the original README that
# accompanied these files, which includes the terms of use and detailed
# description of the files.
#
# The example files bearing the original filenames have not been modified from
# their original form.


def _drop_kwargs(kwargs, *args):
    for arg in args:
        if arg in kwargs:
            kwargs.pop(arg)


class TestSniffer(unittest.TestCase):
    def setUp(self):
        self.positives = [get_data_path(e) for e in [
            'fastq_multi_seq_sanger',
            'fastq_multi_blank_between_records',
            'fastq_multi_ws_lines_between_records',
            'fastq_multi_blank_end_of_file',
            'fastq_multi_ws_lines_end_of_file',
            'fastq_multi_whitespace_stripping',
            'fastq_blank_lines',
            'fastq_whitespace_only_lines',
            'fastq_single_seq_illumina1.3',
            'fastq_wrapping_as_illumina_no_description',
            'fastq_wrapping_as_sanger_no_description',
            'fastq_wrapping_original_sanger_no_description',
            'fastq_writer_illumina1.3_defaults',
            'fastq_writer_sanger_defaults',
            'fastq_writer_sanger_non_defaults',
            'fastq_5_blanks_start_of_file',
            'fastq_5_ws_lines_start_of_file',
            'illumina_full_range_as_illumina.fastq',
            'illumina_full_range_as_sanger.fastq',
            'illumina_full_range_original_illumina.fastq',
            'longreads_as_illumina.fastq',
            'longreads_as_sanger.fastq',
            'longreads_original_sanger.fastq',
            'misc_dna_as_illumina.fastq',
            'misc_dna_as_sanger.fastq',
            'misc_dna_original_sanger.fastq',
            'misc_rna_as_illumina.fastq',
            'misc_rna_as_sanger.fastq',
            'misc_rna_original_sanger.fastq',
            'sanger_full_range_as_illumina.fastq',
            'sanger_full_range_as_sanger.fastq',
            'sanger_full_range_original_sanger.fastq',
            'solexa_full_range_original_solexa.fastq',
            'wrapping_as_illumina.fastq',
            'wrapping_as_sanger.fastq',
            'wrapping_original_sanger.fastq'
        ]]

        self.negatives = [get_data_path(e) for e in [
            'empty',
            'whitespace_only',
            'fastq_multi_blank_start_of_file',
            'fastq_multi_ws_lines_start_of_file',
            'fastq_invalid_blank_after_header',
            'fastq_invalid_blank_after_seq',
            'fastq_invalid_blank_after_plus',
            'fastq_invalid_blank_within_seq',
            'fastq_invalid_blank_within_qual',
            'fastq_invalid_ws_line_after_header',
            'fastq_invalid_ws_line_after_seq',
            'fastq_invalid_ws_line_after_plus',
            'fastq_invalid_ws_line_within_seq',
            'fastq_invalid_ws_line_within_qual',
            'fastq_invalid_missing_header',
            'fastq_invalid_missing_seq_data',
            'error_diff_ids.fastq',
            'error_double_qual.fastq',
            'error_double_seq.fastq',
            'error_long_qual.fastq',
            'error_no_qual.fastq',
            'error_qual_del.fastq',
            'error_qual_escape.fastq',
            'error_qual_null.fastq',
            'error_qual_space.fastq',
            'error_qual_tab.fastq',
            'error_qual_unit_sep.fastq',
            'error_qual_vtab.fastq',
            'error_short_qual.fastq',
            'error_spaces.fastq',
            'error_tabs.fastq',
            'error_trunc_at_seq.fastq',
            'error_trunc_at_plus.fastq',
            'error_trunc_at_qual.fastq',
            'error_trunc_in_title.fastq',
            'error_trunc_in_seq.fastq',
            'error_trunc_in_plus.fastq',
            'error_trunc_in_qual.fastq',
        ]]

    def test_positives(self):
        for fp in self.positives:
            self.assertEqual(_fastq_sniffer(fp), (True, {}))

    def test_negatives(self):
        for fp in self.negatives:
            self.assertEqual(_fastq_sniffer(fp), (False, {}))

    def test_illumina_sniffed(self):
        fp = get_data_path('fastq_single_seq_illumina1.8')
        self.assertEqual(_fastq_sniffer(fp), (True, {'variant':
                                                     'illumina1.8'}))


class TestReaders(unittest.TestCase):
    def setUp(self):
        self.valid_configurations = [
            ([get_data_path('empty'),
              get_data_path('whitespace_only')],
             [{},
              {'variant': 'illumina1.8'},
              {'phred_offset': 33,
               'constructor': DNA}],
             []),

            ([get_data_path('fastq_single_seq_illumina1.3')], [
                {'variant': 'illumina1.3'},
                {'phred_offset': 64},
                {'variant': 'illumina1.3',
                 'constructor': Protein},
            ], [
                ('', 'bar\t baz', 'aCGT', [33, 34, 35, 36])
            ]),

            ([get_data_path('fastq_multi_seq_sanger'),
              get_data_path('fastq_whitespace_only_lines'),
              get_data_path('fastq_blank_lines'),
              get_data_path('fastq_multi_blank_between_records'),
              get_data_path('fastq_multi_ws_lines_between_records'),
              get_data_path('fastq_multi_blank_end_of_file'),
              get_data_path('fastq_multi_ws_lines_end_of_file'),
              get_data_path('fastq_multi_blank_start_of_file'),
              get_data_path('fastq_multi_ws_lines_start_of_file'),
              get_data_path('fastq_multi_whitespace_stripping')], [
                {'variant': 'sanger'},
                {'phred_offset': 33, 'seq_num': 2},
                {'variant': 'sanger',
                 'constructor': partial(RNA, validate=False),
                 'seq_num': 3},
            ], [
                ('foo', 'bar baz', 'AACCGG',
                 [16, 17, 18, 19, 20, 21]),
                ('bar', 'baz foo', 'TTGGCC',
                 [23, 22, 21, 20, 19, 18]),
                ('baz', 'foo bar', 'GATTTC',
                 [20, 21, 22, 23, 24, 18])
            ]),


        ]

        self.invalid_files = [(get_data_path(e[0]), e[1], e[2]) for e in [
            ('fastq_invalid_blank_after_header', FASTQFormatError,
             r'blank or whitespace-only line.*after header.*in FASTQ'),

            ('fastq_invalid_blank_after_seq', FASTQFormatError,
             r"blank or whitespace-only line.*before '\+' in FASTQ"),

            ('fastq_invalid_blank_after_plus', FASTQFormatError,
             r"blank or whitespace-only line.*after '\+'.*in FASTQ"),

            ('fastq_invalid_blank_within_seq', FASTQFormatError,
             r'blank or whitespace-only line.*within sequence.*FASTQ'),

            ('fastq_invalid_blank_within_qual', FASTQFormatError,
             r"blank or whitespace-only line.*within quality scores.*in "
             "FASTQ"),

            ('fastq_invalid_ws_line_after_header', FASTQFormatError,
             r'blank or whitespace-only line.*after header.*in FASTQ'),

            ('fastq_invalid_ws_line_after_seq', FASTQFormatError,
             r"blank or whitespace-only line.*before '\+' in FASTQ"),

            ('fastq_invalid_ws_line_after_plus', FASTQFormatError,
             r"blank or whitespace-only line.*after '\+'.*in FASTQ"),

            ('fastq_invalid_ws_line_within_seq', FASTQFormatError,
             r'blank or whitespace-only line.*within sequence.*FASTQ'),

            ('fastq_invalid_ws_line_within_qual', FASTQFormatError,
             r"blank or whitespace-only line.*within quality scores.*in "
             "FASTQ"),

            ('fastq_invalid_missing_header', FASTQFormatError,
             r"sequence.*header.*start of file: 'seq1 desc1'"),

            ('fastq_invalid_missing_seq_data', FASTQFormatError,
             r'without sequence data'),

            ('error_diff_ids.fastq', FASTQFormatError,
             r"header lines do not match: "
             "'SLXA-B3_649_FC8437_R1_1_1_850_123' != "
             "'SLXA-B3_649_FC8437_R1_1_1_850_124'"),

            ('error_double_qual.fastq', FASTQFormatError,
             r"Extra quality.*'\+SLXA-B3_649_FC8437_R1_1_1_850_123'"),

            ('error_double_seq.fastq', FASTQFormatError,
             r'FASTQ record that is missing a quality \(\+\) header line'),

            ('error_long_qual.fastq', FASTQFormatError, r"Extra quality.*'Y'"),

            ('error_no_qual.fastq', FASTQFormatError,
             r"blank or whitespace-only line.*after '\+'.*in FASTQ"),

            ('error_qual_del.fastq', ValueError,
             r'Decoded Phred score.*out of range'),

            ('error_qual_escape.fastq', ValueError,
             r'Decoded Phred score.*out of range'),

            ('error_qual_null.fastq', ValueError,
             r'Decoded Phred score.*out of range'),

            ('error_qual_space.fastq', ValueError,
             r'Decoded Phred score.*out of range'),

            ('error_qual_tab.fastq', ValueError,
             r'Decoded Phred score.*out of range'),

            ('error_qual_unit_sep.fastq', ValueError,
             r'Decoded Phred score.*out of range'),

            ('error_qual_vtab.fastq', ValueError,
             r'Decoded Phred score.*out of range'),

            ('error_short_qual.fastq', FASTQFormatError,
             r"Extra quality.*'SLXA-B3_649_FC8437_R1_1_1_362_549'"),

            ('error_spaces.fastq', FASTQFormatError,
             r"whitespace.*sequence data: 'GATGTGCAA TACCTTTGTA GAGGAA'"),

            ('error_tabs.fastq', FASTQFormatError,
             r"whitespace.*sequence data: 'GATGTGCAA\\tTACCTTTGTA\\tGAGGAA'"),

            ('error_trunc_at_seq.fastq', FASTQFormatError,
             r'incomplete/truncated.*FASTQ'),

            ('error_trunc_at_plus.fastq', FASTQFormatError,
             r'incomplete/truncated.*FASTQ'),

            ('error_trunc_at_qual.fastq', FASTQFormatError,
             r'incomplete/truncated.*end of file'),

            ('error_trunc_in_title.fastq', FASTQFormatError,
             r'incomplete/truncated.*end of file'),

            ('error_trunc_in_seq.fastq', FASTQFormatError,
             r'incomplete/truncated.*end of file'),

            ('error_trunc_in_plus.fastq', FASTQFormatError,
             r"header lines do not match: "
             "'SLXA-B3_649_FC8437_R1_1_1_183_714' != 'SLXA-B3_649_FC'"),

            ('error_trunc_in_qual.fastq', FASTQFormatError,
             r'incomplete/truncated.*end of file')
        ]]

    def test_fastq_to_generator_valid_files(self):
        for valid_files, kwargs, components in self.valid_configurations:
            for valid in valid_files:
                for observed_kwargs in kwargs:
                    _drop_kwargs(observed_kwargs, 'seq_num')
                    constructor = observed_kwargs.get('constructor', Sequence)

                    expected_kwargs = {}
                    expected_kwargs['lowercase'] = 'introns'
                    observed_kwargs['lowercase'] = 'introns'

                    expected = [constructor(c[2],
                                            metadata={'id': c[0],
                                                      'description': c[1]},
                                positional_metadata={'quality': np.array(c[3],
                                                     dtype=np.uint8)},
                                **expected_kwargs)
                                for c in components]

                    observed = list(_fastq_to_generator(valid,
                                                        **observed_kwargs))
                    self.assertEqual(len(expected), len(observed))
                    for o, e in zip(observed, expected):
                        self.assertEqual(o, e)

    def test_fastq_to_generator_invalid_files_all_variants(self):
        # files that should be invalid for all variants, as well as custom
        # phred offsets
        for fp, error_type, error_msg_regex in self.invalid_files:
            for variant in 'sanger', 'illumina1.3', 'illumina1.8':
                with self.assertRaisesRegex(error_type, error_msg_regex):
                    list(_fastq_to_generator(fp, variant=variant))

            for offset in 33, 64, 40, 77:
                with self.assertRaisesRegex(error_type, error_msg_regex):
                    list(_fastq_to_generator(fp, phred_offset=offset))

    def test_fastq_to_generator_invalid_files_illumina(self):
        # files that should be invalid for illumina1.3 and illumina1.8 variants
        fps = [get_data_path(fp) for fp in
               ['sanger_full_range_original_sanger.fastq',
               'solexa_full_range_original_solexa.fastq']]

        for fp in fps:
            with self.assertRaisesRegex(ValueError, r'out of range \[0, 62\]'):
                list(_fastq_to_generator(fp, variant='illumina1.3'))
            with self.assertRaisesRegex(ValueError, r'out of range \[0, 62\]'):
                list(_fastq_to_generator(fp, variant='illumina1.8'))

    def test_fastq_to_generator_solexa(self):
        # solexa support isn't implemented yet. should raise error even with
        # valid solexa file
        with self.assertRaisesRegex(ValueError, r'Solexa'):
            list(_fastq_to_generator(
                get_data_path('solexa_full_range_original_solexa.fastq'),
                variant='solexa'))

    def test_fastq_to_sequence(self):
        for constructor in [Sequence, DNA, RNA, Protein]:
            for valid_files, kwargs, components in self.valid_configurations:
                for valid in valid_files:
                    # skip empty file case since we cannot read a specific
                    # sequencefrom an empty file
                    if len(components) == 0:
                        continue

                    for observed_kwargs in kwargs:
                        expected_kwargs = {}

                        # TODO:
                        # some of the test files contain characters which are
                        # invalid for RNA, so don't validate for now. Need to
                        # fix this
                        if constructor is RNA:
                            observed_kwargs['validate'] = False
                            expected_kwargs['validate'] = False

                        _drop_kwargs(observed_kwargs, 'constructor')

                        expected_kwargs['lowercase'] = 'introns'
                        observed_kwargs['lowercase'] = 'introns'

                        seq_num = observed_kwargs.get('seq_num', 1)
                        c = components[seq_num - 1]
                        expected = \
                            constructor(
                                c[2], metadata={'id': c[0],
                                                'description': c[1]},
                                positional_metadata={'quality': np.array(c[3],
                                                     dtype=np.uint8)},
                                **expected_kwargs)

                        observed = read(valid, into=constructor,
                                        format='fastq', verify=False,
                                        **observed_kwargs)
                        self.assertEqual(observed, expected)

    def test_fastq_to_tabular_msa(self):
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

        for valid_files, kwargs, components in self.valid_configurations:
            for valid in valid_files:
                for observed_kwargs in kwargs:
                    _drop_kwargs(observed_kwargs, 'seq_num')
                    if 'constructor' not in observed_kwargs:
                        observed_kwargs['constructor'] = CustomSequence
                    constructor = observed_kwargs['constructor']

                    expected_kwargs = {}
                    expected_kwargs['lowercase'] = 'introns'
                    observed_kwargs['lowercase'] = 'introns'

                    expected = TabularMSA(
                        [constructor(
                            c[2], metadata={'id': c[0],
                                            'description': c[1]},
                            positional_metadata={'quality': np.array(c[3],
                                                 dtype=np.uint8)},
                            **expected_kwargs)
                         for c in components])

                    observed = _fastq_to_tabular_msa(valid, **observed_kwargs)
                    self.assertEqual(observed, expected)

    def test_fastq_to_tabular_msa_no_constructor(self):
        with self.assertRaisesRegex(ValueError, r'`constructor`'):
            _fastq_to_tabular_msa(get_data_path('fastq_multi_seq_sanger'))


class TestWriters(unittest.TestCase):
    def setUp(self):
        self.valid_files = [
            ([
                ('f o  o', 'bar\n\nbaz', 'AaCcGg',
                 [16, 17, 18, 19, 20, 21]),
                ('bar', 'baz foo', 'TtGgCc',
                 [23, 22, 21, 20, 19, 18]),
                ('ba\n\t\tz', 'foo bar', 'gAtTtC',
                 [20, 21, 22, 23, 24, 18])
            ], [
                ({'variant': 'sanger'},
                 get_data_path('fastq_writer_sanger_defaults')),
                ({'phred_offset': 33},
                 get_data_path('fastq_writer_sanger_defaults')),
                ({'variant': 'illumina1.8'},
                 get_data_path('fastq_writer_sanger_defaults')),
                ({'variant': 'illumina1.3'},
                 get_data_path('fastq_writer_illumina1.3_defaults')),
                ({'variant': 'sanger', 'id_whitespace_replacement': '%',
                  'description_newline_replacement': '^'},
                 get_data_path('fastq_writer_sanger_non_defaults'))
            ]),
        ]

    def test_generator_to_fastq_kwargs_passed(self):
        for components, kwargs_expected_fp in self.valid_files:
            for kwargs, expected_fp in kwargs_expected_fp:
                def gen():
                    for c in components:
                        yield Sequence(
                            c[2], metadata={'id': c[0], 'description': c[1]},
                            positional_metadata={'quality': c[3]})

                fh = io.StringIO()
                _generator_to_fastq(gen(), fh, **kwargs)
                observed = fh.getvalue()
                fh.close()

                with io.open(expected_fp) as f:
                    expected = f.read()

                self.assertEqual(observed, expected)

    def test_sequence_to_fastq_kwargs_passed(self):
        for constructor in [Sequence, DNA, RNA, Protein]:
            for components, kwargs_expected_fp in self.valid_files:
                for expected_kwargs, expected_fp in kwargs_expected_fp:

                    observed_kwargs = {}
                    # TODO:
                    # some of the test files contain characters which are
                    # invalid for RNA, so don't validate for now. Need to
                    # fix this
                    if constructor is RNA:
                        observed_kwargs['validate'] = False

                    expected_kwargs['lowercase'] = 'introns'
                    observed_kwargs['lowercase'] = 'introns'

                    fh = io.StringIO()
                    for c in components:
                        obj = constructor(
                            c[2],
                            metadata={'id': c[0], 'description': c[1]},
                            positional_metadata={'quality': c[3]},
                            **observed_kwargs)
                        write(obj, into=fh, format='fastq', **expected_kwargs)

                    observed = fh.getvalue()
                    fh.close()

                    with io.open(expected_fp) as f:
                        expected = f.read()

                    self.assertEqual(observed, expected)

    def test_tabular_msa_to_fastq_kwargs_passed(self):
        for components, kwargs_expected_fp in self.valid_files:
            for kwargs, expected_fp in kwargs_expected_fp:
                obj = TabularMSA([
                    Protein(c[2], metadata={'id': c[0], 'description': c[1]},
                            positional_metadata={'quality': c[3]},
                            lowercase='introns')
                    for c in components])

                fh = io.StringIO()
                kwargs['lowercase'] = 'introns'
                _tabular_msa_to_fastq(obj, fh, **kwargs)
                observed = fh.getvalue()
                fh.close()

                with io.open(expected_fp) as f:
                    expected = f.read()

                self.assertEqual(observed, expected)

    def test_generator_to_fastq_no_qual(self):
        def gen():
            yield Sequence('ACGT',
                           metadata={'id': 'foo', 'description': 'bar'},
                           positional_metadata={'quality': range(4)})
            yield Sequence('ACG', metadata={'id': 'foo', 'description': 'bar'})

        with self.assertRaisesRegex(ValueError, r'2nd.*quality scores'):
            _generator_to_fastq(gen(), io.StringIO(), variant='illumina1.8')


class TestConversions(unittest.TestCase):
    def setUp(self):
        self.conversions = [
            (get_data_path('empty'),
             get_data_path('empty'), [
                 ({'variant': 'sanger'}, {'phred_offset': 42}),
            ]),

            (get_data_path('longreads_original_sanger.fastq'),
             get_data_path('longreads_as_sanger.fastq'), [
                 ({'variant': 'sanger'}, {'variant': 'sanger'}),
                 ({'phred_offset': 33}, {'variant': 'sanger'}),
                 ({'variant': 'sanger'}, {'phred_offset': 33})
            ]),
            (get_data_path('longreads_original_sanger.fastq'),
             get_data_path('longreads_as_illumina.fastq'), [
                 ({'variant': 'sanger'}, {'variant': 'illumina1.3'}),
                 ({'phred_offset': 33}, {'variant': 'illumina1.3'}),
                 ({'variant': 'sanger'}, {'phred_offset': 64})
            ]),

            (get_data_path('wrapping_original_sanger.fastq'),
             get_data_path('wrapping_as_sanger.fastq'), [
                 ({'variant': 'sanger'}, {'variant': 'sanger'}),
                 ({'phred_offset': 33}, {'variant': 'sanger'}),
                 ({'variant': 'sanger'}, {'phred_offset': 33})
            ]),
            (get_data_path('wrapping_original_sanger.fastq'),
             get_data_path('wrapping_as_illumina.fastq'), [
                 ({'variant': 'sanger'}, {'variant': 'illumina1.3'}),
                 ({'phred_offset': 33}, {'variant': 'illumina1.3'}),
                 ({'variant': 'sanger'}, {'phred_offset': 64})
            ]),

            (get_data_path('sanger_full_range_original_sanger.fastq'),
             get_data_path('sanger_full_range_as_sanger.fastq'), [
                 ({'variant': 'sanger'}, {'variant': 'sanger'}),
                 ({'phred_offset': 33}, {'variant': 'sanger'}),
                 ({'variant': 'sanger'}, {'phred_offset': 33})
            ]),
            (get_data_path('sanger_full_range_original_sanger.fastq'),
             get_data_path('sanger_full_range_as_illumina.fastq'), [
                 ({'variant': 'sanger'}, {'variant': 'illumina1.3'}),
                 ({'phred_offset': 33}, {'variant': 'illumina1.3'}),
                 ({'variant': 'sanger'}, {'phred_offset': 64})
            ]),

            (get_data_path('illumina_full_range_original_illumina.fastq'),
             get_data_path('illumina_full_range_as_illumina.fastq'), [
                 ({'variant': 'illumina1.3'}, {'variant': 'illumina1.3'}),
                 ({'phred_offset': 64}, {'variant': 'illumina1.3'}),
                 ({'variant': 'illumina1.3'}, {'phred_offset': 64})
            ]),
            (get_data_path('illumina_full_range_original_illumina.fastq'),
             get_data_path('illumina_full_range_as_sanger.fastq'), [
                 ({'variant': 'illumina1.3'}, {'variant': 'sanger'}),
                 ({'phred_offset': 64}, {'variant': 'sanger'}),
                 ({'variant': 'illumina1.3'}, {'phred_offset': 33})
            ]),

            (get_data_path('misc_dna_original_sanger.fastq'),
             get_data_path('misc_dna_as_sanger.fastq'), [
                 ({'variant': 'sanger'}, {'variant': 'sanger'}),
                 ({'phred_offset': 33}, {'variant': 'sanger'}),
                 ({'variant': 'sanger'}, {'phred_offset': 33})
            ]),
            (get_data_path('misc_dna_original_sanger.fastq'),
             get_data_path('misc_dna_as_illumina.fastq'), [
                 ({'variant': 'sanger'}, {'variant': 'illumina1.3'}),
                 ({'phred_offset': 33}, {'variant': 'illumina1.3'}),
                 ({'variant': 'sanger'}, {'phred_offset': 64})
            ]),

            (get_data_path('misc_rna_original_sanger.fastq'),
             get_data_path('misc_rna_as_sanger.fastq'), [
                 ({'variant': 'sanger'}, {'variant': 'sanger'}),
                 ({'phred_offset': 33}, {'variant': 'sanger'}),
                 ({'variant': 'sanger'}, {'phred_offset': 33})
            ]),
            (get_data_path('misc_rna_original_sanger.fastq'),
             get_data_path('misc_rna_as_illumina.fastq'), [
                 ({'variant': 'sanger'}, {'variant': 'illumina1.3'}),
                 ({'phred_offset': 33}, {'variant': 'illumina1.3'}),
                 ({'variant': 'sanger'}, {'phred_offset': 64})
            ]),

            (get_data_path('fastq_wrapping_original_sanger_no_description'),
             get_data_path('fastq_wrapping_as_sanger_no_description'), [
                 ({'variant': 'sanger'}, {'variant': 'sanger'}),
                 ({'phred_offset': 33}, {'variant': 'sanger'}),
                 ({'variant': 'sanger'}, {'phred_offset': 33})
            ]),
            (get_data_path('fastq_wrapping_original_sanger_no_description'),
             get_data_path('fastq_wrapping_as_illumina_no_description'), [
                 ({'variant': 'sanger'}, {'variant': 'illumina1.3'}),
                 ({'phred_offset': 33}, {'variant': 'illumina1.3'}),
                 ({'variant': 'sanger'}, {'phred_offset': 64})
            ]),
        ]

    def test_conversion(self):
        for from_fp, to_fp, kwargs in self.conversions:
            for from_kwargs, to_kwargs in kwargs:
                read_gen = _fastq_to_generator(from_fp, **from_kwargs)
                fh = io.StringIO()

                # will issue warning when truncating quality scores
                with warnings.catch_warnings(record=True):
                    warnings.simplefilter("ignore")
                    _generator_to_fastq(read_gen, fh, **to_kwargs)

                obs = fh.getvalue()
                fh.close()

                with io.open(to_fp) as fh:
                    exp = fh.read()
                self.assertEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()
