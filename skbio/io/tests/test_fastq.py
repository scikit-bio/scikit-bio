# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import zip
from six import StringIO

import unittest
import warnings

from skbio import (read, write, BiologicalSequence, NucleotideSequence,
                   DNASequence, RNASequence, ProteinSequence,
                   SequenceCollection, Alignment)
from skbio.io import FASTQFormatError
from skbio.io.fastq import (
    _fastq_sniffer, _fastq_to_generator, _fastq_to_sequence_collection,
    _fastq_to_alignment, _generator_to_fastq, _biological_sequence_to_fastq,
    _nucleotide_sequence_to_fastq, _dna_sequence_to_fastq,
    _rna_sequence_to_fastq, _protein_sequence_to_fastq,
    _sequence_collection_to_fastq, _alignment_to_fastq)

from skbio.util import get_data_path

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
# The example files have not been modified from their original form.


def _drop_kwargs(kwargs, *args):
    for arg in args:
        if arg in kwargs:
            kwargs.pop(arg)


class TestReaders(unittest.TestCase):
    def setUp(self):
        self.valid_files = [
            (get_data_path('empty'),
             [{},
              {'variant': 'illumina1.8'},
              {'phred_offset': 33, 'constructor': DNASequence}],
             []),

            (get_data_path('fastq_single_seq_illumina1.3'), [
                {'variant': 'illumina1.3'},
                {'phred_offset': 64},
                {'variant': 'illumina1.3', 'constructor': ProteinSequence},
            ], [
                ('', 'bar\t baz', 'ACGT', [33, 34, 35, 36])
            ]),

            (get_data_path('fastq_multi_seq_sanger'), [
                {'variant': 'sanger'},
                {'phred_offset': 33, 'seq_num': 2},
                {'variant': 'sanger', 'constructor': RNASequence,
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
            ('whitespace_only', FASTQFormatError, 'blank line.*FASTQ'),

            ('fastq_invalid_missing_header', FASTQFormatError,
             "sequence.*header.*start of file: 'seq1 desc1'"),

            ('fastq_invalid_missing_seq_data', FASTQFormatError,
             'without sequence data'),

            ('error_diff_ids.fastq', FASTQFormatError,
             "header lines do not match: "
             "'SLXA-B3_649_FC8437_R1_1_1_850_123' != "
             "'SLXA-B3_649_FC8437_R1_1_1_850_124'"),

            ('error_double_qual.fastq', FASTQFormatError,
             "Extra quality.*'\+SLXA-B3_649_FC8437_R1_1_1_850_123'"),

            ('error_double_seq.fastq', FASTQFormatError,
             'FASTQ record that is missing a quality \(\+\) header line'),

            ('error_long_qual.fastq', FASTQFormatError, "Extra quality.*'Y'"),

            ('error_no_qual.fastq', FASTQFormatError,
             'blank line.*FASTQ'),

            ('error_qual_del.fastq', ValueError,
             'Decoded Phred score.*out of range'),

            ('error_qual_escape.fastq', ValueError,
             'Decoded Phred score.*out of range'),

            ('error_qual_null.fastq', ValueError,
             'Decoded Phred score.*out of range'),

            ('error_qual_space.fastq', ValueError,
             'Decoded Phred score.*out of range'),

            ('error_qual_tab.fastq', ValueError,
             'Decoded Phred score.*out of range'),

            ('error_qual_unit_sep.fastq', ValueError,
             'Decoded Phred score.*out of range'),

            ('error_qual_vtab.fastq', ValueError,
             'Decoded Phred score.*out of range'),

            ('error_short_qual.fastq', FASTQFormatError,
             "Extra quality.*'SLXA-B3_649_FC8437_R1_1_1_362_549'"),

            ('error_spaces.fastq', FASTQFormatError,
             "whitespace.*sequence data: 'GATGTGCAA TACCTTTGTA GAGGAA'"),

            ('error_tabs.fastq', FASTQFormatError,
             r"whitespace.*sequence data: 'GATGTGCAA\\tTACCTTTGTA\\tGAGGAA'"),

            ('error_trunc_at_seq.fastq', FASTQFormatError,
             'blank line.*FASTQ'),

            ('error_trunc_at_plus.fastq', FASTQFormatError,
             'blank line.*FASTQ'),

            ('error_trunc_at_qual.fastq', FASTQFormatError,
             'incomplete/truncated.*end of file'),

            ('error_trunc_in_title.fastq', FASTQFormatError,
             'incomplete/truncated.*end of file'),

            ('error_trunc_in_seq.fastq', FASTQFormatError,
             'incomplete/truncated.*end of file'),

            ('error_trunc_in_plus.fastq', FASTQFormatError,
             "header lines do not match: "
             "'SLXA-B3_649_FC8437_R1_1_1_183_714' != 'SLXA-B3_649_FC'"),

            ('error_trunc_in_qual.fastq', FASTQFormatError,
             'incomplete/truncated.*end of file')
        ]]

    def test_fastq_to_generator_valid_files(self):
        for valid, kwargs, components in self.valid_files:
            for kwarg in kwargs:
                _drop_kwargs(kwarg, 'seq_num')
                constructor = kwarg.get('constructor', BiologicalSequence)
                expected = [constructor(c[2], id=c[0], description=c[1],
                            quality=c[3]) for c in components]

                observed = list(_fastq_to_generator(valid, **kwarg))
                self.assertEqual(len(expected), len(observed))
                for o, e in zip(observed, expected):
                    self.assertTrue(o.equals(e))

    def test_fastq_to_generator_invalid_files_all_variants(self):
        # files that should be invalid for all variants, as well as custom
        # phred offsets
        for fp, error_type, error_msg_regex in self.invalid_files:
            for variant in 'sanger', 'illumina1.3', 'illumina1.8':
                with self.assertRaisesRegexp(error_type, error_msg_regex):
                    list(_fastq_to_generator(fp, variant=variant))

            for offset in 33, 64, 40, 77:
                with self.assertRaisesRegexp(error_type, error_msg_regex):
                    list(_fastq_to_generator(fp, phred_offset=offset))

    def test_fastq_to_generator_invalid_files_illumina(self):
        # files that should be invalid for illumina1.3 and illumina1.8 variants
        fps = [get_data_path(fp) for fp in
               'sanger_full_range_original_sanger.fastq',
               'solexa_full_range_original_solexa.fastq']

        for fp in fps:
            with self.assertRaisesRegexp(ValueError, 'out of range \[0, 62\]'):
                list(_fastq_to_generator(fp, variant='illumina1.3'))
            with self.assertRaisesRegexp(ValueError, 'out of range \[0, 62\]'):
                list(_fastq_to_generator(fp, variant='illumina1.8'))

    def test_fastq_to_generator_solexa(self):
        # solexa support isn't implemented yet. should raise error even with
        # valid solexa file
        with self.assertRaises(NotImplementedError):
            list(_fastq_to_generator(
                get_data_path('solexa_full_range_original_solexa.fastq'),
                variant='solexa'))

    def test_fastq_to_sequence(self):
        for constructor in [BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence]:
            for valid, kwargs, components in self.valid_files:
                # skip empty file case since we cannot read a specific sequence
                # from an empty file
                if len(components) == 0:
                    continue

                for kwarg in kwargs:
                    _drop_kwargs(kwarg, 'constructor')

                    seq_num = kwarg.get('seq_num', 1)
                    c = components[seq_num - 1]
                    expected = constructor(c[2], id=c[0], description=c[1],
                                           quality=c[3])

                    observed = read(valid, into=constructor, format='fastq',
                                    verify=False, **kwarg)
                    self.assertTrue(observed.equals(expected))

    def test_fastq_to_sequence_collection(self):
        for valid, kwargs, components in self.valid_files:
            for kwarg in kwargs:
                _drop_kwargs(kwarg, 'seq_num')
                constructor = kwarg.get('constructor', BiologicalSequence)
                expected = SequenceCollection(
                    [constructor(c[2], id=c[0], description=c[1], quality=c[3])
                     for c in components])

                observed = _fastq_to_sequence_collection(valid, **kwarg)
                # TODO remove when #656 is resolved
                self.assertEqual(observed, expected)
                for o, e in zip(observed, expected):
                    self.assertTrue(o.equals(e))

    def test_fastq_to_alignment(self):
        for valid, kwargs, components in self.valid_files:
            for kwarg in kwargs:
                _drop_kwargs(kwarg, 'seq_num')
                constructor = kwarg.get('constructor', BiologicalSequence)
                expected = Alignment(
                    [constructor(c[2], id=c[0], description=c[1], quality=c[3])
                     for c in components])

                observed = _fastq_to_alignment(valid, **kwarg)
                # TODO remove when #656 is resolved
                self.assertEqual(observed, expected)
                for o, e in zip(observed, expected):
                    self.assertTrue(o.equals(e))


class TestWriters(unittest.TestCase):
    def setUp(self):
        self.valid_files = [
            ([
                ('f o  o', 'bar\n\nbaz', 'AACCGG',
                 [16, 17, 18, 19, 20, 21]),
                ('bar', 'baz foo', 'TTGGCC',
                 [23, 22, 21, 20, 19, 18]),
                ('ba\n\t\tz', 'foo bar', 'GATTTC',
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
                        yield BiologicalSequence(c[2], id=c[0],
                                description=c[1], quality=c[3])

                fh = StringIO()
                _generator_to_fastq(gen(), fh, **kwargs)
                observed = fh.getvalue()
                fh.close()

                with open(expected_fp, 'U') as f:
                    expected = f.read()

                self.assertEqual(observed, expected)

    def test_sequence_to_fastq_kwargs_passed(self):
        for constructor in [BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence]:
            for components, kwargs_expected_fp in self.valid_files:
                for kwargs, expected_fp in kwargs_expected_fp:
                    fh = StringIO()
                    for c in components:
                        obj = constructor(c[2], id=c[0], description=c[1],
                                          quality=c[3])
                        write(obj, into=fh, format='fastq', **kwargs)

                    observed = fh.getvalue()
                    fh.close()

                    with open(expected_fp, 'U') as f:
                        expected = f.read()

                    self.assertEqual(observed, expected)

    def test_sequence_collection_to_fastq_kwargs_passed(self):
        for components, kwargs_expected_fp in self.valid_files:
            for kwargs, expected_fp in kwargs_expected_fp:
                obj = SequenceCollection([
                    NucleotideSequence(c[2], id=c[0], description=c[1],
                                       quality=c[3]) for c in components])

                fh = StringIO()
                _sequence_collection_to_fastq(obj, fh, **kwargs)
                observed = fh.getvalue()
                fh.close()

                with open(expected_fp, 'U') as f:
                    expected = f.read()

                self.assertEqual(observed, expected)

    def test_alignment_to_fastq_kwargs_passed(self):
        for components, kwargs_expected_fp in self.valid_files:
            for kwargs, expected_fp in kwargs_expected_fp:
                obj = Alignment([
                    ProteinSequence(c[2], id=c[0], description=c[1],
                                    quality=c[3]) for c in components])

                fh = StringIO()
                _alignment_to_fastq(obj, fh, **kwargs)
                observed = fh.getvalue()
                fh.close()

                with open(expected_fp, 'U') as f:
                    expected = f.read()

                self.assertEqual(observed, expected)


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
                fh = StringIO()

                # will issue warning when truncating quality scores
                with warnings.catch_warnings(record=True):
                    warnings.simplefilter("ignore")
                    _generator_to_fastq(read_gen, fh, **to_kwargs)

                obs = fh.getvalue()
                fh.close()

                with open(to_fp, 'U') as fh:
                    exp = fh.read()
                self.assertEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()
