# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

from skbio import (BiologicalSequence)
from skbio.io import FASTQFormatError
from skbio.io.fastq import (_fastq_to_generator, _fastq_sniffer,
                            ascii_to_phred33, ascii_to_phred64)
from skbio.util import get_data_path


class FASTQSnifferTests(TestCase):
    def setUp(self):
        self.positive_fps_33 = map(get_data_path, [
            'fastq_multi_seq33',
        ])

        self.positive_fps_64 = map(get_data_path, [
            'fastq_single_seq',
            'fastq_multi_seq',
            'fastq_multi_seq64'
        ])

        self.negative_fps = map(get_data_path, [
            'empty',
            'fastq_invalid_missing_header',
            'fastq_invalid_missing_qual_header',
        ])

    def test_positives(self):
        for fp in self.positive_fps_33:
            self.assertEqual(_fastq_sniffer(fp), (True, {'phred_offset':33}))
        for fp in self.positive_fps_64:
            self.assertEqual(_fastq_sniffer(fp), (True, {'phred_offset':64}))

    def test_negatives(self):
        for fp in self.negative_fps:
            self.assertEqual(_fastq_sniffer(fp), (False, {}))


class FASTQReaderTests(TestCase):
    def setUp(self):
        self.empty = ([], {}, map(get_data_path, ['empty']))
        self.single = (
            [BiologicalSequence('AACACCAAACTTCTCCACCACGTGAGCTACAAAAG',
                                id='seq1',
                                quality=ascii_to_phred64(
                                    '````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF'))],
            {'phred_offset': 64},
            map(get_data_path, ['fastq_single_seq'])
        )

        self.multi = (
            [BiologicalSequence('ACGTT', id='k',
                                description='ok',
                                quality=ascii_to_phred64(
                                    r'^\\\Y')),
             BiologicalSequence('AACGGuA', id='desc3',
                                quality=ascii_to_phred64(
                                    'bbbbbbb')),
             BiologicalSequence('AcGtUTu', id='seq1',
                                quality=ascii_to_phred64(
                                    'bbkqwbo')),
             BiologicalSequence('ACGTTGCAccGG',
                                quality=ascii_to_phred64(
                                    'BBBBBBBBBBBB')),
             BiologicalSequence('ACGUU',
                                quality=ascii_to_phred64(
                                    r'^\\\Y'))],
            {'phred_offset': 64},
            map(get_data_path, ['fastq_multi_seq'])
        )

        self.multi33 = (
            [BiologicalSequence('GTTGCTTCTGGCGTGGGTGGGGGGG',
                                id=r'EAS54_6_R1_2_1_443_348',
                                quality=ascii_to_phred33(
                                    r';;;;;;;;;;;9;7;;.7;393333')),
             BiologicalSequence('GATTTGGGGTTCAAAGCAGTATCGATCAAA'
                                'TAGTAAATCCATTTGTTCAACTCACAGTTT',
                                id=r'SEQ_ID',
                                quality=ascii_to_phred33(
                                    r"!''*((((***+))%%%++)(%%%%).1**"
                                    "*-+*''))**55CCF>>>>>>CCCCCCC65"))],
            {'phred_offset': 33},
            map(get_data_path, ['fastq_multi_seq33']))

        self.multi64 = (
            [BiologicalSequence('AACACCAAACTTCTCCACCACGTGAGCTACAAAAG',
                                id=r'GAPC_0015:6:1:1259:10413#0/1',
                                quality=ascii_to_phred64(
                                    r'````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF')),
             BiologicalSequence('TATGTATATATAACATATACATATATACATACATA',
                                id=r'GAPC_0015:6:1:1283:11957#0/1',
                                quality=ascii_to_phred64(
                                    r']KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb')),
             BiologicalSequence('TCAGTTTTCCTCGCCATATTTCACGTCCTAAAGCG',
                                id=r'GAPC_0015:6:1:1284:10484#0/1',
                                quality=ascii_to_phred64(
                                    r'UM_]]U_]Z_Y^\^^``Y]`^SZ]\Ybb`^_LbL_')),
             BiologicalSequence('TGTGCCTATGGAAGCAGTTCTAGGATCCCCTAGAA',
                                id=r'GAPC_0015:6:1:1287:17135#0/1',
                                quality=ascii_to_phred64(
                                    r'^aacccL\ccc\c\cTKTS]KZ\]]I\[Wa^T`^K')),
             BiologicalSequence('AAAGAAAGGAAGAAAAGAAAAAGAAACCCGAGTTA',
                                id=r'GAPC_0015:6:1:1293:3171#0/1',
                                quality=ascii_to_phred64(
                                    r'b`bbbU_[YYcadcda_LbaaabWbaacYcc`a^c')),
             BiologicalSequence('TAATGCCAAAGAAATATTTCCAAACTACATGCTTA',
                                id=r'GAPC_0015:6:1:1297:10729#0/1',
                                quality=ascii_to_phred64(
                                    r'T\ccLbb``bacc]_cacccccLccc\ccTccYL^'))],
            {'phred_offset': 64},
            map(get_data_path, ['fastq_multi_seq64']))

        self.invalid_fps = map(lambda e: (get_data_path(e[0]), e[1]), [
            ('fastq_invalid_bad_id', 'ID mismatch'),
            ('fastq_invalid_missing_header', 'Bad FASTQ format'),
            ('fastq_invalid_bad_qual', 'Failed qual conversion')
        ])

        self.invalid = []

    def test_fastq_to_generator_valid_files(self):

        for exp, kwargs, fps in (self.multi, self.empty, self.single,
                                 self.multi33, self.multi64):
            for fp in fps:
                obs = list(_fastq_to_generator(fp, **kwargs))

                self.assertEqual(len(obs), len(exp))
                for o, e in zip(obs, exp):
                    self.assertTrue(o.equals(e))

    def test_fastq_to_generator_invalid_files(self):
        for fp, error_msg_regex in self.invalid_fps:
            with self.assertRaisesRegexp(FASTQFormatError, error_msg_regex):
                list(_fastq_to_generator(fp, strict=True))


if __name__ == '__main__':
    main()
