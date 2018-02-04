# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from skbio import Sequence, DNA, RNA, Protein
from skbio import read
from skbio.util import get_data_path
from skbio.io import QSeqFormatError
from skbio.io.format.qseq import _qseq_to_generator, _qseq_sniffer
import numpy as np


def _drop_kwargs(kwargs, *args):
    for arg in args:
        if arg in kwargs:
            kwargs.pop(arg)


class TestQSeqBase(unittest.TestCase):
    def setUp(self):
        self.valid_files = [
            (get_data_path('qseq_single_seq_sanger'), [
                {'variant': 'sanger'},
                {'phred_offset': 33},
            ], [
                {'id': 'sanger_1:3:34:-30:30#0/2',
                 'sequence': 'ACGTACGTACGTACGTACGTACGTACTTTTTTTTTTACGTACGTACG'
                             'TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC',
                 'quality': [26, 26, 29, 31, 33, 34, 36, 37, 38, 39, 41, 42,
                             43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
                             55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
                             67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78,
                             79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
                             91, 92, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93,
                             93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93,
                             93, 93, 93, 93, 93, 93, 93, 93, 93, 93],
                 'machine_name': 'sanger',
                 'run_number': 1,
                 'lane_number': 3,
                 'tile_number': 34,
                 'x': -30,
                 'y': 30,
                 'index': 0,
                 'read_number': 2}
            ]),

            (get_data_path('qseq_multi_seq_illumina1.3'), [
                {'variant': 'illumina1.3'},
                {'phred_offset': 64}
            ], [
                {'id': 'illumina_1:3:34:-30:30#0/1',
                 'sequence': 'ACG....ACGTAC',
                 'quality': [50, 53, 2, 2, 2, 2, 50, 2, 3, 5, 6, 7, 8],
                 'machine_name': 'illumina',
                 'run_number': 1,
                 'lane_number': 3,
                 'tile_number': 34,
                 'x': -30,
                 'y': 30,
                 'index': 0,
                 'read_number': 1},
                {'id': 'illumina_1:3:35:-30:30#0/2',
                 'sequence': 'ACGTA.AATAAAC',
                 'quality': [39, 37, 20, 33, 1, 33, 38, 40, 55, 49, 1, 1, 38],
                 'machine_name': 'illumina',
                 'run_number': 1,
                 'lane_number': 3,
                 'tile_number': 35,
                 'x': -30,
                 'y': 30,
                 'index': 0,
                 'read_number': 2}
            ]),

            (get_data_path('qseq_multi_seq_illumina1.3'), [
                {'variant': 'illumina1.3', 'filter': False, 'seq_num': 1},
                {'phred_offset': 64, 'filter': False, 'seq_num': 2},
                {'variant': 'illumina1.3', 'filter': False, 'seq_num': 3,
                 'constructor': Protein},
                {'phred_offset': 64, 'filter': False, 'seq_num': 4,
                 'constructor': DNA},
            ], [
                {'id': 'illumina_1:3:34:-30:30#0/1',
                 'sequence': 'ACG....ACGTAC',
                 'quality': [50, 53, 2, 2, 2, 2, 50, 2, 3, 5, 6, 7, 8],
                 'machine_name': 'illumina',
                 'run_number': 1,
                 'lane_number': 3,
                 'tile_number': 34,
                 'x': -30,
                 'y': 30,
                 'index': 0,
                 'read_number': 1},
                {'id': 'illumina_1:3:34:30:-30#0/1',
                 'sequence': 'CGGGCATTGCA',
                 'quality': [3, 7, 7, 7, 3, 33, 51, 36, 7, 3, 1],
                 'machine_name': 'illumina',
                 'run_number': 1,
                 'lane_number': 3,
                 'tile_number': 34,
                 'x': 30,
                 'y': -30,
                 'index': 0,
                 'read_number': 1},
                {'id': 'illumina_1:3:35:-30:30#0/2',
                 'sequence': 'ACGTA.AATAAAC',
                 'quality': [39, 37, 20, 33, 1, 33, 38, 40, 55, 49, 1, 1, 38],
                 'machine_name': 'illumina',
                 'run_number': 1,
                 'lane_number': 3,
                 'tile_number': 35,
                 'x': -30,
                 'y': 30,
                 'index': 0,
                 'read_number': 2},
                {'id': 'illumina_1:3:35:30:-30#0/3',
                 'sequence': 'CATTTAGGA.TGCA',
                 'quality': [52, 42, 38, 44, 43, 1, 6, 46, 43, 11, 39, 40, 54,
                             13],
                 'machine_name': 'illumina',
                 'run_number': 1,
                 'lane_number': 3,
                 'tile_number': 35,
                 'x': 30,
                 'y': -30,
                 'index': 0,
                 'read_number': 3}
            ])
        ]

        self.invalid_files = [
            (get_data_path('whitespace_only'), [
                {},
                {'variant': 'sanger'}
            ], [
                'blank line',
            ], QSeqFormatError),

            (get_data_path('tsv_10_fields'), [
                {},
                {'variant': 'sanger'},
                {'variant': 'solexa'}
            ], [
                'read',
                '[1, 3]'
            ], QSeqFormatError),

            (get_data_path('tsv_8_fields'), [
                {},
                {'variant': 'sanger'},
                {'variant': 'solexa'}
            ], [
                '8',
                '10 or 11'
            ], QSeqFormatError),


            (get_data_path('qseq_invalid_filter'), [
                {},
                {'phred_offset': 33},
                {'variant': 'solexa'},
                {'variant': 'illumina1.3'},
                {'variant': 'illumina1.8'}
            ], [
                'filter',
                '0 or 1',
            ], QSeqFormatError),

            (get_data_path('qseq_invalid_read'), [
                {},
                {'phred_offset': 33},
                {'variant': 'solexa'},
                {'variant': 'illumina1.3'},
                {'variant': 'illumina1.8'}
            ], [
                'read',
                '[1, 3]',
            ], QSeqFormatError),

            (get_data_path('qseq_invalid_x'), [
                {},
                {'phred_offset': 33},
                {'variant': 'solexa'},
                {'variant': 'illumina1.3'},
                {'variant': 'illumina1.8'}
            ], [
                'x',
                'integer',
            ], QSeqFormatError),

            (get_data_path('qseq_invalid_y'), [
                {},
                {'phred_offset': 33},
                {'variant': 'solexa'},
                {'variant': 'illumina1.3'},
                {'variant': 'illumina1.8'}
            ], [
                'y',
                'integer',
            ], QSeqFormatError),

            (get_data_path('qseq_invalid_lane'), [
                {},
                {'phred_offset': 33},
                {'variant': 'solexa'},
                {'variant': 'illumina1.3'},
                {'variant': 'illumina1.8'}
            ], [
                'lane',
                'positive integer',
            ], QSeqFormatError),

            (get_data_path('qseq_invalid_tile'), [
                {},
                {'phred_offset': 33},
                {'variant': 'solexa'},
                {'variant': 'illumina1.3'},
                {'variant': 'illumina1.8'}
            ], [
                'tile',
                'positive integer',
            ], QSeqFormatError)
        ]


class TestQSeqToGenerator(TestQSeqBase):

    def setUp(self):
        super(TestQSeqToGenerator, self).setUp()
        self.valid_files += [
            (get_data_path('empty'), [{}, {'variant': 'sanger'}], [])
        ]

        self.invalid_files += [
            (get_data_path('qseq_single_seq_sanger'), [
                {'variant': 'illumina1.3'},
                {'variant': 'illumina1.8'}
            ], [
                'out of range',
                '[0, 62]'
            ], ValueError)
        ]

    def test_invalid_files(self):
        for invalid, kwargs, errors, etype in self.invalid_files:
            with self.assertRaises(etype) as cm:
                for kwarg in kwargs:
                    _drop_kwargs(kwarg, 'seq_num')
                    list(_qseq_to_generator(invalid, **kwarg))
            for e in errors:
                self.assertIn(e, str(cm.exception))

    def test_valid_files(self):
        for valid, kwargs, components in self.valid_files:
            for kwarg in kwargs:
                _drop_kwargs(kwarg, 'seq_num')
                constructor = kwarg.get('constructor', Sequence)
                expected = [
                    constructor(
                        c['sequence'],
                        metadata={'id': c['id'],
                                  'machine_name': c['machine_name'],
                                  'run_number': c['run_number'],
                                  'lane_number': c['lane_number'],
                                  'tile_number': c['tile_number'],
                                  'x': c['x'],
                                  'y': c['y'],
                                  'index': c['index'],
                                  'read_number': c['read_number']},
                        positional_metadata={
                            'quality': np.array(c['quality'], dtype=np.uint8)})
                    for c in components]

                observed = list(_qseq_to_generator(valid, **kwarg))
                self.assertEqual(len(expected), len(observed))
                for o, e in zip(observed, expected):
                    self.assertEqual(o, e)


class TestQSeqToSequences(TestQSeqBase):
    def test_invalid_files(self):
        for constructor in [Sequence, DNA, RNA, Protein]:
            for invalid, kwargs, errors, etype in self.invalid_files:
                with self.assertRaises(etype) as cm:
                    for kwarg in kwargs:
                        _drop_kwargs(kwarg, 'constructor', 'filter')

                        read(invalid, format='qseq', verify=False,
                             into=constructor, **kwarg)
                for e in errors:
                    self.assertIn(e, str(cm.exception))

    def test_valid_files(self):
        for constructor in [Sequence, DNA, RNA, Protein]:
            for valid, kwargs, components in self.valid_files:
                for observed_kwargs in kwargs:
                    expected_kwargs = {}
                    # Currently not validating the alphabet for qseq
                    # files that are read in for this test.
                    if hasattr(constructor, 'alphabet'):
                        observed_kwargs['validate'] = False
                        expected_kwargs['validate'] = False
                    _drop_kwargs(observed_kwargs, 'constructor', 'filter')

                    seq_num = observed_kwargs.get('seq_num', 1)
                    c = components[seq_num - 1]
                    expected = constructor(
                        c['sequence'],
                        metadata={'id': c['id'],
                                  'machine_name': c['machine_name'],
                                  'run_number': c['run_number'],
                                  'lane_number': c['lane_number'],
                                  'tile_number': c['tile_number'],
                                  'x': c['x'],
                                  'y': c['y'],
                                  'index': c['index'],
                                  'read_number': c['read_number']},
                        positional_metadata={
                            'quality': np.array(c['quality'], np.uint8)},
                        **expected_kwargs)

                    observed = read(valid, into=constructor,
                                    format='qseq', verify=False,
                                    **observed_kwargs)
                    self.assertEqual(observed, expected)


class TestQSeqSniffer(TestQSeqBase):

    def setUp(self):
        super(TestQSeqSniffer, self).setUp()
        self.invalid_files += [
            (get_data_path('empty'), None, None, None)
        ]

    def test_qseq_sniffer_valid_files(self):
        for valid, _, _ in self.valid_files:
            self.assertEqual(_qseq_sniffer(valid), (True, {}))

    def test_qseq_sniffer_invalid_files(self):
        for invalid, _, _, _ in self.invalid_files:
            self.assertEqual(_qseq_sniffer(invalid), (False, {}))


if __name__ == '__main__':
    unittest.main()
