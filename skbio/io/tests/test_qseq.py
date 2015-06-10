# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from future.builtins import zip

import unittest
from functools import partial

from skbio import SequenceCollection, Sequence, DNA, RNA, Protein
from skbio import read
from skbio.util import get_data_path
from skbio.io import QSeqFormatError
from skbio.io.qseq import (_qseq_to_generator,
                           _qseq_to_sequence_collection, _qseq_sniffer)
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
                ('sanger_1:3:34:-30:30#0/2',
                 'ACGTACGTACGTACGTACGTACGTACTTTTTTTTTTACGTACGTACGTACGT'
                 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC', [
                     26, 26, 29, 31, 33, 34, 36, 37, 38, 39, 41, 42,
                     43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
                     55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
                     67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78,
                     79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
                     91, 92, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93,
                     93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93,
                     93, 93, 93, 93, 93, 93, 93, 93, 93, 93])
            ]),

            (get_data_path('qseq_multi_seq_illumina1.3'), [
                {'variant': 'illumina1.3'},
                {'phred_offset': 64}
            ], [
                ('illumina_1:3:34:-30:30#0/1', 'ACG....ACGTAC', [
                    50, 53, 2, 2, 2, 2, 50, 2, 3, 5, 6, 7, 8]),
                ('illumina_1:3:35:-30:30#0/2', 'ACGTA.AATAAAC', [
                    39, 37, 20, 33, 1, 33, 38, 40, 55, 49, 1, 1, 38])
            ]),

            (get_data_path('qseq_multi_seq_illumina1.3'), [
                {'variant': 'illumina1.3', 'filter': False, 'seq_num': 1},
                {'phred_offset': 64, 'filter': False, 'seq_num': 2},
                {'variant': 'illumina1.3', 'filter': False, 'seq_num': 3,
                 'constructor': partial(Protein, validate=False)},
                {'phred_offset': 64, 'filter': False, 'seq_num': 4,
                 'constructor': partial(DNA, validate=False)},
            ], [
                ('illumina_1:3:34:-30:30#0/1', 'ACG....ACGTAC', [
                    50, 53, 2, 2, 2, 2, 50, 2, 3, 5, 6, 7, 8]),
                ('illumina_1:3:34:30:-30#0/1', 'CGGGCATTGCA', [
                    3, 7, 7, 7, 3, 33, 51, 36, 7, 3, 1]),
                ('illumina_1:3:35:-30:30#0/2', 'ACGTA.AATAAAC', [
                    39, 37, 20, 33, 1, 33, 38, 40, 55, 49, 1, 1, 38]),
                ('illumina_1:3:35:30:-30#0/3', 'CATTTAGGA.TGCA', [
                    52, 42, 38, 44, 43, 1, 6, 46, 43, 11, 39, 40, 54, 13])
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
                        c[1],
                        metadata={'id': c[0]},
                        positional_metadata={
                            'quality': np.array(c[2], dtype=np.uint8)})
                    for c in components]

                observed = list(_qseq_to_generator(valid, **kwarg))
                self.assertEqual(len(expected), len(observed))
                for o, e in zip(observed, expected):
                    self.assertTrue(o.equals(e))


class TestQSeqToSequenceCollection(TestQSeqBase):
    def setUp(self):
        super(TestQSeqToSequenceCollection, self).setUp()
        self.valid_files += [
            (get_data_path('empty'), [{}, {'variant': 'sanger'}],
             SequenceCollection([]))
        ]

    def test_invalid_files(self):
        for invalid, kwargs, errors, etype in self.invalid_files:
            with self.assertRaises(etype) as cm:
                for kwarg in kwargs:
                    _drop_kwargs(kwarg, 'seq_num')
                    _qseq_to_sequence_collection(invalid, **kwarg)
            for e in errors:
                self.assertIn(e, str(cm.exception))

    def test_valid_files(self):
        for valid, kwargs, components in self.valid_files:
            for kwarg in kwargs:
                _drop_kwargs(kwarg, 'seq_num')
                constructor = kwarg.get('constructor', Sequence)
                expected = SequenceCollection([
                    constructor(
                        c[1],
                        metadata={'id': c[0]},
                        positional_metadata={
                            'quality': np.array(c[2], dtype=np.uint8)})
                    for c in components])

                observed = _qseq_to_sequence_collection(valid, **kwarg)
                # TODO remove when #656 is resolved
                self.assertEqual(observed, expected)
                for o, e in zip(observed, expected):
                    self.assertTrue(o.equals(e))


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
        for constructor in [partial(Sequence), partial(DNA, validate=False),
                            partial(RNA, validate=False),
                            partial(Protein, validate=False)]:
            for valid, kwargs, components in self.valid_files:
                for kwarg in kwargs:
                    _drop_kwargs(kwarg, 'constructor', 'filter')

                    seq_num = kwarg.get('seq_num', 1)
                    c = components[seq_num - 1]
                    expected = constructor(
                        c[1],
                        metadata={'id': c[0]},
                        positional_metadata={
                            'quality': np.array(c[2], np.uint8)})

                    observed = read(valid, into=constructor.func,
                                    format='qseq', verify=False, **kwarg)
                    self.assertTrue(observed.equals(expected))


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
