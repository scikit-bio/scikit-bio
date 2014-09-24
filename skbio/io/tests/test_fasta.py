# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils.six import StringIO

from unittest import TestCase, main

from skbio import BiologicalSequence
from skbio.io import FASTAFormatError
from skbio.io.fasta import (_fasta_sniffer, _fasta_to_generator,
                            _generator_to_fasta)
from skbio.util import get_data_path


class FASTATests(TestCase):
    def setUp(self):
        def empty_gen():
            raise StopIteration()
            yield

        bseq1 = BiologicalSequence('ACGT-acgt.', id='seq1',
                                   description='desc1', quality=range(10))
        bseq2 = BiologicalSequence('A', id='seq2')
        bseq3 = BiologicalSequence('AACGGuA', description='desc3')
        blank_seq = BiologicalSequence('')

        def single_seq_gen():
            yield bseq1

        def multi_seq_gen():
            for seq in bseq1, bseq2, bseq3:
                yield seq

        self.objs_fps = map(lambda e: (e[0], get_data_path(e[1])), [
            (empty_gen(), 'empty'),
            (single_seq_gen(), 'fasta_single_seq'),
            (multi_seq_gen(), 'fasta_multi_seq'),
        ])

        def blank_seq_gen():
            for seq in bseq1, blank_seq:
                yield seq

        # generators that cannot be written in fasta format, paired with their
        # expected error message regexps
        self.invalid_objs = [
            (blank_seq_gen(), 'number 2.*empty'),
        ]

    def test_generator_to_fasta(self):
        for obj, fp in self.objs_fps:
            fh = StringIO()
            _generator_to_fasta(obj, fh)
            obs = fh.getvalue()
            fh.close()

            with open(fp, 'U') as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

    def test_generator_to_fasta_invalid_data(self):
        for invalid_obj, error_msg_regexp in self.invalid_objs:
            fh = StringIO()
            with self.assertRaisesRegexp(FASTAFormatError, error_msg_regexp):
                _generator_to_fasta(invalid_obj, fh)
            fh.close()


if __name__ == '__main__':
    main()
