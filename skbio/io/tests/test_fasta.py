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

from skbio import BiologicalSequence, DNA, RNA, Protein
from skbio.io import FASTAFormatError
from skbio.io.fasta import (_fasta_sniffer, _fasta_to_generator,
                            _generator_to_fasta)
from skbio.util import get_data_path


class FASTATests(TestCase):
    def setUp(self):
        def empty_gen():
            raise StopIteration()
            yield

        bio_seq1 = BiologicalSequence('ACGT-acgt.', id='seq1',
                                      description='desc1', quality=range(10))
        bio_seq2 = BiologicalSequence('A', id='seq2')
        bio_seq3 = BiologicalSequence('AACGGuA', description='desc3')
        dna_seq = DNA('ACGTTGCAccGG')
        rna_seq = RNA('ACGUU', quality=[42] * 5)
        prot_seq = Protein('pQqqqPPQQQ', id='proteinseq',
                           description='a very detailed description')

        def single_seq_gen():
            yield bio_seq1

        # multiple sequences of mixed types, lengths, and metadata. lengths are
        # chosen to exercise various splitting cases when testing max_width
        def multi_seq_gen():
            for seq in (bio_seq1, bio_seq2, bio_seq3, dna_seq, rna_seq,
                        prot_seq):
                yield seq

        # store sequence generator to serialize, writer kwargs (if any), and
        # filepath of expected results
        self.objs_fps = map(lambda e: (e[0], e[1], get_data_path(e[2])), [
            (empty_gen(), {}, 'empty'),
            (single_seq_gen(), {}, 'fasta_single_seq'),
            (single_seq_gen(), {'max_width': 1}, 'fasta_max_width_1'),
            (multi_seq_gen(), {}, 'fasta_multi_seq'),
            (multi_seq_gen(), {'max_width': 5}, 'fasta_max_width_5'),
        ])

        blank_seq = BiologicalSequence('')

        def blank_seq_gen():
            for seq in bio_seq1, blank_seq:
                yield seq

        # generators that cannot be written in fasta format, paired with their
        # expected error message regexps
        self.invalid_objs = [
            (blank_seq_gen(), 'number 2.*empty'),
        ]

    def test_generator_to_fasta(self):
        for obj, kwargs, fp in self.objs_fps:
            fh = StringIO()
            _generator_to_fasta(obj, fh, **kwargs)
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
