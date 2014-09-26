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

from skbio import BiologicalSequence, NucleotideSequence, DNA, RNA, Protein
from skbio import SequenceCollection, Alignment
from skbio.io import FASTAFormatError
from skbio.io.fasta import (_generator_to_fasta, _biological_sequence_to_fasta,
                            _nucleotide_sequence_to_fasta,
                            _dna_sequence_to_fasta, _rna_sequence_to_fasta,
                            _protein_sequence_to_fasta,
                            _sequence_collection_to_fasta, _alignment_to_fasta)
from skbio.util import get_data_path


class FASTATests(TestCase):
    def setUp(self):
        def empty_gen():
            raise StopIteration()
            yield

        self.bio_seq1 = BiologicalSequence(
            'ACGT-acgt.', id='seq1', description='desc1', quality=range(10))
        self.bio_seq2 = BiologicalSequence('A', id='seq2')
        self.bio_seq3 = BiologicalSequence('AACGGuA', description='desc3')
        self.nuc_seq = NucleotideSequence('AcGtUTu')
        self.dna_seq = DNA('ACGTTGCAccGG')
        self.rna_seq = RNA('ACGUU', quality=[42] * 5)
        self.prot_seq = Protein(
            'pQqqqPPQQQ', id='proteinseq',
            description='\ndetailed\ndescription \t\twith  new\n\nlines\n\n\n')

        seqs = [self.rna_seq, self.bio_seq1, self.prot_seq]
        self.seq_coll = SequenceCollection(seqs)
        self.align = Alignment(seqs)

        def single_seq_gen():
            yield self.bio_seq1

        def newline_description_gen():
            yield self.prot_seq
            yield DNA('AGGAGAATA', id='foo', description='\n\n\n\n')

        # multiple sequences of mixed types, lengths, and metadata. lengths are
        # chosen to exercise various splitting cases when testing max_width
        def multi_seq_gen():
            for seq in (self.bio_seq1, self.bio_seq2, self.bio_seq3,
                        self.nuc_seq, self.dna_seq, self.rna_seq,
                        self.prot_seq):
                yield seq

        # store sequence generator to serialize, writer kwargs (if any), and
        # filepath of expected results
        self.objs_fps = map(lambda e: (e[0], e[1], get_data_path(e[2])), [
            (empty_gen(), {}, 'empty'),
            (single_seq_gen(), {}, 'fasta_single_seq'),
            (single_seq_gen(), {'max_width': 1}, 'fasta_max_width_1'),
            (multi_seq_gen(), {}, 'fasta_multi_seq'),
            (multi_seq_gen(), {'max_width': 5}, 'fasta_max_width_5'),
            (newline_description_gen(),
             {'description_newline_replacement': ':-)'},
             'fasta_description_newline_replacement_multi_char'),
            (newline_description_gen(),
             {'description_newline_replacement': ''},
             'fasta_description_newline_replacement_empty_str'),
        ])

        self.blank_seq = BiologicalSequence('')

        def blank_seq_gen():
            for seq in self.bio_seq1, self.blank_seq:
                yield seq

        # generators or parameter combos that cannot be written in fasta
        # format, paired with kwargs (if any), error type, and expected error
        # message regexp
        self.invalid_objs = [
            (blank_seq_gen(), {}, FASTAFormatError, '2nd.*empty'),
            (single_seq_gen(),
             {'max_width': 0}, ValueError, 'n=0'),
            (multi_seq_gen(), {'id_whitespace_replacement': '-\n_'},
             FASTAFormatError, 'Newline character'),
            (multi_seq_gen(), {'description_newline_replacement': '-.-\n'},
             FASTAFormatError, 'Newline character'),
        ]

    # extensive tests for generator -> fasta writer since it is used by all
    # other object -> fasta writers

    def test_generator_to_fasta(self):
        for obj, kwargs, fp in self.objs_fps:
            fh = StringIO()
            _generator_to_fasta(obj, fh, **kwargs)
            obs = fh.getvalue()
            fh.close()

            with open(fp, 'U') as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

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
        # Store writer function, sequence object to write, and expected
        # filepaths for each of the invoked keyword arguments (see below).
        test_data = (
            (_biological_sequence_to_fasta, self.bio_seq1,
             ('fasta_single_seq', 'fasta_max_width_1')),
            (_nucleotide_sequence_to_fasta, self.nuc_seq,
             ('fasta_single_nuc_seq', 'fasta_single_nuc_seq_max_width_1')),
            (_dna_sequence_to_fasta, self.dna_seq,
             ('fasta_single_dna_seq', 'fasta_single_dna_seq_max_width_1')),
            (_rna_sequence_to_fasta, self.rna_seq,
             ('fasta_single_rna_seq', 'fasta_single_rna_seq_max_width_1')),
            (_protein_sequence_to_fasta, self.prot_seq,
             ('fasta_single_prot_seq', 'fasta_single_prot_seq_max_width_1')))

        for fn, obj, fps in test_data:
            for kw, fp in zip(({}, {'max_width': 1}), fps):
                fh = StringIO()
                fn(obj, fh, **kw)
                obs = fh.getvalue()
                fh.close()

                with open(get_data_path(fp), 'U') as fh:
                    exp = fh.read()

                self.assertEqual(obs, exp)

    def test_any_sequences_to_fasta(self):
        for fn, obj in ((_sequence_collection_to_fasta, self.seq_coll),
                        (_alignment_to_fasta, self.align)):
            for kw, fp in (({}, 'fasta_3_seqs'),
                           ({'max_width': 3}, 'fasta_3_seqs_max_width_3')):
                fh = StringIO()
                fn(obj, fh, **kw)
                obs = fh.getvalue()
                fh.close()

                with open(get_data_path(fp), 'U') as fh:
                    exp = fh.read()

                self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
