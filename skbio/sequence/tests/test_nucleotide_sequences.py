# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np

from skbio import DNA, RNA, Protein, GeneticCode
from skbio.sequence._nucleotide_mixin import NucleotideMixin
from skbio.sequence import GrammaredSequence
from skbio.util import classproperty
from skbio.metadata import IntervalMetadata


# This file contains tests for functionality of sequence types which implement
# NucleotideMixin. Currently this means DNA and RNA. These types are so
# similar that the testing logic can be shared and parameterized across
# different test data.

class TestNucleotideSequence(unittest.TestCase):
    def setUp(self):
        self.sequence_kinds = frozenset([
            str,
            lambda s: np.frombuffer(s.encode('ascii'), dtype='|S1'),
            lambda s: np.frombuffer(s.encode('ascii'), dtype=np.uint8)])

        dna_str = 'ACGTMRWSYKVHDBN.-'
        dna_comp_str = 'TGCAKYWSRMBDHVN.-'
        dna_rev_comp_str = '-.NVHDBMRSWYKACGT'
        rna_str = 'ACGUMRWSYKVHDBN.-'
        rna_comp_str = 'UGCAKYWSRMBDHVN.-'
        rna_rev_comp_str = '-.NVHDBMRSWYKACGU'
        qual = tuple(range(len(dna_str)))

        self.dna = (DNA, dna_str)
        self.rna = (RNA, rna_str)

        dna_comp = self.dna + (dna_comp_str,)
        rna_comp = self.rna + (rna_comp_str,)

        dna_comp_qual = dna_comp + (qual,)
        rna_comp_qual = rna_comp + (qual,)
        self.all_combos_comp_qual = (dna_comp_qual, rna_comp_qual)

        dna_rev_comp = self.dna + (dna_rev_comp_str,)
        rna_rev_comp = self.rna + (rna_rev_comp_str,)
        self.all_combos_rev_comp = (dna_rev_comp, rna_rev_comp)

        dna_rev_comp_qual = dna_rev_comp + (qual,)
        rna_rev_comp_qual = rna_rev_comp + (qual,)
        self.all_combos_rev_comp_qual = \
            (dna_rev_comp_qual, rna_rev_comp_qual)

    def test_instantiation_with_no_implementation(self):
        class NucleotideSequenceSubclassNoImplementation(NucleotideMixin):
            pass

        with self.assertRaises(TypeError) as cm:
            NucleotideSequenceSubclassNoImplementation()

        self.assertIn("abstract class", str(cm.exception))
        self.assertIn("complement_map", str(cm.exception))

    # TODO: remove when nondegenerate_chars is removed
    def test_nondegenerate_chars(self):
        dna = (DNA, "ACGT")
        rna = (RNA, "ACGU")
        for constructor, nondegenerate in (dna, rna):
            exp = set(nondegenerate)
            self.assertEqual(constructor('').nondegenerate_chars, exp)
            self.assertEqual(constructor.nondegenerate_chars, exp)

    def test_definite_chars(self):
        dna = (DNA, "ACGT")
        rna = (RNA, "ACGU")
        for constructor, definite_char in (dna, rna):
            exp = set(definite_char)
            self.assertEqual(constructor('').definite_chars, exp)
            self.assertEqual(constructor.definite_chars, exp)

    def test_degenerate_map(self):
        dna_exp = (DNA, {
            'B': set(['C', 'T', 'G']), 'D': set(['A', 'T', 'G']),
            'H': set(['A', 'C', 'T']), 'K': set(['T', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'T', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'T']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'T'])
        })

        rna_exp = (RNA, {
            'B': set(['C', 'U', 'G']), 'D': set(['A', 'U', 'G']),
            'H': set(['A', 'C', 'U']), 'K': set(['U', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'U', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'U']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'U'])
        })

        for constructor, degenerate in (dna_exp, rna_exp):
            self.assertEqual(constructor('').degenerate_map, degenerate)
            self.assertEqual(constructor.degenerate_map, degenerate)

    def test_complement_map(self):
        dna_exp = (DNA, {
            '-': '-', '.': '.', 'A': 'T', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'T': 'A', 'W': 'W', 'V': 'B', 'Y': 'R'
        })

        rna_exp = (RNA, {
            '-': '-', '.': '.', 'A': 'U', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'U': 'A', 'W': 'W', 'V': 'B', 'Y': 'R'
        })

        for constructor, comp_map in (dna_exp, rna_exp):
            self.assertEqual(constructor('').complement_map, comp_map)
            self.assertEqual(constructor.complement_map, comp_map)

            # immutable
            constructor.complement_map['A'] = 'X'
            constructor.complement_map['C'] = 'W'
            self.assertEqual(constructor.complement_map, comp_map)
            with self.assertRaises(AttributeError):
                constructor('').complement_map = {'W': 'X'}

    def test_translate_ncbi_table_id(self):
        for seq in RNA('AAAUUUAUGCAU'), DNA('AAATTTATGCAT'):
            # default
            obs = seq.translate()
            self.assertEqual(obs, Protein('KFMH'))

            obs = seq.translate(9)
            self.assertEqual(obs, Protein('NFMH'))

    def test_translate_genetic_code_object(self):
        gc = GeneticCode('M' * 64, '-' * 64)
        for seq in RNA('AAAUUUAUGCAU'), DNA('AAATTTATGCAT'):
            obs = seq.translate(gc)
            self.assertEqual(obs, Protein('MMMM'))

    def test_translate_passes_parameters_through(self):
        exp = Protein('MW')
        for seq in RNA('UAAAUUGUGGUAA'), DNA('TAAATTGTGGTAA'):
            # mix of args and kwargs
            obs = seq.translate(13, reading_frame=2, start='require',
                                stop='require')
            self.assertEqual(obs, exp)

            # kwargs only
            obs = seq.translate(genetic_code=13, reading_frame=2,
                                start='require', stop='require')
            self.assertEqual(obs, exp)

            # args only
            obs = seq.translate(13, 2, 'require', 'require')
            self.assertEqual(obs, exp)

    def test_translate_preserves_metadata(self):
        metadata = {'foo': 'bar', 'baz': 42}
        positional_metadata = {'foo': range(3)}
        for seq in (RNA('AUG', metadata=metadata,
                        positional_metadata=positional_metadata),
                    DNA('ATG', metadata=metadata,
                        positional_metadata=positional_metadata)):
            obs = seq.translate()
            # metadata retained, positional metadata dropped
            self.assertEqual(obs,
                             Protein('M', metadata={'foo': 'bar', 'baz': 42}))

    def test_translate_invalid_id(self):
        for seq in RNA('AUG'), DNA('ATG'):
            with self.assertRaisesRegex(ValueError, r'table_id.*42'):
                seq.translate(42)

    def test_translate_six_frames_ncbi_table_id(self):
        # rc = CAAUUU
        for seq in RNA('AAAUUG'), DNA('AAATTG'):
            # default
            obs = list(seq.translate_six_frames())
            self.assertEqual(obs, [Protein('KL'), Protein('N'), Protein('I'),
                                   Protein('QF'), Protein('N'), Protein('I')])

            obs = list(seq.translate_six_frames(9))
            self.assertEqual(obs, [Protein('NL'), Protein('N'), Protein('I'),
                                   Protein('QF'), Protein('N'), Protein('I')])

    def test_translate_six_frames_genetic_code_object(self):
        gc = GeneticCode('M' * 64, '-' * 64)
        for seq in RNA('AAAUUG'), DNA('AAATTG'):
            obs = list(seq.translate_six_frames(gc))
            self.assertEqual(obs, [Protein('MM'), Protein('M'), Protein('M'),
                                   Protein('MM'), Protein('M'), Protein('M')])

    def test_translate_six_frames_passes_parameters_through(self):
        for seq in RNA('UUUAUGUGGUGA'), DNA('TTTATGTGGTGA'):
            # mix of args and kwargs
            obs = next(seq.translate_six_frames(11, start='require',
                                                stop='require'))
            self.assertEqual(obs, Protein('MW'))

            # kwargs only
            obs = next(seq.translate_six_frames(genetic_code=11,
                                                start='require',
                                                stop='require'))
            self.assertEqual(obs, Protein('MW'))

            # args only
            obs = next(seq.translate_six_frames(11, 'require', 'require'))
            self.assertEqual(obs, Protein('MW'))

    def test_translate_six_frames_preserves_metadata(self):
        metadata = {'foo': 'bar', 'baz': 42}
        positional_metadata = {'foo': range(3)}
        for seq in (RNA('AUG', metadata=metadata,
                        positional_metadata=positional_metadata),
                    DNA('ATG', metadata=metadata,
                        positional_metadata=positional_metadata)):
            obs = list(seq.translate_six_frames())[:2]
            # metadata retained, positional metadata dropped
            self.assertEqual(
                obs,
                [Protein('M', metadata={'foo': 'bar', 'baz': 42}),
                 Protein('', metadata={'foo': 'bar', 'baz': 42})])

    def test_translate_six_frames_invalid_id(self):
        for seq in RNA('AUG'), DNA('ATG'):
            with self.assertRaisesRegex(ValueError, r'table_id.*42'):
                seq.translate_six_frames(42)

    def test_repr(self):
        # basic sanity checks for custom repr stats. more extensive testing is
        # performed on Sequence.__repr__

        for seq in DNA(''), RNA(''):
            obs = repr(seq)
            # obtained from super()
            self.assertIn('has gaps: False', obs)
            # custom to Protein
            self.assertIn('GC-content: 0.00%', obs)

        for seq in DNA('ACGT'), RNA('ACGU'):
            obs = repr(seq)
            self.assertIn('has gaps: False', obs)
            self.assertIn('GC-content: 50.00%', obs)

        for seq in DNA('CST'), RNA('CSU'):
            obs = repr(seq)
            self.assertIn('has gaps: False', obs)
            self.assertIn('GC-content: 66.67%', obs)

        for seq in DNA('GCSSCG'), RNA('GCSSCG'):
            obs = repr(seq)
            self.assertIn('has gaps: False', obs)
            self.assertIn('GC-content: 100.00%', obs)

        for seq in DNA('-GCSSCG.'), RNA('-GCSSCG.'):
            obs = repr(seq)
            self.assertIn('has gaps: True', obs)
            self.assertIn('GC-content: 100.00%', obs)

    def test_complement_without_reverse_empty(self):
        for constructor in (DNA, RNA):
            # without optional attributes
            comp = constructor('').complement()
            self.assertEqual(comp, constructor(''))

            # with optional attributes
            comp = constructor(
                '',
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': []},
                interval_metadata=IntervalMetadata(0)).complement()
            self.assertEqual(
                comp,
                constructor(
                    '',
                    metadata={'id': 'foo', 'description': 'bar'},
                    positional_metadata={'quality': []}))

    def test_complement_without_reverse_non_empty(self):
        for (constructor, seq_str, comp_str,
             qual) in self.all_combos_comp_qual:
            comp = constructor(seq_str).complement()
            self.assertEqual(comp, constructor(comp_str))

            im = IntervalMetadata(len(seq_str))
            im.add([(0, 1)], metadata={'gene': 'p53'})
            comp = constructor(
                seq_str,
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': qual},
                interval_metadata=im).complement()
            self.assertEqual(
                comp,
                constructor(
                    comp_str,
                    metadata={'id': 'foo', 'description': 'bar'},
                    positional_metadata={'quality': qual},
                    interval_metadata=im))

    def test_complement_with_reverse_empty(self):
        for constructor in (DNA, RNA):
            rc = constructor('').complement(reverse=True)
            self.assertEqual(rc, constructor(''))

            rc = constructor(
                '',
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': []},
                interval_metadata=IntervalMetadata(0)).complement(reverse=True)
            self.assertEqual(
                rc,
                constructor(
                    '',
                    metadata={'id': 'foo', 'description': 'bar'},
                    positional_metadata={'quality': []}))

    def test_complement_with_reverse_non_empty(self):
        for (constructor, seq_str, rev_comp_str,
             qual) in self.all_combos_rev_comp_qual:
            rc = constructor(seq_str).complement(reverse=True)
            self.assertEqual(rc, constructor(rev_comp_str))

            length = len(seq_str)
            im = IntervalMetadata(length)
            im.add([(0, 1)], metadata={'gene': 'p53'})
            im_rc = IntervalMetadata(length)
            im_rc.add([(length-1, length)], metadata={'gene': 'p53'})
            original = constructor(
                seq_str,
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={
                    'quality': qual},
                interval_metadata=im)
            rc = original.complement(reverse=True)

            self.assertEqual(
                rc,
                constructor(
                    rev_comp_str,
                    metadata={'id': 'foo', 'description': 'bar'},
                    positional_metadata={'quality':
                                         list(qual)[::-1]},
                    interval_metadata=im_rc))
            # assert the original object is not changed
            self.assertIsNot(original.interval_metadata, im)
            self.assertEqual(original.interval_metadata, im)

    def test_reverse_complement(self):
        # light tests because this just calls
        # NucleotideSequence.complement(reverse=True), which is tested more
        # extensively
        for (constructor, seq_str, rev_comp_str,
             qual) in self.all_combos_rev_comp_qual:
            rc = constructor(
                seq_str,
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': qual}).reverse_complement()
            self.assertEqual(
                rc,
                constructor(
                    rev_comp_str,
                    metadata={'id': 'foo', 'description': 'bar'},
                    positional_metadata={'quality': list(qual)[::-1]}))

    def test_is_reverse_complement_varied_types(self):
        tested = 0
        for constructor, seq_str, rev_comp_str in self.all_combos_rev_comp:
            seq_kinds = self.sequence_kinds.union(frozenset([constructor]))
            for sequence in seq_kinds:
                tested += 1
                seq1 = constructor(seq_str)
                seq2 = sequence(rev_comp_str)

                self.assertTrue(seq1.is_reverse_complement(seq2))

        self.assertEqual(tested, 8)

    def test_is_reverse_complement_empty(self):
        for constructor in (DNA, RNA):
            seq1 = constructor('')
            self.assertTrue(seq1.is_reverse_complement(seq1))

            # optional attributes are ignored, only the sequence is compared
            seq2 = constructor(
                '',
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality':
                                     np.array([], dtype=np.int64)})
            self.assertTrue(seq2.is_reverse_complement(seq2))
            self.assertTrue(seq1.is_reverse_complement(seq2))
            self.assertTrue(seq2.is_reverse_complement(seq1))

    def test_is_reverse_complement_metadata_ignored(self):
        for (constructor, seq_str, rev_comp_str,
             qual) in self.all_combos_rev_comp_qual:
            seq1 = constructor(seq_str)
            seq2 = constructor(
                rev_comp_str,
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': qual})

            self.assertFalse(seq1.is_reverse_complement(seq1))
            self.assertFalse(seq2.is_reverse_complement(seq2))

            self.assertTrue(seq1.is_reverse_complement(seq2))
            self.assertTrue(seq2.is_reverse_complement(seq1))

    def test_is_reverse_complement_non_reverse_complements(self):
        for constructor in (DNA, RNA):
            # same length
            seq1 = constructor('ACAG')
            seq2 = constructor('AAAA')

            self.assertFalse(seq1.is_reverse_complement(seq1))
            self.assertFalse(seq2.is_reverse_complement(seq2))

            self.assertFalse(seq1.is_reverse_complement(seq2))
            self.assertFalse(seq2.is_reverse_complement(seq1))

            # different length
            seq1 = constructor('ACAG')
            seq2 = constructor('AAAAA')

            self.assertFalse(seq1.is_reverse_complement(seq1))
            self.assertFalse(seq2.is_reverse_complement(seq2))

            self.assertFalse(seq1.is_reverse_complement(seq2))
            self.assertFalse(seq2.is_reverse_complement(seq1))

    def test_is_reverse_complement_type_mismatch(self):
        for Class in (DNA, RNA):
            class DifferentSequenceClass(GrammaredSequence):
                @classproperty
                def degenerate_map(cls):
                    return {"X": set("AB")}

                @classproperty
                def definite_chars(cls):
                    return set("ABC")

                @classproperty
                def default_gap_char(cls):
                    return '-'

                @classproperty
                def gap_chars(cls):
                    return set('-.')

            seq1 = Class('ABC')
            seq2 = DifferentSequenceClass('ABC')

            with self.assertRaisesRegex(TypeError,
                                        r"Cannot use.*and "
                                        "DifferentSequenceClass together"):
                seq1.is_reverse_complement(seq2)

    def test_motif_purine_run(self):
        dna = (DNA, "AARC--TCRG", "AA-RC--TCR-G")
        rna = (RNA, "AARC--UCRG", "AA-RC--UCR-G")
        all_sets = (dna, rna)

        for constructor, run1, run2 in all_sets:
            seq = constructor("")
            self.assertEqual(list(seq.find_motifs("purine-run")), [])

            seq = constructor(run1)
            self.assertEqual(list(seq.find_motifs("purine-run")),
                             [slice(0, 3), slice(8, 10)])

            seq = constructor(run2)
            self.assertEqual(list(seq.find_motifs("purine-run", min_length=3,
                                                  ignore=seq.gaps())),
                             [slice(0, 4)])

    def test_motif_pyrimidine_run(self):
        dna = (DNA, "AARC--TCRA", "AA-RC--TCR-A")
        rna = (RNA, "AARC--UCRG", "AA-RC--UCR-G")
        all_sets = (dna, rna)

        for constructor, run1, run2 in all_sets:
            seq = constructor("")
            self.assertEqual(list(seq.find_motifs("pyrimidine-run")), [])

            seq = constructor(run1)
            self.assertEqual(list(seq.find_motifs("pyrimidine-run")),
                             [slice(3, 4), slice(6, 8)])

            seq = constructor(run2)
            self.assertEqual(list(seq.find_motifs("pyrimidine-run",
                                                  min_length=3,
                                                  ignore=seq.gaps())),
                             [slice(4, 9)])

    def test_gc_frequency_and_gc_content(self):
        universal_sets = (('', 0, 0.0), ('ADDDH', 0, 0.0), ('ACGA', 2, 0.5),
                          ('ACGS', 3, 0.75), ('AAAAAAAG', 1, 0.125),
                          ('CCC', 3, 1.0), ('GGG', 3, 1.0), ('SSS', 3, 1.0),
                          ('CGS', 3, 1.0), ('----....', 0, 0.0),
                          ('G--..', 1, 1.0), ('ACGA', 2, 0.5))
        dna = (DNA, universal_sets + (('ATMRWYKVHDBN.-', 0, 0.0),))
        rna = (RNA, universal_sets + (('AUMRWYKVHDBN.-', 0, 0.0),))
        for constructor, current_set in (dna, rna):
            for seq_str, count, ratio in current_set:
                seq = constructor(seq_str)
                self.assertEqual(count, seq.gc_frequency())
                self.assertEqual(count, seq.gc_frequency(relative=False))
                self.assertEqual(ratio, seq.gc_frequency(relative=True))
                self.assertEqual(ratio, seq.gc_content())


if __name__ == "__main__":
    unittest.main()
