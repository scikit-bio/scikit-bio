# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

import numpy as np

from skbio.sequence import DNA, RNA
from skbio.sequence._nucleotide_mixin import NucleotideMixin


# This file contains tests for functionality of sequence types which implement
# NucleotideMixin. Currently this means DNA and RNA. These types are so
# similar that the testing logic can be shared and parameterized across
# different test data.

class TestNucelotideSequence(unittest.TestCase):
    def setUp(self):
        self.sequence_kinds = frozenset([
            str,
            lambda s: np.fromstring(s, dtype='|S1'),
            lambda s: np.fromstring(s, dtype=np.uint8)])

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

    def test_nondegenerate_chars(self):
        dna = (DNA, "ACGT")
        rna = (RNA, "ACGU")
        for constructor, nondegenerate in (dna, rna):
            exp = set(nondegenerate)
            self.assertEqual(constructor('').nondegenerate_chars, exp)
            self.assertEqual(constructor.nondegenerate_chars, exp)

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

    def test_complement_without_reverse_empty(self):
        for constructor in (DNA, RNA):
            # without optional attributes
            comp = constructor('').complement()
            self.assertEqual(comp, constructor(''))

            # with optional attributes
            comp = constructor(
                '',
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': []}).complement()
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

            comp = constructor(
                seq_str,
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': qual}).complement()
            self.assertEqual(
                comp,
                constructor(
                    comp_str,
                    metadata={'id': 'foo', 'description': 'bar'},
                    positional_metadata={'quality': qual}))

    def test_complement_with_reverse_empty(self):
        for constructor in (DNA, RNA):
            rc = constructor('').complement(reverse=True)
            self.assertEqual(rc, constructor(''))

            rc = constructor(
                '',
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': []}).complement(reverse=True)
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

            rc = constructor(
                seq_str,
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={
                    'quality': qual}).complement(reverse=True)
            self.assertEqual(
                rc,
                constructor(
                    rev_comp_str,
                    metadata={'id': 'foo', 'description': 'bar'},
                    positional_metadata={'quality':
                                         list(qual)[::-1]}))

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
            class Subclass(Class):
                pass
            seq1 = Class('ABC')
            seq2 = Subclass('ABC')

            with self.assertRaises(TypeError):
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

if __name__ == "__main__":
    unittest.main()
