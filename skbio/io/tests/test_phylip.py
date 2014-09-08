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

from skbio.io.phylip import _alignment_to_phylip
from skbio import Alignment, DNASequence
from skbio.alignment import SequenceCollectionError
from skbio.util import get_data_path


class AlignmentWriterTests(TestCase):
    def setUp(self):
        d1 = DNASequence('..ACC-GTTGG..', id="d1")
        d2 = DNASequence('TTACCGGT-GGCC', id="d2")
        d3 = DNASequence('.-ACC-GTTGC--', id="d3")
        dna_3_seqs = Alignment([d1, d2, d3])

        self.objs = [dna_3_seqs]
        self.fps = map(get_data_path, ['phylip_dna_3_seqs'])

    def test_write(self):
        for fp, obj in zip(self.fps, self.objs):
            fh = StringIO()
            _alignment_to_phylip(obj, fh)
            obs = fh.getvalue()
            fh.close()

            with open(fp, 'U') as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

    def test_write_unequal_sequence_lengths(self):
        d1 = DNASequence('A-CT', id="d1")
        d2 = DNASequence('TTA', id="d2")
        d3 = DNASequence('.-AC', id="d3")
        obj = Alignment([d1, d2, d3])

        with self.assertRaisesRegexp(SequenceCollectionError, 'equal length'):
            fh = StringIO()
            _alignment_to_phylip(obj, fh)

        # ensure nothing was written to the file before erroring
        obs = fh.getvalue()
        fh.close()
        self.assertEqual(obs, '')

    def test_write_no_sequences(self):
        with self.assertRaisesRegexp(SequenceCollectionError, 'one sequence'):
            fh = StringIO()
            _alignment_to_phylip(Alignment([]), fh)

        obs = fh.getvalue()
        fh.close()
        self.assertEqual(obs, '')

    def test_write_no_positions(self):
        d1 = DNASequence('', id="d1")
        d2 = DNASequence('', id="d2")
        obj = Alignment([d1, d2])

        with self.assertRaisesRegexp(SequenceCollectionError, 'one position'):
            fh = StringIO()
            _alignment_to_phylip(obj, fh)

        obs = fh.getvalue()
        fh.close()
        self.assertEqual(obs, '')


if __name__ == '__main__':
    main()
