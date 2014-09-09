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

        # alignments that can be written in phylip format
        self.objs = [dna_3_seqs]
        self.fps = map(get_data_path, ['phylip_dna_3_seqs'])

        # alignments that cannot be written in phylip format
        self.invalid_objs = [
            # unequal length
            (Alignment([DNASequence('A-CT', id="d1"),
                        DNASequence('TTA', id="d2"),
                        DNASequence('.-AC', id="d3")]), 'equal length'),

            # invalid chars
            (Alignment([DNASequence('ACG', id="d1"),
                        DNASequence('FOO', id="d2")]), 'valid characters'),

            # no seqs
            (Alignment([]), 'one sequence'),

            # no positions
            (Alignment([DNASequence('', id="d1"),
                        DNASequence('', id="d2")]), 'one position')
        ]

    def test_write(self):
        for fp, obj in zip(self.fps, self.objs):
            fh = StringIO()
            _alignment_to_phylip(obj, fh)
            obs = fh.getvalue()
            fh.close()

            with open(fp, 'U') as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

    def test_write_invalid_alignment(self):
        for invalid_obj, error_msg_regexp in self.invalid_objs:
            fh = StringIO()
            with self.assertRaisesRegexp(SequenceCollectionError,
                                         error_msg_regexp):
                _alignment_to_phylip(invalid_obj, fh)

            # ensure nothing was written to the file before the error was
            # thrown
            obs = fh.getvalue()
            fh.close()
            self.assertEqual(obs, '')


if __name__ == '__main__':
    main()
