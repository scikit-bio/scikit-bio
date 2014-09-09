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

from skbio.io import PhylipFormatError
from skbio.io.phylip import _alignment_to_phylip
from skbio import Alignment, DNASequence, RNASequence
from skbio.util import get_data_path


class AlignmentWriterTests(TestCase):
    def setUp(self):
        # ids all same length, seqs longer than 10 chars
        dna_3_seqs = Alignment([
            DNASequence('..ACC-GTTGG..', id="d1"),
            DNASequence('TTACCGGT-GGCC', id="d2"),
            DNASequence('.-ACC-GTTGC--', id="d3")])

        # id lengths from 0 to 10, with mixes of numbers, characters, and
        # spaces. sequence characters are a mix of cases and gap characters
        variable_length_ids = Alignment([
            RNASequence('.-ACGU'),
            RNASequence('UGCA-.', id='a'),
            RNASequence('.ACGU-', id='bb'),
            RNASequence('ugca-.', id='1'),
            RNASequence('AaAaAa', id='abcdefghij'),
            RNASequence('GGGGGG', id='ab def42ij')])

        # alignments that can be written in phylip format
        self.objs = [dna_3_seqs, variable_length_ids]
        self.fps = map(get_data_path,
                       ['phylip_dna_3_seqs', 'phylip_variable_length_ids'])

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
                        DNASequence('', id="d2")]), 'one position'),

            # ids too long
            (Alignment([RNASequence('ACGU', id="foo"),
                        RNASequence('UGCA', id="alongsequenceid")]),
             '10.*alongsequenceid')
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
            with self.assertRaisesRegexp(PhylipFormatError, error_msg_regexp):
                _alignment_to_phylip(invalid_obj, fh)

            # ensure nothing was written to the file before the error was
            # thrown
            obs = fh.getvalue()
            fh.close()
            self.assertEqual(obs, '')


if __name__ == '__main__':
    main()
