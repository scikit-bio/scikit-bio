# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from six import StringIO

from unittest import TestCase, main

from skbio.io import PhylipFormatError
from skbio.io.phylip import _alignment_to_phylip
from skbio import Alignment, DNA, RNA
from skbio.util import get_data_path


class AlignmentWriterTests(TestCase):
    def setUp(self):
        # ids all same length, seqs longer than 10 chars
        dna_3_seqs = Alignment([
            DNA('..ACC-GTTGG..', id="d1"),
            DNA('TTACCGGT-GGCC', id="d2"),
            DNA('.-ACC-GTTGC--', id="d3")])

        # id lengths from 0 to 10, with mixes of numbers, characters, and
        # spaces. sequence characters are a mix of cases and gap characters.
        # sequences are shorter than 10 chars
        variable_length_ids = Alignment([
            RNA('.-ACGU'),
            RNA('UGCA-.', id='a'),
            RNA('.ACGU-', id='bb'),
            RNA('ugca-.', id='1'),
            RNA('AaAaAa', id='abcdefghij'),
            RNA('GGGGGG', id='ab def42ij')])

        # sequences with 20 chars = exactly two chunks of size 10
        two_chunks = Alignment([
            DNA('..ACC-GTTGG..AATGC.C', id='foo'),
            DNA('TTACCGGT-GGCCTA-GCAT', id='bar')])

        # single sequence with more than two chunks
        single_seq_long = Alignment([
            DNA('..ACC-GTTGG..AATGC.C----', id='foo')])

        # single sequence with only a single character (minimal writeable
        # alignment)
        single_seq_short = Alignment([DNA('-')])

        # alignments that can be written in phylip format
        self.objs = [dna_3_seqs, variable_length_ids, two_chunks,
                     single_seq_long, single_seq_short]
        self.fps = map(get_data_path,
                       ['phylip_dna_3_seqs', 'phylip_variable_length_ids',
                        'phylip_two_chunks', 'phylip_single_seq_long',
                        'phylip_single_seq_short'])

        # alignments that cannot be written in phylip format, paired with their
        # expected error message regexps
        self.invalid_objs = [
            # unequal length
            (Alignment([DNA('A-CT', id="d1"),
                        DNA('TTA', id="d2"),
                        DNA('.-AC', id="d3")]), 'equal length'),

            # invalid chars
            (Alignment([DNA('ACG', id="d1"),
                        DNA('FOO', id="d2")]), 'valid characters'),

            # no seqs
            (Alignment([]), 'one sequence'),

            # no positions
            (Alignment([DNA('', id="d1"),
                        DNA('', id="d2")]), 'one position'),

            # ids too long
            (Alignment([RNA('ACGU', id="foo"),
                        RNA('UGCA', id="alongsequenceid")]),
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
            # thrown. TODO remove this check when #674 is resolved
            obs = fh.getvalue()
            fh.close()
            self.assertEqual(obs, '')


if __name__ == '__main__':
    main()
