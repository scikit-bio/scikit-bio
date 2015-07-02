# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import six

import unittest

from skbio import RNA, Protein, GeneticCode


class TestRNA(unittest.TestCase):
    def test_translate_ncbi_table_id(self):
        seq = RNA('AAAUUUAUGCAU')

        # default
        obs = seq.translate()
        self.assertEqual(obs, Protein('KFMH'))

        obs = seq.translate(9)
        self.assertEqual(obs, Protein('NFMH'))

    def test_translate_genetic_code_object(self):
        gc = GeneticCode('M' * 64, '-' * 64)
        obs = RNA('AAAUUUAUGCAU').translate(gc)
        self.assertEqual(obs, Protein('MMMM'))

    def test_translate_passes_parameters_through(self):
        seq = RNA('UAAAUUGUGGUAA')
        obs = seq.translate(13, reading_frame=2, start='require',
                            stop='require')
        self.assertEqual(obs, Protein('MW'))

    def test_translate_invalid_id(self):
        with six.assertRaisesRegex(self, ValueError, 'table_id.*42'):
            RNA('AUG').translate(42)

    def test_translate_six_frames_ncbi_table_id(self):
        seq = RNA('AAAUUG')  # rc = CAAUUU

        # default
        obs = list(seq.translate_six_frames())
        self.assertEqual(obs, [Protein('KL'), Protein('N'), Protein('I'),
                               Protein('QF'), Protein('N'), Protein('I')])

        obs = list(seq.translate_six_frames(9))
        self.assertEqual(obs, [Protein('NL'), Protein('N'), Protein('I'),
                               Protein('QF'), Protein('N'), Protein('I')])

    def test_translate_six_frames_genetic_code_object(self):
        gc = GeneticCode('M' * 64, '-' * 64)
        obs = list(RNA('AAAUUG').translate_six_frames(gc))
        self.assertEqual(obs, [Protein('MM'), Protein('M'), Protein('M'),
                               Protein('MM'), Protein('M'), Protein('M')])

    def test_translate_six_frames_passes_parameters_through(self):
        seq = RNA('UUUAUGUGGUGA')
        obs = next(seq.translate_six_frames(11, start='require',
                                            stop='require'))
        self.assertEqual(obs, Protein('MW'))

    def test_translate_six_frames_invalid_id(self):
        with six.assertRaisesRegex(self, ValueError, 'table_id.*42'):
            RNA('AUG').translate_six_frames(42)


if __name__ == '__main__':
    unittest.main()
