# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
import unittest

import numpy as np
import numpy.testing as npt

from skbio import Sequence, DNA, RNA, Protein, GeneticCode
from skbio.sequence._genetic_code import _ncbi_genetic_codes


class TestGeneticCode(unittest.TestCase):
    def setUp(self):
        self.sgc = GeneticCode.from_ncbi(1)

    def test_from_ncbi_valid_table_ids(self):
        # spot check a few tables
        self.assertEqual(GeneticCode.from_ncbi().name,
                         'Standard')
        self.assertEqual(GeneticCode.from_ncbi(2).name,
                         'Vertebrate Mitochondrial')
        self.assertEqual(GeneticCode.from_ncbi(12).name,
                         'Alternative Yeast Nuclear')
        self.assertEqual(GeneticCode.from_ncbi(25).name,
                         'Candidate Division SR1 and Gracilibacteria')

    def test_from_ncbi_invalid_input(self):
        with self.assertRaisesRegex(ValueError, r'table_id.*7'):
            GeneticCode.from_ncbi(7)
        with self.assertRaisesRegex(ValueError, r'table_id.*42'):
            GeneticCode.from_ncbi(42)

    def test_reading_frames(self):
        exp = [1, 2, 3, -1, -2, -3]
        self.assertEqual(GeneticCode.reading_frames, exp)
        self.assertEqual(self.sgc.reading_frames, exp)

        GeneticCode.reading_frames.append(42)

        self.assertEqual(GeneticCode.reading_frames, exp)
        self.assertEqual(self.sgc.reading_frames, exp)

        with self.assertRaises(AttributeError):
            self.sgc.reading_frames = [1, 2, 42]

    def test_name(self):
        self.assertEqual(self.sgc.name, 'Standard')
        self.assertEqual(GeneticCode('M' * 64, '-' * 64).name, '')
        self.assertEqual(GeneticCode('M' * 64, '-' * 64, 'foo').name, 'foo')

        with self.assertRaises(AttributeError):
            self.sgc.name = 'foo'

    def test_init_varied_equivalent_input(self):
        for args in (('M' * 64, '-' * 64),
                     (Protein('M' * 64), Protein('-' * 64)),
                     (Sequence('M' * 64), Sequence('-' * 64))):
            gc = GeneticCode(*args)
            self.assertEqual(gc.name, '')
            self.assertEqual(gc._amino_acids, Protein('M' * 64))
            self.assertEqual(gc._starts, Protein('-' * 64))
            npt.assert_array_equal(gc._m_character_codon,
                                   np.asarray([0, 0, 0], dtype=np.uint8))
            self.assertEqual(len(gc._start_codons), 0)

    def test_init_invalid_input(self):
        # `amino_acids` invalid protein
        with self.assertRaisesRegex(ValueError, r'Invalid character.*&'):
            GeneticCode('&' * 64, '-' * 64)

        # wrong number of amino acids
        with self.assertRaisesRegex(ValueError, r'amino_acids.*64.*42'):
            GeneticCode('M' * 42, '-' * 64)

        # `amino_acids` missing M
        with self.assertRaisesRegex(ValueError, r'amino_acids.*M.*character'):
            GeneticCode('A' * 64, '-' * 64)

        # `starts` invalid protein
        with self.assertRaisesRegex(ValueError, r'Invalid character.*&'):
            GeneticCode('M' * 64, '&' * 64)

        # wrong number of starts
        with self.assertRaisesRegex(ValueError, r'starts.*64.*42'):
            GeneticCode('M' * 64, '-' * 42)

        # invalid characters in `starts`
        with self.assertRaisesRegex(ValueError, r'starts.*M and - characters'):
            GeneticCode('M' * 64, '-M' * 30 + '*AQR')

    def test_str(self):
        # predefined
        exp = (
            '  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAA'
            'DDEEGGGG\n'
            'Starts = ---M---------------M---------------M--------------------'
            '--------\n'
            'Base1  = UUUUUUUUUUUUUUUUCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGG'
            'GGGGGGGG\n'
            'Base2  = UUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCC'
            'AAAAGGGG\n'
            'Base3  = UCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAG'
            'UCAGUCAG'
        )
        self.assertEqual(str(self.sgc), exp)

        # custom, no name
        obs = str(GeneticCode('M' * 64, '-' * 64))
        self.assertIn('M' * 64, obs)
        self.assertIn('-' * 64, obs)

    def test_repr(self):
        # predefined
        exp = (
            'GeneticCode (Standard)\n'
            '-----------------------------------------------------------------'
            '--------\n'
            '  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAA'
            'DDEEGGGG\n'
            'Starts = ---M---------------M---------------M--------------------'
            '--------\n'
            'Base1  = UUUUUUUUUUUUUUUUCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGG'
            'GGGGGGGG\n'
            'Base2  = UUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCC'
            'AAAAGGGG\n'
            'Base3  = UCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAG'
            'UCAGUCAG'
        )
        self.assertEqual(repr(self.sgc), exp)

        # custom, no name
        obs = repr(GeneticCode('M' * 64, '-' * 64))
        self.assertTrue(obs.startswith('GeneticCode\n'))
        self.assertIn('M' * 64, obs)
        self.assertIn('-' * 64, obs)

    def test_eq(self):
        amino_acids = 'AMPM' * 16
        starts = '--M-' * 16

        equal_gcs = [
            GeneticCode(amino_acids, starts),
            # name should be ignored
            GeneticCode(amino_acids, starts, 'foo'),
            # metadata/positional metadata should be ignored if Sequence
            # subclass is provided
            GeneticCode(
                Protein(amino_acids, metadata={'foo': 'bar'}),
                Protein(starts, positional_metadata={'foo': range(64)}))
        ]

        # every gc should be equal to itself
        for gc in equal_gcs:
            self.assertTrue(gc == gc)
            self.assertFalse(gc != gc)

        # every pair of gcs should be equal. use permutations instead of
        # combinations to test that comparing gc1 to gc2 and gc2 to gc1 are
        # both equal
        for gc1, gc2 in itertools.permutations(equal_gcs, 2):
            self.assertTrue(gc1 == gc2)
            self.assertFalse(gc1 != gc2)

    def test_ne(self):
        class GeneticCodeSubclass(GeneticCode):
            pass

        amino_acids = 'AMPM' * 16
        starts = '--M-' * 16

        unequal_gcs = [
            GeneticCode(amino_acids, starts),
            # type must match
            GeneticCodeSubclass(amino_acids, starts),
            # completely different type
            'foo'
        ]
        # none of the NCBI genetic codes should be equal to each other
        unequal_gcs.extend(_ncbi_genetic_codes.values())

        for gc in unequal_gcs:
            self.assertTrue(gc == gc)
            self.assertFalse(gc != gc)

        for gc1, gc2 in itertools.permutations(unequal_gcs, 2):
            self.assertTrue(gc1 != gc2)
            self.assertFalse(gc1 == gc2)

    def test_translate_preserves_metadata(self):
        obs = self.sgc.translate(
            RNA('AUG', metadata={'foo': 'bar', 'baz': 42},
                positional_metadata={'foo': range(3)}))
        # metadata retained, positional metadata dropped
        self.assertEqual(obs, Protein('M',
                                      metadata={'foo': 'bar', 'baz': 42}))

    def test_translate_default_behavior(self):
        # empty translation
        exp = Protein('')
        for seq in RNA(''), RNA('A'), RNA('AU'):
            obs = self.sgc.translate(seq)
            self.assertEqual(obs, exp)

        # no start or stop codons
        obs = self.sgc.translate(RNA('CCU'))
        self.assertEqual(obs, Protein('P'))

        # multiple alternative start codons, no stop codons, length is multiple
        # of 3
        obs = self.sgc.translate(RNA('CAUUUGCUGAAA'))
        self.assertEqual(obs, Protein('HLLK'))

        # multiple stop codons, length isn't multiple of 3
        obs = self.sgc.translate(RNA('UUUUUUUAAAGUUAAGGGAU'))
        self.assertEqual(obs, Protein('FF*S*G'))

    def test_translate_reading_frame_empty_translation(self):
        exp = Protein('')
        for seq in RNA(''), RNA('A'), RNA('AU'):
            for reading_frame in GeneticCode.reading_frames:
                obs = self.sgc.translate(seq, reading_frame=reading_frame)
                self.assertEqual(obs, exp)

        # reading frames that yield a partial codon
        for reading_frame in 2, 3, -2, -3:
            obs = self.sgc.translate(RNA('AUG'), reading_frame=reading_frame)
            self.assertEqual(obs, exp)

    def test_translate_reading_frame_non_empty_translation(self):
        seq = RNA('AUGGUGGAA')  # rc = UUCCACCAU
        for reading_frame, exp_str in ((1, 'MVE'), (2, 'WW'), (3, 'GG'),
                                       (-1, 'FHH'), (-2, 'ST'), (-3, 'PP')):
            exp = Protein(exp_str)
            obs = self.sgc.translate(seq, reading_frame=reading_frame)
            self.assertEqual(obs, exp)

    def test_translate_start_empty_translation(self):
        exp = Protein('')
        for seq in RNA(''), RNA('A'), RNA('AU'):
            for start in {'optional', 'ignore'}:
                obs = self.sgc.translate(seq, start=start)
                self.assertEqual(obs, exp)

            with self.assertRaisesRegex(ValueError,
                                        r'reading_frame=1.*start=\'require\''):
                self.sgc.translate(seq, start='require')

    def test_translate_start_with_start_codon(self):
        # trim before start codon, replace with M. ensure alternative start
        # codons following the start codon aren't replaced with M. ensure
        # default behavior for handling stop codons is retained
        seq = RNA('CAUUUGCUGAAAUGA')
        exp = Protein('MLK*')
        for start in {'require', 'optional'}:
            obs = self.sgc.translate(seq, start=start)
            self.assertEqual(obs, exp)

        # ignore start codon replacement and trimming; just translate
        exp = Protein('HLLK*')
        obs = self.sgc.translate(seq, start='ignore')
        self.assertEqual(obs, exp)

        # just a start codon, no replacement necessary
        seq = RNA('AUG')
        exp = Protein('M')
        for start in {'require', 'optional', 'ignore'}:
            obs = self.sgc.translate(seq, start=start)
            self.assertEqual(obs, exp)

        # single alternative start codon
        seq = RNA('CUG')
        exp = Protein('M')
        for start in {'require', 'optional'}:
            obs = self.sgc.translate(seq, start=start)
            self.assertEqual(obs, exp)

        exp = Protein('L')
        obs = self.sgc.translate(seq, start='ignore')
        self.assertEqual(obs, exp)

    def test_translate_start_no_start_codon(self):
        seq = RNA('CAACAACAGCAA')
        exp = Protein('QQQQ')
        for start in {'ignore', 'optional'}:
            obs = self.sgc.translate(seq, start=start)
            self.assertEqual(obs, exp)

        with self.assertRaisesRegex(ValueError,
                                    r'reading_frame=1.*start=\'require\''):
            self.sgc.translate(seq, start='require')

        # non-start codon that translates to an AA that start codons also map
        # to. should catch bug if code attempts to search and trim *after*
        # translation -- this must happen *before* translation
        seq = RNA('UUACAA')
        exp = Protein('LQ')
        for start in {'ignore', 'optional'}:
            obs = self.sgc.translate(seq, start=start)
            self.assertEqual(obs, exp)

        with self.assertRaisesRegex(ValueError,
                                    r'reading_frame=1.*start=\'require\''):
            self.sgc.translate(seq, start='require')

    def test_translate_start_no_accidental_mutation(self):
        # `start` mutates a vector in-place that is derived from
        # GeneticCode._offset_table. the current code doesn't perform an
        # explicit copy because numpy's advanced indexing is used, which always
        # returns a copy. test this assumption here in case that behavior
        # changes in the future
        offset_table = self.sgc._offset_table.copy()

        seq = RNA('CAUUUGCUGAAAUGA')
        obs = self.sgc.translate(seq, start='require')
        self.assertEqual(obs, Protein('MLK*'))

        npt.assert_array_equal(self.sgc._offset_table, offset_table)

    def test_translate_stop_empty_translation(self):
        exp = Protein('')
        for seq in RNA(''), RNA('A'), RNA('AU'):
            for stop in {'optional', 'ignore'}:
                obs = self.sgc.translate(seq, stop=stop)
                self.assertEqual(obs, exp)

            with self.assertRaisesRegex(ValueError,
                                        r'reading_frame=1.*stop=\'require\''):
                self.sgc.translate(seq, stop='require')

    def test_translate_stop_with_stop_codon(self):
        # multiple stop codons with trailing codons
        seq = RNA('UGGACUUGAUAUCGUUAGGAU')
        exp = Protein('WT')
        for stop in {'require', 'optional'}:
            obs = self.sgc.translate(seq, stop=stop)
            self.assertEqual(obs, exp)

        # ignore stop codon trimming; just translate
        exp = Protein('WT*YR*D')
        obs = self.sgc.translate(seq, stop='ignore')
        self.assertEqual(obs, exp)

        # ends with single stop codon
        seq = RNA('UGUCUGUAA')
        exp = Protein('CL')
        for stop in {'require', 'optional'}:
            obs = self.sgc.translate(seq, stop=stop)
            self.assertEqual(obs, exp)

        exp = Protein('CL*')
        obs = self.sgc.translate(seq, stop='ignore')
        self.assertEqual(obs, exp)

        # just a stop codon
        seq = RNA('UAG')
        exp = Protein('')
        for stop in {'require', 'optional'}:
            obs = self.sgc.translate(seq, stop=stop)
            self.assertEqual(obs, exp)

        exp = Protein('*')
        obs = self.sgc.translate(seq, stop='ignore')
        self.assertEqual(obs, exp)

    def test_translate_stop_no_stop_codon(self):
        seq = RNA('GAAUCU')
        exp = Protein('ES')
        for stop in {'ignore', 'optional'}:
            obs = self.sgc.translate(seq, stop=stop)
            self.assertEqual(obs, exp)

        with self.assertRaisesRegex(ValueError,
                                    r'reading_frame=1.*stop=\'require\''):
            self.sgc.translate(seq, stop='require')

    def test_translate_trim_to_cds(self):
        seq = RNA('UAAUUGCCUCAUUAAUAACAAUGA')

        # find first start codon, trim all before it, convert alternative start
        # codon to M, finally trim to first stop codon following the start
        # codon
        exp = Protein('MPH')
        for param in {'require', 'optional'}:
            obs = self.sgc.translate(seq, start=param, stop=param)
            self.assertEqual(obs, exp)

        exp = Protein('*LPH**Q*')
        obs = self.sgc.translate(seq, start='ignore', stop='ignore')
        self.assertEqual(obs, exp)

        # alternative reading frame disrupts cds:
        #     AAUUGCCUCAUUAAUAACAAUGA
        #     NCLINNN
        with self.assertRaisesRegex(ValueError,
                                    r'reading_frame=2.*start=\'require\''):
            self.sgc.translate(seq, reading_frame=2, start='require')
        with self.assertRaisesRegex(ValueError,
                                    r'reading_frame=2.*stop=\'require\''):
            self.sgc.translate(seq, reading_frame=2, stop='require')

        exp = Protein('NCLINNN')
        for param in {'ignore', 'optional'}:
            obs = self.sgc.translate(seq, reading_frame=2, start=param,
                                     stop=param)
            self.assertEqual(obs, exp)

    def test_translate_invalid_input(self):
        # invalid sequence type
        with self.assertRaisesRegex(TypeError, r'RNA.*DNA'):
            self.sgc.translate(DNA('ACG'))
        with self.assertRaisesRegex(TypeError, r'RNA.*str'):
            self.sgc.translate('ACG')

        # invalid reading frame
        with self.assertRaisesRegex(ValueError, r'\[1, 2, 3, -1, -2, -3\].*0'):
            self.sgc.translate(RNA('AUG'), reading_frame=0)

        # invalid start
        with self.assertRaisesRegex(ValueError, r'start.*foo'):
            self.sgc.translate(RNA('AUG'), start='foo')

        # invalid stop
        with self.assertRaisesRegex(ValueError, r'stop.*foo'):
            self.sgc.translate(RNA('AUG'), stop='foo')

        # gapped sequence
        with self.assertRaisesRegex(ValueError, r'gapped'):
            self.sgc.translate(RNA('UU-G'))

        # degenerate sequence
        with self.assertRaisesRegex(NotImplementedError, r'degenerate'):
            self.sgc.translate(RNA('RUG'))

    def test_translate_varied_genetic_codes(self):
        # spot check using a few NCBI and custom genetic codes to translate
        seq = RNA('AAUGAUGUGACUAUCAGAAGG')

        # table_id=2
        exp = Protein('NDVTI**')
        obs = GeneticCode.from_ncbi(2).translate(seq)
        self.assertEqual(obs, exp)

        exp = Protein('MTI')
        obs = GeneticCode.from_ncbi(2).translate(seq, start='require',
                                                 stop='require')
        self.assertEqual(obs, exp)

        # table_id=22
        exp = Protein('NDVTIRR')
        obs = GeneticCode.from_ncbi(22).translate(seq)
        self.assertEqual(obs, exp)

        with self.assertRaisesRegex(ValueError,
                                    r'reading_frame=1.*start=\'require\''):
            GeneticCode.from_ncbi(22).translate(seq, start='require',
                                                stop='require')

        # custom, no start codons
        gc = GeneticCode('MWN*' * 16, '-' * 64)
        exp = Protein('MM*MWN*')
        obs = gc.translate(seq)
        self.assertEqual(obs, exp)

        with self.assertRaisesRegex(ValueError,
                                    r'reading_frame=1.*start=\'require\''):
            gc.translate(seq, start='require', stop='require')

    def test_translate_six_frames(self):
        seq = RNA('AUGCUAACAUAAA')  # rc = UUUAUGUUAGCAU

        # test default behavior
        exp = [Protein('MLT*'), Protein('C*HK'), Protein('ANI'),
               Protein('FMLA'), Protein('LC*H'), Protein('YVS')]
        obs = list(self.sgc.translate_six_frames(seq))
        self.assertEqual(obs, exp)

        # test that start/stop are respected
        exp = [Protein('MLT'), Protein('C'), Protein('ANI'),
               Protein('MLA'), Protein('LC'), Protein('YVS')]
        obs = list(self.sgc.translate_six_frames(seq, start='optional',
                                                 stop='optional'))
        self.assertEqual(obs, exp)

    def test_translate_six_frames_preserves_metadata(self):
        seq = RNA('AUG', metadata={'foo': 'bar', 'baz': 42},
                  positional_metadata={'foo': range(3)})
        obs = list(self.sgc.translate_six_frames(seq))[:2]
        # metadata retained, positional metadata dropped
        self.assertEqual(
            obs,
            [Protein('M', metadata={'foo': 'bar', 'baz': 42}),
             Protein('', metadata={'foo': 'bar', 'baz': 42})])


if __name__ == '__main__':
    unittest.main()
