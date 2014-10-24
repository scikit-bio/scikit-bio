# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import range

import unittest

import numpy.testing as npt

from skbio.io._base import (_chunk_str, _decode_qual_to_phred,
                            _encode_phred_to_qual)


class ChunkStrTests(unittest.TestCase):
    def test_even_split(self):
        self.assertEqual(_chunk_str('abcdef', 6, ' '), 'abcdef')
        self.assertEqual(_chunk_str('abcdef', 3, ' '), 'abc def')
        self.assertEqual(_chunk_str('abcdef', 2, ' '), 'ab cd ef')
        self.assertEqual(_chunk_str('abcdef', 1, ' '), 'a b c d e f')
        self.assertEqual(_chunk_str('a', 1, ' '), 'a')
        self.assertEqual(_chunk_str('abcdef', 2, ''), 'abcdef')

    def test_no_split(self):
        self.assertEqual(_chunk_str('', 2, '\n'), '')
        self.assertEqual(_chunk_str('a', 100, '\n'), 'a')
        self.assertEqual(_chunk_str('abcdef', 42, '|'), 'abcdef')

    def test_uneven_split(self):
        self.assertEqual(_chunk_str('abcdef', 5, '|'), 'abcde|f')
        self.assertEqual(_chunk_str('abcdef', 4, '|'), 'abcd|ef')
        self.assertEqual(_chunk_str('abcdefg', 3, ' - '), 'abc - def - g')

    def test_invalid_n(self):
        with self.assertRaisesRegexp(ValueError, 'n=0'):
            _chunk_str('abcdef', 0, ' ')

        with self.assertRaisesRegexp(ValueError, 'n=-42'):
            _chunk_str('abcdef', -42, ' ')


class PhredDecoderTests(unittest.TestCase):
    def test_missing_variant_and_phred_offset(self):
        with self.assertRaises(ValueError) as cm:
            _decode_qual_to_phred('abcd')
        self.assertIn('`variant`', str(cm.exception))
        self.assertIn('`phred_offset`', str(cm.exception))
        self.assertIn('decode', str(cm.exception))

    def test_variant_and_phred_offset_provided(self):
        with self.assertRaises(ValueError) as cm:
            _decode_qual_to_phred('abcd', variant='sanger', phred_offset=64)
        self.assertIn('both', str(cm.exception))
        self.assertIn('`variant`', str(cm.exception))
        self.assertIn('`phred_offset`', str(cm.exception))

    def test_solexa_variant(self):
        with self.assertRaises(NotImplementedError) as cm:
            _decode_qual_to_phred('abcd', variant='solexa')
        self.assertIn('719', str(cm.exception))

    def test_unrecognized_variant(self):
        with self.assertRaises(ValueError) as cm:
            _decode_qual_to_phred('abcd', variant='illumina')
        self.assertIn('variant', str(cm.exception))
        self.assertIn("'illumina'", str(cm.exception))

    def test_empty_qual_str(self):
        self.assertEqual(_decode_qual_to_phred('', variant='sanger'), [])

    def test_sanger_variant(self):
        # test entire range of possible ascii chars for sanger
        all_sanger_ascii = ('!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP'
                            'QRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~')
        obs = _decode_qual_to_phred(all_sanger_ascii, variant='sanger')
        self.assertEqual(obs, list(range(94)))

        with self.assertRaises(ValueError) as cm:
            _decode_qual_to_phred('a b', variant='sanger')
        self.assertIn('-1', str(cm.exception))
        self.assertIn('[0, 93]', str(cm.exception))

    def test_illumina13_variant(self):
        # test entire range of possible ascii chars for illumina1.3
        all_illumina13_ascii = ('@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijk'
                                'lmnopqrstuvwxyz{|}~')
        obs = _decode_qual_to_phred(all_illumina13_ascii,
                                    variant='illumina1.3')
        self.assertEqual(obs, list(range(63)))

        with self.assertRaises(ValueError) as cm:
            _decode_qual_to_phred('a!b', variant='illumina1.3')
        self.assertIn('-31', str(cm.exception))
        self.assertIn('[0, 62]', str(cm.exception))

    def test_illumina18_variant(self):
        # test entire range of possible ascii chars for illumina1.8
        all_illumina18_ascii = ('!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKL'
                                'MNOPQRSTUVWXYZ[\\]^_')
        obs = _decode_qual_to_phred(all_illumina18_ascii,
                                    variant='illumina1.8')
        self.assertEqual(obs, list(range(63)))

        with self.assertRaises(ValueError) as cm:
            _decode_qual_to_phred('AaB', variant='illumina1.8')
        self.assertIn('64', str(cm.exception))
        self.assertIn('[0, 62]', str(cm.exception))

    def test_custom_phred_offset(self):
        ascii_chars = '*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\'
        obs = _decode_qual_to_phred(ascii_chars, phred_offset=42)
        self.assertEqual(obs, list(range(51)))

        with self.assertRaises(ValueError) as cm:
            _decode_qual_to_phred(ascii_chars, phred_offset=43)
        self.assertIn('-1', str(cm.exception))
        self.assertIn('[0, 83]', str(cm.exception))

        with self.assertRaises(ValueError) as cm:
            _decode_qual_to_phred(ascii_chars, phred_offset=0)
        self.assertIn('`phred_offset`', str(cm.exception))
        self.assertIn('0', str(cm.exception))
        self.assertIn('printable', str(cm.exception))

        with self.assertRaises(ValueError) as cm:
            _decode_qual_to_phred(ascii_chars, phred_offset=127)
        self.assertIn('`phred_offset`', str(cm.exception))
        self.assertIn('127', str(cm.exception))
        self.assertIn('printable', str(cm.exception))


class PhredEncoderTests(unittest.TestCase):
    def test_missing_variant_and_phred_offset(self):
        with self.assertRaises(ValueError) as cm:
            _encode_phred_to_qual([1, 2, 3])
        self.assertIn('`variant`', str(cm.exception))
        self.assertIn('`phred_offset`', str(cm.exception))
        self.assertIn('encode', str(cm.exception))

    def test_variant_and_phred_offset_provided(self):
        with self.assertRaises(ValueError) as cm:
            _encode_phred_to_qual([1, 2, 3], variant='sanger', phred_offset=64)
        self.assertIn('both', str(cm.exception))
        self.assertIn('`variant`', str(cm.exception))
        self.assertIn('`phred_offset`', str(cm.exception))

    def test_solexa_variant(self):
        with self.assertRaises(NotImplementedError) as cm:
            _encode_phred_to_qual([1, 2, 3], variant='solexa')
        self.assertIn('719', str(cm.exception))

    def test_unrecognized_variant(self):
        with self.assertRaises(ValueError) as cm:
            _encode_phred_to_qual([1, 2, 3], variant='illumina')
        self.assertIn('variant', str(cm.exception))
        self.assertIn("'illumina'", str(cm.exception))

    def test_no_phred_scores(self):
        self.assertEqual(_encode_phred_to_qual([], variant='sanger'), '')

    def test_sanger_variant(self):
        # test entire range of possible ascii chars for sanger
        all_sanger_ascii = ('!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP'
                            'QRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~')
        obs = _encode_phred_to_qual(list(range(94)), variant='sanger')
        self.assertEqual(obs, all_sanger_ascii)

        with self.assertRaises(ValueError) as cm:
            _encode_phred_to_qual([42, -1, 33], variant='sanger')
        self.assertIn('-1', str(cm.exception))
        self.assertIn('[0, 93]', str(cm.exception))

        obs = npt.assert_warns(UserWarning, _encode_phred_to_qual,
                               [42, 94, 33], variant='sanger')
        self.assertEqual(obs, 'K~B')

    def test_illumina13_variant(self):
        # test entire range of possible ascii chars for illumina1.3
        all_illumina13_ascii = ('@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijk'
                                'lmnopqrstuvwxyz{|}~')
        obs = _encode_phred_to_qual(list(range(63)), variant='illumina1.3')
        self.assertEqual(obs, all_illumina13_ascii)

        with self.assertRaises(ValueError) as cm:
            _encode_phred_to_qual([42, -1, 33], variant='illumina1.3')
        self.assertIn('-1', str(cm.exception))
        self.assertIn('[0, 62]', str(cm.exception))

        obs = npt.assert_warns(UserWarning, _encode_phred_to_qual,
                               [42, 63, 33], variant='illumina1.3')
        self.assertEqual(obs, 'j~a')

    def test_illumina18_variant(self):
        # test entire range of possible ascii chars for illumina1.8
        all_illumina18_ascii = ('!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKL'
                                'MNOPQRSTUVWXYZ[\\]^_')
        obs = _encode_phred_to_qual(list(range(63)), variant='illumina1.8')
        self.assertEqual(obs, all_illumina18_ascii)

        with self.assertRaises(ValueError) as cm:
            _encode_phred_to_qual([42, -1, 33], variant='illumina1.8')
        self.assertIn('-1', str(cm.exception))
        self.assertIn('[0, 62]', str(cm.exception))

        obs = npt.assert_warns(UserWarning, _encode_phred_to_qual,
                               [42, 63, 33], variant='illumina1.8')
        self.assertEqual(obs, 'K_B')

    def test_custom_phred_offset(self):
        ascii_chars = '*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\'
        obs = _encode_phred_to_qual(list(range(51)), phred_offset=42)
        self.assertEqual(obs, ascii_chars)

        with self.assertRaises(ValueError) as cm:
            _encode_phred_to_qual([42, -1, 33], phred_offset=42)
        self.assertIn('-1', str(cm.exception))
        self.assertIn('[0, 84]', str(cm.exception))

        obs = npt.assert_warns(UserWarning, _encode_phred_to_qual,
                               [42, 255, 33], phred_offset=42)
        self.assertEqual(obs, 'T~K')


if __name__ == '__main__':
    unittest.main()
