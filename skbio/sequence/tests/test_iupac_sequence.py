# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio.sequence._iupac_sequence import IUPACSequence
from skbio.util import classproperty

class IUPACSequenceSubclassNoImplementation(IUPACSequence):
    pass

class ExampleIUPACSequence(IUPACSequence):

    @classproperty
    def degenerate_map(self):
        return {"X": set("AB"), "Y": set("BC"), "Z": set("AC")}
    @classproperty
    def nondegenerate_chars(self):
        return set("ABC")

class TestIUPACSequence(TestCase):

    def test_instantiation_with_no_implementation(self):
        with self.assertRaises(TypeError) as cm:
            t = IUPACSequenceSubclassNoImplementation()

        self.assertIn("abstract class", str(cm.exception))
        self.assertIn("nondegenerate_chars", str(cm.exception))
        self.assertIn("degenerate_map", str(cm.exception))

    def test_degenerate_chars(self):
        expected = set("XYZ")
        self.assertIs(type(ExampleIUPACSequence.degenerate_chars), set)
        self.assertEqual(ExampleIUPACSequence.degenerate_chars, expected)

        ExampleIUPACSequence.degenerate_chars.add("W")
        self.assertEqual(ExampleIUPACSequence.degenerate_chars, expected)

        self.assertEqual(ExampleIUPACSequence('').degenerate_chars, expected)

        with self.assertRaises(AttributeError):
            ExampleIUPACSequence('').degenerate_chars = set("BAR")

    def test_nondegenerate_chars(self):
        expected = set("ABC")
        self.assertEqual(ExampleIUPACSequence.nondegenerate_chars, expected)

        ExampleIUPACSequence.degenerate_chars.add("D")
        self.assertEqual(ExampleIUPACSequence.nondegenerate_chars, expected)

        self.assertEqual(ExampleIUPACSequence('').nondegenerate_chars,
                         expected)

        with self.assertRaises(AttributeError):
            ExampleIUPACSequence('').nondegenerate_chars = set("BAR")

    def test_gap_chars(self):
        expected = set(".-")
        self.assertIs(type(ExampleIUPACSequence.gap_chars), set)
        self.assertEqual(ExampleIUPACSequence.gap_chars, expected)

        ExampleIUPACSequence.gap_chars.add("_")
        self.assertEqual(ExampleIUPACSequence.gap_chars, expected)

        self.assertEqual(ExampleIUPACSequence('').gap_chars, expected)

        with self.assertRaises(AttributeError):
            ExampleIUPACSequence('').gap_chars = set("_ =")

    def test_alphabet(self):
        expected = set("ABC.-XYZ")
        self.assertIs(type(ExampleIUPACSequence.alphabet), set)
        self.assertEqual(ExampleIUPACSequence.alphabet, expected)

        ExampleIUPACSequence.alphabet.add("DEF")
        self.assertEqual(ExampleIUPACSequence.alphabet, expected)

        self.assertEqual(ExampleIUPACSequence('').alphabet, expected)

        with self.assertRaises(AttributeError):
            ExampleIUPACSequence('').alphabet = set("ABCDEFG.-WXYZ")

    def test_degenerate_map(self):
        expected = {"X": set("AB"), "Y": set("BC"), "Z": set("AC")}
        self.assertEqual(ExampleIUPACSequence.degenerate_map, expected)

        ExampleIUPACSequence.degenerate_map['W'] = set("ABC")
        ExampleIUPACSequence.degenerate_map['X'] = set("CA")
        self.assertEqual(ExampleIUPACSequence.degenerate_map, expected)

        self.assertEqual(ExampleIUPACSequence('').degenerate_map, expected)

        with self.assertRaises(AttributeError):
            ExampleIUPACSequence('').degenerate_map = {'W': "ABC"}

    def test_gaps(self):
        self.assertIs(type(ExampleIUPACSequence("").gaps()), np.ndarray)
        self.assertIs(ExampleIUPACSequence("").gaps().dtype, np.dtype('bool'))
        npt.assert_equal(ExampleIUPACSequence("ABCXBZYABC").gaps(),
                         np.zeros(10).astype(bool))

        npt.assert_equal(ExampleIUPACSequence(".-.-.").gaps(),
                         np.ones(5).astype(bool))

        npt.assert_equal(ExampleIUPACSequence("A.B-C.X-Y.").gaps(),
                         np.array([0, 1] * 5, dtype=bool))

        npt.assert_equal(ExampleIUPACSequence("AB.AC.XY-").gaps(),
                         np.array([0, 0, 1] * 3, dtype=bool))

        npt.assert_equal(ExampleIUPACSequence("A.BC.-").gaps(),
                         np.array([0, 1, 0, 0, 1, 1], dtype=bool))

    def test_has_gaps(self):
        self.assertIs(type(ExampleIUPACSequence("").has_gaps()), bool)
        self.assertIs(type(ExampleIUPACSequence("-").has_gaps()), bool)

        self.assertFalse(ExampleIUPACSequence("").has_gaps())
        self.assertFalse(ExampleIUPACSequence("ABCXYZ").has_gaps())

        self.assertTrue(ExampleIUPACSequence("-").has_gaps())
        self.assertTrue(ExampleIUPACSequence("ABCXYZ-").has_gaps())

    def test_degenerates(self):
        self.assertIs(type(ExampleIUPACSequence("").degenerates()), np.ndarray)
        self.assertIs(ExampleIUPACSequence("").degenerates().dtype,
                      np.dtype('bool'))

        npt.assert_equal(ExampleIUPACSequence("ABCBC-.AB.").degenerates(),
                         np.zeros(10).astype(bool))

        npt.assert_equal(ExampleIUPACSequence("ZYZYZ").degenerates(),
                         np.ones(5).astype(bool))

        npt.assert_equal(ExampleIUPACSequence("AX.Y-ZBXCZ").degenerates(),
                         np.array([0, 1] * 5, dtype=bool))

        npt.assert_equal(ExampleIUPACSequence("ABXACY.-Z").degenerates(),
                         np.array([0, 0, 1] * 3, dtype=bool))

        npt.assert_equal(ExampleIUPACSequence("AZBCXY").degenerates(),
                         np.array([0, 1, 0, 0, 1, 1], dtype=bool))

    def test_has_degenerates(self):
        self.assertIs(type(ExampleIUPACSequence("").has_degenerates()), bool)
        self.assertIs(type(ExampleIUPACSequence("X").has_degenerates()), bool)

        self.assertFalse(ExampleIUPACSequence("").has_degenerates())
        self.assertFalse(ExampleIUPACSequence("A-.BC").has_degenerates())

        self.assertTrue(ExampleIUPACSequence("Z").has_degenerates())
        self.assertTrue(ExampleIUPACSequence("ABC.XYZ-").has_degenerates())


if __name__ == "__main__":
    main()
