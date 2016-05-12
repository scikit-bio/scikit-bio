# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections.abc
import types
import unittest

import numpy as np

from skbio.stats.distance import PairwiseDistances
from skbio.util._testing import ReallyEqualMixin


class TestMappingType(unittest.TestCase):
    def test_isinstance_empty(self):
        d = PairwiseDistances({})

        self.assertIsInstance(d, collections.abc.Mapping)

    def test_isinstance_non_empty(self):
        d = PairwiseDistances({('a', 'b'): 42.0, ('b', 'c'): 43.0})

        self.assertIsInstance(d, collections.abc.Mapping)

    def test_immutable_mapping(self):
        d = PairwiseDistances({('a', 'b'): 42.0, ('b', 'c'): 43.0})

        self.assertNotIsInstance(d, collections.abc.MutableMapping)
        with self.assertRaises(TypeError):
            d[('a', 'b')] = 42.5


class TestConstructor(unittest.TestCase):
    def test_invalid_input_type(self):
        with self.assertRaisesRegex(TypeError, 'dict.*set'):
            PairwiseDistances({42.0, 43.0, 44.0})

    def test_invalid_dict_key_type(self):
        with self.assertRaisesRegex(TypeError, "key.*tuple.*'a'"):
            PairwiseDistances({'a': 42.0, ('b', 'c'): 43.0})

        with self.assertRaisesRegex(TypeError, "key.*tuple.*('a', 'b', 'c')"):
            PairwiseDistances({('a', 'b', 'c'): 42.0, ('b', 'c'): 43.0})

    def test_invalid_distance_type(self):
        with self.assertRaisesRegex(TypeError, 'distance.*float.*str'):
            PairwiseDistances({('a', 'c'): 'abc', ('b', 'a'): 43.0})

        with self.assertRaisesRegex(TypeError, 'distance.*float.*int'):
            PairwiseDistances({('a', 'c'): 42, ('b', 'a'): 43.0})

    def test_nan_distance(self):
        with self.assertRaisesRegex(ValueError, 'distance.*NaN'):
            PairwiseDistances({('a', 'c'): 42.0, ('b', 'a'): np.nan})

    def test_non_hollow_distances(self):
        with self.assertRaisesRegex(ValueError, "nonzero distance.*'a'.*42.0"):
            PairwiseDistances({('a', 'a'): 42.0, ('b', 'a'): 43.0})

    def test_asymmetric_distances(self):
        with self.assertRaises(ValueError) as cm:
            PairwiseDistances({('a', 'b'): 42.0, ('b', 'a'): 43.0})

        # Not using assertRaisesRegex because order of IDs isn't guaranteed.
        error = str(cm.exception)
        self.assertIn("asymmetric distance", error)
        self.assertIn("'a'", error)
        self.assertIn("'b'", error)
        self.assertIn("42.0", error)
        self.assertIn("43.0", error)

    def test_empty(self):
        d = PairwiseDistances({})

        self.assertEqual(d.ids, set())
        self.assertEqual(len(d), 0)

    def test_single_within_distance(self):
        d = PairwiseDistances({('abc', 'abc'): 0.0})

        self.assertEqual(d.ids, {'abc'})
        self.assertEqual(len(d), 1)
        self.assertEqual(d['abc', 'abc'], 0.0)

    def test_single_between_distance(self):
        d = PairwiseDistances({('abc', 'x'): 0.9})

        self.assertEqual(d.ids, {'abc', 'x'})
        self.assertEqual(len(d), 3)
        self.assertEqual(d['abc', 'abc'], 0.0)
        self.assertEqual(d['x', 'x'], 0.0)
        self.assertEqual(d['abc', 'x'], 0.9)

    def test_multiple_distances(self):
        d = PairwiseDistances({('a', 'b'): 1.0,
                               ('a', 'c'): 42.5,
                               ('c', 'd'): 0.0})

        self.assertEqual(d.ids, {'a', 'b', 'c', 'd'})
        self.assertEqual(len(d), 7)
        self.assertEqual(d['a', 'a'], 0.0)
        self.assertEqual(d['b', 'b'], 0.0)
        self.assertEqual(d['c', 'c'], 0.0)
        self.assertEqual(d['d', 'd'], 0.0)
        self.assertEqual(d['a', 'b'], 1.0)
        self.assertEqual(d['a', 'c'], 42.5)
        self.assertEqual(d['c', 'd'], 0.0)

    def test_redundant_symmetric_distances(self):
        d = PairwiseDistances({('a', 'b'): 42.0, ('b', 'a'): 42.0})

        self.assertEqual(d.ids, {'a', 'b'})
        self.assertEqual(len(d), 3)
        self.assertEqual(d['a', 'a'], 0.0)
        self.assertEqual(d['b', 'b'], 0.0)
        self.assertEqual(d['a', 'b'], 42.0)

    def test_redundant_hollow_distances(self):
        d = PairwiseDistances({('a', 'a'): 0.0,
                               ('b', 'a'): 42.0,
                               ('b', 'b'): 0.0})

        self.assertEqual(d.ids, {'a', 'b'})
        self.assertEqual(len(d), 3)
        self.assertEqual(d['a', 'a'], 0.0)
        self.assertEqual(d['b', 'b'], 0.0)
        self.assertEqual(d['a', 'b'], 42.0)

    def test_constructor_copies_distances(self):
        distances = {('a', 'b'): 12.3, ('c', 'b'): 3.21}
        d = PairwiseDistances(distances)

        # Mutating `distances` dict after constructing PairwiseDistances
        # shouldn't mutate PairwiseDistances.
        distances['c', 'b'] = 100.0
        distances['a', 'd'] = 42.0

        self.assertEqual(d.ids, {'a', 'b', 'c'})
        self.assertEqual(len(d), 5)
        self.assertEqual(d['a', 'a'], 0.0)
        self.assertEqual(d['b', 'b'], 0.0)
        self.assertEqual(d['c', 'c'], 0.0)
        self.assertEqual(d['a', 'b'], 12.3)
        self.assertEqual(d['c', 'b'], 3.21)


class TestIDs(unittest.TestCase):
    def test_return_type(self):
        d = PairwiseDistances({('a', 'b'): 1.2345})

        self.assertIsInstance(d.ids, frozenset)

    def test_read_only_property(self):
        d = PairwiseDistances({('a', 'b'): 1.2345})

        self.assertEqual(d.ids, {'a', 'b'})

        with self.assertRaises(AttributeError):
            d.ids = {'c', 'd'}

        self.assertEqual(d.ids, {'a', 'b'})

    def test_empty(self):
        d = PairwiseDistances({})

        self.assertEqual(d.ids, set())

    def test_single_within_distance(self):
        d = PairwiseDistances({('a', 'a'): 0.0})

        self.assertEqual(d.ids, {'a'})

    def test_single_between_distance(self):
        d = PairwiseDistances({('a', 'b'): 42.0})

        self.assertEqual(d.ids, {'a', 'b'})

    def test_multiple_distances(self):
        d = PairwiseDistances({('a', 'b'): 1.0,
                               ('a', 'c'): 42.5,
                               ('c', 'd'): 0.0,
                               ('e', 'e'): 0.0})

        self.assertEqual(d.ids, {'a', 'b', 'c', 'd', 'e'})


class TestBool(unittest.TestCase):
    def test_empty(self):
        d = PairwiseDistances({})

        self.assertFalse(bool(d))

    def test_single_within_distance(self):
        d = PairwiseDistances({('xyz', 'xyz'): 0.0})

        self.assertTrue(bool(d))

    def test_single_between_distance(self):
        d = PairwiseDistances({('xyz', '123'): 0.5})

        self.assertTrue(bool(d))

    def test_multiple_distances(self):
        d = PairwiseDistances({('a', 'b'): 1.0,
                               ('a', 'c'): 42.5,
                               ('c', 'd'): 0.0,
                               ('e', 'e'): 0.0})

        self.assertTrue(bool(d))


class TestLen(unittest.TestCase):
    def test_empty(self):
        d = PairwiseDistances({})

        self.assertEqual(len(d), 0)

    def test_single_within_distance(self):
        d = PairwiseDistances({('xyz', 'xyz'): 0.0})

        self.assertEqual(len(d), 1)

    def test_single_between_distance(self):
        d = PairwiseDistances({('xyz', '123'): 0.5})

        self.assertEqual(len(d), 3)

    def test_multiple_redundant_distances(self):
        d = PairwiseDistances({('a', 'b'): 1.0,
                               ('a', 'c'): 42.5,
                               ('a', 'a'): 0.0,
                               ('c', 'a'): 42.5,
                               ('c', 'd'): 0.0,
                               ('d', 'c'): 0.0,
                               ('e', 'e'): 0.0})

        self.assertEqual(len(d), 8)


class TestIter(unittest.TestCase):
    def test_return_type(self):
        d = PairwiseDistances({('a', 'c'): 1.0, ('a', 'd'): 2.0})

        self.assertIsInstance(iter(d), types.GeneratorType)

    def test_empty(self):
        d = PairwiseDistances({})

        id_pairs = list(iter(d))
        self.assertEqual(id_pairs, [])

    def test_single_within_distance(self):
        d = PairwiseDistances({('xyz', 'xyz'): 0.0})

        id_pairs = list(iter(d))

        self.assertEqual(id_pairs, [('xyz', 'xyz')])

    def test_single_between_distance(self):
        d = PairwiseDistances({('xyz', '123'): 0.5})

        id_pairs = list(iter(d))

        self.assertEqual(len(id_pairs), 3)

        id_pairs = set(id_pairs)

        self.assertEqual(id_pairs, {('xyz', 'xyz'),
                                    ('123', '123'),
                                    ('xyz', '123')})

    def test_multiple_redundant_distances(self):
        d = PairwiseDistances({('a', 'b'): 1.0,
                               ('a', 'c'): 42.5,
                               ('a', 'a'): 0.0,
                               ('c', 'a'): 42.5,
                               ('c', 'd'): 0.0,
                               ('e', 'e'): 0.0})

        id_pairs = list(iter(d))

        self.assertEqual(len(id_pairs), 8)

        id_pairs = set(id_pairs)

        exp = {
            ('a', 'a'),
            ('b', 'b'),
            ('c', 'c'),
            ('d', 'd'),
            ('e', 'e'),
            ('a', 'b'),
            ('c', 'd'),
        }

        exp1 = exp | {('a', 'c')}
        exp2 = exp | {('c', 'a')}

        self.assertTrue((id_pairs == exp1) or (id_pairs == exp2))


class TestReversed(unittest.TestCase):
    def test_irreversible(self):
        d = PairwiseDistances({('c', 'd'): 1.0, ('d', 'e'): 3.0})

        with self.assertRaises(TypeError):
            reversed(d)


class TestHash(unittest.TestCase):
    def test_unhashable(self):
        d = PairwiseDistances({('c', 'd'): 1.0, ('d', 'e'): 3.0})

        with self.assertRaises(TypeError):
            hash(d)


class TestGetItem(unittest.TestCase):
    def test_non_id_pair(self):
        d = PairwiseDistances({('a', 'b'): 1.0,
                               ('a', 'c'): 42.5,
                               ('c', 'd'): 0.0})

        with self.assertRaisesRegex(TypeError, "id_pair.*tuple.*'a'"):
            d['a']

        with self.assertRaisesRegex(TypeError, 'id_pair.*tuple.*{.*}'):
            d[{'a', 'b'}]

        with self.assertRaisesRegex(TypeError,
                                    "id_pair.*tuple.*('a', 'b', 'c')"):
            d['a', 'b', 'c']

    def test_non_existing_between_distance(self):
        d = PairwiseDistances({('a', 'b'): 1.0,
                               ('a', 'c'): 42.5,
                               ('c', 'd'): 0.0})

        with self.assertRaises(KeyError):
            d['d', 'a']

        with self.assertRaises(KeyError):
            d['a', 'd']

    def test_non_existing_within_distance(self):
        d = PairwiseDistances({('a', 'b'): 1.0,
                               ('a', 'c'): 42.5,
                               ('c', 'd'): 0.0})

        with self.assertRaises(KeyError):
            d['e', 'e']

    def test_empty(self):
        d = PairwiseDistances({})

        with self.assertRaises(KeyError):
            d['', '']

        with self.assertRaises(KeyError):
            d['a', 'b']

    def test_return_type_within_distance(self):
        d = PairwiseDistances({('a', 'b'): 42.5})

        self.assertIsInstance(d['a', 'a'], float)

    def test_return_type_between_distance(self):
        d = PairwiseDistances({('a', 'b'): 42.5})

        self.assertIsInstance(d['a', 'b'], float)
        self.assertIsInstance(d['b', 'a'], float)

    def test_between_distances(self):
        d = PairwiseDistances({('a', 'b'): 1.0,
                               ('a', 'c'): 42.5,
                               ('c', 'd'): 0.0})

        self.assertEqual(d['a', 'b'], 1.0)
        self.assertEqual(d['b', 'a'], 1.0)

        self.assertEqual(d['a', 'c'], 42.5)
        self.assertEqual(d['c', 'a'], 42.5)

        self.assertEqual(d['c', 'd'], 0.0)
        self.assertEqual(d['d', 'c'], 0.0)

    def test_within_distances(self):
        d = PairwiseDistances({('a', 'b'): 1.0,
                               ('a', 'c'): 42.5,
                               ('c', 'd'): 0.0})

        self.assertEqual(d['a', 'a'], 0.0)
        self.assertEqual(d['b', 'b'], 0.0)
        self.assertEqual(d['c', 'c'], 0.0)
        self.assertEqual(d['d', 'd'], 0.0)


class TestEq(unittest.TestCase, ReallyEqualMixin):
    def test_subclass(self):
        class PairwiseDistancesSubclass(PairwiseDistances):
            pass

        d1 = PairwiseDistances({})
        d2 = PairwiseDistancesSubclass({})

        self.assertReallyEqual(d1, d2)

    def test_empty(self):
        d1 = PairwiseDistances({})
        d2 = PairwiseDistances({})

        self.assertReallyEqual(d1, d2)

    def test_single_within_distance(self):
        d1 = PairwiseDistances({('xyz', 'xyz'): 0.0})
        d2 = PairwiseDistances({('xyz', 'xyz'): 0.0})

        self.assertReallyEqual(d1, d2)

    def test_single_between_distance(self):
        d1 = PairwiseDistances({('xyz', '123'): 0.5})
        d2 = PairwiseDistances({('xyz', '123'): 0.5})

        self.assertReallyEqual(d1, d2)

    def test_multiple_distances_partial(self):
        d1 = PairwiseDistances({('a', 'b'): 0.5, ('a', 'c'): 42.0})
        d2 = PairwiseDistances({('a', 'b'): 0.5, ('a', 'c'): 42.0})

        self.assertReallyEqual(d1, d2)

    def test_multiple_distances_full(self):
        d1 = PairwiseDistances({('a', 'b'): 0.5,
                                ('a', 'c'): 42.0,
                                ('b', 'c'): 0.9})
        d2 = PairwiseDistances({('a', 'b'): 0.5,
                                ('a', 'c'): 42.0,
                                ('b', 'c'): 0.9})

        self.assertReallyEqual(d1, d2)

    def test_isomorphic_id_pairs(self):
        d1 = PairwiseDistances({('a', 'b'): 0.5,
                                ('a', 'c'): 42.0,
                                ('b', 'd'): 0.9})
        d2 = PairwiseDistances({('b', 'a'): 0.5,
                                ('c', 'a'): 42.0,
                                ('d', 'b'): 0.9})

        self.assertReallyEqual(d1, d2)

    def test_redundant_isomorphic_distances(self):
        d1 = PairwiseDistances({('a', 'b'): 0.5,
                                ('a', 'c'): 42.0})
        d2 = PairwiseDistances({('a', 'b'): 0.5,
                                ('b', 'a'): 0.5,
                                ('a', 'c'): 42.0,
                                ('c', 'a'): 42.0,
                                ('a', 'a'): 0.0,
                                ('b', 'b'): 0.0,
                                ('c', 'c'): 0.0})

        self.assertReallyEqual(d1, d2)


class TestNe(unittest.TestCase, ReallyEqualMixin):
    def test_different_type(self):
        d = PairwiseDistances({})

        self.assertReallyNotEqual(d, {})

    def test_different_length(self):
        d1 = PairwiseDistances({})
        d2 = PairwiseDistances({('a', 'a'): 0.0})

        self.assertReallyNotEqual(d1, d2)

    def test_different_ids(self):
        d1 = PairwiseDistances({('a', 'a'): 0.0, ('b', 'b'): 0.0})
        d2 = PairwiseDistances({('a', 'a'): 0.0, ('c', 'c'): 0.0})

        self.assertReallyNotEqual(d1, d2)

    def test_different_between_distances(self):
        d1 = PairwiseDistances({('a', 'b'): 42.0, ('a', 'c'): 43.0})
        d2 = PairwiseDistances({('a', 'b'): 42.0, ('a', 'c'): 44.0})

        self.assertReallyNotEqual(d1, d2)

    def test_different_between_id_pairs(self):
        d1 = PairwiseDistances({('a', 'b'): 42.0, ('a', 'c'): 43.0})
        d2 = PairwiseDistances({('a', 'b'): 42.0, ('b', 'c'): 1.0})

        self.assertReallyNotEqual(d1, d2)


# Tests for mixin methods provided by ABC and not overridden by
# PairwiseDistances.

class TestGet(unittest.TestCase):
    def setUp(self):
        self.d = PairwiseDistances({('A', 'B'): 42.0, ('A', 'C'): 43.0})

    def test_non_id_pair(self):
        with self.assertRaisesRegex(TypeError, "`id_pair`.*tuple.*'A'"):
            self.d.get('A')

    def test_existing_within_distance(self):
        self.assertEqual(self.d.get(('C', 'C')), 0.0)

    def test_non_existing_within_distance(self):
        self.assertIsNone(self.d.get(('D', 'D')))

    def test_existing_between_distance(self):
        self.assertEqual(self.d.get(('A', 'C')), 43.0)
        self.assertEqual(self.d.get(('C', 'A')), 43.0)

    def test_non_existing_between_distance(self):
        self.assertIsNone(self.d.get(('B', 'C')))
        self.assertIsNone(self.d.get(('C', 'B')))

    def test_non_existing_distance_non_default(self):
        self.assertEqual(self.d.get(('B', 'C'), default='foo'), 'foo')
        self.assertEqual(self.d.get(('C', 'B'), default='foo'), 'foo')


class TestContains(unittest.TestCase):
    def setUp(self):
        self.d = PairwiseDistances({('A', 'B'): 42.0, ('A', 'C'): 43.0})

    def test_non_id_pair(self):
        with self.assertRaisesRegex(TypeError, "`id_pair`.*tuple.*'A'"):
            'A' in self.d

    def test_exists_within(self):
        self.assertTrue(('A', 'A') in self.d)
        self.assertTrue(('B', 'B') in self.d)
        self.assertTrue(('C', 'C') in self.d)

    def test_not_exists_within(self):
        self.assertFalse(('D', 'D') in self.d)
        self.assertFalse(('E', 'E') in self.d)

    def test_exists_between(self):
        self.assertTrue(('A', 'B') in self.d)
        self.assertTrue(('B', 'A') in self.d)
        self.assertTrue(('A', 'C') in self.d)
        self.assertTrue(('C', 'A') in self.d)

    def test_not_exists_between(self):
        self.assertFalse(('B', 'C') in self.d)
        self.assertFalse(('C', 'B') in self.d)
        self.assertFalse(('X', 'Y') in self.d)


class TestKeys(unittest.TestCase):
    def test_empty(self):
        d = PairwiseDistances({})

        keys = d.keys()

        self.assertIsInstance(keys, collections.abc.KeysView)
        self.assertEqual(list(keys), [])

    def test_non_empty(self):
        d = PairwiseDistances({('A', 'B'): 42.0, ('A', 'C'): 43.0})

        keys = d.keys()

        self.assertIsInstance(keys, collections.abc.KeysView)
        self.assertEqual(len(keys), 5)

        keys = set(keys)
        self.assertEqual(
            keys, {('A', 'B'), ('A', 'C'), ('A', 'A'), ('B', 'B'), ('C', 'C')})


class TestItems(unittest.TestCase):
    def test_empty(self):
        d = PairwiseDistances({})

        items = d.items()

        self.assertIsInstance(items, collections.abc.ItemsView)
        self.assertEqual(list(items), [])

    def test_non_empty(self):
        d = PairwiseDistances({('A', 'B'): 42.0, ('A', 'C'): 43.0})

        items = d.items()

        self.assertIsInstance(items, collections.abc.ItemsView)
        self.assertEqual(len(items), 5)

        items = sorted(items)
        self.assertEqual(items, [(('A', 'A'), 0.0),
                                 (('A', 'B'), 42.0),
                                 (('A', 'C'), 43.0),
                                 (('B', 'B'), 0.0),
                                 (('C', 'C'), 0.0)])


class TestValues(unittest.TestCase):
    def test_empty(self):
        d = PairwiseDistances({})

        values = d.values()

        self.assertIsInstance(values, collections.abc.ValuesView)
        self.assertEqual(list(values), [])

    def test_non_empty(self):
        d = PairwiseDistances({('A', 'B'): 42.0, ('A', 'C'): 43.0})

        values = d.values()

        self.assertIsInstance(values, collections.abc.ValuesView)
        self.assertEqual(len(values), 5)

        values = sorted(values)
        self.assertEqual(values, [0.0, 0.0, 0.0, 42.0, 43.0])


class TestReprAndStr(unittest.TestCase):
    # Basic sanity checks -- more extensive testing of formatting and special
    # cases is performed in PairwiseDistancesReprDoctests below. Here we only
    # test that pieces of the repr are present. These tests also exercise
    # coverage in case doctests stop counting towards coverage in the future.

    def test_empty(self):
        d = PairwiseDistances({})

        obs = repr(d)

        self.assertEqual(obs, str(d))
        self.assertTrue(obs.startswith('PairwiseDistances'))
        self.assertIn('ID count: 0', obs)
        self.assertIn('"Between" distances count: 0', obs)

    def test_single_within_distance(self):
        d = PairwiseDistances({('A', 'A'): 0.0})

        obs = repr(d)

        self.assertEqual(obs, str(d))
        self.assertTrue(obs.startswith('PairwiseDistances'))
        self.assertIn('ID count: 1', obs)
        self.assertIn('"Between" distances count: 0', obs)
        self.assertIn("'A': 0", obs)

    def test_single_between_distance(self):
        d = PairwiseDistances({('A', 'B'): 0.5})

        obs = repr(d)

        self.assertEqual(obs, str(d))
        self.assertTrue(obs.startswith('PairwiseDistances'))
        self.assertIn('ID count: 2', obs)
        self.assertIn('"Between" distances count: 1', obs)
        self.assertIn("'A': 1", obs)
        self.assertIn("'B': 1", obs)

    def test_multiple_distances_truncated(self):
        d = PairwiseDistances({('A', 'B'): 0.5,
                               ('A', 'C'): 0.4,
                               ('A', 'D'): 0.9,
                               ('A', 'E'): 0.9,
                               ('F', 'G'): 0.9,
                               ('H', 'I'): 0.9})

        obs = repr(d)

        self.assertEqual(obs, str(d))
        self.assertTrue(obs.startswith('PairwiseDistances'))
        self.assertIn('ID count: 9', obs)
        self.assertIn('"Between" distances count: 6', obs)
        self.assertIn("'A': 4", obs)
        self.assertIn("...", obs)
        self.assertIn("'I': 1", obs)


# NOTE: this must be a *separate* class for doctests only (no unit tests). nose
# will not run the unit tests otherwise.
#
# These doctests exercise the correct formatting of PairwiseDistances's repr in
# a variety of situations. They are more extensive than the unit tests above
# (TestReprAndStr) but cannot be relied upon for coverage (the unit tests take
# care of this).
class PairwiseDistancesReprDoctests:
    r"""
    >>> from skbio.stats.distance import PairwiseDistances

    Empty (minimal):

    >>> PairwiseDistances({})
    PairwiseDistances
    --------------------------------------
    Stats:
        ID count: 0
        "Between" distances count: 0
        Total "between" distances count: 0
        Percent complete: 100.00%
    --------------------------------------

    Single "within" distance:

    >>> PairwiseDistances({('A', 'A'): 0.0})
    PairwiseDistances
    --------------------------------------
    Stats:
        ID count: 1
        "Between" distances count: 0
        Total "between" distances count: 0
        Percent complete: 100.00%
    --------------------------------------
    'A': 0 "between" distances

    Single "between" distance:

    >>> PairwiseDistances({('A', 'B'): 0.5})
    PairwiseDistances
    --------------------------------------
    Stats:
        ID count: 2
        "Between" distances count: 1
        Total "between" distances count: 1
        Percent complete: 100.00%
    --------------------------------------
    'A': 1 "between" distance
    'B': 1 "between" distance

    Multiple "within" distances (no "between" distances):

    >>> PairwiseDistances({('ABC', 'ABC'): 0.0,
    ...                    ('DEF', 'DEF'): 0.0,
    ...                    ('GHI', 'GHI'): 0.0})
    PairwiseDistances
    --------------------------------------
    Stats:
        ID count: 3
        "Between" distances count: 0
        Total "between" distances count: 3
        Percent complete: 0.00%
    --------------------------------------
    'ABC': 0 "between" distances
    'DEF': 0 "between" distances
    'GHI': 0 "between" distances

    Multiple "between" distances (partial):

    >>> PairwiseDistances({('A', 'B'): 0.5, ('A', 'C'): 0.4, ('D', 'D'): 0.0})
    PairwiseDistances
    --------------------------------------
    Stats:
        ID count: 4
        "Between" distances count: 2
        Total "between" distances count: 6
        Percent complete: 33.33%
    --------------------------------------
    'A': 2 "between" distances
    'B': 1 "between" distance
    'C': 1 "between" distance
    'D': 0 "between" distances

    Multiple "between" distances (full):

    >>> PairwiseDistances({('A', 'B'): 0.5, ('A', 'C'): 0.4, ('B', 'C'): 0.9})
    PairwiseDistances
    --------------------------------------
    Stats:
        ID count: 3
        "Between" distances count: 3
        Total "between" distances count: 3
        Percent complete: 100.00%
    --------------------------------------
    'A': 2 "between" distances
    'B': 2 "between" distances
    'C': 2 "between" distances

    Five IDs (max untruncated output):

    >>> PairwiseDistances({('A', 'B'): 0.5,
    ...                    ('B', 'C'): 0.4,
    ...                    ('C', 'D'): 0.9,
    ...                    ('D', 'E'): 0.9})
    PairwiseDistances
    ---------------------------------------
    Stats:
        ID count: 5
        "Between" distances count: 4
        Total "between" distances count: 10
        Percent complete: 40.00%
    ---------------------------------------
    'A': 1 "between" distance
    'B': 2 "between" distances
    'C': 2 "between" distances
    'D': 2 "between" distances
    'E': 1 "between" distance

    More than five IDs (truncated output):

    >>> PairwiseDistances({('A', 'B'): 0.5,
    ...                    ('A', 'C'): 0.4,
    ...                    ('A', 'D'): 0.9,
    ...                    ('A', 'E'): 0.9,
    ...                    ('F', 'G'): 0.9,
    ...                    ('H', 'I'): 0.9})
    PairwiseDistances
    ---------------------------------------
    Stats:
        ID count: 9
        "Between" distances count: 6
        Total "between" distances count: 36
        Percent complete: 16.67%
    ---------------------------------------
    'A': 4 "between" distances
    'B': 1 "between" distance
    ...
    'H': 1 "between" distance
    'I': 1 "between" distance

    """
    pass


if __name__ == '__main__':
    unittest.main()
