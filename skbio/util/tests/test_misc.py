# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from io import BytesIO
import re

from skbio.util._misc import (
    cardinal_to_ordinal,
    safe_md5,
    find_duplicates,
    MiniRegistry,
    chunk_str,
    resolve_key,
)


class TestMiniRegistry(TestCase):
    def setUp(self):
        self.registry = MiniRegistry()

    def test_decoration(self):
        self.assertNotIn("name1", self.registry)
        self.assertNotIn("name2", self.registry)
        self.n1_called = False
        self.n2_called = False

        @self.registry("name1")
        def some_registration1():
            self.n1_called = True

        @self.registry("name2")
        def some_registration2():
            self.n2_called = True

        self.assertIn("name1", self.registry)
        self.assertEqual(some_registration1, self.registry["name1"])
        self.assertIn("name2", self.registry)
        self.assertEqual(some_registration2, self.registry["name2"])

        self.registry["name1"]()
        self.assertTrue(self.n1_called)
        self.registry["name2"]()
        self.assertTrue(self.n2_called)

    def test_copy(self):
        @self.registry("name")
        def some_registration():
            pass

        new = self.registry.copy()
        self.assertIsNot(new, self.registry)

        @new("other")
        def other_registration():
            pass

        self.assertIn("name", self.registry)
        self.assertNotIn("other", self.registry)

        self.assertIn("other", new)
        self.assertIn("name", new)

    def test_everything(self):
        class SomethingToInterpolate:
            def interpolate_me():
                """First line

                Some description of things, also this:

                Other things are happening now.
                """

            def dont_interpolate_me():
                """First line

                Some description of things, also this:

                Other things are happening now.
                """

        class Subclass(SomethingToInterpolate):
            pass

        @self.registry("a")
        def a():
            """x"""

        @self.registry("b")
        def b():
            """y"""

        @self.registry("c")
        def c():
            """z"""

        subclass_registry = self.registry.copy()

        @subclass_registry("o")
        def o():
            """p"""

        self.registry.interpolate(SomethingToInterpolate, "interpolate_me")
        subclass_registry.interpolate(Subclass, "interpolate_me")

        # Python 3.13+ removes leading whitespaces from a docstring, therefore it is
        # necessary to make this edit.
        def _strip_spaces(s):
            return re.sub(r'^\s+', '', s, flags=re.MULTILINE)

        obs = _strip_spaces(SomethingToInterpolate.interpolate_me.__doc__)
        exp = ("First line\nSome description of things, also this:\n'a'\nx"
               "\n'b'\ny\n'c'\nz\nOther things are happening now.\n")
        self.assertEqual(obs, exp)

        obs = _strip_spaces(SomethingToInterpolate.dont_interpolate_me.__doc__)
        exp = ("First line\nSome description of things, also this:\nOther things "
               "are happening now.\n")
        self.assertEqual(obs, exp)

        obs = _strip_spaces(Subclass.interpolate_me.__doc__)
        exp = ("First line\nSome description of things, also this:\n'a'\nx\n'b'\ny\n"
               "'c'\nz\n'o'\np\nOther things are happening now.\n")
        self.assertEqual(obs, exp)

        obs = _strip_spaces(Subclass.dont_interpolate_me.__doc__)
        exp = ("First line\nSome description of things, also this:\nOther things "
               "are happening now.\n")
        self.assertEqual(obs, exp)


class ResolveKeyTests(TestCase):
    def test_callable(self):
        def func(x):
            return str(x)

        self.assertEqual(resolve_key(1, func), "1")
        self.assertEqual(resolve_key(4, func), "4")

    def test_index(self):
        class MetadataHaver(dict):
            @property
            def metadata(self):
                return self

        obj = MetadataHaver({'foo': 123})
        self.assertEqual(resolve_key(obj, 'foo'), 123)

        obj = MetadataHaver({'foo': 123, 'bar': 'baz'})
        self.assertEqual(resolve_key(obj, 'bar'), 'baz')

    def test_wrong_type(self):
        with self.assertRaises(TypeError):
            resolve_key({'foo': 1}, 'foo')


class ChunkStrTests(TestCase):
    def test_even_split(self):
        self.assertEqual(chunk_str('abcdef', 6, ' '), 'abcdef')
        self.assertEqual(chunk_str('abcdef', 3, ' '), 'abc def')
        self.assertEqual(chunk_str('abcdef', 2, ' '), 'ab cd ef')
        self.assertEqual(chunk_str('abcdef', 1, ' '), 'a b c d e f')
        self.assertEqual(chunk_str('a', 1, ' '), 'a')
        self.assertEqual(chunk_str('abcdef', 2, ''), 'abcdef')

    def test_no_split(self):
        self.assertEqual(chunk_str('', 2, '\n'), '')
        self.assertEqual(chunk_str('a', 100, '\n'), 'a')
        self.assertEqual(chunk_str('abcdef', 42, '|'), 'abcdef')

    def test_uneven_split(self):
        self.assertEqual(chunk_str('abcdef', 5, '|'), 'abcde|f')
        self.assertEqual(chunk_str('abcdef', 4, '|'), 'abcd|ef')
        self.assertEqual(chunk_str('abcdefg', 3, ' - '), 'abc - def - g')

    def test_invalid_n(self):
        with self.assertRaisesRegex(ValueError, r'n=0'):
            chunk_str('abcdef', 0, ' ')

        with self.assertRaisesRegex(ValueError, r'n=-42'):
            chunk_str('abcdef', -42, ' ')


class SafeMD5Tests(TestCase):
    def test_safe_md5(self):
        exp = 'ab07acbb1e496801937adfa772424bf7'

        fd = BytesIO(b'foo bar baz')
        obs = safe_md5(fd)
        self.assertEqual(obs.hexdigest(), exp)

        fd.close()


class CardinalToOrdinalTests(TestCase):
    def test_valid_range(self):
        # taken and modified from http://stackoverflow.com/a/20007730/3776794
        exp = ['0th', '1st', '2nd', '3rd', '4th', '5th', '6th', '7th', '8th',
               '9th', '10th', '11th', '12th', '13th', '14th', '15th', '16th',
               '17th', '18th', '19th', '20th', '21st', '22nd', '23rd', '24th',
               '25th', '26th', '27th', '28th', '29th', '30th', '31st', '32nd',
               '100th', '101st', '42042nd']
        obs = [cardinal_to_ordinal(n) for n in
               list(range(0, 33)) + [100, 101, 42042]]
        self.assertEqual(obs, exp)

    def test_invalid_n(self):
        with self.assertRaisesRegex(ValueError, r'-1'):
            cardinal_to_ordinal(-1)


class TestFindDuplicates(TestCase):
    def test_empty_input(self):
        def empty_gen():
            yield from ()

        for empty in [], (), '', set(), {}, empty_gen():
            self.assertEqual(find_duplicates(empty), set())

    def test_no_duplicates(self):
        self.assertEqual(find_duplicates(['a', 'bc', 'def', 'A']), set())

    def test_one_duplicate(self):
        self.assertEqual(find_duplicates(['a', 'bc', 'def', 'a']), set(['a']))

    def test_many_duplicates(self):
        self.assertEqual(find_duplicates(['a', 'bc', 'bc', 'def', 'a']),
                         set(['a', 'bc']))

    def test_all_duplicates(self):
        self.assertEqual(
            find_duplicates(('a', 'bc', 'bc', 'def', 'a', 'def', 'def')),
            set(['a', 'bc', 'def']))

    def test_mixed_types(self):
        def gen():
            yield from ('a', 1, 'bc', 2, 'a', 2, 2, 3.0)

        self.assertEqual(find_duplicates(gen()), set(['a', 2]))


if __name__ == '__main__':
    main()
