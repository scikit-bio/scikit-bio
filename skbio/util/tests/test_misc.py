# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import range
from six import BytesIO

import unittest
from tempfile import NamedTemporaryFile, mkdtemp
from os.path import exists, join
from shutil import rmtree
from uuid import uuid4

from skbio.util import (cardinal_to_ordinal, safe_md5, remove_files,
                        create_dir, find_duplicates, flatten,
                        is_casava_v180_or_later)
from skbio.util._misc import _handle_error_codes, MiniRegistry, chunk_str


class TestMiniRegistry(unittest.TestCase):
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
        class SomethingToInterpolate(object):
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

        self.assertEqual(SomethingToInterpolate.interpolate_me.__doc__,
                         "First line\n\n                Some description of th"
                         "ings, also this:\n\n\t'a'\n\t  x\n\t'b'\n\t  y\n\t'c"
                         "'\n\t  z\n\n                Other things are happeni"
                         "ng now.\n                ")
        self.assertEqual(SomethingToInterpolate.dont_interpolate_me.__doc__,
                         "First line\n\n                Some description of th"
                         "ings, also this:\n\n                Other things are"
                         " happening now.\n                ")
        self.assertEqual(Subclass.interpolate_me.__doc__,
                         "First line\n\n                Some description of th"
                         "ings, also this:\n\n\t'a'\n\t  x\n\t'b'\n\t  y\n\t'c"
                         "'\n\t  z\n\t'o'\n\t  p\n\n                Other thin"
                         "gs are happening now.\n                ")
        self.assertEqual(Subclass.dont_interpolate_me.__doc__,
                         "First line\n\n                Some description of th"
                         "ings, also this:\n\n                Other things are"
                         " happening now.\n                ")


class ChunkStrTests(unittest.TestCase):
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
        with self.assertRaisesRegexp(ValueError, 'n=0'):
            chunk_str('abcdef', 0, ' ')

        with self.assertRaisesRegexp(ValueError, 'n=-42'):
            chunk_str('abcdef', -42, ' ')


class MiscTests(unittest.TestCase):
    def setUp(self):
        self.dirs_to_remove = []

    def tearDown(self):
        for element in self.dirs_to_remove:
            rmtree(element)

    def test_is_casava_v180_or_later(self):
        self.assertFalse(is_casava_v180_or_later(b'@foo'))
        id_ = b'@M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0'
        self.assertTrue(is_casava_v180_or_later(id_))

        with self.assertRaises(ValueError):
            is_casava_v180_or_later(b'foo')

    def test_safe_md5(self):
        exp = 'ab07acbb1e496801937adfa772424bf7'

        fd = BytesIO(b'foo bar baz')
        obs = safe_md5(fd)
        self.assertEqual(obs.hexdigest(), exp)

        fd.close()

    def test_remove_files(self):
        # create list of temp file paths
        test_fds = [NamedTemporaryFile(delete=False) for i in range(5)]
        test_filepaths = [element.name for element in test_fds]

        # should work just fine
        remove_files(test_filepaths)

        # check that an error is raised on trying to remove the files...
        self.assertRaises(OSError, remove_files, test_filepaths)

        # touch one of the filepaths so it exists
        extra_file = NamedTemporaryFile(delete=False).name
        test_filepaths.append(extra_file)

        # no error is raised on trying to remove the files
        # (although 5 don't exist)...
        remove_files(test_filepaths, error_on_missing=False)
        # ... and the existing file was removed
        self.assertFalse(exists(extra_file))

        # try to remove them with remove_files and verify that an IOError is
        # raises
        self.assertRaises(OSError, remove_files, test_filepaths)

        # now get no error when error_on_missing=False
        remove_files(test_filepaths, error_on_missing=False)

    def test_create_dir(self):
        # create a directory
        tmp_dir_path = mkdtemp()

        # create a random temporary directory name
        tmp_dir_path2 = join(mkdtemp(), str(uuid4()))
        tmp_dir_path3 = join(mkdtemp(), str(uuid4()))

        self.dirs_to_remove += [tmp_dir_path, tmp_dir_path2, tmp_dir_path3]

        # create on existing dir raises OSError if fail_on_exist=True
        self.assertRaises(OSError, create_dir, tmp_dir_path,
                          fail_on_exist=True)
        self.assertEqual(create_dir(tmp_dir_path, fail_on_exist=True,
                                    handle_errors_externally=True), 1)

        # return should be 1 if dir exist and fail_on_exist=False
        self.assertEqual(create_dir(tmp_dir_path, fail_on_exist=False), 1)

        # if dir not there make it and return always 0
        self.assertEqual(create_dir(tmp_dir_path2), 0)
        self.assertEqual(create_dir(tmp_dir_path3, fail_on_exist=True), 0)

    def test_handle_error_codes_no_error(self):
        obs = _handle_error_codes('/foo/bar/baz')
        self.assertEqual(obs, 0)

    def test_flatten(self):
        self.assertEqual(flatten(['aa', 'bb', 'cc']), list('aabbcc'))
        self.assertEqual(flatten([1, [2, 3], [[4, [5]]]]), [1, 2, 3, [4, [5]]])


class CardinalToOrdinalTests(unittest.TestCase):
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
        with self.assertRaisesRegexp(ValueError, '-1'):
            cardinal_to_ordinal(-1)


class TestFindDuplicates(unittest.TestCase):
    def test_empty_input(self):
        def empty_gen():
            raise StopIteration()
            yield

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
            for e in 'a', 1, 'bc', 2, 'a', 2, 2, 3.0:
                yield e

        self.assertEqual(find_duplicates(gen()), set(['a', 2]))


if __name__ == '__main__':
    unittest.main()
