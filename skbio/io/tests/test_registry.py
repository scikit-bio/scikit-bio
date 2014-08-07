#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
try:
    # future >= 0.12
    from future.backports.test.support import import_fresh_module
except ImportError:
    from future.standard_library.test.support import import_fresh_module
from io import StringIO

import unittest
from tempfile import mkstemp

from skbio.io import (DuplicateRegistrationError,
                      FormatIdentificationError, FileFormatError)


class TestClass(object):
    def __init__(self, l):
        self.list = l


class TestClassA(TestClass):
    pass


class TestClassB(TestClass):
    pass


class RegistryTest(unittest.TestCase):
    def setUp(self):
        self.module = import_fresh_module('skbio.io._registry')


class TestRegisterAndReader(RegistryTest):
    def test_get_reader_no_match(self):
        self.assertEqual(None, self.module.get_reader('not_a_format',
                                                      TestClass))

    def test_get_reader_too_many_args(self):
        with self.assertRaises(TypeError) as cm:
            self.module.get_reader('not_a_format', TestClass, 1)
        self.assertTrue('get_reader' in str(cm.exception))

    def test_register_reader_on_generator(self):
        @self.module.register_reader('format1')
        def format1_reader_generator(fh):
            yield

        self.assertEqual(format1_reader_generator,
                         self.module.get_reader('format1'))

        self.assertEqual(format1_reader_generator,
                         self.module.get_reader('format1', None))

        @self.module.register_reader('format2', None)
        def format2_reader_generator(fh):
            yield

        self.assertEqual(format2_reader_generator,
                         self.module.get_reader('format2'))

        self.assertEqual(format2_reader_generator,
                         self.module.get_reader('format2', None))

    def test_register_reader_with_too_many_args(self):
        with self.assertRaises(TypeError) as cm:
            @self.module.register_reader('format1', TestClass, 1)
            def too_many_args(fh):
                return
        self.assertTrue('register_reader' in str(cm.exception))

    def test_register_reader_on_many(self):
        @self.module.register_reader('format1', TestClassA)
        def format1_reader(fh):
            return

        @self.module.register_reader('format1', TestClassB)
        def format1_reader_b(fh):
            return

        @self.module.register_reader('format2', TestClassA)
        def format2_reader(fh):
            return

        @self.module.register_reader('format3', TestClassB)
        def format3_reader(fh):
            return

        self.assertEqual(format1_reader,
                         self.module.get_reader('format1', TestClassA))

        self.assertEqual(format1_reader_b,
                         self.module.get_reader('format1', TestClassB))

        self.assertEqual(format2_reader,
                         self.module.get_reader('format2', TestClassA))

        self.assertEqual(None,
                         self.module.get_reader('format2', TestClassB))

        self.assertEqual(None,
                         self.module.get_reader('format3', TestClassA))

        self.assertEqual(format3_reader,
                         self.module.get_reader('format3', TestClassB))

    def test_register_reader_over_existing(self):
        with self.assertRaises(DuplicateRegistrationError) as cm:
            @self.module.register_reader('format1', TestClassA)
            def format1_reader(fh):
                return

            @self.module.register_reader('format1', TestClassA)
            def duplicate_format1_reader(fh):
                return

        self.assertTrue('format1' in str(cm.exception))
        self.assertTrue('reader' in str(cm.exception))
        self.assertTrue(TestClassA.__name__ in str(cm.exception))


class TestRegisterAndGetWriter(RegistryTest):
    def test_get_writer_no_match(self):
        self.assertEqual(None, self.module.get_writer('not_a_format',
                                                      TestClass))

    def test_get_writer_too_many_args(self):
        with self.assertRaises(TypeError) as cm:
            self.module.get_writer('not_a_format', TestClass, 1)
        self.assertTrue('get_writer' in str(cm.exception))

    def test_register_writer_on_generator(self):
        @self.module.register_writer('format1')
        def format1_writer_generator(obj, fh):
            yield

        self.assertEqual(format1_writer_generator,
                         self.module.get_writer('format1'))

        self.assertEqual(format1_writer_generator,
                         self.module.get_writer('format1', None))

        @self.module.register_writer('format2', None)
        def format2_writer_generator(obj, fh):
            yield

        self.assertEqual(format2_writer_generator,
                         self.module.get_writer('format2'))

        self.assertEqual(format2_writer_generator,
                         self.module.get_writer('format2', None))

    def test_register_writer_with_too_many_args(self):
        with self.assertRaises(TypeError) as cm:
            @self.module.register_writer('format1', TestClass, 1)
            def too_many_args(obj, fh):
                return
        self.assertTrue('register_writer' in str(cm.exception))

    def test_register_writer_on_many(self):
        @self.module.register_writer('format1', TestClassA)
        def format1_writer(obj, fh):
            return

        @self.module.register_writer('format1', TestClassB)
        def format1_writer_b(obj, fh):
            return

        @self.module.register_writer('format2', TestClassA)
        def format2_writer(obj, fh):
            return

        @self.module.register_writer('format3', TestClassB)
        def format3_writer(obj, fh):
            return

        self.assertEqual(format1_writer,
                         self.module.get_writer('format1', TestClassA))

        self.assertEqual(format1_writer_b,
                         self.module.get_writer('format1', TestClassB))

        self.assertEqual(format2_writer,
                         self.module.get_writer('format2', TestClassA))

        self.assertEqual(None,
                         self.module.get_writer('format2', TestClassB))

        self.assertEqual(None,
                         self.module.get_writer('format3', TestClassA))

        self.assertEqual(format3_writer,
                         self.module.get_writer('format3', TestClassB))

    def test_register_writer_over_existing(self):
        with self.assertRaises(DuplicateRegistrationError) as cm:
            @self.module.register_writer('format1', TestClassA)
            def format1_writer(obj, fh):
                return

            @self.module.register_writer('format1', TestClassA)
            def duplicate_format1_writer(obj, fh):
                return

        self.assertTrue('format1' in str(cm.exception))
        self.assertTrue('writer' in str(cm.exception))
        self.assertTrue(TestClassA.__name__ in str(cm.exception))


class TestRegisterAndGetIdentifer(RegistryTest):
    def test_get_identifier_no_match(self):
        self.assertEqual(None, self.module.get_identifier('not_a_format'))

    def test_register_identifier_on_many(self):
        @self.module.register_identifier('format1')
        def format1_identifier(fh):
            return

        @self.module.register_identifier('format2')
        def format2_identifier(fh):
            return

        @self.module.register_identifier('format3')
        def format3_identifier(fh):
            return

        self.assertEqual(format1_identifier,
                         self.module.get_identifier('format1'))

        self.assertEqual(format2_identifier,
                         self.module.get_identifier('format2'))

        self.assertEqual(format3_identifier,
                         self.module.get_identifier('format3'))

    def test_register_identifier_over_existing(self):
        with self.assertRaises(DuplicateRegistrationError) as cm:
            @self.module.register_identifier('format1')
            def format1_identifier(fh):
                return

            @self.module.register_identifier('format1')
            def duplicate_format1_identifier(fh):
                return

        self.assertTrue('format1' in str(cm.exception))


class TestListReadFormats(RegistryTest):
    def test_no_read_formats(self):
        @self.module.register_reader('format1', TestClassA)
        def this_isnt_on_clsB(fh):
            return

        self.assertEqual([], self.module.list_read_formats(TestClassB))

    def test_one_read_format(self):
        @self.module.register_reader('format1', TestClass)
        def format1_cls(fh):
            return

        self.assertEqual(['format1'], self.module.list_read_formats(TestClass))

    def test_many_read_formats(self):
        @self.module.register_reader('format1', TestClassA)
        def format1_clsA(fh):
            return

        @self.module.register_reader('format2', TestClassA)
        def format2_clsA(fh):
            return

        @self.module.register_reader('format3', TestClassA)
        def format3_clsA(fh):
            return

        @self.module.register_reader('format3', TestClassB)
        def format3_clsB(fh):
            return

        @self.module.register_reader('format4', TestClassB)
        def format4_clsB(fh):
            return

        formats = self.module.list_read_formats(TestClassA)

        self.assertTrue('format1' in formats)
        self.assertTrue('format2' in formats)
        self.assertTrue('format3' in formats)
        self.assertTrue('format4' not in formats)


class TestListWriteFormats(RegistryTest):
    def test_no_read_formats(self):
        @self.module.register_writer('format1', TestClassA)
        def this_isnt_on_clsB(fh):
            return

        self.assertEqual([], self.module.list_write_formats(TestClassB))

    def test_one_read_format(self):
        @self.module.register_writer('format1', TestClass)
        def format1_cls(fh):
            return

        self.assertEqual(['format1'],
                         self.module.list_write_formats(TestClass))

    def test_many_read_formats(self):
        @self.module.register_writer('format1', TestClassA)
        def format1_clsA(fh):
            return

        @self.module.register_writer('format2', TestClassA)
        def format2_clsA(fh):
            return

        @self.module.register_writer('format3', TestClassA)
        def format3_clsA(fh):
            return

        @self.module.register_writer('format3', TestClassB)
        def format3_clsB(fh):
            return

        @self.module.register_writer('format4', TestClassB)
        def format4_clsB(fh):
            return

        formats = self.module.list_write_formats(TestClassA)

        self.assertTrue('format1' in formats)
        self.assertTrue('format2' in formats)
        self.assertTrue('format3' in formats)
        self.assertTrue('format4' not in formats)


class TestGuessFormat(RegistryTest):
    def setUp(self):
        super(TestGuessFormat, self).setUp()

        @self.module.register_identifier('format1')
        def format1_identifier(fh):
            return '1' in fh.readline()

        @self.module.register_identifier('format2')
        def format2_identifier(fh):
            return '2' in fh.readline()

        @self.module.register_identifier('format3')
        def format3_identifier(fh):
            return '3' in fh.readline()

        @self.module.register_identifier('format4')
        def format4_identifier(fh):
            return '4' in fh.readline()

        @self.module.register_reader('format3', TestClass)
        def reader3(fh):
            return

        @self.module.register_reader('format4', TestClass)
        def reader4(fh):
            return

    def test_no_matches(self):
        fh = StringIO(u"no matches here")
        with self.assertRaises(FormatIdentificationError) as cm:
            self.module.guess_format(fh)
        self.assertTrue(str(fh) in str(cm.exception))

        with self.assertRaises(FormatIdentificationError) as cm:
            self.module.guess_format(fh, cls=TestClass)

        with self.assertRaises(FormatIdentificationError) as cm:
            self.module.guess_format(fh, cls=TestClassB)

        fh.close()

    def test_one_match(self):
        fh = StringIO(u"contains a 3")
        self.assertEqual('format3', self.module.guess_format(fh))

    def test_many_matches(self):
        fh = StringIO(u"1234 will match all")
        with self.assertRaises(FormatIdentificationError) as cm:
            self.module.guess_format(fh)
        self.assertTrue("format1" in str(cm.exception))
        self.assertTrue("format2" in str(cm.exception))
        self.assertTrue("format3" in str(cm.exception))
        self.assertTrue("format4" in str(cm.exception))
        fh.close()

    def test_no_matches_w_cls(self):
        fh = StringIO(u"no matches here")
        with self.assertRaises(FormatIdentificationError) as cm:
            self.module.guess_format(fh, cls=TestClass)
        self.assertTrue(str(fh) in str(cm.exception))
        fh.close()

    def test_one_match_w_cls(self):
        fh = StringIO(u"contains a 3")
        self.assertEqual('format3',
                         self.module.guess_format(fh, cls=TestClass))

    def test_many_matches_w_cls(self):
        fh = StringIO(u"1234 will only format3 and format4 w/ class")
        with self.assertRaises(FormatIdentificationError) as cm:
            self.module.guess_format(fh, cls=TestClass)
        self.assertTrue("format1" not in str(cm.exception))
        self.assertTrue("format2" not in str(cm.exception))
        # Only format3 and format4 have a definition for the provided class.
        self.assertTrue("format3" in str(cm.exception))
        self.assertTrue("format4" in str(cm.exception))
        fh.close()


class TestRead(RegistryTest):
    def test_format_and_into_are_none(self):
        fh = StringIO()
        with self.assertRaises(ValueError):
            self.module.read(fh)

        fh.close()

    def test_format_is_none(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_identifier('format')
        def identifier(fh):
            return '1' in fh.readline()

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        instance = self.module.read(fh, into=TestClass)
        self.assertEqual([1, 2, 3, 4], instance.list)
        fh.close()

    def test_into_is_none(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_reader('format')
        def reader(fh):
            for value in [int(x) for x in fh.read().split('\n')]:
                yield value

        generator = self.module.read(fh, format='format')
        for a, b in zip([1, 2, 3, 4], generator):
            self.assertEqual(a, b)
        self.assertTrue(not fh.closed)
        fh.close()

    def test_into_is_none_real_file(self):
        _, fp = mkstemp()
        with open(fp, 'w') as fh:
            fh.write('1\n2\n3\n4')

        @self.module.register_reader('format')
        def reader(fh):
            for value in [int(x) for x in fh.read().split('\n')]:
                yield value

        generator = self.module.read(fp, format='format')
        for a, b in zip([1, 2, 3, 4], generator):
            self.assertEqual(a, b)
        self.assertTrue(fh.closed)

    def test_reader_does_not_exist(self):
        with self.assertRaises(FileFormatError) as cm:
            self.module.read(None, format='not_a_format', into=TestClass)

        self.assertTrue(TestClass.__name__ in str(cm.exception))
        self.assertTrue('not_a_format' in str(cm.exception))

        with self.assertRaises(FileFormatError) as cm:
            self.module.read(None, format='not_a_format2')

        self.assertTrue('generator' in str(cm.exception))
        self.assertTrue('not_a_format2' in str(cm.exception))

    def test_reader_exists(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_identifier('format')
        def identifier(fh):
            return '1' in fh.readline()

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        instance = self.module.read(fh, format='format', into=TestClass)
        self.assertEqual([1, 2, 3, 4], instance.list)
        fh.close()

    def test_reader_exists_real_file(self):
        _, fp = mkstemp()
        with open(fp, 'w') as fh:
            fh.write('1\n2\n3\n4')

        @self.module.register_identifier('format')
        def identifier(fh):
            return '1' in fh.readline()

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        instance = self.module.read(fp, format='format', into=TestClass)
        self.assertEqual([1, 2, 3, 4], instance.list)


class TestWrite(RegistryTest):
    def test_format_is_none(self):
        fh = StringIO()
        with self.assertRaises(ValueError):
            self.module.write({}, into=fh)
        fh.close()

    def test_into_is_none(self):
        fh = StringIO()
        with self.assertRaises(ValueError):
            self.module.write({}, format='format')
        fh.close()

    def test_writer_does_not_exist(self):
        fh = StringIO()
        with self.assertRaises(FileFormatError) as cm:
            self.module.write({}, format='not_a_format', into=fh)

        self.assertTrue('not_a_format' in str(cm.exception))
        self.assertTrue(str(fh) in str(cm.exception))
        fh.close()

    def test_writer_exists(self):
        obj = TestClass(['1', '2', '3', '4'])
        fh = StringIO()

        @self.module.register_writer('format', TestClass)
        def writer(obj, fh):
            fh.write(u'\n'.join(obj.list))

        self.module.write(obj, format='format', into=fh)
        fh.seek(0)
        self.assertEqual("1\n2\n3\n4", fh.read())
        fh.close()

    def test_writer_exists_real_file(self):
        obj = TestClass(['1', '2', '3', '4'])
        _, fp = mkstemp()

        @self.module.register_writer('format', TestClass)
        def writer(obj, fh):
            fh.write('\n'.join(obj.list))

        self.module.write(obj, format='format', into=fp)

        with open(fp, 'U') as fh:
            self.assertEqual("1\n2\n3\n4", fh.read())

if __name__ == '__main__':
    unittest.main()
