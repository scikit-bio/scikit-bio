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
import os

import unittest
import warnings
from tempfile import mkstemp

from skbio.io import (DuplicateRegistrationError, FormatIdentificationWarning,
                      InvalidRegistrationError, UnrecognizedFormatError,
                      ArgumentOverrideWarning)
from skbio.io._registry import empty_file_sniffer
from skbio.util import TestingUtilError, get_data_path


class TestClass(object):
    def __init__(self, l):
        self.list = l

    def __eq__(self, other):
        # They are only equal when the class is EXACTLY the same. We don't want
        # readers to return knockoff instances...
        return self.__class__ is other.__class__ and self.list == other.list

    def __repr__(self):
        return "%s(%s)" % (str(self.__class__.__name__), str(self.list))


class TestClassA(TestClass):
    pass


class TestClassB(TestClass):
    pass


class RegistryTest(unittest.TestCase):
    def setUp(self):
        # A fresh module needs to be imported for each test because the
        # registry stores its state in the module which is by default
        # only loaded once.
        self.module = import_fresh_module('skbio.io._registry')
        self.fd1, self.fp1 = mkstemp()
        self.fd2, self.fp2 = mkstemp()

    def tearDown(self):
        os.remove(self.fp1)
        os.close(self.fd1)
        os.remove(self.fp2)
        os.close(self.fd2)


class TestRegisterAndGetReader(RegistryTest):
    def test_get_reader_no_match(self):
        self.assertEqual(None, self.module.get_reader('not_a_format',
                                                      TestClass))

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

    def test_get_reader_when_only_writer_exists(self):
        @self.module.register_writer('format', TestClass)
        def format_reader(fh):
            return

        self.assertEqual(None, self.module.get_reader('format', TestClass))

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

    def test_register_reader_generator_with_not_a_generator(self):
        @self.module.register_reader('format')
        def not_a_generator(fp):
            return 'oops'

        fh = StringIO()
        with self.assertRaises(InvalidRegistrationError):
            next(self.module.get_reader('format')(fh))
        fh.close()


class TestRegisterAndGetWriter(RegistryTest):
    def test_get_writer_no_match(self):
        self.assertEqual(None, self.module.get_writer('not_a_format',
                                                      TestClass))

    def test_get_writer_when_only_reader_exists(self):
        @self.module.register_reader('format', TestClass)
        def format_reader(fh):
            return

        self.assertEqual(None, self.module.get_writer('format', TestClass))

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

    def test_register_writer_over_existing_generator(self):
        with self.assertRaises(DuplicateRegistrationError) as cm:
            @self.module.register_writer('format1')
            def format1_writer(obj, fh):
                return

            @self.module.register_writer('format1')
            def duplicate_format1_writer(obj, fh):
                return

        self.assertTrue('format1' in str(cm.exception))
        self.assertTrue('writer' in str(cm.exception))
        self.assertTrue('generator' in str(cm.exception))


class TestRegisterAndGetSniffer(RegistryTest):
    def test_get_sniffer_no_match(self):
        self.assertEqual(None, self.module.get_sniffer('not_a_format'))

    def test_register_sniffer_on_many(self):
        fh1 = StringIO(u'1')
        fh2 = StringIO(u'2')
        fh3 = StringIO(u'3')

        @self.module.register_sniffer('format1')
        def format1_sniffer(fh):
            return '1' in fh.readline(), {}

        @self.module.register_sniffer('format2')
        def format2_sniffer(fh):
            return '2' in fh.readline(), {}

        @self.module.register_sniffer('format3')
        def format3_sniffer(fh):
            return '3' in fh.readline(), {}

        self.assertEqual(format1_sniffer,
                         self.module.get_sniffer('format1'))

        self.assertEqual(format2_sniffer,
                         self.module.get_sniffer('format2'))

        self.assertEqual(format3_sniffer,
                         self.module.get_sniffer('format3'))

    def test_register_sniffer_over_existing(self):
        with self.assertRaises(DuplicateRegistrationError) as cm:
            @self.module.register_sniffer('format1')
            def format1_sniffer(fh):
                return False, {}

            @self.module.register_sniffer('format1')
            def duplicate_format1_sniffer(fh):
                return False, {}

        self.assertTrue('format1' in str(cm.exception))

    def test_sniffer_warns_on_exception(self):
        @self.module.register_sniffer('format')
        def format_sniffer(fh):
            raise TestingUtilError("Sniffer will return False and warn.")

        fh = StringIO()
        sniffer = self.module.get_sniffer('format')
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            with self.assertRaises(FormatIdentificationWarning):
                sniffer(fh)

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("ignore")
            result, kwargs = sniffer(fh)
            self.assertFalse(result)
            self.assertEqual({}, kwargs)

        fh.close()


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

        @self.module.register_writer('format5', TestClassA)
        def format5_clsA(fh):
            return

        formats = self.module.list_read_formats(TestClassA)
        self.assertTrue('format1' in formats)
        self.assertTrue('format2' in formats)
        self.assertTrue('format3' in formats)
        self.assertTrue('format4' not in formats)
        self.assertTrue('format5' not in formats)


class TestListWriteFormats(RegistryTest):
    def test_no_write_formats(self):
        @self.module.register_writer('format1', TestClassA)
        def this_isnt_on_clsB(fh):
            return

        self.assertEqual([], self.module.list_write_formats(TestClassB))

    def test_one_write_format(self):
        @self.module.register_writer('format1', TestClass)
        def format1_cls(fh):
            return

        self.assertEqual(['format1'],
                         self.module.list_write_formats(TestClass))

    def test_many_write_formats(self):
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

        @self.module.register_reader('format5', TestClassA)
        def format5_clsA(fh):
            return

        formats = self.module.list_write_formats(TestClassA)

        self.assertTrue('format1' in formats)
        self.assertTrue('format2' in formats)
        self.assertTrue('format3' in formats)
        self.assertTrue('format4' not in formats)
        self.assertTrue('format5' not in formats)


class TestSniff(RegistryTest):
    def setUp(self):
        super(TestSniff, self).setUp()

        @self.module.register_sniffer('format1')
        def format1_sniffer(fh):
            return '1' in fh.readline(), {}

        @self.module.register_sniffer('format2')
        def format2_sniffer(fh):
            return '2' in fh.readline(), {}

        @self.module.register_sniffer('format3')
        def format3_sniffer(fh):
            return '3' in fh.readline(), {}

        @self.module.register_sniffer('format4')
        def format4_sniffer(fh):
            return '4' in fh.readline(), {}

        @self.module.register_reader('format3', TestClass)
        def reader3(fh):
            return

        @self.module.register_reader('format4', TestClass)
        def reader4(fh):
            return

    def test_no_matches(self):
        fh = StringIO(u"no matches here")
        fh2 = StringIO(u"no matches here")
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.module.sniff(fh)
        self.assertTrue(str(fh) in str(cm.exception))

        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.module.sniff(fh, cls=TestClass)

        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.module.sniff(fh, cls=TestClassB)

        fh.close()

    def test_one_match(self):
        fh = StringIO(u"contains a 3")
        self.assertEqual('format3', self.module.sniff(fh)[0])

    def test_many_matches(self):
        fh = StringIO(u"1234 will match all")
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.module.sniff(fh)
        self.assertTrue("format1" in str(cm.exception))
        self.assertTrue("format2" in str(cm.exception))
        self.assertTrue("format3" in str(cm.exception))
        self.assertTrue("format4" in str(cm.exception))
        fh.close()

    def test_no_matches_w_cls(self):
        fh = StringIO(u"no matches here")
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.module.sniff(fh, cls=TestClass)
        self.assertTrue(str(fh) in str(cm.exception))
        fh.close()

    def test_one_match_w_cls(self):
        fh = StringIO(u"contains a 3")
        self.assertEqual('format3',
                         self.module.sniff(fh, cls=TestClass)[0])

    def test_many_matches_w_cls(self):
        fh = StringIO(u"1234 will only format3 and format4 w/ class")
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.module.sniff(fh, cls=TestClass)
        self.assertTrue("format1" not in str(cm.exception))
        self.assertTrue("format2" not in str(cm.exception))
        # Only format3 and format4 have a definition for the provided class.
        self.assertTrue("format3" in str(cm.exception))
        self.assertTrue("format4" in str(cm.exception))
        fh.close()

    def test_that_mode_is_used(self):
        fp = self.fp1
        with open(fp, 'w') as fh:
            fh.write('@\n#\n')

        @self.module.register_sniffer('format')
        def sniffer(fh):
            self.assertEqual(self.expected_mode, fh.mode)
            return '@' in fh.readline(), {}

        self.expected_mode = 'U'
        self.module.sniff(fp)

        self.expected_mode = 'r'
        self.module.sniff(fp, mode='r')

    def test_position_not_mutated_real_file(self):
        @self.module.register_sniffer('format')
        def sniffer(fh):
            return True, {}

        with open(get_data_path('real_file')) as fh:
            fh.seek(2)
            self.module.sniff(fh)
            self.assertEqual('b\n', next(fh))

    def test_position_not_mutated_fileish(self):
        @self.module.register_sniffer('format')
        def sniffer(fh):
            return True, {}

        fh = StringIO(u'a\nb\nc\nd\n')
        fh.seek(2)
        self.module.sniff(fh)
        self.assertEqual('b\n', next(fh))


class TestRead(RegistryTest):
    def test_format_and_into_are_none(self):
        fh = StringIO()
        with self.assertRaises(ValueError):
            self.module.read(fh)

        fh.close()

    def test_format_is_none(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_sniffer('format')
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        instance = self.module.read(fh, into=TestClass)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)
        fh.close()

    def test_into_is_none(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_reader('format')
        def reader(fh):
            for value in [int(x) for x in fh.read().split('\n')]:
                yield value

        generator = self.module.read(fh, format='format')
        first_run = True
        for a, b in zip(generator, [1, 2, 3, 4]):
            if first_run:
                fh.seek(3)
                first_run = False
            self.assertEqual(a, b)
            self.assertEqual(3, fh.tell())
        fh.close()

    def test_into_is_none_real_file(self):
        fp = self.fp1
        with open(fp, 'w') as fh:
            fh.write('1\n2\n3\n4')

        self._test_fh = None

        @self.module.register_reader('format')
        def reader(fh):
            self._test_fh = fh
            for value in [int(x) for x in fh.read().split('\n')]:
                yield value

        generator = self.module.read(fp, format='format')
        for a, b in zip(generator, [1, 2, 3, 4]):
            self.assertEqual(a, b)
        self.assertTrue(self._test_fh.closed)

    def test_reader_does_not_exist(self):
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.module.read(None, format='not_a_format', into=TestClass)

        self.assertTrue(TestClass.__name__ in str(cm.exception))
        self.assertTrue('not_a_format' in str(cm.exception))

        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.module.read(None, format='not_a_format2')

        self.assertTrue('generator' in str(cm.exception))
        self.assertTrue('not_a_format2' in str(cm.exception))

    def test_reader_is_not_generator(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_sniffer('format')
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @self.module.register_reader('format')
        def reader(fh):
            # Not a generator!
            return TestClass([int(x) for x in fh.read().split('\n')])

        with self.assertRaises(InvalidRegistrationError):
            next(self.module.read(fh, format='format'))

        fh.close()

    def test_reader_empty_file(self):
        fh = StringIO()

        @self.module.register_sniffer('format')
        def sniffer(fh):
            return False, {}

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            return

        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.module.read(fh, into=TestClass)
        self.assertIn('<emptyfile>', str(cm.exception))

        fh.close()

    def test_reader_exists_with_verify_true(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_sniffer('format')
        def sniffer(fh):
            self.was_verified = True
            return '1' in fh.readline(), {}

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        self.was_verified = False
        instance = self.module.read(fh, format='format', into=TestClass,
                                    verify=True)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)
        self.assertTrue(self.was_verified)

        # Remove if read-context management is support in the future.
        fh.seek(0)

        self.was_verified = False
        instance = self.module.read(fh, format='format', into=TestClass)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)
        self.assertTrue(self.was_verified)

        fh.close()

    def test_warning_raised(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_sniffer('format')
        def sniffer(fh):
            self.was_verified = True
            return False, {}

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            with self.assertRaises(FormatIdentificationWarning):
                self.was_verified = False
                instance = self.module.read(fh, format='format',
                                            into=TestClass, verify=True)
                self.assertEqual(TestClass([1, 2, 3, 4]), instance)
                self.assertTrue(self.was_verified)

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            with self.assertRaises(FormatIdentificationWarning):
                self.was_verified = False
                instance = self.module.read(fh, format='format',
                                            into=TestClass)
                self.assertEqual(TestClass([1, 2, 3, 4]), instance)
                self.assertTrue(self.was_verified)

        fh.close()

    def test_reader_exists_with_verify_false(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_sniffer('format')
        def sniffer(fh):
            self.was_verified = True
            return '1' in fh.readline(), {}

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        self.was_verified = False
        instance = self.module.read(fh, format='format', into=TestClass,
                                    verify=False)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)
        self.assertFalse(self.was_verified)
        fh.close()

    def test_reader_exists_real_file(self):
        fp = self.fp1
        with open(fp, 'w') as fh:
            fh.write('1\n2\n3\n4')

        @self.module.register_sniffer('format')
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        instance = self.module.read(fp, format='format', into=TestClass)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)

    def test_read_kwargs_passed_generator(self):
        @self.module.register_sniffer('format')
        def sniffer(fh):
            return True, {'arg1': 15, 'arg2': 'abc'}

        @self.module.register_reader('format')
        def reader(fh, **kwargs):
            self.assertEqual(kwargs['arg1'], 15)
            self.assertEqual(kwargs['arg2'], 'abc')
            self.assertEqual(kwargs['arg3'], [1])
            yield

        next(self.module.read(StringIO(), format='format', arg3=[1]))

    def test_read_kwargs_passed_and_override(self):
        @self.module.register_sniffer('format')
        def sniffer(fh):
            return True, {'arg1': 15, 'arg2': 'abc', 'override': 30}

        @self.module.register_reader('format', TestClass)
        def reader(fh, **kwargs):
            self.assertEqual(kwargs['arg1'], 15)
            self.assertEqual(kwargs['arg2'], 'abc')
            self.assertEqual(kwargs['arg3'], [1])
            return

        self.module.read(StringIO(u'notempty'), into=TestClass, arg3=[1])

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            # Should raise no warning and thus no error.
            self.module.read(StringIO(u'notempty'), into=TestClass, arg3=[1],
                             override=30)
            # Should raise a warning and thus an error.
            with self.assertRaises(ArgumentOverrideWarning):
                self.module.read(StringIO(u'notempty'), into=TestClass,
                                 arg3=[1], override=100)

    def test_that_mode_is_used(self):
        fp = self.fp1
        with open(fp, 'w') as fh:
            fh.write('1\n2\n3\n4')

        @self.module.register_sniffer('format')
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            self.assertEqual(self.expected_mode, fh.mode)
            return TestClass([int(x) for x in fh.read().split('\n')])

        self.expected_mode = 'U'
        instance = self.module.read(fp, format='format', into=TestClass)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)

        self.expected_mode = 'r'
        instance = self.module.read(fp, format='format', into=TestClass,
                                    mode='r')
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)


class TestWrite(RegistryTest):
    def test_writer_does_not_exist(self):
        fh = StringIO()
        with self.assertRaises(UnrecognizedFormatError) as cm:
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
        fp = self.fp1

        @self.module.register_writer('format', TestClass)
        def writer(obj, fh):
            fh.write('\n'.join(obj.list))

        self.module.write(obj, format='format', into=fp)

        with open(fp, 'U') as fh:
            self.assertEqual("1\n2\n3\n4", fh.read())

    def test_writer_passed_kwargs(self):
        @self.module.register_reader('format')
        def reader(fh):
            yield

        @self.module.register_writer('format')
        def writer(obj, fh, **kwargs):
            self.assertEqual(kwargs['passed'], True)

        generator = self.module.get_reader('format')(None)
        self.module.write(generator, format='format',
                          into=StringIO(), passed=True)

    def test_that_mode_is_used(self):
        obj = TestClass(['1', '2', '3', '4'])
        fp = self.fp1

        @self.module.register_writer('format', TestClass)
        def writer(obj, fh):
            fh.write('\n'.join(obj.list))
            self.assertEqual(self.expected_mode, fh.mode)

        self.expected_mode = 'w'
        self.module.write(obj, format='format', into=fp)

        with open(fp, 'U') as fh:
            self.assertEqual("1\n2\n3\n4", fh.read())

        fp = self.fp2
        self.expected_mode = 'a'
        self.module.write(obj, format='format', into=fp, mode='a')

        with open(fp, 'U') as fh:
            self.assertEqual("1\n2\n3\n4", fh.read())


class TestInitializeOOPInterface(RegistryTest):
    def setUp(self):
        super(TestInitializeOOPInterface, self).setUp()

        class UnassumingClass(object):
            pass

        class ClassWithDefault(object):
            default_write_format = 'favfmt'

        self.unassuming_class = UnassumingClass
        self.class_with_default = ClassWithDefault

    def test_no_readers_writers(self):
        self.module.initialize_oop_interface()
        self.assertFalse(hasattr(self.unassuming_class, 'read'))
        self.assertFalse(hasattr(self.unassuming_class, 'write'))
        self.assertFalse(hasattr(self.class_with_default, 'read'))
        self.assertFalse(hasattr(self.class_with_default, 'write'))

    def test_readers_only(self):
        @self.module.register_reader('favfmt', self.unassuming_class)
        def fvfmt_to_unasumming_class(fh):
            return

        @self.module.register_reader('favfmt')
        def fvfmt_to_gen(fh):
            yield

        @self.module.register_reader('favfmt2', self.unassuming_class)
        def fvfmt2_to_unasumming_class(fh):
            return

        self.module.initialize_oop_interface()

        self.assertTrue(hasattr(self.unassuming_class, 'read'))
        self.assertFalse(hasattr(self.unassuming_class, 'write'))
        self.assertFalse(hasattr(self.class_with_default, 'read'))
        self.assertFalse(hasattr(self.class_with_default, 'write'))

        self.assertIn('favfmt', self.unassuming_class.read.__doc__)
        self.assertIn('favfmt2', self.unassuming_class.read.__doc__)

    def test_writers_only(self):
        @self.module.register_writer('favfmt', self.class_with_default)
        def favfmt(fh):
            pass

        @self.module.register_writer('favfmt')
        def gen_to_favfmt(fh):
            pass

        @self.module.register_writer('favfmt2', self.class_with_default)
        def favfmt2(fh):
            pass

        self.module.initialize_oop_interface()

        self.assertFalse(hasattr(self.unassuming_class, 'read'))
        self.assertFalse(hasattr(self.unassuming_class, 'write'))
        self.assertFalse(hasattr(self.class_with_default, 'read'))
        self.assertTrue(hasattr(self.class_with_default, 'write'))

        self.assertIn('favfmt', self.class_with_default.write.__doc__)
        self.assertIn('favfmt2', self.class_with_default.write.__doc__)

    def test_writers_no_default_format(self):
        @self.module.register_writer('favfmt', self.unassuming_class)
        def favfmt(fh):
            pass

        @self.module.register_writer('favfmt')
        def gen_to_favfmt(fh):
            pass

        @self.module.register_writer('favfmt2', self.unassuming_class)
        def favfmt2(fh):
            pass
        with self.assertRaises(NotImplementedError) as cm:
            self.module.initialize_oop_interface()

        self.assertIn('default_write_format', str(cm.exception))

    def test_readers_writers(self):
        @self.module.register_reader('favfmt', self.unassuming_class)
        def fvfmt_to_unasumming_class(fh):
            return

        @self.module.register_reader('favfmt', self.class_with_default)
        def fvfmt_to_class_w_default(fh):
            return

        @self.module.register_reader('favfmt')
        def fvfmt_to_gen(fh):
            yield

        @self.module.register_reader('favfmt2', self.unassuming_class)
        def fvfmt2_to_unasumming_class(fh):
            return

        @self.module.register_reader('favfmt2', self.class_with_default)
        def fvfmt2_to_class_w_default(fh):
            return

        @self.module.register_writer('favfmt', self.class_with_default)
        def favfmt(fh):
            pass

        @self.module.register_writer('favfmt')
        def gen_to_favfmt(fh):
            pass

        @self.module.register_writer('favfmt2', self.class_with_default)
        def favfmt2(fh):
            pass

        self.module.initialize_oop_interface()

        self.assertTrue(hasattr(self.unassuming_class, 'read'))
        self.assertFalse(hasattr(self.unassuming_class, 'write'))

        self.assertTrue(hasattr(self.class_with_default, 'read'))
        self.assertTrue(hasattr(self.class_with_default, 'write'))

        self.assertIn('favfmt', self.unassuming_class.read.__doc__)
        self.assertIn('favfmt2', self.unassuming_class.read.__doc__)

        self.assertIn('favfmt', self.class_with_default.read.__doc__)
        self.assertIn('favfmt2', self.class_with_default.read.__doc__)

        self.assertIn('favfmt', self.class_with_default.write.__doc__)
        self.assertIn('favfmt2', self.class_with_default.write.__doc__)

    def test_read_kwargs_passed(self):
        self.was_called = False

        @self.module.register_sniffer('favfmt')
        def fvfmt_sniffer(fh):
            return True, {}

        @self.module.register_reader('favfmt', self.class_with_default)
        def fvfmt_to_class_w_default(fh, **kwargs):
            self.assertEqual('a', kwargs['a'])
            self.assertEqual(123, kwargs['b'])
            self.was_called = True

        self.module.initialize_oop_interface()
        fh = StringIO(u'notempty')
        self.class_with_default.read(fh, a='a', b=123)

        self.assertTrue(self.was_called)
        fh.close()

    def test_write_kwargs_passed(self):
        self.was_called = False

        @self.module.register_writer('favfmt', self.class_with_default)
        def favfmt(obj, fh, **kwargs):
            self.assertEqual('a', kwargs['a'])
            self.assertEqual(123, kwargs['b'])
            self.was_called = True

        self.module.initialize_oop_interface()
        fh = StringIO()
        self.class_with_default().write(fh, a='a', b=123)

        self.assertTrue(self.was_called)
        fh.close()


class TestEmptyFileSniffer(unittest.TestCase):
    def test_blank_file(self):
        fh = StringIO()
        self.assertTrue(empty_file_sniffer(fh)[0])
        fh.close()

    def test_whitespace_file(self):
        fh = StringIO(u' ')
        self.assertTrue(empty_file_sniffer(fh)[0])
        fh.close()
        fh = StringIO(u'\n')
        self.assertTrue(empty_file_sniffer(fh)[0])
        fh.close()
        fh = StringIO(u'\t')
        self.assertTrue(empty_file_sniffer(fh)[0])
        fh.close()

    def test_mixed_whitespace_file(self):
        fh = StringIO(u'\n\n\t\n \t \t \n \n \n\n')
        self.assertTrue(empty_file_sniffer(fh)[0])
        fh.close()

    def test_not_empty_file(self):
        fh = StringIO(u'\n\n\t\n a\t \t \n \n \n\n')
        self.assertFalse(empty_file_sniffer(fh)[0])
        fh.close()

if __name__ == '__main__':
    unittest.main()
