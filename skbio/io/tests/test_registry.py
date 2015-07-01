# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from six.moves import zip_longest

from io import StringIO
import io
import os
import unittest
import warnings
from tempfile import mkstemp

from skbio.io import (DuplicateRegistrationError, FormatIdentificationWarning,
                      UnrecognizedFormatError, ArgumentOverrideWarning)
from skbio.io.registry import IORegistry, FileSentinel
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
        self.registry = IORegistry()
        self.fd1, self.fp1 = mkstemp()
        self.fd2, self.fp2 = mkstemp()

    def tearDown(self):
        os.remove(self.fp1)
        os.close(self.fd1)
        os.remove(self.fp2)
        os.close(self.fd2)


class TestRegisterAndGetReader(RegistryTest):
    def test_get_reader_no_match(self):
        self.assertIs(None, self.registry.get_reader('not_a_format',
                                                     TestClass))

    def test_get_reader_when_only_writer_exists(self):
        format1 = self.registry.create_format('format1')

        @format1.writer(TestClass)
        def format_writer(fh):
            return

        self.assertEqual(None, self.registry.get_reader('format', TestClass))

    def test_register_reader_on_many(self):
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3')

        @format1.reader(TestClassA)
        def format1_reader(fh):
            return

        @format1.reader(TestClassB)
        def format1_reader_b(fh):
            return

        @format2.reader(TestClassA)
        def format2_reader(fh):
            return

        @format3.reader(TestClassB)
        def format3_reader(fh):
            return

        self.assertIs(format1_reader,
                      self.registry.get_reader('format1', TestClassA))

        self.assertIs(format1_reader_b,
                      self.registry.get_reader('format1', TestClassB))

        self.assertIs(format2_reader,
                      self.registry.get_reader('format2', TestClassA))

        self.assertIs(None, self.registry.get_reader('format2', TestClassB))

        self.assertIs(None, self.registry.get_reader('format3', TestClassA))

        self.assertIs(format3_reader,
                      self.registry.get_reader('format3', TestClassB))

    def test_register_reader_over_existing(self):
        format1 = self.registry.create_format('format1')
        with self.assertRaises(DuplicateRegistrationError) as cm:
            @format1.reader(TestClassA)
            def format1_reader(fh):
                return

            @format1.reader(TestClassA)
            def duplicate_format1_reader(fh):
                return

        self.assertTrue('format1' in str(cm.exception))
        self.assertTrue('reader' in str(cm.exception))
        self.assertTrue(TestClassA.__name__ in str(cm.exception))

    def test_register_reader_over_existing_override(self):
        format1 = self.registry.create_format('format1')

        @format1.reader(TestClassA)
        def format1_reader(fh):
            return

        self.assertIs(format1_reader,
                      self.registry.get_reader('format1', TestClassA))

        @format1.reader(TestClassA, override=True)
        def duplicate_format1_reader(fh):
            return

        self.assertIs(duplicate_format1_reader,
                      self.registry.get_reader('format1', TestClassA))


class TestRegisterAndGetWriter(RegistryTest):
    def test_get_writer_no_match(self):
        self.assertEqual(None, self.registry.get_writer('not_a_format',
                                                        TestClass))

    def test_get_writer_when_only_reader_exists(self):
        format = self.registry.create_format('format')

        @format.reader(TestClass)
        def format_reader(fh):
            return

        self.assertEqual(None, self.registry.get_writer('format', TestClass))

    def test_register_writer_on_many(self):
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3')

        @format1.writer(TestClassA)
        def format1_writer(obj, fh):
            return

        @format1.writer(TestClassB)
        def format1_writer_b(obj, fh):
            return

        @format2.writer(TestClassA)
        def format2_writer(obj, fh):
            return

        @format3.writer(TestClassB)
        def format3_writer(obj, fh):
            return

        self.assertEqual(format1_writer,
                         self.registry.get_writer('format1', TestClassA))

        self.assertEqual(format1_writer_b,
                         self.registry.get_writer('format1', TestClassB))

        self.assertEqual(format2_writer,
                         self.registry.get_writer('format2', TestClassA))

        self.assertEqual(None,
                         self.registry.get_writer('format2', TestClassB))

        self.assertEqual(None,
                         self.registry.get_writer('format3', TestClassA))

        self.assertEqual(format3_writer,
                         self.registry.get_writer('format3', TestClassB))

    def test_register_writer_over_existing(self):
        format1 = self.registry.create_format('format1')
        with self.assertRaises(DuplicateRegistrationError) as cm:
            @format1.writer(TestClassA)
            def format1_writer(obj, fh):
                return

            @format1.writer(TestClassA)
            def duplicate_format1_writer(obj, fh):
                return

        self.assertTrue('format1' in str(cm.exception))
        self.assertTrue('writer' in str(cm.exception))
        self.assertTrue(TestClassA.__name__ in str(cm.exception))

    def test_register_writer_over_existing_override(self):
        format1 = self.registry.create_format('format1')

        @format1.writer(TestClassA)
        def format1_writer(obj, fh):
            return

        self.assertIs(format1_writer,
                      self.registry.get_writer('format1', TestClassA))

        @format1.writer(TestClassA, override=True)
        def duplicate_format1_writer(obj, fh):
            return

        self.assertIs(duplicate_format1_writer,
                      self.registry.get_writer('format1', TestClassA))


class TestRegisterAndGetSniffer(RegistryTest):
    def test_get_sniffer_no_match(self):
        self.assertEqual(None, self.registry.get_sniffer('not_a_format'))

    def test_register_sniffer_on_many(self):
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3')

        @format1.sniffer()
        def format1_sniffer(fh):
            return '1' in fh.readline(), {}

        @format2.sniffer()
        def format2_sniffer(fh):
            return '2' in fh.readline(), {}

        @format3.sniffer()
        def format3_sniffer(fh):
            return '3' in fh.readline(), {}

        self.assertEqual(format1_sniffer,
                         self.registry.get_sniffer('format1'))

        self.assertEqual(format2_sniffer,
                         self.registry.get_sniffer('format2'))

        self.assertEqual(format3_sniffer,
                         self.registry.get_sniffer('format3'))

    def test_register_sniffer_over_existing(self):
        format1 = self.registry.create_format('format1')

        with self.assertRaises(DuplicateRegistrationError) as cm:
            @format1.sniffer()
            def format1_sniffer(fh):
                return False, {}

            @format1.sniffer()
            def duplicate_format1_sniffer(fh):
                return False, {}

        self.assertTrue('format1' in str(cm.exception))

    def test_register_sniffer_over_existing_override(self):
        format1 = self.registry.create_format('format1')

        @format1.sniffer()
        def format1_sniffer(fh):
            return False, {}

        self.assertIs(self.registry.get_sniffer('format1'), format1_sniffer)

        @format1.sniffer(override=True)
        def duplicate_format1_sniffer(fh):
            return False, {}

        self.assertIs(self.registry.get_sniffer('format1'),
                      duplicate_format1_sniffer)

    def test_sniffer_warns_on_exception(self):
        format = self.registry.create_format('format')

        @format.sniffer()
        def format_sniffer(fh):
            raise TestingUtilError("Sniffer will return False and warn.")

        fh = StringIO()
        sniffer = self.registry.get_sniffer('format')
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
        format1 = self.registry.create_format('format1')

        @format1.reader(TestClassA)
        def this_isnt_on_clsB(fh):
            return

        self.assertEqual([], self.registry.list_read_formats(TestClassB))

    def test_one_read_format(self):
        format1 = self.registry.create_format('format1')

        @format1.reader(TestClass)
        def format1_cls(fh):
            return

        self.assertEqual(['format1'],
                         self.registry.list_read_formats(TestClass))

    def test_many_read_formats(self):
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3')
        format4 = self.registry.create_format('format4')
        format5 = self.registry.create_format('format5')

        @format1.reader(TestClassA)
        def format1_clsA(fh):
            return

        @format2.reader(TestClassA)
        def format2_clsA(fh):
            return

        @format3.reader(TestClassA)
        def format3_clsA(fh):
            return

        @format3.reader(TestClassB)
        def format3_clsB(fh):
            return

        @format4.reader(TestClassB)
        def format4_clsB(fh):
            return

        @format5.writer(TestClassA)
        def format5_clsA(fh):
            return

        formats = self.registry.list_read_formats(TestClassA)
        self.assertTrue('format1' in formats)
        self.assertTrue('format2' in formats)
        self.assertTrue('format3' in formats)
        self.assertTrue('format4' not in formats)
        self.assertTrue('format5' not in formats)


class TestListWriteFormats(RegistryTest):
    def test_no_write_formats(self):
        format1 = self.registry.create_format('format1')

        @format1.writer(TestClassA)
        def this_isnt_on_clsB(fh):
            return

        self.assertEqual([], self.registry.list_write_formats(TestClassB))

    def test_one_write_format(self):
        format1 = self.registry.create_format('format1')

        @format1.writer(TestClass)
        def format1_cls(fh):
            return

        self.assertEqual(['format1'],
                         self.registry.list_write_formats(TestClass))

    def test_many_write_formats(self):
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3')
        format4 = self.registry.create_format('format4')
        format5 = self.registry.create_format('format5')

        @format1.writer(TestClassA)
        def format1_clsA(fh):
            return

        @format2.writer(TestClassA)
        def format2_clsA(fh):
            return

        @format3.writer(TestClassA)
        def format3_clsA(fh):
            return

        @format3.writer(TestClassB)
        def format3_clsB(fh):
            return

        @format4.writer(TestClassB)
        def format4_clsB(fh):
            return

        @format5.reader(TestClassA)
        def format5_clsA(fh):
            return

        formats = self.registry.list_write_formats(TestClassA)

        self.assertTrue('format1' in formats)
        self.assertTrue('format2' in formats)
        self.assertTrue('format3' in formats)
        self.assertTrue('format4' not in formats)
        self.assertTrue('format5' not in formats)


class TestSniff(RegistryTest):
    def setUp(self):
        super(TestSniff, self).setUp()
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3')
        format4 = self.registry.create_format('format4')

        @format1.sniffer()
        def format1_sniffer(fh):
            return '1' in fh.readline(), {}

        @format2.sniffer()
        def format2_sniffer(fh):
            return '2' in fh.readline(), {}

        @format3.sniffer()
        def format3_sniffer(fh):
            return '3' in fh.readline(), {}

        @format4.sniffer()
        def format4_sniffer(fh):
            return '4' in fh.readline(), {}

        @format3.reader(TestClass)
        def reader3(fh):
            return

        @format4.reader(TestClass)
        def reader4(fh):
            return

    def test_no_matches(self):
        fh = StringIO(u"no matches here")
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.registry.sniff(fh)
        self.assertTrue(str(fh) in str(cm.exception))

        fh.close()

    def test_one_match(self):
        fh = StringIO(u"contains a 3")
        self.assertEqual('format3', self.registry.sniff(fh)[0])

    def test_many_matches(self):
        fh = StringIO(u"1234 will match all")
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.registry.sniff(fh)
        self.assertTrue("format1" in str(cm.exception))
        self.assertTrue("format2" in str(cm.exception))
        self.assertTrue("format3" in str(cm.exception))
        self.assertTrue("format4" in str(cm.exception))
        fh.close()

#    def test_that_mode_is_used(self):
#        fp = self.fp1
#        with open(fp, 'w') as fh:
#            fh.write('@\n#\n')
#
#        @self.registry.register_sniffer('format')
#        def sniffer(fh):
#            self.assertEqual(self.expected_mode, fh.mode)
#            return '@' in fh.readline(), {}
#
#        self.expected_mode = 'U'
#        self.registry.sniff(fp)
#
#        self.expected_mode = 'r'
#        self.registry.sniff(fp, mode='r')

    def test_position_not_mutated_real_file(self):
        formatx = self.registry.create_format('formatx')

        @formatx.sniffer()
        def sniffer(fh):
            return True, {}

        with io.open(get_data_path('real_file')) as fh:
            fh.seek(2)
            self.registry.sniff(fh)
            self.assertEqual(fh.tell(), 2)
            self.assertEqual('b\n', fh.readline())

    def test_position_not_mutated_fileish(self):
        formatx = self.registry.create_format('formatx')

        @formatx.sniffer()
        def sniffer(fh):
            return True, {}

        fh = StringIO(u'a\nb\nc\nd\n')
        fh.seek(2)
        self.registry.sniff(fh)
        self.assertEqual('b\n', fh.readline())


class TestRead(RegistryTest):
    def test_format_and_into_are_none(self):
        fh = StringIO()
        with self.assertRaises(ValueError):
            self.registry.read(fh)

        fh.close()

    def test_format_is_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO(u'1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        instance = self.registry.read(fh, into=TestClass)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)
        fh.close()

    def test_into_is_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO(u'1\n2\n3\n4')

        @format1.reader(None)
        def reader(fh):
            for value in [int(x) for x in fh.read().split('\n')]:
                yield value

        generator = self.registry.read(fh, format='format1')
        first_run = True
        for a, b in zip(generator, [1, 2, 3, 4]):
            if first_run:
                fh.seek(3)
                first_run = False
            self.assertEqual(a, b)
            self.assertEqual(3, fh.tell())
        fh.close()

    def test_into_is_none_real_file(self):
        format1 = self.registry.create_format('format1')

        fp = self.fp1
        with open(fp, 'w') as fh:
            fh.write('1\n2\n3\n4')

        self._test_fh = None

        @format1.reader(None)
        def reader(fh):
            self._test_fh = fh
            for value in [int(x) for x in fh.read().split('\n')]:
                yield value

        generator = self.registry.read(fp, format='format1')
        for a, b in zip_longest(generator, [1, 2, 3, 4]):
            self.assertEqual(a, b)
        self.assertTrue(self._test_fh.closed)

    def test_reader_does_not_exist(self):
        fh = StringIO()
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.registry.read(fh, format='not_a_format', into=TestClass)

        self.assertTrue(TestClass.__name__ in str(cm.exception))
        self.assertTrue('not_a_format' in str(cm.exception))

        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.registry.read(fh, format='not_a_format2')

        self.assertTrue('generator' in str(cm.exception))
        self.assertTrue('not_a_format2' in str(cm.exception))

#    def test_reader_is_not_generator(self):
#        fh = StringIO(u'1\n2\n3\n4')
#
#        @self.registry.register_sniffer('format')
#        def sniffer(fh):
#            return '1' in fh.readline(), {}
#
#        @self.registry.register_reader('format')
#        def reader(fh):
#            # Not a generator!
#            return TestClass([int(x) for x in fh.read().split('\n')])
#
#        with self.assertRaises(InvalidRegistrationError):
#            next(self.registry.read(fh, format='format'))
#
#        fh.close()

    def test_reader_exists_with_verify_true(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO(u'1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            self.was_verified = True
            return '1' in fh.readline(), {}

        @format1.reader(TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        self.was_verified = False
        instance = self.registry.read(fh, format='format1', into=TestClass,
                                      verify=True)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)
        self.assertTrue(self.was_verified)

        # Remove if read-context management is support in the future.
        fh.seek(0)

        self.was_verified = False
        instance = self.registry.read(fh, format='format1', into=TestClass)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)
        self.assertTrue(self.was_verified)

        fh.close()

    def test_warning_raised(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO(u'1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            self.was_verified = True
            return False, {}

        @format1.reader(TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            with self.assertRaises(FormatIdentificationWarning):
                self.was_verified = False
                instance = self.registry.read(fh, format='format1',
                                              into=TestClass, verify=True)
                self.assertEqual(TestClass([1, 2, 3, 4]), instance)
                self.assertTrue(self.was_verified)

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            with self.assertRaises(FormatIdentificationWarning):
                self.was_verified = False
                instance = self.registry.read(fh, format='format1',
                                              into=TestClass)
                self.assertEqual(TestClass([1, 2, 3, 4]), instance)
                self.assertTrue(self.was_verified)

        fh.close()

    def test_reader_exists_with_verify_false(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO(u'1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            self.was_verified = True
            return '1' in fh.readline(), {}

        @format1.reader(TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        self.was_verified = False
        instance = self.registry.read(fh, format='format1', into=TestClass,
                                      verify=False)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)
        self.assertFalse(self.was_verified)
        fh.close()

    def test_reader_exists_real_file(self):
        format1 = self.registry.create_format('format1')

        fp = self.fp1
        with open(fp, 'w') as fh:
            fh.write('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        instance = self.registry.read(fp, format='format1', into=TestClass)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)

    def test_read_kwargs_passed_generator(self):
        format1 = self.registry.create_format('format1')

        @format1.sniffer()
        def sniffer(fh):
            return True, {'arg1': 15, 'arg2': 'abc'}

        @format1.reader(None)
        def reader(fh, **kwargs):
            self.assertEqual(kwargs['arg1'], 15)
            self.assertEqual(kwargs['arg2'], 'abc')
            self.assertEqual(kwargs['arg3'], [1])
            yield

        next(self.registry.read(StringIO(), format='format1', arg3=[1]))

    def test_read_kwargs_passed_and_override(self):
        format1 = self.registry.create_format('format1')

        @format1.sniffer()
        def sniffer(fh):
            return True, {'arg1': 15, 'arg2': 'abc', 'override': 30}

        @format1.reader(TestClass)
        def reader(fh, **kwargs):
            self.assertEqual(kwargs['arg1'], 15)
            self.assertEqual(kwargs['arg2'], 'abc')
            self.assertEqual(kwargs['arg3'], [1])
            return

        self.registry.read(StringIO(u'notempty'), into=TestClass, arg3=[1])

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            # Should raise no warning and thus no error.
            self.registry.read(StringIO(u'notempty'), into=TestClass, arg3=[1],
                               override=30)
            # Should raise a warning and thus an error.
            with self.assertRaises(ArgumentOverrideWarning):
                self.registry.read(StringIO(u'notempty'), into=TestClass,
                                   arg3=[1], override=100)

    def test_that_encoding_is_used(self):
        format1 = self.registry.create_format('format1')

        fp = get_data_path('big5_file')

        @format1.sniffer()
        def sniffer(fh):
            return u'\u4f60' in fh.readline(), {}

        @format1.reader(TestClass)
        def reader(fh):
            self.assertEqual(self._expected_enc, fh.encoding)
            return TestClass(fh.readlines())

        self._expected_enc = 'big5'
        instance = self.registry.read(fp, into=TestClass, encoding='big5')
        self.assertEqual(TestClass([u'\u4f60\u597d\n']), instance)

    def test_file_sentinel_many(self):
        format1 = self.registry.create_format('format1')

        extra = get_data_path('real_file')
        extra_2 = get_data_path('real_file_2')
        fh = StringIO(u'1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(TestClass)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertEqual('a\nb\nc\nd\ne\n', extra.read())
            self.assertEqual('!\n@\n#\n$\n%\nThe realest.\n', extra_2.read())
            return TestClass([int(x) for x in fh.read().split('\n')])

        instance = self.registry.read(fh, format='format1', into=TestClass,
                                      extra=extra, extra_2=extra_2)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)

        fh.close()

    def test_file_sentinel_converted_to_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO(u'1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(TestClass)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)
            return TestClass([int(x) for x in fh.read().split('\n')])

        instance = self.registry.read(fh, format='format1', into=TestClass)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)

        fh.close()

    def test_file_sentinel_pass_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO(u'1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(TestClass)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)
            return TestClass([int(x) for x in fh.read().split('\n')])

        instance = self.registry.read(fh, format='format1', into=TestClass,
                                      extra=None)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)

        fh.close()

    def test_file_sentinel_generator_many(self):
        format1 = self.registry.create_format('format1')

        extra = get_data_path('real_file')
        extra_2 = get_data_path('real_file_2')
        fh = StringIO(u'1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(None)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertEqual('a\nb\nc\nd\ne\n', extra.read())
            self.assertEqual('!\n@\n#\n$\n%\nThe realest.\n', extra_2.read())
            yield TestClass([int(x) for x in fh.read().split('\n')])

        gen = self.registry.read(fh, format='format1', extra=extra,
                                 extra_2=extra_2)
        self.assertEqual(TestClass([1, 2, 3, 4]), next(gen))

        fh.close()

    def test_file_sentinel_converted_to_none_generator(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO(u'1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(None)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)
            yield TestClass([int(x) for x in fh.read().split('\n')])

        gen = self.registry.read(fh, format='format1')
        self.assertEqual(TestClass([1, 2, 3, 4]), next(gen))

        fh.close()

    def test_file_sentinel_pass_none_generator(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO(u'1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(None)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)
            yield TestClass([int(x) for x in fh.read().split('\n')])

        gen = self.registry.read(fh, format='format1', extra=None)
        self.assertEqual(TestClass([1, 2, 3, 4]), next(gen))

        fh.close()


class TestWrite(RegistryTest):
    def test_writer_does_not_exist(self):
        fh = StringIO()
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.registry.write({}, format='not_a_format', into=fh)

        self.assertTrue('not_a_format' in str(cm.exception))
        self.assertTrue(str(fh) in str(cm.exception))
        fh.close()

    def test_writer_exists(self):
        format1 = self.registry.create_format('format1')

        obj = TestClass(['1', '2', '3', '4'])
        fh = StringIO()

        @format1.writer(TestClass)
        def writer(obj, fh):
            fh.write(u'\n'.join(obj.list))

        self.registry.write(obj, format='format1', into=fh)
        fh.seek(0)
        self.assertEqual("1\n2\n3\n4", fh.read())
        fh.close()

    def test_writer_exists_real_file(self):
        format1 = self.registry.create_format('format1')

        obj = TestClass(['1', '2', '3', '4'])
        fp = self.fp1

        @format1.writer(TestClass)
        def writer(obj, fh):
            fh.write(u'\n'.join(obj.list))

        self.registry.write(obj, format='format1', into=fp)

        with io.open(fp) as fh:
            self.assertEqual(u"1\n2\n3\n4", fh.read())

    def test_writer_passed_kwargs(self):
        format1 = self.registry.create_format('format1')

        @format1.reader(None)
        def reader(fh):
            yield

        @format1.writer(None)
        def writer(obj, fh, **kwargs):
            self.assertEqual(kwargs['passed'], True)

        generator = self.registry.get_reader('format1', None)([])
        self.registry.write(generator, format='format1',
                            into=StringIO(), passed=True)

    def test_that_encoding_is_used(self):
        format1 = self.registry.create_format('format1')

        obj = TestClass([u'\u4f60\u597d\n'])  # Ni Hau
        fp = self.fp1

        @format1.writer(TestClass)
        def writer(obj, fh):
            fh.write(u''.join(obj.list))
            self.assertEqual(self._expected_encoding, fh.encoding)

        self._expected_encoding = 'big5'
        self.registry.write(obj, format='format1', into=fp, encoding='big5')

        with io.open(fp, mode='rb') as fh:
            # This would have been b'\xe4\xbd\xa0\xe5\xa5\xbd\n' in utf8
            self.assertEqual(b'\xa7A\xa6n\n', fh.read())

    def test_file_sentinel_many(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO()

        @format1.writer(TestClass)
        def writer(obj, fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            extra.write(u'oh yeah...')
            extra_2.write(u'oh no...')

        self.registry.write(TestClass([]), format='format1', into=fh,
                            extra=self.fp1, extra_2=self.fp2)
        with open(self.fp1) as f1:
            self.assertEqual('oh yeah...', f1.read())

        with open(self.fp2) as f2:
            self.assertEqual('oh no...', f2.read())

        fh.close()

    def test_file_sentinel_converted_to_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO()

        @format1.writer(TestClass)
        def writer(obj, fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)

        self.registry.write(TestClass([]), format='format1', into=fh)

        fh.close()

    def test_file_sentinel_pass_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO()

        @format1.writer(TestClass)
        def writer(obj, fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)

        self.registry.write(TestClass([]), format='format1', into=fh,
                            extra=None)

        fh.close()


# class TestInitializeOOPInterface(RegistryTest):
#    def setUp(self):
#        super(TestInitializeOOPInterface, self).setUp()
#
#        class UnassumingClass(object):
#            pass
#
#        class ClassWithDefault(object):
#            default_write_format = 'favfmt'
#
#        self.unassuming_class = UnassumingClass
#        self.class_with_default = ClassWithDefault
#
#    def test_no_readers_writers(self):
#        self.registry.initialize_oop_interface()
#        self.assertFalse(hasattr(self.unassuming_class, 'read'))
#        self.assertFalse(hasattr(self.unassuming_class, 'write'))
#        self.assertFalse(hasattr(self.class_with_default, 'read'))
#        self.assertFalse(hasattr(self.class_with_default, 'write'))
#
#    def test_readers_only(self):
#        @self.registry.register_reader('favfmt', self.unassuming_class)
#        def fvfmt_to_unasumming_class(fh):
#            return
#
#        @self.registry.register_reader('favfmt')
#        def fvfmt_to_gen(fh):
#            yield
#
#        @self.registry.register_reader('favfmt2', self.unassuming_class)
#        def fvfmt2_to_unasumming_class(fh):
#            return
#
#        self.registry.initialize_oop_interface()
#
#        self.assertTrue(hasattr(self.unassuming_class, 'read'))
#        self.assertFalse(hasattr(self.unassuming_class, 'write'))
#        self.assertFalse(hasattr(self.class_with_default, 'read'))
#        self.assertFalse(hasattr(self.class_with_default, 'write'))
#
#        self.assertIn('favfmt', self.unassuming_class.read.__doc__)
#        self.assertIn('favfmt2', self.unassuming_class.read.__doc__)
#
#    def test_writers_only(self):
#        @self.registry.register_writer('favfmt', self.class_with_default)
#        def favfmt(fh):
#            pass
#
#        @self.registry.register_writer('favfmt')
#        def gen_to_favfmt(fh):
#            pass
#
#        @self.registry.register_writer('favfmt2', self.class_with_default)
#        def favfmt2(fh):
#            pass
#
#        self.registry.initialize_oop_interface()
#
#        self.assertFalse(hasattr(self.unassuming_class, 'read'))
#        self.assertFalse(hasattr(self.unassuming_class, 'write'))
#        self.assertFalse(hasattr(self.class_with_default, 'read'))
#        self.assertTrue(hasattr(self.class_with_default, 'write'))
#
#        self.assertIn('favfmt', self.class_with_default.write.__doc__)
#        self.assertIn('favfmt2', self.class_with_default.write.__doc__)
#
#    def test_writers_no_default_format(self):
#        @self.registry.register_writer('favfmt', self.unassuming_class)
#        def favfmt(fh):
#            pass
#
#        @self.registry.register_writer('favfmt')
#        def gen_to_favfmt(fh):
#            pass
#
#        @self.registry.register_writer('favfmt2', self.unassuming_class)
#        def favfmt2(fh):
#            pass
#        with self.assertRaises(NotImplementedError) as cm:
#            self.registry.initialize_oop_interface()
#
#        self.assertIn('default_write_format', str(cm.exception))
#
#    def test_readers_writers(self):
#        @self.registry.register_reader('favfmt', self.unassuming_class)
#        def fvfmt_to_unasumming_class(fh):
#            return
#
#        @self.registry.register_reader('favfmt', self.class_with_default)
#        def fvfmt_to_class_w_default(fh):
#            return
#
#        @self.registry.register_reader('favfmt')
#        def fvfmt_to_gen(fh):
#            yield
#
#        @self.registry.register_reader('favfmt2', self.unassuming_class)
#        def fvfmt2_to_unasumming_class(fh):
#            return
#
#        @self.registry.register_reader('favfmt2', self.class_with_default)
#        def fvfmt2_to_class_w_default(fh):
#            return
#
#        @self.registry.register_writer('favfmt', self.class_with_default)
#        def favfmt(fh):
#            pass
#
#        @self.registry.register_writer('favfmt')
#        def gen_to_favfmt(fh):
#            pass
#
#        @self.registry.register_writer('favfmt2', self.class_with_default)
#        def favfmt2(fh):
#            pass
#
#        self.registry.initialize_oop_interface()
#
#        self.assertTrue(hasattr(self.unassuming_class, 'read'))
#        self.assertFalse(hasattr(self.unassuming_class, 'write'))
#
#        self.assertTrue(hasattr(self.class_with_default, 'read'))
#        self.assertTrue(hasattr(self.class_with_default, 'write'))
#
#        self.assertIn('favfmt', self.unassuming_class.read.__doc__)
#        self.assertIn('favfmt2', self.unassuming_class.read.__doc__)
#
#        self.assertIn('favfmt', self.class_with_default.read.__doc__)
#        self.assertIn('favfmt2', self.class_with_default.read.__doc__)
#
#        self.assertIn('favfmt', self.class_with_default.write.__doc__)
#        self.assertIn('favfmt2', self.class_with_default.write.__doc__)
#
#    def test_read_kwargs_passed(self):
#        self.was_called = False
#
#        @self.registry.register_sniffer('favfmt')
#        def fvfmt_sniffer(fh):
#            return True, {}
#
#        @self.registry.register_reader('favfmt', self.class_with_default)
#        def fvfmt_to_class_w_default(fh, **kwargs):
#            self.assertEqual('a', kwargs['a'])
#            self.assertEqual(123, kwargs['b'])
#            self.was_called = True
#
#        self.registry.initialize_oop_interface()
#        fh = StringIO(u'notempty')
#        self.class_with_default.read(fh, a='a', b=123)
#
#        self.assertTrue(self.was_called)
#        fh.close()
#
#    def test_write_kwargs_passed(self):
#        self.was_called = False
#
#        @self.registry.register_writer('favfmt', self.class_with_default)
#        def favfmt(obj, fh, **kwargs):
#            self.assertEqual('a', kwargs['a'])
#            self.assertEqual(123, kwargs['b'])
#            self.was_called = True
#
#        self.registry.initialize_oop_interface()
#        fh = StringIO()
#        self.class_with_default().write(fh, a='a', b=123)
#
#        self.assertTrue(self.was_called)
#        fh.close()


if __name__ == '__main__':
    unittest.main()
