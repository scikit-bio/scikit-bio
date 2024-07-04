# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from io import StringIO
import io
import itertools
import os
import unittest
import warnings
import types
from tempfile import mkstemp
from sys import platform

from skbio.io import (FormatIdentificationWarning, UnrecognizedFormatError,
                      ArgumentOverrideWarning, io_registry, sniff,
                      create_format)
from skbio.io.registry import (IORegistry, FileSentinel, Format,
                               DuplicateRegistrationError,
                               InvalidRegistrationError)
from skbio.util import get_data_path
from skbio.util._exception import TestingUtilError
from skbio import DNA, read, write


class MockClass:
    def __init__(self, list_):
        self.list = list_

    def __eq__(self, other):
        # They are only equal when the class is EXACTLY the same. We don't want
        # readers to return knockoff instances...
        return self.__class__ is other.__class__ and self.list == other.list

    def __repr__(self):
        return "%s(%s)" % (str(self.__class__.__name__), str(self.list))


class MockClassA(MockClass):
    pass


class MockClassB(MockClass):
    pass


class TestFormatAndIORegistry(unittest.TestCase):
    def test_add_duplicate_format(self):
        f = Format('Example')
        r = IORegistry()
        r.add_format(f)
        with self.assertRaises(DuplicateRegistrationError):
            r.add_format(Format('Example'))


class RegistryTest(unittest.TestCase):
    def setUp(self):
        self.registry = IORegistry()
        self.fd1, self.fp1 = mkstemp()
        self.fd2, self.fp2 = mkstemp()

    def tearDown(self):
        os.close(self.fd1)
        os.remove(self.fp1)
        os.close(self.fd2)
        os.remove(self.fp2)


class TestRegisterAndGetReader(RegistryTest):
    def test_get_reader_no_match(self):
        self.assertIs(None, self.registry.get_reader('not_a_format',
                                                     MockClass))

    def test_get_reader_when_only_writer_exists(self):
        format1 = self.registry.create_format('format1')

        @format1.writer(MockClass)
        def format_writer(fh):
            return

        self.assertEqual(None, self.registry.get_reader('format', MockClass))

    def test_register_reader_on_many(self):
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3')
        format4 = self.registry.create_format('format4', encoding='binary')
        format5 = self.registry.create_format('format5', encoding='binary')

        @format1.reader(MockClassA)
        def format1_reader(fh):
            return

        @format1.reader(MockClassB)
        def format1_reader_b(fh):
            return

        @format2.reader(MockClassA)
        def format2_reader(fh):
            return

        @format3.reader(MockClassB)
        def format3_reader(fh):
            return

        @format4.reader(MockClassA)
        def format4_reader(fh):
            return

        @format4.reader(MockClassB)
        def format4_reader_b(fh):
            return

        @format5.reader(None)
        def format5_reader(fh):
            return

        self.assertIs(format1_reader,
                      self.registry.get_reader('format1', MockClassA))

        self.assertIs(format1_reader_b,
                      self.registry.get_reader('format1', MockClassB))

        self.assertIs(format2_reader,
                      self.registry.get_reader('format2', MockClassA))

        self.assertIs(None, self.registry.get_reader('format2', MockClassB))

        self.assertIs(None, self.registry.get_reader('format3', MockClassA))

        self.assertIs(format3_reader,
                      self.registry.get_reader('format3', MockClassB))

        self.assertIs(format4_reader,
                      self.registry.get_reader('format4', MockClassA))

        self.assertIs(format4_reader_b,
                      self.registry.get_reader('format4', MockClassB))

        self.assertIs(format5_reader,
                      self.registry.get_reader('format5', None))

        self.assertIs(None, self.registry.get_reader('format5', MockClassA))

        self.assertIs(None, self.registry.get_reader('format5', MockClassB))

    def test_register_reader_over_existing(self):
        format1 = self.registry.create_format('format1')
        with self.assertRaises(DuplicateRegistrationError) as cm:
            @format1.reader(MockClassA)
            def format1_reader(fh):
                return

            @format1.reader(MockClassA)
            def duplicate_format1_reader(fh):
                return

        self.assertTrue('format1' in str(cm.exception))
        self.assertTrue('reader' in str(cm.exception))
        self.assertTrue(MockClassA.__name__ in str(cm.exception))

    def test_register_reader_over_existing_override(self):
        format1 = self.registry.create_format('format1')

        @format1.reader(MockClassA)
        def format1_reader(fh):
            return

        self.assertIs(format1_reader,
                      self.registry.get_reader('format1', MockClassA))

        @format1.reader(MockClassA, override=True)
        def duplicate_format1_reader(fh):
            return

        self.assertIs(duplicate_format1_reader,
                      self.registry.get_reader('format1', MockClassA))

    def test_mistype_reader_registration(self):
        format1 = self.registry.create_format('format1')

        with self.assertRaises(InvalidRegistrationError):
            @format1.reader
            def left_out_parens(fh):
                return


class TestRegisterAndGetWriter(RegistryTest):
    def test_get_writer_no_match(self):
        self.assertEqual(None, self.registry.get_writer('not_a_format',
                                                        MockClass))

    def test_get_writer_when_only_reader_exists(self):
        format = self.registry.create_format('format')

        @format.reader(MockClass)
        def format_reader(fh):
            return

        self.assertEqual(None, self.registry.get_writer('format', MockClass))

    def test_register_writer_on_many(self):
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3')
        format4 = self.registry.create_format('format4', encoding='binary')
        format5 = self.registry.create_format('format5', encoding='binary')

        @format1.writer(MockClassA)
        def format1_writer(obj, fh):
            return

        @format1.writer(MockClassB)
        def format1_writer_b(obj, fh):
            return

        @format2.writer(MockClassA)
        def format2_writer(obj, fh):
            return

        @format3.writer(MockClassB)
        def format3_writer(obj, fh):
            return

        @format4.writer(MockClassA)
        def format4_writer(fh):
            return

        @format4.writer(MockClassB)
        def format4_writer_b(fh):
            return

        @format5.writer(None)
        def format5_writer(fh):
            return

        self.assertEqual(format1_writer,
                         self.registry.get_writer('format1', MockClassA))

        self.assertEqual(format1_writer_b,
                         self.registry.get_writer('format1', MockClassB))

        self.assertEqual(format2_writer,
                         self.registry.get_writer('format2', MockClassA))

        self.assertEqual(None,
                         self.registry.get_writer('format2', MockClassB))

        self.assertEqual(None,
                         self.registry.get_writer('format3', MockClassA))

        self.assertEqual(format3_writer,
                         self.registry.get_writer('format3', MockClassB))

        self.assertIs(format4_writer,
                      self.registry.get_writer('format4', MockClassA))

        self.assertIs(format4_writer_b,
                      self.registry.get_writer('format4', MockClassB))

        self.assertIs(format5_writer,
                      self.registry.get_writer('format5', None))

        self.assertIs(None, self.registry.get_writer('format5', MockClassA))

        self.assertIs(None, self.registry.get_writer('format5', MockClassB))

    def test_register_writer_over_existing(self):
        format1 = self.registry.create_format('format1')
        with self.assertRaises(DuplicateRegistrationError) as cm:
            @format1.writer(MockClassA)
            def format1_writer(obj, fh):
                return

            @format1.writer(MockClassA)
            def duplicate_format1_writer(obj, fh):
                return

        self.assertTrue('format1' in str(cm.exception))
        self.assertTrue('writer' in str(cm.exception))
        self.assertTrue(MockClassA.__name__ in str(cm.exception))

    def test_register_writer_over_existing_override(self):
        format1 = self.registry.create_format('format1')

        @format1.writer(MockClassA)
        def format1_writer(obj, fh):
            return

        self.assertIs(format1_writer,
                      self.registry.get_writer('format1', MockClassA))

        @format1.writer(MockClassA, override=True)
        def duplicate_format1_writer(obj, fh):
            return

        self.assertIs(duplicate_format1_writer,
                      self.registry.get_writer('format1', MockClassA))

    def test_mistype_writer_registration(self):
        format1 = self.registry.create_format('format1')

        with self.assertRaises(InvalidRegistrationError):
            @format1.writer
            def left_out_parens(fh):
                return


class TestRegisterAndGetSniffer(RegistryTest):
    def test_get_sniffer_no_match(self):
        self.assertEqual(None, self.registry.get_sniffer('not_a_format'))

    def test_register_sniffer_on_many(self):
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3', encoding='binary')

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

    def test_mistype_sniffer_registration(self):
        format1 = self.registry.create_format('format1')

        with self.assertRaises(InvalidRegistrationError):
            @format1.sniffer
            def left_out_parens(fh):
                return


class TestListReadFormats(RegistryTest):
    def test_no_read_formats(self):
        format1 = self.registry.create_format('format1')

        @format1.reader(MockClassA)
        def this_isnt_on_clsB(fh):
            return

        self.assertEqual([], self.registry.list_read_formats(MockClassB))

    def test_one_read_format(self):
        format1 = self.registry.create_format('format1')

        @format1.reader(MockClass)
        def format1_cls(fh):
            return

        self.assertEqual(['format1'],
                         self.registry.list_read_formats(MockClass))

    def test_many_read_formats(self):
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3', encoding='binary')
        format4 = self.registry.create_format('format4')
        format5 = self.registry.create_format('format5', encoding='binary')

        @format1.reader(MockClassA)
        def format1_clsA(fh):
            return

        @format2.reader(MockClassA)
        def format2_clsA(fh):
            return

        @format3.reader(MockClassA)
        def format3_clsA(fh):
            return

        @format3.reader(MockClassB)
        def format3_clsB(fh):
            return

        @format4.reader(MockClassB)
        def format4_clsB(fh):
            return

        @format5.writer(MockClassA)
        def format5_clsA(fh):
            return

        formats = self.registry.list_read_formats(MockClassA)
        self.assertTrue('format1' in formats)
        self.assertTrue('format2' in formats)
        self.assertTrue('format3' in formats)
        self.assertTrue('format4' not in formats)
        self.assertTrue('format5' not in formats)


class TestListWriteFormats(RegistryTest):
    def test_no_write_formats(self):
        format1 = self.registry.create_format('format1')

        @format1.writer(MockClassA)
        def this_isnt_on_clsB(fh):
            return

        self.assertEqual([], self.registry.list_write_formats(MockClassB))

    def test_one_write_format(self):
        format1 = self.registry.create_format('format1')

        @format1.writer(MockClass)
        def format1_cls(fh):
            return

        self.assertEqual(['format1'],
                         self.registry.list_write_formats(MockClass))

    def test_many_write_formats(self):
        format1 = self.registry.create_format('format1')
        format2 = self.registry.create_format('format2')
        format3 = self.registry.create_format('format3', encoding='binary')
        format4 = self.registry.create_format('format4')
        format5 = self.registry.create_format('format5', encoding='binary')

        @format1.writer(MockClassA)
        def format1_clsA(fh):
            return

        @format2.writer(MockClassA)
        def format2_clsA(fh):
            return

        @format3.writer(MockClassA)
        def format3_clsA(fh):
            return

        @format3.writer(MockClassB)
        def format3_clsB(fh):
            return

        @format4.writer(MockClassB)
        def format4_clsB(fh):
            return

        @format5.reader(MockClassA)
        def format5_clsA(fh):
            return

        formats = self.registry.list_write_formats(MockClassA)

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
        # No sniffer for this format:
        self.registry.create_format('format5')

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

        @format3.reader(MockClass)
        def reader3(fh):
            return

        @format4.reader(MockClass)
        def reader4(fh):
            return

    def test_no_matches(self):
        fh = StringIO("no matches here")
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.registry.sniff(fh)
        self.assertTrue(str(fh) in str(cm.exception))

        fh.close()

    def test_one_match(self):
        fh = StringIO("contains a 3")
        self.assertEqual('format3', self.registry.sniff(fh)[0])

    def test_many_matches(self):
        fh = StringIO("1234 will match all")
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.registry.sniff(fh)
        self.assertTrue("format1" in str(cm.exception))
        self.assertTrue("format2" in str(cm.exception))
        self.assertTrue("format3" in str(cm.exception))
        self.assertTrue("format4" in str(cm.exception))
        fh.close()

    def test_that_encoding_is_used(self):
        formatx = self.registry.create_format('formatx')

        fp = get_data_path('big5_file')

        @formatx.sniffer()
        def sniffer(fh):
            self.assertEqual('big5', fh.encoding)
            return True, {}

        fmt, _ = self.registry.sniff(fp, encoding='big5')
        self.assertEqual(fmt, 'formatx')

    def test_passing_newline_raises_error(self):
        formatx = self.registry.create_format('formatx')

        fp = get_data_path('real_file')

        @formatx.sniffer()
        def sniffer(fh):
            return True, {}

        with self.assertRaisesRegex(TypeError, r'`newline`'):
            self.registry.sniff(fp, newline='\r')

    def test_non_default_encoding(self):
        big5_format = self.registry.create_format('big5_format',
                                                  encoding='big5')

        @big5_format.sniffer()
        def sniffer(fh):
            self.assertEqual(self._expected_encoding, fh.encoding)
            return True, {}

        self._expected_encoding = 'big5'
        fmt, _ = self.registry.sniff(self.fp1)
        self.assertEqual(fmt, 'big5_format')

        self._expected_encoding = 'UTF-8'
        fmt, _ = self.registry.sniff(self.fp1, encoding='UTF-8')
        self.assertEqual(fmt, 'big5_format')

    def test_non_default_newline(self):
        formatx = self.registry.create_format('formatx', newline='\r')

        fp = get_data_path('real_file')

        @formatx.sniffer()
        def sniffer(fh):
            with io.open(fp, newline=None) as f:
                content = f.read()
            self.assertEqual(content, 'a\nb\nc\nd\ne\n')
            return True, {}

        fmt, _ = self.registry.sniff(fp)
        self.assertEqual(fmt, 'formatx')

    def test_position_not_mutated_real_file(self):
        formatx = self.registry.create_format('formatx')

        @formatx.sniffer()
        def sniffer(fh):
            return True, {}

        with io.open(get_data_path('real_file')) as fh:
            if platform == "win32":
                fh.seek(3)
                self.registry.sniff(fh)
                self.assertEqual(fh.tell(), 3)
                self.assertEqual('b\n', fh.readline())
            else:
                fh.seek(2)
                self.registry.sniff(fh)
                self.assertEqual(fh.tell(), 2)
                self.assertEqual('b\n', fh.readline())

    def test_position_not_mutated_fileish(self):
        formatx = self.registry.create_format('formatx')

        @formatx.sniffer()
        def sniffer(fh):
            return True, {}

        fh = StringIO('a\nb\nc\nd\n')
        fh.seek(2)
        self.registry.sniff(fh)
        self.assertEqual('b\n', fh.readline())

    def test_sniff_with_errors_in_sniffer(self):
        formatx = self.registry.create_format('formatx', encoding='ascii')

        @formatx.sniffer()
        def sniffer(fh):
            raise Exception("OH NO!")

        fp = get_data_path('big5_file')
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('error')
            with self.assertRaises(FormatIdentificationWarning):
                fmt, _ = self.registry.sniff(fp)

    def test_sniff_with_encoding_errors(self):
        formatx = self.registry.create_format('formatx', encoding='ascii')

        @formatx.sniffer()
        def sniffer(fh):
            fh.read()
            return True, {}

        fp = get_data_path('big5_file')
        with self.assertRaises(UnrecognizedFormatError):
            fmt, _ = self.registry.sniff(fp, errors='strict')
        # errors is set to ignore by default, so our sniffer will return
        # true even though read() didn't entirely work for ascii
        fmt, _ = self.registry.sniff(fp)
        self.assertEqual(fmt, 'formatx')

    def test_binary_sniffer(self):
        binf = self.registry.create_format('binf', encoding='binary')

        @binf.sniffer()
        def sniffer(fh):
            self.assertIsInstance(fh, (io.BufferedReader, io.BufferedRandom))
            return True, {}

        fmt, _ = self.registry.sniff(self.fp1)
        self.assertEqual(fmt, 'binf')

    def test_text_sniffer(self):
        textf = self.registry.create_format('textf', encoding=None)

        @textf.sniffer()
        def sniffer(fh):
            self.assertIsInstance(fh, io.TextIOBase)
            return True, {}

        fmt, _ = self.registry.sniff(self.fp1)
        self.assertEqual(fmt, 'textf')

    def test_sniff_with_illegal_encoding(self):
        binf = self.registry.create_format('binf', encoding='binary')
        textf = self.registry.create_format('textf', encoding=None)

        @binf.sniffer()
        def binf_sniffer(fh):
            return True, {}

        @textf.sniffer()
        def textf_sniffer(fh):
            return True, {}

        # Should skip binary sniffers
        fmt, _ = self.registry.sniff(self.fp1, encoding=None)
        self.assertEqual(fmt, 'textf')
        # Should skip text sniffers
        fmt, _ = self.registry.sniff(self.fp1, encoding='binary')
        self.assertEqual(fmt, 'binf')

        with self.assertRaises(ValueError):
            self.registry.sniff(['some content\n'], encoding='binary')

        with self.assertRaises(ValueError):
            binf_sniffer(self.fp1, encoding=None)

        with self.assertRaises(ValueError):
            textf_sniffer(self.fp1, encoding='binary')

    def test_binary_fall_through(self):
        binf = self.registry.create_format('binf', encoding='binary')
        textf = self.registry.create_format('textf', encoding=None)

        @binf.sniffer()
        def binf_sniffer(fh):
            self._check_binf = True
            return False, {}

        @textf.sniffer()
        def textf_sniffer(fh):
            self._check_textf = True
            return True, {}

        self._check_binf = False
        self._check_textf = False

        fmt, _ = self.registry.sniff(self.fp1)
        self.assertEqual(fmt, 'textf')

        self.assertTrue(self._check_binf)
        self.assertTrue(self._check_textf)

    def test_sniff_gzip(self):
        expected = "This is some content\nIt occurs on more than one line\n"

        formata = self.registry.create_format('formata', encoding='binary')
        formatb = self.registry.create_format('formatb')
        formatc = self.registry.create_format('formatc')

        @formata.sniffer()
        def formata_sniffer(fh):
            self._check_f1 = True
            self.assertEqual(fh.read(), expected.encode('ascii'))
            return False, {}

        @formatb.sniffer()
        def formatb_sniffer(fh):
            self._check_f2 = True
            self.assertEqual(fh.read(), expected)
            return True, {}

        @formatc.sniffer()
        def formatc_sniffer(fh):
            self._check_f3 = True
            self.assertEqual(fh.read(), expected)
            return False, {}

        self._check_f1 = False
        self._check_f2 = False
        self._check_f3 = False
        self.registry.sniff(get_data_path('example_file.gz'))
        self.assertTrue(self._check_f1)
        self.assertTrue(self._check_f2)
        self.assertTrue(self._check_f3)

    def test_text_skip_binary(self):
        binf = self.registry.create_format('binf', encoding='binary')
        textf = self.registry.create_format('textf', encoding=None)

        @binf.sniffer()
        def binf_sniffer(fh):
            self._check_binf = True
            return True, {}

        @textf.sniffer()
        def textf_sniffer(fh):
            self._check_textf = True
            return True, {}

        self._check_binf = False
        self._check_textf = False

        fmt, _ = self.registry.sniff(['text'])
        self.assertEqual(fmt, 'textf')

        self.assertFalse(self._check_binf)
        self.assertTrue(self._check_textf)

        self._check_binf = False
        self._check_textf = False

        fmt, _ = self.registry.sniff(self.fp1, encoding=None)
        self.assertEqual(fmt, 'textf')

        self.assertFalse(self._check_binf)
        self.assertTrue(self._check_textf)

    def test_text_skip_text(self):
        binf = self.registry.create_format('binf', encoding='binary')
        textf = self.registry.create_format('textf', encoding=None)

        @binf.sniffer()
        def binf_sniffer(fh):
            self._check_binf = True
            return True, {}

        @textf.sniffer()
        def textf_sniffer(fh):
            self._check_textf = True
            return True, {}

        self._check_binf = False
        self._check_textf = False

        fmt, _ = self.registry.sniff(self.fp1, encoding='binary')
        self.assertEqual(fmt, 'binf')

        self.assertTrue(self._check_binf)
        self.assertFalse(self._check_textf)


class TestRead(RegistryTest):
    def test_format_and_into_are_none(self):
        fh = StringIO()
        with self.assertRaises(ValueError):
            self.registry.read(fh)

        fh.close()

    def test_format_is_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(MockClass)
        def reader(fh):
            self.assertIsInstance(fh, io.TextIOBase)
            return MockClass([int(x) for x in fh.read().split('\n')])

        instance = self.registry.read(fh, into=MockClass)
        self.assertEqual(MockClass([1, 2, 3, 4]), instance)
        fh.close()

    def test_into_is_none_and_no_generator_reader(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO('1\n2\n3\n4')

        @format1.reader(MockClass)
        def reader(fh):
            self.assertIsInstance(fh, io.TextIOBase)
            return

        with self.assertRaisesRegex(UnrecognizedFormatError,
                                    r"Cannot read 'format1'.*Possible.*include"
                                    ": MockClass"):
            self.registry.read(fh, format='format1')

    def test_into_is_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO('1\n2\n3\n4')

        @format1.reader(None)
        def reader(fh):
            self.assertIsInstance(fh, io.TextIOBase)
            yield from [int(x) for x in fh.read().split('\n')]

        generator = self.registry.read(fh, format='format1')
        self.assertIsInstance(generator, types.GeneratorType)
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
            yield from [int(x) for x in fh.read().split('\n')]

        generator = self.registry.read(fp, format='format1')
        for a, b in itertools.zip_longest(generator, [1, 2, 3, 4]):
            self.assertEqual(a, b)
        self.assertTrue(self._test_fh.closed)

    def test_reader_does_not_exist(self):
        fh = StringIO()
        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.registry.read(fh, format='not_a_format', into=MockClass)

        self.assertTrue(MockClass.__name__ in str(cm.exception))
        self.assertTrue('not_a_format' in str(cm.exception))

        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.registry.read(fh, format='not_a_format2')

        self.assertTrue('generator' in str(cm.exception))
        self.assertTrue('not_a_format2' in str(cm.exception))

    def test_reader_exists_with_verify_true(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            self.was_verified = True
            return '1' in fh.readline(), {}

        @format1.reader(MockClass)
        def reader(fh):
            return MockClass([int(x) for x in fh.read().split('\n')])

        self.was_verified = False
        instance = self.registry.read(fh, format='format1', into=MockClass,
                                      verify=True)
        self.assertEqual(MockClass([1, 2, 3, 4]), instance)
        self.assertTrue(self.was_verified)

        # Remove if read-context management is support in the future.
        fh.seek(0)

        self.was_verified = False
        instance = self.registry.read(fh, format='format1', into=MockClass)
        self.assertEqual(MockClass([1, 2, 3, 4]), instance)
        self.assertTrue(self.was_verified)

        fh.close()

    def test_warning_raised(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            self.was_verified = True
            return False, {}

        @format1.reader(MockClass)
        def reader(fh):
            return MockClass([int(x) for x in fh.read().split('\n')])

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            with self.assertRaises(FormatIdentificationWarning):
                self.was_verified = False
                instance = self.registry.read(fh, format='format1',
                                              into=MockClass, verify=True)
                self.assertEqual(MockClass([1, 2, 3, 4]), instance)
                self.assertTrue(self.was_verified)

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            with self.assertRaises(FormatIdentificationWarning):
                self.was_verified = False
                instance = self.registry.read(fh, format='format1',
                                              into=MockClass)
                self.assertEqual(MockClass([1, 2, 3, 4]), instance)
                self.assertTrue(self.was_verified)

        fh.close()

    def test_reader_exists_with_verify_false(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            self.was_verified = True
            return '1' in fh.readline(), {}

        @format1.reader(MockClass)
        def reader(fh):
            return MockClass([int(x) for x in fh.read().split('\n')])

        self.was_verified = False
        instance = self.registry.read(fh, format='format1', into=MockClass,
                                      verify=False)
        self.assertEqual(MockClass([1, 2, 3, 4]), instance)
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

        @format1.reader(MockClass)
        def reader(fh):
            return MockClass([int(x) for x in fh.read().split('\n')])

        instance = self.registry.read(fp, format='format1', into=MockClass)
        self.assertEqual(MockClass([1, 2, 3, 4]), instance)

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

        @format1.reader(MockClass)
        def reader(fh, **kwargs):
            self.assertEqual(kwargs['arg1'], 15)
            self.assertEqual(kwargs['arg2'], 'abc')
            self.assertEqual(kwargs['arg3'], [1])
            return

        self.registry.read(StringIO('notempty'), into=MockClass, arg3=[1])

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            # Should raise no warning and thus no error.
            self.registry.read(StringIO('notempty'), into=MockClass, arg3=[1],
                               override=30)
            # Should raise a warning and thus an error.
            with self.assertRaises(ArgumentOverrideWarning):
                self.registry.read(StringIO('notempty'), into=MockClass,
                                   arg3=[1], override=100)

    def test_that_encoding_is_used(self):
        format1 = self.registry.create_format('format1')

        fp = get_data_path('big5_file')

        @format1.sniffer()
        def sniffer(fh):
            return '\u4f60' in fh.readline(), {}

        @format1.reader(MockClass)
        def reader(fh):
            self.assertEqual(self._expected_enc, fh.encoding)
            return MockClass(fh.readlines())

        @format1.reader(None)
        def reader_gen(fh):
            self.assertEqual(self._expected_enc, fh.encoding)
            yield MockClass(fh.readlines())

        self._expected_enc = 'big5'
        instance = self.registry.read(fp, into=MockClass, encoding='big5')
        self.assertEqual(MockClass(['\u4f60\u597d\n']), instance)

        self._expected_enc = 'big5'
        gen = self.registry.read(fp, format='format1', encoding='big5')
        self.assertEqual(MockClass(['\u4f60\u597d\n']), next(gen))

    def test_non_default_encoding(self):
        format1 = self.registry.create_format('format1', encoding='big5')

        fp = get_data_path('big5_file')

        @format1.sniffer()
        def sniffer(fh):
            return True, {}

        @format1.reader(MockClass)
        def reader(fh):
            self.assertEqual(self._expected_enc, fh.encoding)
            return MockClass(fh.readlines())

        @format1.reader(None)
        def reader_gen(fh):
            self.assertEqual(self._expected_enc, fh.encoding)
            yield MockClass(fh.readlines())

        self._expected_enc = 'big5'
        instance = self.registry.read(fp, into=MockClass)
        self.assertEqual(MockClass(['\u4f60\u597d\n']), instance)

        gen = self.registry.read(fp, format='format1')
        self.assertEqual(MockClass(['\u4f60\u597d\n']), next(gen))
        gen.close()

        self._expected_enc = 'utf8'
        with self.assertRaises(UnicodeDecodeError):
            self.registry.read(fp, into=MockClass, encoding='utf8')

        with self.assertRaises(UnicodeDecodeError):
            self.registry.read(fp, format='format1', encoding='utf8')

    def test_passing_newline_raises_error(self):
        formatx = self.registry.create_format('formatx')

        fp = get_data_path('real_file')

        @formatx.sniffer()
        def sniffer(fh):
            return True, {}

        @formatx.reader(MockClass)
        def reader(fh):
            return MockClass(fh.readlines())

        @formatx.reader(None)
        def reader_gen(fh):
            yield MockClass(fh.readlines())

        with self.assertRaisesRegex(TypeError, r'`newline`'):
            self.registry.read(fp, into=MockClass, newline='\r')

        with self.assertRaisesRegex(TypeError, r'`newline`'):
            self.registry.read(fp, format='formatx', newline='\r')

    def test_non_default_newline(self):
        formatx = self.registry.create_format('formatx', newline='\r')

        fp = get_data_path('real_file')

        @formatx.sniffer()
        def sniffer(fh):
            return True, {}

        @formatx.reader(MockClass)
        def reader(fh):
            with io.open(fp, newline=None) as f:
                content = [''.join(f.readlines())]
            return MockClass(content)

        @formatx.reader(None)
        def reader_gen(fh):
            with io.open(fp, newline=None) as f:
                content = [''.join(f.readlines())]
            yield MockClass(content)

        instance = self.registry.read(fp, into=MockClass)
        self.assertEqual(instance, MockClass(['a\nb\nc\nd\ne\n']))

        gen = self.registry.read(fp, format='formatx')
        self.assertEqual(next(gen), MockClass(['a\nb\nc\nd\ne\n']))
        gen.close()

    def test_file_sentinel_many(self):
        format1 = self.registry.create_format('format1')

        extra = get_data_path('real_file')
        extra_2 = get_data_path('real_file_2')
        fh = StringIO('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(MockClass)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertEqual('a\nb\nc\nd\ne\n', extra.read())
            self.assertEqual('!\n@\n#\n$\n%\nThe realest.\n', extra_2.read())
            return MockClass([int(x) for x in fh.read().split('\n')])

        instance = self.registry.read(fh, format='format1', into=MockClass,
                                      extra=extra, extra_2=extra_2)
        self.assertEqual(MockClass([1, 2, 3, 4]), instance)

        fh.close()

    def test_file_sentinel_converted_to_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(MockClass)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)
            return MockClass([int(x) for x in fh.read().split('\n')])

        instance = self.registry.read(fh, format='format1', into=MockClass)
        self.assertEqual(MockClass([1, 2, 3, 4]), instance)

        fh.close()

    def test_file_sentinel_pass_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(MockClass)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)
            return MockClass([int(x) for x in fh.read().split('\n')])

        instance = self.registry.read(fh, format='format1', into=MockClass,
                                      extra=None)
        self.assertEqual(MockClass([1, 2, 3, 4]), instance)

        fh.close()

    def test_file_sentinel_generator_many(self):
        format1 = self.registry.create_format('format1')

        extra = get_data_path('real_file')
        extra_2 = get_data_path('real_file_2')
        fh = StringIO('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(None)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertEqual('a\nb\nc\nd\ne\n', extra.read())
            self.assertEqual('!\n@\n#\n$\n%\nThe realest.\n', extra_2.read())
            yield MockClass([int(x) for x in fh.read().split('\n')])

        gen = self.registry.read(fh, format='format1', extra=extra,
                                 extra_2=extra_2)
        self.assertEqual(MockClass([1, 2, 3, 4]), next(gen))

        fh.close()

    def test_file_sentinel_converted_to_none_generator(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(None)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)
            yield MockClass([int(x) for x in fh.read().split('\n')])

        gen = self.registry.read(fh, format='format1')
        self.assertEqual(MockClass([1, 2, 3, 4]), next(gen))

        fh.close()

    def test_file_sentinel_pass_none_generator(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO('1\n2\n3\n4')

        @format1.sniffer()
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @format1.reader(None)
        def reader(fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)
            yield MockClass([int(x) for x in fh.read().split('\n')])

        gen = self.registry.read(fh, format='format1', extra=None)
        self.assertEqual(MockClass([1, 2, 3, 4]), next(gen))

        fh.close()

    def test_read_with_illegal_encoding(self):
        binf = self.registry.create_format('binf', encoding='binary')
        textf = self.registry.create_format('textf', encoding=None)

        @binf.sniffer()
        def binf_sniffer(fh):
            return True, {}

        @binf.reader(MockClass)
        def binf_reader(fh):
            return MockClass(['bin'])

        @binf.reader(None)
        def binf_reader_gen(fh):
            yield MockClass(['bin'])

        @textf.sniffer()
        def textf_sniffer(fh):
            return True, {}

        @textf.reader(MockClass)
        def textf_reader(fh):
            return MockClass(['text'])

        @textf.reader(None)
        def textf_reader_gen(fh):
            yield MockClass(['text'])

        # Should skip binary sniffers
        instance = self.registry.read(self.fp1, encoding=None, into=MockClass)
        self.assertEqual(instance, MockClass(['text']))
        gen = self.registry.read(self.fp1, encoding=None, format='textf')
        self.assertEqual(next(gen), MockClass(['text']))
        gen.close()
        # Should skip text sniffers
        instance = self.registry.read(self.fp1, encoding='binary',
                                      into=MockClass)
        self.assertEqual(instance, MockClass(['bin']))
        gen = self.registry.read(self.fp1, encoding='binary', format='binf')
        self.assertEqual(next(gen), MockClass(['bin']))
        gen.close()

        with self.assertRaises(ValueError):
            self.registry.read(['some content\n'], encoding='binary',
                               into=MockClass)

        with self.assertRaises(ValueError):
            self.registry.read(['some content\n'], format='textf',
                               encoding='binary', into=MockClass)

        with self.assertRaises(ValueError):
            self.registry.read(['some content\n'], format='textf',
                               encoding='binary', verify=False, into=MockClass)

        with self.assertRaises(ValueError):
            self.registry.read(['some content\n'], format='textf',
                               encoding='binary')

        with self.assertRaises(ValueError):
            self.registry.read(['some content\n'], format='textf',
                               encoding='binary', verify=False)

        with self.assertRaises(ValueError):
            self.registry.read(self.fp1, format='binf',
                               encoding=None, into=MockClass)

        with self.assertRaises(ValueError):
            self.registry.read(self.fp1, format='binf',
                               encoding=None, verify=False, into=MockClass)

        with self.assertRaises(ValueError):
            self.registry.read(self.fp1, format='binf',
                               encoding=None)

        with self.assertRaises(ValueError):
            self.registry.read(self.fp1, format='binf',
                               encoding=None, verify=False)

    def test_read_with_binary_encoding(self):
        binf = self.registry.create_format('binf', encoding='binary')

        @binf.reader(MockClass)
        def reader1(fh):
            self.assertIsInstance(fh, (io.BufferedReader, io.BufferedRandom))
            return MockClass(['woo'])

        @binf.reader(None)
        def reader2(fh):
            self.assertIsInstance(fh, (io.BufferedReader, io.BufferedRandom))
            yield MockClass(['woo'])

        instance = self.registry.read(self.fp1, format='binf', verify=False,
                                      into=MockClass)
        self.assertEqual(MockClass(['woo']), instance)

        gen = self.registry.read(self.fp1, format='binf', verify=False,
                                 into=None)
        self.assertEqual(MockClass(['woo']), next(gen))
        gen.close()

    def test_io_kwargs_passed(self):
        format1 = self.registry.create_format('format1')

        @format1.sniffer()
        def sniffer(fh):
            return True, {}

        @format1.reader(MockClass)
        def reader1(fh):
            self.assertEqual(fh.errors, 'replace')
            return MockClass(['woo'])

        @format1.reader(None)
        def reader1_gen(fh):
            self.assertEqual(fh.errors, 'replace')
            yield MockClass(['woo'])

        obj = self.registry.read(self.fp1, into=MockClass, errors='replace')
        self.assertEqual(obj, MockClass(['woo']))
        gen = self.registry.read(self.fp1, format='format1', errors='replace')
        self.assertEqual(next(gen), MockClass(['woo']))
        gen.close()

    def test_read_empty_file_gen_with_format(self):
        format1 = self.registry.create_format('format1')

        @format1.sniffer()
        def sniffer(fh):
            return True, {}

        @format1.reader(None)
        def reader1(fh):
            return
            yield

        with io.StringIO("") as fh:
            gen = self.registry.read(fh, format='format1')

        self.assertEqual(list(gen), [])


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

        obj = MockClass(['1', '2', '3', '4'])
        fh = StringIO()

        @format1.writer(MockClass)
        def writer(obj, fh):
            self.assertIsInstance(fh, io.TextIOBase)
            fh.write('\n'.join(obj.list))

        self.registry.write(obj, format='format1', into=fh)
        fh.seek(0)
        self.assertEqual("1\n2\n3\n4", fh.read())
        fh.close()

    def test_writer_exists_real_file(self):
        format1 = self.registry.create_format('format1')

        obj = MockClass(['1', '2', '3', '4'])
        fp = self.fp1

        @format1.writer(MockClass)
        def writer(obj, fh):
            self.assertIsInstance(fh, io.TextIOBase)
            fh.write('\n'.join(obj.list))

        self.registry.write(obj, format='format1', into=fp)

        with io.open(fp) as fh:
            self.assertEqual("1\n2\n3\n4", fh.read())

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

        obj = MockClass(['\u4f60\u597d\n'])  # Ni Hau
        fp = self.fp1

        @format1.writer(MockClass)
        def writer(obj, fh):
            fh.write(''.join(obj.list))
            self.assertEqual(self._expected_encoding, fh.encoding)

        self._expected_encoding = 'big5'
        self.registry.write(obj, format='format1', into=fp, encoding='big5')

        with io.open(fp, mode='rb') as fh:
            # This would have been b'\xe4\xbd\xa0\xe5\xa5\xbd\n' in utf8
            self.assertEqual(b"\xa7A\xa6n\n", fh.read().replace(b"\r\n", b"\n"))

    def test_non_default_encoding(self):
        format1 = self.registry.create_format('format1', encoding='big5')

        obj = MockClass(['\u4f60\u597d\n'])  # Ni Hau
        fp = self.fp1

        @format1.writer(MockClass)
        def writer(obj, fh):
            fh.write(''.join(obj.list))
            self.assertEqual(self._expected_encoding, fh.encoding)

        self._expected_encoding = 'big5'
        self.registry.write(obj, format='format1', into=fp)

        with io.open(fp, mode='rb') as fh:
            self.assertEqual(b"\xa7A\xa6n\n", fh.read().replace(b"\r\n", b"\n"))

        self._expected_encoding = 'utf8'
        self.registry.write(obj, format='format1', into=fp, encoding='utf8')

        with io.open(fp, mode='rb') as fh:
            self.assertEqual(b"\xe4\xbd\xa0\xe5\xa5\xbd\n",
                             fh.read().replace(b"\r\n", b"\n"))

    def test_that_newline_is_used(self):
        format1 = self.registry.create_format('format1')

        obj = MockClass(['a\n', 'b\n', 'c\n'])
        fp = self.fp1

        @format1.writer(MockClass)
        def writer(obj, fh):
            fh.write(''.join(obj.list))

        self.registry.write(obj, format='format1', into=fp, newline='\r')

        with io.open(fp, mode='rb') as fh:
            self.assertEqual(b'a\rb\rc\r', fh.read())

    def test_non_default_newline(self):
        format1 = self.registry.create_format('format1', newline='\r')

        obj = MockClass(['a\n', 'b\n', 'c\n'])
        fp = self.fp1

        @format1.writer(MockClass)
        def writer(obj, fh):
            fh.write(''.join(obj.list))

        self.registry.write(obj, format='format1', into=fp)

        with io.open(fp, mode='rb') as fh:
            self.assertEqual(b'a\rb\rc\r', fh.read())

        self.registry.write(obj, format='format1', into=fp, newline='\n')

        with io.open(fp, mode='rb') as fh:
            self.assertEqual(b'a\nb\nc\n', fh.read())

    def test_file_sentinel_many(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO()

        @format1.writer(MockClass)
        def writer(obj, fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            extra.write('oh yeah...')
            extra_2.write('oh no...')

        self.registry.write(MockClass([]), format='format1', into=fh,
                            extra=self.fp1, extra_2=self.fp2)
        with open(self.fp1) as f1:
            self.assertEqual('oh yeah...', f1.read())

        with open(self.fp2) as f2:
            self.assertEqual('oh no...', f2.read())

        fh.close()

    def test_file_sentinel_converted_to_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO()

        @format1.writer(MockClass)
        def writer(obj, fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)

        self.registry.write(MockClass([]), format='format1', into=fh)

        fh.close()

    def test_file_sentinel_pass_none(self):
        format1 = self.registry.create_format('format1')

        fh = StringIO()

        @format1.writer(MockClass)
        def writer(obj, fh, extra=FileSentinel, other=2, extra_2=FileSentinel):
            self.assertIsNone(extra)
            self.assertIsNone(extra_2)

        self.registry.write(MockClass([]), format='format1', into=fh,
                            extra=None)

        fh.close()

    def test_write_with_illegal_encoding(self):
        binf = self.registry.create_format('binf', encoding='binary')
        textf = self.registry.create_format('textf', encoding=None)

        @binf.writer(MockClass)
        def writer(obj, fh):
            pass

        @textf.writer(MockClass)
        def writer2(obj, fh):
            pass

        with self.assertRaises(ValueError):
            self.registry.write(MockClass([]), into=self.fp1, format='binf',
                                encoding=None)

        with self.assertRaises(ValueError):
            self.registry.write(MockClass([]), into=self.fp1, format='textf',
                                encoding='binary')

    def test_write_binary_format(self):
        format1 = self.registry.create_format('format1', encoding='binary')

        obj = MockClass([b'a\n', b'b\n', b'c\n'])
        fp = self.fp1

        @format1.writer(MockClass)
        def writer(obj, fh):
            self.assertIsInstance(fh, (io.BufferedWriter, io.BufferedRandom))
            fh.write(b''.join(obj.list))

        self.registry.write(obj, format='format1', into=fp)

        with io.open(fp, mode='rb') as fh:
            self.assertEqual(b'a\nb\nc\n', fh.read())

    def test_io_kwargs_passed(self):
        format1 = self.registry.create_format('format1', encoding='ascii')

        obj = MockClass(['a\n', 'b\n', 'c\n'])
        fp = self.fp1
        f = io.BytesIO()

        @format1.writer(MockClass)
        def writer(obj, fh):
            iterator = iter(obj.list)
            fh.write(next(iterator))
            fh.flush()  # Flush should be a noop for bz2

            for line in iterator:
                fh.write(line)

        self.registry.write(obj, format='format1', into=fp, compression='bz2', newline='\n')
        self.registry.write(obj, format='format1', into=f, compression='bz2', newline='\n')
        expected = (
            b'BZh91AY&SY\x03\x89\x0c\xa6\x00\x00\x01\xc1\x00\x00\x108\x00 \x00'
            b'!\x9ah3M\x1c\xb7\x8b\xb9"\x9c(H\x01\xc4\x86S\x00')

        with io.open(fp, mode='rb') as fh:
            self.assertEqual(expected, fh.read())

        self.assertEqual(expected, f.getvalue())


class TestMonkeyPatch(RegistryTest):
    def setUp(self):
        super(TestMonkeyPatch, self).setUp()

        class UnassumingClass:
            pass

        class ClassWithDefault:
            default_write_format = 'favfmt'

        class NoMonkeySee:
            pass

        self.unassuming_class = UnassumingClass
        self.class_with_default = ClassWithDefault
        self.no_monkey_see = NoMonkeySee

    def test_no_readers_writers(self):
        self.registry.monkey_patch()
        self.assertFalse(hasattr(self.unassuming_class, 'read'))
        self.assertFalse(hasattr(self.unassuming_class, 'write'))
        self.assertFalse(hasattr(self.class_with_default, 'read'))
        self.assertFalse(hasattr(self.class_with_default, 'write'))

    def test_readers_only(self):
        favfmt = self.registry.create_format('favfmt')
        favfmt2 = self.registry.create_format('favfmt2')

        @favfmt.reader(self.unassuming_class)
        def fvfmt_to_unasumming_class(fh):
            return

        @favfmt.reader(None)
        def fvfmt_to_gen(fh):
            yield

        @favfmt2.reader(self.unassuming_class)
        def fvfmt2_to_unasumming_class(fh):
            return

        self.registry.monkey_patch()

        self.assertTrue(hasattr(self.unassuming_class, 'read'))
        self.assertFalse(hasattr(self.unassuming_class, 'write'))
        self.assertFalse(hasattr(self.class_with_default, 'read'))
        self.assertFalse(hasattr(self.class_with_default, 'write'))

        self.assertIn('favfmt', self.unassuming_class.read.__doc__)
        self.assertIn('favfmt2', self.unassuming_class.read.__doc__)

    def test_writers_only(self):
        favfmt = self.registry.create_format('favfmt')
        favfmt2 = self.registry.create_format('favfmt2')

        @favfmt.writer(self.class_with_default)
        def favfmt_writer(fh):
            pass

        @favfmt.writer(None)
        def gen_to_favfmt(fh):
            pass

        @favfmt2.writer(self.class_with_default)
        def favfmt2_writer(fh):
            pass

        self.registry.monkey_patch()

        self.assertFalse(hasattr(self.unassuming_class, 'read'))
        self.assertFalse(hasattr(self.unassuming_class, 'write'))
        self.assertFalse(hasattr(self.class_with_default, 'read'))
        self.assertTrue(hasattr(self.class_with_default, 'write'))

        self.assertIn('favfmt', self.class_with_default.write.__doc__)
        self.assertIn('favfmt2', self.class_with_default.write.__doc__)

    def test_writers_no_default_format(self):
        favfmt = self.registry.create_format('favfmt')
        favfmt2 = self.registry.create_format('favfmt2')

        @favfmt.writer(self.unassuming_class)
        def favfmt_writer(fh):
            pass

        @favfmt.writer(None)
        def gen_to_favfmt(fh):
            pass

        @favfmt2.writer(self.unassuming_class)
        def favfmt2_writer(fh):
            pass
        with self.assertRaises(NotImplementedError) as cm:
            self.registry.monkey_patch()

        self.assertIn('default_write_format', str(cm.exception))

    def test_readers_writers(self):
        favfmt = self.registry.create_format('favfmt')
        favfmt2 = self.registry.create_format('favfmt2')

        @favfmt.reader(self.unassuming_class)
        def fvfmt_to_unasumming_class(fh):
            return

        @favfmt.reader(self.class_with_default)
        def fvfmt_to_class_w_default(fh):
            return

        @favfmt.reader(None)
        def fvfmt_to_gen(fh):
            yield

        @favfmt2.reader(self.unassuming_class)
        def fvfmt2_to_unasumming_class(fh):
            return

        @favfmt2.reader(self.class_with_default)
        def fvfmt2_to_class_w_default(fh):
            return

        @favfmt.writer(self.class_with_default)
        def favfmt_writer(fh):
            pass

        @favfmt.writer(None)
        def gen_to_favfmt(fh):
            pass

        @favfmt2.writer(self.class_with_default)
        def favfmt2_writer(fh):
            pass

        @favfmt2.reader(self.no_monkey_see, monkey_patch=True)
        def favfmt2_to_monkey(fh):
            pass

        @favfmt2.writer(self.no_monkey_see, monkey_patch=False)
        def monkey_to_favfmt2(fh):
            pass

        self.registry.monkey_patch()

        self.assertTrue(hasattr(self.unassuming_class, 'read'))
        self.assertFalse(hasattr(self.unassuming_class, 'write'))

        self.assertTrue(hasattr(self.class_with_default, 'read'))
        self.assertTrue(hasattr(self.class_with_default, 'write'))

        self.assertTrue(hasattr(self.no_monkey_see, 'read'))
        self.assertFalse(hasattr(self.no_monkey_see, 'write'))

        self.assertIn('favfmt', self.unassuming_class.read.__doc__)
        self.assertIn('favfmt2', self.unassuming_class.read.__doc__)

        self.assertIn('favfmt', self.class_with_default.read.__doc__)
        self.assertIn('favfmt2', self.class_with_default.read.__doc__)

        self.assertIn('favfmt', self.class_with_default.write.__doc__)
        self.assertIn('favfmt2', self.class_with_default.write.__doc__)

        self.assertIn('favfmt2', self.no_monkey_see.read.__doc__)

    def test_read_kwargs_passed(self):
        favfmt = self.registry.create_format('favfmt')
        self.was_called = False

        @favfmt.sniffer()
        def fvfmt_sniffer(fh):
            return True, {}

        @favfmt.reader(self.class_with_default)
        def fvfmt_to_class_w_default(fh, **kwargs):
            self.assertEqual('a', kwargs['a'])
            self.assertEqual(123, kwargs['b'])
            self.was_called = True

        self.registry.monkey_patch()
        fh = StringIO('notempty')
        self.class_with_default.read(fh, a='a', b=123)

        self.assertTrue(self.was_called)
        fh.close()

    def test_write_kwargs_passed(self):
        favfmt = self.registry.create_format('favfmt')
        self.was_called = False

        @favfmt.writer(self.class_with_default)
        def favfmt_writer(obj, fh, **kwargs):
            self.assertEqual('a', kwargs['a'])
            self.assertEqual(123, kwargs['b'])
            self.was_called = True

        self.registry.monkey_patch()
        fh = StringIO()
        self.class_with_default().write(fh, a='a', b=123)

        self.assertTrue(self.was_called)
        fh.close()


class TestModuleFunctions(unittest.TestCase):

    def test_sniff_matches(self):
        exp = io_registry.sniff(['(a, b);'])
        result = sniff(['(a, b);'])
        self.assertEqual(exp, result)
        self.assertEqual('newick', exp[0])
        self.assertEqual({}, exp[1])

    def test_read_matches(self):
        input = ['>\n', 'ACGT\n']
        exp = io_registry.read(input, into=DNA)
        result = read(input, into=DNA)
        self.assertEqual(exp, result)
        self.assertEqual(exp, DNA('ACGT', metadata={'id': '',
                                                    'description': ''}))

    def test_write_matches(self):
        input = DNA('ACGT')
        exp = io_registry.write(input, format='fasta', into=[])
        result = write(input, format='fasta', into=[])
        self.assertEqual(exp, result)
        self.assertEqual(exp, ['>\n', 'ACGT\n'])

    def test_create_format_matches(self):
        with self.assertRaises(DuplicateRegistrationError):
            io_registry.create_format('fasta')

        with self.assertRaises(DuplicateRegistrationError):
            create_format('fasta')


if __name__ == '__main__':
    unittest.main()
