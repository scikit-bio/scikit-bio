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

from skbio.io import (DuplicateRegistrationError, UnprovenFormatWarning,
                      UnrecognizedFormatError, ArgumentOverrideWarning)
from skbio.io._registry import empty_file_sniffer


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

        self.assertEqual(None, self.module.get_reader(['not_a_format',
                                                       'still_not'],
                                                      TestClass))

        self.assertEqual(None, self.module.get_reader('Nope, Sorry',
                                                      TestClass))

    def test_register_reader_on_generator(self):
        @self.module.register_reader('format1')
        def format1_reader_generator(fh):
            yield

        @self.module.register_reader(['compound', 'format'])
        def compound_format_reader(fh1, fh2):
            yield

        self.assertEqual(format1_reader_generator,
                         self.module.get_reader('format1'))

        self.assertEqual(format1_reader_generator,
                         self.module.get_reader('format1', None))

        self.assertEqual(compound_format_reader,
                         self.module.get_reader(['compound', 'format']))

        self.assertEqual(compound_format_reader,
                         self.module.get_reader('compound, format'))

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


class TestRegisterAndGetWriter(RegistryTest):
    def test_get_writer_no_match(self):
        self.assertEqual(None, self.module.get_writer('not_a_format',
                                                      TestClass))

        self.assertEqual(None, self.module.get_writer(['not_a_format',
                                                       'still_not'],
                                                      TestClass))

        self.assertEqual(None, self.module.get_writer('Nope, Sorry',
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

        @self.module.register_writer(['compound', 'format'])
        def compound_format_writer(fh1, fh2):
            yield

        self.assertEqual(format1_writer_generator,
                         self.module.get_writer('format1'))

        self.assertEqual(format1_writer_generator,
                         self.module.get_writer('format1', None))

        self.assertEqual(compound_format_writer,
                         self.module.get_writer(['compound', 'format']))

        self.assertEqual(compound_format_writer,
                         self.module.get_writer('compound, format'))

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

    def test_register_compound_sniffer(self):
        with self.assertRaises(ValueError):
            @self.module.register_sniffer(['f1', 'f2'])
            def this_wont_work(fh1, fh2):
                return False, {}

        with self.assertRaises(ValueError):
            @self.module.register_sniffer('f1, f2')
            def this_still_wont_work(fh1, fh2):
                return False, {}

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

        compound_sniffer = self.module.get_sniffer(['format3',
                                                    'format2'])

        compound_sniffer2 = self.module.get_sniffer('format3, format2')

        dont_do_this_ever = self.module.get_sniffer([['format2', 'format1'],
                                                     'format3'])

        self.assertTrue(compound_sniffer([fh3, fh2])[0])
        self.assertTrue(compound_sniffer2([fh3, fh2])[0])
        self.assertTrue(not compound_sniffer2([fh2, fh3])[0])
        self.assertTrue(dont_do_this_ever([[fh2, fh1], fh3])[0])

        with self.assertRaises(ValueError):
            dont_do_this_ever(fh1)

        with self.assertRaises(ValueError):
            dont_do_this_ever([fh1, fh2, fh3])

    def test_register_sniffer_over_existing(self):
        with self.assertRaises(DuplicateRegistrationError) as cm:
            @self.module.register_sniffer('format1')
            def format1_sniffer(fh):
                return False, {}

            @self.module.register_sniffer('format1')
            def duplicate_format1_sniffer(fh):
                return False, {}

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

        @self.module.register_writer('format5', TestClassA)
        def format5_clsA(fh):
            return

        @self.module.register_reader('formatB, formatA', TestClassA)
        def formatAB_clsA(fh):
            return

        @self.module.register_reader(['formatX', 'formatY'], TestClassA)
        def formatXY_clsA(fh):
            return

        @self.module.register_reader(['formatM', 'formatN'], TestClassB)
        def formatMN_clsB(fh):
            return

        formats = self.module.list_read_formats(TestClassA)
        self.assertTrue('format1' in formats)
        self.assertTrue('format2' in formats)
        self.assertTrue('format3' in formats)
        self.assertTrue('format4' not in formats)
        self.assertTrue('format5' not in formats)
        self.assertTrue('formatB, formatA' in formats)
        self.assertTrue('formatX, formatY' in formats)
        self.assertTrue('formatM, formatN' not in formats)


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

        @self.module.register_reader('format5', TestClassA)
        def format5_clsA(fh):
            return

        @self.module.register_writer('formatB, formatA', TestClassA)
        def formatAB_clsA(fh):
            return

        @self.module.register_writer(['formatX', 'formatY'], TestClassA)
        def formatXY_clsA(fh):
            return

        @self.module.register_writer(['formatM', 'formatN'], TestClassB)
        def formatMN_clsB(fh):
            return

        formats = self.module.list_write_formats(TestClassA)

        self.assertTrue('format1' in formats)
        self.assertTrue('format2' in formats)
        self.assertTrue('format3' in formats)
        self.assertTrue('format4' not in formats)
        self.assertTrue('format5' not in formats)
        self.assertTrue('formatB, formatA' in formats)
        self.assertTrue('formatX, formatY' in formats)
        self.assertTrue('formatM, formatN' not in formats)


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

        with self.assertRaises(UnrecognizedFormatError) as cm:
            self.module.sniff([fh, fh2], cls=TestClassB)

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

    def test_compound_matches(self):
        fh = StringIO(u'1')
        fh2 = StringIO(u'2')

        self.assertEqual('format1, format2', self.module.sniff([fh, fh2])[0])

        self.assertEqual('format2, format1', self.module.sniff([fh2, fh])[0])

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

    def test_compound_matches_w_cls(self):
        fh = StringIO(u'3')
        fh2 = StringIO(u'4')

        self.assertEqual('format3, format4',
                         self.module.sniff([fh, fh2], cls=TestClass)[0])

        self.assertEqual('format4, format3',
                         self.module.sniff([fh2, fh], cls=TestClass)[0])

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

    def test_that_mode_is_used_compound(self):
        fp1 = self.fp1
        with open(fp1, 'w') as fh:
            fh.write('@\n#\n')

        fp2 = self.fp2
        with open(fp2, 'w') as fh:
            fh.write('!\n#\n')

        @self.module.register_sniffer('c1')
        def c1_sniffer(fh):
            self.assertEqual(self.expected_mode, fh.mode)
            return '@' in fh.readline(), {}

        @self.module.register_sniffer('c2')
        def c2_sniffer(fh):
            self.assertEqual(self.expected_mode, fh.mode)
            return '!' in fh.readline(), {}

        self.expected_mode = 'U'
        self.module.sniff([fp1, fp2])

        self.expected_mode = 'r'
        self.module.sniff([fp1, fp2], mode='r')


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

    def test_into_is_none_compound_format(self):
        fh = StringIO(u'1\n3')
        fh2 = StringIO(u'2\n4')

        @self.module.register_reader(['odd', 'even'])
        def reader(odd, even):
            for o, e in zip(odd, even):
                yield int(o.rstrip('\n'))
                yield int(e.rstrip('\n'))

        generator = self.module.read([fh, fh2], format='odd, even')
        first_run = True
        for a, b in zip(generator, [1, 2, 3, 4]):
            if first_run:
                fh.seek(3)
                fh2.seek(2)
                first_run = False
            self.assertEqual(a, b)
            self.assertEqual(3, fh.tell())
            self.assertEqual(2, fh2.tell())

        fh2.seek(0)
        fh.seek(0)

        generator = self.module.read([fh2, fh], format='even, odd')
        first_run = True
        for a, b in zip(generator, [1, 2, 3, 4]):
            if first_run:
                fh.seek(5)
                fh2.seek(1)
                first_run = False
            self.assertEqual(a, b)
            self.assertEqual(5, fh.tell())
            self.assertEqual(1, fh2.tell())

        with self.assertRaises(ValueError):
            self.module.read([fh], format='even, odd')

        fh.close()
        fh2.close()

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

        self.was_verified = False
        instance = self.module.read(fh, format='format', into=TestClass)
        self.assertEqual(TestClass([1, 2, 3, 4]), instance)
        self.assertTrue(self.was_verified)

        fh.close()

    def test_reader_compound_format_w_verify(self):
        fh = StringIO(u'1\n3')
        fh2 = StringIO(u'2\n4')
        fh3 = StringIO(u'9\n9')

        @self.module.register_reader('odd, even, nines', TestClass)
        def reader(odd, even, nines):
            state = []
            for o, e, n in zip(odd, even, nines):
                state.append(int(o.rstrip('\n')))
                state.append(int(e.rstrip('\n')))
                state.append(int(n.rstrip('\n')))
            return TestClass(state)

        @self.module.register_sniffer('nines')
        def nines_sniffer(fh):
            self.was_verified_n = True
            return '9' in fh.readline(), {}

        @self.module.register_sniffer('odd')
        def odd_sniffer(fh):
            self.was_verified_o = True
            return '1' in fh.readline(), {}

        @self.module.register_sniffer('even')
        def even_sniffer(fh):
            self.was_verified_e = True
            return '2' in fh.readline(), {}

        self.was_verified_e = False
        self.was_verified_o = False
        self.was_verified_n = False
        instance = self.module.read([fh2, fh3, fh], format=['even', 'nines',
                                                            'odd'],
                                    into=TestClass, verify=True)
        self.assertEqual(TestClass([1, 2, 9, 3, 4, 9]), instance)
        self.assertTrue(self.was_verified_e)
        self.assertTrue(self.was_verified_o)
        self.assertTrue(self.was_verified_n)

        self.was_verified_e = False
        self.was_verified_o = False
        self.was_verified_n = False
        instance = self.module.read([fh, fh2, fh3], format=['odd',
                                                            'even', 'nines'],
                                    into=TestClass)
        self.assertEqual(TestClass([1, 2, 9, 3, 4, 9]), instance)
        self.assertTrue(self.was_verified_e)
        self.assertTrue(self.was_verified_o)
        self.assertTrue(self.was_verified_n)

        fh.close()
        fh2.close()
        fh3.close()

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
            with self.assertRaises(UnprovenFormatWarning):
                self.was_verified = False
                instance = self.module.read(fh, format='format',
                                            into=TestClass, verify=True)
                self.assertEqual(TestClass([1, 2, 3, 4]), instance)
                self.assertTrue(self.was_verified)

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            with self.assertRaises(UnprovenFormatWarning):
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

    def test_reader_into_none_w_mutate_fh(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_sniffer('format')
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @self.module.register_reader('format')
        def reader(fh):
            for x in fh.read().split('\n'):
                yield int(x)

        fh.seek(0)
        generator = self.module.read(fh, format='format', mutate_fh=True)
        for a, b in zip(generator, [1, 2, 3, 4]):
            self.assertEqual(a, b)
        self.assertNotEqual(0, fh.tell())
        fh.close()

    def test_reader_w_mutate_fh(self):
        fh = StringIO(u'1\n2\n3\n4')

        @self.module.register_sniffer('format')
        def sniffer(fh):
            return '1' in fh.readline(), {}

        @self.module.register_reader('format', TestClass)
        def reader(fh):
            return TestClass([int(x) for x in fh.read().split('\n')])

        fh.seek(0)
        instance = self.module.read(fh, format='format', into=TestClass,
                                    mutate_fh=True)
        self.assertNotEqual(0, fh.tell())

        self.assertEqual(TestClass([1, 2, 3, 4]), instance)
        fh.close()

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

        self.module.read(StringIO(), into=TestClass, arg3=[1])

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("error")
            # Should raise no warning and thus no error.
            self.module.read(StringIO(), into=TestClass, arg3=[1],
                             override=30)
            # Should raise a warning and thus an error.
            with self.assertRaises(ArgumentOverrideWarning):
                self.module.read(StringIO(), into=TestClass, arg3=[1],
                                 override=100)

    def test_read_kwargs_passed_w_compound_format(self):
        @self.module.register_sniffer('format1')
        def format1_sniffer(fh):
            return True, {'partial': 1, 'overlap': 1, 'list': 0}

        @self.module.register_sniffer('format2')
        def format2_sniffer(fh):
            return True, {'arg1': 'a', 'overlap': 2, 'list': [1, 2]}

        @self.module.register_sniffer('format3')
        def format3_sniffer(fh):
            return True, {'partial': 3, 'overlap': 3, 'list': [3, 4]}

        @self.module.register_reader(['format1', 'format2', 'format3'])
        def reader(f1, f2, f3, **kwargs):
            self.assertEqual(kwargs['partial'], [1, None, 3])
            self.assertEqual(kwargs['overlap'], [1, 2, 3])
            self.assertEqual(kwargs['arg1'], [None, 'a', None])
            self.assertEqual(kwargs['provided'], [True, True, True])
            self.assertEqual(kwargs['list'], [0, [1, 2], [3, 4]])
            return

        fh1 = StringIO()
        fh2 = StringIO()
        fh3 = StringIO()
        self.module.read([fh1, fh2, fh3],
                         format='format1, format2, format3',
                         provided=True)
        self.module.read([fh3, fh1, fh2],
                         format='format3, format1, format2',
                         provided=True)
        fh1.close()
        fh2.close()
        fh3.close()

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

    def test_writer_compound_format(self):
        fh = StringIO()
        fh2 = StringIO()
        fh3 = StringIO()
        fhb = StringIO()
        fh2b = StringIO()
        fh3b = StringIO()

        @self.module.register_writer('odd, even, nines', TestClass)
        def writer(obj, odd, even, nines):
            i = 0
            while i+2 < len(obj.list):
                odd.write(str(obj.list[i]) + u"\n")
                even.write(str(obj.list[i+1]) + u"\n")
                nines.write(str(obj.list[i+2]) + u"\n")
                i = i+3

        self.module.write(TestClass([1, 2, 9, 3, 4, 9]), into=[fh2, fh3, fh],
                          format=['even', 'nines', 'odd'])

        self.module.write(TestClass([1, 2, 9, 3, 4, 9]), into=[fhb, fh2b,
                                                               fh3b],
                          format='odd, even, nines')

        fh.seek(0)
        fh2.seek(0)
        fh3.seek(0)
        fhb.seek(0)
        fh2b.seek(0)
        fh3b.seek(0)

        self.assertEqual(u'1\n3\n', fh.read())
        self.assertEqual(u'1\n3\n', fhb.read())

        self.assertEqual(u'2\n4\n', fh2.read())
        self.assertEqual(u'2\n4\n', fh2b.read())

        self.assertEqual(u'9\n9\n', fh3.read())
        self.assertEqual(u'9\n9\n', fh3b.read())

        fh.close()
        fh2.close()
        fh3.close()
        fhb.close()
        fh2b.close()
        fh3b.close()

        with self.assertRaises(ValueError):
            self.module.write(TestClass([1, 2, 9]), into=[fh],
                              format='even, odd, nines')

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
