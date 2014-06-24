# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.core.exception import RecordError
from skbio.parse.record_finder import (DelimitedRecordFinder,
                                       LabeledRecordFinder, LineGrouper,
                                       TailedRecordFinder)


class TailedRecordFinderTests(TestCase):
    def setUp(self):
        self.endswith_period = lambda x: x.endswith('.')
        self.period_tail_finder = TailedRecordFinder(self.endswith_period)

    def test_parsers(self):
        lines = '>abc\ndef\nz.\n>efg\nz.'.split()
        fl = self.period_tail_finder
        self.assertEqual(list(fl(lines)),
                         [['>abc', 'def', 'z.'], ['>efg', 'z.']])

    def test_parsers_empty(self):
        fl = self.period_tail_finder
        self.assertEqual(list(fl(['  ', '\n'])), [])
        self.assertEqual(list(fl([])), [])

    def test_parsers_strip(self):
        fl = self.period_tail_finder
        lines = '>abc  \n \t def\n  z. \t\n>efg \nz.'.split('\n')
        self.assertEqual(list(fl(lines)),
                         [['>abc', ' \t def', '  z.'], ['>efg', 'z.']])

    def test_parsers_leftover(self):
        f = self.period_tail_finder
        good = ['abc  \n',
                'def\n',
                '.\n',
                'ghi \n',
                'j.',
                ]
        blank = ['', '   ', '\t    \t\n\n']
        bad = ['abc']

        result = [['abc', 'def', '.'], ['ghi', 'j.']]

        self.assertEqual(list(f(good)), result)
        self.assertEqual(list(f(good + blank)), result)
        self.assertRaises(RecordError, list, f(good + bad))

        f2 = TailedRecordFinder(self.endswith_period, strict=False)
        self.assertEqual(list(f2(good + bad)), result + [['abc']])

    def test_parsers_ignore(self):
        def never(line):
            return False

        def ignore_labels(line):
            return (not line) or line.isspace() or line.startswith('#')

        lines = ['abc', '\n', '1.', 'def', '#ignore', '2.']
        self.assertEqual(list(TailedRecordFinder(self.endswith_period)(lines)),
                         [['abc', '1.'], ['def', '#ignore', '2.']])
        self.assertEqual(list(TailedRecordFinder(self.endswith_period,
                                                 ignore=never)(lines)),
                         [['abc', '', '1.'], ['def', '#ignore', '2.']])
        self.assertEqual(list(TailedRecordFinder(self.endswith_period,
                                                 ignore=ignore_labels)(lines)),
                         [['abc', '1.'], ['def', '2.']])


class DelimitedRecordFinderTests(TestCase):
    def test_parsers(self):
        lines = 'abc\ndef\n//\nefg\n//'.split()
        self.assertEqual(list(DelimitedRecordFinder('//')(lines)),
                         [['abc', 'def', '//'], ['efg', '//']])
        self.assertEqual(list(DelimitedRecordFinder('//', keep_delimiter=False)
                              (lines)),
                         [['abc', 'def'], ['efg']])

    def test_parsers_empty(self):
        self.assertEqual(list(DelimitedRecordFinder('//')(['  ', '\n'])), [])
        self.assertEqual(list(DelimitedRecordFinder('//')([])), [])

    def test_parsers_strip(self):
        lines = '  \t   abc  \n \t   def\n  // \t\n\t\t efg \n//'.split('\n')
        self.assertEqual(list(DelimitedRecordFinder('//')(lines)),
                         [['abc', 'def', '//'], ['efg', '//']])

    def test_parsers_error(self):
        good = ['  \t   abc  \n',
                '\t   def\n',
                '// \t\n',
                '\t\n',
                '\t efg \n',
                '\t\t//\n',
                ]
        blank = ['', '   ', '\t    \t\n\n']
        bad = ['abc']

        result = [['abc', 'def', '//'], ['efg', '//']]
        r = DelimitedRecordFinder('//')

        self.assertEqual(list(r(good)), result)
        self.assertEqual(list(r(good + blank)), result)
        try:
            list(r(good + bad))
        except RecordError:
            pass
        else:
            raise AssertionError("Parser failed to raise error on bad data")

        r = DelimitedRecordFinder('//', strict=False)
        self.assertEqual(list(r(good + bad)), result + [['abc']])

    def test_parsers_ignore(self):
        def never(line):
            return False

        def ignore_labels(line):
            return (not line) or line.isspace() or line.startswith('#')

        lines = ['>abc', '\n', '1', '$$', '>def', '#ignore', '2', '$$']
        self.assertEqual(list(DelimitedRecordFinder('$$')(lines)),
                         [['>abc', '1', '$$'], ['>def', '#ignore', '2', '$$']])
        self.assertEqual(list(DelimitedRecordFinder('$$',
                                                    ignore=never)(lines)),
                         [['>abc', '', '1', '$$'],
                          ['>def', '#ignore', '2', '$$']])
        self.assertEqual(
            list(DelimitedRecordFinder('$$', ignore=ignore_labels)(lines)),
            [['>abc', '1', '$$'], ['>def', '2', '$$']])


class LabeledRecordFinderTests(TestCase):
    def setUp(self):
        self.FastaLike = LabeledRecordFinder(lambda x: x.startswith('>'))

    def test_parsers(self):
        lines = '>abc\ndef\n//\n>efg\n//'.split()
        fl = self.FastaLike
        self.assertEqual(list(fl(lines)),
                         [['>abc', 'def', '//'], ['>efg', '//']])

    def test_parsers_empty(self):
        fl = self.FastaLike
        self.assertEqual(list(fl(['  ', '\n'])), [])
        self.assertEqual(list(fl([])), [])

    def test_parsers_strip(self):
        fl = self.FastaLike
        lines = '  \t   >abc  \n \t   def\n  // \t\n\t\t >efg \n//'.split('\n')
        self.assertEqual(list(fl(lines)),
                         [['>abc', 'def', '//'], ['>efg', '//']])

    def test_parsers_leftover(self):
        fl = self.FastaLike
        good = ['  \t   >abc  \n',
                '\t   def\n',
                '\t\n',
                '\t >efg \n',
                'ghi',
                ]
        blank = ['', '   ', '\t    \t\n\n']
        bad = ['>abc']

        result = [['>abc', 'def'], ['>efg', 'ghi']]

        self.assertEqual(list(fl(good)), result)
        self.assertEqual(list(fl(good + blank)), result)
        self.assertEqual(list(fl(good + bad)), result + [['>abc']])

    def test_parsers_ignore(self):
        def never(line):
            return False

        def ignore_labels(line):
            return (not line) or line.isspace() or line.startswith('#')

        def is_start(line):
            return line.startswith('>')

        lines = ['>abc', '\n', '1', '>def', '#ignore', '2']
        self.assertEqual(list(LabeledRecordFinder(is_start)(lines)),
                         [['>abc', '1'], ['>def', '#ignore', '2']])
        self.assertEqual(list(LabeledRecordFinder(is_start,
                                                  ignore=never)(lines)),
                         [['>abc', '', '1'], ['>def', '#ignore', '2']])
        self.assertEqual(list(LabeledRecordFinder(is_start,
                                                  ignore=ignore_labels)(
            lines)),
            [['>abc', '1'], ['>def', '2']])

    def test_constructor_is_none(self):
        lrf = LabeledRecordFinder(lambda x: x.strip().startswith('>'),
                                  constructor=None)
        lines = '  \t   >abc  \n \t   def\n  // \t\n\t\t >efg \n//'.split('\n')

        obs = list(lrf(lines))
        exp = [['  \t   >abc  ', ' \t   def', '  // \t'], ['\t\t >efg ', '//']]
        self.assertEqual(obs, exp)


class LineGrouperTests(TestCase):
    def test_parser(self):
        good = ['  \t   >abc  \n',
                '\t   def\n',
                '\t\n',
                '\t >efg \n',
                'ghi',
                ]
        c = LineGrouper(2)
        self.assertEqual(list(c(good)), [['>abc', 'def'], ['>efg', 'ghi']])
        c = LineGrouper(1)
        self.assertEqual(list(c(good)), [['>abc'], ['def'], ['>efg'], ['ghi']])
        c = LineGrouper(4)
        self.assertEqual(list(c(good)), [['>abc', 'def', '>efg', 'ghi']])
        # shouldn't work if not evenly divisible
        c = LineGrouper(3)
        self.assertRaises(RecordError, list, c(good))

    def test_parser_ignore(self):
        def never(line):
            return False

        def ignore_labels(line):
            return (not line) or line.isspace() or line.startswith('#')

        lines = ['abc', '\n', '1', 'def', '#ignore', '2']
        self.assertEqual(list(LineGrouper(1)(lines)),
                         [['abc'], ['1'], ['def'], ['#ignore'], ['2']])
        self.assertEqual(list(LineGrouper(1, ignore=never)(lines)),
                         [[i.strip()] for i in lines])
        self.assertEqual(list(LineGrouper(2, ignore=ignore_labels)(lines)),
                         [['abc', '1'], ['def', '2']])

    def test_constructor_is_none(self):
        lines = ['abc', ' def   ', ' ghi', 'jkl  ']

        # should strip
        exp = [['abc', 'def'], ['ghi', 'jkl']]
        obs = list(LineGrouper(2)(lines))
        self.assertEqual(obs, exp)

        # should not strip
        exp = [['abc', ' def   '], [' ghi', 'jkl  ']]
        obs = list(LineGrouper(2, constructor=None)(lines))
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
