#!/usr/bin/env python
"""Unit tests for recordfinders: parsers that group the lines for a record.
"""


from bipy.parse.record import RecordError
from bipy.parse.record_finder import (DelimitedRecordFinder, 
                    LabeledRecordFinder, LineGrouper,TailedRecordFinder)
from bipy.util.unit_test import TestCase, main

class TailedRecordFinderTests(TestCase):
    """Tests of the TailedRecordFinder factory function."""
    def setUp(self):
        """Define a standard TailedRecordFinder"""
        self.endswith_period = lambda x: x.endswith('.')
        self.period_tail_finder = TailedRecordFinder(self.endswith_period)
        
    def test_parsers(self):
        """TailedRecordFinder should split records into lines correctly"""
        lines = '>abc\ndef\nz.\n>efg\nz.'.split()
        fl = self.period_tail_finder
        self.assertEqual(list(fl(lines)), \
            [['>abc', 'def', 'z.'], ['>efg','z.']])

    def test_parsers_empty(self):
        """TailedRecordFinder should return empty list on empty lines"""
        fl = self.period_tail_finder
        self.assertEqual(list(fl(['  ','\n'])), [])
        self.assertEqual(list(fl([])), [])

    def test_parsers_strip(self):
        """TailedRecordFinder should trim each line correctly"""
        fl = self.period_tail_finder
        lines = '>abc  \n \t def\n  z. \t\n>efg \nz.'.split('\n')
        self.assertEqual(list(fl(lines)), \
            [['>abc', ' \t def', '  z.'], ['>efg','z.']])
   
    def test_parsers_leftover(self):
        """TailedRecordFinder should raise error or yield leftover"""
        f = self.period_tail_finder
        good =  [   'abc  \n',
                    'def\n',
                    '.\n',
                    'ghi \n',
                    'j.',
                ]
        blank = ['', '   ', '\t    \t\n\n']
        bad = ['abc']

        result = [['abc', 'def','.'], ['ghi','j.']]
        
        self.assertEqual(list(f(good)), result)
        self.assertEqual(list(f(good+blank)), result)
        self.assertRaises(RecordError, list, f(good+bad))
        
        f2 = TailedRecordFinder(self.endswith_period, strict=False)
        self.assertEqual(list(f2(good+bad)), result + [['abc']])

    def test_parsers_ignore(self):
        """TailedRecordFinder should skip lines to ignore."""
        def never(line):
            return False
        
        def ignore_labels(line):
            return (not line) or line.isspace() or line.startswith('#')
        
        lines = ['abc','\n','1.','def','#ignore','2.']
        self.assertEqual(list(TailedRecordFinder(self.endswith_period)(lines)), 
            [['abc', '1.'],['def','#ignore','2.']])
        self.assertEqual(list(TailedRecordFinder(self.endswith_period, 
            ignore=never)(lines)),
            [['abc', '', '1.'],['def','#ignore','2.']])
        self.assertEqual(list(TailedRecordFinder(self.endswith_period, 
            ignore=ignore_labels)(lines)),
            [['abc','1.'],['def','2.']])

class DelimitedRecordFinderTests(TestCase):
    """Tests of the DelimitedRecordFinder factory function."""
    def test_parsers(self):
        """DelimitedRecordFinder should split records into lines correctly"""
        lines = 'abc\ndef\n//\nefg\n//'.split()
        self.assertEqual(list(DelimitedRecordFinder('//')(lines)), \
            [['abc', 'def', '//'], ['efg','//']])
        self.assertEqual(list(DelimitedRecordFinder('//', keep_delimiter=False)
            (lines)), \
            [['abc', 'def'], ['efg']])

    def test_parsers_empty(self):
        """DelimitedRecordFinder should return empty list on empty lines"""
        self.assertEqual(list(DelimitedRecordFinder('//')(['  ','\n'])), [])
        self.assertEqual(list(DelimitedRecordFinder('//')([])), [])

    def test_parsers_strip(self):
        """DelimitedRecordFinder should trim each line correctly"""
        lines = '  \t   abc  \n \t   def\n  // \t\n\t\t efg \n//'.split('\n')
        self.assertEqual(list(DelimitedRecordFinder('//')(lines)), \
            [['abc', 'def', '//'], ['efg','//']])
   
    def test_parsers_error(self):
        """DelimitedRecordFinder should raise RecordError if trailing data"""
        good =  [   '  \t   abc  \n',
                    '\t   def\n',
                    '// \t\n',
                    '\t\n',
                    '\t efg \n',
                    '\t\t//\n',
                ]
        blank = ['', '   ', '\t    \t\n\n']
        bad = ['abc']

        result = [['abc', 'def', '//'], ['efg','//']]
        r = DelimitedRecordFinder('//')
        
        self.assertEqual(list(r(good)), result)
        self.assertEqual(list(r(good+blank)), result)
        try:
            list(r(good+bad))
        except RecordError:
            pass
        else:
            raise AssertionError, "Parser failed to raise error on bad data"

        r = DelimitedRecordFinder('//', strict=False)
        self.assertEqual(list(r(good+bad)), result + [['abc']])

    def test_parsers_ignore(self):
        """DelimitedRecordFinder should skip lines to ignore."""
        def never(line):
            return False
        
        def ignore_labels(line):
            return (not line) or line.isspace() or line.startswith('#')
        
        lines = ['>abc','\n','1', '$$', '>def','#ignore','2', '$$']
        self.assertEqual(list(DelimitedRecordFinder('$$')(lines)), 
            [['>abc', '1', '$$'],['>def','#ignore','2', '$$']])
        self.assertEqual(list(DelimitedRecordFinder('$$', 
            ignore=never)(lines)),
            [['>abc', '', '1', '$$'],['>def','#ignore','2','$$']])
        self.assertEqual(list(DelimitedRecordFinder('$$', 
            ignore=ignore_labels)(lines)),
            [['>abc','1','$$'],['>def','2','$$']])
 
class LabeledRecordFinderTests(TestCase):
    """Tests of the LabeledRecordFinder factory function."""
    def setUp(self):
        """Define a standard LabeledRecordFinder"""
        self.FastaLike = LabeledRecordFinder(lambda x: x.startswith('>'))
        
    def test_parsers(self):
        """LabeledRecordFinder should split records into lines correctly"""
        lines = '>abc\ndef\n//\n>efg\n//'.split()
        fl = self.FastaLike
        self.assertEqual(list(fl(lines)), \
            [['>abc', 'def', '//'], ['>efg','//']])

    def test_parsers_empty(self):
        """LabeledRecordFinder should return empty list on empty lines"""
        fl = self.FastaLike
        self.assertEqual(list(fl(['  ','\n'])), [])
        self.assertEqual(list(fl([])), [])

    def test_parsers_strip(self):
        """LabeledRecordFinder should trim each line correctly"""
        fl = self.FastaLike
        lines = '  \t   >abc  \n \t   def\n  // \t\n\t\t >efg \n//'.split('\n')
        self.assertEqual(list(fl(lines)), \
            [['>abc', 'def', '//'], ['>efg','//']])
   
    def test_parsers_leftover(self):
        """LabeledRecordFinder should not raise RecordError if last line label"""
        fl = self.FastaLike
        good =  [   '  \t   >abc  \n',
                    '\t   def\n',
                    '\t\n',
                    '\t >efg \n',
                    'ghi',
                ]
        blank = ['', '   ', '\t    \t\n\n']
        bad = ['>abc']

        result = [['>abc', 'def'], ['>efg','ghi']]
        
        self.assertEqual(list(fl(good)), result)
        self.assertEqual(list(fl(good+blank)), result)
        self.assertEqual(list(fl(good+bad)), result + [['>abc']])

    def test_parsers_ignore(self):
        """LabeledRecordFinder should skip lines to ignore."""
        def never(line):
            return False
        
        def ignore_labels(line):
            return (not line) or line.isspace() or line.startswith('#')
        
        def is_start(line):
            return line.startswith('>')
        
        lines = ['>abc','\n','1','>def','#ignore','2']
        self.assertEqual(list(LabeledRecordFinder(is_start)(lines)), 
            [['>abc', '1'],['>def','#ignore','2']])
        self.assertEqual(list(LabeledRecordFinder(is_start, 
            ignore=never)(lines)),
            [['>abc', '', '1'],['>def','#ignore','2']])
        self.assertEqual(list(LabeledRecordFinder(is_start, 
            ignore=ignore_labels)(lines)),
            [['>abc','1'],['>def','2']])
 
class LineGrouperTests(TestCase):
    """Tests of the LineGrouper class."""
    def test_parser(self):
        """LineGrouper should return n non-blank lines at a time"""
        good =  [   '  \t   >abc  \n',
                    '\t   def\n',
                    '\t\n',
                    '\t >efg \n',
                    'ghi',
                ]
        c = LineGrouper(2)
        self.assertEqual(list(c(good)), [['>abc', 'def'],['>efg','ghi']])
        c = LineGrouper(1)
        self.assertEqual(list(c(good)), [['>abc'], ['def'],['>efg'],['ghi']])
        c = LineGrouper(4)
        self.assertEqual(list(c(good)), [['>abc', 'def','>efg','ghi']])
        #shouldn't work if not evenly divisible
        c = LineGrouper(3)
        self.assertRaises(RecordError, list, c(good))

    def test_parser_ignore(self):
        """LineGrouper should skip lines to ignore."""
        def never(line):
            return False
        
        def ignore_labels(line):
            return (not line) or line.isspace() or line.startswith('#')
        
        lines = ['abc','\n','1','def','#ignore','2']
        self.assertEqual(list(LineGrouper(1)(lines)), 
            [['abc'], ['1'],['def'],['#ignore'],['2']])
        self.assertEqual(list(LineGrouper(1, ignore=never)(lines)),
            [[i.strip()] for i in lines])
        self.assertEqual(list(LineGrouper(2, ignore=ignore_labels)(lines)),
            [['abc','1'],['def','2']])
        
         
if __name__ == '__main__':
    main()
