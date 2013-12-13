#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""Tests of transformation and composition functions .
"""

from bipy.util.transform import (apply_each, bools, bool_each,
    conjoin, all, both,
    disjoin, any, either, negate, none, neither, compose, compose_many,
    per_shortest, per_longest, for_seq, 
    has_field, extract_field, test_field, index, test_container, 
    trans_except, trans_all, make_trans, find_any, find_no, find_all,
    keep_chars,exclude_chars, reorder, reorder_inplace,
    float_from_string, first, last, first_in_set, last_in_set, 
    first_not_in_set, last_not_in_set, first_index, last_index, 
    first_index_in_set, last_index_in_set, first_index_not_in_set,
    last_index_not_in_set, perm, comb, cross_comb, _increment_comb, identity,
    select)

from bipy.util.unit_test import TestCase, main

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Sandra Smit", "Zongzhi Liu"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

class has_x(object):
    #convenience class for has_field and related functions
    def __init__(self, x):
        self.x = x
    def __hash__(self):
        return hash(self.x)
    def __str__(self):
        return str(self.x)

class has_y(object):
    #convenience class for has_field and related functions
    def __init__(self, y):
        self.y = y
    def __hash__(self):
        return hash(self.y)
    def __str__(self):
        return str(self.y)


class metafunctionsTests(TestCase):
    """Tests of standalone functions."""
    def setUp(self):
        """Define some standard functions and data."""
        self.Numbers = range(20)
        self.SmallNumbers = range(3)
        self.SmallNumbersRepeated = range(5) * 4
        self.Letters = 'abcde'
        self.Mixed = list(self.Letters) + range(5)
        self.firsts = 'ab2'
        self.seconds = '0bc'

        self.is_char = lambda x: isinstance(x, str) and len(x) == 1
        self.is_vowel = lambda x: x in 'aeiou'
        self.is_consonant = lambda x: x not in 'aeiuo'
        self.is_number = lambda x: isinstance(x, int)
        self.is_odd_number = lambda x: x%2
        self.is_odd_letter = lambda x: x in 'acegikmoqs'
        self.is_zero = lambda x: x == 0
        self.is_small = lambda x: x < 3
        self.double = lambda x: x * 2
        self.minusone = lambda x: x - 1

        #function to test *args, **kwargs)
        self.is_alpha_digit = lambda first, second: \
            first.isalpha() and second.isdigit()
        self.is_digit_alpha = lambda first, second: \
            first.isdigit() and second.isalpha()

    def test_identity(self):
        """should return same object"""
        foo = [1,'a',lambda x: x]
        exp = id(foo)
        self.assertEqual(id(identity(foo)), exp)

    def test_select_sequence(self):
        """select should work on a sequence with a list of indices"""
        chars = 'abcdefghij'
        strings = list(chars)

        tests = {   (0,):['a'],
                    (-1,):['j'],
                    (0, 2, 4): ['a', 'c', 'e'],
                    (9,8,7,6,5,4,3,2,1,0):list('jihgfedcba'),
                    (-8, 8): ['c', 'i'],
                    ():[],
                }
        for test, result in tests.items():
            self.assertEqual(select(test, chars), result)
            self.assertEqual(select(test, strings), result)

    def test_select_empty(self):
        """select should raise error if indexing into empty sequence"""
        self.assertRaises(IndexError, select, [1], [])

    def test_select_mapping(self):
        """select should return the values corresponding to a list of keys"""
        values = {'a':5, 'b':2, 'c':4, 'd':6, 'e':7}
        self.assertEqual(select('abc', values), [5,2,4])
        self.assertEqual(select(['e','e','e'], values), [7,7,7])
        self.assertEqual(select(('e', 'b', 'a'), values), [7, 2, 5])
        #check that it raises KeyError on anything out of range
        self.assertRaises(KeyError, select, 'abx', values)

    def test_apply_each(self):
        """apply_each should apply each function to args, kwargs"""
        self.assertEqual(apply_each( \
            [self.is_char, self.is_vowel, self.is_consonant, self.is_number], \
            self.Letters[0]), [True, True, False, False])
            
        self.assertEqual(apply_each( \
            [self.is_char, self.is_vowel, self.is_consonant, self.is_number], \
            self.Letters[1]), [True, False, True, False])

        self.assertEqual(apply_each( \
            [self.double, self.minusone], self.SmallNumbers[0]), [0, -1])

        self.assertEqual(apply_each( \
            [self.double, self.minusone], self.SmallNumbers[1]), [2, 0])

        expects = [[True, False], [False, False], [False, True]]
        for i in range(len(expects)):
            self.assertEqual(apply_each( \
                [self.is_alpha_digit, self.is_digit_alpha], 
                self.firsts[i], self.seconds[i]), expects[i])
            self.assertEqual(apply_each( \
                [self.is_alpha_digit, self.is_digit_alpha], 
                self.firsts[i], second=self.seconds[i]), expects[i])
            self.assertEqual(apply_each( \
                [self.is_alpha_digit, self.is_digit_alpha], 
                second=self.seconds[i], first=self.firsts[i]), expects[i])

    def test_bools(self):
        """bools should convert items to True or False."""
        self.assertEqual(bools(self.Letters), [True]*5)
        self.assertEqual(bools(self.Numbers), [False] + [True]*19)

    def test_bool_each(self):
        """bool_each should return boolean version of applying each f to args"""
        self.assertEqual(bool_each([self.double, self.minusone], \
            self.SmallNumbers[0]), [False, True])

        self.assertEqual(bool_each([self.double, self.minusone], \
            self.SmallNumbers[1]), [True, False])

    def test_conjoin(self):
        """conjoin should return True if all components True"""
        self.assertEqual(conjoin([self.is_odd_letter,self.is_vowel],'a'), True)
        self.assertEqual(conjoin([self.is_odd_letter,self.is_vowel], x='b'), 
            False)
        self.assertEqual(conjoin([self.is_odd_letter,self.is_vowel],'c'), False)
        self.assertEqual(conjoin([self.is_odd_letter,self.is_vowel],'e'), True)
        #technically, this one should be true as well, but I left it off to 
        #have an even vowel test case...
        self.assertEqual(conjoin([self.is_odd_letter,self.is_vowel],'u'), False)
        #should short-circuit, i.e. not evaluate later cases after False
        self.assertEqual(conjoin([self.is_odd_letter, self.fail], 'b'), False)
        self.assertRaises(AssertionError, conjoin, \
            [self.is_odd_letter, self.fail], 'a')
        
    def test_all(self):
        """all should return a function returning True if all components True"""
        odd_vowel = all([self.is_odd_letter, self.is_vowel, self.is_char])
        self.assertEqual(odd_vowel('a'), True)
        self.assertEqual(map(odd_vowel, 'abceu'), 
            [True,False,False,True,False])
        odd_number = all([self.is_odd_number, self.is_number])
        self.assertEqual(map(odd_number, range(5)), [False,True]*2+[False])
        #should short-circuit, i.e. not evaluate later cases after False
        self.assertEqual(all([self.is_odd_letter, self.fail])('b'), False)
        self.assertRaises(AssertionError, all([self.is_odd_letter,self.fail]),\
            'a')

    def test_both(self):
        """both should return True if both components True"""
        odd_vowel = both(self.is_odd_letter, self.is_vowel)
        self.assertEqual(map(odd_vowel, 'abcu'), [True,False,False,False])
        #should short-circuit
        self.assertEqual(both(self.is_odd_letter, self.fail)('b'), False)
        self.assertRaises(AssertionError, both(self.is_odd_letter, self.fail),\
            'a')

    def test_disjoin(self):
        """disjoin should return True if any component True"""
        self.assertEqual(disjoin([self.is_odd_letter,self.is_vowel], 'a'), True)
        self.assertEqual(disjoin([self.is_odd_letter,self.is_vowel], 'b'),False)
        self.assertEqual(disjoin([self.is_odd_letter,self.is_vowel], 'c'), True)
        self.assertEqual(disjoin([self.is_odd_letter,self.is_vowel], x='u'),
            True)
        #should short-circuit after first True
        self.assertEqual(disjoin([self.is_odd_letter, self.fail], 'a'), True)
        self.assertRaises(AssertionError, \
            disjoin, [self.is_odd_letter, self.fail], 'b')
        
    def test_any(self):
        """any should return a function returning True if any component True"""
        odd_vowel = any([self.is_odd_letter, self.is_vowel])
        self.assertEqual(odd_vowel('a'), True)
        self.assertEqual(map(odd_vowel, 'abceu'), [True,False,True,True,True])
        odd = any([self.is_odd_number, self.is_small])
        self.assertEqual(map(odd, range(5)), [True]*4+[False])
        #should short-circuit after first True
        self.assertEqual(any([self.is_odd_letter, self.fail])(x='a'), True)
        self.assertRaises(AssertionError, any([self.is_odd_letter,self.fail]),\
            'b')

    def test_either(self):
        """either should return function returning True if either component True"""
        odd_vowel = either(self.is_odd_letter, self.is_vowel)
        self.assertEqual(map(odd_vowel, 'abcu'), [True,False,True,True])
        #should short-circuit
        self.assertEqual(either(self.is_odd_letter, self.fail)(x='a'), True)
        self.assertRaises(AssertionError, \
            either(self.is_odd_letter, self.fail), 'b')
        
    def test_negate(self):
        """negate should return True if no component True"""
        self.assertEqual(negate([self.is_odd_letter,self.is_vowel], 'a'), False)
        self.assertEqual(negate([self.is_odd_letter,self.is_vowel], 'b'), True)
        self.assertEqual(negate([self.is_odd_letter,self.is_vowel], 'c'), False)
        self.assertEqual(negate([self.is_odd_letter,self.is_vowel], 'u'), False)
        #should short-circuit after first True
        self.assertEqual(negate([self.is_odd_letter, self.fail], x='a'), False)
        self.assertRaises(AssertionError, \
            negate, [self.is_odd_letter, self.fail], 'b')
        
    def test_none(self):
        """none should return a function returning True if no component True"""
        odd_vowel = none([self.is_odd_letter, self.is_vowel])
        self.assertEqual(odd_vowel('a'), False)
        self.assertEqual(map(odd_vowel, 'abceu'), [False,True] + [False]*3)
        odd = none([self.is_odd_number, self.is_small])
        self.assertEqual(map(odd, range(5)), [False]*4+[True])
        #should short-circuit after first True
        self.assertEqual(none([self.is_odd_letter, self.fail])(x='a'), False)
        self.assertRaises(AssertionError, none([self.is_odd_letter,self.fail]),\
            'b')

    def test_neither(self):
        """neither should return function returning True if each component False"""
        odd_vowel = neither(self.is_odd_letter, self.is_vowel)
        self.assertEqual(map(odd_vowel, 'abcu'), [False,True,False,False])
        #should short-circuit
        self.assertEqual(neither(self.is_odd_letter, self.fail)(x='a'), False)
        self.assertRaises(AssertionError, \
            neither(self.is_odd_letter, self.fail), 'b')
 
    def test_compose(self):
        """compose should return function returning f(g(x))"""
        ds = compose(self.double, self.minusone)
        sd = compose(self.minusone, self.double)
        self.assertEqual(ds(5), 8)
        self.assertEqual(sd(x=5), 9)
        #check that it works when arg lists are different
        commafy = compose(','.join, list)
        self.assertEqual(commafy('abc'), 'a,b,c')
        self.assertEqual(commafy(''), '')
        self.assertEqual(commafy('a'), 'a')

    def test_compose_many(self):
        """compose_many should return composition of all args"""
        from numpy import arange
        def to_strings(x):
            return map(str, x)
        printable_range = compose_many(''.join, to_strings, range)
        printable_arange = compose_many(''.join, to_strings, arange)
        
        self.assertEqual(printable_range(3), '012')
        self.assertEqual(printable_range(0), '')
        self.assertEqual(printable_range(5), '01234')

        self.assertEqual(printable_arange(stop=51, start=10, step=10),
            '1020304050')

    def test_identity(self):
        """identity should return x"""
        for i in ['a', 'abc', None, '', [], [1], 1, 2**50, 0.3e-50, {'a':3}]:
            assert identity(i) is i
            
    def test_has_field(self):
        """has_field should return True if specified field exists."""
        x = has_x(1)
        y = has_y(1)
        check_x = has_field('x')
        self.assertEqual(check_x(x), True)
        self.assertEqual(check_x(y), False)
        check_y = has_field('y')
        self.assertEqual(check_y(x), False)
        self.assertEqual(check_y(y), True)
        del y.y
        self.assertEqual(check_y(y), False)
        y.x = 3
        self.assertEqual(check_x(y), True)

    def test_extract_field(self):
        """extract_field should apply constructor to field, or return None"""
        num = has_x('1')
        alpha = has_x('x')
        y = has_y('1')
        extractor = extract_field('x')
        self.assertEqual(extractor(num), '1')
        self.assertEqual(extractor(alpha), 'x')
        self.assertEqual(extractor(y), None)

        int_extractor = extract_field('x', int)
        self.assertEqual(int_extractor(num), 1)
        self.assertEqual(int_extractor(alpha), None)
        self.assertEqual(int_extractor(y), None)

    def test_test_field(self):
        """test_field should return boolean result of applying constructor"""
        num = has_x('5')
        alpha = has_x('x')
        zero = has_x(0)
        y = has_y('5')
        
        tester = test_field('x')
        self.assertEqual(tester(num), True)
        self.assertEqual(tester(alpha), True)
        self.assertEqual(tester(y), False)

        int_tester = test_field('x', int)
        self.assertEqual(int_tester(num), True)
        self.assertEqual(int_tester(alpha), False)
        self.assertEqual(int_tester(y), False)
        self.assertEqual(int_tester(zero), False)

    def test_index(self):
        """index should index objects by specified field or identity"""
        num = has_x(5)
        let = has_x('5')
        zer = has_x('0')
        non = has_x(None)
        y = has_y(3)
        items = [num, let, zer, non, y]
        duplicates = items * 3

        basic_indexer = index()
        i = basic_indexer(items)
        self.assertEqual(i, {num:[num], let:[let], zer:[zer], non:[non], y:[y]})
        #test reusability
        i = basic_indexer([3,3,4])
        self.assertEqual(i, {3:[3, 3], 4:[4]})
        #test duplicates
        d = basic_indexer(duplicates)
        self.assertEqual(d, {num:[num]*3, let:[let]*3, zer:[zer]*3, \
            non:[non]*3, y:[y]*3})
        #test with constructor
        str_indexer = index(str)
        i = str_indexer(items)
        self.assertEqual(i, {'5':[num,let], '0':[zer], 'None':[non], '3':[y]})
        #test order correct in duplicates
        i = str_indexer(duplicates)
        self.assertEqual(i, {'5':[num,let,num,let,num,let], '0':[zer,zer,zer],
            'None':[non,non,non], '3':[y,y,y]})
        #test with squashing
        overwriter = index(str, overwrite=True)
        i = overwriter(duplicates)
        self.assertEqual(i, {'5':let, '0':zer, 'None':non, '3':y})

    def test_test_container(self):
        """test_container should return True or False in a typesafe way."""
        test_dict = test_container({'a':1})
        test_list = test_container([1,2,3])
        test_str = test_container('438hfanvr438')

        for item in (1, 2, 3):
            assert test_list(item)
            assert not test_dict(item)
            assert not test_str(item)

        assert test_dict('a')
        assert not test_list('a')
        assert test_str('a')

        for item in ('4', 'h', 'fan'):
            assert not test_dict(item)
            assert not test_list(item)
            assert test_str(item)

        for item in (['x','y'],{},{'a':3},'@#@',('a','b'),None,False):
            assert not test_dict(item)
            assert not test_list(item)
            assert not test_str(item)

class SequenceFunctionsTests(TestCase):
    """Tests of standalone functions for dealing with sequences."""
    def test_per_shortest(self):
        """per_shortest should divide by min(len(x), len(y))"""
        self.assertEqual(per_shortest(20, 'aaaaaa', 'bbbb'), 5)
        self.assertEqual(per_shortest(20, 'aaaaaa', 'b'), 20)
        self.assertEqual(per_shortest(20, 'a', 'bbbbb'), 20)
        self.assertEqual(per_shortest(20, '', 'b'), 0)
        self.assertEqual(per_shortest(20, '', ''), 0)
        #check that it does it in floating-point
        self.assertEqual(per_shortest(1, 'aaaaaa', 'bbbb'), 0.25)
        #check that it raises TypeError on non-seq
        self.assertRaises(TypeError, per_shortest, 1, 2, 3)

    def test_per_longest(self):
        """per_longest should divide by max(len(x), len(y))"""
        self.assertEqual(per_longest(20, 'aaaaaa', 'bbbb'), 20/6.0)
        self.assertEqual(per_longest(20, 'aaaaaa', 'b'), 20/6.0)
        self.assertEqual(per_longest(20, 'a', 'bbbbb'), 20/5.0)
        self.assertEqual(per_longest(20, '', 'b'), 20)
        self.assertEqual(per_longest(20, '', ''), 0)
        #check that it does it in floating-point
        self.assertEqual(per_longest(1, 'aaaaaa', 'bbbb'), 1/6.0)
        #check that it raises TypeError on non-seq
        self.assertRaises(TypeError, per_longest, 1, 2, 3)

    def test_for_seq(self):
        """for_seq should return the correct function"""
        is_eq = lambda x,y: x == y
        is_ne = lambda x,y: x != y
        lt_5 =  lambda x,y: x + y < 5
        diff =  lambda x,y: x - y

        sumsq = lambda x: sum([i*i for i in x])

        long_norm = lambda s, x, y: (s + 0.0) / max(len(x), len(y))
        times_two = lambda s, x, y: 2*s

        empty = []
        s1 = [1,2,3,4,5]
        s2 = [1,3,2,4,5]
        s3 = [1,1,1,1,1]
        s4 = [5,5,5,5,5]
        s5 = [3,3,3,3,3]
        short = [1]

        #test behavior with default aggregator and normalizer
        f = for_seq(is_eq)
        self.assertAlmostEqual(f(s1, s1), 1.0)
        self.assertAlmostEqual(f(s1, short), 1.0)
        self.assertAlmostEqual(f(short, s1), 1.0)
        self.assertAlmostEqual(f(short, s4), 0.0)
        self.assertAlmostEqual(f(s4, short), 0.0)
        self.assertAlmostEqual(f(s1,s2), 0.6)
        
        f = for_seq(is_ne)
        self.assertAlmostEqual(f(s1, s1), 0.0)
        self.assertAlmostEqual(f(s1, short), 0.0)
        self.assertAlmostEqual(f(short, s1), 0.0)
        self.assertAlmostEqual(f(short, s4), 1.0)
        self.assertAlmostEqual(f(s4, short), 1.0)
        self.assertAlmostEqual(f(s1, s2), 0.4)
         
        f = for_seq(lt_5)
        self.assertAlmostEqual(f(s3,s3), 1.0)
        self.assertAlmostEqual(f(s3,s4), 0.0)
        self.assertAlmostEqual(f(s2,s3), 0.6)

        f = for_seq(diff)
        self.assertAlmostEqual(f(s1,s1), 0.0)
        self.assertAlmostEqual(f(s4,s1), 2.0)
        self.assertAlmostEqual(f(s1,s4), -2.0)

        #test behavior with different aggregator
        f = for_seq(diff)
        self.assertAlmostEqual(f(s1,s5), 0)
        f = for_seq(diff, aggregator=sum)
        self.assertAlmostEqual(f(s1,s5), 0)
        f = for_seq(diff, aggregator=sumsq)
        self.assertAlmostEqual(f(s1,s5), 2.0)

        #test behavior with different normalizer
        f = for_seq(diff, aggregator=sumsq, normalizer=None)
        self.assertAlmostEqual(f(s1,s5), 10)
        f = for_seq(diff, aggregator=sumsq)
        self.assertAlmostEqual(f(s1,s5), 2.0)
        f = for_seq(diff, aggregator=sumsq, normalizer=times_two)
        self.assertAlmostEqual(f(s1,s5), 20)
        f = for_seq(diff, aggregator=sumsq)
        self.assertAlmostEqual(f(s5,short), 4)
        f = for_seq(diff, aggregator=sumsq, normalizer=long_norm)
        self.assertAlmostEqual(f(s5,short), 0.8)
        


class Filter_Criteria_Tests(TestCase):
    """Tests of standalone functions used as filter criteria"""

    def test_trans_except(self):
        """trans_except should return trans table mapping non-good chars to x"""
        a = trans_except('Aa', '-')
        none = trans_except('', '*')
        some = trans_except('zxcvbnm,.zxcvbnm,.', 'V')

        self.assertEqual('abcABA'.translate(a), 'a--A-A')
        self.assertEqual(''.translate(a), '')
        self.assertEqual('12345678'.translate(a), '--------')

        self.assertEqual(''.translate(none), '')
        self.assertEqual('abcdeEFGHI12345&*(!@'.translate(none), '*'*20)

        self.assertEqual('qazwsxedcrfv'.translate(some),'VVzVVxVVcVVv') 

    def test_trans_all(self):
        """trans_all should return trans table mapping all bad chars to x"""
        a = trans_all('Aa', '-')
        none = trans_all('', '*')
        some = trans_all('zxcvbnm,.zxcvbnm,.', 'V')

        self.assertEqual('abcABA'.translate(a), '-bc-B-')
        self.assertEqual(''.translate(a), '')
        self.assertEqual('12345678'.translate(a), '12345678')

        self.assertEqual(''.translate(none), '')
        self.assertEqual('abcdeEFGHI12345&*(!@'.translate(none), \
            'abcdeEFGHI12345&*(!@')

        self.assertEqual('qazwsxedcrfv'.translate(some),'qaVwsVedVrfV') 

    def test_make_trans(self):
        """make_trans should return trans table mapping chars to default"""
        a = make_trans()
        self.assertEqual('abc123'.translate(a), 'abc123')
        a = make_trans('a', 'x')
        self.assertEqual('abc123'.translate(a), 'xbc123')
        a = make_trans('ac', 'xa')
        self.assertEqual('abc123'.translate(a), 'xba123')
        a = make_trans('ac', 'xa', '.')
        self.assertEqual('abc123'.translate(a), 'x.a...')
        self.assertRaises(ValueError, make_trans, 'ac', 'xa', 'av')

    def test_find_any(self):
        """find_any should be True if one of the words is in the string"""
        
        f = find_any('ab')
        self.assertEqual(f(''),0) #empty
        self.assertRaises(AttributeError,f,None) # none
        self.assertEqual(f('cde'),0) #none of the elements
        self.assertEqual(f('axxx'),1) #one of the elements
        self.assertEqual(f('bxxx'),1) #one of the elements
        self.assertEqual(f('axxxb'),1) #all elements
        self.assertEqual(f('aaaa'),1) #repeated element

        # works on any sequence
        f = find_any(['foo','bar'])
        self.assertEqual(f("joe"),0)
        self.assertEqual(f("only foo"),1)
        self.assertEqual(f("bar and foo"),1)
        
        # does NOT work on numbers
        
    def test_find_no(self):
        """find_no should be True if none of the words in the string"""
       
        f = find_no('ab')
        self.assertEqual(f(''),1) #empty
        self.assertRaises(AttributeError,f,None) # none
        self.assertEqual(f('cde'),1) #none of the elements
        self.assertEqual(f('axxx'),0) #one of the elements
        self.assertEqual(f('bxxx'),0) #one of the elements
        self.assertEqual(f('axxxb'),0) #all elements
        self.assertEqual(f('aaaa'),0) #repeated element

        # works on any sequence
        f = find_no(['foo','bar'])
        self.assertEqual(f("joe"),1)
        self.assertEqual(f("only foo"),0)
        self.assertEqual(f("bar and foo"),0)
        
        # does NOT work on numbers

    def test_find_all(self):
        """find_all should be True if all words appear in the string"""

        f = find_all('ab')
        self.assertEqual(f(''),0) #empty
        self.assertRaises(AttributeError,f,None) # none
        self.assertEqual(f('cde'),0) #none of the elements
        self.assertEqual(f('axxx'),0) #one of the elements
        self.assertEqual(f('bxxx'),0) #one of the elements
        self.assertEqual(f('axxxb'),1) #all elements
        self.assertEqual(f('aaaa'),0) #repeated element

        # works on any sequence
        f = find_all(['foo','bar'])
        self.assertEqual(f("joe"),0)
        self.assertEqual(f("only foo"),0)
        self.assertEqual(f("bar and foo"),1)
        
        # does NOT work on numbers

    def test_keep_chars(self):
        """keep_chars returns a string containing only chars in keep"""
        f = keep_chars('ab c3*[')
        self.assertEqual(f(''),'') #empty
        self.assertRaises(AttributeError,f,None) #None
        
        #one character, case sensitive
        self.assertEqual(f('b'),'b')
        self.assertEqual(f('g'),'')
        self.assertEqual(f('xyz123'),'3')
        self.assertEqual(f('xyz  123'),'  3')
        
        #more characters, case sensitive
        self.assertEqual(f('kjbwherzcagebcujrkcs'),'bcabcc')
        self.assertEqual(f('f[ffff*ff*fff3fff'),'[**3')

        # case insensitive
        f = keep_chars('AbC',False)
        self.assertEqual(f('abcdef'),'abc')
        self.assertEqual(f('ABCDEF'),'ABC')
        self.assertEqual(f('aBcDeF'),'aBc')

    def test_exclude_chars(self):
        """exclude_chars returns string containing only chars not in exclude"""
        
        f = exclude_chars('ab c3*[')
        self.assertEqual(f(''),'') #empty
        self.assertRaises(AttributeError,f,None) #None
        
        #one character, case sensitive
        self.assertEqual(f('b'),'')
        self.assertEqual(f('g'),'g')
        self.assertEqual(f('xyz123'),'xyz12')
        self.assertEqual(f('xyz  123'),'xyz12')
        
        #more characters, case sensitive
        self.assertEqual(f('axxxbxxxcxxx'),'xxxxxxxxx')

        # case insensitive
        f = exclude_chars('AbC',False)
        self.assertEqual(f('abcdef'),'def')
        self.assertEqual(f('ABCDEF'),'DEF')
        self.assertEqual(f('aBcDeF'),'DeF')

    def test_reorder(self):
        """reorder should always use the same order when invoked"""
        list_test = reorder([3,2,1])
        dict_test = reorder(['x','y','z'])
        multi_test = reorder([3,2,2])
        null_test = reorder([])

        first_seq = 'abcde'
        second_seq = [3,4,5,6,7]
        empty_list = []
        empty_dict = {}
        full_dict = {'a':3, 'c':5, 'x':'abc','y':'234','z':'qaz'}

        for i in (first_seq, second_seq, empty_list, empty_dict):
            self.assertEqual(null_test(i), [])

        self.assertEqual(list_test(first_seq), ['d','c','b'])
        self.assertEqual(list_test(second_seq), [6,5,4])
        self.assertEqual(multi_test(first_seq), ['d','c','c'])
        self.assertEqual(dict_test(full_dict), ['abc','234','qaz'])

        self.assertRaises(KeyError, dict_test, empty_dict)
        self.assertRaises(IndexError, list_test, empty_list)
        
    def test_reorder_inplace(self):
        """reorder_inplace should replace object's data with new order"""
        attr_test = reorder_inplace([3,2,1], 'Data')
        obj_test = reorder_inplace([3,2,2])

        seq = [3,4,5,6,7]

        class obj(object):
            pass

        o = obj()
        o.XYZ = [9, 7, 5]
        o.Data = ['a','b','c','d','e']
        orig_data = o.Data

        self.assertEqual(obj_test(seq), [6,5,5])
        self.assertEqual(seq, [6,5,5])
        assert attr_test(o) is o
        self.assertEqual(o.XYZ, [9,7,5])
        self.assertEqual(o.Data, ['d','c','b'])
        assert orig_data is o.Data
    
    def test_float_from_string(self):
        """float_from_string should ignore funny chars"""
        ffs = float_from_string
        self.assertEqual(ffs('3.5'), 3.5)
        self.assertEqual(ffs('     -3.45e-10   '), float('     -3.45e-10   '))
        self.assertEqual(ffs('jsdjhsdf[]()0.001IVUNZSDFl]]['), 0.001)

    def test_first_index(self):
        """first_index should return index of first occurrence where f(s)"""
        vowels = 'aeiou'
        is_vowel = lambda x: x in vowels
        s1 = 'ebcua'
        s2 = 'bcbae'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(first_index(s1, is_vowel), 0)
        self.assertEqual(first_index(s2, is_vowel), 3)
        self.assertEqual(first_index(s3, is_vowel), None)
        self.assertEqual(first_index(s4, is_vowel), None)
          
    def test_last_index(self):
        """last_index should return index of last occurrence where f(s)"""
        vowels = 'aeiou'
        is_vowel = lambda x: x in vowels
        s1 = 'ebcua'
        s2 = 'bcbaef'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(last_index(s1, is_vowel), 4)
        self.assertEqual(last_index(s2, is_vowel), 4)
        self.assertEqual(last_index(s3, is_vowel), None)
        self.assertEqual(last_index(s4, is_vowel), None)

    def test_first_index_in_set(self):
        """first_index_in_set should return index of first occurrence """
        vowels = 'aeiou'
        s1 = 'ebcua'
        s2 = 'bcbae'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(first_index_in_set(s1, vowels), 0)
        self.assertEqual(first_index_in_set(s2, vowels), 3)
        self.assertEqual(first_index_in_set(s3, vowels), None)
        self.assertEqual(first_index_in_set(s4, vowels), None)
          
    def test_last_index_in_set(self):
        """last_index_in_set should return index of last occurrence"""
        vowels = 'aeiou'
        s1 = 'ebcua'
        s2 = 'bcbaef'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(last_index_in_set(s1, vowels), 4)
        self.assertEqual(last_index_in_set(s2, vowels), 4)
        self.assertEqual(last_index_in_set(s3, vowels), None)
        self.assertEqual(last_index_in_set(s4, vowels), None)

    def test_first_index_not_in_set(self):
        """first_index_not_in_set should return index of first occurrence """
        vowels = 'aeiou'
        s1 = 'ebcua'
        s2 = 'bcbae'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(first_index_not_in_set(s1, vowels), 1)
        self.assertEqual(first_index_not_in_set(s2, vowels), 0)
        self.assertEqual(first_index_not_in_set(s3, vowels), None)
        self.assertEqual(first_index_not_in_set(s4, vowels), 0)
          
    def test_last_index_not_in_set(self):
        """last_index_not_in_set should return index of last occurrence"""
        vowels = 'aeiou'
        s1 = 'ebcua'
        s2 = 'bcbaef'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(last_index_not_in_set(s1, vowels), 2)
        self.assertEqual(last_index_not_in_set(s2, vowels), 5)
        self.assertEqual(last_index_not_in_set(s3, vowels), None)
        self.assertEqual(last_index_not_in_set(s4, vowels), 2)

    def test_first(self):
        """first should return first occurrence where f(s)"""
        vowels = 'aeiou'
        is_vowel = lambda x: x in vowels
        s1 = 'ebcua'
        s2 = 'bcbae'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(first(s1, is_vowel), 'e')
        self.assertEqual(first(s2, is_vowel), 'a')
        self.assertEqual(first(s3, is_vowel), None)
        self.assertEqual(first(s4, is_vowel), None)
          
    def test_last(self):
        """last should return last occurrence where f(s)"""
        vowels = 'aeiou'
        is_vowel = lambda x: x in vowels
        s1 = 'ebcua'
        s2 = 'bcbaef'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(last(s1, is_vowel), 'a')
        self.assertEqual(last(s2, is_vowel), 'e')
        self.assertEqual(last(s3, is_vowel), None)
        self.assertEqual(last(s4, is_vowel), None)

    def test_first_in_set(self):
        """first_in_set should return first occurrence """
        vowels = 'aeiou'
        s1 = 'ebcua'
        s2 = 'bcbae'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(first_in_set(s1, vowels), 'e')
        self.assertEqual(first_in_set(s2, vowels), 'a')
        self.assertEqual(first_in_set(s3, vowels), None)
        self.assertEqual(first_in_set(s4, vowels), None)
          
    def test_last_in_set(self):
        """last_in_set should return last occurrence"""
        vowels = 'aeiou'
        s1 = 'ebcua'
        s2 = 'bcbaef'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(last_in_set(s1, vowels), 'a')
        self.assertEqual(last_in_set(s2, vowels), 'e')
        self.assertEqual(last_in_set(s3, vowels), None)
        self.assertEqual(last_in_set(s4, vowels), None)

    def test_first_not_in_set(self):
        """first_not_in_set should return first occurrence """
        vowels = 'aeiou'
        s1 = 'ebcua'
        s2 = 'bcbae'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(first_not_in_set(s1, vowels), 'b')
        self.assertEqual(first_not_in_set(s2, vowels), 'b')
        self.assertEqual(first_not_in_set(s3, vowels), None)
        self.assertEqual(first_not_in_set(s4, vowels), 'c')
          
    def test_last_not_in_set(self):
        """last_not_in_set should return last occurrence"""
        vowels = 'aeiou'
        s1 = 'ebcua'
        s2 = 'bcbaef'
        s3 = ''
        s4 = 'cbd'
        self.assertEqual(last_not_in_set(s1, vowels), 'c')
        self.assertEqual(last_not_in_set(s2, vowels), 'f')
        self.assertEqual(last_not_in_set(s3, vowels), None)
        self.assertEqual(last_not_in_set(s4, vowels), 'd')

    def test_perm(self):
        """perm should return correct permutations"""
        self.assertEqual(list(perm('abc')), ['abc','acb','bac','bca','cab','cba']) 

    def test_comb(self):
        """comb should return correct combinations"""
        self.assertEqual(list(comb(range(5), 0)),
                [])
        self.assertEqual(list(comb(range(5), 1)),
                [[0], [1], [2], [3], [4]])
        self.assertEqual(list(comb(range(5), 2)),
                [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4], [2, 3],
                [2, 4], [3, 4]])
        self.assertEqual(list(comb(range(5), 3)),
               [[0, 1, 2], [0, 1, 3], [0, 1, 4], [0, 2, 3], [0, 2, 4], [0, 3, 4],
                [1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])
        self.assertEqual(list(comb(range(5), 4)),
                [[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 3, 4], [0, 2, 3, 4], [1, 2, 3, 4]])
        self.assertEqual(list(comb(range(5), 5)),
                [[0, 1, 2, 3, 4]])


    def test_cross_comb(self):
        """cross_comb should produce correct combinations"""
        v1 = range(2)
        v2 = range(3)
        v3 = list('abc')
        vv1 = ([e] for e in v1)
        v1_x_v2 = [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2]]
        v1v2v3 = [[0, 0, 'a'], [0, 0, 'b'], [0, 0, 'c'], [0, 1, 'a'],
                  [0, 1, 'b'], [0, 1, 'c'], [0, 2, 'a'], [0, 2, 'b'],
                  [0, 2, 'c'], [1, 0, 'a'], [1, 0, 'b'], [1, 0, 'c'],
                  [1, 1, 'a'], [1, 1, 'b'], [1, 1, 'c'], [1, 2, 'a'],
                  [1, 2, 'b'], [1, 2, 'c']]
        self.assertEqual(list( _increment_comb(vv1, v2)), v1_x_v2)
        self.assertEqual(list( cross_comb([v1, v2])), v1_x_v2)
        self.assertEqual(list(cross_comb([v1, v2, v3])), v1v2v3)



 #run tests if invoked from the commandline
if __name__ == '__main__':
    main()
