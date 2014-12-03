# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.parse.record import (DelimitedSplitter,
                                GenericRecord, MappedRecord, TypeSetter,
                                list_adder, dict_adder,
                                LineOrientedConstructor, int_setter,
                                bool_setter, string_and_strip, FieldWrapper,
                                StrictFieldWrapper, raise_unknown_field,
                                FieldMorpher)
from skbio.io import FieldError


class recordsTests(TestCase):

    """Tests of top-level functionality in records."""

    def test_string_and_strip(self):
        """string_and_strip should convert all items to strings and strip them
        """
        self.assertEqual(string_and_strip(), [])
        self.assertEqual(string_and_strip('\t', ' ', '\n\t'), ['', '', ''])
        self.assertEqual(string_and_strip('\ta\tb', 3, '   cde   e', None),
                         ['a\tb', '3', 'cde   e', 'None'])

    def test_raise_unknown_field(self):
        """raise_unknown_field should always raise FieldError"""
        self.assertRaises(FieldError, raise_unknown_field, 'xyz', 123)


class DelimitedSplitterTests(TestCase):

    """Tests of the DelimitedSplitter factory function."""

    def test_parsers(self):
        """DelimitedSplitter should return function with correct behavior"""
        empty = DelimitedSplitter()
        space = DelimitedSplitter(None)
        semicolon = DelimitedSplitter(';')
        twosplits = DelimitedSplitter(';', 2)
        allsplits = DelimitedSplitter(';', None)
        lastone = DelimitedSplitter(';', -1)
        lasttwo = DelimitedSplitter(';', -2)

        self.assertEqual(empty('a   b  c'), ['a', 'b  c'])
        self.assertEqual(empty('abc'), ['abc'])
        self.assertEqual(empty('   '), [])

        self.assertEqual(empty('a  b  c'), space('a  b  c'))
        self.assertEqual(semicolon('  a  ; b   ;  c  d'), ['a', 'b   ;  c  d'])
        self.assertEqual(twosplits('  a  ; b   ;  c  d'), ['a', 'b', 'c  d'])
        self.assertEqual(allsplits(' a ;  b  ; c;;d;e  ;'),
                         ['a', 'b', 'c', '', 'd', 'e', ''])
        self.assertEqual(lastone(' a ;  b  ; c;;d;e  ;'),
                         ['a ;  b  ; c;;d;e', ''])
        self.assertEqual(lasttwo(' a ;  b  ; c;;d;e  ;'),
                         ['a ;  b  ; c;;d', 'e', ''])
        self.assertEqual(lasttwo(''), [])
        self.assertEqual(lasttwo('x'), ['x'])
        self.assertEqual(lasttwo('x;'), ['x', ''])


class GenericRecordTests(TestCase):

    """Tests of the GenericRecord class"""
    class gr(GenericRecord):
        Required = {'a': 'x', 'b': [], 'c': {}}

    def test_init(self):
        """GenericRecord init should work OK empty or with data"""
        self.assertEqual(GenericRecord(), {})
        self.assertEqual(GenericRecord({'a': 1}), {'a': 1})
        assert isinstance(GenericRecord(), GenericRecord)

    def test_init_subclass(self):
        """GenericRecord subclass init should include required data"""
        self.assertEqual(self.gr(), {'a': 'x', 'b': [], 'c': {}})
        self.assertEqual(self.gr({'a': []}), {'a': [], 'b': [], 'c': {}})
        assert isinstance(self.gr(), self.gr)
        assert isinstance(self.gr(), GenericRecord)

    def test_delitem(self):
        """GenericRecord delitem should fail if item required"""
        g = self.gr()
        g['d'] = 3
        self.assertEqual(g, {'a': 'x', 'b': [], 'c': {}, 'd': 3})
        del g['d']
        self.assertEqual(g, {'a': 'x', 'b': [], 'c': {}})
        self.assertRaises(AttributeError, g.__delitem__, 'a')
        g['c'][3] = 4
        self.assertEqual(g['c'], {3: 4})

    def test_copy(self):
        """GenericRecord copy should include attributes and set correct class
        """
        g = self.gr()
        g['a'] = 'abc'
        g.X = 'y'
        h = g.copy()
        self.assertEqual(g, h)
        assert isinstance(h, self.gr)
        self.assertEqual(h.X, 'y')
        self.assertEqual(h, {'a': 'abc', 'b': [], 'c': {}})


class MappedRecordTests(TestCase):

    """Tests of the MappedRecord class"""

    def setUp(self):
        """Define a few standard MappedRecords"""
        self.empty = MappedRecord()
        self.single = MappedRecord({'a': 3})
        self.several = MappedRecord(a=4, b=5, c='a', d=[1, 2, 3])

    def test_init_empty(self):
        """MappedRecord empty init should work OK"""
        g = MappedRecord()
        self.assertEqual(g, {})

    def test_init_data(self):
        """MappedRecord should work like normal dict init"""
        exp = {'a': 3, 'b': 4}
        self.assertEqual(MappedRecord({'a': 3, 'b': 4}), exp)
        self.assertEqual(MappedRecord(a=3, b=4), exp)
        self.assertEqual(MappedRecord([['a', 3], ['b', 4]]), exp)

    def test_init_subclass(self):
        """MappedRecord subclasses should behave as expected"""
        class rec(MappedRecord):
            Required = {'a': {}, 'b': 'xyz', 'c': 3}
            Aliases = {'B': 'b'}

        r = rec()
        self.assertEqual(r, {'a': {}, 'b': 'xyz', 'c': 3})
        # test that subclassing is correct
        s = r.copy()
        assert isinstance(s, rec)
        # test Aliases
        s.B = 0
        self.assertEqual(s, {'a': {}, 'b': 0, 'c': 3})
        # test Required
        try:
            del s.B
        except AttributeError:
            pass
        else:
            raise AssertionError("Subclass failed to catch requirement")

    def test_getattr(self):
        """MappedRecord getattr should look in dict after real attrs"""
        s = self.several
        self.assertEqual(s.Aliases, {})
        self.assertEqual(s.a, 4)
        self.assertEqual(s.d, [1, 2, 3])
        for key in s:
            self.assertEqual(getattr(s, key), s[key])
        assert 'xyz' not in s
        self.assertEqual(s.xyz, None)
        self.assertEqual(s['xyz'], None)
        s.Aliases = {'xyz': 'a'}
        self.assertEqual(s['xyz'], 4)

    def test_setattr(self):
        """MappedRecord setattr should add to dict"""
        s = self.single
        # check that we haven't screwed up normal attribute setting
        assert 'Aliases' not in s
        s.Aliases = {'x': 'y'}
        assert 'Aliases' not in s
        self.assertEqual(s.Aliases, {'x': 'y'})
        s.x = 5
        assert 'x' in s
        self.assertEqual(s['x'], 5)
        self.assertEqual(s.x, 5)
        s.Aliases = {'XYZ': 'b'}
        s.XYZ = 3
        self.assertEqual(s.b, 3)

    def test_delattr(self):
        """MappedRecord delattr should work for 'normal' and other attributes
        """
        s = self.single
        s.__dict__['x'] = 'y'
        assert 'x' not in s
        self.assertEqual(s.x, 'y')
        del s.x
        self.assertEqual(s.x, None)
        self.assertEqual(s, {'a': 3})
        # try it for an internal attribute: check it doesn't delete anything
        # else
        s.b = 4
        self.assertEqual(s, {'a': 3, 'b': 4})
        del s.a
        self.assertEqual(s, {'b': 4})
        del s.abc
        self.assertEqual(s, {'b': 4})
        s.Required = {'b': True}
        try:
            del s.b
        except AttributeError:
            pass
        else:
            raise AssertionError("Allowed deletion of required attribute""")
        s.a = 3
        self.assertEqual(s.a, 3)
        s.Aliases = {'xyz': 'a'}
        del s.xyz
        self.assertEqual(s.a, None)

    def test_getitem(self):
        """MappedRecord getitem should work only for keys, not attributes"""
        s = self.single
        self.assertEqual(s['Required'], None)
        self.assertEqual(s['a'], 3)
        self.assertEqual(s['xyz'], None)
        self.assertEqual(s[list('abc')], None)
        s.Aliases = {'xyz': 'a'}
        self.assertEqual(s['xyz'], 3)

    def test_setitem(self):
        """MappedRecord setitem should work only for keys, not attributes"""
        s = self.single
        s['Required'] = None
        self.assertEqual(s, {'a': 3, 'Required': None})
        self.assertEqual(s.Required, {})
        self.assertNotEqual(s.Required, None)
        s['c'] = 5
        self.assertEqual(s, {'a': 3, 'c': 5, 'Required': None})
        # still not allowed unhashable objects as keys
        self.assertRaises(TypeError, s.__setitem__, range(3))
        s.Aliases = {'C': 'c'}
        s['C'] = 3
        self.assertEqual(s, {'a': 3, 'c': 3, 'Required': None})

    def test_delitem(self):
        """MappedRecord delitem should only work for keys, not attributes"""
        s = self.single
        del s['Required']
        self.assertEqual(s.Required, {})
        s.Required = {'a': True}
        try:
            del s['a']
        except AttributeError:
            pass
        else:
            raise AssertionError("Allowed deletion of required item")
        s.Aliases = {'B': 'b'}
        s.b = 5
        self.assertEqual(s.b, 5)
        del s.B
        self.assertEqual(s.b, None)

    def test_contains(self):
        """MappedRecord contains should use aliases, but not apply to attrs"""
        s = self.single
        assert 'a' in s
        assert 'b' not in s
        s.b = 5
        assert 'b' in s
        assert 'Required' not in s
        assert 'A' not in s
        s.Aliases = {'A': 'a'}
        assert 'A' in s

    def test_get(self):
        """MappedRecord get should be typesafe against unhashables"""
        s = self.single
        self.assertEqual(s.get(1, 6), 6)
        self.assertEqual(s.get('a', 'xyz'), 3)
        self.assertEqual(s.get('ABC', 'xyz'), 'xyz')
        s.Aliases = {'ABC': 'a'}
        self.assertEqual(s.get('ABC', 'xyz'), 3)
        self.assertEqual(s.get([1, 2, 3], 'x'), 'x')

    def test_setdefault(self):
        """MappedRecord setdefault should not be typesafe against unhashables
        """
        s = self.single
        x = s.setdefault('X', 'xyz')
        self.assertEqual(x, 'xyz')
        self.assertEqual(s, {'a': 3, 'X': 'xyz'})
        self.assertRaises(TypeError, s.setdefault, ['a', 'b'], 'xyz')

    def test_update(self):
        """MappedRecord update should transparently convert keys"""
        s = self.single
        s.b = 999
        s.Aliases = {'XYZ': 'x', 'ABC': 'a'}
        d = {'ABC': 111, 'CVB': 222}
        s.update(d)
        self.assertEqual(s, {'a': 111, 'b': 999, 'CVB': 222})

    def test_copy(self):
        """MappedRecord copy should return correct class"""
        s = self.single
        t = s.copy()
        assert isinstance(t, MappedRecord)
        s.Aliases = {'XYZ': 'x'}
        u = s.copy()
        u.Aliases['ABC'] = 'a'
        self.assertEqual(s.Aliases, {'XYZ': 'x'})
        self.assertEqual(t.Aliases, {})
        self.assertEqual(u.Aliases, {'XYZ': 'x', 'ABC': 'a'})

    def test_subclass(self):
        """MappedRecord subclassing should work correctly"""
        class ret3(MappedRecord):
            DefaultValue = 3
            ClassData = 'xyz'

        x = ret3({'ABC': 777, 'DEF': '999'})
        self.assertEqual(x.ZZZ, 3)
        self.assertEqual(x.ABC, 777)
        self.assertEqual(x.DEF, '999')
        self.assertEqual(x.ClassData, 'xyz')
        x.ZZZ = 6
        self.assertEqual(x.ZZZ, 6)
        self.assertEqual(x.ZZ, 3)
        x.ClassData = 'qwe'
        self.assertEqual(x.ClassData, 'qwe')
        self.assertEqual(ret3.ClassData, 'xyz')

    def test_DefaultValue(self):
        """MappedRecord DefaultValue should give new copy when requested"""
        class m(MappedRecord):
            DefaultValue = []

        a = m()
        b = m()
        assert a['abc'] is not b['abc']
        assert a['abc'] == b['abc']


class dummy(object):

    """Do-nothing class whose attributes can be freely abused."""
    pass


class TypeSetterTests(TestCase):

    """Tests of the TypeSetter class"""

    def test_setter_empty(self):
        """TypeSetter should set attrs to vals on empty init"""
        d = dummy()
        ident = TypeSetter()
        ident(d, 'x', 'abc')
        self.assertEqual(d.x, 'abc')
        ident(d, 'y', 3)
        self.assertEqual(d.y, 3)
        ident(d, 'x', 2)
        self.assertEqual(d.x, 2)

    def test_setter_typed(self):
        """TypeSetter should set attrs to constructor(val) when specified"""
        d = dummy()
        i = TypeSetter(int)
        i(d, 'zz', 3)
        self.assertEqual(d.zz, 3)
        i(d, 'xx', '456')
        self.assertEqual(d.xx, 456)


class TypeSetterLikeTests(TestCase):

    """Tests of the functions that behave similarly to TypeSetter products"""

    def test_list_adder(self):
        """list_adder should add items to list, creating if necessary"""
        d = dummy()
        list_adder(d, 'x', 3)
        self.assertEqual(d.x, [3])
        list_adder(d, 'x', 'abc')
        self.assertEqual(d.x, [3, 'abc'])
        list_adder(d, 'y', [2, 3])
        self.assertEqual(d.x, [3, 'abc'])
        self.assertEqual(d.y, [[2, 3]])

    def test_dict_adder(self):
        """dict_adder should add items to dict, creating if necessary"""
        d = dummy()
        dict_adder(d, 'x', 3)
        self.assertEqual(d.x, {3: None})
        dict_adder(d, 'x', 'ab')
        self.assertEqual(d.x, {3: None, 'a': 'b'})
        dict_adder(d, 'x', ['a', 0])
        self.assertEqual(d.x, {3: None, 'a': 0})
        dict_adder(d, 'y', None)
        self.assertEqual(d.x, {3: None, 'a': 0})
        self.assertEqual(d.y, {None: None})


class LineOrientedConstructorTests(TestCase):

    """Tests of the LineOrientedConstructor class"""

    def test_init_empty(self):
        """LOC empty init should succeed with expected defaults"""
        l = LineOrientedConstructor()
        self.assertEqual(l.Lines, [])
        self.assertEqual(l.LabelSplitter(' ab  cd  '), ['ab', 'cd'])
        self.assertEqual(l.FieldMap, {})
        self.assertEqual(l.Constructor, MappedRecord)
        self.assertEqual(l.Strict, False)

    def test_empty_LOC(self):
        """LOC empty should fail if strict, fill fields if not strict"""
        data = ["abc   def", "3  n", "\t  abc   \txyz\n\n", "fgh   "]
        l = LineOrientedConstructor()
        result = l()
        self.assertEqual(result, {})
        result = l([])
        self.assertEqual(result, {})
        result = l(['   ', '\n\t   '])
        self.assertEqual(result, {})
        result = l(data)
        self.assertEqual(result, {'abc': 'xyz', '3': 'n', 'fgh': None})

    def test_full_LOC(self):
        """LOC should behave as expected when initialized with rich data"""
        data = ["abc\t def", " 3 \t n", "  abc   \txyz\n\n", "x\t5", "fgh   ",
                "x\t3    "]

        class rec(MappedRecord):
            Required = {'abc': []}
        maps = {'abc': list_adder, 'x': int_setter, 'fgh': bool_setter}
        label_splitter = DelimitedSplitter('\t')
        constructor = rec
        strict = True
        loc_bad = LineOrientedConstructor(data, label_splitter, maps,
                                          constructor, strict)
        self.assertRaises(FieldError, loc_bad)
        strict = False
        loc_good = LineOrientedConstructor(data, label_splitter, maps,
                                           constructor, strict)
        result = loc_good()
        assert isinstance(result, rec)
        self.assertEqual(result,
                         {'abc': ['def', 'xyz'], '3': 'n',
                          'fgh': False, 'x': 3})


class fake_dict(dict):

    """Test that constructors return the correct subclass"""
    pass


class FieldWrapperTests(TestCase):

    """Tests of the FieldWrapper factory function"""

    def test_default(self):
        """Default FieldWrapper should wrap fields and labels"""
        fields = list('abcde')
        f = FieldWrapper(fields)
        self.assertEqual(f(''), {})
        self.assertEqual(f('xy za '), {'a': 'xy', 'b': 'za'})
        self.assertEqual(f('1   2\t\t 3  \n4 5 6'),
                         {'a': '1', 'b': '2', 'c': '3', 'd': '4', 'e': '5'})

    def test_splitter(self):
        """FieldWrapper with splitter should use that splitter"""
        fields = ['label', 'count']
        splitter = DelimitedSplitter(':', -1)
        f = FieldWrapper(fields, splitter)
        self.assertEqual(f(''), {})
        self.assertEqual(f('nknasd:'), {'label': 'nknasd', 'count': ''})
        self.assertEqual(
            f('n:k:n:a:sd  '),
            {'label': 'n:k:n:a',
             'count': 'sd'})

    def test_constructor(self):
        """FieldWrapper with constructor should use that constructor"""
        fields = list('abc')
        f = FieldWrapper(fields, constructor=fake_dict)
        self.assertEqual(f('x y'), {'a': 'x', 'b': 'y'})
        assert isinstance(f('x y'), fake_dict)


class StrictFieldWrapperTests(TestCase):

    """Tests of the StrictFieldWrapper factory function"""

    def test_default(self):
        """Default StrictFieldWrapper should wrap fields if count correct"""
        fields = list('abcde')
        f = StrictFieldWrapper(fields)
        self.assertEqual(f('1   2\t\t 3  \n4 5 '),
                         {'a': '1', 'b': '2', 'c': '3', 'd': '4', 'e': '5'})
        self.assertRaises(FieldError, f, '')
        self.assertRaises(FieldError, f, 'xy za ')

    def test_splitter(self):
        """StrictFieldWrapper with splitter should use that splitter"""
        fields = ['label', 'count']
        splitter = DelimitedSplitter(':', -1)
        f = StrictFieldWrapper(fields, splitter)
        self.assertEqual(
            f('n:k:n:a:sd  '),
            {'label': 'n:k:n:a',
             'count': 'sd'})
        self.assertEqual(f('nknasd:'), {'label': 'nknasd', 'count': ''})
        self.assertRaises(FieldError, f, '')

    def test_constructor(self):
        """StrictFieldWrapper with constructor should use that constructor"""
        fields = list('ab')
        f = StrictFieldWrapper(fields, constructor=fake_dict)
        self.assertEqual(f('x y'), {'a': 'x', 'b': 'y'})
        assert isinstance(f('x y'), fake_dict)


class FieldMorpherTests(TestCase):

    """Tests of the FieldMorpher class."""

    def test_default(self):
        """FieldMorpher default should use correct constructors"""
        fm = FieldMorpher({'a': int, 'b': str})
        self.assertEqual(fm({'a': '3', 'b': 456}), {'a': 3, 'b': '456'})

    def test_default_error(self):
        """FieldMorpher default should raise FieldError on unknown fields"""
        fm = FieldMorpher({'a': int, 'b': str})
        self.assertRaises(FieldError, fm, {'a': '3', 'b': 456, 'c': '4'})

    def test_altered_default(self):
        """FieldMorpher with default set should apply it"""
        func = lambda x, y: (str(x), float(y) - 0.5)
        fm = FieldMorpher({'3': str, 4: int}, func)
        # check that recognized values aren't tampered with
        self.assertEqual(fm({3: 3, 4: '4'}), {'3': '3', 4: 4})
        # check that unrecognized values get the appropriate conversion
        self.assertEqual(fm({3: 3, 5: '5'}), {'3': '3', '5': 4.5})

if __name__ == '__main__':
    main()
