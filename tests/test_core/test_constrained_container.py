#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from copy import copy, deepcopy

from bipy.util.unit_test import TestCase, main
from bipy.core.constrained_container import (FunctionWrapper,
    ConstraintError, ConstrainedContainer, ClassChecker,
    ConstrainedString, ConstrainedList, ConstrainedDict,
    MappedString, MappedList, MappedDict, iterable)

class TopLevelTests(TestCase):
    """Tests of top-level functions """
    
    def test_iterable(self):
        """iterable(x) should return x or [x], always an iterable result"""
        self.assertEqual(iterable('x'), 'x')
        self.assertEqual(iterable(''), '')
        self.assertEqual(iterable(3), [3])
        self.assertEqual(iterable(None), [None])
        self.assertEqual(iterable({'a':1}), {'a':1})
        self.assertEqual(iterable(['a','b','c']), ['a', 'b', 'c'])

class _my_dict(dict):
    """Used for testing subclass behavior of ClassChecker"""
    pass

class ClassCheckerTests(TestCase):
    """Unit tests for the ClassChecker class."""

    def setUp(self):
        """define a few standard checkers"""
        self.strcheck = ClassChecker(str)
        self.intcheck = ClassChecker(int, long)
        self.numcheck = ClassChecker(float, int, long)
        self.emptycheck = ClassChecker()
        self.dictcheck = ClassChecker(dict)
        self.mydictcheck = ClassChecker(_my_dict)
        
    def test_init_good(self):
        """ClassChecker should init OK when initialized with classes"""
        self.assertEqual(self.strcheck.Classes, [str])
        self.assertEqual(self.numcheck.Classes, [float, int, long])
        self.assertEqual(self.emptycheck.Classes, [])

    def test_init_bad(self):
        """ClassChecker should raise TypeError if initialized with non-class"""
        self.assertRaises(TypeError, ClassChecker, 'x')
        self.assertRaises(TypeError, ClassChecker, str, None)

    def test_contains(self):
        """ClassChecker should return True only if given instance of class"""
        self.assertEqual(self.strcheck.__contains__('3'), True)
        self.assertEqual(self.strcheck.__contains__('ahsdahisad'), True)
        self.assertEqual(self.strcheck.__contains__(3), False)
        self.assertEqual(self.strcheck.__contains__({3:'c'}), False)

        self.assertEqual(self.intcheck.__contains__('ahsdahisad'), False)
        self.assertEqual(self.intcheck.__contains__('3'), False)
        self.assertEqual(self.intcheck.__contains__(3.0), False)
        self.assertEqual(self.intcheck.__contains__(3), True)
        self.assertEqual(self.intcheck.__contains__(4**60), True)
        self.assertEqual(self.intcheck.__contains__(4**60 * -1), True)

        d = _my_dict()
        self.assertEqual(self.dictcheck.__contains__(d), True)
        self.assertEqual(self.dictcheck.__contains__({'d':1}), True)
        self.assertEqual(self.mydictcheck.__contains__(d), True)
        self.assertEqual(self.mydictcheck.__contains__({'d':1}), False)

        self.assertEqual(self.emptycheck.__contains__('d'), False)

        self.assertEqual(self.numcheck.__contains__(3), True)
        self.assertEqual(self.numcheck.__contains__(3.0), True)
        self.assertEqual(self.numcheck.__contains__(-3), True)
        self.assertEqual(self.numcheck.__contains__(-3.0), True)
        self.assertEqual(self.numcheck.__contains__(3e-300), True)
        self.assertEqual(self.numcheck.__contains__(0), True)
        self.assertEqual(self.numcheck.__contains__(4**1000), True)
        self.assertEqual(self.numcheck.__contains__('4**1000'), False)

    def test_str(self):
        """ClassChecker str should be the same as str(self.Classes)"""
        for c in [self.strcheck, self.intcheck, self.numcheck, self.emptycheck,
            self.dictcheck, self.mydictcheck]:
            self.assertEqual(str(c), str(c.Classes))

    def test_copy(self):
        """copy.copy should work correctly on ClassChecker"""
        c = copy(self.strcheck)
        assert c is not self.strcheck
        assert '3' in c
        assert 3 not in c
        assert c.Classes is self.strcheck.Classes

    def test_deepcopy(self):
        """copy.deepcopy should work correctly on ClassChecker"""
        c = deepcopy(self.strcheck)
        assert c is not self.strcheck
        assert '3' in c
        assert 3 not in c
        assert c.Classes is not self.strcheck.Classes

class _simple_container(object):
    """example of a container to constrain"""
    def __init__(self, data):
        self._data = list(data)
    def __getitem__(self, item):
        return self._data.__getitem__(item)

class _constrained_simple_container(_simple_container, ConstrainedContainer):
    """constrained version of _simple_container"""
    def __init__(self, data):
        _simple_container.__init__(self, data)
        ConstrainedContainer.__init__(self, None)

class ConstrainedContainerTests(TestCase):
    """Tests of the generic ConstrainedContainer interface."""
    def setUp(self):
        """Make a couple of standard containers"""
        self.alphabet = _constrained_simple_container('abc')
        self.numbers = _constrained_simple_container([1,2,3])
        self.alphacontainer = 'abcdef'
        self.numbercontainer = ClassChecker(int)
        
    def test_matchesConstraint(self):
        """ConstrainedContainer matchesConstraint should return true if items ok"""
        self.assertEqual(self.alphabet.matchesConstraint(self.alphacontainer), \
            True)
        self.assertEqual(self.alphabet.matchesConstraint(self.numbercontainer),\
            False)
        self.assertEqual(self.numbers.matchesConstraint(self.alphacontainer), \
            False)
        self.assertEqual(self.numbers.matchesConstraint(self.numbercontainer),\
            True)

    def test_otherIsValid(self):
        """ConstrainedContainer should use constraint for checking other"""
        self.assertEqual(self.alphabet.otherIsValid('12d8jc'), True)
        self.alphabet.Constraint = self.alphacontainer
        self.assertEqual(self.alphabet.otherIsValid('12d8jc'), False)
        self.alphabet.Constraint = list('abcdefghijkl12345678')
        self.assertEqual(self.alphabet.otherIsValid('12d8jc'), True)
        self.assertEqual(self.alphabet.otherIsValid('z'), False)

    def test_itemIsValid(self):
        """ConstrainedContainer should use constraint for checking item"""
        self.assertEqual(self.alphabet.itemIsValid(3), True)
        self.alphabet.Constraint = self.alphacontainer
        self.assertEqual(self.alphabet.itemIsValid(3), False)
        self.assertEqual(self.alphabet.itemIsValid('a'), True)

    def test_sequenceIsValid(self):
        """ConstrainedContainer should use constraint for checking sequence"""
        self.assertEqual(self.alphabet.sequenceIsValid('12d8jc'), True)
        self.alphabet.Constraint = self.alphacontainer
        self.assertEqual(self.alphabet.sequenceIsValid('12d8jc'), False)
        self.alphabet.Constraint = list('abcdefghijkl12345678')
        self.assertEqual(self.alphabet.sequenceIsValid('12d8jc'), True)
        self.assertEqual(self.alphabet.sequenceIsValid('z'), False)

    def test_Constraint(self):
        """ConstrainedContainer should only allow valid constraints to be set"""
        try:
            self.alphabet.Constraint = self.numbers
        except ConstraintError:
            pass
        else:
            raise AssertionError, \
            "Failed to raise ConstraintError with invalid constraint."
        self.alphabet.Constraint = 'abcdefghi'
        self.alphabet.Constraint = ['a','b', 'c', 1, 2, 3]
        self.numbers.Constraint = range(20)
        self.numbers.Constraint = xrange(20)
        self.numbers.Constraint = [5,1,3,7,2]
        self.numbers.Constraint = {1:'a',2:'b',3:'c'}
        self.assertRaises(ConstraintError, setattr, self.numbers, \
            'Constraint', '1')
            
class ConstrainedStringTests(TestCase):
    """Tests that ConstrainedString can only contain allowed items."""
    
    def test_init_good_data(self):
        """ConstrainedString should init OK if string matches constraint"""
        self.assertEqual(ConstrainedString('abc', 'abcd'), 'abc')
        self.assertEqual(ConstrainedString('', 'abcd'), '')
        items = [1,2,3.2234, (['a'], ['b'],), 'xyz']
        #should accept anything str() does if no constraint is passed
        self.assertEqual(ConstrainedString(items), str(items))
        self.assertEqual(ConstrainedString(items, None), str(items))
        self.assertEqual(ConstrainedString('12345'), str(12345))
        self.assertEqual(ConstrainedString(12345, '1234567890'), str(12345))
        #check that list is formatted correctly and chars are all there
        test_list = [1,2,3,4,5]
        self.assertEqual(ConstrainedString(test_list, '][, 12345'), str(test_list))

    def test_init_bad_data(self):
        """ConstrainedString should fail init if unknown chars in string"""
        self.assertRaises(ConstraintError, ConstrainedString, 1234, '123')
        self.assertRaises(ConstraintError, ConstrainedString, '1234', '123')
        self.assertRaises(ConstraintError, ConstrainedString, [1,2,3], '123')

    def test_add_prevents_bad_data(self):
        """ConstrainedString should allow addition only of compliant string"""
        a = ConstrainedString('123', '12345')
        b = ConstrainedString('444', '4')
        c = ConstrainedString('45', '12345')
        d = ConstrainedString('x')
        self.assertEqual(a + b, '123444')
        self.assertEqual(a + c, '12345')
        self.assertRaises(ConstraintError, b.__add__, c)
        self.assertRaises(ConstraintError, c.__add__, d)
        #should be OK if constraint removed
        b.Constraint = None
        self.assertEqual(b + c, '44445')
        self.assertEqual(b + d, '444x')
        #should fail if we add the constraint back
        b.Constraint = '4x'
        self.assertEqual(b + d, '444x')
        self.assertRaises(ConstraintError, b.__add__, c)
        #check that added strings retain constraint
        self.assertRaises(ConstraintError, (a+b).__add__, d)
   
    def test_mul(self):
        """ConstrainedString mul amd rmul should retain constraint"""
        a = ConstrainedString('123', '12345')
        b = 3*a
        c = b*2
        self.assertEqual(b, '123123123')
        self.assertEqual(c, '123123123123123123')
        self.assertRaises(ConstraintError, b.__add__, 'x')
        self.assertRaises(ConstraintError, c.__add__, 'x')
   
    def test_getslice(self):
        """ConstrainedString getslice should remember constraint"""
        a = ConstrainedString('123333', '12345')
        b = a[2:4]
        self.assertEqual(b, '33')
        self.assertEqual(b.Constraint, '12345')

    def test_getitem(self):
        """ConstrainedString getitem should handle slice objects"""
        a = ConstrainedString('7890543', '1234567890')
        self.assertEqual(a[0], '7')
        self.assertEqual(a[1], '8')
        self.assertEqual(a[-1], '3')
        self.assertRaises(AttributeError, getattr, a[1], 'Alphabet')
        self.assertEqual(a[1:6:2], '804')
        self.assertEqual(a[1:6:2].Constraint, '1234567890')

    def test_init_masks(self):
        """ConstrainedString should init OK with masks"""
        def mask(x):
            return str(int(x) + 3)
        a = ConstrainedString('12333', '45678', mask)
        self.assertEqual(a, '45666')
        assert 'x' not in a
        self.assertRaises(TypeError, a.__contains__, 1)

class MappedStringTests(TestCase):
    """MappedString should behave like ConstrainedString, but should map items."""
    def test_init_masks(self):
        """MappedString should init OK with masks"""
        def mask(x):
            return str(int(x) + 3)
        a = MappedString('12333', '45678', mask)
        self.assertEqual(a, '45666')
        assert 1 in a
        assert 'x' not in a

    

class ConstrainedListTests(TestCase):
    """Tests that bad data can't sneak into ConstrainedLists."""
    def test_init_good_data(self):
        """ConstrainedList should init OK if list matches constraint"""
        self.assertEqual(ConstrainedList('abc', 'abcd'), list('abc'))
        self.assertEqual(ConstrainedList('', 'abcd'), list(''))
        items = [1,2,3.2234, (['a'], ['b'],), list('xyz')]
        #should accept anything str() does if no constraint is passed
        self.assertEqual(ConstrainedList(items), items)
        self.assertEqual(ConstrainedList(items, None), items)
        self.assertEqual(ConstrainedList('12345'), list('12345'))
        #check that list is formatted correctly and chars are all there
        test_list = list('12345')
        self.assertEqual(ConstrainedList(test_list, '12345'), test_list)

    def test_init_bad_data(self):
        """ConstrainedList should fail init with items not in constraint"""
        self.assertRaises(ConstraintError, ConstrainedList, '1234', '123')
        self.assertRaises(ConstraintError,ConstrainedList,[1,2,3],['1','2','3'])

    def test_add_prevents_bad_data(self):
        """ConstrainedList should allow addition only of compliant data"""
        a = ConstrainedList('123', '12345')
        b = ConstrainedList('444', '4')
        c = ConstrainedList('45', '12345')
        d = ConstrainedList('x')
        self.assertEqual(a + b, list('123444'))
        self.assertEqual(a + c, list('12345'))
        self.assertRaises(ConstraintError, b.__add__, c)
        self.assertRaises(ConstraintError, c.__add__, d)
        #should be OK if constraint removed
        b.Constraint = None
        self.assertEqual(b + c, list('44445'))
        self.assertEqual(b + d, list('444x'))
        #should fail if we add the constraint back
        b.Constraint = {'4':1, 5:2}
        self.assertRaises(ConstraintError, b.__add__, c)
                
    
    def test_iadd_prevents_bad_data(self):
        """ConstrainedList should allow in-place addition only of compliant data"""
        a = ConstrainedList('12', '123')
        a += '2'
        self.assertEqual(a, list('122'))
        self.assertEqual(a.Constraint, '123')
        self.assertRaises(ConstraintError, a.__iadd__, '4')
    
    def test_imul(self):
        """ConstrainedList imul should preserve constraint"""
        a = ConstrainedList('12', '123')
        a *= 3
        self.assertEqual(a, list('121212'))
        self.assertEqual(a.Constraint, '123')

    def test_mul(self):
        """ConstrainedList mul should preserve constraint"""
        a = ConstrainedList('12', '123')
        b = a * 3
        self.assertEqual(b, list('121212'))
        self.assertEqual(b.Constraint, '123')

    def test_rmul(self):
        """ConstrainedList rmul should preserve constraint"""
        a = ConstrainedList('12', '123')
        b = 3 * a
        self.assertEqual(b, list('121212'))
        self.assertEqual(b.Constraint, '123')

    def test_setitem(self):
        """ConstrainedList setitem should work only if item in constraint"""
        a = ConstrainedList('12', '123')
        a[0] = '3'
        self.assertEqual(a, list('32'))
        self.assertRaises(ConstraintError, a.__setitem__, 0, 3)
        a = ConstrainedList('1'*20, '123')
        self.assertRaises(ConstraintError, a.__setitem__, slice(0,1,1), [3])
        self.assertRaises(ConstraintError, a.__setitem__, slice(0,1,1), ['111'])
        a[2:9:2] = '3333'
        self.assertEqual(a, list('11313131311111111111'))

    def test_append(self):
        """ConstrainedList append should work only if item in constraint"""
        a = ConstrainedList('12', '123')
        a.append('3')
        self.assertEqual(a, list('123'))
        self.assertRaises(ConstraintError, a.append, 3)

    def test_extend(self):
        """ConstrainedList extend should work only if all items in constraint"""
        a = ConstrainedList('12', '123')
        a.extend('321')
        self.assertEqual(a, list('12321'))
        self.assertRaises(ConstraintError, a.extend, ['1','2', 3])

    def test_insert(self):
        """ConstrainedList insert should work only if item in constraint"""
        a = ConstrainedList('12', '123')
        a.insert(0, '2')
        self.assertEqual(a, list('212'))
        self.assertRaises(ConstraintError, a.insert, 0, [2])

    def test_getslice(self):
        """ConstrainedList getslice should remember constraint"""
        a = ConstrainedList('123333', '12345')
        b = a[2:4]
        self.assertEqual(b, list('33'))
        self.assertEqual(b.Constraint, '12345')

    def test_setslice(self):
        """ConstrainedList setslice should fail if slice has invalid chars"""
        a = ConstrainedList('123333', '12345')
        a[2:4] = ['2','2']
        self.assertEqual(a, list('122233'))
        self.assertRaises(ConstraintError, a.__setslice__, 2,4, [2,2])
        a[:] = []
        self.assertEqual(a, [])
        self.assertEqual(a.Constraint, '12345')

    def test_setitem_masks(self):
        """ConstrainedList setitem with masks should transform input"""
        a = ConstrainedList('12333', range(5), lambda x: int(x) + 1)
        self.assertEqual(a, [2,3,4,4,4])
        self.assertRaises(ConstraintError, a.append, 4)
        b = a[1:3]
        assert b.Mask is a.Mask
        assert '1' not in a
        assert '2' not in a
        assert 2 in a
        assert 'x' not in a

class MappedListTests(TestCase):
    """MappedList should behave like ConstrainedList, but map items."""
    def test_setitem_masks(self):
        """MappedList setitem with masks should transform input"""
        a = MappedList('12333', range(5), lambda x: int(x) + 1)
        self.assertEqual(a, [2,3,4,4,4])
        self.assertRaises(ConstraintError, a.append, 4)
        b = a[1:3]
        assert b.Mask is a.Mask
        assert '1' in a
        assert 'x' not in a

class ConstrainedDictTests(TestCase):
    """Tests that bad data can't sneak into ConstrainedDicts."""
    def test_init_good_data(self):
        """ConstrainedDict should init OK if list matches constraint"""
        self.assertEqual(ConstrainedDict(dict.fromkeys('abc'), 'abcd'), \
            dict.fromkeys('abc'))
        self.assertEqual(ConstrainedDict('', 'abcd'), dict(''))
        items = [1,2,3.2234, tuple('xyz')]
        #should accept anything dict() does if no constraint is passed
        self.assertEqual(ConstrainedDict(dict.fromkeys(items)), \
            dict.fromkeys(items))
        self.assertEqual(ConstrainedDict(dict.fromkeys(items), None), \
            dict.fromkeys(items))
        self.assertEqual(ConstrainedDict([(x,1) for x in '12345']), \
            dict.fromkeys('12345', 1))
        #check that list is formatted correctly and chars are all there
        test_dict = dict.fromkeys('12345')
        self.assertEqual(ConstrainedDict(test_dict, '12345'), test_dict)

    def test_init_sequence(self):
        """ConstrainedDict should init from sequence, unlike normal dict"""
        self.assertEqual(ConstrainedDict('abcda'), {'a':2,'b':1,'c':1,'d':1})

    def test_init_bad_data(self):
        """ConstrainedDict should fail init with items not in constraint"""
        self.assertRaises(ConstraintError, ConstrainedDict, \
            dict.fromkeys('1234'), '123')
        self.assertRaises(ConstraintError,ConstrainedDict, \
        dict.fromkeys([1,2,3]),['1','2','3'])
   
    def test_setitem(self):
        """ConstrainedDict setitem should work only if key in constraint"""
        a = ConstrainedDict(dict.fromkeys('12'), '123')
        a['1'] = '3'
        self.assertEqual(a, {'1':'3','2':None})
        self.assertRaises(ConstraintError, a.__setitem__, 1, '3')

    def test_copy(self):
        """ConstrainedDict copy should retain constraint"""
        a = ConstrainedDict(dict.fromkeys('12'), '123')
        b = a.copy()
        self.assertEqual(a.Constraint, b.Constraint)
        self.assertRaises(ConstraintError, a.__setitem__, 1, '3')
        self.assertRaises(ConstraintError, b.__setitem__, 1, '3')

    def test_fromkeys(self):
        """ConstrainedDict instance fromkeys should retain constraint"""
        a = ConstrainedDict(dict.fromkeys('12'), '123')
        b = a.fromkeys('23')
        self.assertEqual(a.Constraint, b.Constraint)
        self.assertRaises(ConstraintError, a.__setitem__, 1, '3')
        self.assertRaises(ConstraintError, b.__setitem__, 1, '3')
        b['2'] = 5
        self.assertEqual(b, {'2':5, '3':None})

    def test_setdefault(self):
        """ConstrainedDict setdefault shouldn't allow bad keys"""
        a = ConstrainedDict({'1':None, '2': 'xyz'}, '123')
        self.assertEqual(a.setdefault('2', None), 'xyz')
        self.assertEqual(a.setdefault('1', None), None)
        self.assertRaises(ConstraintError, a.setdefault, 'x', 3)
        a.setdefault('3', 12345)
        self.assertEqual(a, {'1':None, '2':'xyz', '3': 12345})

    def test_update(self):
        """ConstrainedDict should allow update only of compliant data"""
        a = ConstrainedDict(dict.fromkeys('123'), '12345')
        b = ConstrainedDict(dict.fromkeys('444'), '4')
        c = ConstrainedDict(dict.fromkeys('45'), '12345')
        d = ConstrainedDict([['x','y']])
        a.update(b)
        self.assertEqual(a, dict.fromkeys('1234'))
        a.update(c)
        self.assertEqual(a, dict.fromkeys('12345'))
        self.assertRaises(ConstraintError, b.update, c)
        self.assertRaises(ConstraintError, c.update, d)
        #should be OK if constraint removed
        b.Constraint = None
        b.update(c)
        self.assertEqual(b, dict.fromkeys('45'))
        b.update(d)
        self.assertEqual(b, {'4':None, '5':None, 'x':'y'})
        #should fail if we add the constraint back
        b.Constraint = {'4':1, 5:2, '5':1, 'x':1}
        self.assertRaises(ConstraintError, b.update, {4:1})
        b.update({5:1})
        self.assertEqual(b, {'4':None, '5':None, 'x':'y', 5:1})
   
    def test_setitem_masks(self):
        """ConstrainedDict setitem should work only if key in constraint"""
        key_mask = str
        val_mask = lambda x: int(x) + 3
        d = ConstrainedDict({1:4, 2:6}, '123', key_mask, val_mask)
        d[1] = '456'
        self.assertEqual(d, {'1':459,'2':9,})
        d['1'] = 234
        self.assertEqual(d, {'1':237,'2':9,})
        self.assertRaises(ConstraintError, d.__setitem__, 4, '3')
        e = d.copy()
        assert e.Mask is d.Mask
        assert '1' in d
        assert not 1 in d

class MappedDictTests(TestCase):
    """MappedDict should work like ConstrainedDict, but map keys."""
    
    def test_setitem_masks(self):
        """MappedDict setitem should work only if key in constraint"""
        key_mask = str
        val_mask = lambda x: int(x) + 3
        d = MappedDict({1:4, 2:6}, '123', key_mask, val_mask)
        d[1] = '456'
        self.assertEqual(d, {'1':459,'2':9,})
        d['1'] = 234
        self.assertEqual(d, {'1':237,'2':9,})
        self.assertRaises(ConstraintError, d.__setitem__, 4, '3')
        e = d.copy()
        assert e.Mask is d.Mask
        assert '1' in d
        assert 1 in d
        assert 1 not in d.keys()
        assert 'x' not in d.keys()

    def test_getitem(self):
        """MappedDict getitem should automatically map key."""
        key_mask = str
        d = MappedDict({}, '123', key_mask)
        self.assertEqual(d, {})
        d['1'] = 5
        self.assertEqual(d, {'1':5})
        self.assertEqual(d[1], 5)

    def test_get(self):
        """MappedDict get should automatically map key."""
        key_mask = str
        d = MappedDict({}, '123', key_mask)
        self.assertEqual(d, {})
        d['1'] = 5
        self.assertEqual(d, {'1':5})
        self.assertEqual(d.get(1, 'x'), 5)
        self.assertEqual(d.get(5, 'x'), 'x')

    def test_has_key(self):
        """MappedDict has_key should automatically map key."""
        key_mask = str
        d = MappedDict({}, '123', key_mask)
        self.assertEqual(d, {})
        d['1'] = 5
        assert d.has_key('1')
        assert d.has_key(1)
        assert not d.has_key('5')


if __name__ == "__main__":
    main()