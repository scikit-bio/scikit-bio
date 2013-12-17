#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from bipy.util.transform import identity

def iterable(item):
    """If item is iterable, returns item. Otherwise, returns [item].

    Useful for guaranteeing a result that can be iterated over.
    """
    try:
        iter(item)
        return item
    except TypeError:
        return [item]

class FunctionWrapper(object):
    """Wraps a function to hide it from a class so that it isn't a method."""
    def __init__(self, Function):
        self.Function = Function
    def __call__(self, *args, **kwargs):
        return self.Function(*args, **kwargs)

class ClassChecker(object):
    """Container for classes: 'if t in x == True' if t is the right class."""
    def __init__(self, *Classes):
        """Returns a new ClassChecker that accepts specified classes."""
        type_type = type(str)
        for c in Classes:
            if type(c) != type_type:
                raise TypeError, \
                "ClassChecker found non-type object '%s' in parameter list." % c
        self.Classes = list(Classes)

    def __contains__(self, item):
        """Returns True if item is a subclass of one of the classes in self."""
        for c in self.Classes:
            if isinstance(item, c):
                return True
        return False

    def __str__(self):
        """Informal string representation: returns list"""
        return str(self.Classes)

class ConstraintError(Exception):
    """Raised when constraint on a container is violated."""
    pass

class ConstrainedContainer(object):
    """Mixin class providing constraint checking to a container.

    Container should have a Constraint property that __contains__ the items
    that will be allowed in the container. Can also have a Mask property that
    contains a function that will be applied to each item (a) on checking the
    item for validity, and (b) on inserting the item in the container. 

    WARNING: Because the Mask is evaluated both when the item is checked and 
    when it is inserted, any side-effects it has are applied _twice_. This 
    means that any Mask that mutates the object or changes global variables
    is unlikely to do what you want!
    """
    _constraint = None
    Mask = FunctionWrapper(identity)

    def _mask_for_new(self):
        """Returns self.Mask only if different from class data."""
        if self.Mask is not self.__class__.Mask:
            return self.Mask
        else:
            return None

    def __init__(self, Constraint=None, Mask=None):
        """Returns new ConstrainedContainer, incorporating constraint.
        
        WARNING: Does not perform validation. It is the subclass's 
        responsibility to perform validation during __init__ or __new__!
        """
        if Constraint is not None:
            self._constraint = Constraint
        if Mask is not None:
            self.Mask = Mask

    def matchesConstraint(self, constraint):
        """Returns True if all items in self are allowed."""
        #First checks if constraints are compatible. If not, or if the current
        #sequence has no constraint, does item by item search.

        #bail out if self or constraint is empty
        if not constraint or not self:
            return True
        #try checking constraints for compatibility
        if self.Constraint:
            try:
                constraint_ok = True
                for c in self.Constraint:
                    if c not in constraint:
                            constraint_ok = False
                            break
                if constraint_ok:
                    return True
            except TypeError:
                pass #e.g. tried to check wrong type item in string alphabet

        #get here if either self.Constraint is empty, or if we found an item
        #in self.Constraint that wasn't in the other constraint. In either case,
        #we need to check self item by item.
        if self:
            try:
                for i in self:
                    if i not in constraint:
                        return False
            except TypeError:   #e.g. tried to check int in string alphabet
                return False
        return True

    def otherIsValid(self, other):
        """Returns True if other has only items allowed in self.Constraint."""
        #First, checks other.Constrant for compatibility.
        #If other.Constraint is incompatible, checks items in other.
        mask = self.Mask
        constraint = self.Constraint
        if not constraint or not other:
            return True     #bail out if empty
        try:
            #if other has a constraint, check whether it's compatible
            other_constraint = other.Constraint
            if other_constraint:
                for c in map(mask, other_constraint):
                    if c not in constraint:
                        raise ConstraintError
                return True
        except (ConstraintError, AttributeError, TypeError):
            pass
        #get here if other doesn't have a constraint or if other's constraint
        #isn't valid on self's constraint.
        try:
            for item in map(mask, other):
                if item not in constraint:
                    return False
        except TypeError:
            return False    #e.g. tried to check int in str alphabet
        return True

    def itemIsValid(self, item):
        """Returns True if single item is in self.Constraint."""
        try:
            if (not self.Constraint) or self.Mask(item) in self.Constraint:
                return True
            else:
                return False
        except (TypeError, ConstraintError):  #wrong type or not allowed
            return False

    def sequenceIsValid(self, sequence):
        """Returns True if all items in sequence are in self.Constraint."""
        is_valid = self.itemIsValid
        for i in map(self.Mask, sequence):
            if not is_valid(i):
                return False
        return True

    def _get_constraint(self):
        """Accessor for constraint."""
        return self._constraint

    def _set_constraint(self, constraint):
        """Mutator for constraint."""
        if self.matchesConstraint(constraint):
            self._constraint = constraint
        else:
            raise ConstraintError, \
            "Sequence '%s' incompatible with constraint '%s'" % (self, constraint)

    Constraint = property(_get_constraint, _set_constraint)


class ConstrainedString(str, ConstrainedContainer):
    """String that is always valid on a specified constraint."""
    def __new__(cls, data, Constraint=None, Mask=None):
        """Constructor class method for validated ConstrainedString."""
        mask = Mask or cls.Mask
        if data == '':
            pass    #map can't handle an empty sequence, sadly...
        elif isinstance(data, str):
            data = ''.join(map(mask, data))
        else:
            try:
                data = str(map(mask, data))
            except (TypeError, ValueError):
                data = str(mask(data))
        new_string = str.__new__(cls, data)
        curr_constraint = Constraint or new_string.Constraint
        if curr_constraint and new_string:
            for c in new_string:
                try:
                    is_valid = c in curr_constraint
                except TypeError:
                    is_valid = False
                if not is_valid:
                    raise ConstraintError, \
                    "Character '%s' not in constraint '%s'" % (c, curr_constraint)
        return new_string

    def __init__(self, data, Constraint=None, Mask=None):
        """Constructor for validated ConstrainedString."""
        ConstrainedContainer.__init__(self, Constraint, Mask)

    def __add__(self, other):
        """Returns copy of self added to copy of other if constraint correct."""
        if not self.otherIsValid(other):
            raise ConstraintError, \
            "Sequence '%s' doesn't meet constraint" % other
        result = self.__class__(str(self) + ''.join(map(self.Mask, other)), \
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __mul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(str.__mul__(self, multiplier),
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __rmul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(str.__rmul__(self, multiplier),
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __getslice__(self, *args, **kwargs):
        """Make sure slice remembers the constraint."""
        result = self.__class__(str.__getslice__(self, *args, **kwargs), \
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __getitem__(self, index):
        """Make sure extended slice remembers the constraint."""
        items = str.__getitem__(self, index)
        if isinstance(index, slice):
            mask = self._mask_for_new()
            result = self.__class__(items, Constraint=self.Constraint)
            if mask:
                result.Mask = mask
            return result
        else:
            return items

class MappedString(ConstrainedString):
    """As for ConstrainedString, but maps __contained__ and __getitem__."""

    def __contains__(self, item):
        """Ensure that contains applies the mask."""
        try:
            return super(MappedString, self).__contains__(self.Mask(item))
        except (TypeError, ValueError):
            return False


class ConstrainedList(ConstrainedContainer, list):
    """List that is always valid on a specified constraint."""

    def __init__(self, data=None, Constraint=None, Mask=None):
        """Constructor for validated ConstrainedList."""
        ConstrainedContainer.__init__(self, Constraint, Mask)
        if data:
            self.extend(data)

    def __add__(self, other):
        """Returns copy of self added to copy of other if constraint correct."""
        result = self.__class__(list(self) + map(self.Mask, other) , \
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __iadd__(self, other):
        """Adds other to self if constraint correct."""
        other = map(self.Mask, other)
        if self.otherIsValid(other):
            return list.__iadd__(self, other)
        else:
            raise ConstraintError, \
            "Sequence '%s' has items not in constraint '%s'" \
                % (other, self.Constraint)

    def __mul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(list(self) * multiplier,
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __rmul__(self, multiplier):
        """Returns copy of self multiplied by multiplier."""
        result = self.__class__(list(self) * multiplier,
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __setitem__(self, index, item):
        """Sets self[index] to item if item in constraint. Handles slices"""
        if isinstance(index, slice):
            if not self.otherIsValid(item):
                raise ConstraintError, \
                "Sequence '%s' contains items not in constraint '%s'." % \
                (item, self.Constraint)
            item = map(self.Mask, item)
        else:
            if not self.itemIsValid(item):
                raise ConstraintError, "Item '%s' not in constraint '%s'" % \
                    (item, self.Constraint)
            item = self.Mask(item)
        list.__setitem__(self, index, item)

    def append(self, item):
        """Appends item to self."""
        if not self.itemIsValid(item):
            raise ConstraintError, "Item '%s' not in constraint '%s'" % \
                (item, self.Constraint)
        list.append(self, self.Mask(item))

    def extend(self, sequence):
        """Appends sequence to self."""
        if self.otherIsValid(sequence):
            list.extend(self, map(self.Mask, sequence))
        else:
            raise ConstraintError, "Some items in '%s' not in constraint '%s'"\
                % (sequence, self.Constraint)

    def insert(self, position, item):
        """Inserts item at position in self."""
        if not self.itemIsValid(item):
            raise ConstraintError, "Item '%s' not in constraint '%s'" % \
                (item, self.Constraint)
        list.insert(self, position, self.Mask(item))

    def __getslice__(self, *args, **kwargs):
        """Make sure slice remembers the constraint."""
        result = self.__class__(list.__getslice__(self, *args, **kwargs), \
            Constraint=self.Constraint)
        mask = self._mask_for_new()
        if mask:
            result.Mask = mask
        return result

    def __setslice__(self, start, end, sequence):
        """Make sure invalid data can't get into slice."""
        if self.otherIsValid(sequence):
            list.__setslice__(self, start, end, map(self.Mask, sequence))
        else:
            raise ConstraintError, \
                "Sequence '%s' has items not in constraint '%s'"\
                % (sequence, self.Constraint)

class MappedList(ConstrainedList):
    """As for ConstrainedList, but maps items on contains and getitem."""

    def __contains__(self, item):
        """Ensure that contains applies the mask."""
        try:
            return super(MappedList, self).__contains__(self.Mask(item))
        except (TypeError, ValueError):
            return False

class ConstrainedDict(ConstrainedContainer, dict):
    """Dict containing only keys that are valid on a specified constraint.
    
    Default behavior when fed a sequence is to store counts of the items in
    that sequence, which is not the standard dict interface (should raise a
    ValueError instead) but which is surprisingly useful in practice.
    """
    ValueMask = FunctionWrapper(identity)

    def _get_mask_and_valmask(self):
        """Helper method to check whether Mask and ValueMask were set."""
        if self.Mask is self.__class__.Mask:
            mask = None
        else:
            mask = self.Mask

        if self.ValueMask is self.__class__.ValueMask:
            valmask = None
        else:
            valmask = self.ValueMask
        return mask, valmask

    def __init__(self, data=None, Constraint=None, Mask=None, ValueMask=None):
        """Constructor for validated ConstrainedDict."""
        ConstrainedContainer.__init__(self, Constraint, Mask)
        if ValueMask is not None:
            self.ValueMask = ValueMask
        if data:
            try:
                self.update(data)
            except (ValueError, TypeError):
                for d in map(self.Mask, iterable(data)):
                    curr = self.get(d, 0)
                    self[d] = curr + 1

    def __setitem__(self, key, value):
        """Sets self[key] to value if value in constraint."""
        if not self.itemIsValid(key):
            raise ConstraintError, "Item '%s' not in constraint '%s'" % \
                (key, self.Constraint)
        key, value = self.Mask(key), self.ValueMask(value)
        dict.__setitem__(self, key, value)

    def copy(self):
        """Should return copy of self, including constraint."""
        mask, valmask = self._get_mask_and_valmask()
        return self.__class__(self, Constraint=self.Constraint, Mask=mask,
                ValueMask=valmask)

    def fromkeys(self, keys, value=None):
        """Returns new dictionary with same constraint as self."""
        mask, valmask = self._get_mask_and_valmask()
        return self.__class__(dict.fromkeys(keys, value),
            Constraint=self.Constraint, Mask=mask, ValueMask=valmask)

    def setdefault(self, key, default=None):
        """Returns self[key], setting self[key]=default if absent."""
        key, default = self.Mask(key), self.ValueMask(default)
        if key not in self:
            self[key] = default
        return self[key]

    def update(self, other):
        """Updates self with items in other.
        
        Implementation note: currently uses __setitem__, so no need to apply
        masks in this method.
        """
        if not hasattr(other, 'keys'):
            other = dict(other)
        for key in other:
            self[key] = other[key]

class MappedDict(ConstrainedDict):
    """As for ConstrainedDict, but maps keys on contains and getitem."""

    def __contains__(self, item):
        """Ensure that contains applies the mask."""
        try:
            return super(MappedDict, self).__contains__(self.Mask(item))
        except (TypeError, ValueError):
            return False

    def __getitem__(self, item):
        """Ensure that getitem applies the mask."""
        return super(MappedDict, self).__getitem__(self.Mask(item))

    def get(self, item, default=None):
        """Ensure that get applies the mask."""
        return super(MappedDict, self).get(self.Mask(item), default)

    def has_key(self, item):
        """Ensure that has_key applies the mask."""
        return super(MappedDict, self).has_key(self.Mask(item))
