#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from future.builtins import int
from future.utils import viewitems

from copy import deepcopy
from skbio.core.exception import FieldError


def string_and_strip(*items):
    """Converts items to strings and strips them."""
    return [str(i).strip() for i in items]


def DelimitedSplitter(delimiter=None, max_splits=1):
    """Returns function that returns stripped fields split by delimiter.

    Unlike the default behavior of split, max_splits can be negative, in
    which case it counts from the end instead of the start (i.e. splits
    at the _last_ delimiter, last two delimiters, etc. for -1, -2, etc.)
    However, if the delimiter is None (the default) and max_splits is
    negative, will not preserve internal spaces.

    Note: leaves empty fields in place.
    """
    is_int = isinstance(max_splits, int)
    if is_int and (max_splits > 0):
        def parser(line):
            return [i.strip() for i in line.split(delimiter, max_splits)]
    elif is_int and (max_splits < 0):
        def parser(line):
            to_insert = delimiter or ' '  # re-join fields w/ space if None
            fields = line.split(delimiter)
            if (fields == []) or (fields == ['']):
                return []  # empty string or only delimiter: return nothing
            # if not enough fields, count from the start, not the end
            if len(fields) < max_splits:
                first_fields = fields[0]
                last_fields = fields[1:]
            # otherwise, count off the last n fields and join the remainder
            else:
                first_fields = fields[:max_splits]
                last_fields = fields[max_splits:]
            pieces = []
            # if first_fields is empty, don't make up an extra empty string
            if first_fields:
                pieces.append(to_insert.join(first_fields))
            pieces.extend(last_fields)
            return [i.strip() for i in pieces]

    else:  # ignore max_splits if it was 0
        def parser(line):
            return [i.strip() for i in line.split(delimiter)]
    return parser

# The following provide examples of the kinds of functions DelimitedSplitter
# returns.
semi_splitter = DelimitedSplitter(';', None)
space_pairs = DelimitedSplitter(None)
equal_pairs = DelimitedSplitter('=')
last_colon = DelimitedSplitter(':', -1)


class GenericRecord(dict):

    """Holds data for a generic field ->: value mapping.

    Override Required with {name:prototype} mapping. Each required name
    will get a deepcopy of its prototype. For example, use an empty list to
    guarantee that each instance has its own list for a particular field to
    which items can be appended.

    Raises AttributeError on attempt to delete required item, but does not
    raise an exception on attempt to delete absent item.

    This class explicitly does _not_ override __getitem__ or __setitem__ for
    performance reasons: if you need to transform keys on get/set or if you
    need to access items as attributes and vice versa, use MappedRecord
    instead.
    """
    Required = {}

    def __init__(self, *args, **kwargs):
        """Reads kwargs as properties of self."""
        # perform init on temp dict to preserve interface: will then translate
        # aliased keys when loading into self
        temp = {}
        dict.__init__(temp, *args, **kwargs)
        self.update(temp)
        for name, prototype in viewitems(self.Required):
            if name not in self:
                self[name] = deepcopy(prototype)

    def __delitem__(self, item):
        """Deletes item or raises exception if item required.

        Note: Fails silently if item absent.
        """
        if item in self.Required:
            raise AttributeError("%s is a required item" % (item,))
        try:
            super(GenericRecord, self).__delitem__(item)
        except KeyError:
            pass

    def copy(self):
        """Coerces copy to correct type"""
        temp = self.__class__(super(GenericRecord, self).copy())
        # don't forget to copy attributes!
        for attr, val in viewitems(self.__dict__):
            temp.__dict__[attr] = deepcopy(val)
        return temp


class MappedRecord(GenericRecord):

    """GenericRecord that maps names of fields onto standardized names.

    Override Aliases in subclass for new mapping of OldName->NewName. Each
    OldName can have only one NewName, but it's OK if several OldNames map
    to the same NewName.

    Note: can access fields either as items or as attributes. In addition,
    can access either using nonstandard names or using standard names.

    Implementation note: currently, just a dict with appropriate get/set
    overrides and ability to access items as attributes. Attribute access
    is about 10x slower than in GenericRecord, so make sure you need the
    additional capabilities if you use MappedRecord instead of GenericRecord.

    WARNING: MappedRecord pretends to have every attribute, so will never raise
    AttributeError when trying to find an unknown attribute. This feature can
    cause surprising interactions when a Delegator is delegating its
    attributes to a MappedRecord, since any attributes defined in __init__ will
    be set in the MappedRecord and not in the object itself. The solution is
    to use the self.__dict__['AttributeName'] = foo syntax to force the
    attributes to be set in the object and not the MappedRecord to which it
    forwards.
    """
    Aliases = {}

    DefaultValue = None

    def _copy(self, prototype):
        """Returns a copy of item."""
        if hasattr(prototype, 'copy'):
            return prototype.copy()
        elif isinstance(prototype, list):
            return prototype[:]
        elif isinstance(prototype, str) or isinstance(prototype, int) or\
                isinstance(prototype, long) or isinstance(prototype, tuple)\
                or isinstance(prototype, complex) or prototype is None:
            return prototype  # immutable type: use directly
        else:
            return deepcopy(prototype)

    def __init__(self, *args, **kwargs):
        """Reads kwargs as properties of self."""
        # perform init on temp dict to preserve interface: will then translate
        # aliased keys when loading into self
        temp = {}
        unalias = self.unalias
        dict.__init__(temp, *args, **kwargs)
        for key, val in viewitems(temp):
            self[unalias(key)] = val
        for name, prototype in viewitems(self.Required):
            new_name = unalias(name)
            if new_name not in self:
                self[new_name] = self._copy(prototype)

    def unalias(self, key):
        """Returns dealiased name for key, or key if not in alias."""
        try:
            return self.Aliases.get(key, key)
        except TypeError:
            return key

    def __getattr__(self, attr):
        """Returns None if field is absent, rather than raising exception."""
        if attr in self:
            return self[attr]
        elif attr in self.__dict__:
            return self.__dict__[attr]
        elif attr.startswith('__'):  # don't retrieve private class attrs
            raise AttributeError
        elif hasattr(self.__class__, attr):
            return getattr(self.__class__, attr)
        else:
            return self._copy(self.DefaultValue)

    def __setattr__(self, attr, value):
        """Sets attribute in self if absent, converting name if necessary."""
        normal_attr = self.unalias(attr)
        # we overrode __getattr__, so have to simulate getattr(self, attr) by
        # calling superclass method and checking for AttributeError.
        # BEWARE: dict defines __getattribute__, not __getattr__!
        try:
            super(MappedRecord, self).__getattribute__(normal_attr)
            super(MappedRecord, self).__setattr__(normal_attr, value)
        except AttributeError:
            self[normal_attr] = value

    def __delattr__(self, attr):
        """Deletes attribute, converting name if necessary. Fails silently."""
        normal_attr = self.unalias(attr)
        if normal_attr in self.Required:
            raise AttributeError("%s is a required attribute" % (attr,))
        else:
            try:
                super(MappedRecord, self).__delattr__(normal_attr)
            except AttributeError:
                del self[normal_attr]

    def __getitem__(self, item):
        """Returns default if item is absent, rather than raising exception."""
        normal_item = self.unalias(item)
        return self.get(normal_item, self._copy(self.DefaultValue))

    def __setitem__(self, item, val):
        """Sets item, converting name if necessary."""
        super(MappedRecord, self).__setitem__(self.unalias(item), val)

    def __delitem__(self, item):
        """Deletes item, converting name if necessary. Fails silently."""
        normal_item = self.unalias(item)
        super(MappedRecord, self).__delitem__(normal_item)

    def __contains__(self, item):
        """Tests membership, converting name if necessary."""
        return super(MappedRecord, self).__contains__(self.unalias(item))

    def get(self, item, default):
        """Returns self[item] or default if not present. Silent on unhashable.
        """
        try:
            return super(MappedRecord, self).get(self.unalias(item), default)
        except TypeError:
            return default

    def setdefault(self, key, default=None):
        """Returns self[key] or default (and sets self[key]=default)"""
        return super(MappedRecord, self).setdefault(self.unalias(key), default)

    def update(self, *args, **kwargs):
        """Updates self with items in other"""
        temp = {}
        unalias = self.unalias
        temp.update(*args, **kwargs)
        for key, val in viewitems(temp):
            self[unalias(key)] = val

# The following methods are useful for handling particular types of fields in
# line-oriented parsers


def TypeSetter(constructor=None):
    """Returns function that takes obj, field, val and sets obj.field = val.

    constructor can be any callable that returns an object.
    """
    if constructor:
        def setter(obj, field, val):
            setattr(obj, field, constructor(val))
    else:
        def setter(obj, field, val):
            setattr(obj, field, val)
    return setter

int_setter = TypeSetter(int)
str_setter = TypeSetter(str)
list_setter = TypeSetter(list)
tuple_setter = TypeSetter(tuple)
dict_setter = TypeSetter(dict)
float_setter = TypeSetter(float)
complex_setter = TypeSetter(complex)
bool_setter = TypeSetter(bool)
identity_setter = TypeSetter()


def list_adder(obj, field, val):
    """Adds val to list in obj.field, creating list if necessary."""
    try:
        getattr(obj, field).append(val)
    except AttributeError:
        setattr(obj, field, [val])


def dict_adder(obj, field, val):
    """If val is a sequence, adds key/value pair in obj.field: else adds key.
    """
    try:
        key, value = val
    except (ValueError, TypeError):
        key, value = val, None
    try:
        getattr(obj, field)[key] = value
    except AttributeError:
        setattr(obj, field, {key: value})


class LineOrientedConstructor(object):

    """Constructs a MappedRecord from a sequence of lines."""

    def __init__(self, Lines=None, LabelSplitter=space_pairs, FieldMap=None,
                 Constructor=MappedRecord, Strict=False):
        """Returns new LineOrientedConstructor.

        Fields:
            Lines: set of lines to construct record from (for convenience).
            Default is None.

            LabelSplitter: function that returns (label, data) tuple.
            Default is to split on first space and strip components.

            FieldMap: dict of {fieldname:handler} functions. Each function
            has the signature (obj, field, val) and performs an inplace
            action like setting field to val or appending val to field.
            Default is empty dict.

            Constructor: constructor for the resulting object.
            Default is MappedRecord: beware of using constructors that don't
            subclass MappedRecord.

            Strict: boolean controlling whether to raise error on unrecognized
            field. Default is False.
        """
        self.Lines = Lines or []
        self.LabelSplitter = LabelSplitter
        self.FieldMap = FieldMap or {}
        self.Constructor = Constructor
        self.Strict = Strict

    def __call__(self, Lines=None):
        """Returns the record constructed from Lines, or self.Lines"""
        if Lines is None:
            Lines = self.Lines
        result = self.Constructor()
        fieldmap = self.FieldMap
        aka = result.unalias

        splitter = self.LabelSplitter
        for line in Lines:
            # find out how many items we got, setting key and val appropiately
            items = list(splitter(line))
            num_items = len(items)
            if num_items == 2:  # typical case: key-value pair
                raw_field, val = items
            elif num_items > 2:
                raw_field = items[0]
                val = items[1:]
            elif len(items) == 1:
                raw_field, val = items[0], None
            elif not items:  # presumably had line with just a delimiter?
                continue
            # figure out if we know the field under its original name or as
            # an alias
            if raw_field in fieldmap:
                field, mapper = raw_field, fieldmap[raw_field]
            else:
                new_field = aka(raw_field)
                if new_field in fieldmap:
                    field, mapper = new_field, fieldmap[new_field]
                else:
                    if self.Strict:
                        raise FieldError(
                            "Got unrecognized field %s" %
                            (raw_field,))
                    else:
                        identity_setter(result, raw_field, val)
                    continue
            # if we found the field in the fieldmap, apply the correct function
            try:
                mapper(result, field, val)
            except:  # Warning: this is a catchall for _any_ exception,
                        # and may mask what's actually going wrong.
                if self.Strict:
                    raise FieldError("Could not handle line %s" % (line,))
        return result


def FieldWrapper(fields, splitter=None, constructor=None):
    """Returns dict containing field->val mapping, one level.

    fields should be list of fields, in order.
    splitter should be something like a DelimitedSplitter that converts the
        line into a sequence of fields.
    constructor is a callable applied to the dict after construction.

    Call result on a _single_ line, not a list of lines.

    Note that the constructor should take a dict and return an object of some
    useful type. Additionally, it is the _constructor's_ responsibility to
    complain if there are not enough fields, since zip will silently truncate
    at the shorter sequence. This is actually useful in the case where many of
    the later fields are optional.
    """
    if splitter is None:
        splitter = DelimitedSplitter(None, None)
    if constructor:
        def parser(line):
            return constructor(dict(zip(fields, splitter(line))))
    else:
        def parser(line):
            return dict(zip(fields, splitter(line)))
    return parser


def StrictFieldWrapper(fields, splitter=None, constructor=None):
    """Returns dict containing field->val mapping, one level.

    fields should be list of fields, in order.
    splitter should be something like a DelimitedSplitter that converts the
        line into a sequence of fields.
    constructor is a callable applied to the dict after construction.

    Call result on a _single_ line, not a list of lines.

    Note that the constructor should take a dict and return an object of some
    useful type. Raises RecordError if the wrong number of fields are returned
    from the split.
    """
    if splitter is None:
        splitter = DelimitedSplitter(None, None)
    if constructor:
        def parser(line):
            items = splitter(line)
            if len(items) != len(fields):
                raise FieldError("Expected %s items but got %s: %s" %
                                 (len(fields), len(items), items))
            return constructor(dict(zip(fields, items)))
    else:
        def parser(line):
            items = splitter(line)
            if len(items) != len(fields):
                raise FieldError("Expected %s items but got %s: %s" %
                                 (len(fields), len(items), items))
            return dict(zip(fields, items))
    return parser


def raise_unknown_field(field, data):
    """Raises a FieldError, displaying the offending field and data."""
    raise FieldError("Got unknown field %s with data %s" % (field, data))


class FieldMorpher(object):

    """When called, applies appropriate constructors to each value of dict.

    Initialize using a dict of fieldname:constructor pairs.
    """

    def __init__(self, Constructors, Default=raise_unknown_field):
        """Returns a new FieldMorpher, using appropriate constructors.

        If a field is unknown, will try to set key and value to the results
        of Default(key, value): in other words, the signature of Default should
        take a key and a value and should return a key and a value. The
        built-in value of Default raises a FieldError instead, but it will
        often be useful to do things like return the key/value pair unchanged,
        or to strip the key and the value and then add them.
        """
        self.Constructors = Constructors
        self.Default = Default

    def __call__(self, data):
        """Returns a new dict containing information converted from data."""
        result = {}
        default = self.Default
        cons = self.Constructors
        for key, val in viewitems(data):
            if key in cons:
                result[key] = cons[key](val)
            else:
                new_key, new_val = default(key, val)
                # if we now recognize the key, use its constructor on the old
                # val
                if new_key in cons:
                    result[new_key] = cons[new_key](val)
                # otherwise, enter the new key and the new val
                else:
                    result[new_key] = new_val
        return result
