#!/usr/bin/env python
"""Provides transformations of functions and other objects.

Includes:

Standard combinatorial higher-order functions adapted from David Mertz (2003),
"Text Processing in Python", Chapter 1.

Functions for performing complex tests on strings, e.g. includes_any or
includes_all.

Functions for generating combinations, permutations, or cartesian products
of lists.
"""
from __future__ import division
from operator import add, and_, or_
from numpy import logical_and, logical_or, logical_not
from string import maketrans
from cogent.maths.stats.util import Freqs
from cogent.util.misc import identity, select

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Sandra Smit", "Rob Knight","Zongzhi Liu"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Production"

#standard combinatorial HOF's from Mertz
def apply_each(functions, *args, **kwargs):
    """Returns list containing result of applying each function to args."""
    return [f(*args, **kwargs) for f in functions]

def bools(items):
    """Returns list of booleans: reflects state of each item."""
    return map(bool, items)

def bool_each(functions, *args, **kwargs):
    """Returns list of booleans: results of applying each function to args."""
    return bools(apply_each(functions,*args, **kwargs))

def conjoin(functions, *args, **kwargs):
    """Returns True if all functions return True when applied to args."""
    for f in functions:
        if not f(*args, **kwargs):
            return False
    return True

def all(functions):
    """Returns function that returns True when all components return True."""
    def apply_to(*args, **kwargs):
        return conjoin(functions, *args, **kwargs)
    return apply_to

def both(f,g):
    """Returns function that returns True when functions f and g return True."""
    def apply_to(*args, **kwargs):
        #use operator.__and__ to make it compatible to numpy array operation
        #return logical_and(f(*args, **kwargs), g(*args, **kwargs))
        return f(*args, **kwargs) and g(*args, **kwargs)
    return apply_to

def disjoin(functions, *args, **kwargs):
    """Returns True if any of the component functions return True."""
    for f in functions:
        if f(*args, **kwargs):
            return True
    return False

def any(functions):
    """Returns a function that returns True if any component returns True."""
    def apply_to(*args, **kwargs):
        return disjoin(functions, *args, **kwargs)
    return apply_to

def either(f,g):
    """Returns a function that returns True if either f or g returns True."""
    def apply_to(*args, **kwargs):
        return f(*args, **kwargs) or g(*args, **kwargs)
    return apply_to

def negate(functions, *args, **kwargs):
    """Returns True if all functions return False."""
    for f in functions:
        if f(*args, **kwargs):
            return False
    return True

def none(functions):
    """Returns a function that returns True if all components return False."""
    def apply_to(*args, **kwargs):
        return negate(functions, *args, **kwargs)
    return apply_to

def neither(f,g):
    """Returns a function that returns True if neither f not g returns True."""
    def apply_to(*args, **kwargs):
        return not(f(*args, **kwargs)) and not(g(*args, **kwargs))
    return apply_to

def compose(f,g):
    """Returns a function that returns the result of applying f to g(x)."""
    def apply_to(*args, **kwargs):
        return f(g(*args, **kwargs))
    return apply_to

def compose_many(*functions):
    """Returns a function that composes all input functions."""
    funs = list(functions)
    funs.reverse()
    def apply_to(*args, **kwargs):
        result = funs[0](*args, **kwargs)
        for f in funs[1:]:
            next = f(result)
            result = next
        return result
    return apply_to

#factory for making functions that apply to sequences

def per_shortest(total, x, y):
    """Divides total by min(len(x), len(y)).
    
    Useful for normalizing per-item results from sequences that are zipped
    together. Always returns 0 if one of the sequences is empty (to 
    avoid divide by zero error).
    """
    shortest = min(len(x), len(y))
    if not shortest:
        return 0
    else:
        return total/shortest

def per_longest(total, x, y):
    """Divides total by max(len(x), len(y)).
    
    Useful for normalizing per-item results from sequences that are zipped
    together. Always returns 0 if one of the sequences is empty (to 
    avoid divide by zero error).
    """
    longest = max(len(x), len(y))
    if not longest:
        return 0
    else:
        return total/longest

class for_seq(object):
    """Returns function that applies f(i,j) to i,j in zip(first, second).

    f: f(i,j) applying to elements of the sequence.

    aggregator: method to reduce the list of results to a scalar. Default: sum.
    
    normalizer: f(total, i, j) that normalizes the total as a function of
    i and j. Default is length_normalizer (divides by the length of the shorter
    of i and j). If normalizer is None, no normalization is performed.

    Will always truncate to length of the shorter sequence (because of the use
    of zip).
    """
    def __init__(self, f, aggregator=sum, normalizer=per_shortest):
        self.f = f
        self.aggregator = aggregator
        self.normalizer = normalizer

    def __call__(self, first, second):
        f = self.f
        if self.normalizer is None:
            return self.aggregator([f(i,j) for i,j in zip(first, second)])
        else:
            return self.normalizer(self.aggregator(\
                [f(i,j) for i,j in zip(first,second)]), first, second)

#convenience functions for modifying objects

def has_field(field_name):
    """Returns a function that returns True if the obj has the field_name."""
    def field_checker(obj):
        return hasattr(obj, field_name)
    return field_checker

def extract_field(field_name, constructor=None):
    """Returns a function that returns the value of the specified field.

    If set, the constructor function will be applied to obj.field_name.
    Returns None if the constructor fails or the attribute doesn't exist.
    """
    f = constructor or identity
    def result(x):
        try:
            return f(getattr(x, field_name))
        except:
            return None
    return result

def test_field(field_name, constructor=None):
    """Returns True if obj.field_name is True. False otherwise.

    If set, the constructor will be applied to the obj.field_name.
    If accessing the field raises an exception, returns False.
    """
    extractor = extract_field(field_name, constructor)
    def result(x):
        return bool(extractor(x))
    return result

def index(constructor=None, overwrite=False):
    """Returns a function that constructs a dict mapping constructor to object.
   
    If constructor is None, tries to make the objects the keys (only works if
    the objects defines __hash__.

    If overwrite is True, returns single item for each value. Otherwise, always
    returns a list (even if one item) for each value.
    """
    f = constructor or identity
    if overwrite:
        def result(items):
            return dict([(f(i), i) for i in items])
        return result
    else:
        def result(items):
            index = {}
            for i in items:
                curr = f(i)
                if curr in index:
                    index[curr].append(i)
                else:
                    index[curr] = [i]
            return index
        return result

def test_container(container):
    """Returns function that tests safely if item in container.

    Does not raise TypeError (e.g. for dict checks.)
    """
    def result(item):
        try:
            return item in container
        except TypeError:
            return False
    return result

allchars = maketrans('','')

def trans_except(good_chars, default):
    """Returns translation table mapping all but the 'good chars' to default."""
    all_list = list(allchars)
    for i, char in enumerate(all_list):
        if char not in good_chars:
            all_list[i] = default
    return ''.join(all_list)
    
def trans_all(bad_chars, default):
    """Returns translation table mapping all the 'bad chars' to default."""
    all_list = list(allchars)
    for i, char in enumerate(all_list):
        if char in bad_chars:
            all_list[i] = default
    return ''.join(all_list)

def make_trans(frm='', to='', default=None):
    """Like built-in maketrans, but sets all unspecified chars to default."""
    if default is None:
        all_list = list(allchars)
    else:
        if len(default) != 1:
            raise ValueError, 'make_trans default must be single char: got %s' \
            % default
        all_list = [default] * 256
    for orig, new in zip(frm, to):
        all_list[ord(orig)] = new
    return ''.join(all_list)
    

def find_any(words, case_sens = False):
    """Tests if any of the given words occurs in the given string.
    
    This filter is case INsensitive by default.
    """
    if not case_sens:
        used_words = [w.lower() for w in words]
    else:
        used_words = words
    
    def apply_to(s):
        if not case_sens:
            used_s = s.lower()
        else:
            used_s = s
        for w in used_words:
            if used_s.find(w) > -1:
                return True
        return False
    return apply_to

def find_no(words, case_sens = False):
    """Returns True if none of the words appears in s.
    
    This filter is case INsensitive by default.
    """
    f=find_any(words,case_sens)
    def apply_to(s):
        return not f(s)
    return apply_to


def find_all(words, case_sens=False):
    """Returns True if all given words appear in s.
    
    This filter is case INsensitive by default.
    """
    if not case_sens:
        used_words = [w.lower() for w in words]
    else:
        used_words = words

    def apply_to(s):
        if not case_sens:
            used_s = s.lower()
        else:
            used_s = s
        for w in used_words:
            # if w doesn't appear in lc_s
            if not used_s.find(w) > -1:
                return False
        return True
    return apply_to


def keep_if_more(items,x,case_sens=False):
    """Returns True if #items in s > x. False otherwise.
    
    This filter is case INsensitive by default.
    """
    
    x = int(x)
    if x < 0:
        raise IndexError, "x should be >= 0"
    
    if not case_sens:
        used_items = [str(item).lower() for item in items]
    else:
        used_items = items
    
    def find_more_good(s):
        if s and not case_sens:
            for i in s:
                used_s = [str(item).lower() for item in s]
        else:
            used_s = s
        fd = Freqs(used_s)
        value_list = [fd[i] for i in fd if i in used_items]
        if value_list:
            count = reduce(add, value_list)
            return count > x
        else:
            return False
    return find_more_good

def exclude_if_more(items,x,case_sens=False):
    """Returns True if #items in s < x.
    
    This filter is case INsensitive by default.
    """
    f = keep_if_more(items,x,case_sens)
    def apply_to(s):
        return not f(s)
    return apply_to

def keep_if_more_other(items,x,case_sens=False):
    """Returns True if #items in s other than those in items > x. 
    
    This filter is case INsensitive by default.
    """
    x = int(x)
    if x < 0:
        raise IndexError, "x should be >= 0"

    if not case_sens:
        used_items = [str(item).lower() for item in items]
    else:
        used_items = items
        
    def apply_to(s):
        if s and not case_sens:
            used_s = [str(item).lower() for item in s]
        else:
            used_s = s
        fd = Freqs(used_s)
        value_list = [fd[i] for i in fd if i not in used_items]
        if value_list:
            count = reduce(add, value_list)
            return count > x
        else:
            return False
    return apply_to

def exclude_if_more_other(items,x,case_sens=False):
    """Returns True if #items other than in items in s < x.
    
    This filter is case INsensitive by default.
    """ 
    f = keep_if_more_other(items,x,case_sens)
    def apply_to(s):
        return not f(s)
    return apply_to

'''
def keep_chars(keep, case_sens=True):
    """Returns a filter function f(s) that returns a filtered string.

    Specifically, strips out everything in s that is not in keep.
    This filter is case sensitive by default.
    """
    allchars = maketrans('', '')
    if not case_sens:
        low = keep.lower()
        up = keep.upper()
        keep = low + up
    delchars = ''.join([c for c in allchars if c not in keep])
    #make the filter function, capturing allchars and delchars in closure
    def filter_function(s, a=allchars, d=delchars):
        return s.translate(a, d)
    #return the filter function
    return filter_function
'''

class keep_chars(object):
    """Returns a filter object o(s): call to return a filtered string.

    Specifically, strips out everything in s that is not in keep.
    This filter is case sensitive by default.
    """
    allchars = maketrans('','')
    def __init__(self, keep, case_sens=True):
        """Returns a new keep_chars object, based on string keep"""
        if not case_sens:
            low = keep.lower()
            up = keep.upper()
            keep = low + up
        self.delchars = ''.join([c for c in self.allchars if c not in keep])
    
    def __call__(self, s, a=None, d=None):
        """f(s) -> s, translates using self.allchars and self.delchars"""
        if a is None: a = self.allchars
        if d is None: d = self.delchars
        return s.translate(a,d)

def exclude_chars(exclude,case_sens=True):
    """Returns a filter function f(s) that returns a filtered string.

    Specifically, strips out everything is s that is in exlude.
    This filter is case sensitive by default.
    """
    allchars = maketrans('','')
    if not case_sens:
        low = exclude.lower()
        up = exclude.upper()
        exclude = low + up
    delchars = ''.join([c for c in allchars if c in exclude])

    def filter_function(s,a=allchars, d=delchars):
        return s.translate(a,d)
    return filter_function

def reorder(order):
    """Returns a function that rearranges sequence into specified order.

    order should be a sequence of indices (for lists, tuples, strings, etc.)
    or keys (for dicts, etc.).

    Always returns a list.

    Will raise the appropriate IndexError or KeyError if items does not contain
    any position reqested by the order.
    """
    def result(items):
        return select(order, items)
    return result

def reorder_inplace(order, attr=None):
    """Returns a function that rearranges the items in attr, in place.

    If attr is None (the default), reorders the object itself.

    Uses slice assignment, so, unlike the original reorder, will only work on
    mutable sequences (e.g. lists, but not tuples or strings or dicts).
    """
    if attr:
        def result(obj):
            curr = getattr(obj, attr)
            curr[:] = select(order, curr)
            return obj
    else:
        def result(obj):
            obj[:] = select(order, obj)
            return obj
    return result

maybe_number = keep_chars('0123456789.+-eE')

def float_from_string(data):
    """Extracts a floating point number from string in data, if possible."""
    return float(maybe_number(data))

   
def first_index(seq, f):
    """Returns index of first item in seq where f(item) is True, or None.
    
    To invert the function, use lambda f: not f
    """
    for i, s in enumerate(seq):
        if f(s):
            return i

def last_index(seq, f):
    """Returns index of last item in seq where f(item) is True, or None.
    
    To invert the function, use lambda f: not f
    
    NOTE: We could do this slightly more efficiently by iterating over s in
    reverse order, but then it wouldn't work on generators that can't be
    reversed.
    """
    found = None
    for i, s in enumerate(seq):
        if f(s):
            found = i
    return found

def first_index_in_set(seq, items):
    """Returns index of first occurrence of any of items in seq, or None."""
    for i, s in enumerate(seq):
        if s in items:
            return i

def last_index_in_set(seq, items):
    """Returns index of last occurrence of any of items in seq, or None.
    
    NOTE: We could do this slightly more efficiently by iterating over s in
    reverse order, but then it wouldn't work on generators that can't be
    reversed.
    """
    found = None
    for i, s in enumerate(seq):
        if s in items:
            found = i
    return found

def first_index_not_in_set(seq, items):
    """Returns index of first occurrence of any of items in seq, or None."""
    for i, s in enumerate(seq):
        if not s in items:
            return i

def last_index_not_in_set(seq, items):
    """Returns index of last occurrence of any of items in seq, or None.
    
    NOTE: We could do this slightly more efficiently by iterating over s in
    reverse order, but then it wouldn't work on generators that can't be
    reversed.
    """
    found = None
    for i, s in enumerate(seq):
        if s not in items:
            found = i
    return found

def first(seq, f):
    """Returns first item in seq where f(item) is True, or None.
    
    To invert the function, use lambda f: not f
    """
    for s in seq:
        if f(s):
            return s

def last(seq, f):
    """Returns last item in seq where f(item) is True, or None.
    
    To invert the function, use lambda f: not f
    
    NOTE: We could do this slightly more efficiently by iterating over s in
    reverse order, but then it wouldn't work on generators that can't be
    reversed.
    """
    found = None
    for s in seq:
        if f(s):
            found = s
    return found

def first_in_set(seq, items):
    """Returns first occurrence of any of items in seq, or None."""
    for s in seq:
        if s in items:
            return s

def last_in_set(seq, items):
    """Returns index of last occurrence of any of items in seq, or None.
    
    NOTE: We could do this slightly more efficiently by iterating over s in
    reverse order, but then it wouldn't work on generators that can't be
    reversed.
    """
    found = None
    for s in seq:
        if s in items:
            found = s
    return found

def first_not_in_set(seq, items):
    """Returns first occurrence of any of items in seq, or None."""
    for s in seq:
        if not s in items:
            return s

def last_not_in_set(seq, items):
    """Returns last occurrence of any of items in seq, or None.
    
    NOTE: We could do this slightly more efficiently by iterating over s in
    reverse order, but then it wouldn't work on generators that can't be
    reversed.
    """
    found = None
    for s in seq:
        if s not in items:
            found = s
    return found

def perm(items, n=None):
    """Yields each successive permutation of items.

    This version from Raymond Hettinger, 2006/03/23
    """
    if n is None:
        n = len(items)
    for i in range(len(items)):
        v = items[i:i+1]
        if n == 1:
            yield v
        else:
            rest = items[:i] + items[i+1:]
            for p in perm(rest, n-1):
                yield v + p

def comb(items, n=None):
    """Yields each successive combination of n items.

    items: a slicable sequence.
    n: number of items in each combination
    This version from Raymond Hettinger, 2006/03/23
    """
    if n is None:
        n = len(items)
    for i in range(len(items)):
        v = items[i:i+1]
        if n == 1:
            yield v
        else:
            rest = items[i+1:]
            for c in comb(rest, n-1):
                yield v + c

def _increment_comb(outcomes, vector):
    """Yields each new outcome as an expansion of existing outcomes.
    """
    for outcome in outcomes:
        for e in vector:
            yield outcome + [e]

def cross_comb(vectors):
    """Yields each cross combination of a sequence of sequences (e.g. lists).

    i.e. returns the Cartesian product of the sequences.

    vectors: must be a seq of sequences.

    Recipe from the Python cookbook.

    Profiling shows that this is slightly slower than the previous
    implementation of cartesian_product for long lists, but faster for short
    lists. The speed penalty for long lists is outweighed by not having to
    keep the lists in memory.
    """
    result = ([],)
    for vector in vectors:
        result = _increment_comb(result, vector)
    return result

cartesian_product = cross_comb  #standard, but obscure, name for this function
